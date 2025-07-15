import os
import datetime
from concurrent.futures import ProcessPoolExecutor
import re
import ast
from abc import ABC, abstractmethod
from typing import Dict
import pandas as pd
import numpy as np
from pathlib import Path
from sqlalchemy import create_engine
from sqlalchemy import text
from sqlalchemy.pool import NullPool
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise_distances
from Bio.Align import PairwiseAligner, substitution_matrices
import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.express as px
import snf
from scipy.linalg import eigh
from itertools import product

from app_components.function_dump import (
    extractor,
    parser,
    calculate_cdr_chars,
    calculate_antigen_chars,
    to_list,
    prefix,
    clear_dir,
    hbdscan_cluster,
    DBCV,
    UMAP,
    merger_func
)

from app_components.config import INTERNAL_FILES_DIR, COMPUTATION_DIR, SQL_DIR

#Separation of all application functions into steps, which are brought together in a recipe like way, to accomplish different tasks; aka chem, clean, or extract from version 1 (=version lame) of this program.
#Object oriented approach, where dfs are the key object being traded between classes/functions, enables more accesible operation with different file formats. I cannot train my database solely with a former parser that reads only NAstructural database csvs.
#CLI is connected and feeds a list of steps depending on the type of recipe the user wants to use: ex. complete to do a full treatment of teh raw data, or simple for testing purposes. Inspect clauses are handled independently. I preferred inspection with manually created functions, because pandas is nicer than dunder methods. 

#New implementation, introducing a logger so that the streamlit app can work with multithreading. This class assumes the function of a visitor, logging the former print statements and storing them, so that the GUI can read them...
class FileLogger:
    def __init__(self, job_id: str, log_dir: Path = INTERNAL_FILES_DIR):
        log_dir.mkdir(exist_ok=True)
        self.filename = log_dir / f"{job_id}.log"
        self.file_handle = open(self.filename, 'w', encoding='utf-8')
    def log(self, message: str):
        log_entry = f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}\n"
        print(log_entry, end='')
        if self.file_handle and not self.file_handle.closed:
            self.file_handle.write(log_entry)
            self.file_handle.flush()
    def get_filename(self) -> str:
        return self.filename
    def close(self):
        if self.file_handle and not self.file_handle.closed:
            self.file_handle.close()

class Pipeline: #Receives a "recipe"/list of Step subclasses needed to produce a desired output
    def __init__(self, steps, logger: FileLogger):
        self.steps = steps
        self.logger = logger
    def run(self) -> Dict[str, pd.DataFrame]:
        data: Dict[str, pd.DataFrame] = {}
        for step in self.steps:
            data = step.process(data)
        return data

class Step(ABC): #Parent class to be method overrriden by specialised children
    def __init__(self, logger: FileLogger):
        self.logger = logger
    @abstractmethod
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        pass

class Concatenation(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        all_cdr_frames = [pd.read_csv(f) for f in INTERNAL_FILES_DIR.glob('*_cdr.csv')]
        all_antigen_frames = [pd.read_csv(f) for f in INTERNAL_FILES_DIR.glob('*_antigen.csv')]

        final_data: Dict[str, pd.DataFrame] = {}

        if all_cdr_frames:
            final_data['cdr'] = pd.concat(all_cdr_frames, ignore_index=True)
            self.logger.log(f"Concatenation → Final CDR DataFrame has {len(final_data['cdr'])} rows.")
        else:
            final_data['cdr'] = pd.DataFrame()
            self.logger.log("Concatenation → No CDR data found.")
        
        if all_antigen_frames:
            final_data['antigen'] = pd.concat(all_antigen_frames, ignore_index=True)
            self.logger.log(f"Concatenation → Final Antigen DataFrame has {len(final_data['antigen'])} rows.")
        else:
            final_data['antigen'] = pd.DataFrame()
            self.logger.log("Concatenation → No Antigen data found.")
            
        return final_data

class Write(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        COMPUTATION_DIR.mkdir(exist_ok=True)
        for key, df in data.items():
            filename = f"{key}.csv"
            filepath = COMPUTATION_DIR / filename
            df.to_csv(filepath, index=False)
            self.logger.log(f"Wrote DataFrame '{key}' to {filepath}")
        self.logger.log(f"Write → wrote {len(data)} tables: {list(data.keys())}")
        return {}

#CHEMICAL CHARACTERISTICS

class CDRComputation(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data = data.copy()
        if 'cdr' in new_data and not new_data['cdr'].empty:
            self.logger.log("CDRComputation → Starting...")
            new_data['cdr'] = calculate_cdr_chars(new_data['cdr'])
            self.logger.log(f"CDRComputation → Processed 'cdr', final shape: {new_data['cdr'].shape}")
        return new_data

class AntigenComputation(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data = data.copy()
        if 'antigen' in new_data and not new_data['antigen'].empty:
            self.logger.log("AntigenComputation → Starting...")
            new_data['antigen'] = calculate_antigen_chars(new_data['antigen'])
            self.logger.log(f"AntigenComputation → Processed 'antigen', final shape: {new_data['antigen'].shape}")
        return new_data
#CHEMICAL CHARACTERISTICS

# Multiple antibodies can bind to a single antigen, many antigens can sature a single antibody, this is reflected in the source data too...
class FlattenDuplicates(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        if 'cdr' in data:
            cdr_df = data['cdr'].copy()
            if not cdr_df['h3_chain'].duplicated().any():
                self.logger.log(f"FlattenDuplicates → no CDR duplicates found, keeping {len(cdr_df)} rows")
                new_data['cdr'] = cdr_df
            else:
                search_dict = {
                    col: 'first'
                    for col in cdr_df.columns
                    if col not in ('h3_chain', 'pdb_id')
                }
                search_dict['pdb_id'] = lambda ids: list(ids.unique())
                flat_cdr = (
                    cdr_df
                    .groupby('h3_chain')
                    .agg(search_dict)
                    .reset_index()
                )
                new_data['cdr'] = flat_cdr
        if 'antigen' in data:
            df = data['antigen'].copy()
            if not df['antigen_seq'].duplicated().any():
                self.logger.log(f"FlattenDuplicates → no Antigen duplicates found, keeping {len(df)} rows")
                new_data['antigen'] = df
            else:
                agg_dict = {
                    col: 'first'
                    for col in df.columns
                    if col not in ('antigen_seq', 'corresponding_pdb_antibody')
                }
                agg_dict['corresponding_pdb_antibody'] = lambda ids: list(ids.unique())
                flat_df = (
                    df
                    .groupby('antigen_seq')
                    .agg(agg_dict)
                    .reset_index()
                )
                new_data['antigen'] = flat_df
        if 'antigen' in new_data:
            self.logger.log(f"FlattenDuplicates → flattened antigen, rows: {new_data['antigen'].shape[0]}")
        if 'cdr' in new_data:
            self.logger.log(f"FlattenDuplicates → flattened CDR, rows: {new_data['cdr'].shape[0]}")
        return new_data

class CleanUp(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        current_log_file = os.path.basename(self.logger.get_filename())
        clear_dir(str(INTERNAL_FILES_DIR), exception=current_log_file)
        clear_dir(str(COMPUTATION_DIR))
        return {}

class PreWalker(Step):
    def __init__(self, input_path: str, logger: FileLogger):
        super().__init__(logger)
        self.input_path = Path(input_path)

    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        if not self.input_path.is_dir():
            self.logger.log(f"PreWalker → Input path '{self.input_path}' is not a valid directory.")
            data["paths"] = []
            return data

        processed_stems = {f.stem.replace('_antigen', '').replace('_cdr', '') for f in INTERNAL_FILES_DIR.glob('*.csv')}
        self.logger.log(f"PreWalker → Found {len(processed_stems)} stems from previously processed files.")

        all_csv_files_in_path = {f.resolve() for f in self.input_path.rglob("*.csv")}

        new_source_files = []
        for file in all_csv_files_in_path:

            if file.name.endswith(('_cdr.csv', '_antigen.csv')):
                continue

            if file.stem in processed_stems:
                continue

            new_source_files.append(str(file))

        data["paths"] = new_source_files
        self.logger.log(f"PreWalker → Found {len(new_source_files)} new source files to process.")
        return data

#Re wrote the following piece of shit, it was causing me headaches before
class Walker(Step):
    @staticmethod
    def _process_file(args):
        filepath, idx = args
        df_raw = extractor(filepath)
        antigen_df, cdr_df = parser(df_raw)
        return filepath, idx, antigen_df, cdr_df
        
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
            tasks = []
            INTERNAL_FILES_DIR.mkdir(exist_ok=True)
            for idx, filepath_str in enumerate(data.get("paths", []), start=1):
                if filepath_str.endswith('_cdr.csv') or filepath_str.endswith('_antigen.csv'):
                    self.logger.log(f"Walker → Skipping already-parsed file: {os.path.basename(filepath_str)}")
                    continue
                tasks.append((filepath_str, idx))
            if not tasks:
                self.logger.log("Walker → No new files to process, HEADACHE")
                return {} 
            with ProcessPoolExecutor(max_workers=2) as executor:
                for filepath, idx, antigen_df, cdr_df in executor.map(self._process_file, tasks):
                    base_key = Path(filepath).stem
                    if not antigen_df.empty:
                        filename_agn = f"{base_key}_antigen.csv"
                        filepath_agn = INTERNAL_FILES_DIR / filename_agn
                        antigen_df.to_csv(filepath_agn, index=False)
                        self.logger.log(f"Wrote Antigen DataFrame from '{base_key}' to {filepath_agn}")
                    if not cdr_df.empty:
                        filename_abd = f"{base_key}_cdr.csv"
                        filepath_abd = INTERNAL_FILES_DIR / filename_abd
                        cdr_df.to_csv(filepath_abd, index=False)
                        self.logger.log(f"Wrote CDR DataFrame from '{base_key}' to {filepath_abd}")
                    self.logger.log(
                        f"Walker → processed {filepath!r}: "
                        f"antigen rows={antigen_df.shape[0]}, "
                        f"cdr rows={cdr_df.shape[0]}"
                    )
            return {}

#Offering a simple to implement solution to a complicated issue: remove purification tags from antigenic sequences no matter if they are alone or in groups or whatever.
class RmPurificationTags(Step):
    def __init__(self, logger: FileLogger):
        super().__init__(logger)
        self.motifs = {
            "c-Myc": "EQKLISEEDL",
            "6-His": "HHHHHH",
            "5-His": "HHHHH",
            "FLAG":  "DYKDDDDK",
            "V5":    "GKPIPNPLLGLDST",
        }
        self.max_gap = 6
        self.alteration = "|".join(map(re.escape, self.motifs.values()))
        self.clusters = re.compile(
            rf'(?:{self.alteration}(?:[A-Z]{{1,{self.max_gap}}}{self.alteration})*)'
        )
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        if "cdr" in data:
            new_data["cdr"] = data["cdr"].copy()
        if "antigen" in data:
            df = data["antigen"].copy()
            for idx, seq in df["antigen_seq"].items():
                df.at[idx, "antigen_seq"] = self.clusters.sub("", seq)
            new_data["antigen"] = df
            self.logger.log(f"RmPurificationTags → cleaned antigenic sequences in the following number of rows: {new_data['antigen'].shape[0]}")
        return new_data

#Remember, this step is ONLY intended for assigning unique identifiers, this is uninformative encoding, all prior information is lost dramatically; the sequence cannot possibly be discerned from the computed ID. 
class AssignIDs(Step):
    def __init__(self, logger: FileLogger):
        super().__init__(logger)
        self.amino_acid_rubric = {
        'A': 1,  'C': 2,  'D': 3,  'E': 4,  'F': 5,
        'G': 6,  'H': 7,  'I': 8,  'K': 9,  'L':10,
        'M':11,  'N':12,  'P':13,  'Q':14,  'R':15,
        'S':16,  'T':17,  'V':18,  'W':19,  'Y':20,}
        self.re_pattern = re.compile(r'[ACDEFGHIKLMNPQRSTVWY]')
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        if 'cdr' in data:
            cdr_df = data["cdr"].copy()
            cdr_df["cdr_computed_id"] = pd.NA
            for idx, seq in cdr_df["h3_chain"].items():
                if re.search(r'[^ACDEFGHIKLMNPQRSTVWY]', seq):
                    raise ValueError(f"Illegal residue in {seq!r} at row {idx}")
                else:
                    cdr_df.at[idx, "cdr_computed_id"] = self.re_pattern.sub(lambda m: str(self.amino_acid_rubric[m.group(0)]), seq)
            new_data["cdr"] = cdr_df
            self.logger.log(f"AssignIDs → assigned {new_data['cdr'].shape[0]} sequence based IDs")
        if "antigen" in data:
            ant_df = data["antigen"].copy()
            ant_df["antigen_computed_id"] = pd.NA
            for idx, seq in ant_df["antigen_seq"].items():
                if re.search(r'[^ACDEFGHIKLMNPQRSTVWY]', seq):
                    raise ValueError(f"Illegal residue in {seq!r} at row {idx}")
                else:
                    ant_df.at[idx, "antigen_computed_id"] = self.re_pattern.sub(lambda m: str(self.amino_acid_rubric[m.group(0)]), seq)
            new_data["antigen"] = ant_df
            self.logger.log(f"AssignIDs → assigned {new_data['antigen'].shape[0]} sequence based IDs")
        return new_data


#Bear in mind this algorithm prefers the use of dictionaries to compact data; might not be entirely sound..
#Importantly, one could use pandas to make use of the C based hash join?! However I am more comfortable with dictionaries and feel that an inner join with pandas might miss out on the multiplicity of connections between different antigen/cdr objects... this is entirely based on vibes
#Algorithmic complexity is therefore higher than pandas but whatever...preliminary trials showed that to parse throuh the entire NAStructuralDB deduplicated dataset and find 679 relationships, it takes 5 seconds on M1 Macintosh without connection to power supply!
class ComputeRelationships(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        antigens = data["antigen"].copy()
        cdrs = data["cdr"].copy()
        antigens["pdb_list"] = antigens["corresponding_pdb_antibody"].apply(to_list)
        cdrs["pdb_list"]     = cdrs["pdb_id"].apply(to_list)
        ag_exp  = antigens.explode("pdb_list")[["antigen_computed_id", "pdb_list"]]
        cdr_exp = cdrs.explode("pdb_list")[["cdr_computed_id",     "pdb_list"]]
        pairs = (
            pd.merge(ag_exp, cdr_exp, on="pdb_list", how="inner")
              [["antigen_computed_id", "cdr_computed_id"]]
              .drop_duplicates()
              .reset_index(drop=True)
        )
        new_data["antigen"]       = antigens.drop(columns=["pdb_list"])
        new_data["cdr"]           = cdrs.drop(columns=["pdb_list"])
        new_data["relationships"] = pairs
        return new_data
    
class WorkWithDatabase(Step):
    def __init__(self, logger: FileLogger):
        super().__init__(logger)
        self.connection_string = "postgresql://chrismitsacopoulos:password@localhost/pasteurdb"
        self.ddl_dir = SQL_DIR
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        ddl_paths = sorted(self.ddl_dir.glob("*.sql"), key=prefix)
        if not ddl_paths:
            raise RuntimeError(f"No DDL files found in {self.ddl_dir}")
        engine = create_engine(self.connection_string, echo=True, poolclass=NullPool, pool_pre_ping=True)
        with engine.begin() as conn:
            first = ddl_paths[0]
            self.logger.log(f"WorkWithDatabase → creating staging tables from {first.name}")
            for stmt in first.read_text(encoding="utf-8").split(";"):
                stmt = stmt.strip()
                if stmt and not re.fullmatch(r"(?i)(BEGIN|COMMIT|ROLLBACK)", stmt):
                    conn.execute(text(stmt))
            self.logger.log("WorkWithDatabase → inserting DataFrames into staging tables...")
            for dict_key, dict_df in data.items():
                if not dict_df.empty:
                    stg_table = f"staging_{dict_key}"
                    self.logger.log(f" • loading {len(dict_df)} rows into {stg_table}")
                    dict_df.to_sql(stg_table, conn, if_exists="append", index=False)
            self.logger.log("WorkWithDatabase → altering staging_antigen.antigen_is_incomplete to BOOLEAN...")
            conn.execute(text(      """
                ALTER TABLE staging_antigen
                  ALTER COLUMN antigen_is_incomplete TYPE BOOLEAN
                  USING (antigen_is_incomplete = 1);
                """)
            )
            self.logger.log("WorkWithDatabase → altering staging_cdr.h3_is_incomplete to BOOLEAN...")
            conn.execute(text(         """
                ALTER TABLE staging_cdr
                  ALTER COLUMN h3_is_incomplete TYPE BOOLEAN
                  USING (h3_is_incomplete = 1);
                """)
            )
            self.logger.log("WorkWithDatabase → altering staging_cdr.l3_is_incomplete to BOOLEAN...")
            conn.execute(text( """
                ALTER TABLE staging_cdr
                  ALTER COLUMN l3_is_incomplete TYPE BOOLEAN
                  USING (l3_is_incomplete = 1);
                """)
            )
            for ddl_path in ddl_paths[1:]:
                self.logger.log(f"WorkWithDatabase → applying post-load DDL {ddl_path.name}")
                for stmt in ddl_path.read_text(encoding="utf-8").split(";"):
                    stmt = stmt.strip()
                    if stmt and not re.fullmatch(r"(?i)(BEGIN|COMMIT|ROLLBACK)", stmt):
                        conn.execute(text(stmt))
        self.logger.log("WorkWithDatabase → all DDL applied and data loaded.")
        return data
    
class ComputeDistanceMatrices(Step):
    """
    Explanation: BLOSUM45 due to hypervariability of CDR3H sequences, this will be tweaked with BLOSUM62 as a substitute to observe differences in analysis outcome; regardless studies demonstrate preference for BLOSUM45 when handling CDR. Additionally, to make this pairwise (BLOSUM) similarity matrix useful, just like with CTD stats, values are normalised to eliminate strong outlier biases but for BLOSUM, the max normalised score of the matrix is subtracted against by all smaller values. As such, max similarity will have a distance of zero, while all others...since BLOSUM pairwise distances are computed pairwise in an NxN matrix, there will always be  of pairwise comparison on the same same sequence for this methodology to be just. 
    """
    def _calculate_blosum_score(self, seq1: str, seq2: str) -> float:
        aligner = PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM45")
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -1.0
        try:
            return aligner.score(seq1, seq2)
        except Exception:
            return 0.0
    """
    Im keeping the same data dict retrieval despite its non-existence when starting up the recipe, purely because I know it works; bad practice nonetheless...
    """
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        self.logger.log("ComputeDistanceMatrices → Starting...")
        data = {}
        path_to_csv = COMPUTATION_DIR / "cdr.csv"
        data["cdr"] = pd.read_csv(path_to_csv)
        if 'cdr' not in data or data['cdr'].empty:
            self.logger.log("ComputeDistanceMatrices → SKIPPING: 'cdr' DataFrame not found or is empty.")
            return data

        cdr_df = data['cdr']
        ctd_path = INTERNAL_FILES_DIR / 'ctd_dist_matrix.npz'
        blosum_path = INTERNAL_FILES_DIR / 'blosum_dist_matrix.npz'

        if 'h3_ctd' not in cdr_df.columns:
            self.logger.log("ComputeDistanceMatrices → SKIPPING CTD: 'h3_ctd' column not found.")
        elif ctd_path.exists():
            self.logger.log(f"ComputeDistanceMatrices → Found cached CTD matrix. Loading from {ctd_path}.")
            loaded_data = np.load(ctd_path)
            data['ctd_dist_matrix'] = loaded_data['matrix']
        else:
            self.logger.log("ComputeDistanceMatrices → No CTD cache. Processing distances...")
            try:
                sequences = cdr_df["h3_chain"].to_numpy(dtype=object)
                ctd_vectors = cdr_df['h3_ctd'].apply(ast.literal_eval).tolist()
                feature_matrix = np.array(ctd_vectors)
                scaler = StandardScaler()
                scaled_features = scaler.fit_transform(feature_matrix)
                ctd_dist_matrix = pairwise_distances(scaled_features, metric='manhattan')
                np.savez(ctd_path, matrix=ctd_dist_matrix, sequences=sequences)
                self.logger.log(f"ComputeDistanceMatrices → Saved new CTD distance matrix to {ctd_path}.")
                data['ctd_dist_matrix'] = ctd_dist_matrix
            except Exception as e:
                self.logger.log(f"ComputeDistanceMatrices → ERROR processing CTD distances: {e}")
        if blosum_path.exists():
            self.logger.log(f"ComputeDistanceMatrices → Found cached BLOSUM matrix. Loading from {blosum_path}.")
            loaded_data = np.load(blosum_path)
            data['blosum_dist_matrix'] = loaded_data['matrix']
        else:
            self.logger.log("ComputeDistanceMatrices → No BLOSUM cache. Processing distances...")
            try:
                sequences = cdr_df['h3_chain'].tolist()
                num_sequences = len(sequences)
                blosum_similarity_matrix = np.zeros((num_sequences, num_sequences))
                for i in range(num_sequences):
                    for j in range(i, num_sequences):
                        score = self._calculate_blosum_score(sequences[i], sequences[j])
                        blosum_similarity_matrix[i, j] = blosum_similarity_matrix[j, i] = score
                
                max_blosum_score = np.max(blosum_similarity_matrix)
                blosum_dist_matrix = max_blosum_score - blosum_similarity_matrix
                np.savez(blosum_path, matrix=blosum_dist_matrix, sequences=np.array(sequences, dtype=object))
                self.logger.log(f"ComputeDistanceMatrices → Saved new BLOSUM distance matrix to {blosum_path}.")
                data['blosum_dist_matrix'] = blosum_dist_matrix
            except Exception as e:
                self.logger.log(f"ComputeDistanceMatrices → ERROR processing BLOSUM distances: {e}")

        self.logger.log("ComputeDistanceMatrices → Finished.")
        return data

class TuneSNFParameters(Step):
    def __init__(self, logger: FileLogger, k_range: list = [10, 20, 30], mu_range: list = [0.3, 0.5, 0.8]):
        super().__init__(logger)
        self.k_range = k_range
        self.mu_range = mu_range

    def _calculate_eigengap(self, affinity_matrix: np.ndarray, num_eig: int = 10) -> float:
        diag = np.sum(affinity_matrix, axis=1)
        diag[diag == 0] = 1e-9 
        D_inv_sqrt = np.diag(1.0 / np.sqrt(diag))
        L = np.eye(affinity_matrix.shape[0]) - D_inv_sqrt @ affinity_matrix @ D_inv_sqrt
        eigenvalues = eigh(L, eigvals_only=True, subset_by_index=[0, num_eig-1])
        eigenvalues.sort()
        return np.max(np.diff(eigenvalues))

    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        self.logger.log("TuneSNFParameters → Starting hyperparameter tuning for SNF.")
        if 'ctd_dist_matrix' not in data or 'blosum_dist_matrix' not in data:
            self.logger.log("TuneSNFParameters → SKIPPING: Required distance matrices not found.")
            return data
        best_score = -1
        best_params = {}
        param_grid = list(product(self.k_range, self.mu_range))
        self.logger.log(f"TuneSNFParameters → Testing {len(param_grid)} parameter combinations...")
        for k, mu in param_grid:
            affinity_ctd = snf.compute.affinity_matrix(data['ctd_dist_matrix'], K=k, mu=mu)
            affinity_blosum = snf.compute.affinity_matrix(data['blosum_dist_matrix'], K=k, mu=mu)
            fused_matrix = snf.snf([affinity_ctd, affinity_blosum], K=k)
            score = self._calculate_eigengap(fused_matrix)
            self.logger.log(f"TuneSNFParameters → Tested K={k}, mu={mu}. Score (Eigen-gap): {score:.4f}")
            if score > best_score:
                best_score = score
                best_params = {'K': k, 'mu': mu}
        
        self.logger.log(f"TuneSNFParameters → Tuning complete. Best parameters found: {best_params} with score {best_score:.4f}")
        data['best_snf_params'] = best_params
        return data
    
class FuseAndProject(Step):
    def __init__(self, logger: FileLogger, k_neighbors: int = 20, mu_param: float = 0.5):
        super().__init__(logger)
        self.k_neighbors = k_neighbors
        self.mu_param = mu_param
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        self.logger.log("FuseAndProject → Starting network fusion and projection.")

        if 'best_snf_params' in data:
            self.logger.log(f"FuseAndProject → Using tuned SNF parameters: {data['best_snf_params']}")
            self.k_neighbors = data['best_snf_params']['K']
            self.mu_param = data['best_snf_params']['mu']
        else:
            self.logger.log(f"FuseAndProject → Using default SNF parameters: K={self.k_neighbors}, mu={self.mu_param}")

        if 'ctd_dist_matrix' not in data or 'blosum_dist_matrix' not in data:
            self.logger.log("FuseAndProject → SKIPPING: Required distance matrices not found in data dictionary.")
            return data
        fused_matrix_path = INTERNAL_FILES_DIR / 'fused_affinity_matrix.npz'
        plot_filepath = INTERNAL_FILES_DIR / 'snf_umap_projection.html'

        if fused_matrix_path.exists():
            self.logger.log(f"FuseAndProject → Found cached fused matrix. Loading from {fused_matrix_path}.")
            loaded_data = np.load(fused_matrix_path)
            fused_affinity_matrix = loaded_data['matrix']
        else:
            self.logger.log("FuseAndProject → No cache found. Creating and fusing affinity matrices...")
            ctd_dist = data['ctd_dist_matrix']
            blosum_dist = data['blosum_dist_matrix']

            affinity_ctd = snf.compute.affinity_matrix(ctd_dist, K=self.k_neighbors, mu=self.mu_param)
            affinity_blosum = snf.compute.affinity_matrix(blosum_dist, K=self.k_neighbors, mu=self.mu_param)
            
            fused_affinity_matrix = snf.snf([affinity_ctd, affinity_blosum], K=self.k_neighbors)
            np.savez(fused_matrix_path, matrix=fused_affinity_matrix)
            self.logger.log(f"FuseAndProject → Saved new fused affinity matrix to {fused_matrix_path}.")
    
        data['fused_affinity_matrix'] = fused_affinity_matrix
        self.logger.log("FuseAndProject → Added 'fused_affinity_matrix' to data dictionary for clustering.")

        self.logger.log("FuseAndProject → Generating UMAP embedding from fused network...")
        fused_distance_matrix = 1 - fused_affinity_matrix
        coords_2d = UMAP(fused_distance_matrix) 

        embedding_df = pd.DataFrame({
            'h3_chain': data['cdr']['h3_chain'],
            'umap_x': coords_2d[:, 0],
            'umap_y': coords_2d[:, 1]
        })

        merged_df = merger_func(embedding_df, data["cdr"]) 
        data['snf_embedding'] = merged_df
        self.logger.log(f"FuseAndProject → Complete. Added 'snf_embedding' DataFrame with {len(merged_df)} rows.")

        self.logger.log("FuseAndProject → Generating interactive plot...")
        fig = px.scatter(
            merged_df, 
            x='umap_x',
            y='umap_y',
            hover_name='h3_chain',
            hover_data=['pdb_id', 'heavy_host_organism_name'], 
            title='UMAP Projection of Fused CTD and BLOSUM Networks'
        )
        fig.update_traces(marker=dict(size=5, opacity=0.8))
        fig.update_layout(xaxis_title="UMAP Dimension 1", yaxis_title="UMAP Dimension 2")
        fig.write_html(plot_filepath)
        self.logger.log(f"FuseAndProject → Saved interactive plot to {plot_filepath}")
        
        return data
    
class GetDensity(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        self.logger.log("GetDensity → Starting cluster stability analysis.")
        path_to_csv = INTERNAL_FILES_DIR / "fused_affinity_matrix.npz"
        npz_file = np.load(path_to_csv)

        self.logger.log("GetDensity → Analysing stability on the high-dimensional fused network...")
        fused_affinity_matrix = npz_file['matrix']
        fused_distance_matrix = 1 - fused_affinity_matrix
        stability_plot_path = INTERNAL_FILES_DIR / 'stability_plot.html'

        stability_df = DBCV(fused_distance_matrix)
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=stability_df['min_cluster_size'],
            y=stability_df['silhouette_score'], 
            mode='lines+markers',
            name='Silhouette Score'
        ))
        fig.update_layout(
            title='HDBSCAN Cluster Stability (Higher is Better)',
            xaxis_title='Minimum Cluster Size',
            yaxis_title='Silhouette Score' 
        )
        fig.write_html(stability_plot_path)
        self.logger.log(f"GetDensity → Saved interactive stability plot to {stability_plot_path}")

        best_size = stability_df.loc[stability_df['silhouette_score'].idxmax()]
        self.logger.log(f"GetDensity → Analysis complete. Optimal min_cluster_size found: {int(best_size['min_cluster_size'])}")
        
        data['cluster_stability'] = stability_df
        data['best_min_cluster_size'] = int(best_size['min_cluster_size'])
        return data

class HDBSCAN(Step):
    def __init__(self, logger: FileLogger, min_cluster_size: int):
        super().__init__(logger)
        if not isinstance(min_cluster_size, int) or min_cluster_size < 1:
            raise ValueError("min_cluster_size must be an integer greater than 1.")
        self.min_cluster_size = min_cluster_size

    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data = data.copy()
        if 'umap_embedding' not in new_data:
            self.logger.log("HDBSCAN → Skipping: 'umap_embedding' DataFrame not found.")
            return new_data

        self.logger.log(f"HDBSCAN → Clustering with min_cluster_size={self.min_cluster_size}")
        embedding_df = new_data['umap_embedding'].copy()
        coords_2d = embedding_df[['umap_x', 'umap_y']].to_numpy()

        cluster_labels = hbdscan_cluster(coords_2d, self.min_cluster_size)
        embedding_df['cluster_id'] = cluster_labels
        new_data['umap_clusters'] = embedding_df

        interactive_plot_path = os.path.join(INTERNAL_FILES_DIR, 'cluster_plot.html')
        
        fig = px.scatter(
            embedding_df,
            x='umap_x',
            y='umap_y',
            color=embedding_df['cluster_id'].astype(str),
            hover_name='h3_chain',
            hover_data= ['pdb_id', 'heavy_host_organism_name'],
            color_discrete_map={'-1': 'lightgrey'},
            title=f'HDBSCAN Clustering (min_cluster_size={self.min_cluster_size})'
        )
        fig.update_traces(marker=dict(size=5))
        fig.update_layout(xaxis_title="UMAP 1", yaxis_title="UMAP 2")
        fig.write_html(interactive_plot_path)

        self.logger.log(f"HDBSCAN → Saved interactive plot to {interactive_plot_path}")
        n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
        self.logger.log(f"HDBSCAN → Clustering complete. Found {n_clusters} clusters.")
        
        if 'umap_embedding' in new_data:
             del new_data['umap_embedding']

        return new_data