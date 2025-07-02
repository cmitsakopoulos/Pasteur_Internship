import os
import datetime
import time
from concurrent.futures import ProcessPoolExecutor
import re
from abc import ABC, abstractmethod
from typing import Dict
import pandas as pd
import numpy as np
from pathlib import Path
from sqlalchemy import create_engine
from sqlalchemy import text
from sqlalchemy.pool import NullPool
from rapidfuzz import process, distance
import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.express as px

from function_dump import (
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

#Separation of all application functions into steps, which are brought together in a recipe like way, to accomplish different tasks; aka chem, clean, or extract from version 1 (=version lame) of this program.
#Object oriented approach, where dfs are the key object being traded between classes/functions, enables more accesible operation with different file formats. I cannot train my database solely with a former parser that reads only NAstructural database csvs.
#CLI is connected and feeds a list of steps depending on the type of recipe the user wants to use: ex. complete to do a full treatment of teh raw data, or simple for testing purposes. Inspect clauses are handled independently. I preferred inspection with manually created functions, because pandas is nicer than dunder methods. 

#New implementation, introducing a logger so that the streamlit app can work with multithreading. This class assumes the function of a visitor, logging the former print statements and storing them, so that the GUI can read them...
class FileLogger:
    def __init__(self, job_id: str, log_dir: str = "Internal_Files"):
        os.makedirs(log_dir, exist_ok=True)
        self.filename = os.path.abspath(os.path.join(log_dir, f"{job_id}.log"))
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
        cdr_frames = []
        antigen_frames = []
        for key, df in data.items():
            if key.startswith('cdr_'):
                cdr_frames.append(df)
            elif key.startswith('antigen_'):
                antigen_frames.append(df)
        result: Dict[str, pd.DataFrame] = {}
        if cdr_frames:
            result['cdr'] = pd.concat(cdr_frames, ignore_index=True)
        else:
            result['cdr'] = pd.DataFrame()
            self.logger.log("No CDRS, moving on")
        if antigen_frames:
            result['antigen'] = pd.concat(antigen_frames, ignore_index=True)
        else:
            result['antigen'] = pd.DataFrame()
            self.logger.log("No Antigens, moving on")
        self.logger.log(f"Concatenation → output keys: {list(result.keys())}")
        return result

class Write(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        output_dir = os.path.join(os.path.dirname(__file__), "Computation_Deposit")
        os.makedirs(output_dir, exist_ok=True)
        for key, df in data.items():
            filename = f"{key}.csv"
            filepath = os.path.join(output_dir, filename)
            df.to_csv(filepath, index=False)
            self.logger.log(f"Wrote DataFrame '{key}' to {filepath}")
        self.logger.log(f"Write → wrote {len(data)} tables: {list(data.keys())}")
        return {}

#CHEMICAL CHARACTERISTICS
class CDRComputation(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = data.copy()
        if 'cdr' in new_data:
            computation_cdr_df = new_data["cdr"]
            calculate_cdr_chars(computation_cdr_df)
            new_data['cdr'] = computation_cdr_df
        if 'cdr' in new_data:
            self.logger.log(f"CDRComputation → processed cdr, rows: {new_data['cdr'].shape[0]}")
        return new_data

class AntigenComputation(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = data.copy()
        if 'antigen' in data:
            computation_antigen_df = new_data["antigen"]
            calculate_antigen_chars(computation_antigen_df)
            new_data['antigen'] = computation_antigen_df
        if 'antigen' in new_data:
            self.logger.log(f"AntigenComputation → processed antigen, rows: {new_data['antigen'].shape[0]}")
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

class CleanUp(Step): #Connected to function in function dump, just cleans up the code here
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        current_log_file = self.logger.get_filename()
        exception =os.path.basename(current_log_file)
        internal_dir = os.path.join(os.path.dirname(__file__), "Internal_Files")
        output_dir = os.path.join(os.path.dirname(__file__), "Computation_Deposit")
        clear_dir(internal_dir, exception)
        clear_dir(output_dir)
        data = {}
        return data

class PreWalker(Step): #Not incredibly efficient but very flexible method to check if files have already been parsed/extracted...it will also read files within the common directory if needed.
    def __init__(self, input_path: str, logger: FileLogger):
        super().__init__(logger) 
        self.input_path = Path(input_path)
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        internal_dir = Path(__file__).parent / "Internal_Files"
        input_files = list(self.input_path.rglob("*.csv"))
        if internal_dir.is_dir():
            processed_stems = {p.stem for p in internal_dir.rglob("*.csv")}
            tasks_to_be = [
                str(f.resolve())
                for f in input_files
                if f.stem not in processed_stems
            ]
            dropped = len(list(internal_dir.glob('*.csv')))
            data["paths"] = tasks_to_be
            self.logger.log(f"Dropped {dropped} already-processed files; {len(tasks_to_be)} remain.")
            for p in internal_dir.rglob("*.csv"):
                stem = p.stem.lower()
                df = pd.read_csv(p)
                if "antigen" in stem:
                    data[f"antigen_{stem}"] = df
                elif "cdr" in stem:
                    data[f"cdr_{stem}"] = df
        else:
            data["paths"] = [str(f.resolve()) for f in input_files]
        return data
    
#Named after the os import function "walk", traverses user provided directory to build our shared dict, which is then parsed further...
class Walker(Step):
    @staticmethod
    def _process_file(args):
        filepath, idx = args
        df_raw = extractor(filepath)
        antigen_df, cdr_df = parser(df_raw)
        return filepath, idx, antigen_df, cdr_df
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        tasks = []
        output_dir_pf = os.path.join(os.path.dirname(__file__), "Internal_Files")
        os.makedirs(output_dir_pf, exist_ok=True)
        for idx, filepath in enumerate(data["paths"], start=1):
            tasks.append((filepath, idx))
        with ProcessPoolExecutor(max_workers=2) as executor:
            for filepath, idx, antigen_df, cdr_df in executor.map(self._process_file, tasks):
                base_key = Path(filepath).stem
                if not antigen_df.empty:
                    new_data[f"antigen_{base_key}"] = antigen_df
                    filename_agn = f"{base_key}_antigen.csv"
                    filepath_agn = os.path.join(output_dir_pf, filename_agn)
                    antigen_df.to_csv(filepath_agn, index=False)
                    self.logger.log(f"Wrote Antigen DataFrame from '{base_key}' to {filepath_agn}")
                if not cdr_df.empty:
                    new_data[f"cdr_{base_key}"] = cdr_df
                    filename_abd = f"{base_key}_cdr.csv"
                    filepath_abd = os.path.join(output_dir_pf, filename_abd)
                    cdr_df.to_csv(filepath_abd, index=False)
                    self.logger.log(f"Wrote CDR DataFrame from '{base_key}' to {filepath_abd}")
                self.logger.log(
                    f"Walker → processed {filepath!r}: "
                    f"antigen rows={antigen_df.shape[0]}, "
                    f"cdr rows={cdr_df.shape[0]}"
                )
        for key, df in data.items(): #Without it the dataframes were being nuked...
            if key != 'paths' and key not in new_data:
                new_data[key] = df
        return new_data

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
        self.ddl_dir = Path(__file__).resolve().parent / 'sql_files'
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
    
class LevenshteinDistance(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        output_dir = os.path.join(os.path.dirname(__file__), "Internal_Files")
        os.makedirs(output_dir, exist_ok=True)
        output_filepath = os.path.join(output_dir, 'h3_chain_dist_matrix.npz')
        dataframe = data["cdr"]
        h3_sequences = dataframe["h3_chain"].tolist()
        self.logger.log(f"Extracted {len(h3_sequences)} unique sequences from 'h3_chain' column.")
        self.logger.log("Calculating pairwise distance matrix...")
        start_time = time.time()
        distance_matrix = process.cdist(
            h3_sequences,
            h3_sequences,
            scorer=distance.Levenshtein.distance,
            workers=-1
        )
        end_time = time.time()
        self.logger.log(f"Matrix calculation finished in {end_time - start_time:.4f} seconds.")
        self.logger.log(f"Shape of the distance matrix: {distance_matrix.shape}")
        self.logger.log(f"Saving matrix to '{output_filepath}'...")
        np.savez(output_filepath, matrix=distance_matrix, sequences=np.array(h3_sequences, dtype=object))
        self.logger.log("Save complete.")
        return data

#Note I store the UMAP embedding on the distance matrix as parquet to avoid PreWalker conflicts with csv.
#I prefer the user has to issue different recipes for each of the separate functions, as keep in mind I am working with 1.2k sequences, someone with 100k is ducked, the system might crash if it keeps resources too high for too long (running one step after the other nonstop).
class GenerateUMAP(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data = data.copy()
        internal_dir = os.path.join(os.path.dirname(__file__), "Internal_Files")
        storage_dir = os.path.join(os.path.dirname(__file__), "Computation_Deposit")
        embedding_filepath = os.path.join(internal_dir, 'umap_embedding.parquet')
        plot_filepath = os.path.join(internal_dir, 'unclustered_umap_plot.html')
        dist_matrix_filepath = os.path.join(internal_dir, 'h3_chain_dist_matrix.npz')

        embedding_df = None

        for filename in os.listdir(storage_dir):
            path_full_scanning = os.path.join(storage_dir, filename)
            try:
                file_key = Path(filename).stem
                new_data[file_key] = pd.read_csv(path_full_scanning)
            except Exception as e:
                self.logger.log("Computation_Deposit is polluted, how on earth did you manage that? Delete trash and RERUN")
                raise e

        if os.path.exists(embedding_filepath):
            self.logger.log(f"GenerateUMAP → Found cached embedding at {embedding_filepath}. Loading from file.")
            embedding_df = pd.read_parquet(embedding_filepath)
        else:
            if not os.path.exists(dist_matrix_filepath):
                self.logger.log(f"GenerateUMAP → Problem: Distance matrix not found at expected: {dist_matrix_filepath}")
                raise FileNotFoundError(f"ERROR: Ensure to compute distances, nothing found at expected: {dist_matrix_filepath}")

            self.logger.log(f"GenerateUMAP → Loading distance matrix from {dist_matrix_filepath}")
            loaded_data = np.load(dist_matrix_filepath, allow_pickle=True)
            distance_matrix, sequences = loaded_data['matrix'], loaded_data['sequences']

            self.logger.log("GenerateUMAP → Generating UMAP embedding...")
            coords_2d = UMAP(distance_matrix)
            
            embedding_df = pd.DataFrame({
                'h3_chain': sequences,
                'umap_x': coords_2d[:, 0],
                'umap_y': coords_2d[:, 1]
            })

            embedding_df.to_parquet(embedding_filepath, index=False)
            self.logger.log(f"GenerateUMAP → Saved new embedding to {embedding_filepath}")
    
        merged_df = merger_func(embedding_df, new_data["cdr"])
        new_data['umap_embedding'] = merged_df
        self.logger.log(f"GenerateUMAP → Complete. Added 'umap_embedding' DataFrame with {len(embedding_df)} rows.")
        
        self.logger.log("GenerateUMAP → Generating interactive unclustered plot...")
        fig = px.scatter(
            merged_df, 
            x='umap_x',
            y='umap_y',
            hover_name='h3_chain',
            hover_data=['pdb_id', 'heavy_host_organism_name'], 
            title='Unclustered UMAP Projection of Sequences')
        fig.update_traces(marker=dict(size=5, color='blue', opacity=0.7))
        fig.update_layout(xaxis_title="UMAP Dimension 1", yaxis_title="UMAP Dimension 2")
        fig.write_html(plot_filepath)
        self.logger.log(f"GenerateUMAP → Saved interactive plot to {plot_filepath}")
        
        return new_data

#HDBSCAN still needs the user to declare the minimum size for a cluster, but you cannot be sure of what the correct number should be if you dont know your data well; hyperpatametre tuning therefore is done here...
class GetDensity(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data = data.copy()
        if 'umap_embedding' not in new_data:
            self.logger.log("GetDensity → Skipping: 'umap_embedding' data not found.")
            return new_data

        self.logger.log("GetDensity → Analysing HDBSCAN cluster stability...")
        embedding_df = new_data['umap_embedding']
        coords_2d = embedding_df[['umap_x', 'umap_y']].to_numpy()
        
        internal_dir = os.path.join(os.path.dirname(__file__), "Internal_Files")
        stability_plot_path = os.path.join(internal_dir, 'stability_plot.html')

        stability_df = DBCV(coords_2d)
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=stability_df['min_cluster_size'],
            y=stability_df['dbcv_score'],
            mode='lines+markers',
            name='DBCV Score'
        ))
        fig.update_layout(
            title='HDBSCAN Cluster Stability (Higher is Better)',
            xaxis_title='Minimum Cluster Size (i.e. Stability)',
            yaxis_title='Relative Validity (DBCV Score)'
        )
        fig.write_html(stability_plot_path)
        self.logger.log(f"GetDensity → Saved interactive stability plot to {stability_plot_path}")

        best_size = stability_df.loc[stability_df['dbcv_score'].idxmax()]
        self.logger.log(f"GetDensity → Analysis complete. Best score found:\n{best_size}")
        
        new_data['cluster_stability'] = stability_df
        return new_data

#Observe cluster aveues in form of dendrogram.
class Dendrogram(Step):
    def __init__(self, logger: FileLogger, min_cluster_size: int):
        super().__init__(logger)
        self.min_cluster_size = min_cluster_size

    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data = data.copy()
        if 'umap_embedding' not in new_data:
            self.logger.log("Dendrogram → Skipping: 'umap_embedding' data not found.")
            return new_data

        self.logger.log(f"Dendrogram → Generating dendrogram with min_cluster_size={self.min_cluster_size}...")
        embedding_df = new_data['umap_embedding']
        coords_2d = embedding_df[['umap_x', 'umap_y']].to_numpy()
        
        internal_dir = os.path.join(os.path.dirname(__file__), "Internal_Files")
        tree_plot_path = os.path.join(internal_dir, 'dendrogram_plot.html')
        
        fig = ff.create_dendrogram(coords_2d, labels=embedding_df['sequence'].to_list())
        fig.update_layout(title='Sequence Cluster Dendrogram')
        fig.write_html(tree_plot_path)
        self.logger.log(f"Dendrogram → Saved interactive dendrogram to {tree_plot_path}")
        
        return new_data

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

        internal_dir = os.path.join(os.path.dirname(__file__), "Internal_Files")
        interactive_plot_path = os.path.join(internal_dir, 'cluster_plot.html')
        
        fig = px.scatter(
            embedding_df,
            x='umap_x',
            y='umap_y',
            color=embedding_df['cluster_id'].astype(str),
            hover_name='sequence',
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