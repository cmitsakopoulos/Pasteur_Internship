BEGIN;
DROP TABLE IF EXISTS staging_cdr;
DROP TABLE IF EXISTS staging_antigen;
DROP TABLE IF EXISTS staging_relationships;
CREATE TABLE staging_cdr (
  h3_chain                      TEXT,
  database_origin               TEXT,
  heavy_taxonomy_id             TEXT,
  heavy_host_organism_name      TEXT,
  resolution                    DOUBLE PRECISION,
  method                        TEXT,
  last_update                   TIMESTAMPTZ,
  l3_chain                      TEXT,
  pdb_id                        TEXT,      
  h3_geary_hydrophobicity                      DOUBLE PRECISION[],
  h3_pi                        DOUBLE PRECISION,
  h3_inflamed_geary_charge        DOUBLE PRECISION[],
  h3_blood_geary_charge          DOUBLE PRECISION[],
  l3_geary_hydrophobicity                      DOUBLE PRECISION[],
  l3_pi                         DOUBLE PRECISION,
  l3_inflamed_geary_charge        DOUBLE PRECISION[],
  l3_blood_geary_charge          DOUBLE PRECISION[],
  h3_is_incomplete              FLOAT, 
  h3_n_glycosilation_sites      FLOAT,
  h3_o_glycosilation_sites      FLOAT,
  l3_is_incomplete              FLOAT,
  l3_n_glycosilation_sites      FLOAT,
  l3_o_glycosilation_sites      FLOAT,
  cdr_computed_id               INT
);
CREATE TABLE staging_antigen (
  antigen_seq                   TEXT,
  database_origin               TEXT,
  antigen_organism_name         TEXT,
  antigen_host_organism         TEXT,
  antigen_taxonomy_id           TEXT,
  resolution                    DOUBLE PRECISION,
  method                        TEXT,
  last_update                   TIMESTAMPTZ,
  corresponding_pdb_antibody    TEXT,      
  antigen_geary_hydrophobicity                 DOUBLE PRECISION[],
  antigen_pi      DOUBLE PRECISION,
  antigen_inflamed_geary_charge   DOUBLE PRECISION[],
  antigen_blood_geary_charge     DOUBLE PRECISION[],
  antigen_is_incomplete         FLOAT,
  antigen_computed_id           INT
);
CREATE TABLE staging_relationships (
  antigen_computed_id  INT,
  cdr_computed_id      INT
);
COMMIT;