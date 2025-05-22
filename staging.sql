-- MAKE STAGING TABLE FOR ANTIGENS
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
  antigen_gravy                 DOUBLE PRECISION,
  antigen_pI                    DOUBLE PRECISION,
  antigen_net_charge_inflamed   DOUBLE PRECISION,
  antigen_net_charge_normal     DOUBLE PRECISION,
  antigen_is_incomplete         FLOAT,
  antigen_computed_id           INT
);

-- COMMIT THE CSV INTO SQL
\copy staging_antigen
  FROM '/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Computation_Deposit/antigen_20250522_125947.csv'
  WITH (
    FORMAT csv,
    HEADER,
    NULL ''            
  );

--Because of how I set up the int as a boolean flag within the pipeline, we will have to stick with this method
ALTER TABLE staging_antigen
  ALTER COLUMN antigen_is_incomplete TYPE BOOLEAN
  USING (antigen_is_incomplete = 1);

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
  h3_gravy                      DOUBLE PRECISION,
  h3_pI                         DOUBLE PRECISION,
  h3_net_charge_inflamed        DOUBLE PRECISION,
  h3_net_charge_normal          DOUBLE PRECISION,
  l3_gravy                      DOUBLE PRECISION,
  l3_pI                         DOUBLE PRECISION,
  l3_net_charge_inflamed        DOUBLE PRECISION,
  l3_net_charge_normal          DOUBLE PRECISION,
  h3_is_incomplete              FLOAT, 
  h3_N_gylcosylation_sites      FLOAT,
  h3_O_gylcosylation_sites      FLOAT,
  l3_is_incomplete              FLOAT,
  l3_N_gylcosylation_sites      FLOAT,
  l3_O_gylcosylation_sites      FLOAT,
  cdr_computed_id               INT
);

\copy staging_cdr
  FROM '/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Computation_Deposit/cdr_20250522_125947.csv'
  WITH (
    FORMAT csv,
    HEADER,
    NULL ''            
  );

ALTER TABLE staging_cdr ALTER COLUMN cdr_is_incomplete TYPE BOOLEAN USING (cdr_is_incomplete = 1);

CREATE TABLE relationships_staging (
  antigen_computed_id  INT,
  cdr_computed_id      INT
);

\copy relationships_staging
  FROM '/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Computation_Deposit/relationships_20250522_125947.csv'
  WITH (
    FORMAT csv,
    HEADER,
    NULL ''            
  );

  