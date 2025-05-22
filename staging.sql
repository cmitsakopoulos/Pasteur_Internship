-- MAKE STAGING TABLE FOR ANTIGENS
CREATE TABLE staging_antigen (
  antigen_seq                   TEXT,
  database_origin               TEXT,
  antigen_organism_name         TEXT,
  antigen_host_organism         TEXT,
  antigen_taxonomy_id           INT,
  resolution                    DOUBLE PRECISION,
  method                        TEXT,
  last_update                   TIMESTAMPTZ,
  corresponding_pdb_antibody    TEXT,      
  antigen_gravy                 DOUBLE PRECISION,
  antigen_pI                    DOUBLE PRECISION,
  antigen_net_charge_inflamed   DOUBLE PRECISION,
  antigen_net_charge_normal     DOUBLE PRECISION,
  antigen_is_incomplete         SMALLINT,
  antigen_computed_id           INT
);

--Because of how I set up the int as a boolean flag within the pipeline, we will have to stick with this method
ALTER TABLE staging_antigen
  ALTER COLUMN antigen_is_incomplete TYPE BOOLEAN
  USING (antigen_is_incomplete = 1);