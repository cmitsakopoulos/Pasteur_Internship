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
  h3_gravy                      DOUBLE PRECISION,
  h3_pi                        DOUBLE PRECISION,
  h3_net_charge_inflamed        DOUBLE PRECISION,
  h3_net_charge_normal          DOUBLE PRECISION,
  l3_gravy                      DOUBLE PRECISION,
  l3_pi                         DOUBLE PRECISION,
  l3_net_charge_inflamed        DOUBLE PRECISION,
  l3_net_charge_normal          DOUBLE PRECISION,
  h3_is_incomplete              FLOAT, 
  h3_n_gylcosylation_sites      FLOAT,
  h3_o_gylcosylation_sites      FLOAT,
  l3_is_incomplete              FLOAT,
  l3_n_gylcosylation_sites      FLOAT,
  l3_o_gylcosylation_sites      FLOAT,
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
  antigen_gravy                 DOUBLE PRECISION,
  antigen_pi      DOUBLE PRECISION,
  antigen_net_charge_inflamed   DOUBLE PRECISION,
  antigen_net_charge_normal     DOUBLE PRECISION,
  antigen_is_incomplete         FLOAT,
  antigen_computed_id           INT
);

CREATE TABLE staging_relationships (
  antigen_computed_id  INT,
  cdr_computed_id      INT
);
COMMIT;