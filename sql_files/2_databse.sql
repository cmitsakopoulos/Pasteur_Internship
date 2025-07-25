BEGIN;
DROP TABLE IF EXISTS antigen_antibody_binding;
DROP TABLE IF EXISTS antigen_primary;
DROP TABLE IF EXISTS cdr3_primary;
DROP TABLE IF EXISTS antigen_central;
DROP TABLE IF EXISTS cdr3_central;

CREATE TABLE antigen_central (
    antigen_id          INT, 
    organism_name         TEXT,
    host_organism         TEXT,
    method              TEXT,
    taxonomy_id         TEXT,
    database_origin TEXT,
    corresponding_pdb_antibody    TEXT, 
    is_incomplete BOOLEAN,
    resolution_angstrom DOUBLE PRECISION,
    last_update         TIMESTAMPTZ,
    CONSTRAINT pk_antigen_id_main PRIMARY KEY (antigen_id)
);

CREATE TABLE antigen_primary (
    antigen_id INT PRIMARY KEY,
    sequence   TEXT,
    isoelectric             DOUBLE PRECISION,
    geary_hydrophobicity            TEXT,
    blood_geary_charge       TEXT,
    inflamed_geary_charge TEXT,
    CONSTRAINT fk_antigen_primary_to_central
      FOREIGN KEY (antigen_id)
      REFERENCES antigen_central(antigen_id) ON DELETE CASCADE
);

CREATE INDEX index_antigen_with_cdr3
  ON antigen_primary (antigen_id); 

CREATE TABLE cdr3_central (
    cdr3_id             INT,
    method              TEXT,
    database_origin TEXT,
    heavy_taxonomy_id         TEXT,
    heavy_host_organism_name      TEXT,
    pdb_id                        TEXT,
    resolution_angstrom DOUBLE PRECISION,
    h3_is_incomplete BOOLEAN,
    species             TEXT,
    taxonomy_id         TEXT,
    l3_is_incomplete BOOLEAN,
    last_update         TIMESTAMPTZ,
    CONSTRAINT pk_cdr3_id_main PRIMARY KEY (cdr3_id));

CREATE TABLE cdr3_primary (
    cdr3_id     INT PRIMARY KEY,
    h3_chain TEXT,
    l3_chain    TEXT,
    h3_geary_hydrophobicity              TEXT,
    l3_geary_hydrophobicity                      TEXT,
    h3_pi        DOUBLE PRECISION,
    l3_pi        DOUBLE PRECISION,
    h3_n_glycosilation_sites      FLOAT,
    h3_o_glycosilation_sites      FLOAT,
    l3_n_glycosilation_sites      FLOAT,
    l3_o_glycosilation_sites      FLOAT,
    l3_inflamed_geary_charge  TEXT,
    l3_blood_geary_charge TEXT,
    h3_inflamed_geary_charge      TEXT,
    h3_blood_geary_charge      TEXT,
    CONSTRAINT fk_cdr3_primary_to_central
      FOREIGN KEY (cdr3_id)
      REFERENCES cdr3_central(cdr3_id) ON DELETE CASCADE
);

CREATE TABLE antigen_antibody_binding (
  cdr3_id INT,
  antigen_id INT,
  PRIMARY KEY (cdr3_id, antigen_id),
  FOREIGN KEY (cdr3_id) REFERENCES cdr3_central(cdr3_id),
  FOREIGN KEY (antigen_id) REFERENCES antigen_central(antigen_id)
);

CREATE INDEX idx_binding_by_antigen 
  ON antigen_antibody_binding (antigen_id);

COMMIT;