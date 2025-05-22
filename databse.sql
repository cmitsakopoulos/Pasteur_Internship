BEGIN;
DROP TABLE IF EXISTS cdr3_primary;
DROP TABLE IF EXISTS cdr3_central;
DROP TABLE IF EXISTS antigen_central;
DROP TABLE IF EXISTS antigen_primary;
DROP TABLE IF EXISTS antigen_antibody_binding;

CREATE TABLE antigen_central (
    antigen_id          INT      PRIMARY KEY, 
    name                TEXT,
    species             TEXT,
    method              TEXT,
    taxonomy_id         INT,
    resolution_angstrom DOUBLE PRECISION,
    last_update         TIMESTAMPTZ
);

CREATE TABLE antigen_primary (
    antigen_id INT PRIMARY KEY,
    sequence   TEXT,
    isoelectric             DOUBLE PRECISION,
    gravy            DOUBLE PRECISION,
    net_charge_normal       DOUBLE PRECISION,
    net_charge_inflammation DOUBLE PRECISION,
    CONSTRAINT fk_antigen_primary_to_central
      FOREIGN KEY (antigen_id)
      REFERENCES antigen_central(antigen_id)
);

CREATE INDEX index_antigen_with_cdr3
  ON antigen_primary (antigen_id); --Indexing should help with lookup speed especially when automating look up commands for training an AI model....but might need to remove if overhead is too much; if too much ram usage/space usage

CREATE TABLE cdr3_central (
    cdr3_id             INT    PRIMARY KEY,
    method              TEXT,
    resolution_angstrom DOUBLE PRECISION,
    species             TEXT,
    taxonomy_id         INT,
    last_update         TIMESTAMPTZ);

CREATE TABLE cdr3_primary (
    cdr3_id     INT PRIMARY KEY,
    sequence    TEXT,
    gravy              DOUBLE PRECISION,
    isoelectric        DOUBLE PRECISION,
    net_charge_normal  DOUBLE PRECISION,
    net_charge_inflamed DOUBLE PRECISION,
    CONSTRAINT fk_cdr3_primary_to_central
      FOREIGN KEY (cdr3_id)
      REFERENCES cdr3_central(cdr3_id)
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

