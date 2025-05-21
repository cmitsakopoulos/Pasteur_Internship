DROP TABLE IF EXISTS cdr3_chemical;
DROP TABLE IF EXISTS cdr3_primary;
DROP TABLE IF EXISTS cdr3_central;
DROP TABLE IF EXISTS antigen_central;

CREATE TABLE antigen_central (
    antigen_id          VARCHAR(12)      PRIMARY KEY, 
    name                TEXT,
    species             TEXT,
    method              TEXT,
    taxonomy_id         INT,
    resolution_angstrom NUMERIC(5,2),
    last_update         TIMESTAMPTZ
);

CREATE TABLE antigen_primary (
    antigen_id VARCHAR(12) PRIMARY KEY,
    sequence   TEXT,
    CONSTRAINT fk_antigen_primary_to_central
      FOREIGN KEY (antigen_id)
      REFERENCES antigen_central(antigen_id)
);
CREATE INDEX idx_antigen_primary_antigen_id
  ON antigen_primary (antigen_id);

CREATE TABLE antigen_chemical (
    antigen_id              VARCHAR(12) PRIMARY KEY,
    isoelectric             NUMERIC(6,4),
    gravy_score             NUMERIC(6,4),
    net_charge_normal       NUMERIC(6,4),
    net_charge_inflammation NUMERIC(6,4),
    CONSTRAINT fk_antigen_chemical_to_central
      FOREIGN KEY (antigen_id)
      REFERENCES antigen_central(antigen_id)
);
CREATE INDEX idx_antigen_chemical_antigen_id
  ON antigen_chemical (antigen_id);

CREATE TABLE cdr3_central (
    cdr3_id             VARCHAR(4)    PRIMARY KEY,
    antigen_id          VARCHAR(12)   NOT NULL
      REFERENCES antigen_central(antigen_id),
    method              TEXT,
    resolution_angstrom NUMERIC(5,2),
    species             TEXT,
    taxonomy_id         INT,
    last_update         TIMESTAMPTZ,
    CONSTRAINT fk_cdr3_central_to_antigen
      FOREIGN KEY (antigen_id)
      REFERENCES antigen_central(antigen_id)
);
CREATE INDEX idx_cdr3_central_antigen_id
  ON cdr3_central (antigen_id);

CREATE TABLE cdr3_primary (
    cdr3_id     VARCHAR(4) PRIMARY KEY,
    sequence    TEXT,
    CONSTRAINT fk_cdr3_primary_to_central
      FOREIGN KEY (cdr3_id)
      REFERENCES cdr3_central(cdr3_id)
);

CREATE TABLE cdr3_chemical (
    cdr3_id            VARCHAR(4)   PRIMARY KEY,
    gravy              NUMERIC(6,4),
    isoelectric        NUMERIC(6,4),
    net_charge_normal  NUMERIC(6,4),
    net_charge_inflamed NUMERIC(6,4),
    CONSTRAINT fk_cdr3_chemical_to_central
      FOREIGN KEY (cdr3_id)
      REFERENCES cdr3_central(cdr3_id)
);