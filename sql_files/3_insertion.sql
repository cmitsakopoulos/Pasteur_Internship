BEGIN;
------------------------------------------------------------
-- 1.  ANTIGEN  (parent)
------------------------------------------------------------
INSERT INTO antigen_central (
    antigen_id,
    organism_name,
    host_organism,
    method,
    taxonomy_id,
    database_origin,
    corresponding_pdb_antibody,
    is_incomplete,
    resolution_angstrom,
    last_update
)
SELECT
    antigen_computed_id,
    antigen_organism_name,
    antigen_host_organism,
    method,
    antigen_taxonomy_id,
    database_origin,
    corresponding_pdb_antibody,
    COALESCE(antigen_is_incomplete, FALSE)            AS is_incomplete,
    resolution                                         AS resolution_angstrom,
    last_update
FROM staging_antigen
WHERE antigen_computed_id IS NOT NULL
ON CONFLICT (antigen_id) DO NOTHING;



------------------------------------------------------------
-- 2.  ANTIGEN PRIMARY (child 1-to-1 with antigen_central)
------------------------------------------------------------
INSERT INTO antigen_primary (
    antigen_id,
    sequence,
    isoelectric,
    geary_hydrophobicity,
    blood_geary_charge,
    inflamed_geary_charge
)
SELECT
    antigen_computed_id,
    antigen_seq,
    antigen_pi,
    antigen_geary_hydrophobicity,
    antigen_blood_geary_charge,
    antigen_inflamed_geary_charge
FROM staging_antigen
WHERE antigen_computed_id IS NOT NULL
ON CONFLICT (antigen_id) DO NOTHING;



------------------------------------------------------------
-- 3.  CDR3 (parent)
------------------------------------------------------------
INSERT INTO cdr3_central (
    cdr3_id,
    method,
    database_origin,
    heavy_taxonomy_id,
    heavy_host_organism_name,
    pdb_id,
    resolution_angstrom,
    species,
    taxonomy_id,
    last_update,
    h3_is_incomplete,
    l3_is_incomplete
)
SELECT
    cdr_computed_id,
    method,
    database_origin,
    heavy_taxonomy_id,
    heavy_host_organism_name,
    pdb_id,
    resolution                                  AS resolution_angstrom,
    NULL::text                                  AS species,        -- <-- add if available
    NULL::text                             AS taxonomy_id,    -- <-- add if available
    last_update,
    h3_is_incomplete,
    l3_is_incomplete
FROM staging_cdr
WHERE cdr_computed_id IS NOT NULL
ON CONFLICT (cdr3_id) DO NOTHING;



------------------------------------------------------------
-- 4.  CDR3 PRIMARY  (child 1-to-1 with cdr3_central)
------------------------------------------------------------
INSERT INTO cdr3_primary (
    cdr3_id,
    h3_chain,
    l3_chain,
    h3_geary_hydrophobicity,
    l3_geary_hydrophobicity,
    h3_pi,
    l3_pi,
    h3_n_glycosilation_sites,
    h3_o_glycosilation_sites,
    l3_n_glycosilation_sites,
    l3_o_glycosilation_sites,
    l3_inflamed_geary_charge,
    l3_blood_geary_charge,
    h3_inflamed_geary_charge,
    h3_blood_geary_charge
)
SELECT
    cdr_computed_id,
    h3_chain,
    l3_chain,
    h3_geary_hydrophobicity,
    l3_geary_hydrophobicity,
    h3_pi,
    l3_pi,
    h3_n_glycosilation_sites,
    h3_o_glycosilation_sites,
    l3_n_glycosilation_sites,
    l3_o_glycosilation_sites,
    h3_blood_geary_charge,       
    h3_inflamed_geary_charge,
    l3_inflamed_geary_charge,
    l3_blood_geary_charge
FROM staging_cdr
WHERE cdr_computed_id IS NOT NULL
ON CONFLICT (cdr3_id) DO NOTHING;



------------------------------------------------------------
-- 5.  ANTIGEN ⟷ CDR3  (junction)
------------------------------------------------------------
INSERT INTO antigen_antibody_binding (
    cdr3_id,
    antigen_id
)
SELECT
    cdr_computed_id,
    antigen_computed_id
FROM staging_relationships
WHERE cdr_computed_id    IS NOT NULL
  AND antigen_computed_id IS NOT NULL
ON CONFLICT DO NOTHING;

COMMIT;