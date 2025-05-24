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
    gravy,
    net_charge_normal,
    net_charge_inflammation
)
SELECT
    antigen_computed_id,
    antigen_seq,
    antigen_pI,
    antigen_gravy,
    antigen_net_charge_normal,
    antigen_net_charge_inflamed
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
    h3_gravy,
    l3_gravy,
    h3_pI,
    l3_pI,
    h3_N_gylcosylation_sites,
    h3_O_gylcosylation_sites,
    l3_N_gylcosylation_sites,
    l3_O_gylcosylation_sites,
    l3_net_charge_inflamed,
    l3_net_charge_normal,
    h3_net_charge_inflamed,
    h3_net_charge_normal
)
SELECT
    cdr_computed_id,
    h3_chain,
    l3_chain,
    h3_gravy,
    l3_gravy,
    h3_pI,
    l3_pI,
    h3_N_gylcosylation_sites,
    h3_O_gylcosylation_sites,
    l3_N_gylcosylation_sites,
    l3_O_gylcosylation_sites,
    h3_net_charge_normal,       
    h3_net_charge_inflamed,
    l3_net_charge_inflamed,
    l3_net_charge_normal
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
FROM relationships_staging
WHERE cdr_computed_id    IS NOT NULL
  AND antigen_computed_id IS NOT NULL
ON CONFLICT DO NOTHING;

COMMIT;