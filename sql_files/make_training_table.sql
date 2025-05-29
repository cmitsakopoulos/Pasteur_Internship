-- The idea behind merging the training table content in SQL, lies in reducing performance overhead that is required when using SQLAlchemy with pandas and having to undergo additional processing within Python, before we even get to sourcing the information we need.
-- IMPORTANTLY: I can export the database image for you and give it to you, as I used it, minimising reproducibility complexity; obviously...https://snapshooter.com/learn/import-export-postgresql-database
BEGIN;
DROP TABLE training_dataset IF EXISTS;
CREATE TABLE training_dataset AS
SELECT
  relationship.cdr3_id,
  relationship.antigen_id,
  cdrp.h3_chain,
  cdrp.l3_chain,
  cdrp.h3_gravy,
  cdrp.l3_gravy,
  cdrp.h3_pI,
  cdrp.l3_pI,
  cdrp.h3_N_gylcosylation_sites,
  cdrp.h3_O_gylcosylation_sites,
  cdrp.l3_N_gylcosylation_sites,
  cdrp.l3_O_gylcosylation_sites,
  cdrp.l3_net_charge_inflamed,
  cdrp.l3_net_charge_normal,
  cdrp.h3_net_charge_inflamed,
  cdrp.h3_net_charge_normal,
  apr.sequence                     AS antigen_sequence,
  apr.isoelectric                  AS antigen_isoelectric,
  apr.gravy                        AS antigen_gravy,
  apr.net_charge_normal            AS antigen_net_charge_normal,
  apr.net_charge_inflammation      AS antigen_net_charge_inflamed
FROM antigen_antibody_binding AS relationship
JOIN cdr3_primary    AS cdrp
  ON relationship.cdr3_id    = cdrp.cdr3_id
JOIN antigen_primary AS apr
  ON relationship.antigen_id = apr.antigen_id
;
COMMIT;