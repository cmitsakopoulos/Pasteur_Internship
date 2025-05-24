-- Making the training table is preferable as select commands wont have to be run anew to make ex. a dataframe in Python, then training a model on that. Instead, the table is ready and the only computiational step is information acquisition from a pre-computed table
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
  cdrp.net_charge_normal   AS cdr_net_charge_normal,
  cdrp.net_charge_inflamed AS cdr_net_charge_inflamed,
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