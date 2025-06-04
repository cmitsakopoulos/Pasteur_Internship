BEGIN;
DROP TABLE IF EXISTS training_dataset;
CREATE TABLE training_dataset AS
SELECT
  relationship.cdr3_id,
  relationship.antigen_id,
  cdrp.h3_chain,
  cdrp.l3_chain,
  replace(replace(cdrp.h3_geary_hydrophobicity, '[', '{'), ']', '}')::double precision[] AS h3_geary_hydrophobicity,
  replace(replace(cdrp.l3_geary_hydrophobicity, '[', '{'), ']', '}')::double precision[] AS l3_geary_hydrophobicity,
  cdrp.h3_pi,
  cdrp.l3_pi,
  cdrp.h3_n_glycosilation_sites,
  cdrp.h3_o_glycosilation_sites,
  cdrp.l3_n_glycosilation_sites,
  cdrp.l3_o_glycosilation_sites,
  replace(replace(cdrp.l3_inflamed_geary_charge, '[', '{'), ']', '}')::double precision[] AS l3_inflamed_geary_charge,
  replace(replace(cdrp.l3_blood_geary_charge, '[', '{'), ']', '}')::double precision[] AS l3_blood_geary_charge,
  replace(replace(cdrp.h3_inflamed_geary_charge, '[', '{'), ']', '}')::double precision[] AS h3_inflamed_geary_charge,
  replace(replace(cdrp.h3_blood_geary_charge, '[', '{'), ']', '}')::double precision[] AS h3_blood_geary_charge,
  apr.sequence                     AS antigen_sequence,
  apr.isoelectric                  AS antigen_isoelectric,
  replace(replace(apr.geary_hydrophobicity, '[', '{'), ']', '}')::double precision[] AS antigen_geary_hydrophobicity,
  replace(replace(apr.blood_geary_charge, '[', '{'), ']', '}')::double precision[] AS antigen_blood_geary_charge,
  replace(replace(apr.inflamed_geary_charge, '[', '{'), ']', '}')::double precision[] AS antigen_inflamed_geary_charge
FROM antigen_antibody_binding AS relationship
JOIN cdr3_primary    AS cdrp
  ON relationship.cdr3_id    = cdrp.cdr3_id
JOIN antigen_primary AS apr
  ON relationship.antigen_id = apr.antigen_id
;
COMMIT;