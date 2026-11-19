-- Cohort extraction for: Early net fluid balance and extubation success in
-- mechanically ventilated adults with acute kidney injury (MIMIC-IV v3.1).
--
-- Platform : Google BigQuery (GoogleSQL), PhysioNet MIMIC-IV v3.1 tables.
-- Before running, replace YOUR_PROJECT.YOUR_DATASET with a dataset you can write to.
-- Run the whole script as one query job. The final SELECT returns one row per ICU stay;
-- export it as CSV and use it as input to analysis/run_analysis.py.

DECLARE WEAN_WINDOW_HOURS INT64 DEFAULT 48;

-- Selected intravenous fluids and blood products
DECLARE INTAKE_ITEMIDS ARRAY<INT64> DEFAULT [
  225158,  -- NaCl 0.9%
  225159,  -- NaCl 0.45%
  225828,  -- Lactated Ringer's
  220949,  -- Dextrose 5%
  225823,  -- D5 1/2NS
  225825,  -- D5NS
  220864,  -- Albumin 5%
  225168,  -- Packed red blood cells
  220970   -- Fresh frozen plasma
];

DECLARE VASOPRESSOR_ITEMIDS ARRAY<INT64> DEFAULT [
  221906,  -- Norepinephrine
  221289,  -- Epinephrine
  221749,  -- Phenylephrine
  221662,  -- Dopamine
  222315,  -- Vasopressin
  221653   -- Dobutamine
];

CREATE TEMP TABLE cohort_all AS

WITH icu AS (
  SELECT
    subject_id, hadm_id, stay_id,
    CAST(intime AS TIMESTAMP) AS intime,
    CAST(outtime AS TIMESTAMP) AS outtime,
    ROW_NUMBER() OVER (PARTITION BY subject_id ORDER BY intime, stay_id) AS icu_stay_seq
  FROM `physionet-data.mimiciv_3_1_icu.icustays`
),

adm AS (
  SELECT
    subject_id, hadm_id,
    CAST(admittime AS TIMESTAMP) AS admittime,
    CAST(dischtime AS TIMESTAMP) AS dischtime,
    CAST(deathtime AS TIMESTAMP) AS deathtime,
    hospital_expire_flag
  FROM `physionet-data.mimiciv_3_1_hosp.admissions`
),

-- KDIGO stage >= 1 within 48 h of ICU admission
aki_48h AS (
  SELECT k.stay_id, MAX(k.aki_stage) AS max_aki_stage_48h
  FROM `physionet-data.mimiciv_3_1_derived.kdigo_stages` k
  JOIN icu i ON k.stay_id = i.stay_id
  WHERE CAST(k.charttime AS TIMESTAMP) >= i.intime
    AND CAST(k.charttime AS TIMESTAMP) < TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR)
    AND k.aki_stage >= 1
  GROUP BY k.stay_id
),

-- Invasive ventilation episodes (exact status; non-invasive ventilation excluded)
vent_invasive AS (
  SELECT DISTINCT
    v.stay_id, i.subject_id, i.hadm_id,
    CAST(v.starttime AS TIMESTAMP) AS starttime,
    CAST(v.endtime AS TIMESTAMP) AS endtime
  FROM `physionet-data.mimiciv_3_1_derived.ventilation` v
  JOIN icu i ON v.stay_id = i.stay_id
  WHERE v.ventilation_status = 'InvasiveVent'
    AND v.starttime IS NOT NULL
),

-- First invasive episode per ICU stay (a later episode is never substituted)
first_vent AS (
  SELECT stay_id, starttime AS first_vent_start, endtime AS first_vent_end
  FROM vent_invasive
  QUALIFY ROW_NUMBER() OVER (PARTITION BY stay_id ORDER BY starttime, endtime) = 1
),

-- New invasive episodes within 48 h after the end of the first episode
revent AS (
  SELECT
    fv.stay_id,
    COUNT(DISTINCT IF(
      v.starttime > fv.first_vent_end
      AND v.starttime <= TIMESTAMP_ADD(fv.first_vent_end, INTERVAL WEAN_WINDOW_HOURS HOUR)
      AND v.starttime >= a.admittime
      AND (a.dischtime IS NULL OR v.starttime <= a.dischtime)
      AND (a.deathtime IS NULL OR v.starttime < a.deathtime),
      v.starttime, NULL)) AS revent_cnt_48h,                -- any ICU stay of the admission
    COUNT(DISTINCT IF(
      v.stay_id = fv.stay_id
      AND v.starttime > fv.first_vent_end
      AND v.starttime <= TIMESTAMP_ADD(fv.first_vent_end, INTERVAL WEAN_WINDOW_HOURS HOUR)
      AND v.starttime >= a.admittime
      AND (a.dischtime IS NULL OR v.starttime <= a.dischtime)
      AND (a.deathtime IS NULL OR v.starttime < a.deathtime),
      v.starttime, NULL)) AS revent_cnt_48h_same_stay,      -- index ICU stay only (primary outcome)
    COUNT(DISTINCT IF(v.stay_id != fv.stay_id AND v.starttime = fv.first_vent_end,
                      v.stay_id, NULL)) AS invasive_other_stay_at_index_end
  FROM first_vent fv
  JOIN icu i ON fv.stay_id = i.stay_id
  LEFT JOIN adm a ON i.hadm_id = a.hadm_id
  LEFT JOIN vent_invasive v ON i.subject_id = v.subject_id AND i.hadm_id = v.hadm_id
  GROUP BY fv.stay_id
),

windows AS (
  SELECT
    i.subject_id, i.hadm_id, i.stay_id, i.intime, i.outtime, i.icu_stay_seq,
    fv.first_vent_start, fv.first_vent_end,
    LEAST(TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR), fv.first_vent_end) AS pre_extub_or_48h_end
  FROM icu i
  JOIN first_vent fv ON i.stay_id = fv.stay_id
),

-- Selected IV intake, 0-48 h (assigned by event start time)
intake_0_48 AS (
  SELECT
    ie.stay_id,
    SUM(CASE
          WHEN LOWER(TRIM(COALESCE(ie.amountuom, ''))) = 'ml' THEN ie.amount
          WHEN LOWER(TRIM(COALESCE(ie.amountuom, ''))) = 'l' THEN ie.amount * 1000
        END) AS intake_ml_0_48,
    COUNTIF(ie.amount IS NOT NULL AND LOWER(TRIM(COALESCE(ie.amountuom, ''))) NOT IN ('ml', 'l')) AS intake_unknown_unit_events,
    COUNTIF(CAST(ie.endtime AS TIMESTAMP) > TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR)) AS intake_crosses_48h_events
  FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
  JOIN icu i ON ie.stay_id = i.stay_id
  WHERE CAST(ie.starttime AS TIMESTAMP) >= i.intime
    AND CAST(ie.starttime AS TIMESTAMP) < TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR)
    AND ie.itemid IN UNNEST(INTAKE_ITEMIDS)
    AND LOWER(COALESCE(ie.statusdescription, '')) != 'rewritten'
  GROUP BY ie.stay_id
),

-- Selected IV intake, ICU admission to the earlier of extubation or 48 h
intake_pre_extub_or_48 AS (
  SELECT
    ie.stay_id,
    SUM(CASE
          WHEN LOWER(TRIM(COALESCE(ie.amountuom, ''))) = 'ml' THEN ie.amount
          WHEN LOWER(TRIM(COALESCE(ie.amountuom, ''))) = 'l' THEN ie.amount * 1000
        END) AS intake_ml_pre_extub_or_48,
    COUNTIF(ie.amount IS NOT NULL AND LOWER(TRIM(COALESCE(ie.amountuom, ''))) NOT IN ('ml', 'l')) AS intake_pre_unknown_unit_events,
    COUNTIF(CAST(ie.endtime AS TIMESTAMP) > w.pre_extub_or_48h_end) AS intake_pre_crosses_end_events
  FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
  JOIN windows w ON ie.stay_id = w.stay_id
  WHERE CAST(ie.starttime AS TIMESTAMP) >= w.intime
    AND CAST(ie.starttime AS TIMESTAMP) < w.pre_extub_or_48h_end
    AND ie.itemid IN UNNEST(INTAKE_ITEMIDS)
    AND LOWER(COALESCE(ie.statusdescription, '')) != 'rewritten'
  GROUP BY ie.stay_id
),

-- All recorded inputevents volumes, 0-48 h (comprehensive-balance sensitivity analysis)
all_input_0_48 AS (
  SELECT
    ie.stay_id,
    SUM(CASE
          WHEN LOWER(COALESCE(ie.amountuom, '')) = 'ml' THEN COALESCE(ie.amount, 0)
          WHEN LOWER(COALESCE(ie.amountuom, '')) = 'l' THEN COALESCE(ie.amount, 0) * 1000
          ELSE 0
        END) AS all_input_ml_0_48
  FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
  JOIN icu i ON ie.stay_id = i.stay_id
  WHERE CAST(ie.starttime AS TIMESTAMP) >= i.intime
    AND CAST(ie.starttime AS TIMESTAMP) < TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR)
    AND LOWER(COALESCE(ie.statusdescription, '')) != 'rewritten'
  GROUP BY ie.stay_id
),

urine_0_48 AS (
  SELECT u.stay_id, SUM(COALESCE(u.urineoutput, 0)) AS urine_ml_0_48
  FROM `physionet-data.mimiciv_3_1_derived.urine_output` u
  JOIN icu i ON u.stay_id = i.stay_id
  WHERE CAST(u.charttime AS TIMESTAMP) >= i.intime
    AND CAST(u.charttime AS TIMESTAMP) < TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR)
  GROUP BY u.stay_id
),

urine_pre_extub_or_48 AS (
  SELECT u.stay_id, SUM(COALESCE(u.urineoutput, 0)) AS urine_ml_pre_extub_or_48
  FROM `physionet-data.mimiciv_3_1_derived.urine_output` u
  JOIN windows w ON u.stay_id = w.stay_id
  WHERE CAST(u.charttime AS TIMESTAMP) >= w.intime
    AND CAST(u.charttime AS TIMESTAMP) < w.pre_extub_or_48h_end
  GROUP BY u.stay_id
),

-- RRT-related volumes in outputevents (none are captured; ultrafiltration is charted in chartevents)
rrt_output_0_48 AS (
  SELECT
    oe.stay_id,
    SUM(oe.value) AS rrt_output_ml_0_48_candidate,
    COUNT(*) AS rrt_output_candidate_records
  FROM `physionet-data.mimiciv_3_1_icu.outputevents` oe
  JOIN `physionet-data.mimiciv_3_1_icu.d_items` di ON oe.itemid = di.itemid
  JOIN icu i ON oe.stay_id = i.stay_id
  WHERE CAST(oe.charttime AS TIMESTAMP) >= i.intime
    AND CAST(oe.charttime AS TIMESTAMP) < TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR)
    AND LOWER(COALESCE(oe.valueuom, '')) = 'ml'
    AND (REGEXP_CONTAINS(LOWER(COALESCE(di.label, '')), r'ultra.?filtration|dialysis|dialysate|crrt')
         OR REGEXP_CONTAINS(LOWER(COALESCE(di.category, '')), r'dialysis|crrt'))
    AND NOT REGEXP_CONTAINS(LOWER(COALESCE(di.label, '')), r'urine|foley|void|irrigat')
  GROUP BY oe.stay_id
),

-- Non-urine, non-RRT outputevents volumes, 0-48 h
nonurine_nonrrt_output_0_48 AS (
  SELECT oe.stay_id, SUM(COALESCE(oe.value, 0)) AS nonurine_nonrrt_output_ml_0_48
  FROM `physionet-data.mimiciv_3_1_icu.outputevents` oe
  JOIN `physionet-data.mimiciv_3_1_icu.d_items` di ON oe.itemid = di.itemid
  JOIN icu i ON oe.stay_id = i.stay_id
  WHERE CAST(oe.charttime AS TIMESTAMP) >= i.intime
    AND CAST(oe.charttime AS TIMESTAMP) < TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR)
    AND LOWER(COALESCE(oe.valueuom, '')) = 'ml'
    AND NOT REGEXP_CONTAINS(LOWER(COALESCE(di.label, '')), r'urine|foley|void|irrigat')
    AND NOT (REGEXP_CONTAINS(LOWER(COALESCE(di.label, '')), r'ultra.?filtration|dialysis|dialysate|crrt')
             OR REGEXP_CONTAINS(LOWER(COALESCE(di.category, '')), r'dialysis|crrt'))
  GROUP BY oe.stay_id
),

rrt_48h AS (
  SELECT r.stay_id, MAX(r.dialysis_present) AS dialysis_present_48h
  FROM `physionet-data.mimiciv_3_1_derived.rrt` r
  JOIN icu i ON r.stay_id = i.stay_id
  WHERE CAST(r.charttime AS TIMESTAMP) >= i.intime
    AND CAST(r.charttime AS TIMESTAMP) < TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR)
  GROUP BY r.stay_id
),

detail AS (
  SELECT stay_id, admission_age, gender
  FROM `physionet-data.mimiciv_3_1_derived.icustay_detail`
),

sofa AS (
  SELECT stay_id, sofa AS first_day_sofa
  FROM `physionet-data.mimiciv_3_1_derived.first_day_sofa`
),

wt AS (
  SELECT stay_id, weight
  FROM `physionet-data.mimiciv_3_1_derived.first_day_weight`
),

comorb AS (
  SELECT
    hadm_id,
    MAX(CASE
          WHEN icd_version = 9 AND (REGEXP_CONTAINS(icd_code, r'^428')
                                    OR REGEXP_CONTAINS(icd_code, r'^402[019]1')
                                    OR REGEXP_CONTAINS(icd_code, r'^404[019][13]')) THEN 1
          WHEN icd_version = 10 AND (REGEXP_CONTAINS(icd_code, r'^I50')
                                     OR icd_code IN ('I110', 'I130', 'I132')) THEN 1
          ELSE 0 END) AS chf_icd,
    MAX(CASE
          WHEN icd_version = 9 AND REGEXP_CONTAINS(icd_code, r'^41[0-4]') THEN 1
          WHEN icd_version = 10 AND REGEXP_CONTAINS(icd_code, r'^I2[0-5]') THEN 1
          ELSE 0 END) AS cad_icd,
    MAX(CASE
          WHEN icd_version = 9 AND REGEXP_CONTAINS(icd_code, r'^585') THEN 1
          WHEN icd_version = 10 AND REGEXP_CONTAINS(icd_code, r'^N18') THEN 1
          ELSE 0 END) AS ckd_icd,
    MAX(CASE
          WHEN icd_version = 9 AND (REGEXP_CONTAINS(icd_code, r'^038')
                                    OR icd_code IN ('99591', '99592', '78552')) THEN 1
          WHEN icd_version = 10 AND (REGEXP_CONTAINS(icd_code, r'^A4[01]')
                                     OR REGEXP_CONTAINS(icd_code, r'^R652')
                                     OR icd_code = 'R572') THEN 1
          ELSE 0 END) AS sepsis_icd
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
  GROUP BY hadm_id
),

sepsis3 AS (
  SELECT stay_id, MAX(CASE WHEN sepsis3 THEN 1 ELSE 0 END) AS sepsis3_derived
  FROM `physionet-data.mimiciv_3_1_derived.sepsis3`
  GROUP BY stay_id
),

vasopressor_48h AS (
  SELECT ie.stay_id, 1 AS vasopressor_use_48h
  FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
  JOIN icu i ON ie.stay_id = i.stay_id
  WHERE CAST(ie.starttime AS TIMESTAMP) >= i.intime
    AND CAST(ie.starttime AS TIMESTAMP) < TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR)
    AND ie.itemid IN UNNEST(VASOPRESSOR_ITEMIDS)
  GROUP BY ie.stay_id
),

diuretic_48h AS (
  SELECT i.stay_id, 1 AS diuretic_use_48h
  FROM icu i
  JOIN `physionet-data.mimiciv_3_1_hosp.prescriptions` p ON i.hadm_id = p.hadm_id
  WHERE CAST(p.starttime AS TIMESTAMP) < TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR)
    AND COALESCE(CAST(p.stoptime AS TIMESTAMP), CAST(p.starttime AS TIMESTAMP)) >= i.intime
    AND REGEXP_CONTAINS(LOWER(COALESCE(p.drug, '')),
          r'furosemide|lasix|bumetanide|bumex|torsemide|metolazone|chlorothiazide|hydrochlorothiazide|spironolactone')
  GROUP BY i.stay_id
),

lactate_0_48 AS (
  SELECT
    i.stay_id,
    ARRAY_AGG(l.valuenum IGNORE NULLS ORDER BY CAST(l.charttime AS TIMESTAMP) ASC LIMIT 1)[SAFE_OFFSET(0)] AS first_lactate_0_48,
    MIN(l.valuenum) AS min_lactate_0_48,
    MAX(l.valuenum) AS max_lactate_0_48
  FROM icu i
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l ON i.subject_id = l.subject_id AND i.hadm_id = l.hadm_id
  JOIN `physionet-data.mimiciv_3_1_hosp.d_labitems` dl ON l.itemid = dl.itemid
  WHERE CAST(l.charttime AS TIMESTAMP) >= i.intime
    AND CAST(l.charttime AS TIMESTAMP) < TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR)
    AND l.valuenum IS NOT NULL
    AND LOWER(dl.label) = 'lactate'
  GROUP BY i.stay_id
)

SELECT
  w.subject_id, w.hadm_id, w.stay_id,
  w.intime AS icu_intime, w.outtime AS icu_outtime,
  a.admittime AS hospital_admittime, a.dischtime AS hospital_dischtime,
  a.deathtime AS hospital_deathtime, a.hospital_expire_flag,
  w.icu_stay_seq,
  CASE WHEN w.icu_stay_seq = 1 THEN 1 ELSE 0 END AS is_first_icu_stay,
  d.admission_age, d.gender, s.first_day_sofa, wt.weight,
  aki.max_aki_stage_48h,
  w.first_vent_start, w.first_vent_end,

  -- Exclusion reasons are assigned hierarchically in the order listed
  CASE
    WHEN w.first_vent_end IS NULL THEN 'missing_episode_end'
    WHEN w.first_vent_end <= w.first_vent_start THEN 'invalid_episode_duration'
    WHEN w.first_vent_end < w.intime THEN 'episode_end_before_ICU'
    WHEN w.outtime IS NULL THEN 'missing_ICU_outtime'
    WHEN w.first_vent_start > w.outtime OR w.first_vent_end > w.outtime THEN 'episode_end_after_ICU'
    WHEN d.admission_age IS NULL OR d.admission_age < 18 THEN 'not_confirmed_adult'
    WHEN a.admittime IS NULL OR a.dischtime IS NULL THEN 'missing_hospital_times'
    WHEN w.first_vent_end < a.admittime OR w.first_vent_end > a.dischtime THEN 'episode_end_outside_hospital'
    WHEN a.hospital_expire_flag = 1 AND a.deathtime IS NULL THEN 'hospital_death_time_missing'
    WHEN a.deathtime IS NOT NULL AND a.deathtime <= w.first_vent_end THEN 'death_at_or_before_episode_end'
    ELSE 'eligible'
  END AS eligibility_status,

  TIMESTAMP_DIFF(w.first_vent_end, w.intime, HOUR) AS hours_icu_to_extubation,
  CASE WHEN w.first_vent_end <= TIMESTAMP_ADD(w.intime, INTERVAL 48 HOUR) THEN 1 ELSE 0 END AS extubated_within_48h_of_icu,
  CASE WHEN w.first_vent_end > TIMESTAMP_ADD(w.intime, INTERVAL 48 HOUR) THEN 1 ELSE 0 END AS extubated_after_48h_of_icu,

  COALESCE(in48.intake_ml_0_48, 0) AS intake_ml_0_48,
  COALESCE(u48.urine_ml_0_48, 0) AS urine_ml_0_48,
  COALESCE(in48.intake_ml_0_48, 0) - COALESCE(u48.urine_ml_0_48, 0) AS net_fluid_ml_0_48,
  LEAST(COALESCE(in48.intake_ml_0_48, 0) - COALESCE(u48.urine_ml_0_48, 0), 60000) AS net_fluid_ml_0_48_capped,
  SAFE_DIVIDE(LEAST(COALESCE(in48.intake_ml_0_48, 0) - COALESCE(u48.urine_ml_0_48, 0), 60000), wt.weight) AS net_fluid_ml_per_kg,

  COALESCE(inpre.intake_ml_pre_extub_or_48, 0) AS intake_ml_pre_extub_or_48,
  COALESCE(upre.urine_ml_pre_extub_or_48, 0) AS urine_ml_pre_extub_or_48,
  COALESCE(inpre.intake_ml_pre_extub_or_48, 0) - COALESCE(upre.urine_ml_pre_extub_or_48, 0) AS net_fluid_ml_pre_extub_or_48,

  COALESCE(allin.all_input_ml_0_48, 0) AS all_input_ml_0_48,
  COALESCE(nonout.nonurine_nonrrt_output_ml_0_48, 0) AS nonurine_nonrrt_output_ml_0_48,
  rrtout.rrt_output_ml_0_48_candidate,
  COALESCE(rrtout.rrt_output_candidate_records, 0) AS rrt_output_candidate_records,
  0 AS uf_capture_validated,
  COALESCE(in48.intake_unknown_unit_events, 0) AS intake_unknown_unit_events,
  COALESCE(inpre.intake_pre_unknown_unit_events, 0) AS intake_pre_unknown_unit_events,
  COALESCE(in48.intake_crosses_48h_events, 0) AS intake_crosses_48h_events,
  COALESCE(inpre.intake_pre_crosses_end_events, 0) AS intake_pre_crosses_end_events,
  CAST(NULL AS FLOAT64) AS net_fluid_comprehensive_ml_0_48_candidate,
  COALESCE(allin.all_input_ml_0_48, 0) - COALESCE(u48.urine_ml_0_48, 0) AS all_recorded_input_minus_urine_ml_0_48,

  COALESCE(r48.dialysis_present_48h, 0) AS dialysis_present_48h,

  COALESCE(rv.revent_cnt_48h, 0) AS revent_cnt_48h,
  COALESCE(rv.revent_cnt_48h_same_stay, 0) AS revent_cnt_48h_same_stay,
  COALESCE(rv.invasive_other_stay_at_index_end, 0) AS invasive_other_stay_at_index_end,
  CASE WHEN COALESCE(rv.revent_cnt_48h_same_stay, 0) = 0 THEN 1 ELSE 0 END AS weaning_success_48h_same_stay,
  CASE WHEN w.first_vent_end > TIMESTAMP_ADD(w.intime, INTERVAL 48 HOUR) THEN 1 ELSE 0 END AS extubated_strictly_after_48h_of_icu,
  CASE WHEN COALESCE(rv.revent_cnt_48h, 0) = 0 THEN 1 ELSE 0 END AS weaning_success_48h_original,

  CASE WHEN a.deathtime > w.first_vent_end
        AND a.deathtime <= TIMESTAMP_ADD(w.first_vent_end, INTERVAL WEAN_WINDOW_HOURS HOUR)
       THEN 1 ELSE 0 END AS death_within_48h_after_extubation,

  -- A missing death time is treated as no recorded death
  CASE WHEN COALESCE(rv.revent_cnt_48h, 0) = 0
        AND (a.deathtime IS NULL
             OR a.deathtime > TIMESTAMP_ADD(w.first_vent_end, INTERVAL WEAN_WINDOW_HOURS HOUR))
       THEN 1 ELSE 0 END AS weaning_success_48h_alive,

  CASE WHEN a.dischtime < TIMESTAMP_ADD(w.first_vent_end, INTERVAL WEAN_WINDOW_HOURS HOUR)
        AND (a.deathtime IS NULL OR a.deathtime > TIMESTAMP_ADD(w.first_vent_end, INTERVAL WEAN_WINDOW_HOURS HOUR))
       THEN 1 ELSE 0 END AS discharged_alive_before_48h_followup,

  CASE WHEN a.dischtime >= TIMESTAMP_ADD(w.first_vent_end, INTERVAL WEAN_WINDOW_HOURS HOUR)
        OR (a.deathtime > w.first_vent_end
            AND a.deathtime <= TIMESTAMP_ADD(w.first_vent_end, INTERVAL WEAN_WINDOW_HOURS HOUR))
       THEN 1 ELSE 0 END AS complete_inhospital_48h_followup,

  COALESCE(c.chf_icd, 0) AS chf_icd,
  COALESCE(c.cad_icd, 0) AS cad_icd,
  COALESCE(c.ckd_icd, 0) AS ckd_icd,
  COALESCE(c.sepsis_icd, 0) AS sepsis_icd,
  COALESCE(sep.sepsis3_derived, 0) AS sepsis3_derived,
  COALESCE(vp.vasopressor_use_48h, 0) AS vasopressor_use_48h,
  COALESCE(diur.diuretic_use_48h, 0) AS diuretic_use_48h,
  lac.first_lactate_0_48, lac.min_lactate_0_48, lac.max_lactate_0_48

FROM windows w
JOIN aki_48h aki ON w.stay_id = aki.stay_id
LEFT JOIN adm a ON w.hadm_id = a.hadm_id
LEFT JOIN revent rv ON w.stay_id = rv.stay_id
LEFT JOIN intake_0_48 in48 ON w.stay_id = in48.stay_id
LEFT JOIN urine_0_48 u48 ON w.stay_id = u48.stay_id
LEFT JOIN intake_pre_extub_or_48 inpre ON w.stay_id = inpre.stay_id
LEFT JOIN urine_pre_extub_or_48 upre ON w.stay_id = upre.stay_id
LEFT JOIN all_input_0_48 allin ON w.stay_id = allin.stay_id
LEFT JOIN rrt_output_0_48 rrtout ON w.stay_id = rrtout.stay_id
LEFT JOIN nonurine_nonrrt_output_0_48 nonout ON w.stay_id = nonout.stay_id
LEFT JOIN rrt_48h r48 ON w.stay_id = r48.stay_id
LEFT JOIN detail d ON w.stay_id = d.stay_id
LEFT JOIN sofa s ON w.stay_id = s.stay_id
LEFT JOIN wt ON w.stay_id = wt.stay_id
LEFT JOIN comorb c ON w.hadm_id = c.hadm_id
LEFT JOIN sepsis3 sep ON w.stay_id = sep.stay_id
LEFT JOIN vasopressor_48h vp ON w.stay_id = vp.stay_id
LEFT JOIN diuretic_48h diur ON w.stay_id = diur.stay_id
LEFT JOIN lactate_0_48 lac ON w.stay_id = lac.stay_id;

CREATE OR REPLACE TABLE `YOUR_PROJECT.YOUR_DATASET.aki_fluid_wean_exclusions` AS
SELECT * FROM cohort_all WHERE eligibility_status != 'eligible';

CREATE OR REPLACE TABLE `YOUR_PROJECT.YOUR_DATASET.aki_fluid_wean_cohort` AS
SELECT * FROM cohort_all WHERE eligibility_status = 'eligible';

-- Integrity checks
ASSERT (SELECT COUNT(*) = COUNT(DISTINCT stay_id) FROM `YOUR_PROJECT.YOUR_DATASET.aki_fluid_wean_cohort`)
  AS 'Duplicate stay_id in cohort';
ASSERT (
  SELECT COUNTIF(first_vent_end < icu_intime OR first_vent_end > icu_outtime
                 OR hospital_deathtime <= first_vent_end
                 OR weaning_success_48h_alive > weaning_success_48h_original
                 OR (death_within_48h_after_extubation = 1 AND weaning_success_48h_alive = 1)
                 OR (hospital_deathtime IS NULL AND revent_cnt_48h = 0 AND weaning_success_48h_alive != 1)) = 0
  FROM `YOUR_PROJECT.YOUR_DATASET.aki_fluid_wean_cohort`)
  AS 'Cohort failed timing or outcome checks';

-- Export this result (all rows) as CSV
SELECT * FROM `YOUR_PROJECT.YOUR_DATASET.aki_fluid_wean_cohort` ORDER BY stay_id;
