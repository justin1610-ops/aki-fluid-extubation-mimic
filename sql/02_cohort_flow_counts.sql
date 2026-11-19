-- Counts for the cohort flow diagram (Figure 1).
-- Run after 01_extract_cohort.sql. Export the result as CSV (columns: step, n)
-- and pass it to analysis/run_analysis.py with --flow.

WITH icu AS (
  SELECT stay_id, CAST(intime AS TIMESTAMP) AS intime
  FROM `physionet-data.mimiciv_3_1_icu.icustays`
),
aki AS (
  SELECT DISTINCT k.stay_id
  FROM `physionet-data.mimiciv_3_1_derived.kdigo_stages` k
  JOIN icu i USING (stay_id)
  WHERE CAST(k.charttime AS TIMESTAMP) >= i.intime
    AND CAST(k.charttime AS TIMESTAMP) < TIMESTAMP_ADD(i.intime, INTERVAL 48 HOUR)
    AND k.aki_stage >= 1
),
imv AS (
  SELECT DISTINCT stay_id
  FROM `physionet-data.mimiciv_3_1_derived.ventilation`
  WHERE ventilation_status = 'InvasiveVent' AND starttime IS NOT NULL
)
SELECT '1_all_icu_stays' AS step, COUNT(*) AS n FROM icu
UNION ALL SELECT '2_aki_within_48h', COUNT(*) FROM aki
UNION ALL SELECT '3_aki_and_invasivevent', COUNT(*) FROM aki JOIN imv USING (stay_id)
UNION ALL
SELECT CONCAT('4_excluded_', eligibility_status), COUNT(*)
FROM `YOUR_PROJECT.YOUR_DATASET.aki_fluid_wean_exclusions`
GROUP BY eligibility_status
UNION ALL
SELECT '5_final_eligible', COUNT(*) FROM `YOUR_PROJECT.YOUR_DATASET.aki_fluid_wean_cohort`
ORDER BY step;
