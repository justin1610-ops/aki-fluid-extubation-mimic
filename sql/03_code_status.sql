-- Code status documentation for the analytic cohort.
-- Sources: ICU charting (chartevents, items labelled "Code Status") and provider
-- orders (poe / poe_detail, field "Code status").
-- Run after 01_extract_cohort.sql.
--   Result 1: inventory of recorded values (to check the text patterns below).
--   Result 2: one row per ICU stay; export as CSV and pass to
--             analysis/run_analysis.py with --code-status.

CREATE TEMP TABLE code_status_events AS
WITH cohort AS (
  SELECT stay_id, hadm_id, icu_intime, first_vent_end
  FROM `YOUR_PROJECT.YOUR_DATASET.aki_fluid_wean_cohort`
),
chart AS (
  SELECT c.stay_id, CAST(ce.charttime AS TIMESTAMP) AS event_time,
         'chartevents' AS source, ce.value AS value
  FROM cohort c
  JOIN `physionet-data.mimiciv_3_1_icu.chartevents` ce ON ce.stay_id = c.stay_id
  JOIN `physionet-data.mimiciv_3_1_icu.d_items` di ON di.itemid = ce.itemid
  WHERE LOWER(di.label) LIKE '%code status%'
    AND ce.value IS NOT NULL
),
orders AS (
  SELECT c.stay_id, CAST(p.ordertime AS TIMESTAMP) AS event_time,
         'poe' AS source, pd.field_value AS value
  FROM cohort c
  JOIN `physionet-data.mimiciv_3_1_hosp.poe` p ON p.hadm_id = c.hadm_id
  JOIN `physionet-data.mimiciv_3_1_hosp.poe_detail` pd ON pd.poe_id = p.poe_id
  WHERE LOWER(pd.field_name) LIKE '%code status%'
    AND pd.field_value IS NOT NULL
)
SELECT e.stay_id, e.event_time, e.source, e.value,
       REGEXP_CONTAINS(LOWER(e.value), r'comfort') AS is_cmo,
       REGEXP_CONTAINS(LOWER(e.value), r'dnr|dni|do not|no cpr|cpr not') AS is_limitation
FROM (SELECT * FROM chart UNION ALL SELECT * FROM orders) e;

-- Result 1: value inventory
SELECT source, value, is_cmo, is_limitation, COUNT(*) AS records, COUNT(DISTINCT stay_id) AS stays
FROM code_status_events
GROUP BY source, value, is_cmo, is_limitation
ORDER BY stays DESC;

-- Result 2: per-stay flags
-- Documentation within 6 h after ventilator discontinuation is counted as
-- "by extubation" to allow for charting delay.
SELECT
  c.stay_id,
  IF(COUNTIF(e.stay_id IS NOT NULL) > 0, 1, 0) AS any_code_status_documented,
  IF(COUNTIF(e.is_cmo AND e.event_time <= TIMESTAMP_ADD(c.first_vent_end, INTERVAL 6 HOUR)) > 0, 1, 0) AS cmo_by_extubation,
  IF(COUNTIF(e.is_cmo AND e.event_time <= TIMESTAMP_ADD(c.first_vent_end, INTERVAL 48 HOUR)) > 0, 1, 0) AS cmo_by_48h_after,
  IF(COUNTIF((e.is_cmo OR e.is_limitation) AND e.event_time <= TIMESTAMP_ADD(c.first_vent_end, INTERVAL 6 HOUR)) > 0, 1, 0) AS dnr_dni_or_cmo_by_extubation
FROM `YOUR_PROJECT.YOUR_DATASET.aki_fluid_wean_cohort` c
LEFT JOIN code_status_events e ON e.stay_id = c.stay_id
GROUP BY c.stay_id
ORDER BY c.stay_id;
