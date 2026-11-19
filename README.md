# Early net fluid balance and extubation success in ventilated adults with AKI (MIMIC-IV)

Code for the cohort extraction and statistical analyses of:

> Kim M, Lee H-C, Lee B. Early Net Fluid Balance and Extubation Success in Mechanically Ventilated Adults with Acute Kidney Injury: A Retrospective Cohort Study Using MIMIC-IV. *Journal of Clinical Medicine* (manuscript jcm-4580114).

## Data

The study uses MIMIC-IV v3.1. The data are not included in this repository. Access requires credentialed approval on PhysioNet (https://physionet.org/content/mimiciv/3.1/), completion of the required training, and signing of the data use agreement. The SQL queries run on the MIMIC-IV tables hosted on Google BigQuery (`physionet-data.mimiciv_3_1_*`), including the derived concept tables (`mimiciv_3_1_derived`).

## Repository structure

```
sql/
  01_extract_cohort.sql       Cohort selection, exposure, outcome, and covariate extraction
  02_cohort_flow_counts.sql   Counts for the cohort flow diagram (Figure 1)
  03_code_status.sql          Code status / comfort-measures-only documentation
analysis/
  run_analysis.py             All statistical analyses, tables, and figures
requirements.txt              Python package versions used for the published results
```

## How to reproduce

1. In BigQuery, replace `YOUR_PROJECT.YOUR_DATASET` in the SQL files with a dataset you can write to.
2. Run `sql/01_extract_cohort.sql` as a single query job. Export the final result (all rows) as `cohort.csv`.
3. Run `sql/02_cohort_flow_counts.sql` and export the result as `flow_counts.csv`.
4. Run `sql/03_code_status.sql` and export the second result as `code_status.csv`.
5. Run the analysis (Python 3.13):

```
pip install -r requirements.txt
python analysis/run_analysis.py --data cohort.csv --flow flow_counts.csv \
    --code-status code_status.csv --out results
```

`results/` will contain `results.json` (all estimates reported in the manuscript), `table1.csv`, and Figures 1, 2, and S1 (PNG and TIFF).

Expected cohort size: 25,199 ICU stays from 21,907 patients (25,012 with complete covariate data).

## Definitions

- **Cohort:** adult ICU stays with KDIGO stage ≥1 AKI within 48 h of ICU admission and a first invasive ventilation episode (`ventilation_status = 'InvasiveVent'`) ending within the index ICU stay and hospital admission, without in-hospital death at or before the end of that episode.
- **Exposure:** selected intravenous fluids and blood products minus urine output during the first 48 h after ICU admission, capped at +60 L. Infusion volumes are assigned by start time. Quartile cutoffs are computed in the full cohort and applied unchanged in all subgroup analyses.
- **Primary outcome:** no new invasive ventilation episode within 48 h after the end of the first episode in the same ICU stay.
- **Secondary outcome:** primary outcome and alive at 48 h after extubation.
- **Models:** logistic regression adjusted for age, sex, first-day SOFA, body weight, maximum AKI stage within 48 h, and RRT within 48 h; complete-case analysis.

## Notes

- The analyses in the original submission were run in R. All results in the revised and published manuscript were produced with the Python code in this repository.
- Ultrafiltration volumes are not included in the fluid balance (see manuscript).
- `first_lactate_0_48` can differ between runs when two lactate values share the same timestamp; it is not used in the analyses.

## License

MIT License (see `LICENSE`).
