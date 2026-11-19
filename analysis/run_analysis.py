"""Reproduce all analyses, tables, and figures of the manuscript.

Input : CSV exported from sql/01_extract_cohort.sql (one row per ICU stay).
Optional inputs:
  --flow  CSV exported from sql/02_cohort_flow_counts.sql (columns: step, n)
  --code-status  CSV exported from sql/03_code_status.sql (one row per ICU stay)
Output: tables (CSV), figures (PNG/TIFF), and results.json in --out.

Usage:
  python analysis/run_analysis.py --data cohort.csv --out results \
      [--flow flow_counts.csv] [--code-status code_status.csv]
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from matplotlib.patches import FancyBboxPatch
from patsy import dmatrix
from scipy import stats

plt.rcParams["font.family"] = "DejaVu Sans"

QS = ["Q1", "Q2", "Q3", "Q4"]
QT = 'C(Q, Treatment("Q1"))'
BASE = "admission_age + male + first_day_sofa + weight + max_aki_stage_48h"
EXPANDED = (" + chf_icd + cad_icd + ckd_icd + sepsis3_derived + vasopressor_use_48h"
            " + diuretic_use_48h + lactate + lactate_missing")


# ---------------------------------------------------------------- data
def load(path):
    d = pd.read_csv(path)
    for c in ["icu_intime", "first_vent_end", "hospital_deathtime"]:
        d[c] = pd.to_datetime(d[c].astype("string").str.replace(" UTC", ""), errors="coerce")
    d["male"] = (d["gender"] == "M").astype(int)
    d["bal"] = d["net_fluid_ml_0_48_capped"] / 1000            # primary exposure, L
    d["bal_uncapped"] = d["net_fluid_ml_0_48"] / 1000
    d["bal_pre"] = d["net_fluid_ml_pre_extub_or_48"] / 1000     # window ending at extubation or 48 h
    d["bal_comp"] = (d["all_input_ml_0_48"] - d["urine_ml_0_48"]
                     - d["nonurine_nonrrt_output_ml_0_48"]) / 1000
    d["bal_kg"] = d["net_fluid_ml_per_kg"]
    # Quartile cutoffs from the full analytic cohort; applied unchanged in all subgroups.
    cuts = d["bal"].quantile([0.25, 0.5, 0.75]).values
    d["Q"] = pd.cut(d["bal"], [-np.inf, *cuts, np.inf], labels=QS).astype(str)
    # Outcomes
    d["Y"] = d["weaning_success_48h_same_stay"]                 # primary: no reventilation, same ICU stay
    d["D48"] = d["death_within_48h_after_extubation"]
    d["YD"] = ((d["Y"] == 1) & (d["D48"] == 0)).astype(int)     # death within 48 h counted as failure
    d["Y_hosp"] = d["weaning_success_48h_original"]             # reventilation in any ICU stay of the admission
    d["late"] = d["extubated_strictly_after_48h_of_icu"].astype(int)
    d["lactate_missing"] = d["max_lactate_0_48"].isna().astype(int)
    return d, cuts


def complete_cases(d):
    cc = d.dropna(subset=["weight"]).copy()
    cc["lactate"] = cc["max_lactate_0_48"].fillna(cc["max_lactate_0_48"].median())
    return cc


# ---------------------------------------------------------------- models
def logit(df, y, x, rrt=True, extra="", cluster=False):
    f = f"{y} ~ {x} + {BASE}" + (" + dialysis_present_48h" if rrt else "") + extra
    kw = dict(cov_type="cluster", cov_kwds={"groups": df["subject_id"]}) if cluster else {}
    return smf.logit(f, df).fit(disp=0, **kw)


def odds_ratio(m, term):
    lo, hi = np.exp(m.conf_int().loc[term])
    return {"or": float(np.exp(m.params[term])), "lo": float(lo), "hi": float(hi), "p": float(m.pvalues[term])}


def quartile_ors(m):
    return {q: odds_ratio(m, f'C(Q, Treatment("Q1"))[T.{q}]') for q in QS[1:]}


def standardized(m, df):
    """Standardized adjusted probabilities by quartile and differences vs Q1 (delta method)."""
    beta, V = m.params.values, m.cov_params().values
    rhs = m.model.formula.split("~")[1]

    def design(q):
        return np.asarray(dmatrix(rhs, df.assign(Q=pd.Categorical([q] * len(df), categories=QS)),
                                  return_type="dataframe"))
    x1 = design("Q1"); p1 = 1 / (1 + np.exp(-x1 @ beta)); g1 = (x1 * (p1 * (1 - p1))[:, None]).mean(0)
    out = {"Q1": {"prob": float(p1.mean() * 100)}}
    for q in QS[1:]:
        x = design(q); p = 1 / (1 + np.exp(-x @ beta)); g = (x * (p * (1 - p))[:, None]).mean(0) - g1
        diff = (p.mean() - p1.mean()) * 100; se = np.sqrt(g @ V @ g) * 100
        out[q] = {"prob": float(p.mean() * 100), "diff": float(diff),
                  "lo": float(diff - 1.96 * se), "hi": float(diff + 1.96 * se)}
    return out


# ---------------------------------------------------------------- splines
def rcs_basis(x, knots):
    p = lambda u: np.maximum(u, 0) ** 3
    k = knots; n = len(k); cols = [x]
    for j in range(n - 2):
        cols.append((p(x - k[j]) - p(x - k[n - 2]) * (k[n - 1] - k[j]) / (k[n - 1] - k[n - 2])
                     + p(x - k[n - 1]) * (k[n - 2] - k[j]) / (k[n - 1] - k[n - 2])) / (k[n - 1] - k[0]) ** 2)
    return np.column_stack(cols)


def spline_figure(df, rrt, title, path):
    knots = np.percentile(df["bal"], [5, 35, 65, 95])
    cov = ["admission_age", "male", "first_day_sofa", "weight", "max_aki_stage_48h"] + (["dialysis_present_48h"] if rrt else [])
    X = sm.add_constant(np.column_stack([rcs_basis(df["bal"].values, knots), df[cov].values.astype(float)]))
    m = sm.Logit(df["Y"].values, X).fit(disp=0); b, V = m.params, m.cov_params()
    nl = [2, 3]
    p_nonlinear = float(stats.chi2.sf(b[nl] @ np.linalg.solve(V[np.ix_(nl, nl)], b[nl]), 2))
    ref = {"admission_age": round(df["admission_age"].median()), "male": 1,
           "first_day_sofa": round(df["first_day_sofa"].median()), "weight": round(df["weight"].median()),
           "max_aki_stage_48h": 2, "dialysis_present_48h": 0}
    lo, hi = np.percentile(df["bal"], [1, 99]); xs = np.linspace(lo, hi, 200)
    Xg = np.column_stack([np.ones(len(xs)), rcs_basis(xs, knots), np.tile([ref[c] for c in cov], (len(xs), 1))])
    eta = Xg @ b; se = np.sqrt(np.einsum("ij,jk,ik->i", Xg, V, Xg)); f = lambda z: 1 / (1 + np.exp(-z))
    fig, ax = plt.subplots(figsize=(8, 5.7))
    ax.fill_between(xs, f(eta - 1.96 * se), f(eta + 1.96 * se), color="#A6A6A6", lw=0)
    ax.plot(xs, f(eta), color="black", lw=2.5)
    ax.set_title(title, loc="left", fontsize=17)
    ax.set_xlabel("Net fluid balance in first 48 hours (L)", fontsize=14)
    ax.set_ylabel("Predicted probability of successful extubation", fontsize=14)
    ax.grid(color="#EBEBEB"); ax.set_axisbelow(True)
    adj = " ".join(f"{'gender' if c == 'male' else c}={'M' if c == 'male' else ref[c]}" for c in cov)
    fig.text(0.985, 0.012, "Adjusted to: " + adj, ha="right", fontsize=8.3)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path.with_suffix(".png"), dpi=200)
    fig.savefig(path.with_suffix(".tiff"), dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return {"p_nonlinear": p_nonlinear, "knots": knots.tolist()}


def flow_figure(counts, path):
    c = dict(zip(counts["step"], counts["n"]))
    ex = {k[len("4_excluded_"):]: v for k, v in c.items() if k.startswith("4_excluded_")}
    outside_icu = ex.get("episode_end_after_ICU", 0) + ex.get("episode_end_before_ICU", 0)
    lines = ["Excluded (n = {:,}):".format(sum(ex.values())),
             f"• Episode ending outside the index ICU stay (n = {outside_icu})",
             f"• Missing ICU discharge time (n = {ex.get('missing_ICU_outtime', 0)})",
             f"• Episode ending outside the hospital admission (n = {ex.get('episode_end_outside_hospital', 0)})",
             f"• In-hospital death without recorded death time (n = {ex.get('hospital_death_time_missing', 0)})",
             f"• Death at or before the episode end (n = {ex.get('death_at_or_before_episode_end', 0)})"]
    blue = dict(fc="#D6E4F0", ec="#2E75B6"); orange = dict(fc="#FBE5D6", ec="#C55A11")
    fig, ax = plt.subplots(figsize=(11, 7.6)); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    def box(x, y, w, h, text, style, ha="center", fs=10.5):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1", lw=2, **style))
        ax.text(x + w / 2 if ha == "center" else x + 1.5, y + h / 2, "\n".join(text), ha=ha, va="center", fontsize=fs, linespacing=1.45)
    n1, n2, n3, n4 = c["1_all_icu_stays"], c["2_aki_within_48h"], c["3_aki_and_invasivevent"], c["5_final_eligible"]
    mains = [(89, 10, ["Total ICU admissions in MIMIC-IV v3.1", f"n = {n1:,}"]),
             (67, 10, ["ICU admissions with AKI within 48 h", f"n = {n2:,}"]),
             (41, 13, ["ICU admissions with AKI and a qualifying invasive", "mechanical ventilation episode (InvasiveVent)", f"n = {n3:,}"]),
             (1, 16, ["Final analytic cohort", "Early AKI, invasive mechanical ventilation,", "and documented first extubation attempt", f"n = {n4:,}"])]
    for y, h, t in mains: box(1, y, 49, h, t, blue)
    for (y1, _, _), (y2, h2, _) in zip(mains, mains[1:]):
        ax.annotate("", xy=(25.5, y2 + h2), xytext=(25.5, y1), arrowprops=dict(arrowstyle="-|>", color="#2E75B6", lw=2))
    excl = [(78, 9, ["Excluded: No AKI within 48 h", f"n = {n1 - n2:,}"], "center"),
            (52, 16, ["Excluded: No qualifying invasive mechanical", "ventilation episode (InvasiveVent)", "with documented termination", f"n = {n2 - n3:,}"], "center"),
            (17, 22, lines, "left")]
    for y, h, t, ha in excl:
        box(52, y, 47.5, h, t, orange, ha=ha, fs=9.3 if ha == "left" else 10.5)
        ax.annotate("", xy=(52, y + h / 2), xytext=(25.5, y + h / 2), arrowprops=dict(arrowstyle="-|>", color="#C55A11", lw=2))
    fig.savefig(path.with_suffix(".png"), dpi=200, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(path.with_suffix(".tiff"), dpi=300, bbox_inches="tight", pad_inches=0.05, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


# ---------------------------------------------------------------- tables
def table1(d):
    def med(s, fmt="{:,.0f}", sep="–"):
        txt = f"{fmt.format(s.median())} ({fmt.format(s.quantile(.25))}{sep}{fmt.format(s.quantile(.75))})"
        return txt.replace("-", "−")

    def npct(s):
        return f"{int(s.sum()):,} ({s.mean() * 100:.1f}%)"
    kw = lambda c: stats.kruskal(*[g[c].dropna() for _, g in d.groupby("Q")]).pvalue
    chi = lambda c: stats.chi2_contingency(pd.crosstab(d["Q"], d[c]))[1]
    G = {q: g for q, g in d.groupby("Q")}
    d = d.assign(female=(d["gender"] == "F").astype(int))
    G = {q: g for q, g in d.groupby("Q")}
    rows = [["Age, years"] + [med(G[q]["admission_age"]) for q in QS] + [kw("admission_age")],
            ["Female sex, n (%)"] + [npct(G[q]["female"]) for q in QS] + [chi("female")],
            ["SOFA score (day 1)"] + [med(G[q]["first_day_sofa"], "{:.1f}") for q in QS] + [kw("first_day_sofa")],
            ["Body weight, kg"] + [med(G[q]["weight"]) for q in QS] + [kw("weight")]]
    for lab, col in [("Congestive heart failure", "chf_icd"), ("Coronary artery disease", "cad_icd"),
                     ("Chronic kidney disease", "ckd_icd"), ("Sepsis-3", "sepsis3_derived"),
                     ("Vasopressor use within 48 h", "vasopressor_use_48h"), ("Diuretic use within 48 h", "diuretic_use_48h"),
                     ("RRT within 48 h", "dialysis_present_48h")]:
        rows.append([lab + ", n (%)"] + [npct(G[q][col]) for q in QS] + [chi(col)])
    rows.append(["Max AKI stage within 48 h, n (%)", "", "", "", "", chi("max_aki_stage_48h")])
    for s in (1, 2, 3):
        rows.append([f"  Stage {s}"] + [npct(G[q]["max_aki_stage_48h"] == s) for q in QS] + [""])
    rows.append(["Total IV intake, mL (0–48 h)"] + [med(G[q]["intake_ml_0_48"]) for q in QS] + [kw("intake_ml_0_48")])
    rows.append(["Total urine output, mL (0–48 h)"] + [med(G[q]["urine_ml_0_48"]) for q in QS] + [kw("urine_ml_0_48")])
    rows.append(["Net fluid balance, L (0–48 h)"] + [med(G[q]["bal"], "{:.2f}", " to ") for q in QS] + [kw("bal")])
    cols = ["Characteristic"] + [f"{q} (n={len(G[q]):,})" for q in QS] + ["P"]
    t = pd.DataFrame(rows, columns=cols)
    t["P"] = t["P"].map(lambda v: v if v == "" else ("<0.001" if v < 0.001 else f"{v:.3f}"))
    return t


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--flow")
    ap.add_argument("--code-status")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    d, cuts = load(a.data)
    cc = complete_cases(d)
    R = {"n_cohort": len(d), "n_subjects": int(d["subject_id"].nunique()), "n_complete_case": len(cc),
         "quartile_cutoffs_L": cuts.tolist(), "n_success": int(d["Y"].sum()),
         "missing_weight": int(d["weight"].isna().sum()), "missing_lactate": int(d["max_lactate_0_48"].isna().sum()),
         "hours_to_extubation": {"median": float(d["hours_icu_to_extubation"].median()),
                                 "q1": float(d["hours_icu_to_extubation"].quantile(.25)),
                                 "q3": float(d["hours_icu_to_extubation"].quantile(.75)),
                                 "mean": float(d["hours_icu_to_extubation"].mean())},
         "extubated_within_48h": int(d["extubated_within_48h_of_icu"].sum()),
         "extubated_after_48h_by_quartile_pct": d.groupby("Q")["late"].mean().mul(100).to_dict(),
         "crude_success_pct": d.groupby("Q")["Y"].mean().mul(100).to_dict(),
         "crude_success_death_as_failure_pct": d.groupby("Q")["YD"].mean().mul(100).to_dict(),
         "rrt_marker_full_cohort": int(d["dialysis_present_48h"].sum()),
         "rrt_marker_complete_case": int(cc["dialysis_present_48h"].sum())}

    nr = cc[cc["dialysis_present_48h"] == 0]
    late = cc[cc["late"] == 1]
    early = cc[cc["late"] == 0]
    first = cc[cc["is_first_icu_stay"] == 1]
    R["n_no_rrt_full"] = int((d["dialysis_present_48h"] == 0).sum())
    R["n_no_rrt"] = len(nr); R["n_late_full"] = int(d["late"].sum()); R["n_late"] = len(late)
    R["n_late_no_rrt"] = int((late["dialysis_present_48h"] == 0).sum()); R["n_first_stay"] = len(first)

    # Primary and death-as-failure models
    for y in ("Y", "YD"):
        m = logit(cc, y, QT)
        R[f"{y}_quartiles"] = quartile_ors(m)
        R[f"{y}_standardized"] = standardized(m, cc)
        R[f"{y}_per_L"] = odds_ratio(logit(cc, y, "bal"), "bal")
        m = logit(nr, y, QT, rrt=False)
        R[f"{y}_no_rrt_quartiles"] = quartile_ors(m)
        R[f"{y}_no_rrt_standardized"] = standardized(m, nr)
        R[f"{y}_no_rrt_per_L"] = odds_ratio(logit(nr, y, "bal", rrt=False), "bal")
        R[f"{y}_late_quartiles"] = quartile_ors(logit(late, y, QT))
        R[f"{y}_late_per_L"] = odds_ratio(logit(late, y, "bal"), "bal")
        R[f"{y}_late_no_rrt_quartiles"] = quartile_ors(logit(late[late["dialysis_present_48h"] == 0], y, QT, rrt=False))
        R[f"{y}_first_stay_quartiles"] = quartile_ors(logit(first, y, QT))
        R[f"{y}_expanded_quartiles"] = quartile_ors(logit(cc, y, QT, extra=EXPANDED))
    R["unadjusted_quartiles"] = quartile_ors(smf.logit(f"Y ~ {QT}", cc).fit(disp=0))

    # Sensitivity analyses (Supplementary Table S1)
    R["S1"] = {
        "uncapped_per_L": odds_ratio(logit(cc, "Y", "bal_uncapped"), "bal_uncapped"),
        "weight_indexed_per_mL_kg": odds_ratio(logit(cc, "Y", "bal_kg"), "bal_kg"),
        "rrt_interaction": odds_ratio(logit(cc, "Y", "bal * dialysis_present_48h", rrt=False), "bal:dialysis_present_48h"),
        "cluster_robust_per_L": odds_ratio(logit(cc, "Y", "bal", cluster=True), "bal"),
        "first_stay_per_L": odds_ratio(logit(first, "Y", "bal"), "bal"),
        "after_48h_per_L": R["Y_late_per_L"],
        "expanded_per_L": odds_ratio(logit(cc, "Y", "bal", extra=EXPANDED), "bal"),
        "pre_extubation_window_per_L": odds_ratio(logit(cc, "Y", "bal_pre"), "bal_pre"),
        "comprehensive_per_L": odds_ratio(logit(cc, "Y", "bal_comp"), "bal_comp"),
        "hospital_wide_reventilation_per_L": odds_ratio(logit(cc, "Y_hosp", "bal"), "bal"),
        "no_rrt_per_L": R["Y_no_rrt_per_L"],
    }
    R["hospital_wide_reventilation_quartiles"] = quartile_ors(logit(cc, "Y_hosp", QT))
    R["additional_failures_hospital_wide"] = int(((d["Y"] == 1) & (d["Y_hosp"] == 0)).sum())

    # Extubation timing
    R["within_48h_per_L"] = odds_ratio(logit(early, "Y", "bal"), "bal")
    m = logit(cc, "Y", "bal * late")
    b, V = m.params, m.cov_params()

    def lincomb(w):
        e = sum(b[k] * v for k, v in w.items())
        se = np.sqrt(sum(V.loc[i, j] * w[i] * w[j] for i in w for j in w))
        return {"or": float(np.exp(e)), "lo": float(np.exp(e - 1.96 * se)), "hi": float(np.exp(e + 1.96 * se))}
    R["timing_common_covariates"] = {"within_48h": lincomb({"bal": 1}), "after_48h": lincomb({"bal": 1, "bal:late": 1}),
                                     "p_interaction": float(m.pvalues["bal:late"])}
    full = smf.logit(f"Y ~ late * (bal + {BASE} + dialysis_present_48h)", cc).fit(disp=0)
    R["timing_all_terms_interacted_p"] = float(full.pvalues["late:bal"])

    # Early post-extubation deaths
    hrs = (d["hospital_deathtime"] - d["first_vent_end"]).dt.total_seconds() / 3600
    q4dead = (d["Q"] == "Q4") & (d["D48"] == 1)
    R["early_deaths"] = {
        "n": int(d["D48"].sum()),
        "by_quartile": {q: [int(g["D48"].sum()), len(g)] for q, g in d.groupby("Q")},
        "q4_within_6h_pct": float((hrs[q4dead] <= 6).mean() * 100),
        "without_reventilation": int(((d["D48"] == 1) & (d["revent_cnt_48h_same_stay"] == 0)).sum()),
        "discharged_alive_before_48h": int(d["discharged_alive_before_48h_followup"].sum())}
    R["excluding_early_deaths_quartiles"] = quartile_ors(logit(cc[cc["D48"] == 0], "Y", QT))

    # Code status (optional)
    if a.code_status:
        cs = pd.read_csv(a.code_status)
        dd = d.merge(cs, on="stay_id", how="left")
        for c in [c for c in cs.columns if c != "stay_id"]:
            if dd[c].dtype != object:
                dd[c] = dd[c].fillna(0)
        dead = dd["D48"] == 1
        summ = {}
        for flag in [c for c in cs.columns if c.startswith(("cmo_", "dnr_", "any_"))]:
            summ[flag] = {
                "all_stays_pct": float(dd[flag].mean() * 100),
                "early_deaths": {q: [int(g.loc[g["D48"] == 1, flag].sum()), int((g["D48"] == 1).sum())] for q, g in dd.groupby("Q")},
                "survivors_pct": float(dd.loc[~dead, flag].mean() * 100)}
            summ[flag]["all_stays_by_quartile_pct"] = {q: float(g[flag].mean() * 100) for q, g in dd.groupby("Q")}
        R["code_status"] = summ
        for flag in ["cmo_by_extubation", "dnr_dni_or_cmo_by_extubation"]:
            if flag in dd.columns:
                keep = dd.loc[dd[flag] == 0, "stay_id"]
                sub = cc[cc["stay_id"].isin(keep)]
                R[f"excluding_{flag}"] = {
                    "n": len(sub),
                    "early_death_pct_by_quartile": {q: float(g["D48"].mean() * 100) for q, g in sub.groupby("Q")},
                    "Y_quartiles": quartile_ors(logit(sub, "Y", QT)),
                    "YD_quartiles": quartile_ors(logit(sub, "YD", QT))}

    # Tables and figures
    table1(d).to_csv(out / "table1.csv", index=False)
    R["spline_full"] = spline_figure(cc, True, "Restricted cubic spline analysis", out / "Figure2")
    R["spline_no_rrt"] = spline_figure(nr, False, "Sensitivity analysis excluding patients receiving RRT", out / "FigureS1")
    if a.flow:
        flow_figure(pd.read_csv(a.flow), out / "Figure1")
    (out / "results.json").write_text(json.dumps(R, indent=2, default=float))
    print(f"Done. Results written to {out}/")


if __name__ == "__main__":
    main()
