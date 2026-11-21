"""Risk-based quality management: key risk indicators by site and centralized
statistical monitoring.

Key risk indicators (KRIs), each compared with the rest of the study:
  * adverse events per subject-month of exposure (under-reporting risk)
  * open queries per 100 data points (data quality)
  * protocol deviations per subject
  * median days from visit to data entry (data currency)
  * lab results in a non-standard unit

Centralized statistical monitoring looks for data that are too regular to
be real: terminal-digit preference in blood pressure and unusually low
within-subject variability of blood pressure and pulse across visits.

Status: red when a site differs from the rest of the study beyond a
statistical threshold and a practical one; amber for smaller signals.
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def kris(raw, queries, deviations):
    dm, ae, ex, sv, vs, lb = raw["dm"], raw["ae"], raw["ex"], raw["sv"], raw["vs"], raw["lb"]
    site = dm.set_index("SUBJID").SITEID
    ex = ex.assign(SITEID=ex.SUBJID.map(site),
                   months=(pd.to_datetime(ex.EXENDAT) - pd.to_datetime(ex.EXSTDAT)).dt.days / 30.44)
    exposure = ex.groupby("SITEID").months.sum()
    n_ae = ae.SUBJID.map(site).value_counts().reindex(exposure.index, fill_value=0)
    datapoints = (vs.SUBJID.map(site).value_counts() * 4 + lb.SUBJID.map(site).value_counts()).reindex(exposure.index)
    nq = queries.SITEID.astype(str).str.zfill(2).value_counts().reindex(exposure.index, fill_value=0)
    nd = deviations.SITEID.astype(str).str.zfill(2).value_counts().reindex(exposure.index, fill_value=0)
    subjects = dm.SITEID.value_counts().reindex(exposure.index)
    lag = (pd.to_datetime(sv.ENTRYDAT) - pd.to_datetime(sv.SVSTDAT)).dt.days
    lag_med = lag.groupby(sv.SUBJID.map(site)).median()
    nonstd = lb[(lb.LBTEST == "CREAT") & (lb.LBORRESU == "umol/L")].SUBJID.map(site).value_counts().reindex(
        exposure.index, fill_value=0)

    k = pd.DataFrame({"subjects": subjects, "exposure_months": exposure.round(1), "adverse_events": n_ae,
                      "ae_per_subject_month": (n_ae / exposure).round(3),
                      "queries_per_100_datapoints": (100 * nq / datapoints).round(2),
                      "deviations_per_subject": (nd / subjects).round(2),
                      "median_entry_lag_days": lag_med, "nonstandard_lab_units": nonstd})
    # AE under-reporting: Poisson test of the site's count against the rate at all other sites.
    p_ae, rr = [], []
    for s in k.index:
        rest = k.drop(s)
        rate = rest.adverse_events.sum() / rest.exposure_months.sum()
        expected = rate * k.loc[s, "exposure_months"]
        rr.append(k.loc[s, "adverse_events"] / expected)
        p_ae.append(stats.poisson.cdf(k.loc[s, "adverse_events"], expected))
    k["ae_rate_ratio_vs_others"] = np.round(rr, 2)
    k["ae_underreporting_p"] = p_ae
    other_lag = [lag[sv.SUBJID.map(site) != s].median() for s in k.index]
    k["entry_lag_vs_others"] = (k.median_entry_lag_days / np.array(other_lag)).round(2)
    return k


def csm(raw):
    dm, vs = raw["dm"], raw["vs"].drop_duplicates(["SUBJID", "VISIT"])
    vs = vs[(vs.SYSBP.between(70, 220)) & (vs.SYSBP > vs.DIABP)]       # exclude values already under query
    site = dm.set_index("SUBJID").SITEID
    vs = vs.assign(SITEID=vs.SUBJID.map(site))
    rows = []
    sd = vs.groupby("SUBJID").SYSBP.std().rename("sd").to_frame().assign(SITEID=lambda d: d.index.map(site))
    sdp = vs.groupby("SUBJID").PULSE.std().rename("sd").to_frame().assign(SITEID=lambda d: d.index.map(site))
    for s, g in vs.groupby("SITEID"):
        digits = pd.concat([g.SYSBP, g.DIABP]).astype(int) % 10
        obs = np.bincount(digits, minlength=10)
        chi = stats.chisquare(obs)
        mine, others = sd[sd.SITEID == s].sd.dropna(), sd[sd.SITEID != s].sd.dropna()
        mw = stats.mannwhitneyu(mine, others, alternative="less")
        pm, po = sdp[sdp.SITEID == s].sd.dropna(), sdp[sdp.SITEID != s].sd.dropna()
        mwp = stats.mannwhitneyu(pm, po, alternative="less")
        rows.append({"SITEID": s, "bp_values": int(obs.sum()), "share_ending_in_0": round(obs[0] / obs.sum(), 3),
                     "digit_preference_p": chi.pvalue,
                     "median_within_subject_sbp_sd": round(mine.median(), 2),
                     "others_median_sd": round(others.median(), 2), "low_variability_p": mw.pvalue,
                     "median_within_subject_pulse_sd": round(pm.median(), 2),
                     "others_median_pulse_sd": round(po.median(), 2), "low_pulse_variability_p": mwp.pvalue})
    return pd.DataFrame(rows).set_index("SITEID")


def status(k, c, n_sites):
    alpha = 0.05 / n_sites
    st = pd.DataFrame(index=k.index)
    st["AE reporting"] = np.select([(k.ae_underreporting_p < alpha) & (k.ae_rate_ratio_vs_others < 0.5),
                                    k.ae_rate_ratio_vs_others < 0.7], ["red", "amber"], "green")
    st["Data entry lag"] = np.select([k.entry_lag_vs_others > 2, k.entry_lag_vs_others > 1.5], ["red", "amber"], "green")
    st["Query rate"] = np.select([k.queries_per_100_datapoints > 2 * k.queries_per_100_datapoints.median(),
                                  k.queries_per_100_datapoints > 1.5 * k.queries_per_100_datapoints.median()],
                                 ["red", "amber"], "green")
    st["Deviations"] = np.select([k.deviations_per_subject > 2 * k.deviations_per_subject.median(),
                                  k.deviations_per_subject > 1.5 * k.deviations_per_subject.median()],
                                 ["red", "amber"], "green")
    st["Lab units"] = np.where(k.nonstandard_lab_units > 0, "amber", "green")
    st["Digit preference"] = np.where((c.digit_preference_p < alpha) & (c.share_ending_in_0 > 0.3), "red", "green")
    # Either blood pressure or pulse varying across visits at less than half the rate seen elsewhere.
    low_sbp = (c.low_variability_p < alpha) & (c.median_within_subject_sbp_sd < 0.5 * c.others_median_sd)
    low_hr = (c.low_pulse_variability_p < alpha) & (c.median_within_subject_pulse_sd < 0.5 * c.others_median_pulse_sd)
    st["Low variability"] = np.where(low_sbp | low_hr, "red", "green")
    return st


def heatmap(st):
    colors = {"green": 0, "amber": 1, "red": 2}
    m = st.replace(colors).astype(int).values
    fig, ax = plt.subplots(figsize=(8, 4))
    cmap = matplotlib.colors.ListedColormap(["#4caf50", "#ffb300", "#e53935"])
    ax.imshow(m, cmap=cmap, vmin=0, vmax=2, aspect="auto")
    ax.set_xticks(range(st.shape[1]), st.columns, rotation=30, ha="right")
    ax.set_yticks(range(st.shape[0]), [f"Site {s}" for s in st.index])
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            ax.text(j, i, st.iat[i, j], ha="center", va="center", fontsize=7, color="white" if m[i, j] != 1 else "black")
    ax.set_title("Risk indicators by site")
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "results", "figures", "site_risk_heatmap.png"), dpi=150)
    plt.close(fig)


def main(raw, queries, deviations):
    k = kris(raw, queries, deviations)
    c = csm(raw)
    st = status(k, c, len(k))
    k.to_csv(os.path.join(ROOT, "results", "tables", "site_kris.csv"))
    c.to_csv(os.path.join(ROOT, "results", "tables", "central_statistical_monitoring.csv"))
    st.to_csv(os.path.join(ROOT, "results", "tables", "site_risk_status.csv"))
    heatmap(st)
    print(k.to_string())
    print(c.to_string())
    print(st.to_string())
    return k, c, st
