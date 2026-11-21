"""Programmed edit checks that raise data queries, as an EDC system would.

The check catalog is in specs/edit_checks.csv. Each check returns the
records it fails; together they form the query listing that data management
sends to sites. Checks are then scored against the errors injected into the
simulated data, so we know how many real problems they catch (recall) and how
many queries they raise on clean data (false positives).
"""

import json
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_raw(data):
    raw = os.path.join(data, "raw")
    read = lambda n: pd.read_csv(os.path.join(raw, f"{n}.csv"), dtype={"SITEID": str})
    return {n: read(n) for n in ["dm", "ie", "sv", "vs", "lb", "ae", "ex", "ds"]}


def run_checks(d):
    vs, ae, lb, sv, dm, ex = d["vs"], d["ae"], d["lb"], d["sv"], d["dm"], d["ex"]
    q = []

    def add(check, df, form, field, value_col, msg):
        for _, r in df.iterrows():
            q.append({"check_id": check, "SUBJID": r.SUBJID, "form": form,
                      "visit_or_record": r.get("VISIT", r.get("AETERM", "")), "field": field,
                      "value": r[value_col] if value_col else "", "query_text": msg})

    add("EC01", vs[(vs.SYSBP < 70) | (vs.SYSBP > 220)], "VS", "SYSBP", "SYSBP",
        "Systolic BP outside 70-220 mmHg. Please verify against source.")
    add("EC02", vs[(vs.DIABP < 40) | (vs.DIABP > 130)], "VS", "DIABP", "DIABP",
        "Diastolic BP outside 40-130 mmHg. Please verify against source.")
    add("EC03", vs[vs.SYSBP <= vs.DIABP], "VS", "SYSBP/DIABP", "SYSBP",
        "Systolic BP is not greater than diastolic BP. Please check for transposition.")
    add("EC04", vs[(vs.PULSE < 35) | (vs.PULSE > 150)], "VS", "PULSE", "PULSE",
        "Pulse outside 35-150 beats/min. Please verify.")
    add("EC05", vs[(vs.TEMP < 34.5) | (vs.TEMP > 40.5)], "VS", "TEMP", "TEMP",
        "Temperature outside 34.5-40.5 C. Please verify.")
    dup = vs[vs.duplicated(["SUBJID", "VISIT"], keep="first")]
    add("EC06", dup, "VS", "record", None, "Duplicate vital signs record for this visit. Please delete one.")
    svd = sv.merge(dm[["SUBJID", "RFICDAT"]], on="SUBJID")
    add("EC07", svd[pd.to_datetime(svd.SVSTDAT) < pd.to_datetime(svd.RFICDAT)], "SV", "SVSTDAT", "SVSTDAT",
        "Visit date is before the informed consent date.")
    add("EC08", ae[pd.to_datetime(ae.AEENDAT) < pd.to_datetime(ae.AESTDAT)], "AE", "AEENDAT", "AEENDAT",
        "AE end date is before the start date.")
    add("EC09", lb[lb.LBORRESU.isna() | (lb.LBORRESU.astype(str).str.strip() == "")], "LB", "LBORRESU", "LBTEST",
        "Lab unit is missing.")
    aex = ae.merge(ex[["SUBJID", "EXSTDAT"]], on="SUBJID")
    add("EC10", aex[pd.to_datetime(aex.AESTDAT) < pd.to_datetime(aex.EXSTDAT)], "AE", "AESTDAT", "AESTDAT",
        "AE started before first dose; confirm it should be recorded as an adverse event, not medical history.")
    add("EC11", ae[(ae.AESER == "Y") & (ae.AESEV == "MILD")], "AE", "AESER/AESEV", "AESEV",
        "Serious AE recorded as mild severity. Please confirm seriousness criteria.")
    ql = pd.DataFrame(q)
    ql.insert(0, "query_id", [f"Q{i:04d}" for i in range(1, len(ql) + 1)])
    ql["SITEID"] = ql.SUBJID.str[:2]
    ql["status"] = "OPEN"
    return ql


def score(ql, injected):
    """Compare the queries with the injected errors."""
    keys = {
        "EC01": {(s, v) for s, v, f in injected["vs_out_of_range"] if f == "SYSBP"},
        "EC03": {(s, v) for s, v, f in injected["vs_out_of_range"] if f == "SYSBP<DIABP"},
        "EC06": {tuple(x) for x in injected["duplicate_vs"]},
        "EC07": {tuple(x) for x in injected["date_before_consent"]},
        "EC08": {tuple(x) for x in injected["ae_end_before_start"]},
        "EC09": {tuple(x) for x in injected["lb_missing_unit"]},
    }
    rows = []
    for cid, truth in keys.items():
        sub = ql[ql.check_id == cid]
        if cid == "EC09":
            raised = set(zip(sub.SUBJID, sub.visit_or_record, sub.value))
            truth = {(s, v, t) for s, v, t in truth}
        else:
            raised = set(zip(sub.SUBJID, sub.visit_or_record))
        tp = len(raised & truth)
        rows.append({"check_id": cid, "injected_errors": len(truth), "queries_raised": len(raised),
                     "errors_caught": tp, "recall": round(tp / max(len(truth), 1), 3),
                     "queries_on_clean_data": len(raised - truth)})
    return pd.DataFrame(rows)


def main(data="data"):
    d = load_raw(data)
    ql = run_checks(d)
    ql.to_csv(os.path.join(ROOT, "results", "tables", "query_listing.csv"), index=False)
    injected = json.load(open(os.path.join(data, "injected.json")))
    sc = score(ql, injected)
    sc.to_csv(os.path.join(ROOT, "results", "tables", "edit_check_performance.csv"), index=False)
    summary = ql.groupby("check_id").size().rename("queries").reset_index()
    print(summary.to_string(index=False))
    print(sc.to_string(index=False))
    return ql


if __name__ == "__main__":
    main()
