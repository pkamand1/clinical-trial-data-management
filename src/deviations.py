"""Data-driven protocol deviations.

Rules (from a fictional protocol):
  PD01 visit outside its window around the target study day        minor
  PD02 scheduled visit missed while the subject remained on study     minor
  PD03 study drug compliance below 80%                               major
  PD04 study procedure performed before informed consent             major

Detected deviations are compared with those injected into the simulation.
"""

import json
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WINDOWS = {"WEEK 2": (15, 3, 3), "WEEK 4": (29, 3, 3), "WEEK 8": (57, 5, 5), "WEEK 12": (85, 7, 7)}
ORDER = ["SCREENING", "DAY 1", "WEEK 2", "WEEK 4", "WEEK 8", "WEEK 12"]


def detect(raw):
    sv, ex, dm = raw["sv"], raw["ex"], raw["dm"]
    out = []
    day1 = sv[sv.VISIT == "DAY 1"].set_index("SUBJID").SVSTDAT.map(pd.Timestamp)
    for _, r in sv[sv.VISIT.isin(WINDOWS)].iterrows():
        if r.SUBJID not in day1:
            continue
        target, before, after = WINDOWS[r.VISIT]
        study_day = (pd.Timestamp(r.SVSTDAT) - day1[r.SUBJID]).days + 1
        if not (target - before <= study_day <= target + after):
            out.append({"deviation_id": "PD01", "SUBJID": r.SUBJID, "visit": r.VISIT, "severity": "minor",
                        "detail": f"study day {study_day}, window {target - before} to {target + after}"})
    completed = set(raw["ds"].loc[raw["ds"].DSDECOD == "COMPLETED", "SUBJID"])
    for subj, g in sv.groupby("SUBJID"):
        done = [ORDER.index(v) for v in g.VISIT]
        # Completers owe every visit; subjects who withdrew owe those before their last attended visit.
        last_expected = len(ORDER) if subj in completed else max(done)
        for k in range(2, last_expected):
            if k not in done:
                out.append({"deviation_id": "PD02", "SUBJID": subj, "visit": ORDER[k], "severity": "minor",
                            "detail": "visit not performed while the subject was on study"})
    ex = ex.assign(compliance=ex.DOSES_TAKEN / ex.DOSES_PLANNED)
    for _, r in ex[ex.compliance < 0.80].iterrows():
        out.append({"deviation_id": "PD03", "SUBJID": r.SUBJID, "visit": "", "severity": "major",
                    "detail": f"compliance {r.compliance:.0%}"})
    svd = sv.merge(dm[["SUBJID", "RFICDAT"]], on="SUBJID")
    for _, r in svd[pd.to_datetime(svd.SVSTDAT) < pd.to_datetime(svd.RFICDAT)].iterrows():
        out.append({"deviation_id": "PD04", "SUBJID": r.SUBJID, "visit": r.VISIT, "severity": "major",
                    "detail": f"visit {r.SVSTDAT} before consent {r.RFICDAT}"})
    dv = pd.DataFrame(out)
    dv["SITEID"] = dv.SUBJID.str[:2]
    return dv


def score(dv, inj):
    pairs = {"PD01": ({tuple(x) for x in inj["visit_out_of_window"]}, lambda d: set(zip(d.SUBJID, d.visit))),
             "PD02": ({tuple(x) for x in inj["missed_visit"]}, lambda d: set(zip(d.SUBJID, d.visit))),
             "PD03": ({x[0] for x in inj["dose_noncompliance"]}, lambda d: set(d.SUBJID)),
             "PD04": ({x[0] for x in inj["date_before_consent"]}, lambda d: set(d.SUBJID))}
    rows = []
    for k, (truth, key) in pairs.items():
        found = key(dv[dv.deviation_id == k])
        rows.append({"deviation_id": k, "injected": len(truth), "detected": len(found),
                     "injected_detected": len(found & truth), "recall": round(len(found & truth) / max(len(truth), 1), 3),
                     "detected_not_injected": len(found - truth)})
    return pd.DataFrame(rows)


def main(data="data"):
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from edit_checks import load_raw
    raw = load_raw(data)
    dv = detect(raw)
    dv.to_csv(os.path.join(ROOT, "results", "tables", "protocol_deviations.csv"), index=False)
    sc = score(dv, json.load(open(os.path.join(data, "injected.json"))))
    sc.to_csv(os.path.join(ROOT, "results", "tables", "deviation_detection.csv"), index=False)
    print(sc.to_string(index=False))


if __name__ == "__main__":
    main()
