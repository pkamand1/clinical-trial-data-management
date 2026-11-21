"""Simulate raw EDC exports for a fictional phase 2 trial.

Study XYZ-201: randomized, double-blind, placebo-controlled, 240 subjects at
8 sites, visits at screening, day 1, and weeks 2, 4, 8, and 12. Raw forms use
CDASH-style variable names: DM, IE, SV (visit dates), VS, LB, AE, EX, DS.

Errors and site behaviors are injected deliberately and recorded in
data/injected.json, so the edit checks, deviation rules, and monitoring can be
scored against what is actually wrong:

  * out-of-range and transposed values, impossible dates, an AE that ends
    before it starts, duplicate records, visits before informed consent
  * visits outside their protocol windows and missed visits
  * one site reporting creatinine in umol/L instead of mg/dL
  * Site 07: vital signs with too little variability and every blood pressure
    ending in 0 (a classic sign of data not measured as recorded)
  * Site 03: adverse events reported at about a third of the rate of other sites
"""

import json
import os
import sys
from datetime import date, timedelta

import numpy as np
import pandas as pd

SEED = 201
STUDY = "XYZ-201"
SITES = ["01", "02", "03", "04", "05", "06", "07", "08"]
VISITS = [("SCREENING", -14, 0, 7), ("DAY 1", 1, 0, 0), ("WEEK 2", 15, 3, 3), ("WEEK 4", 29, 3, 3),
          ("WEEK 8", 57, 5, 5), ("WEEK 12", 85, 7, 7)]          # name, target study day, window before, after
AE_TERMS = ["headache", "Headache", "nausea", "Nausea ", "fatigue", "tiredness", "dizziness", "rash",
            "diarrhoea", "diarrhea", "upper respiratory infection", "URI", "back pain", "insomnia",
            "elevated ALT", "ALT increased"]


def main(out="data"):
    rng = np.random.default_rng(SEED)
    raw = os.path.join(out, "raw")
    os.makedirs(raw, exist_ok=True)
    inj = {k: [] for k in ["vs_out_of_range", "date_before_consent", "ae_end_before_start", "duplicate_vs",
                           "visit_out_of_window", "missed_visit", "lb_missing_unit", "dose_noncompliance"]}
    dm, ie, sv, vs, lb, ae, ex, ds = ([] for _ in range(8))
    subj_n = 0
    for site in SITES:
        n_site = 30
        for _ in range(n_site):
            subj_n += 1
            subj = f"{site}-{subj_n:03d}"
            consent = date(2025, 2, 1) + timedelta(days=int(rng.integers(0, 150)))
            screen = consent + timedelta(days=int(rng.integers(0, 4)))
            arm = rng.choice(["XYZ 10 mg", "Placebo"])
            age = int(rng.integers(25, 76))
            dm.append({"SUBJID": subj, "SITEID": site, "BRTHYR": consent.year - age, "SEX": rng.choice(["M", "F"]),
                       "RACE": rng.choice(["WHITE", "BLACK OR AFRICAN AMERICAN", "ASIAN", "OTHER"], p=[.6, .2, .15, .05]),
                       "ETHNIC": rng.choice(["HISPANIC OR LATINO", "NOT HISPANIC OR LATINO"], p=[.18, .82]),
                       "RFICDAT": consent.isoformat(), "ARM": arm})
            ie.append({"SUBJID": subj, "IEYN": "Y", "IEDAT": screen.isoformat()})
            day1 = screen + timedelta(days=int(rng.integers(7, 15)))
            withdrew_after = rng.choice(len(VISITS), p=[0, 0, .02, .03, .04, .91])   # last visit index completed
            base_sbp = rng.normal(132, 14)
            base_dbp, base_hr = 16 + 0.5 * base_sbp + rng.normal(0, 6), rng.normal(72, 9)
            effect = -6 if arm == "XYZ 10 mg" else 0
            base_creat, base_alt = rng.lognormal(np.log(0.9), 0.2), rng.lognormal(np.log(24), 0.35)
            last_visit_date = day1
            for vi, (vname, tday, wb, wa) in enumerate(VISITS):
                if vi > withdrew_after:
                    break
                if vname == "SCREENING":
                    vdate = screen
                else:
                    if vi >= 2 and rng.random() < 0.025:                     # missed visit
                        inj["missed_visit"].append([subj, vname])
                        continue
                    shift = int(rng.integers(-wb, wa + 1)) if wb or wa else 0
                    if vi >= 2 and rng.random() < 0.04:                     # out of window
                        shift = int(rng.choice([-1, 1]) * (wa + rng.integers(2, 6)))
                        inj["visit_out_of_window"].append([subj, vname])
                    vdate = day1 + timedelta(days=tday - 1 + shift)
                if vi == 1 and rng.random() < 0.01:                          # visit dated before consent
                    vdate = consent - timedelta(days=3)
                    inj["date_before_consent"].append([subj, vname])
                entry_lag = rng.gamma(2, 2.5 if site != "05" else 9)        # site 05 enters data slowly
                sv.append({"SUBJID": subj, "VISIT": vname, "SVSTDAT": vdate.isoformat(),
                           "ENTRYDAT": (vdate + timedelta(days=int(entry_lag))).isoformat()})
                last_visit_date = vdate
                weeks = max(0, (tday - 1) / 7)
                trt = effect * min(1, weeks / 4)
                if site == "07":                                            # low variability, rounded
                    sbp = 10 * round((base_sbp + trt + rng.normal(0, 1.5)) / 10)
                    dbp = 10 * round((base_dbp + rng.normal(0, 1.0)) / 10)
                    hr = round(base_hr + rng.normal(0, 1.0))
                else:
                    sbp = round(base_sbp + trt + rng.normal(0, 8))
                    dbp = round(base_dbp + trt / 2 + rng.normal(0, 6))
                    hr = round(base_hr + rng.normal(0, 6))
                temp = round(rng.normal(36.7, 0.3), 1)
                rec = {"SUBJID": subj, "VISIT": vname, "VSDAT": vdate.isoformat(), "SYSBP": sbp, "DIABP": dbp,
                       "PULSE": hr, "TEMP": temp, "SYSBP_U": "mmHg", "DIABP_U": "mmHg", "PULSE_U": "beats/min",
                       "TEMP_U": "C"}
                if rng.random() < 0.012:                                    # typing error
                    rec["SYSBP"] = rec["SYSBP"] * 10
                    inj["vs_out_of_range"].append([subj, vname, "SYSBP"])
                elif rng.random() < 0.008:                                  # systolic and diastolic transposed
                    rec["SYSBP"], rec["DIABP"] = rec["DIABP"], rec["SYSBP"]
                    inj["vs_out_of_range"].append([subj, vname, "SYSBP<DIABP"])
                vs.append(rec)
                if rng.random() < 0.01:
                    vs.append(dict(rec))
                    inj["duplicate_vs"].append([subj, vname])
                if vname in ("SCREENING", "WEEK 4", "WEEK 12"):
                    creat = base_creat * rng.lognormal(0, 0.06)
                    unit = "mg/dL"
                    if site == "04":
                        creat, unit = creat * 88.42, "umol/L"
                    for test, val, u in [("CREAT", round(creat, 2 if unit == "mg/dL" else 0), unit),
                                         ("ALT", round(base_alt * rng.lognormal(0, 0.15)), "U/L"),
                                         ("HGB", round(rng.normal(13.6, 1.3), 1), "g/dL")]:
                        if rng.random() < 0.006:
                            u = ""
                            inj["lb_missing_unit"].append([subj, vname, test])
                        lb.append({"SUBJID": subj, "VISIT": vname, "LBDAT": vdate.isoformat(), "LBTEST": test,
                                   "LBORRES": val, "LBORRESU": u})
            # Exposure: daily dosing between day 1 and the last visit, with compliance
            last = last_visit_date
            planned = max(1, (last - day1).days)
            compliance = min(1.0, rng.beta(18, 1.2))
            if rng.random() < 0.06:
                compliance = rng.uniform(0.55, 0.78)
                inj["dose_noncompliance"].append([subj])
            ex.append({"SUBJID": subj, "EXTRT": arm.replace("XYZ 10 mg", "XYZ").replace("Placebo", "PLACEBO"),
                       "EXDOSE": 10 if arm != "Placebo" else 0, "EXDOSU": "mg", "EXSTDAT": day1.isoformat(),
                       "EXENDAT": last.isoformat(), "DOSES_PLANNED": planned,
                       "DOSES_TAKEN": int(round(planned * compliance))})
            # Adverse events
            ae_rate = 0.9 * (1.3 if arm != "Placebo" else 1.0) * (0.33 if site == "03" else 1.0)
            for _ in range(rng.poisson(ae_rate * planned / 84)):
                st = day1 + timedelta(days=int(rng.integers(0, planned + 1)))
                en = st + timedelta(days=int(rng.integers(1, 15)))
                rec = {"SUBJID": subj, "AETERM": rng.choice(AE_TERMS), "AESTDAT": st.isoformat(),
                       "AEENDAT": en.isoformat(), "AESEV": rng.choice(["MILD", "MODERATE", "SEVERE"], p=[.6, .32, .08]),
                       "AESER": "Y" if rng.random() < 0.03 else "N",
                       "AEREL": rng.choice(["RELATED", "NOT RELATED"], p=[.35, .65])}
                if rng.random() < 0.015:
                    rec["AEENDAT"] = (st - timedelta(days=3)).isoformat()
                    inj["ae_end_before_start"].append([subj, rec["AETERM"]])
                ae.append(rec)
            completed = withdrew_after == len(VISITS) - 1
            ds.append({"SUBJID": subj, "DSDECOD": "COMPLETED" if completed else
                       rng.choice(["ADVERSE EVENT", "WITHDRAWAL BY SUBJECT", "LOST TO FOLLOW-UP"]),
                       "DSSTDAT": last.isoformat()})

    for name, rows in [("dm", dm), ("ie", ie), ("sv", sv), ("vs", vs), ("lb", lb), ("ae", ae), ("ex", ex), ("ds", ds)]:
        pd.DataFrame(rows).to_csv(os.path.join(raw, f"{name}.csv"), index=False)
    inj["sites"] = {"low_variability_digit_preference": "07", "ae_underreporting": "03",
                    "creatinine_in_umol_per_l": "04", "slow_data_entry": "05"}
    json.dump(inj, open(os.path.join(out, "injected.json"), "w"), indent=1)
    print(f"{len(dm)} subjects, {len(sv)} visits, {len(vs)} VS rows, {len(lb)} LB rows, {len(ae)} AEs")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
