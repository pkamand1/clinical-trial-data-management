"""Conformance checks on the SDTM datasets, in the spirit of the validation
rules sponsors run before submission (a small, transparent subset)."""

import os
import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REQUIRED = {"DM": ["STUDYID", "DOMAIN", "USUBJID", "SUBJID", "RFSTDTC", "SITEID", "SEX", "ARMCD", "ARM", "COUNTRY"],
            "SV": ["STUDYID", "DOMAIN", "USUBJID", "VISITNUM"],
            "VS": ["STUDYID", "DOMAIN", "USUBJID", "VSSEQ", "VSTESTCD", "VSTEST"],
            "LB": ["STUDYID", "DOMAIN", "USUBJID", "LBSEQ", "LBTESTCD", "LBTEST"],
            "AE": ["STUDYID", "DOMAIN", "USUBJID", "AESEQ", "AETERM", "AEDECOD"],
            "EX": ["STUDYID", "DOMAIN", "USUBJID", "EXSEQ", "EXTRT"],
            "DS": ["STUDYID", "DOMAIN", "USUBJID", "DSSEQ", "DSTERM", "DSDECOD"]}
CT = {("DM", "SEX"): {"M", "F", "U", "UNDIFFERENTIATED"}, ("AE", "AESEV"): {"MILD", "MODERATE", "SEVERE"},
      ("AE", "AESER"): {"Y", "N"}, ("VS", "VSTESTCD"): {"SYSBP", "DIABP", "PULSE", "TEMP"},
      ("LB", "LBTESTCD"): {"CREAT", "ALT", "HGB"}, ("LB", "LBNRIND"): {"LOW", "NORMAL", "HIGH", ""}}
ISO = re.compile(r"^\d{4}(-\d{2}(-\d{2}(T\d{2}:\d{2}(:\d{2})?)?)?)?$")


def run(domains):
    res = []

    def add(rule, dom, desc, n):
        res.append({"rule": rule, "domain": dom, "description": desc, "records_failing": int(n),
                    "status": "pass" if n == 0 else "issue"})

    dm_ids = set(domains["DM"].USUBJID)
    for dom, df in domains.items():
        missing = [v for v in REQUIRED[dom] if v not in df.columns or df[v].isna().all()]
        add("SD01", dom, "required variables present and populated", len(missing))
        add("SD02", dom, "every USUBJID exists in DM", (~df.USUBJID.isin(dm_ids)).sum())
        seqc = f"{dom}SEQ"
        if seqc in df:
            add("SD03", dom, "--SEQ unique within USUBJID", df.duplicated(["USUBJID", seqc]).sum())
        for c in [c for c in df.columns if c.endswith("DY")]:
            add("SD04", dom, f"{c} is never 0", (df[c] == 0).sum())
        for c in [c for c in df.columns if c.endswith("DTC")]:
            add("SD05", dom, f"{c} is ISO 8601", (~df[c].dropna().astype(str).str.match(ISO)).sum())
        for (d2, var), allowed in CT.items():
            if d2 == dom:
                add("SD06", dom, f"{var} uses controlled terminology", (~df[var].fillna("").isin(allowed)).sum())
    lb = domains["LB"]
    add("SD07", "LB", "one standard unit per test", (lb.groupby("LBTESTCD").LBSTRESU.nunique() > 1).sum())
    add("SD08", "LB", "standard result present when original result present",
        (lb.LBORRES.ne("") & lb.LBSTRESN.isna()).sum())
    add("SD09", "AE", "AEDECOD coded for every AETERM", domains["AE"].AEDECOD.isna().sum())
    ae = domains["AE"]
    add("SD10", "AE", "AEENDTC not before AESTDTC", (pd.to_datetime(ae.AEENDTC) < pd.to_datetime(ae.AESTDTC)).sum())
    vs = domains["VS"]
    add("SD11", "VS", "one result per subject, visit, and test", vs.duplicated(["USUBJID", "VISITNUM", "VSTESTCD"]).sum())
    return pd.DataFrame(res)


def main(domains):
    r = run(domains)
    r.to_csv(os.path.join(ROOT, "results", "tables", "sdtm_conformance.csv"), index=False)
    print(r[r.status == "issue"].to_string(index=False))
    return r
