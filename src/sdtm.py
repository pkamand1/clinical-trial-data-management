"""Map the raw EDC data to CDISC SDTM domains (DM, SV, VS, LB, AE, EX, DS).

Follows SDTM Implementation Guide conventions: USUBJID built from STUDYID and
SUBJID, ISO 8601 dates, study days relative to RFSTDTC with no day 0,
--SEQ within subject, original and standardized results (--ORRES/--STRESN with
units converted to one standard unit per test), controlled terminology for
test codes, VISITNUM and EPOCH, and normal range indicators for labs.

Outputs SAS transport files (.xpt, the format submitted to the FDA) and CSV
copies in results/sdtm/, and a variable-level specification in
specs/sdtm_spec.csv generated from the same metadata the code uses.

Values are carried as collected. Open queries are resolved at the source in
EDC, not edited here, so records with open queries are mapped as they stand.
AEDECOD uses a small illustrative dictionary in place of MedDRA, which is
licensed.
"""

import os

import numpy as np
import pandas as pd
import pyreadstat

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STUDY = "XYZ-201"
VISITNUM = {"SCREENING": 1, "DAY 1": 2, "WEEK 2": 3, "WEEK 4": 4, "WEEK 8": 5, "WEEK 12": 6}
VS_TESTS = {"SYSBP": ("Systolic Blood Pressure", "mmHg"), "DIABP": ("Diastolic Blood Pressure", "mmHg"),
            "PULSE": ("Pulse Rate", "beats/min"), "TEMP": ("Temperature", "C")}
LB_TESTS = {"CREAT": ("Creatinine", "mg/dL", 0.6, 1.3), "ALT": ("Alanine Aminotransferase", "U/L", 7, 56),
            "HGB": ("Hemoglobin", "g/dL", 12.0, 17.5)}
UNIT_CONVERSION = {("CREAT", "umol/L"): 1 / 88.42}       # to mg/dL
AE_DICTIONARY = {"headache": "Headache", "nausea": "Nausea", "fatigue": "Fatigue", "tiredness": "Fatigue",
                 "dizziness": "Dizziness", "rash": "Rash", "diarrhoea": "Diarrhoea", "diarrhea": "Diarrhoea",
                 "upper respiratory infection": "Upper respiratory tract infection",
                 "uri": "Upper respiratory tract infection", "back pain": "Back pain", "insomnia": "Insomnia",
                 "elevated alt": "Alanine aminotransferase increased",
                 "alt increased": "Alanine aminotransferase increased"}

LABELS = {
    "STUDYID": "Study Identifier", "DOMAIN": "Domain Abbreviation", "USUBJID": "Unique Subject Identifier",
    "SUBJID": "Subject Identifier for the Study", "SITEID": "Study Site Identifier",
    "RFSTDTC": "Subject Reference Start Date/Time", "RFENDTC": "Subject Reference End Date/Time",
    "RFICDTC": "Date/Time of Informed Consent", "AGE": "Age", "AGEU": "Age Units", "SEX": "Sex", "RACE": "Race",
    "ETHNIC": "Ethnicity", "ARMCD": "Planned Arm Code", "ARM": "Description of Planned Arm",
    "ACTARMCD": "Actual Arm Code", "ACTARM": "Description of Actual Arm", "COUNTRY": "Country",
    "VISITNUM": "Visit Number", "VISIT": "Visit Name", "EPOCH": "Epoch", "SVSTDTC": "Start Date/Time of Visit",
    "SVSTDY": "Study Day of Start of Visit",
    "VSSEQ": "Sequence Number", "VSTESTCD": "Vital Signs Test Short Name", "VSTEST": "Vital Signs Test Name",
    "VSORRES": "Result or Finding in Original Units", "VSORRESU": "Original Units",
    "VSSTRESC": "Character Result/Finding in Std Format", "VSSTRESN": "Numeric Result/Finding in Standard Units",
    "VSSTRESU": "Standard Units", "VSDTC": "Date/Time of Measurements", "VSDY": "Study Day of Vital Signs",
    "LBSEQ": "Sequence Number", "LBTESTCD": "Lab Test or Examination Short Name", "LBTEST": "Lab Test or Examination Name",
    "LBORRES": "Result or Finding in Original Units", "LBORRESU": "Original Units",
    "LBSTRESC": "Character Result/Finding in Std Format", "LBSTRESN": "Numeric Result/Finding in Standard Units",
    "LBSTRESU": "Standard Units", "LBSTNRLO": "Normal Range Lower Limit-Std Units",
    "LBSTNRHI": "Normal Range Upper Limit-Std Units", "LBNRIND": "Reference Range Indicator",
    "LBDTC": "Date/Time of Specimen Collection", "LBDY": "Study Day of Specimen Collection",
    "AESEQ": "Sequence Number", "AETERM": "Reported Term for the Adverse Event",
    "AEDECOD": "Dictionary-Derived Term", "AESEV": "Severity/Intensity", "AESER": "Serious Event",
    "AEREL": "Causality", "AESTDTC": "Start Date/Time of Adverse Event", "AEENDTC": "End Date/Time of Adverse Event",
    "AESTDY": "Study Day of Start of Adverse Event", "AEENDY": "Study Day of End of Adverse Event",
    "EXSEQ": "Sequence Number", "EXTRT": "Name of Treatment", "EXDOSE": "Dose", "EXDOSU": "Dose Units",
    "EXDOSFRQ": "Dosing Frequency per Interval", "EXROUTE": "Route of Administration",
    "EXSTDTC": "Start Date/Time of Treatment", "EXENDTC": "End Date/Time of Treatment",
    "EXSTDY": "Study Day of Start of Treatment", "EXENDY": "Study Day of End of Treatment",
    "DSSEQ": "Sequence Number", "DSTERM": "Reported Term for the Disposition Event",
    "DSDECOD": "Standardized Disposition Term", "DSCAT": "Category for Disposition Event",
    "DSSTDTC": "Start Date/Time of Disposition Event", "DSSTDY": "Study Day of Start of Disposition Event",
}
DERIVATION = {
    "USUBJID": "STUDYID || '-' || SUBJID", "RFSTDTC": "EX.EXSTDAT (first dose)", "RFENDTC": "EX.EXENDAT (last dose)",
    "AGE": "year of RFICDTC minus DM.BRTHYR", "ARMCD": "XYZ10 or PBO from DM.ARM",
    "EPOCH": "SCREENING before RFSTDTC, otherwise TREATMENT", "--DY": "date − RFSTDTC + 1 on or after RFSTDTC, date − RFSTDTC before (no day 0)",
    "--SEQ": "sequence within USUBJID, ordered by date and test", "VSSTRESN": "numeric VSORRES (units already standard)",
    "LBSTRESN": "LBORRES converted to standard units (creatinine umol/L ÷ 88.42 → mg/dL)",
    "LBNRIND": "LOW, NORMAL, or HIGH against LBSTNRLO and LBSTNRHI", "AEDECOD": "illustrative dictionary lookup of AETERM",
}


def study_day(dtc, ref):
    d = (pd.to_datetime(dtc) - pd.to_datetime(ref)).dt.days
    return np.where(d >= 0, d + 1, d)


def seq(df, by):
    return df.groupby("USUBJID").cumcount() + 1 if by is None else \
        df.sort_values(["USUBJID"] + by).groupby("USUBJID").cumcount() + 1


def build(raw):
    dmr, ex = raw["dm"], raw["ex"]
    first = ex.set_index("SUBJID")
    dm = pd.DataFrame({
        "STUDYID": STUDY, "DOMAIN": "DM", "USUBJID": STUDY + "-" + dmr.SUBJID, "SUBJID": dmr.SUBJID,
        "RFSTDTC": dmr.SUBJID.map(first.EXSTDAT), "RFENDTC": dmr.SUBJID.map(first.EXENDAT),
        "RFICDTC": dmr.RFICDAT, "SITEID": dmr.SITEID,
        "AGE": pd.to_datetime(dmr.RFICDAT).dt.year - dmr.BRTHYR, "AGEU": "YEARS", "SEX": dmr.SEX, "RACE": dmr.RACE,
        "ETHNIC": dmr.ETHNIC, "ARMCD": dmr.ARM.map({"XYZ 10 mg": "XYZ10", "Placebo": "PBO"}), "ARM": dmr.ARM,
        "ACTARMCD": dmr.ARM.map({"XYZ 10 mg": "XYZ10", "Placebo": "PBO"}), "ACTARM": dmr.ARM, "COUNTRY": "USA"})
    ref = dm.set_index("SUBJID").RFSTDTC
    usub = lambda s: STUDY + "-" + s
    epoch = lambda dtc, subj: np.where(pd.to_datetime(dtc) < pd.to_datetime(subj.map(ref)), "SCREENING", "TREATMENT")

    svr = raw["sv"]
    sv = pd.DataFrame({"STUDYID": STUDY, "DOMAIN": "SV", "USUBJID": usub(svr.SUBJID), "VISITNUM": svr.VISIT.map(VISITNUM),
                       "VISIT": svr.VISIT, "EPOCH": epoch(svr.SVSTDAT, svr.SUBJID), "SVSTDTC": svr.SVSTDAT,
                       "SVSTDY": study_day(svr.SVSTDAT, svr.SUBJID.map(ref))})

    vsr = raw["vs"].melt(id_vars=["SUBJID", "VISIT", "VSDAT"], value_vars=list(VS_TESTS), var_name="VSTESTCD",
                         value_name="VSORRES")
    vs = pd.DataFrame({"STUDYID": STUDY, "DOMAIN": "VS", "USUBJID": usub(vsr.SUBJID), "VSTESTCD": vsr.VSTESTCD,
                       "VSTEST": vsr.VSTESTCD.map(lambda t: VS_TESTS[t][0]), "VSORRES": vsr.VSORRES.map(lambda v: f"{v:g}"),
                       "VSORRESU": vsr.VSTESTCD.map(lambda t: VS_TESTS[t][1]), "VSSTRESC": vsr.VSORRES.map(lambda v: f"{v:g}"),
                       "VSSTRESN": vsr.VSORRES.astype(float), "VSSTRESU": vsr.VSTESTCD.map(lambda t: VS_TESTS[t][1]),
                       "VISITNUM": vsr.VISIT.map(VISITNUM), "VISIT": vsr.VISIT, "EPOCH": epoch(vsr.VSDAT, vsr.SUBJID),
                       "VSDTC": vsr.VSDAT, "VSDY": study_day(vsr.VSDAT, vsr.SUBJID.map(ref))})

    lbr = raw["lb"].copy()
    conv = [UNIT_CONVERSION.get((t, u), 1.0) for t, u in zip(lbr.LBTEST, lbr.LBORRESU.fillna(""))]
    stdn = (lbr.LBORRES.astype(float) * conv).round(2)
    lo = lbr.LBTEST.map(lambda t: LB_TESTS[t][2])
    hi = lbr.LBTEST.map(lambda t: LB_TESTS[t][3])
    lb = pd.DataFrame({"STUDYID": STUDY, "DOMAIN": "LB", "USUBJID": usub(lbr.SUBJID), "LBTESTCD": lbr.LBTEST,
                       "LBTEST": lbr.LBTEST.map(lambda t: LB_TESTS[t][0]), "LBORRES": lbr.LBORRES.map(lambda v: f"{v:g}"),
                       "LBORRESU": lbr.LBORRESU.fillna(""),
                       # A result without a unit cannot be standardized: leave the standard result empty.
                       "LBSTRESC": np.where(lbr.LBORRESU.isna(), "", stdn.map(lambda v: f"{v:g}")),
                       "LBSTRESN": np.where(lbr.LBORRESU.isna(), np.nan, stdn),
                       "LBSTRESU": lbr.LBTEST.map(lambda t: LB_TESTS[t][1]), "LBSTNRLO": lo, "LBSTNRHI": hi,
                       "VISITNUM": lbr.VISIT.map(VISITNUM), "VISIT": lbr.VISIT, "EPOCH": epoch(lbr.LBDAT, lbr.SUBJID),
                       "LBDTC": lbr.LBDAT, "LBDY": study_day(lbr.LBDAT, lbr.SUBJID.map(ref))})
    lb["LBNRIND"] = np.select([lb.LBSTRESN < lb.LBSTNRLO, lb.LBSTRESN > lb.LBSTNRHI, lb.LBSTRESN.notna()],
                              ["LOW", "HIGH", "NORMAL"], "")

    aer = raw["ae"]
    ae = pd.DataFrame({"STUDYID": STUDY, "DOMAIN": "AE", "USUBJID": usub(aer.SUBJID), "AETERM": aer.AETERM.str.strip(),
                       "AEDECOD": aer.AETERM.str.strip().str.lower().map(AE_DICTIONARY), "AESEV": aer.AESEV,
                       "AESER": aer.AESER, "AEREL": aer.AEREL, "EPOCH": "TREATMENT", "AESTDTC": aer.AESTDAT,
                       "AEENDTC": aer.AEENDAT, "AESTDY": study_day(aer.AESTDAT, aer.SUBJID.map(ref)),
                       "AEENDY": study_day(aer.AEENDAT, aer.SUBJID.map(ref))})

    ex_ = pd.DataFrame({"STUDYID": STUDY, "DOMAIN": "EX", "USUBJID": usub(ex.SUBJID), "EXTRT": ex.EXTRT,
                        "EXDOSE": ex.EXDOSE, "EXDOSU": ex.EXDOSU, "EXDOSFRQ": "QD", "EXROUTE": "ORAL",
                        "EPOCH": "TREATMENT", "EXSTDTC": ex.EXSTDAT, "EXENDTC": ex.EXENDAT,
                        "EXSTDY": study_day(ex.EXSTDAT, ex.SUBJID.map(ref)), "EXENDY": study_day(ex.EXENDAT, ex.SUBJID.map(ref))})
    dsr = raw["ds"]
    ds = pd.DataFrame({"STUDYID": STUDY, "DOMAIN": "DS", "USUBJID": usub(dsr.SUBJID), "DSTERM": dsr.DSDECOD,
                       "DSDECOD": dsr.DSDECOD, "DSCAT": "DISPOSITION EVENT", "EPOCH": "TREATMENT",
                       "DSSTDTC": dsr.DSSTDAT, "DSSTDY": study_day(dsr.DSSTDAT, dsr.SUBJID.map(ref))})

    out = {"DM": dm.drop(columns=[]), "SV": sv, "VS": vs, "LB": lb, "AE": ae, "EX": ex_, "DS": ds}
    sort_keys = {"VS": ["VSDTC", "VSTESTCD"], "LB": ["LBDTC", "LBTESTCD"], "AE": ["AESTDTC", "AETERM"],
                 "EX": ["EXSTDTC"], "DS": ["DSSTDTC"]}
    for dom, keys in sort_keys.items():
        df = out[dom].sort_values(["USUBJID"] + keys).reset_index(drop=True)
        df.insert(3, f"{dom}SEQ", df.groupby("USUBJID").cumcount() + 1)
        out[dom] = df
    out["SV"] = out["SV"].sort_values(["USUBJID", "VISITNUM"]).reset_index(drop=True)
    return out


def write(domains):
    d = os.path.join(ROOT, "results", "sdtm")
    os.makedirs(d, exist_ok=True)
    spec = []
    for dom, df in domains.items():
        df.to_csv(os.path.join(d, f"{dom.lower()}.csv"), index=False)
        labels = [LABELS.get(c, c) for c in df.columns]
        x = df.copy()
        for c in x.columns:
            if x[c].dtype == object:
                x[c] = x[c].fillna("").astype(str)
        pyreadstat.write_xport(x, os.path.join(d, f"{dom.lower()}.xpt"), file_format_version=5,
                               table_name=dom, column_labels=labels)
        for c in df.columns:
            key = "--SEQ" if c.endswith("SEQ") else ("--DY" if c.endswith("DY") and c != "DOMAIN" else c)
            spec.append({"domain": dom, "variable": c, "label": LABELS.get(c, ""),
                         "type": "Num" if pd.api.types.is_numeric_dtype(df[c]) else "Char",
                         "origin": "Derived" if key in DERIVATION else ("Assigned" if c in ("STUDYID", "DOMAIN", "EPOCH")
                                                                       else "Collected"),
                         "derivation": DERIVATION.get(key, "")})
    pd.DataFrame(spec).to_csv(os.path.join(ROOT, "specs", "sdtm_spec.csv"), index=False)


def main(data="data"):
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from edit_checks import load_raw
    domains = build(load_raw(data))
    write(domains)
    for k, v in domains.items():
        print(f"{k}: {len(v):,} records, {v.shape[1]} variables")
    return domains


if __name__ == "__main__":
    main()
