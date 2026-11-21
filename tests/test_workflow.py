"""Score the workflow against the errors and site behaviors injected into the data."""

import json
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
inj = json.load(open(os.path.join(ROOT, "data", "injected.json")))


def test_edit_checks_catch_every_injected_error_without_false_queries():
    sc = pd.read_csv(os.path.join(TAB, "edit_check_performance.csv"))
    assert (sc.recall == 1.0).all()
    assert (sc.queries_on_clean_data == 0).all()


def test_injected_deviations_detected():
    sc = pd.read_csv(os.path.join(TAB, "deviation_detection.csv"))
    assert (sc.recall == 1.0).all()


def test_red_flags_are_exactly_the_planted_site_problems():
    st = pd.read_csv(os.path.join(TAB, "site_risk_status.csv"), dtype={"Unnamed: 0": str}).set_index("Unnamed: 0")
    red = {(s, col) for col in st.columns for s in st.index if st.loc[s, col] == "red"}
    sites = inj["sites"]
    assert red == {(sites["ae_underreporting"], "AE reporting"), (sites["slow_data_entry"], "Data entry lag"),
                   (sites["low_variability_digit_preference"], "Digit preference"),
                   (sites["low_variability_digit_preference"], "Low variability")}
    assert st.loc[sites["creatinine_in_umol_per_l"], "Lab units"] == "amber"


def test_sdtm_conformance_issues_trace_to_open_queries():
    conf = pd.read_csv(os.path.join(TAB, "sdtm_conformance.csv")).set_index("rule")
    issues = conf[conf.status == "issue"]
    assert set(issues.index) == {"SD08", "SD10", "SD11"}
    assert issues.loc["SD08", "records_failing"] == len(inj["lb_missing_unit"])
    assert issues.loc["SD10", "records_failing"] == len(inj["ae_end_before_start"])
    assert issues.loc["SD11", "records_failing"] == 4 * len(inj["duplicate_vs"])


def test_sdtm_study_day_and_units():
    lb = pd.read_csv(os.path.join(ROOT, "results", "sdtm", "lb.csv"))
    vs = pd.read_csv(os.path.join(ROOT, "results", "sdtm", "vs.csv"))
    assert (vs.VSDY != 0).all() and (lb.LBDY != 0).all()
    creat = lb[(lb.LBTESTCD == "CREAT") & lb.LBSTRESN.notna()]
    assert creat.LBSTRESU.eq("mg/dL").all() and creat.LBSTRESN.between(0.3, 3).all()
