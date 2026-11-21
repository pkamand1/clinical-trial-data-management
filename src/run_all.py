"""Run the whole data management workflow from the raw EDC export.

Usage: python src/run_all.py [--data data]
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deviations  # noqa: E402
import edit_checks  # noqa: E402
import rbqm  # noqa: E402
import sdtm  # noqa: E402
import sdtm_checks  # noqa: E402
import simulate_edc  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(ROOT, "data"))
    a = ap.parse_args()
    for d in ("tables", "figures", "sdtm"):
        os.makedirs(os.path.join(ROOT, "results", d), exist_ok=True)
    simulate_edc.main(a.data)
    queries = edit_checks.main(a.data)
    deviations.main(a.data)
    raw = edit_checks.load_raw(a.data)
    dv = pd.read_csv(os.path.join(ROOT, "results", "tables", "protocol_deviations.csv"), dtype={"SITEID": str})
    domains = sdtm.main(a.data)
    conf = sdtm_checks.main(domains)
    rbqm.main(raw, queries, dv)

    # Database lock readiness: what must be resolved, by site.
    site = raw["dm"].set_index("SUBJID").SITEID
    ready = pd.DataFrame({
        "subjects": raw["dm"].SITEID.value_counts(),
        "open_queries": queries.SUBJID.map(site).value_counts(),
        "major_deviations_to_review": dv[dv.severity == "major"].SUBJID.map(site).value_counts(),
        "minor_deviations": dv[dv.severity == "minor"].SUBJID.map(site).value_counts(),
    }).fillna(0).astype(int).sort_index()
    ready.index.name = "SITEID"
    ready.loc["Total"] = ready.sum()
    ready.to_csv(os.path.join(ROOT, "results", "tables", "lock_readiness.csv"))
    print(ready.to_string())
    print("SDTM conformance issues:", int(conf[conf.status == "issue"].records_failing.sum()),
          "records, all traceable to open queries")


if __name__ == "__main__":
    main()
