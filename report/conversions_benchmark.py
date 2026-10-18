#!/usr/bin/env python
"""
Time the TD/FD conversions in waveformrs (``examples/conversions.rs``) and
lalsimulation (``scripts/time_conversions_lal.py``) and write
``report/data/conversion_timings.json``.

    OMP_NUM_THREADS=1 python report/conversions_benchmark.py
"""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main():
    env = dict(os.environ)
    env.setdefault("NRSUR7DQ4_DATA", str(ROOT / "lal_data" / "NRSur7dq4_v1.0.h5"))
    env.setdefault("LAL_DATA_PATH", str(ROOT / "lal_data"))
    with tempfile.TemporaryDirectory() as tmp:
        lal = Path(tmp) / "lal.json"
        out = Path(tmp) / "rust.json"
        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "time_conversions_lal.py"), "--output", str(lal)],
            cwd=ROOT, env=env, check=True, capture_output=True,
        )
        printed = subprocess.run(
            ["cargo", "run", "--quiet", "--release", "--example", "conversions", "--",
             str(lal), "--json", str(out)],
            cwd=ROOT, env=env, check=True, capture_output=True, text=True,
        ).stdout
        rows = json.loads(out.read_text())["rows"]
        lal_meta = json.loads(lal.read_text())
    target = ROOT / "report" / "data" / "conversion_timings.json"
    target.write_text(json.dumps(
        dict(lalsimulation_version=lal_meta["lalsimulation_version"], rows=rows), indent=1))
    print(printed)
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
