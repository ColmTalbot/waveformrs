#!/usr/bin/env python
"""
Time NRSur7dq4 in waveformrs (``examples/nrsur7dq4.rs``) and lalsimulation
(``scripts/time_nrsur7dq4_lal.py``) on the same configurations and write
``report/data/nrsur_timings.json``.

    OMP_NUM_THREADS=1 python report/nrsur_benchmark.py

Needs ``lal_data/NRSur7dq4_v1.0.h5`` (or ``NRSUR7DQ4_DATA`` / ``LAL_DATA_PATH``).
"""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main():
    data = ROOT / "lal_data"
    env = dict(os.environ)
    env.setdefault("NRSUR7DQ4_DATA", str(data / "NRSur7dq4_v1.0.h5"))
    env.setdefault("LAL_DATA_PATH", str(data))
    rust = json.loads(subprocess.run(
        ["cargo", "run", "--quiet", "--release", "--example", "nrsur7dq4", "--", "--json"],
        cwd=ROOT, env=env, capture_output=True, text=True, check=True,
    ).stdout)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "lal.json"
        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "time_nrsur7dq4_lal.py"), "--json", str(path)],
            cwd=ROOT, env=env, check=True, capture_output=True,
        )
        lal_times = json.loads(path.read_text())
    out = dict(rust=rust, lal=lal_times)
    target = ROOT / "report" / "data" / "nrsur_timings.json"
    target.write_text(json.dumps(out, indent=1))
    for r, l in zip(rust["configurations"], lal_times["configurations"]):
        print(f"M={r['total_mass']:>5} rate={r['sampling_rate']:>6} f_min={r['f_min']:>4}: "
              f"Rust {r['fresh'] * 1e3:.2f} ms (workspace {r['reused'] * 1e3:.2f}, extrinsic "
              f"{r['extrinsic'] * 1e3:.3f}), LAL {l['polarizations'] * 1e3:.2f} ms")
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
