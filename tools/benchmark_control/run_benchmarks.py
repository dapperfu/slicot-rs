#!/usr/bin/env python3
"""
Benchmark common control routines: run the Rust bench_csv binary (up to 30s per size),
then plot results with seaborn. Run from the repository root:

    python3 tools/benchmark_control/run_benchmarks.py

Requires: pandas, seaborn, matplotlib
"""

import csv
import subprocess
import sys
from pathlib import Path

def main():
    repo_root = Path(__file__).resolve().parent.parent.parent
    cmd = ["cargo", "run", "--release", "--bin", "bench_csv"]
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,
            cwd=repo_root,
        )
    except FileNotFoundError:
        print("error: cargo not found", file=sys.stderr)
        sys.exit(1)
    except subprocess.TimeoutError:
        print("error: bench_csv timed out (max 600s total)", file=sys.stderr)
        sys.exit(1)

    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    lines = result.stdout.strip().splitlines()
    if not lines:
        print("error: no CSV output", file=sys.stderr)
        sys.exit(1)

    reader = csv.DictReader(lines)
    rows = list(reader)
    if not rows:
        print("error: no data rows", file=sys.stderr)
        sys.exit(1)

    try:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
    except ImportError as e:
        print("error: install pandas, seaborn, matplotlib:", e, file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(rows)
    df["n"] = df["n"].astype(int)
    df["time_us"] = df["time_us"].astype(float)

    sns.set_theme(style="darkgrid")
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.lineplot(data=df, x="n", y="time_us", hue="routine", marker="o", ax=ax)
    ax.set_xlabel("Matrix size n")
    ax.set_ylabel("Time per call (µs)")
    ax.set_title("Control routines benchmark (up to 30s per size)")
    out = Path(__file__).parent / "benchmark_control.png"
    fig.savefig(out, dpi=150)
    print(f"Plot saved to {out}")
    plt.close()


if __name__ == "__main__":
    main()
