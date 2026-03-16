#!/usr/bin/env python3
"""
Benchmark common control routines: run Rust bench_csv (up to 30s per size) and the
Fortran SLICOT benchmark, then plot Rust vs Fortran with seaborn. Run from repo root:

    python3 tools/benchmark_control/run_benchmarks.py

Requires: pandas, seaborn, matplotlib
"""

import csv
import subprocess
import sys
from pathlib import Path


def run_rust_bench(repo_root: Path, timeout: int = 600):
    """Run Rust bench_csv; return list of dicts with routine, n, time_us, impl=Rust."""
    cmd = ["cargo", "run", "--release", "--bin", "bench_csv"]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=repo_root,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr or result.stdout or "bench_csv failed")
    lines = result.stdout.strip().splitlines()
    if not lines or lines[0].strip().lower() != "routine,n,time_us":
        raise RuntimeError("bench_csv did not output expected CSV header")
    reader = csv.DictReader(lines)
    rows = [dict(r) for r in reader]
    for r in rows:
        r["impl"] = "Rust"
    return rows


def run_fortran_bench(repo_root: Path, timeout: int = 300):
    """Run Fortran benchmark script; return list of dicts with routine, n, time_us, impl=Fortran."""
    script = repo_root / "tools" / "slicot-fortran" / "run_fortran_benchmarks.sh"
    if not script.is_file():
        return []
    result = subprocess.run(
        [str(script)],
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=repo_root,
        shell=True,
    )
    if result.returncode != 0:
        return []
    lines = [s.strip() for s in result.stdout.strip().splitlines()]
    # Parse CSV data lines (routine,n,time_us); ignore build log and header
    rows = []
    for line in lines:
        parts = line.split(",", 2)
        if len(parts) != 3:
            continue
        try:
            routine, n_str, time_str = [p.strip() for p in parts]
            n = int(n_str)
            time_us = float(time_str)
            if n <= 0 or time_us < 0:
                continue
            rows.append({"routine": routine, "n": n, "time_us": time_us, "impl": "Fortran"})
        except (ValueError, TypeError):
            continue
    return rows


def main():
    repo_root = Path(__file__).resolve().parent.parent.parent

    print("Running Rust benchmarks (bench_csv)...")
    try:
        rust_rows = run_rust_bench(repo_root)
    except FileNotFoundError:
        print("error: cargo not found", file=sys.stderr)
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print("error: bench_csv timed out", file=sys.stderr)
        sys.exit(1)
    except RuntimeError as e:
        print("error:", e, file=sys.stderr)
        sys.exit(1)

    if not rust_rows:
        print("error: no Rust benchmark data", file=sys.stderr)
        sys.exit(1)

    print("Running Fortran benchmarks...")
    fortran_rows = run_fortran_bench(repo_root)
    if not fortran_rows:
        print("warning: no Fortran data (build/run may have failed); plotting Rust only.")

    try:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
    except ImportError as e:
        print("error: install pandas, seaborn, matplotlib:", e, file=sys.stderr)
        sys.exit(1)

    df_rust = pd.DataFrame(rust_rows)
    df_rust["n"] = df_rust["n"].astype(int)
    df_rust["time_us"] = df_rust["time_us"].astype(float)

    all_dfs = [df_rust]
    if fortran_rows:
        df_fortran = pd.DataFrame(fortran_rows)
        df_fortran["n"] = df_fortran["n"].astype(int)
        df_fortran["time_us"] = df_fortran["time_us"].astype(float)
        all_dfs.append(df_fortran)

    df = pd.concat(all_dfs, ignore_index=True)
    df = df[df["n"] <= 1024]

    sns.set_theme(style="darkgrid")
    fig, ax = plt.subplots(figsize=(10, 6))
    has_fortran = fortran_rows and "Fortran" in df["impl"].values
    if has_fortran:
        sns.lineplot(
            data=df,
            x="n",
            y="time_us",
            hue="routine",
            style="impl",
            markers=True,
            ax=ax,
        )
    else:
        sns.lineplot(
            data=df,
            x="n",
            y="time_us",
            hue="routine",
            marker="o",
            ax=ax,
        )
    ax.set_xlabel("Matrix size n")
    ax.set_ylabel("Time per call (µs)")
    ax.set_title("Control routines: Rust vs Fortran (up to 30s per size for Rust)")
    out = Path(__file__).parent / "benchmark_control.png"
    fig.savefig(out, dpi=150)
    print(f"Plot saved to {out}")
    plt.close()


if __name__ == "__main__":
    main()
