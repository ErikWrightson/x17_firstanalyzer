#!/usr/bin/env python3
"""
@author Rafayel Paremuzyan <rafopar@jlab.org>
Adapted by Erik Wrightson for use in this project to run batch jobs.

Run_AnalysisBatch.py

Usage: python3 Run_AnalysisBatch.py <RUN>

Workflow:
  1. Scan /volatile/hallb/prad/x17_replay/3cl_full_data/ for prad_<RUN>_recon_<FileNo>.root files.
  2. Make sure outfiles/<RUN>/ exists.
  3. For every file, launch  ../x17_firstanalysis <RUN> <FileNo>
     in parallel, capped at MAX_PARALLEL=18 concurrent jobs.
  4. When all jobs finish, hadd outfiles/<RUN>/x17_Ana_<RUN>_*.root
     into outfiles/<RUN>/x17_Ana_<RUN>_All.root.
  5. Delete the individual per-file root files.

Filename convention:
    prad_024014_recon_015.root   (RUN=24014, FileNo=15)
The 6-digit RUN and 5-digit FileNo widths are preserved when forming output
filenames (e.g. x17_Ana_024014_00015.root).
"""

import argparse
import glob
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# ── Configuration ──────────────────────────────────────────────────────────
MAX_PARALLEL = 18
EXECUTABLE   = "/work/hallb/prad/wrightso/x17_firstanalyzer/x17_firstanalyzer"
DATA_DIR     = "/volatile/hallb/prad/x17_replay/3cl_full_data/"
OUT_BASE     = "/work/hallb/prad/wrightso/x17_firstanalyzer/outfiles"
RUN_WIDTH    = 6        # zero-pad width for run numbers
FILE_WIDTH   = 3        # zero-pad width for file numbers
# ─────────────────────────────────────────────────────────────────────────


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def fmt_seconds(s):
    if s < 60:
        return f"{s:.1f}s"
    m, s = divmod(s, 60)
    if m < 60:
        return f"{int(m)}m{int(s):02d}s"
    h, m = divmod(m, 60)
    return f"{int(h)}h{int(m):02d}m{int(s):02d}s"


def print_table(title, headers, rows):
    col_widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            col_widths[i] = max(col_widths[i], len(str(cell)))
    sep = "+" + "+".join("-" * (w + 2) for w in col_widths) + "+"
    def fmt(cells):
        return "|" + "|".join(
            f" {str(c):<{col_widths[i]}} " for i, c in enumerate(cells)
        ) + "|"
    print(title)
    print(sep)
    print(fmt(headers))
    print(sep)
    for row in rows:
        print(fmt(row))
    print(sep)


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------

def run_analysis(run, file_no, loc):
    """
    Launch  ../x17_firstanalyzer <RUN> <FileNo>  and wait.
    Returns (file_no, returncode, elapsed_seconds, stderr_tail).
    """

    t0 = time.time()
    try:
        file_name = OUT_BASE +"/" + str(run) + "/" + f"x17_Ana_0{run}_" + str(file_no)
        #print(file_name)
        command = [EXECUTABLE, "-f", str(loc), "-b", str(run), "-o", str(file_name)]#, str(run), str(file_no)]
        
        result = subprocess.run(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
        rc  = result.returncode
        err = result.stderr or ""
    except FileNotFoundError:
        return file_no, -1, 0.0, f"executable not found: {EXECUTABLE}"
    except OSError as exc:
        return file_no, -1, 0.0, f"launch failed: {exc}"

    elapsed = time.time() - t0
    tail = "\n".join(err.strip().splitlines()[-5:])
    return file_no, rc, elapsed, tail


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Run x17_firstanalyzer in parallel for all recon files of a run, then hadd."
    )
    parser.add_argument("run", type=int, help="Run number (e.g. 24014)")
    parser.add_argument("--max-parallel", type=int, default=MAX_PARALLEL,
                        help=f"Max concurrent jobs (default: {MAX_PARALLEL}).")
    args = parser.parse_args()

    run         = args.run
    run_str     = f"{run:0{RUN_WIDTH}d}"          # e.g. "024014"
    n_workers   = args.max_parallel

    sp = os.path.join(DATA_DIR, "*"+str(run)+"*")
    raw_dir     = glob.glob(sp)[0]
    out_dir     = os.path.join(OUT_BASE, str(run))

    # --- Preflight checks --------------------------------------------------
    if not os.path.isdir(raw_dir):
        print(f"Error: directory not found: {raw_dir}")
        sys.exit(1)

    if not os.path.isfile(EXECUTABLE) or not os.access(EXECUTABLE, os.X_OK):
        print(f"Error: '{EXECUTABLE}' is not an executable file (cwd: {os.getcwd()}).")
        sys.exit(1)

    os.makedirs(out_dir, exist_ok=True)
    print(f"Output directory ready: {out_dir}")

    # --- Discover input files ----------------------------------------------
    pattern = os.path.join(raw_dir, f"prad_{run_str}_recon_[0-9][0-9][0-9].root")
    files   = sorted(glob.glob(pattern))

    if not files:
        print(f"No files matching {pattern}")
        sys.exit(1)

    file_pat = re.compile(rf"^prad_{run_str}_recon_(\d+)\.root$")
    file_nos = []
    for f in files:
        m = file_pat.match(os.path.basename(f))
        if m:
            file_nos.append(int(m.group(1)))
        else:
            print(f"  WARNING: skipping unexpected filename: {os.path.basename(f)}")

    if not file_nos:
        print("No matching file numbers extracted; nothing to do.")
        sys.exit(1)

    file_nos = sorted(set(file_nos))
    n_total  = len(file_nos)

    print(f"Found {n_total} input file(s) for run {run_str}")
    print(f"Launching '{EXECUTABLE}' with up to {n_workers} concurrent job(s)")
    print("=" * 70)

    # ------------------------------------------------------------------ #
    # 1. Run x17_firstanalysis in parallel                                         #
    # ------------------------------------------------------------------ #
    results  = {}     # {file_no: (rc, elapsed, stderr_tail)}
    t_start  = time.time()
    width    = len(str(n_total))
    done     = 0
    
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        fut_to_fn = {pool.submit(run_analysis, run, fn, loc): (fn, loc)
                     for fn, loc in zip(file_nos, files)}
        for fut in as_completed(fut_to_fn):
            fn, rc, elapsed, tail = fut.result()
            results[fn] = (rc, elapsed, tail)
            done += 1
            status = "OK" if rc == 0 else f"FAIL(rc={rc})"
            print(f"[{done:>{width}}/{n_total}]  File {fn:0{FILE_WIDTH}d}  "
                  f"{status:>10}   {fmt_seconds(elapsed)}")

    total_elapsed = time.time() - t_start

    n_ok   = sum(1 for rc, _, _ in results.values() if rc == 0)
    n_fail = n_total - n_ok

    print()
    print(f"All {n_total} job(s) finished in {fmt_seconds(total_elapsed)}  "
          f"(OK: {n_ok}, FAIL: {n_fail})")

    # ------------------------------------------------------------------ #
    # 2. hadd outputs                                                      #
    # ------------------------------------------------------------------ #
    out_pattern    = os.path.join(out_dir, f"x17_Ana_{run_str}_[0-9]*.root")
    individual_outputs = sorted(glob.glob(out_pattern))
    # Exclude the merged file from the inputs in case of a re-run
    merged_path    = os.path.join(out_dir, f"x17_Ana_{run_str}_All.root")
    individual_outputs = [f for f in individual_outputs if f != merged_path]

    if not individual_outputs:
        print("\nNo per-file root outputs found to merge. Skipping hadd.")
    else:
        print()
        print("=" * 70)
        print(f"Merging {len(individual_outputs)} root file(s) -> "
              f"{os.path.basename(merged_path)}")
        print("=" * 70)

        cmd = ["hadd", "-f", merged_path] + individual_outputs
        # Print a short version so we don't spam the terminal with hundreds of paths
        print(f"  hadd -f {merged_path}  [{len(individual_outputs)} input files]")

        t0 = time.time()
        hadd_result = subprocess.run(cmd)
        t_hadd = time.time() - t0

        if hadd_result.returncode == 0:
            print(f"\nhadd finished in {fmt_seconds(t_hadd)} -> {merged_path}")
            # ---------------------------------------------------------- #
            # 3. Delete individual files                                   #
            # ---------------------------------------------------------- #
            print(f"\nRemoving {len(individual_outputs)} individual root file(s) ...")
            n_removed = 0
            for f in individual_outputs:
                try:
                    os.remove(f)
                    n_removed += 1
                except OSError as exc:
                    print(f"  WARNING: could not remove {f}: {exc}")
            print(f"Removed {n_removed} file(s).")
        else:
            print(f"\nERROR: hadd exited with code {hadd_result.returncode}.")
            print("Individual root files have NOT been deleted.")

    # ------------------------------------------------------------------ #
    # 4. Final summary                                                     #
    # ------------------------------------------------------------------ #
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Run                 : {run_str}")
    print(f"  Input files         : {n_total}")
    print(f"  Jobs successful     : {n_ok}")
    print(f"  Jobs failed         : {n_fail}")
    print(f"  Wall-clock time     : {fmt_seconds(total_elapsed)}")
    if individual_outputs:
        print(f"  Merged output       : {merged_path}")

    failed = [(fn, rc, el, tail)
              for fn, (rc, el, tail) in results.items() if rc != 0]
    if failed:
        failed.sort(key=lambda x: x[0])
        print()
        print_table(
            "Failed jobs:",
            ["FileNo", "Exit code", "Elapsed"],
            [(f"{fn:0{FILE_WIDTH}d}", rc, fmt_seconds(el))
             for fn, rc, el, _ in failed],
        )
        print()
        print("Stderr tail for failed jobs:")
        for fn, rc, _, tail in failed:
            print(f"  File {fn:0{FILE_WIDTH}d}:")
            if tail:
                for line in tail.splitlines():
                    print(f"    {line}")
            else:
                print("    (no stderr output)")


if __name__ == "__main__":
    main()
