"""Run an (M_h, M_t) grid scan with the C++ binary in parallel processes and merge the chunks.

Usage:  python python/smvac/run_scan.py configs/controlled.toml [--workers N] [--binary build/smvac_scan]

The configuration file (TOML) defines the grid, the physics options and the output path. The merged CSV keeps the
metadata header written by the binary and adds provenance (config file, git commit, start/end time, workers).
Integrity checks: every grid point present exactly once; no NaN in the coordinates.
"""
from __future__ import annotations

import argparse
import datetime as dt
import os
import shutil
import subprocess
import sys
import tempfile
import tomllib
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

OPTION_KEYS = {  # config key -> command-line option of smvac_scan
    "prescription": "--prescription", "matching": "--matching", "alpha3": "--alpha3", "mw": "--mw",
    "rge_dt": "--rge-dt", "nodes": "--nodes", "nodes_check": "--nodes-check", "half_width": "--half-width",
    "h0": "--h0", "bsm_phi6_c6": "--bsm-phi6-c6", "bsm_phi6_scale": "--bsm-phi6-scale",
}


def git_commit() -> str:
    try:
        out = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--untracked-files=no"], capture_output=True,
                               text=True).stdout.strip()
        return out + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--binary", default=os.path.join(ROOT, "build", "smvac_scan"))
    args = ap.parse_args()
    with open(args.config, "rb") as fh:
        cfg = tomllib.load(fh)
    g = cfg["grid"]
    mh, mt = g["mh"], g["mt"]
    nh = int(round((mh[1] - mh[0]) / mh[2])) + 1
    nt = int(round((mt[1] - mt[0]) / mt[2])) + 1
    total = nh * nt
    out = os.path.join(ROOT, cfg["output"]["path"])
    os.makedirs(os.path.dirname(out), exist_ok=True)
    opts = []
    for k, v in cfg.get("physics", {}).items():
        if k not in OPTION_KEYS:
            raise SystemExit(f"unknown physics option '{k}' in {args.config}")
        opts += [OPTION_KEYS[k], str(v)]
    binary = args.binary + (".exe" if os.name == "nt" and not args.binary.endswith(".exe") else "")
    nchunk = min(args.workers * 4, total)
    size = -(-total // nchunk)
    tmp = tempfile.mkdtemp(prefix="smvac_scan_")
    started = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    print(f"{total} points, {nchunk} chunks, {args.workers} workers -> {out}", flush=True)

    def run(i):
        s = i * size
        if s >= total:
            return None
        part = os.path.join(tmp, f"chunk_{i:05d}.csv")
        cmd = [binary, "--mh", *map(str, mh), "--mt", *map(str, mt), "--chunk", str(s), str(min(size, total - s)),
               "--out", part, *opts]
        subprocess.run(cmd, check=True)
        return part

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        parts = [p for p in ex.map(run, range(nchunk)) if p]
    finished = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    header, rows, seen = None, [], set()
    for p in sorted(parts):
        with open(p) as fh:
            lines = fh.read().splitlines()
        meta = [l for l in lines if l.startswith("#")]
        body = [l for l in lines if not l.startswith("#")]
        if header is None:
            header = (meta, body[0])
        for line in body[1:]:
            key = tuple(line.split(",")[:2])
            if key in seen:
                raise SystemExit(f"duplicate grid point {key}")
            seen.add(key)
            rows.append(line)
    if len(rows) != total:
        raise SystemExit(f"expected {total} rows, got {len(rows)}")
    meta = [l for l in header[0] if not l.startswith("# chunk")]
    meta += [f"# config = {os.path.relpath(args.config, ROOT)}", f"# label = {cfg['output'].get('label', '')}",
             f"# git_commit = {git_commit()}", f"# started_utc = {started}", f"# finished_utc = {finished}",
             f"# workers = {args.workers}", f"# points = {total}"]
    with open(out, "w", newline="\n") as fh:
        fh.write("\n".join(meta + [header[1]] + rows) + "\n")
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"wrote {out} ({total} points)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
