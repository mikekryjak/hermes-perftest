#!/usr/bin/env python3

"""
make_window.py - Create a window test from a parent run.

A window is a slice of a parent run: restarted from the parent's state at
time a and run to time b. It is named by its plasma time on a 0.1 ms grid,
e.g. `test2_2.0-3.0ms` is test 2 from 2.0 ms to 3.0 ms.

The restart files for the start time (the seed) come from a seed library at
<seeds>/<hermes_sha>/<parent case>/<a>ms/. A missing seed is cut from the parent's
dumps with sdtools' make_restart.py.

The new case gets the parent's BOUT.inp, so its physics settings match the
state the seed came from, with 50 outputs spanning the window. The seed is
copied into base/ and into the case, so reset_test.py works on it. Applying a
recipe, recording and launching are left to the usual run steps; the window
runs with -restart.

Usage:
    make_window.py test2_2.0-3.0ms <parent case> <new case> [--seeds <dir>]
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

from boutdata.collect import collect

NOUT = 50
ONE_MS = 95788  # 1 ms in normalised units, 1e-3 * e * Bnorm / Mp at Bnorm = 1 T
SEED_TOLERANCE_MS = 0.01

WINDOW_RE = re.compile(r"^(test\d+)_(\d+\.\d)-(\d+\.\d)ms$")
PARENT_RE = re.compile(r"^(test\d+)_")
HERMES_SHA_RE = re.compile(r"Git Version of Hermes:\s*([0-9a-f]{40})")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a window test from a parent run."
    )
    parser.add_argument("window", help="Window name, e.g. test2_2.0-3.0ms")
    parser.add_argument("parent", help="Parent run case directory")
    parser.add_argument("case", help="New case directory to create")
    parser.add_argument(
        "--seeds",
        help="Seed library directory. Default: $SOLVEROPT_SEEDS, else seeds/ "
        "beside the parent run",
    )
    parser.add_argument(
        "--nout",
        type=int,
        default=NOUT,
        help=f"Number of outputs across the window (default {NOUT}). A long "
        "window needs more, or its output interval outgrows a stall limit "
        "set for the parent's cadence.",
    )
    return parser.parse_args()


def parse_window(window):
    """Return (test, start, end) as (str, Decimal, Decimal)."""
    m = WINDOW_RE.match(window)
    if not m:
        raise SystemExit(
            f"Bad window name {window!r}: expected e.g. test2_2.0-3.0ms, "
            "times with one decimal on a 0.1 ms grid"
        )
    test, start, end = m.group(1), Decimal(m.group(2)), Decimal(m.group(3))
    if end <= start:
        raise SystemExit(f"Window {window} ends before it starts")
    return test, start, end


def hermes_sha(parent):
    """Return the Hermes-3 commit recorded in the parent's BOUT.log.0."""
    log = parent / "BOUT.log.0"
    if not log.exists():
        raise SystemExit(f"No BOUT.log.0 in parent {parent}")
    m = HERMES_SHA_RE.search(log.read_text(errors="replace"))
    if not m:
        raise SystemExit(f"No Hermes-3 commit found in {log}")
    return m.group(1)


def seed_time_ms(seed):
    """Return the simulation time of the restart files in seed, in ms."""
    tt = float(collect("tt", path=str(seed), prefix="BOUT.restart", info=False))
    omega_ci = float(
        collect("Omega_ci", path=str(seed), prefix="BOUT.restart", info=False)
    )
    return tt / omega_ci * 1e3


def get_seed(parent, seeds, sha, start):
    """Return the parent's seed directory at start, cutting it if missing."""
    seed = seeds / sha / parent.name / f"{start}ms"
    if not list(seed.glob("BOUT.restart.*.nc")):
        make_restart = shutil.which("make_restart.py")
        if make_restart is None:
            raise SystemExit("make_restart.py (sdtools) not found on PATH")
        subprocess.run(
            [sys.executable, make_restart, str(parent), str(seed), str(start)],
            check=True,
        )

    t_ms = seed_time_ms(seed)
    if abs(t_ms - float(start)) > SEED_TOLERANCE_MS:
        raise SystemExit(
            f"Seed {seed} is at {t_ms:.6g} ms, more than {SEED_TOLERANCE_MS} ms "
            f"from the requested {start} ms. The parent has no output there."
        )
    return seed, t_ms


def write_input(parent, case, start, end, nout=NOUT):
    """Copy the parent's BOUT.inp, setting nout and timestep for the window."""
    step = (end - start) / nout
    lines = (parent / "BOUT.inp").read_text().splitlines(keepends=True)
    done = set()
    for i, line in enumerate(lines):
        if line.lstrip().startswith("["):
            break  # nout and timestep are top-level options
        key = line.split("=")[0].strip()
        if key == "nout":
            lines[i] = f"nout = {nout}   # window {start}-{end} ms\n"
            done.add(key)
        elif key == "timestep":
            lines[i] = f"timestep = {ONE_MS} * {step}   # {ONE_MS} = 1ms in normalised units\n"
            done.add(key)
    if done != {"nout", "timestep"}:
        raise SystemExit(f"Parent BOUT.inp lacks top-level nout or timestep")
    (case / "BOUT.inp").write_text("".join(lines))
    return step


def main():
    args = parse_args()
    test, start, end = parse_window(args.window)

    parent = Path(args.parent).expanduser().resolve()
    case = Path(args.case).expanduser().resolve()
    if not parent.is_dir():
        raise SystemExit(f"Parent case not found: {parent}")
    m = PARENT_RE.match(parent.name)
    if not m or m.group(1) != test:
        raise SystemExit(f"Parent {parent.name} is not a {test} run")
    if case.exists():
        raise SystemExit(f"Case already exists: {case}")

    seeds = Path(
        args.seeds or os.environ.get("SOLVEROPT_SEEDS") or parent.parent / "seeds"
    ).expanduser().resolve()

    sha = hermes_sha(parent)
    seed, t_ms = get_seed(parent, seeds, sha, start)

    case.mkdir(parents=True)
    step = write_input(parent, case, start, end, args.nout)
    (case / "base").mkdir()
    for f in sorted(seed.glob("BOUT.restart.*.nc")):
        shutil.copy2(f, case / "base" / f.name)
        shutil.copy2(f, case / f.name)

    print(f"Created {case.name}")
    print(f"  parent  {parent.name} (Hermes-3 {sha[:12]})")
    print(f"  seed    {seed} at {t_ms:.6g} ms")
    print(f"  output  nout = {args.nout}, timestep = {ONE_MS} * {step} ({start}-{end} ms)")
    print("  next    apply a recipe, open the row, launch with -restart")


if __name__ == "__main__":
    main()
