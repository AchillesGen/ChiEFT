#!/usr/bin/env python3
"""Free-nucleon cross sections from Sherpa with the ChiEFT UFO (constant Delta width).

    python3 sherpa_scan.py photo     gamma N -> pi N, 4 channels      -> output/sherpa_photo.csv
    python3 sherpa_scan.py cc        nu_mu N -> mu- N' pi, W < 1.4 GeV -> output/sherpa_cc.csv

Needs the ChiEFT_UFO Sherpa model built and installed (see the top-level
README).  The Sherpa executable is taken from $SHERPA, default `Sherpa`.
Each point runs in its own directory under output/sherpa_runs/.  Points
already in the csv are skipped, so an interrupted scan can be resumed.
"""
import argparse
import csv
import glob
import os
import shutil
import string
import subprocess
import sys
import zipfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from chieft import PARAM_CARD  # noqa: E402

OUT = os.path.join(HERE, "output")
SHERPA = os.environ.get("SHERPA", "Sherpa")
MN = 0.938272
# Comix does not integrate gamma N with the nucleon exactly at rest, so the
# photon cards give it E = 0.938273 GeV (pz = 1.4 MeV, head-on).  The photon
# energy recorded is the equivalent one on a nucleon at rest.
E_NUCLEON_PHOTO = 0.938273
REST_FACTOR = (E_NUCLEON_PHOTO + np.sqrt(E_NUCLEON_PHOTO ** 2 - MN ** 2)) / MN
MPI0, MPIC = 0.1349768, 0.1395704

PHOTO = {  # channel: (target, process, pion mass)
    "ppi0": (2212, "2212 22 -> 2212 111", MPI0),
    "npip": (2212, "2212 22 -> 2112 211", MPIC),
    "npi0": (2112, "2112 22 -> 2112 111", MPI0),
    "ppim": (2112, "2112 22 -> 2212 -211", MPIC),
}
CC = {  # channel: (target, process, pion mass)
    "ppip": (2212, "14 2212 -> 13 2212 211", MPIC),
    "ppi0": (2112, "14 2112 -> 13 2212 111", MPI0),
    "npip": (2112, "14 2112 -> 13 2112 211", MPIC),
}
WCUT = 1.4
SEEDS = (1, 2)        # CC points are the average of two seeds
MMU = 0.10566


def photon_grid():
    thr = ((MN + MPI0) ** 2 - MN ** 2) / (2 * MN)
    return np.unique(np.concatenate([
        np.arange(thr + 0.002, 0.20, 0.010),
        np.arange(0.200, 0.500, 0.005),
        np.arange(0.500, 0.800, 0.020),
        np.arange(0.800, 1.501, 0.050)]).round(6))


def neutrino_grid():
    """Contains the ANL/BNL bin centres; dense where sigma rises steeply."""
    pts = {0.4, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5}
    pts |= {round(0.300 + 0.025 * i, 4) for i in range(15)}
    pts |= {round(0.650 + 0.050 * i, 4) for i in range(8)}
    pts |= {round(1.000 + 0.100 * i, 4) for i in range(7)}
    return np.array(sorted(pts))


def run_sherpa(card, rundir):
    """Total cross section in ub.  Charged-pion channels can crash Comix on the
    first (cold) build; the retry then reuses the process map it left behind."""
    os.makedirs(rundir, exist_ok=True)
    with open(os.path.join(rundir, "Run.yaml"), "w") as f:
        f.write(card)
    for attempt in range(2):
        with open(os.path.join(rundir, "sherpa.log"), "w") as log:
            rc = subprocess.run([SHERPA, "Run.yaml"], cwd=rundir, stdout=log,
                                stderr=subprocess.STDOUT).returncode
        if rc == 0:
            break
    else:
        raise RuntimeError(f"Sherpa failed in {rundir} (exit {rc}); see sherpa.log")
    res = os.path.join(rundir, "Results")
    with zipfile.ZipFile(os.path.join(rundir, "Results.zip")) as z:
        z.extractall(res)
    stats = glob.glob(os.path.join(res, "*", "*", "Statistics.dat"))
    last = open(stats[0]).readlines()[-1].split()
    shutil.rmtree(res)
    os.remove(os.path.join(rundir, "Results.zip"))
    return float(last[1]) * 1e-6, float(last[2]) * 1e-6           # pb -> ub


def scan(kind):
    tmpl = string.Template(open(os.path.join(HERE, "sherpa", f"{kind}.yaml")).read())
    path = os.path.join(OUT, f"sherpa_{kind}.csv")
    done = {}
    if os.path.exists(path):
        done = {(r["channel"], r["E"]): r for r in csv.DictReader(open(path))}
    channels, grid = (PHOTO, photon_grid()) if kind == "photo" else (CC, neutrino_grid())
    for ch, (target, process, mpi) in channels.items():
        lepton = 0.0 if kind == "photo" else MMU
        thr = ((MN + mpi + lepton) ** 2 - MN ** 2) / (2 * MN)
        for E in grid[grid * (REST_FACTOR if kind == "photo" else 1.0) > thr + 1e-4]:
            key = (ch, "%.6f" % E)
            if key in done or any(r["channel"] == ch and r.get("E_beam") == key[1]
                                  for r in done.values()):
                continue
            sub = dict(param_card=PARAM_CARD, nucleon=target, energy="%.6f" % E,
                       mass=MN if kind == "cc" else E_NUCLEON_PHOTO,
                       process=process, wcut=WCUT)
            if kind == "cc":
                h = [int(x) for x in process.split("->")[1].split() if int(x) != 13]
                runs = [run_sherpa(tmpl.substitute(sub, seed=s, hadron1=h[0], hadron2=h[1]),
                                   os.path.join(OUT, "sherpa_runs", f"cc_{ch}_{E:.4f}_s{s}"))
                        for s in SEEDS]
                xs = np.mean([r[0] for r in runs])
                err = np.std([r[0] for r in runs], ddof=1) / np.sqrt(len(runs))
            else:
                xs, err = run_sherpa(tmpl.substitute(sub),
                                     os.path.join(OUT, "sherpa_runs", f"photo_{ch}_{E:.4f}"))
                key = (ch, "%.6f" % (E * REST_FACTOR))
            done[key] = dict(channel=ch, E=key[1], E_beam="%.6f" % E,
                             sigma_ub="%.8g" % xs, err_ub="%.3g" % err)
            print(f"{kind} {ch} E = {E:.4f} GeV: {xs:.6g} ub", flush=True)
            os.makedirs(OUT, exist_ok=True)
            with open(path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=["channel", "E", "E_beam", "sigma_ub", "err_ub"])
                w.writeheader()
                w.writerows(done.values())
    print("wrote", path)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["photo", "cc"])
    scan(ap.parse_args().kind)
