#!/usr/bin/env python3
"""Semi-analytic tables from the ChiPT UFO (package ../chipt), for what Sherpa
cannot do: the running Delta width and the rho/omega exchange.

    python3 tables.py photo      photoproduction variants        -> output/photo_curves.csv
    python3 tables.py cc         CC running-width factor w(E_nu) -> output/running_width_cc.csv
    python3 tables.py nc         NC running-width factor w(E_nu) -> output/running_width_nc.csv
    python3 tables.py electron   e N -> e' N pi running-width responses
                                                                 -> output/running_width_electron.csv

Runs in parallel on half the cores (--jobs to change).  cc and nc save each
point as it finishes and skip saved points when rerun.

photo curves (sigma in ub vs E_gamma):
    base         this work (as Sherpa, to ~1%)
    no_ff        base with the piNDelta vertex form factor set to 1
    vec          this work + rho/omega exchange
    bw           vec with the Breit-Wigner Delta parameters and the running width
    delta_const  Delta diagrams only, constant width (gamma p -> p pi0)
    delta_k3     Delta diagrams only, running width at the pole parameters
    bw_formula   unitary Breit-Wigner cross section at the pole parameters

cc, nc: w(E_nu) = sigma(Breit-Wigner, running width) / sigma(pole, constant
width); CC with and without W < 1.4 GeV, NC without a cut for neutrinos and
antineutrinos (channel name + "~").  Multiplying a Sherpa constant-width
prediction by w (event by event, at the event's E_nu) gives the running-width
prediction.

electron: the virtual-photon responses S_T, S_L of both models on a (W, Q^2)
grid; the event weight is (S_T^bw + eps S_L^bw) / (S_T^pole + eps S_L^pole).
"""
import argparse
import csv
import os
import sys
from multiprocessing import Pool

from dataclasses import replace

# One BLAS thread per worker (the matrices are 4x4), so --jobs is the core count.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from chipt import model, photo, weak  # noqa: E402

OUT = os.path.join(HERE, "output")
_models = {}


def get(delta, running=None, no_ff=False):
    key = (delta, running, no_ff)
    if key not in _models:
        m = model.load(delta, running)
        if no_ff:
            m = replace(m, ff=dict(m.ff, FDeltaLam=lambda p2: 1.0))
        _models[key] = m
    return _models[key]


# ---------------------------------------------------------------------- photo
def photo_point(args):
    ch, eg = args
    rows = [(eg, ch, "base", photo.sigma(ch, eg, get("pole"))),
            (eg, ch, "no_ff", photo.sigma(ch, eg, get("pole", no_ff=True))),
            (eg, ch, "vec", photo.sigma(ch, eg, get("pole"), vec=True)),
            (eg, ch, "bw", photo.sigma(ch, eg, get("bw"), vec=True))]
    if ch == "ppi0":
        rows += [(eg, ch, "delta_const", photo.sigma(ch, eg, get("pole"), born=False)),
                 (eg, ch, "delta_k3", photo.sigma(ch, eg, get("pole", True), born=False)),
                 (eg, ch, "bw_formula", photo.breit_wigner(eg, get("pole")))]
    return rows


def photo_table(pool):
    egs = np.round(np.arange(0.150, 0.6001, 0.005), 4)
    jobs = [(ch, eg) for ch in photo.CHANNELS for eg in egs]
    rows = [r for rs in pool.map(photo_point, jobs) for r in rs]
    write("photo_curves.csv", ["Egamma", "channel", "curve", "sigma_ub"], rows)


# ---------------------------------------------------------------------- width
ENU = (0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0)


def width_point(args):
    ch, E, wmax, nubar = args
    # A finer grid without the W cut above 1 GeV, where the W range is wide.
    grid = (16, 10, 10, 8) if (wmax is None and E >= 1.0) else (10, 6, 6, 5)
    pole = weak.sigma(ch, E, get("pole"), wmax=wmax, grid=grid, nubar=nubar)
    bw = weak.sigma(ch, E, get("bw"), wmax=wmax, grid=grid, nubar=nubar)
    return (ch + ("~" if nubar else ""), "none" if wmax is None else "%g" % wmax,
            "%.4f" % E, "%.6f" % (bw / pole))


def width_table(pool, name, jobs):
    """Each point is appended to the csv as it finishes and points already in it
    are skipped, so an interrupted table can be resumed."""
    header = ["channel", "wmax", "E_nu", "w"]
    path = os.path.join(OUT, name)
    key = lambda ch, E, wmax, nubar: (ch + ("~" if nubar else ""),
                                      "none" if wmax is None else "%g" % wmax, "%.4f" % E)
    done = {}
    if os.path.exists(path):
        done = {tuple(r[:3]): r for r in list(csv.reader(open(path)))[1:]}
    todo = [j for j in jobs if key(*j) not in done]
    # The fine-grid points are the slow ones; start them first.
    todo.sort(key=lambda j: (j[2] is not None, -j[1]))
    os.makedirs(OUT, exist_ok=True)
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(header)
        for row in pool.imap_unordered(width_point, todo, chunksize=1):
            done[tuple(row[:3])] = row
            w.writerow(row)
            f.flush()
    write(name, header, [done[key(*j)] for j in jobs])


def cc_table(pool):
    jobs = [(ch, E, wm, False) for wm in (1.4, None) for ch in weak.CC_CHANNELS for E in ENU]
    width_table(pool, "running_width_cc.csv", jobs)


def nc_table(pool):
    jobs = [(ch, E, None, nb) for nb in (False, True) for ch in weak.NC_CHANNELS for E in ENU]
    width_table(pool, "running_width_nc.csv", jobs)


# ------------------------------------------------------------------- electron
W_GRID = np.round(np.arange(1.08, 1.8001, 0.01), 6)
Q2_GRID = (0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.55, 0.7, 0.85, 1.0)


def electron_point(args):
    ch, W, Q2 = args
    tp, lp = photo.responses(ch, W, Q2, get("pole"))
    tb, lb = photo.responses(ch, W, Q2, get("bw"))
    return ch, W, Q2, tp, lp, tb, lb


def electron_table(pool):
    jobs = [(ch, W, Q2) for ch in photo.CHANNELS for Q2 in Q2_GRID for W in W_GRID]
    write("running_width_electron.csv",
          ["channel", "W", "Q2", "ST_pole", "SL_pole", "ST_bw", "SL_bw"],
          pool.map(electron_point, jobs, chunksize=8))


def write(name, header, rows):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print("wrote", path)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", choices=["photo", "cc", "nc", "electron"])
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    a = ap.parse_args()
    with Pool(a.jobs) as pool:
        {"photo": photo_table, "cc": cc_table, "nc": nc_table,
         "electron": electron_table}[a.table](pool)
