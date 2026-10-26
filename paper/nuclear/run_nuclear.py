#!/usr/bin/env python3
"""Achilles + NUISANCE predictions for the nuclear and electron-scattering samples.

    python3 run_nuclear.py generate [--experiment NAME ...] [--events N] [--jobs J]
        Achilles events, one file per experiment of analyses.yml
        -> work/events/<experiment>.hepmc.gz.  Needs $ACHILLES_BUILD (an Achilles
        build with bin/achilles and data/) and the ChiPT_UFO Sherpa model.

    python3 run_nuclear.py analyze [--experiment NAME ...] [--docker]
        Bins the events with NUISANCE, nominal and with the running-width weights,
        -> results/<sample>.json.  Needs pyNUISANCE; --docker runs this stage in
        the achilles-physval:chipt image (docker/Dockerfile).
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
WORK = os.path.join(HERE, "work")
RESULTS = os.path.join(HERE, "results")
IMAGE = "achilles-physval:chipt"


class _Include(str):
    pass


class _Loader(yaml.SafeLoader):
    pass


class _Dumper(yaml.SafeDumper):
    pass


_Loader.add_constructor("!include", lambda l, n: _Include(l.construct_scalar(n)))
_Dumper.add_representer(_Include, lambda d, v: d.represent_scalar("!include", str(v)))


def event_file(exp):
    return os.path.join(WORK, "events", exp["name"] + ".hepmc.gz")


# ------------------------------------------------------------------- generate
def generate(exp, events, seed):
    build = os.environ["ACHILLES_BUILD"]
    rundir = os.path.join(WORK, "runs", exp["name"])
    os.makedirs(rundir, exist_ok=True)
    for name, target in (("data", os.path.join(build, "data")),
                         ("flux", os.path.join(HERE, "flux")),
                         ("FormFactors.yml", os.path.join(HERE, "FormFactors.yml"))):
        if not os.path.lexists(os.path.join(rundir, name)):
            os.symlink(target, os.path.join(rundir, name))
    shutil.copy(os.path.join(REPO, "param_card.dat"), rundir)

    card = yaml.load(open(os.path.join(HERE, exp["card"])), Loader=_Loader)
    card["Main"]["NEvents"] = events
    fifo = os.path.join(rundir, "events.hepmc")
    card["Main"]["Output"] = {"Format": "NuHepMC", "Name": fifo, "Zipped": False}
    card["Options"] = yaml.load(open(os.path.join(rundir, card["Options"])), Loader=_Loader)
    card["Options"].setdefault("Initialize", {})["Seed"] = seed
    with open(os.path.join(rundir, "run.yml"), "w") as f:
        yaml.dump(card, f, Dumper=_Dumper, sort_keys=False)

    # Achilles writes plain text into a FIFO that gzip compresses on the fly.
    out = event_file(exp)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if os.path.lexists(fifo):
        os.remove(fifo)
    os.mkfifo(fifo)
    rfd = os.open(fifo, os.O_RDONLY | os.O_NONBLOCK)
    hold = os.open(fifo, os.O_WRONLY)         # gzip sees EOF only once this closes
    os.set_blocking(rfd, True)
    with open(out + ".part", "wb") as gz_out:
        gz = subprocess.Popen(["gzip", "-1", "-c"], stdin=rfd, stdout=gz_out)
    os.close(rfd)
    print(f"[run ] {exp['name']}: {events} events in {rundir}", flush=True)
    with open(os.path.join(rundir, "achilles.log"), "w") as log:
        subprocess.run([os.path.join(build, "bin", "achilles"), "run.yml"], cwd=rundir,
                       stdout=log, stderr=subprocess.STDOUT)
    os.close(hold)
    gz.wait()
    os.remove(fifo)
    ok = "Event Run Concluded - Success" in open(os.path.join(rundir, "achilles.log")).read()
    if ok:
        os.replace(out + ".part", out)
    print(f"[{'done' if ok else 'FAIL'}] {exp['name']}", flush=True)
    return ok


# -------------------------------------------------------------------- analyze
def binned(sample, w, scale, smear):
    vals = np.bincount(sample.bins, weights=w, minlength=sample.nbins)[:sample.nbins] * scale
    sumw2 = np.bincount(sample.bins, weights=w ** 2, minlength=sample.nbins)[:sample.nbins] \
        * scale ** 2
    if smear is not None:
        vals, sumw2 = smear @ vals, np.clip(smear @ sumw2, 0.0, None)
    return {"values": vals.tolist(), "stat_err": np.sqrt(sumw2).tolist()}


def analyze(cfg, exp):
    import nuisance
    import running_width as rw
    path = event_file(exp)
    if not os.path.exists(path):
        print(f"[miss] {path}")
        return
    tables = {k: os.path.join(HERE, v) for k, v in cfg["running_width"].items()}
    weights = rw.Weights(tables["cc"], tables["nc"], tables["electron"])
    divisor = float(exp["target_nucleons"]) if exp["data_per"] == "nucleon" else 1.0
    samples = nuisance.histogram(path, exp["samples"], divisor, rw.add_columns, rw.COLUMNS)
    os.makedirs(RESULTS, exist_ok=True)
    for name in exp["samples"]:
        s = samples[name]
        values, cov, edges = nuisance.data_table(name)
        f = float(cfg.get("data_scale", {}).get(name, 1.0))
        values, cov = values * f, cov * f ** 2
        scale = np.ones(s.nbins)
        if name in cfg.get("bin_edges", {}):
            edges = np.asarray(cfg["bin_edges"][name], float)
            scale /= np.diff(edges)
        if name in cfg.get("solid_angle", {}):
            scale /= float(cfg["solid_angle"][name])
        smear = None
        if name in cfg.get("smearing", {}):
            smear = np.loadtxt(os.path.join(HERE, cfg["smearing"][name]), delimiter=",", ndmin=2)
        res = {"measurement": name, "experiment": exp["name"],
               "target_label": exp["target_label"],
               "edges": None if edges is None else list(edges),
               "data": values.tolist(), "data_cov": cov.tolist(),
               "predictions": {
                   "nominal": binned(s, s.weights, scale, smear),
                   "running_width": binned(s, s.weights * weights(s.extra), scale, smear)}}
        with open(os.path.join(RESULTS, name + ".json"), "w") as f:
            json.dump(res, f, indent=1)
        print(f"[bin ] {name}: {s.bins.size} events", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["generate", "analyze"])
    ap.add_argument("--experiment", action="append", help="default: all")
    ap.add_argument("--events", type=int, default=500000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--jobs", type=int, default=1,
                    help="parallel Achilles runs (each needs several GB)")
    ap.add_argument("--docker", action="store_true", help="analyze inside " + IMAGE)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(os.path.join(HERE, "analyses.yml")))
    exps = [e for e in cfg["experiments"] if not a.experiment or e["name"] in a.experiment]

    if a.stage == "generate":
        with ThreadPoolExecutor(a.jobs) as pool:
            ok = list(pool.map(lambda e: generate(e, a.events, a.seed), exps))
        sys.exit(0 if all(ok) else 1)
    if a.docker:
        args = [x for x in sys.argv[1:] if x != "--docker"]
        subprocess.run(["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}",
                        "-e", "HOME=/tmp", "-v", f"{REPO}:{REPO}", "-w", HERE, IMAGE,
                        "python3", os.path.join(HERE, "run_nuclear.py"), *args], check=True)
        return
    for exp in exps:
        analyze(cfg, exp)


if __name__ == "__main__":
    main()
