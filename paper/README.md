# Reproducing the paper figures

Everything needed to regenerate the figures of the ChiPT paper from the model in
this repository (`../ChiPT_UFO`, `../param_card.dat`). Generated files go to
`figures/` and `figures_appendix/`, with the same file names as in the paper.

| Figure | Script | Inputs |
|---|---|---|
| `figures/photoprod_channels` | `free_nucleon/figures.py` | Sherpa photo scan, SAID data |
| `figures/photoprod_pi0p` | `free_nucleon/figures.py` | Sherpa photo scan, `tables.py photo`, SAID data |
| `figures/cc1pi_channels` | `free_nucleon/figures.py` | Sherpa CC scan, ANL/BNL data |
| `figures/nuclear_minerva_pi0` | `nuclear/figures.py` | Achilles + NUISANCE |
| `figures/nuclear_t2k_pip` | `nuclear/figures.py` | Achilles + NUISANCE |
| `figures/electron_carbon` | `nuclear/figures.py` | Achilles + NUISANCE |
| `figures/nuclear_microboone_ncpi0` | `nuclear/figures.py` | Achilles + NUISANCE |
| `figures_appendix/lineshape` | `free_nucleon/figures.py` | `tables.py photo` |
| `figures_appendix/photoprod_channels_appendix` | `free_nucleon/figures.py` | Sherpa photo scan, `tables.py photo` |
| `figures_appendix/photoprod_piNDelta_ff` | `free_nucleon/figures.py` | `tables.py photo` |
| `figures_appendix/cc1pi_appendix` | `free_nucleon/figures.py` | Sherpa CC scan, `tables.py cc` |
| `figures_appendix/nuclear_appendix` | `nuclear/figures.py` | Achilles + NUISANCE |

`nuclear/all_samples.py` also plots every other NUISANCE sample in
`nuclear/analyses.yml` (MicroBooNE CC, MINERvA and T2K STV, ¹⁶O electron
scattering, ...), including the ones not shown in the paper.

## Pipeline

```sh
make tables      # semi-analytic tables: running width, rho/omega, line shape
make sherpa      # gamma N -> pi N and nu_mu N -> mu N' pi with Sherpa
make events      # Achilles events for the nuclear and electron samples
make analysis    # NUISANCE binning of the events (in docker)
make figures
```

Each step writes its results to disk and later steps read them. `make figures`
draws whatever its inputs allow.

1. **`free_nucleon/tables.py`** evaluates the UFO directly (package `chipt/`).
   It covers what Sherpa cannot do: the running Δ width, ρ/ω exchange, and the
   Δ-only line shape. It needs Python only (numpy, matplotlib, pyyaml) and runs
   in parallel on half the cores (`--jobs` to change), one BLAS thread per
   job, so `--jobs N` uses N cores. `cc` and `nc` take
   about half a day each on 7 cores; the others take minutes.
   - `photo`: photoproduction variants.
   - `cc`, `nc`: the CC and NC (ν and ν̄) running-width factors w(E_ν).
   - `electron`: the virtual-photon responses used to reweight electron events.
2. **`free_nucleon/sherpa_scan.py`** computes the free-nucleon cross sections
   with Sherpa at constant Δ width. These are the "this work" curves. The
   `$SHERPA` variable selects the executable.
3. **`nuclear/run_nuclear.py generate`** runs Achilles with the cards in
   `nuclear/cards/`. Set `$ACHILLES_BUILD` to an Achilles build directory
   containing `bin/achilles` and `data/`.
4. **`nuclear/run_nuclear.py analyze --docker`** bins the events with
   pyNUISANCE. The image `docker/Dockerfile` contains NUISANCE with the fixes in
   `nuclear/patches/`. The running-width prediction reweights each CC 1π, NC 1π
   and e 1π hard process with the step-1 tables (`nuclear/running_width.py`).

Both generators must load the Sherpa model built from `../ChiPT_UFO` with the
form-factor plugin (top-level README). Both read `../param_card.dat`, which has
C4AHat = 0 as in the paper.

## Versions used for the paper

| Component | Version |
|---|---|
| FeynRules | 2.3.49 + `../feynrules/feynrules-2.3.49-chipt.patch`, Mathematica 14.3 |
| Sherpa | `570-update-achilles-interface` branch, commit 1cd6f079f |
| Achilles | `sherpa3-0` branch, commit ff313b48 |
| NUISANCE | NUISANCE2 f5f79db6 + `nuclear/patches/`, in `ghcr.io/achillesgen/achilles-physval@sha256:5e1d35ea…` (see `nuclear/docker/Dockerfile`) |

## Data

- `free_nucleon/data/said_*.csv`: total photoproduction cross sections from the
  SAID database. Sources are listed in each file header.
- `free_nucleon/data/gn_npi0_integrated.csv`: γn → nπ⁰, for which SAID has no
  total cross section. We integrated the SAID angular distributions, as
  described in its header.
- `free_nucleon/data/anl_bnl_cc1pip.csv`: ANL and BNL CC1π, as published and
  as re-analysed by Wilkinson et al., PRD 90, 112017 (2014).
- `nuclear/flux/`: the experiments' neutrino fluxes.
- `nuclear/smearing/`: the MicroBooNE NCπ⁰ Wiener-SVD A_C matrices.
- All other measurements come from NUISANCE.
