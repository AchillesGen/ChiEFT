#!/usr/bin/env python3
"""Nuclear and electron-scattering figures of the paper, from results/*.json
(run_nuclear.py), written to ../figures/ and ../figures_appendix/.

    figures/nuclear_minerva_pi0.pdf        MINERvA CC1pi0
    figures/nuclear_t2k_pip.pdf            T2K CC1pi+ and CC1pi+Np
    figures/electron_carbon.pdf            inclusive e-12C
    figures/nuclear_microboone_ncpi0.pdf   MicroBooNE NCpi0
    figures_appendix/nuclear_appendix.pdf  MiniBooNE and MINERvA CC1pi

Thick line: this work (constant width); thin: with the running width.
Annotations and legends are placed where they cover no curve or data point,
extending the y range if necessary.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from chieft import style as ps  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

PAPER = os.path.dirname(HERE)
FTS = ps.FT - 2
NU_CH = r"$\nu_\mu\,{\rm CH} \to \mu^-\pi^+ X$"
ELECTRON = "$E = %s$ GeV, $\\theta = %s^\\circ$"

# (sample, x label, y unit, x range, annotation, legend here?)
FIGS = {
    "figures/nuclear_minerva_pi0": (ps.BW_LABEL, [
        ("MINERvA_CC1pi0_XSec_1DWexp_nu", r"$W_{\rm exp}$ (GeV)", "GeV", (0.6, 1.8),
         r"$\nu_\mu\,{\rm CH} \to \mu^-\pi^0 X$", False),
        ("MINERvA_CC1pi0_XSec_1DPPi0MassDelta_nu", r"$M_{p\pi^0}$ (GeV)", "GeV", (1.0, 1.6),
         r"$W_{\rm exp} < 1.4$ GeV", False),
        ("MINERvA_CC1pi0_XSec_1DTpi_nu", r"$T_\pi$ (GeV)", "GeV", (0.0, 1.0), None, True),
        ("MINERvA_CC1pi0_XSec_1DQ2_nu", r"$Q^2$ (GeV$^2$)", "GeV$^2$", (0.0, 2.0), None, False)]),
    "figures/nuclear_t2k_pip": (ps.BW_LABEL, [
        ("T2K_CC1pip_CH_XSec_1Dppi_nu", r"$p_\pi$ (GeV)", "GeV", (0.0, 2.0), NU_CH, False),
        ("T2K_CC1pip_CH_XSec_1DQ2_nu", r"$Q^2$ (GeV$^2$)", "GeV$^2$", (0.0, 2.0), None, True),
        ("T2K_CC1pipNp_CH_XSec_1DdpTT_nu", r"$\delta p_{TT}$ (MeV)", "MeV", (-700, 700),
         r"$\mu^-\pi^+ p$", False),
        ("T2K_CC1pipNp_CH_XSec_1DpN_nu", r"$p_N$ (MeV)", "MeV", (0, 1500),
         r"$\mu^-\pi^+ p$", False)]),
    "figures_appendix/nuclear_appendix": (ps.BW_SHORT, [
        ("MiniBooNE_CC1pip_XSec_1DTpi_nu", r"$T_\pi$ (MeV)", "MeV", (0, 400),
         r"MiniBooNE $\mu^-\pi^+$", False),
        ("MiniBooNE_CC1pi0_XSec_1Dppi0_nu", r"$p_\pi$ (GeV)", "GeV", (0.0, 1.0),
         r"MiniBooNE $\mu^-\pi^0$", True),
        ("MINERvA_CC1pip_XSec_1DTpi_nu_2017", r"$T_\pi$ (MeV)", "MeV", (0, 350),
         r"MINERvA $\mu^-\pi^\pm$, $W < 1.4$ GeV", False),
        ("MINERvA_CC1pip_XSec_1DQ2_nu_2017", r"$Q^2$ (GeV$^2$)", "GeV$^2$", (0.0, 2.0),
         r"MINERvA $\mu^-\pi^\pm$, $W < 1.4$ GeV", False)]),
    "figures/electron_carbon": (ps.BW_SHORT, [
        ("ElectronData_6_12_0.560_36.000_Barreau:1983ht", r"$\omega$ (GeV)", "GeV",
         (0.0, 0.40), ELECTRON % ("0.560", "36"), True),
        ("ElectronData_6_12_0.560_60.000_Barreau:1983ht", r"$\omega$ (GeV)", "GeV",
         (0.0, 0.46), ELECTRON % ("0.560", "60"), False),
        ("ElectronData_6_12_0.961_37.500_Sealock:1989nx", r"$\omega$ (GeV)", "GeV",
         (0.05, 0.68), ELECTRON % ("0.961", "37.5"), False),
        ("ElectronData_6_12_1.299_37.500_Sealock:1989nx", r"$\omega$ (GeV)", "GeV",
         (0.15, 0.82), ELECTRON % ("1.299", "37.5"), False)]),
    "figures/nuclear_microboone_ncpi0": (ps.BW_SHORT, [
        ("MicroBooNE_NCpi0_XSec_XpPpi0_nu", r"$p_{\pi^0}$ (GeV)", "GeV", (0.0, 1.2),
         r"$\nu\,{\rm Ar} \to \nu\,\pi^0 X$", True),
        ("MicroBooNE_NCpi0_XSec_XpCosThetaPi0_nu", r"$\cos\theta_{\pi^0}$", "", (-1.0, 1.0),
         None, False)]),
}
ANNOT_LOCS = [(0.95, 0.95, "right", "top"), (0.05, 0.95, "left", "top"),
              (0.50, 0.95, "center", "top"), (0.95, 0.80, "right", "top")]
LEG_LOCS = ["upper right", "upper left", "center right", "upper center",
            "center left", "lower right", "lower left"]


def stair_points(ax, edges, vals):
    xs, ys = [], []
    for i, v in enumerate(vals):
        for f in np.linspace(0, 1, 16):
            xs.append(edges[i] + f * (edges[i + 1] - edges[i])); ys.append(v)
    for i in range(1, len(vals)):
        for f in np.linspace(0, 1, 10):
            xs.append(edges[i]); ys.append(vals[i - 1] + f * (vals[i] - vals[i - 1]))
    return ax.transData.transform(np.column_stack([xs, ys]))


def data_boxes(ax, x, y, dy, lo, hi):
    out = []
    for xi, yi, di, a, b in zip(x, y, dy, lo, hi):
        p0 = ax.transData.transform((a, yi - di)); p1 = ax.transData.transform((b, yi + di))
        out.append((min(p0[0], p1[0]) - 3, min(p0[1], p1[1]) - 3,
                    max(p0[0], p1[0]) + 3, max(p0[1], p1[1]) + 3))
    return out


def hits(bb, pts, boxes, pad=2):
    x0, y0, x1, y1 = bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad
    if len(pts) and np.any((pts[:, 0] > x0) & (pts[:, 0] < x1) &
                           (pts[:, 1] > y0) & (pts[:, 1] < y1)):
        return True
    return any(not (b[2] < x0 or b[0] > x1 or b[3] < y0 or b[1] > y1) for b in boxes)


def panel(fig, ax, res, xlab, yunit, xlim, annot, want_legend, legend_rw):
    raw = np.asarray(res["data"], float)
    exp = int(np.floor(np.log10(np.max(np.abs(raw)))))
    sc = 10.0 ** exp
    data = raw / sc
    derr = np.sqrt(np.clip(np.diag(np.asarray(res["data_cov"], float)), 0, None)) / sc
    edges = np.asarray(res["edges"], float)
    x = 0.5 * (edges[:-1] + edges[1:])
    ours = np.asarray(res["predictions"]["nominal"]["values"]) / sc
    rw = np.asarray(res["predictions"]["running_width"]["values"]) / sc
    ax.stairs(ours, edges, baseline=None, color=ps.C_OURS, lw=2.0, zorder=4)
    curves = [ours]
    if not np.allclose(rw, ours):
        ax.stairs(rw, edges, baseline=None, color=ps.C_OURS, lw=1.0, zorder=3)
        curves.append(rw)
    ax.errorbar(x, data, yerr=derr, xerr=[x - edges[:-1], edges[1:] - x], fmt="o",
                ms=3.5, mfc="white", mec=ps.C_DATA, mew=0.9, ecolor=ps.C_DATA,
                elinewidth=0.8, capsize=1.5, capthick=0.8, zorder=6)
    sym = xlab.split(" (")[0]
    per_x = f"/{yunit}" if yunit else ""
    ax.set_xlabel(xlab, fontsize=ps.FT)
    ax.set_ylabel(f"$d\\sigma/d${sym} ($10^{{{exp}}}$ cm$^2${per_x}/{res['target_label']})",
                  fontsize=ps.FT)
    xlim = (max(xlim[0], edges[0]), min(xlim[1], edges[-1]))
    ax.set_xlim(*xlim)
    inside = (x >= xlim[0]) & (x <= xlim[1])
    shown = (edges[:-1] < xlim[1]) & (edges[1:] > xlim[0])
    top = max(float(np.max((data + derr)[inside])), max(float(np.max(c[shown])) for c in curves))
    handles = [Line2D([], [], color=ps.C_OURS, lw=2.0, label="this work")]
    if len(curves) > 1:
        handles.append(Line2D([], [], color=ps.C_OURS, lw=1.0, label=legend_rw))
    handles.append(Line2D([], [], color=ps.C_DATA, marker="o", ls="none", ms=3.5,
                          mfc="white", label="data"))
    for step in range(6):
        ax.set_ylim(0, top * 1.1 * 1.2 ** step)
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        pts = np.vstack([stair_points(ax, edges, c) for c in curves])
        boxes = data_boxes(ax, x, data, derr, edges[:-1], edges[1:])
        placed, ok = [], not annot
        for xa, ya, ha, va in ANNOT_LOCS if annot else []:
            t = ax.text(xa, ya, annot, transform=ax.transAxes, ha=ha, va=va, fontsize=FTS)
            fig.canvas.draw()
            bb = t.get_window_extent(r)
            if not hits(bb, pts, boxes):
                placed.append((bb.x0, bb.y0, bb.x1, bb.y1)); ok = True
                break
            t.remove()
        if ok and want_legend:
            ok = False
            for loc in LEG_LOCS:
                leg = ax.legend(handles=handles, loc=loc, fontsize=FTS, handlelength=2.0,
                                labelspacing=0.3, borderaxespad=0.6, frameon=False)
                fig.canvas.draw()
                if not hits(leg.get_window_extent(r), pts, boxes + placed):
                    ok = True
                    break
                leg.remove()
        if ok:
            return
        for t in list(ax.texts):
            t.remove()
        if ax.get_legend():
            ax.get_legend().remove()
    print(f"  warning: labels overlap in {res['measurement']}")


def make(stem, legend_rw, panels):
    if len(panels) == 4:
        fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.85))
        fig.subplots_adjust(left=0.115, right=0.965, bottom=0.075, top=0.985,
                            wspace=0.40, hspace=0.24)
    else:
        H = 3.885 + 0.664 + 0.133
        fig, axes = plt.subplots(1, 2, figsize=(9.0, H))
        fig.subplots_adjust(left=0.115, right=0.965, bottom=0.664 / H, top=1 - 0.133 / H,
                            wspace=0.40)
    for ax, (name, xlab, yunit, xlim, annot, leg) in zip(axes.flat, panels):
        ps.style_ax(ax)
        res = json.load(open(os.path.join(HERE, "results", name + ".json")))
        panel(fig, ax, res, xlab, yunit, xlim, annot, leg, legend_rw)
        ax.yaxis.set_label_coords(-0.21, 0.5)
        ax.xaxis.set_label_coords(0.5, -0.11)
    path = os.path.join(PAPER, stem + ".pdf")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print("wrote", path)


if __name__ == "__main__":
    for stem, (legend_rw, panels) in FIGS.items():
        if all(os.path.exists(os.path.join(HERE, "results", p[0] + ".json")) for p in panels):
            make(stem, legend_rw, panels)
        else:
            print(f"skipping {stem}: results missing")
