#!/usr/bin/env python3
"""Free-nucleon figures of the paper, from output/ (sherpa_scan.py, tables.py)
and data/, written to ../figures/ and ../figures_appendix/ as in the paper.

    figures/photoprod_channels.pdf                 gamma N -> pi N, all channels
    figures/photoprod_pi0p.pdf                     gamma p -> p pi0 with the variants
    figures/cc1pi_channels.pdf                     nu_mu N -> mu- N' pi vs ANL/BNL
    figures_appendix/lineshape.pdf                 Delta line shape
    figures_appendix/photoprod_channels_appendix.pdf   all channels with the variants
    figures_appendix/photoprod_piNDelta_ff.pdf     effect of the piNDelta form factor
    figures_appendix/cc1pi_appendix.pdf            ANL/BNL with the running width
"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from chieft import style as ps  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

OUT = os.path.join(HERE, "output")
DATA = os.path.join(HERE, "data")
PAPER = os.path.dirname(HERE)

# 2x2 layout: proton target left, neutron right; charged pions on top
PANELS = [("ppim", r"$\gamma n \to p\pi^-$", "said_gn_ppiminus_SGT.csv"),
          ("npip", r"$\gamma p \to n\pi^+$", "said_gp_npiplus_SGT.csv"),
          ("ppi0", r"$\gamma p \to p\pi^0$", "said_gp_ppi0_SGT.csv"),
          ("npi0", r"$\gamma n \to n\pi^0$", "gn_npi0_integrated.csv")]
MARKERS = ["o", "s", "^", "v", "D", "p", "*", "<", ">", "h", "X", "d"]
GREYS = ["#4d4d4d", "#7f7f7f", "#9e9e9e", "#5d5d5d", "#8f8f8f", "#6f6f6f", "#a8a8a8",
         "#606060", "#909090", "#787878", "#525252", "#a0a0a0"]
C_VEC, C_BW = "#1B9E77", "#C62828"
LS_VEC, LS_BW = (0, (1.0, 1.3)), (0, (6.0, 2.0))
H_OURS = Line2D([], [], color=ps.C_OURS, lw=2.2, label="this work")
H_VEC = Line2D([], [], color=C_VEC, lw=2.4, ls=LS_VEC, label=r"$+\ \omega,\rho$ exchange")
H_BW = Line2D([], [], color=C_BW, lw=2.4, ls=LS_BW,
              label=r"$+$ Breit-Wigner $M_\Delta,\Gamma_0$, $\Gamma_\Delta \propto k^3$")
H_DATA = Line2D([], [], color="#7f7f7f", marker="o", markerfacecolor="none", ls="none",
                ms=3.5, label=r"data (refs.\ in caption)")


def read(path):
    with open(path) as f:
        return list(csv.DictReader(r for r in f if not r.startswith("#")))


def curve(rows, ch, key="channel", x="E", y="sigma_ub", **match):
    pts = sorted((float(r[x]), float(r[y])) for r in rows
                 if r[key] == ch and all(r[k] == v for k, v in match.items())
                 and float(r[y]) > 0)
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])


def photo_data(ax, datafile):
    pts = defaultdict(lambda: ([], [], []))
    for r in read(os.path.join(DATA, datafile)):
        pts[r["experiment"]][0].append(float(r["E_gamma_MeV"]) / 1e3)
        pts[r["experiment"]][1].append(float(r["sigma_exp_ub"]))
        pts[r["experiment"]][2].append(float(r["err_ub"]))
    for ((E, s, e), m, c) in zip((v for _, v in sorted(pts.items())), MARKERS, GREYS):
        ax.errorbar(E, s, yerr=e, fmt=m, ms=3.5, capsize=1.5, lw=0.8, ls="none",
                    markerfacecolor="none", color=c, zorder=2)


LS_NOFF = (0, (5.0, 1.5, 1.0, 1.5))


def photo_grid(sherpa=None, variants=None, curves=()):
    """2x2 photoproduction panels: data, `curves` from photo_curves.csv as
    (name, plot kwargs), and the Sherpa result."""
    fig, axes = ps.grid22()
    for ax, (ch, title, datafile) in zip(axes.flat, PANELS):
        photo_data(ax, datafile)
        for name, kw in curves:
            ax.plot(*curve(variants, ch, x="Egamma", curve=name), **kw)
        if sherpa is not None:
            ax.plot(*curve(sherpa, ch), color=ps.C_OURS, lw=2.2, zorder=4)
        ax.text(0.95, 0.92, title, transform=ax.transAxes, ha="right", va="top",
                fontsize=ps.FT)
    axes[0, 0].set_xlim(0.14, 0.6)
    axes[0, 0].set_ylim(0, 419)
    for ax in axes[1, :]:
        ax.set_xlabel(r"$E_\gamma$ (GeV)", fontsize=ps.FT)
        ax.xaxis.set_label_coords(0.5, -0.10)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\sigma$ ($\mu$b)", fontsize=ps.FT)
        ax.yaxis.set_label_coords(-0.16, 0.5)
    return fig, axes


def fig_photoprod_channels(sherpa):
    fig, axes = photo_grid(sherpa)
    axes[0, 1].legend(handles=[H_OURS, H_DATA], fontsize=ps.FT - 2, loc="upper left",
                      frameon=False, bbox_to_anchor=(0.03, 0.97))
    return fig


def legend_strip(fig, handles, ncol):
    """Legend in a strip above the panels: grow the figure, keep the panels square."""
    W0, H0 = fig.get_size_inches()
    H1 = H0 + 0.80
    panel = (0.98 - 0.095) * W0 / 2.05
    bottom = 0.08 * H0 / H1
    fig.set_size_inches(W0, H1)
    fig.subplots_adjust(bottom=bottom, top=bottom + panel * 2.05 / H1)
    fig.legend(handles=handles, fontsize=ps.FT - 2, loc="upper center", ncol=ncol,
               frameon=False, bbox_to_anchor=(0.5375, 0.997), columnspacing=1.6,
               handlelength=2.6)


def fig_photoprod_channels_appendix(sherpa, variants):
    fig, _ = photo_grid(sherpa, variants, [
        ("vec", dict(color=C_VEC, lw=2.4, ls=LS_VEC, zorder=3)),
        ("bw", dict(color=C_BW, lw=2.4, ls=LS_BW, zorder=5))])
    legend_strip(fig, [H_OURS, H_BW, H_VEC, H_DATA], ncol=2)
    return fig


def fig_photoprod_piNDelta_ff(variants):
    """Both curves semi-analytic, so their difference is the form factor alone."""
    fig, _ = photo_grid(None, variants, [
        ("no_ff", dict(color=ps.C_OURS, lw=1.8, ls=LS_NOFF, zorder=3)),
        ("base", dict(color=ps.C_OURS, lw=2.2, zorder=4))])
    legend_strip(fig, [H_OURS, Line2D([], [], color=ps.C_OURS, lw=1.8, ls=LS_NOFF,
                                      label=r"$F_{\pi N\Delta} = 1$"), H_DATA], ncol=3)
    return fig


def fig_photoprod_pi0p(sherpa, variants):
    fig, ax = ps.newfig()
    photo_data(ax, "said_gp_ppi0_SGT.csv")
    ax.plot(*curve(variants, "ppi0", x="Egamma", curve="vec"), color=C_VEC, lw=2.4,
            ls=LS_VEC, zorder=3)
    ax.plot(*curve(variants, "ppi0", x="Egamma", curve="bw"), color=C_BW, lw=2.4,
            ls=LS_BW, zorder=5)
    ax.plot(*curve(sherpa, "ppi0"), color=ps.C_OURS, lw=2.2, zorder=4)
    ax.set_xlim(0.14, 0.6)
    ax.set_ylim(0, 580)
    ax.legend(handles=[H_OURS, H_VEC, H_BW], fontsize=ps.FT - 2, loc="upper right",
              frameon=False, handlelength=2.6, labelspacing=0.3, borderaxespad=0.4)
    ps.label_axes(ax, r"$E_\gamma$ (GeV)", r"$\sigma(\gamma p \to p\pi^0)$ ($\mu$b)")
    return fig


def fig_lineshape(variants):
    MN, MD = 0.938272, 1.209
    fig, ax = ps.newfig()
    egM = (MD ** 2 - MN ** 2) / (2 * MN)
    ax.axvline(egM, color=ps.C_REF, lw=0.8, zorder=1)
    ax.text(egM - 0.004, 25, r"$W = m_\Delta$", rotation=90, ha="right", va="bottom",
            color=ps.C_REF, fontsize=ps.FT - 4)
    styles = [("bw_formula", "k", 1.6, (0, (3, 2)),
               r"Breit–Wigner"),
              ("delta_k3", C_BW, 2.2, (0, (6, 2)), r"$\Gamma_\Delta \propto k^3$"),
              ("delta_const", ps.C_OURS, 2.2, "-",
               r"constant $\Gamma_\Delta$ (this work)")]
    for z, (name, c, lw, ls, _) in enumerate(styles):
        ax.plot(*curve(variants, "ppi0", x="Egamma", curve=name), color=c, lw=lw, ls=ls,
                zorder=3 + z)
    ax.set_xlim(0.14, 0.6)
    ax.set_ylim(0, 580)
    h = [Line2D([], [], color=c, lw=lw, ls=ls, label=lab) for _, c, lw, ls, lab in styles]
    ax.legend(handles=[h[2], h[0], h[1]], fontsize=ps.FT - 2, loc="upper right",
              frameon=False, handlelength=2.4, labelspacing=0.3, borderaxespad=0.4)
    ps.label_axes(ax, r"$E_\gamma$ (GeV)", r"$\sigma_\Delta(\gamma p \to p\pi^0)$ ($\mu$b)")
    return fig


# ------------------------------------------------------------------- ANL/BNL
CC = {"ppip": r"$\nu_\mu\, p \to \mu^-\, p\, \pi^+$",
      "ppi0": r"$\nu_\mu\, n \to \mu^-\, p\, \pi^0$",
      "npip": r"$\nu_\mu\, n \to \mu^-\, n\, \pi^+$"}
UB = 1e8                                   # ub -> 1e-38 cm^2


def anl_bnl():
    """{(channel, experiment): (E_mid, E_lo, E_hi, sigma, err)}, W < 1.4 GeV,
    Wilkinson et al. re-analysis, in 1e-38 cm^2."""
    out = defaultdict(list)
    for r in read(os.path.join(DATA, "anl_bnl_cc1pip.csv")):
        if r["wcut"] == "1.4" and r["correction"] == "wilkinson" \
                and 0.25 <= float(r["E_mid"]) <= 1.65:
            out[(r["channel"], r["experiment"])].append(
                [float(r[k]) for k in ("E_mid", "E_lo", "E_hi")]
                + [float(r["sigma_ub"]) * UB, float(r["dsigma_ub"]) * UB])
    return {k: np.array(sorted(v)).T for k, v in out.items()}


def fig_cc1pi(sherpa, width=None):
    """width: the running_width_cc.csv rows, to add the thin running-width curves."""
    fig, ax = ps.newfig()
    model, variant = {}, {}
    for ch in CC:
        E, S = curve(sherpa, ch)
        model[ch] = (E, S * UB)
        ax.plot(E, S * UB, color=ps.CC_COLOR[ch], lw=2.0, zorder=4, solid_capstyle="round")
        if width is not None:
            Ew, w = curve(width, ch, x="E_nu", y="w", wmax="1.4")
            variant[ch] = (E, S * UB * np.interp(E, Ew, w))
            ax.plot(*variant[ch], color=ps.CC_COLOR[ch], lw=1.1, alpha=0.75, zorder=4)
    for (ch, exp), (Em, Elo, Ehi, S, dS) in anl_bnl().items():
        c = ps.CC_COLOR[ch]
        ax.errorbar(Em, S, xerr=[Em - Elo, Ehi - Em], yerr=dS,
                    fmt="o" if exp == "ANL" else "s", ms=4.5,
                    mfc=c if exp == "ANL" else "white", mec=c, mew=0.9, ecolor=c,
                    elinewidth=0.9, capsize=2.0, capthick=0.9, alpha=0.85, zorder=5)
    top = max(float(S.max()) for _, S in list(model.values()) + list(variant.values()))
    ax.set_xlim(0.3, 1.6)
    ax.set_ylim(0.0, max(1.0, np.ceil(top * 1.15 * 10) / 10))
    ps.label_axes(ax, r"$E_\nu$ (GeV)", r"$\sigma$ ($10^{-38}\,\mathrm{cm}^2$)")

    # direct channel labels, placed clear of the curves and the data
    for ch, lab in CC.items():
        if ch == "ppip" and variant:
            x, y = 1.12, 0.5 * (np.interp(1.3, *model[ch]) + np.interp(1.3, *variant[ch])) + 0.03
        elif ch == "ppip":
            x, y = 1.05, 0.5 * (np.interp(1.25, *model[ch]) + 0.63) + 0.03
        elif ch == "npip":
            x, y = 0.95, 0.062
        else:
            x = 1.0
            span = np.linspace(x, x + 0.45, 20)
            y = max(np.interp(span, *c).max() for c in [model[ch]]
                    + ([variant[ch]] if variant else [])) + 0.085
        ax.text(x, y, lab, color=ps.CC_COLOR[ch], ha="left", va="top", fontsize=ps.FT - 2)
    ax.text(0.33, 0.60, r"$W < 1.4$ GeV", color=ps.C_REF, ha="left", va="center",
            fontsize=ps.FT - 2)
    ink = "#3a3a3a"
    h = [Line2D([], [], color=ink, lw=2.0,
                label="constant width" if variant else "this work")]
    if variant:
        h.append(Line2D([], [], color=ink, lw=1.1, alpha=0.75, label=ps.BW_LABEL))
    h += [Line2D([], [], color=ink, marker="o", ls="none", ms=4.5, mfc=ink,
                 label="ANL (re-analyzed)"),
          Line2D([], [], color=ink, marker="s", ls="none", ms=4.5, mfc="white",
                 label="BNL (re-analyzed)")]
    ax.legend(handles=h, loc="upper left", handlelength=2.0, labelspacing=0.35,
              borderaxespad=0.4, frameon=False, fontsize=ps.FT - 2)
    return fig


def main():
    def have(name):
        path = os.path.join(OUT, name)
        if not os.path.exists(path):
            print(f"skipping figures that need output/{name}")
            return None
        return read(path)

    photo, variants = have("sherpa_photo.csv"), have("photo_curves.csv")
    cc, width = have("sherpa_cc.csv"), have("running_width_cc.csv")
    figs = {}
    if photo:
        figs["figures/photoprod_channels"] = lambda: fig_photoprod_channels(photo)
    if variants:
        figs["figures_appendix/lineshape"] = lambda: fig_lineshape(variants)
        figs["figures_appendix/photoprod_piNDelta_ff"] = \
            lambda: fig_photoprod_piNDelta_ff(variants)
    if photo and variants:
        figs["figures/photoprod_pi0p"] = lambda: fig_photoprod_pi0p(photo, variants)
        figs["figures_appendix/photoprod_channels_appendix"] = \
            lambda: fig_photoprod_channels_appendix(photo, variants)
    if cc:
        figs["figures/cc1pi_channels"] = lambda: fig_cc1pi(cc)
    if cc and width:
        figs["figures_appendix/cc1pi_appendix"] = lambda: fig_cc1pi(cc, width)
    for name, make in figs.items():
        path = os.path.join(PAPER, name + ".pdf")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        make().savefig(path)
        print("wrote", path)


if __name__ == "__main__":
    main()
