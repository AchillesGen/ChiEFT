#!/usr/bin/env python3
"""Every NUISANCE sample of analyses.yml (results/*.json) as a single panel,
including those not shown in the paper -> plots/<sample>.pdf.

Thick line: this work; grey dashed: with the running width; open markers: data.
"""
import json
import os
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from chieft import style as ps  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


def draw(res, labels):
    raw = np.asarray(res["data"], float)
    exponent = int(np.floor(np.log10(np.max(np.abs(raw))))) if np.any(raw) else 0
    sc = 10.0 ** exponent
    data = raw / sc
    derr = np.sqrt(np.clip(np.diag(np.asarray(res["data_cov"], float)), 0, None)) / sc
    edges = np.asarray(res["edges"], float) if res["edges"] else np.arange(data.size + 1.0)
    x = 0.5 * (edges[:-1] + edges[1:])

    fig, ax = ps.newfig()
    curves = [("nominal", "this work", dict(color=ps.C_OURS, lw=2.0, ls="-", zorder=4)),
              ("running_width", r"+ running $\Gamma_\Delta$",
               dict(color=ps.C_REF, lw=1.4, ls=(0, (5, 2.5)), zorder=3))]
    vals = []
    for key, _, st in curves:
        v = np.asarray(res["predictions"][key]["values"], float) / sc
        ax.stairs(v, edges, baseline=None, **st)
        vals.append(v)
    ax.errorbar(x, data, yerr=derr, xerr=[x - edges[:-1], edges[1:] - x], fmt="o", ms=4.5,
                mfc="white", mec=ps.C_DATA, mew=1.0, ecolor=ps.C_DATA, elinewidth=0.9,
                capsize=2.0, capthick=0.9, zorder=6)
    sym, unit = labels.get(res["measurement"], ("x", ""))
    per_x = f"/{unit}" if unit else ""
    ps.label_axes(ax, f"{sym} ({unit})" if unit else sym,
                  f"$d\\sigma/d${sym} ($10^{{{exponent}}}$ cm$^2${per_x}/{res['target_label']})")
    ax.set_xlim(edges[0], edges[-1])
    top = max(float(np.nanmax(data + derr)), max(float(np.nanmax(v)) for v in vals))
    ax.set_ylim(min(0.0, float(np.nanmin(data - derr))), top * 1.3)
    handles = [Line2D([], [], label=lab, **{k: st[k] for k in ("color", "lw", "ls")})
               for _, lab, st in curves]
    handles.append(Line2D([], [], color=ps.C_DATA, marker="o", ls="none", ms=4.5,
                          mfc="white", label=res["experiment"].split("_")[0]))
    ax.legend(handles=handles, loc="best", fontsize=ps.FT - 2, frameon=False,
              handlelength=2.2, labelspacing=0.3)
    return fig


def main():
    cfg = yaml.safe_load(open(os.path.join(HERE, "analyses.yml")))
    labels = {k: tuple(v) for k, v in cfg["labels"].items()}
    out = os.path.join(HERE, "plots")
    os.makedirs(out, exist_ok=True)
    for name in [s for e in cfg["experiments"] for s in e["samples"]]:
        p = os.path.join(HERE, "results", name + ".json")
        if not os.path.exists(p):
            continue
        fig = draw(json.load(open(p)), labels)
        fig.savefig(os.path.join(out, name.replace(":", "_") + ".pdf"))
        plt.close(fig)
        print("plotted", name)


if __name__ == "__main__":
    main()
