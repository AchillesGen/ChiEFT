"""The paper's figure style: LaTeX + Palatino, 15 pt labels, 13 pt ticks and
legends, a 5 x 4.75 in single panel with a 1:1 frame, and a 2x2 grid of square
panels.  "This work" is the one coloured line (C_OURS); data are open grey
markers.  The CC1pi figures colour each channel instead.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                     "font.serif": "Palatino"})

FT = 15
C_OURS = "#0057B7"
C_REF = "#999999"
C_DATA = "#555555"


def style_ax(ax, ft=FT):
    ax.tick_params("both", length=6, width=1, which="major", labelsize=ft - 2)
    ax.tick_params("both", length=4, width=1, which="minor")


def newfig():
    """Single panel with a 1:1 frame."""
    fig = plt.figure(figsize=(5, 4.75))
    ax = fig.add_axes([0.172, 0.137, 0.785, 0.826316])
    style_ax(ax)
    return fig, ax


def label_axes(ax, xlab, ylab):
    ax.set_xlabel(xlab, fontsize=FT)
    ax.set_ylabel(ylab, fontsize=FT)
    ax.xaxis.set_label_coords(0.5, -0.10)
    ax.yaxis.set_label_coords(-0.14, 0.5)


def grid22(sharex=True, sharey=True):
    """2x2 grid of 3.885 in square panels."""
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.85), sharex=sharex, sharey=sharey)
    fig.subplots_adjust(left=0.095, right=0.98, bottom=0.08, top=0.98,
                        wspace=0.05, hspace=0.05)
    for ax in axes.flat:
        style_ax(ax)
    return fig, axes


# CC1pi channels
CC_COLOR = {"ppip": "#0057b7", "ppi0": "#93264f", "npip": "#e8820c"}
BW_LABEL = r"+ Breit--Wigner $M_\Delta, \Gamma_0$, $\Gamma_\Delta \propto k^3$"
BW_SHORT = r"+ BW, $\Gamma_\Delta \propto k^3$"
