"""Per-event running-width weight for the Achilles samples.

Sherpa carries a constant Delta width.  The running-width prediction is the
same events, each nu_mu CC, (anti)neutrino NC or e single-pion hard process
weighted by

    nu_mu CC 1pi   w_ch(E_nu)               from running_width_cc.csv
    nu NC 1pi      w_ch(E_nu)               from running_width_nc.csv (any flavour)
    e 1pi          (S_T^bw + eps S_L^bw) / (S_T^pole + eps S_L^pole)
                   at the event's (W, Q^2, eps), from running_width_electron.csv

(all made by ../free_nucleon/tables.py), and 1 for every other event.  The
channel and kinematics are read off the primary vertex of the event graph,
i.e. the hard process before final-state interactions.
"""
import csv

import numpy as np

COLUMNS = ("rw.enu_GeV", "rw.channel", "rw.W_GeV", "rw.Q2_GeV2", "rw.eps")
PRIMARY = 1
NUCLEONS = (2212, 2112)
NEUTRINOS = (12, 14, 16, -12, -14, -16)
PIONS = (111, 211, -211)

# channel codes: 1-3 nu_mu CC 1pi, 4-7 e 1pi, 8-11 nu NC 1pi, 12-15 nubar NC
# 1pi, 0 anything else
CC_CHANNELS = {1: "ppip", 2: "ppi0", 3: "npip"}
E_CHANNELS = {4: "ppi0", 5: "npip", 6: "npi0", 7: "ppim"}
NC_CHANNELS = {8: "nc_ppi0", 9: "nc_npip", 10: "nc_npi0", 11: "nc_ppim",
               12: "nc_ppi0~", 13: "nc_npip~", 14: "nc_npi0~", 15: "nc_ppim~"}
_CC = {(2212, (211, 2212)): 1, (2112, (111, 2212)): 2, (2112, (211, 2112)): 3}
_NC = {(2212, (111, 2212)): 8, (2212, (211, 2112)): 9,
       (2112, (111, 2112)): 10, (2112, (-211, 2212)): 11}
_E = {(2212, (111, 2212)): 4, (2212, (211, 2112)): 5,
      (2112, (111, 2112)): 6, (2112, (-211, 2212)): 7}
NONE = (-1.0, 0, -1.0, -1.0, -1.0)


def _p4(p, unit):
    m = p.momentum()
    return np.array([m.e(), m.px(), m.py(), m.pz()]) * unit


def _mdot(a, b):
    return a[0] * b[0] - a[1:] @ b[1:]


def hard_process(ev):
    """(E_nu, channel code, W, Q^2, eps) of the event's primary vertex."""
    v = next((v for v in ev.vertices() if v.status() == PRIMARY), None)
    if v is None:
        return NONE
    unit = 1e-3 if "MEV" in str(ev.momentum_unit()).upper() else 1.0
    pin, out = v.particles_in(), v.particles_out()
    nin = [p.pid() for p in pin if p.pid() in NUCLEONS]
    hadp = [p for p in out if p.pid() in NUCLEONS or p.pid() in PIONS]
    had = tuple(sorted(p.pid() for p in hadp))
    if len(nin) != 1:
        return NONE
    nu = [p for p in pin if p.pid() in NEUTRINOS]
    if len(nu) == 1:
        lep = [p.pid() for p in out if abs(p.pid()) in (11, 12, 13, 14, 15, 16)]
        code = 0
        if lep == [nu[0].pid()]:                            # neutral current
            code = _NC.get((nin[0], had), 0)
            code += 4 if code and nu[0].pid() < 0 else 0
        elif nu[0].pid() == 14 and lep == [13]:
            code = _CC.get((nin[0], had), 0)
        return _p4(nu[0], unit)[0], code, -1.0, -1.0, -1.0
    ein = [p for p in pin if p.pid() == 11]
    eout = [p for p in out if p.pid() == 11]
    code = _E.get((nin[0], had), 0)
    if len(ein) != 1 or len(eout) != 1 or not code:
        return NONE
    k, kp = _p4(ein[0], unit), _p4(eout[0], unit)
    q = k - kp
    Q2 = -_mdot(q, q)
    X = sum(_p4(p, unit) for p in hadp)
    W = np.sqrt(max(_mdot(X, X), 0.0))
    # virtual-photon polarisation in the frame the event is written in (nucleus at rest)
    c = (k[1:] @ kp[1:]) / (np.linalg.norm(k[1:]) * np.linalg.norm(kp[1:]))
    tan2 = (1 - c) / (1 + c) if c > -1 else 1e30
    eps = 1.0 / (1.0 + 2.0 * (q[1:] @ q[1:]) / max(Q2, 1e-12) * tan2)
    return -1.0, code, W, Q2, eps


_cache = {"event": None, "value": NONE}


def _column(i):
    def get(ev):
        if ev.event_number() != _cache["event"]:
            _cache["event"], _cache["value"] = ev.event_number(), hard_process(ev)
        return _cache["value"][i]
    return get


def add_columns(fg):
    """Register COLUMNS on a pyNUISANCE EventFrameGen."""
    for i, name in enumerate(COLUMNS):
        (fg.add_int_column if name == "rw.channel" else fg.add_double_column)(name, _column(i))


class Weights:
    """Per-event running-width weights from the three tables."""

    def __init__(self, cc_table, nc_table, electron_table):
        self.nu = {}
        for path in (cc_table, nc_table):
            for r in csv.DictReader(open(path)):
                if r["wmax"] == "none":                 # nuclear samples have no W cut
                    self.nu.setdefault(r["channel"], []).append((float(r["E_nu"]),
                                                                 float(r["w"])))
        self.nu = {ch: np.array(sorted(v)).T for ch, v in self.nu.items()}
        rows = {}
        for r in csv.DictReader(open(electron_table)):
            rows.setdefault(r["channel"], []).append(
                [float(r[k]) for k in ("W", "Q2", "ST_pole", "SL_pole", "ST_bw", "SL_bw")])
        self.e = {}
        for ch, rs in rows.items():
            a = np.array(rs)
            Ws, Qs = np.unique(a[:, 0]), np.unique(a[:, 1])
            cube = np.zeros((4, Ws.size, Qs.size))
            for row in a:
                cube[:, np.searchsorted(Ws, row[0]), np.searchsorted(Qs, row[1])] = row[2:]
            self.e[ch] = (Ws, Qs, cube)

    @staticmethod
    def _bilinear(Ws, Qs, cube, W, Q2):
        W, Q2 = np.clip(W, Ws[0], Ws[-1]), np.clip(Q2, Qs[0], Qs[-1])
        i = np.clip(np.searchsorted(Ws, W) - 1, 0, Ws.size - 2)
        j = np.clip(np.searchsorted(Qs, Q2) - 1, 0, Qs.size - 2)
        tw, tq = (W - Ws[i]) / (Ws[i + 1] - Ws[i]), (Q2 - Qs[j]) / (Qs[j + 1] - Qs[j])
        return ((1 - tw) * (1 - tq) * cube[:, i, j] + tw * (1 - tq) * cube[:, i + 1, j]
                + (1 - tw) * tq * cube[:, i, j + 1] + tw * tq * cube[:, i + 1, j + 1])

    def __call__(self, cols):
        """Weights for the selected events, from their COLUMNS values."""
        code = np.asarray(cols["rw.channel"], dtype=int)
        w = np.ones(code.shape)
        for c, ch in {**CC_CHANNELS, **NC_CHANNELS}.items():
            m = code == c
            if m.any():
                E, R = self.nu[ch]
                w[m] = np.interp(np.asarray(cols["rw.enu_GeV"])[m], E, R)
        for c, ch in E_CHANNELS.items():
            m = code == c
            if m.any():
                tp, lp, tb, lb = self._bilinear(*self.e[ch], np.asarray(cols["rw.W_GeV"])[m],
                                                np.asarray(cols["rw.Q2_GeV2"])[m])
                eps = np.asarray(cols["rw.eps"])[m]
                den = tp + eps * lp
                w[m] = np.where(den > 0, (tb + eps * lb) / np.where(den > 0, den, 1), 1.0)
        return w
