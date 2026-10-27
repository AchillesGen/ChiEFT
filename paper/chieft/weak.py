"""nu_mu N -> mu- N' pi (CC) and nu N -> nu N' pi (NC) from the UFO: the
leptonic current, the hadronic current and the total cross section integrated
over three-body phase space.

Diagrams of W+(q) N -> pi N', as Comix builds them: s- and u-channel nucleon,
t-channel pion, the four-point contact, s- and u-channel Delta and, where the
nucleon keeps its charge, W+ -> pi+ gamma*/Z* absorbed by the nucleon.  The NC
channels (nc_*) exchange a Z instead; a pi0 then has neither pion-in-flight
nor contact term (no Z pi0 pi0 vertex), and there are no seagulls.  The
neutrino is purely left-handed, so the spin average is over the nucleon only.
"""
import numpy as np

from .ufo import (G, HBARC2, METRIC, build, build_bosonic, delta_line, dot,
                  kmom, prop_fermion, prop_scalar, prop_vector, ubar, uspinor)

# Legs named as in the UFO.  N_s/N_u None: no nucleon of that charge exists for
# the s/u channel.  pi_t is the exchanged pion of the pion-in-flight diagram.
CC_CHANNELS = {
    "ppip": dict(Ni="p", Nf="p", pion="Pi+", pi_t="Pi0",
                 N_s=None, N_u="n", D_s="Delta++", D_u="Delta0"),
    "ppi0": dict(Ni="n", Nf="p", pion="Pi0", pi_t="Pi+",
                 N_s="p", N_u="n", D_s="Delta+", D_u="Delta0"),
    "npip": dict(Ni="n", Nf="n", pion="Pi+", pi_t="Pi0",
                 N_s="p", N_u=None, D_s="Delta+", D_u="Delta-"),
}
# Neutral current, Z(q) N -> pi N'.  The Z couples to every neutrino flavour
# alike, so one channel serves nu_e and nu_mu.
NC_CHANNELS = {
    "nc_ppi0": dict(Ni="p", Nf="p", pion="Pi0", pi_t=None, contact=False,
                    N_s="p", N_u="p", D_s="Delta+", D_u="Delta+"),
    "nc_npip": dict(Ni="p", Nf="n", pion="Pi+", pi_t="Pi-",
                    N_s="p", N_u="n", D_s="Delta+", D_u="Delta0"),
    "nc_npi0": dict(Ni="n", Nf="n", pion="Pi0", pi_t=None, contact=False,
                    N_s="n", N_u="n", D_s="Delta0", D_u="Delta0"),
    "nc_ppim": dict(Ni="n", Nf="p", pion="Pi-", pi_t="Pi+",
                    N_s="n", N_u="p", D_s="Delta0", D_u="Delta+"),
}
for _c in NC_CHANNELS.values():
    _c["boson"] = "Z"
CHANNELS = {**CC_CHANNELS, **NC_CHANNELS}
LABEL = {"ppip": r"$\nu_\mu\, p \to \mu^-\, p\, \pi^+$",
         "ppi0": r"$\nu_\mu\, n \to \mu^-\, p\, \pi^0$",
         "npip": r"$\nu_\mu\, n \to \mu^-\, n\, \pi^+$"}
MMU = 0.10566


def anti(name):
    return {"Pi0": "Pi0", "Pi+": "Pi-", "Pi-": "Pi+"}.get(name, name + "~")


def _mpi(model, pion):
    return model["MPi0" if pion == "Pi0" else "MPiP"]


def leptonic_current(k1, k2, model, mlep, nc=False, nubar=False):
    """L[mu, s_lep, s_nu] = ubar(k2) Gamma^mu u(k1) from the UFO lepton vertex.

    nc: the outgoing lepton is the neutrino and the boson the Z.  nubar (NC
    only): the antineutrino current vbar(k1) c gamma^mu P_L v(k2) equals
    ubar(k2) c gamma^mu P_R u(k1) up to a phase, and c gamma^mu P_R =
    g^{mu mu} gamma^0 (c gamma^mu P_L) gamma^0 in the Dirac representation."""
    if nc:
        Vl = build(model.V, {"vm~": -k1, "vm": k2, "Z": k1 - k2}, model.ns, model.ff)
        if nubar:
            g0 = np.diag([1.0, 1.0, -1.0, -1.0]).astype(complex)
            Vl = np.array([METRIC[m, m] * g0 @ Vl[m] @ g0 for m in range(4)])
    else:
        Vl = build(model.V, {"vm~": -k1, "mu-": k2, "W+": k1 - k2}, model.ns, model.ff)
    out = np.zeros((4, 2, 2), dtype=complex)
    for sl in range(2):
        ub = ubar(uspinor(k2, mlep, sl))
        for sn in range(2):
            un = uspinor(k1, 0.0, sn)
            for m in range(4):
                out[m, sl, sn] = ub @ Vl[m] @ un
    return out


def hadronic_current(ch, q, P, k, pp, model):
    """H[mu, row, col] for W+(q) N(P) -> pi(k) N'(p'), or Z(q) N(P) for NC."""
    c = CHANNELS[ch]
    WB = c.get("boson", "W-")            # an incoming W+ is the outgoing W- leg
    Ni, Nf, pion, pi_t = c["Ni"], c["Nf"], c["pion"], c["pi_t"]
    aNi = anti(Ni)
    V, ns, ff = model.V, model.ns, model.ff
    MN, MD = model["MNucleon"], model["MDelta"]
    vs, vu, vt = P + q, P - k, q - k
    H = np.zeros((4, 4, 4), dtype=complex)

    if c["N_s"]:
        N = c["N_s"]
        Vw = build(V, {aNi: -P, N: vs, WB: -q}, ns, ff)
        Vp = build(V, {anti(N): -vs, Nf: pp, pion: k}, ns, ff)
        S = prop_fermion(vs, MN)
        for m in range(4):
            H[m] += Vp @ S @ Vw[m]
    if c["N_u"]:
        N = c["N_u"]
        Vp = build(V, {aNi: -P, N: vu, pion: k}, ns, ff)
        Vw = build(V, {anti(N): -vu, Nf: pp, WB: -q}, ns, ff)
        S = prop_fermion(vu, MN)
        for m in range(4):
            H[m] += Vw[m] @ S @ Vp
    if pi_t:
        Vwpp = build_bosonic(V, {WB: -q, pion: k, pi_t: vt}, ns, ff, [WB])
        Vpi = build(V, {aNi: -P, Nf: pp, anti(pi_t): -vt}, ns, ff)
        D = prop_scalar(vt, _mpi(model, pi_t))
        for m in range(4):
            H[m] += Vwpp[m] * D * Vpi
    if c.get("contact", True):
        H += build(V, {aNi: -P, Nf: pp, WB: -q, pion: k}, ns, ff)
    # W+ -> pi+ gamma*/Z*, the boson absorbed by the nucleon (photon in Feynman
    # gauge, Z in unitary gauge)
    if Ni == Nf and WB == "W-":
        for B, MB, GB in (("A", 0.0, 0.0), ("Z", model["MZ"], model["WZ"])):
            Vwb = build_bosonic(V, {"W-": -q, pion: k, B: vt}, ns, ff, ["W-", B])
            VN = build(V, {aNi: -P, Nf: pp, B: -vt}, ns, ff)
            C = G[:, None] * prop_vector(vt, MB, GB) * G[None, :]
            for m in range(4):
                for a in range(4):
                    for b in range(4):
                        H[m] += Vwb[m, a] * C[a, b] * VN[b]

    s, u = dot(vs, vs), dot(vu, vu)
    Ds = s - MD ** 2 + 1j * MD * model.delta_width(s.real)
    Du = u - MD ** 2 + 1j * MD * model.delta_width(u.real)
    D = c["D_s"]
    VwD = build(V, {aNi: -P, D: vs, WB: -q}, ns, ff)                # [nu, mu, r, c]
    VpD = build(V, {anti(D): -vs, Nf: pp, pion: k}, ns, ff)         # [eta, r, c]
    H += delta_line(VpD, VwD, vs, MD, Ds)
    D = c["D_u"]
    VpD = build(V, {aNi: -P, D: vu, pion: k}, ns, ff)
    VwD = build(V, {anti(D): -vu, Nf: pp, WB: -q}, ns, ff)
    H += delta_line(VwD, VpD, vu, MD, Du)
    return H


def _boost(X, W, a):
    """Four-vector a, given in the rest frame of X (mass W), in X's frame."""
    Xv, av = X[1:], a[1:]
    xa = Xv @ av
    return np.concatenate([[(X[0] * a[0] + xa) / W],
                           av + Xv * (a[0] / W + xa / (W * (X[0] + W)))])


def _lepton_momentum(s, W, m):
    """|k2| of the outgoing lepton (mass m) in the overall CM, for hadronic mass W."""
    if m == 0.0:
        return (s - W ** 2) / (2 * np.sqrt(s))
    lam = (s - (m + W) ** 2) * (s - (m - W) ** 2)
    return np.sqrt(max(0.0, lam)) / (2 * np.sqrt(s))


def _msq(ch, s, W, ctl, cth, phi, model, mlep, nubar):
    """Spin-summed |M|^2 at one phase-space point, overall CM frame."""
    MN = model["MNucleon"]
    mpi = _mpi(model, CHANNELS[ch]["pion"])
    rs = np.sqrt(s)
    E1 = (s - MN ** 2) / (2 * rs)
    k1 = np.array([E1, 0, 0, E1], dtype=complex)
    P = np.array([rs - E1, 0, 0, -E1], dtype=complex)
    ql = _lepton_momentum(s, W, mlep)
    stl = np.sqrt(max(0.0, 1.0 - ctl ** 2))
    k2 = np.array([np.sqrt(ql ** 2 + mlep ** 2), ql * stl, 0, ql * ctl], dtype=complex)
    X = k1 + P - k2
    ks = kmom(W ** 2, mpi, MN)
    sth = np.sqrt(max(0.0, 1.0 - cth ** 2))
    krest = np.array([np.sqrt(ks ** 2 + mpi ** 2), ks * sth * np.cos(phi),
                      ks * sth * np.sin(phi), ks * cth])
    k = _boost(X.real, W, krest).astype(complex)
    pp = X - k

    nc = ch in NC_CHANNELS
    L = leptonic_current(k1, k2, model, mlep, nc, nubar)
    H = hadronic_current(ch, k1 - k2, P, k, pp, model)
    MB, GB = (model["MZ"], model["WZ"]) if nc else (model["MW"], model["WW"])
    C = G[:, None] * prop_vector(k1 - k2, MB, GB) * G[None, :]
    tot = 0.0
    for sp in range(2):
        ui = uspinor(P, MN, sp)
        for sq in range(2):
            ub = ubar(uspinor(pp, MN, sq))
            Hn = np.array([ub @ H[n] @ ui for n in range(4)])
            for sl in range(2):
                for sn in range(2):          # the wrong-helicity term is zero
                    tot += abs(L[:, sl, sn] @ C @ Hn) ** 2
    return tot


def sigma(ch, Enu, model, wmax=None, grid=(10, 6, 6, 5), nubar=False):
    """Total cross section in ub at neutrino energy Enu, nucleon at rest.

    wmax: upper cut on the pi N invariant mass.  grid: Gauss-Legendre nodes in
    (W, cos theta_mu, cos theta*_pi) and trapezoid nodes in phi*_pi.  nubar:
    antineutrino beam (NC channels only)."""
    if nubar and ch not in NC_CHANNELS:
        raise ValueError("nubar is only implemented for the NC channels")
    mlep = 0.0 if ch in NC_CHANNELS else MMU
    MN = model["MNucleon"]
    mpi = _mpi(model, CHANNELS[ch]["pion"])
    s = MN ** 2 + 2 * MN * Enu
    rs = np.sqrt(s)
    Wlo, Whi = MN + mpi, rs - mlep
    if wmax:
        Whi = min(Whi, float(wmax))
    if Whi <= Wlo:
        return 0.0
    nW, ncl, ncth, nphi = grid
    xW, wW = np.polynomial.legendre.leggauss(nW)
    xl, wl = np.polynomial.legendre.leggauss(ncl)
    xh, wh = np.polynomial.legendre.leggauss(ncth)
    xp = 2 * np.pi * (np.arange(nphi) + 0.5) / nphi
    wp = 2 * np.pi / nphi
    # W = Wlo + (Whi - Wlo) u^2 absorbs the sqrt(W - Wlo) threshold behaviour
    us = 0.5 * (xW + 1.0)
    tot = 0.0
    for wi, u in zip(wW, us):
        W = Wlo + (Whi - Wlo) * u ** 2
        pref = _lepton_momentum(s, W, mlep) * kmom(W ** 2, mpi, MN) / (128 * np.pi ** 4 * rs)
        if pref <= 0:
            continue
        acc = sum(wl[i] * wh[j] * wp * _msq(ch, s, W, xl[i], xh[j], phi, model, mlep, nubar)
                  for i in range(ncl) for j in range(ncth) for phi in xp)
        tot += wi * (Whi - Wlo) * u * pref * acc
    return HBARC2 * tot / (4 * (s - MN ** 2))
