"""gamma(*) N -> pi N from the UFO: hadronic current, photoproduction cross
section and the virtual-photon responses used to reweight electron events.

Diagrams: s- and u-channel nucleon, t-channel pion and Kroll-Ruderman contact
(charged pions), s- and u-channel Delta.  For real photons, optionally the
t-channel rho/omega exchange of vecmeson.py, which is not part of the model.
"""
import numpy as np

from . import vecmeson
from .ufo import (HBARC2, METRIC, POL, build, build_bosonic, delta_line, dot,
                  kmom, prop_fermion, prop_scalar, ubar, uspinor)

# gamma N_i -> N_f pi, with the Delta charge states of the s and u channels
CHANNELS = {
    "ppi0": dict(Ni="p", Nf="p", pion="Pi0", D_s="Delta+", D_u="Delta+"),
    "npip": dict(Ni="p", Nf="n", pion="Pi+", D_s="Delta+", D_u="Delta0"),
    "npi0": dict(Ni="n", Nf="n", pion="Pi0", D_s="Delta0", D_u="Delta0"),
    "ppim": dict(Ni="n", Nf="p", pion="Pi-", D_s="Delta0", D_u="Delta+"),
}
LABEL = {"ppi0": r"$\gamma p \to p\pi^0$", "npip": r"$\gamma p \to n\pi^+$",
         "npi0": r"$\gamma n \to n\pi^0$", "ppim": r"$\gamma n \to p\pi^-$"}


def anti(name):
    return {"Pi0": "Pi0", "Pi+": "Pi-", "Pi-": "Pi+"}.get(name, name + "~")


def pion_mass(model, ch):
    return model["MPi0" if CHANNELS[ch]["pion"] == "Pi0" else "MPiP"]


def current(ch, q, P, k, pp, model, born=True, delta=True, vec=False):
    """H[mu, row, col] for gamma*(q) N(P) -> pi(k) N(p'), photon index up."""
    c = CHANNELS[ch]
    Ni, Nf, pion = c["Ni"], c["Nf"], c["pion"]
    aNi, aNf = anti(Ni), anti(Nf)
    V, ns, ff = model.V, model.ns, model.ff
    MN, MD = model["MNucleon"], model["MDelta"]
    vs, vu, vt = P + q, P - k, q - k
    H = np.zeros((4, 4, 4), dtype=complex)

    if born:
        Vg = build(V, {aNi: -P, Ni: vs, "A": -q}, ns, ff)
        Vp = build(V, {aNi: -vs, Nf: pp, pion: k}, ns, ff)
        S = prop_fermion(vs, MN)
        for m in range(4):
            H[m] += Vp @ S @ Vg[m]
        Vp = build(V, {aNi: -P, Nf: vu, pion: k}, ns, ff)
        Vg = build(V, {aNf: -vu, Nf: pp, "A": -q}, ns, ff)
        S = prop_fermion(vu, MN)
        for m in range(4):
            H[m] += Vg[m] @ S @ Vp
        if pion != "Pi0":
            Vgpp = build_bosonic(V, {"A": -q, pion: k, anti(pion): vt}, ns, ff, ["A"])
            Vpi = build(V, {aNi: -P, Nf: pp, pion: -vt}, ns, ff)
            D = prop_scalar(vt, pion_mass(model, ch))
            for m in range(4):
                H[m] += Vgpp[m] * D * Vpi
            H += build(V, {aNi: -P, Nf: pp, "A": -q, pion: k}, ns, ff)

    if delta:
        s, u = dot(vs, vs), dot(vu, vu)
        Ds = s - MD ** 2 + 1j * MD * model.delta_width(s.real)
        Du = u - MD ** 2 + 1j * MD * model.delta_width(u.real)
        Dsn, Dun = c["D_s"], c["D_u"]
        VgD = build(V, {aNi: -P, Dsn: vs, "A": -q}, ns, ff)            # [nu, mu, r, c]
        VpD = build(V, {anti(Dsn): -vs, Nf: pp, pion: k}, ns, ff)      # [eta, r, c]
        H += delta_line(VpD, VgD, vs, MD, Ds)
        VpD = build(V, {aNi: -P, Dun: vu, pion: k}, ns, ff)
        VgD = build(V, {anti(Dun): -vu, Nf: pp, "A": -q}, ns, ff)
        H += delta_line(VgD, VpD, vu, MD, Du)

    if vec:
        # Relative sign fixed against the ELA amplitude of nucl-th/0509020, which
        # vecmeson.py's couplings are fitted with, in this UFO's phase convention.
        H += vecmeson.current(ch, q, k)
    return H


def _cm_kinematics(model, ch, W, Q2):
    """(q, P, |k*|) in the hadronic CM, photon along +z; None below threshold."""
    MN, mpi = model["MNucleon"], pion_mass(model, ch)
    if W <= MN + mpi:
        return None
    nu = (W ** 2 - MN ** 2 - Q2) / (2 * W)
    qa = np.sqrt(nu ** 2 + Q2)
    q = np.array([nu, 0, 0, qa], dtype=complex)
    P = np.array([W - nu, 0, 0, -qa], dtype=complex)
    return q, P, kmom(W ** 2, mpi, MN)


def _angles(model, ch, W, ks, nth):
    """Pion momenta on a Gauss-Legendre grid in cos(theta*), with weights."""
    Epi = np.sqrt(ks ** 2 + pion_mass(model, ch) ** 2)
    for x, w in zip(*np.polynomial.legendre.leggauss(nth)):
        yield np.array([Epi, ks * np.sqrt(1 - x ** 2), 0, ks * x], dtype=complex), w


def _spin_currents(H, P, pp, MN):
    """ubar(p') H^mu u(P) for the four nucleon spin combinations."""
    for si in range(2):
        ui = uspinor(P, MN, si)
        for sf in range(2):
            ub = ubar(uspinor(pp, MN, sf))
            yield np.array([ub @ H[m] @ ui for m in range(4)])


def _transverse(J):
    return sum(abs(sum(METRIC[m, m] * e[m] * J[m] for m in range(4))) ** 2 for e in POL)


def sigma(ch, Egamma, model, nth=16, **opts):
    """Total photoproduction cross section in ub at lab photon energy Egamma.

    opts go to current(): born, delta, vec."""
    MN = model["MNucleon"]
    W = np.sqrt(MN ** 2 + 2 * MN * Egamma)
    kin = _cm_kinematics(model, ch, W, 0.0)
    if kin is None:
        return 0.0
    q, P, ks = kin
    tot = 0.0
    for k, w in _angles(model, ch, W, ks, nth):
        pp = P + q - k
        H = current(ch, q, P, k, pp, model, **opts)
        tot += w * sum(_transverse(J) for J in _spin_currents(H, P, pp, MN)) / 4.0
    qs = q[0].real
    return float((2 * np.pi * HBARC2 * ks / (64 * np.pi ** 2 * W ** 2 * qs) * tot).real)


def responses(ch, W, Q2, model, nth=16):
    """Angle-integrated responses (S_T, S_L) times k*/W at (W, Q^2).

    sigma_T + eps sigma_L is proportional to S_T + eps S_L with a factor that
    does not depend on the Delta, so their ratio between two models is the
    event weight at (W, Q^2, eps)."""
    kin = _cm_kinematics(model, ch, W, Q2)
    if kin is None:
        return 0.0, 0.0
    q, P, ks = kin
    MN, qa = model["MNucleon"], q[3].real
    ST = SL = 0.0
    for k, w in _angles(model, ch, W, ks, nth):
        pp = P + q - k
        for J in _spin_currents(current(ch, q, P, k, pp, model), P, pp, MN):
            ST += w * 0.5 * _transverse(J)
            SL += w * (Q2 / qa ** 2) * abs(J[0]) ** 2
    return ST * ks / W, SL * ks / W


def breit_wigner(Egamma, model):
    """Unitary Breit-Wigner sigma(gamma p -> Delta+ -> p pi0) in ub, with
    Gamma(N pi) ~ k^3 and Gamma(N gamma) ~ q^3 at the model's pole parameters."""
    MN, MD, mpi = model["MNucleon"], model["MDelta"], model["MPi0"]
    W = np.sqrt(MN ** 2 + 2 * MN * Egamma)
    if W <= MN + mpi:
        return 0.0
    q, q0 = (W ** 2 - MN ** 2) / (2 * W), (MD ** 2 - MN ** 2) / (2 * MD)
    k, k0 = kmom(W ** 2, mpi, MN), kmom(MD ** 2, mpi, MN)
    Gpi, Gg = model["GamNPi"] * (k / k0) ** 3, model["GamNGam"] * (q / q0) ** 3
    spin = 4.0 / 4.0                   # (2J + 1) / ((2 s_gamma) (2 s_N + 1)), J = 3/2
    iso = 2.0 / 3.0                    # Delta+ -> p pi0
    return float(4 * np.pi / q ** 2 * spin * (Gg * Gpi / 4)
                 / ((W - MD) ** 2 + (Gpi + Gg) ** 2 / 4) * iso * HBARC2)
