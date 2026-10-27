"""t-channel rho/omega exchange in gamma N -> pi N (photoproduction).

Not part of the ChiEFT model: the couplings are the fit values of Fernandez-
Ramirez, Moya de Guerra and Udias, nucl-th/0509020, Eqs. (15)-(16).

    omega   G = 0.2804   F = 21.11          K = -0.17   m = 0.78257
    rho0    G = 0.1161   F = 2.6            K = 6.3     m = 0.7685
    rho+-   G = 0.0906   F = sqrt(2) 2.6    K = 6.3     m = 0.7665

pi0 p: omega + rho0;  pi0 n: omega + rho0 with F -> -F;  pi+ n, pi- p: rho+-.
Vertices: V pi gamma (i e G / m_pi) eps^{mu nu a b} q_a (-k)_b and
V N N -i F gamma^a + (i/4)(F K / M_N)[gamma^a, Pslash], with the form factor
Lambda^4/(Lambda^4 + (t - m_V^2)^2), Lambda = 1.121 GeV.
"""
import numpy as np

from .ufo import EPS, GAMMA, METRIC

EE = np.sqrt(4 * np.pi / 137.035999)
MN = 0.938272
MPI = 0.1349768
LAMBDA = 1.121
OMEGA = dict(G=0.2804, F=21.11, K=-0.17, m=0.78257)
RHO0 = dict(G=0.1161, F=2.6, K=6.3, m=0.7685)
RHOC = dict(G=0.0906, F=np.sqrt(2) * 2.6, K=6.3, m=0.7665)
EXCHANGES = {"ppi0": [OMEGA, RHO0], "npi0": [OMEGA, dict(RHO0, F=-2.6)],
             "npip": [RHOC], "ppim": [RHOC]}

_G = np.diag(METRIC).real
_GAMLOW = np.einsum("ab,bij->aij", np.diag(_G), GAMMA)


def current(ch, q, k):
    """M[mu, row, col] of the exchange, photon index up; q in, k (pion) out."""
    q, k = np.asarray(q).real, np.asarray(k).real
    vt = k - q
    t = vt @ np.diag(_G) @ vt
    sP = np.einsum("a,aij->ij", -vt, _GAMLOW)
    lc = np.einsum("mnab,a,b->mn", EPS.real, _G * q, _G * (-k))
    out = np.zeros((4, 4, 4), dtype=complex)
    for ex in EXCHANGES[ch]:
        dV = t - ex["m"] ** 2
        ffV = 1.0 / (1.0 + dV ** 2 / LAMBDA ** 4)
        vpg = (1j * EE * ex["G"] / MPI) * lc
        vnn = np.array([-1j * ex["F"] * GAMMA[a]
                        + 0.25j * (ex["F"] * ex["K"] / MN) * (GAMMA[a] @ sP - sP @ GAMMA[a])
                        for a in range(4)])
        for mu in range(4):
            acc = np.zeros((4, 4), dtype=complex)
            for al in range(4):
                acc += _G[al] * vpg[mu, al] * (-1j / dV) * vnn[al]
            out[mu] -= ffV * acc
    return out
