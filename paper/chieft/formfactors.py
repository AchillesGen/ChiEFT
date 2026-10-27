"""The form factors of sherpa/FormFactors.C, evaluated from the same param card.

Vector: Kelly fit (block CHIEFTKELLY).  Axial: z-expansion (CHIEFTZEXP, with
t_cut and t_0 in CHIEFTFF).  Most functions are ratios normalised to 1 at
Q^2 = 0, so the model keeps its own couplings (g_A, C5A = sqrt(2/3) h_D) and
only acquires the Q^2 shape.  The exceptions are the weak-magnetism terms
FW2, FZ2p, FZ2n, which carry absolute normalisation from CVC, and FZ1p, which
also moves sin^2(theta_W) from the model's on-shell value to the low-Q value
Sw2FF.
"""
import math

from .ufo import re_


def _zexpand(A, z):
    r = 0.0
    for a in A[:0:-1]:
        r = (r + a) * z
    return r + A[0]


def _z(Q2, tcut, t0):
    a, b = math.sqrt(tcut + Q2), math.sqrt(tcut - t0)
    return (a - b) / (a + b)


class FormFactors:
    """Nucleon and N -> Delta form factors from a UFO parameter namespace."""

    def __init__(self, ns):
        self.p = {k: re_(v) for k, v in ns.items() if isinstance(v, (int, float, complex))}
        self.MP2 = self.p["MNucleon"] ** 2
        self.sw2_model = re_(ns["sw"] ** 2)        # the model's own (on-shell) value
        self.sw2 = self.p["Sw2FF"]                 # the low-Q value
        self.FA0 = self.FA_raw(0.0)

    # vector
    def dirac_pauli(self, Q2):
        """(F1p, F1n, F2p, F2n) from the Kelly Sachs form factors."""
        p, tau = self.p, Q2 / 4 / self.MP2

        def kelly(a1, b1, b2, b3):
            return (1 + a1 * tau) / (1 + b1 * tau + b2 * tau ** 2 + b3 * tau ** 3)
        Gep = kelly(p["KEp1"], p["KEp2"], p["KEp3"], p["KEp4"])
        Gen = (1.0 / (1 + Q2 / p["KLamSq"]) ** 2
               * p["KEn1"] * tau / (1 + p["KEn2"] * tau))
        Gmp = p["MuP"] * kelly(p["KMp1"], p["KMp2"], p["KMp3"], p["KMp4"])
        Gmn = p["MuN"] * kelly(p["KMn1"], p["KMn2"], p["KMn3"], p["KMn4"])
        return ((Gep + tau * Gmp) / (1 + tau), (Gen + tau * Gmn) / (1 + tau),
                (Gmp - Gep) / (1 + tau), (Gmn - Gen) / (1 + tau))

    # axial
    def FA_raw(self, Q2):
        p = self.p
        return _zexpand([p["AxA%d" % i] for i in range(9)], _z(Q2, p["AxTcut"], p["AxT0"]))

    def FA_ratio(self, Q2):
        return self.FA_raw(Q2) / self.FA0

    def C5A_ratio(self, Q2):
        """C5A(Q^2)/C5A(0); CA5Norm > 0 replaces the model's C5A(0) = sqrt(2/3) h_D."""
        M2, norm = self.p["MADeltaSq"], self.p["CA5Norm"]
        scale = norm / (math.sqrt(2.0 / 3.0) * self.p["hD"]) if norm > 0.0 else 1.0
        return scale / (1 + Q2 / M2) ** 2 / (1 + Q2 / 3 / M2)

    def FrhoCT(self, Q2):
        m2 = self.p["MRhoCT"] ** 2
        return m2 / (m2 + Q2)

    # neutral-current charge combinations, J^Z = J^3 - 2 sw^2 J^EM
    def _FZ(self, Q2, pauli):
        F1p, F1n, F2p, F2n = self.dirac_pauli(Q2)
        a, b = (F2p, F2n) if pauli else (F1p, F1n)
        return 0.5 * (a - b) - 2 * self.sw2 * a, 0.5 * (b - a) - 2 * self.sw2 * b

    def namespace(self):
        """{name: function} for every form factor the UFO's Lorentz structures call."""
        mpic, mpi0 = self.p["MPiP"], self.p["MPi0"]
        L4, MD2 = self.p["LambdaFF"] ** 4, self.p["MDelta"] ** 2
        M2V = self.p["MV2ND"]
        F = {
            "F1p": lambda Q2: self.dirac_pauli(Q2)[0],
            "F1n": lambda Q2: self.dirac_pauli(Q2)[1],
            "F2p": lambda Q2: self.dirac_pauli(Q2)[2],
            "F2n": lambda Q2: self.dirac_pauli(Q2)[3],
            # charged current: CVC isovector combinations, axial ratio
            "FW1": lambda Q2: self.dirac_pauli(Q2)[0] - self.dirac_pauli(Q2)[1],
            "FW2": lambda Q2: self.dirac_pauli(Q2)[2] - self.dirac_pauli(Q2)[3],
            "FWA": self.FA_ratio,
            "FPP": lambda Q2: self.FA_ratio(Q2) / (Q2 + mpic ** 2),
            # neutral current
            "FZ1p": lambda Q2: self._FZ(Q2, False)[0] / (0.5 - 2 * self.sw2_model),
            "FZ1n": lambda Q2: self._FZ(Q2, False)[1] / (-0.5),
            "FZ2p": lambda Q2: self._FZ(Q2, True)[0],
            "FZ2n": lambda Q2: self._FZ(Q2, True)[1],
            "FZA": self.FA_ratio,
            "FPPZ": lambda Q2: self.FA_ratio(Q2) / (Q2 + mpi0 ** 2),
            "FZ1pi": lambda Q2: ((self.dirac_pauli(Q2)[0] - self.dirac_pauli(Q2)[1])
                                 * (1.0 - 2.0 * self.sw2) / (1.0 - 2.0 * self.sw2_model)),
            # contact term and pion pole
            "FrhoCT": self.FrhoCT,
            "FPPpi": lambda Q2: self.FrhoCT(Q2) / (Q2 + mpic ** 2),
            "FPPZpi": lambda Q2: self.FrhoCT(Q2) / (Q2 + mpi0 ** 2),
            # N -> Delta
            "FCA5": self.C5A_ratio,
            "FC6A": lambda Q2: self.C5A_ratio(Q2) / (Q2 + mpic ** 2),
            "FVND": lambda Q2: 1.0 / (1 + Q2 / M2V) ** 2 / (1 + Q2 / (4 * M2V)),
            # strong piNDelta vertex; its argument is p_Delta^2, not Q^2
            "FDeltaLam": lambda p2: L4 / (L4 + (p2 - MD2) ** 2),
        }
        return {k: (lambda f: lambda x: f(re_(x)))(f) for k, f in F.items()}
