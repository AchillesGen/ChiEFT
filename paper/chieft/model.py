"""The ChiEFT model as the semi-analytic code evaluates it.

Two parametrisations of the Delta(1232):

    pole   the param card as is: M = 1.209 GeV, Gamma = 102 MeV (the pole
           position), constant width in the s- and u-channel propagators.
           This is what Sherpa computes.
    bw     the Breit-Wigner parameters M = 1.232 GeV, Gamma(N pi) = 117 MeV,
           Gamma(N gamma) = 0.66 MeV.  h_D, G1 and G2 are internal parameters
           derived from these, so they follow.

and, independently, a running width

    Gamma(s) = Gamma_0 (M/sqrt(s)) (k(s)/k(M^2))^3,  0 below the pi N threshold

(so 0 in the u channel), which Sherpa cannot do.  The running-width variant
of the paper is `bw` with the running width.
"""
from dataclasses import dataclass

import numpy as np

from . import PARAM_CARD, UFO
from .formfactors import FormFactors
from .ufo import kmom, load_ufo, parameter_namespace, re_

BREIT_WIGNER = {"MDelta": 1.232, "GamNPi": 0.117, "GamNGam": 0.00066}


@dataclass
class Model:
    V: object          # the UFO vertices module
    ns: dict           # parameters, resolved
    ff: dict           # form-factor functions
    running: bool      # running Delta width instead of the constant one

    def __getitem__(self, name):
        return re_(self.ns[name])

    def delta_width(self, s):
        """Delta width at invariant mass squared s."""
        G0 = self["WDeltaP"]
        if not self.running:
            return G0
        MN, MD, mpi = self["MNucleon"], self["MDelta"], self["MPi0"]
        if s < (MN + mpi) ** 2:
            return 0.0
        return G0 * (MD / np.sqrt(s)) * (kmom(s, mpi, MN) / kmom(MD ** 2, mpi, MN)) ** 3


def load(delta="pole", running=None, ufo=UFO, param_card=PARAM_CARD, overrides=None):
    """Model with the `delta` parametrisation ('pole' or 'bw').

    running defaults to True for 'bw' and False for 'pole'.  overrides sets
    further external parameters by name, e.g. {'LambdaFF': 1.0}.
    """
    if delta not in ("pole", "bw"):
        raise ValueError("delta must be 'pole' or 'bw', got %r" % delta)
    P, V, L = load_ufo(ufo)
    over = dict(BREIT_WIGNER) if delta == "bw" else {}
    over.update(overrides or {})
    ns = parameter_namespace(P, param_card, over)
    # The width in the propagator is the sum of the partial widths that fix the
    # couplings, so vertex and propagator describe the same resonance.
    ns["WDeltaP"] = re_(ns["GamNPi"]) + re_(ns["GamNGam"])
    ff = FormFactors(ns).namespace()
    return Model(V, ns, ff, delta == "bw" if running is None else running)
