"""Semi-analytic evaluation of the ChiPT UFO, for the paper's figures.

Every vertex, coupling and parameter is read from the UFO (and its param card);
this package only adds what a UFO does not contain: propagators, spinors and
phase space.

    ufo          UFO loading, Dirac algebra and the Lorentz-structure evaluator
    formfactors  the form factors of sherpa/FormFactors.C, in Python
    model        load(): one UFO + param card, constant or running Delta width
    photo        gamma(*) N -> pi N: hadronic current, sigma, virtual-photon responses
    weak         nu_mu N -> mu- N' pi: sigma
    vecmeson     t-channel rho/omega exchange (not part of the model)
    style        the paper's figure style
"""
import os

PAPER = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(PAPER)
UFO = os.path.join(REPO, "ChiPT_UFO")
PARAM_CARD = os.path.join(REPO, "param_card.dat")
