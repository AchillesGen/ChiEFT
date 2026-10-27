# pyNUISANCE for run_nuclear.py analyze: the achilles-physval image with its
# NUISANCE2 rebuilt at the same commit plus the fixes in patches/ (without them
# the MicroBooNE NCpi0, MicroBooNE CC1Mu1p/CC1MuNp and ElectronData samples
# return empty predictions).
#
# Build from paper/nuclear:
#   docker build -f docker/Dockerfile -t achilles-physval:chieft .
ARG BASE=ghcr.io/achillesgen/achilles-physval@sha256:5e1d35ea5c2f8f9af9b1c62cdfd04412a74ca4b7c8f7706b70f0e954f369aa60
FROM ${BASE}

ARG NUISANCE2_REPO=https://github.com/NUISANCEMC/nuisance.git
ARG NUISANCE2_SHA=f5f79db6c7b80b496e249f938982d43291f91de7
ARG JOBS=8

COPY patches/ /tmp/patches/

RUN set -eux; \
    git clone "${NUISANCE2_REPO}" /tmp/nuisance2; \
    git -C /tmp/nuisance2 checkout "${NUISANCE2_SHA}"; \
    for p in /tmp/patches/*.patch; do git -C /tmp/nuisance2 apply "$p"; done; \
    cmake -S /tmp/nuisance2 -B /tmp/n2build -DCMAKE_INSTALL_PREFIX=/opt/nuisance2 \
          -DCMAKE_BUILD_TYPE=Release -DNuHepMC_ENABLED=ON; \
    cmake --build /tmp/n2build -j"${JOBS}"; \
    cmake --install /tmp/n2build; \
    rm -rf /tmp/nuisance2 /tmp/n2build /tmp/patches

RUN python3 -c "import pyNUISANCE as pn; r = pn.RecordFactory().make_record({'type': 'nuisance2'}); r.analysis('MicroBooNE_NCpi0_XSec_XpPpi0_nu'); r.analysis('ElectronData_6_12_0.560_36.000_Barreau:1983ht')"
