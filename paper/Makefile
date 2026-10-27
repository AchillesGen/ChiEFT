# Reproduce the figures of the ChiEFT paper.  See README.md for requirements.
#
#   make tables          semi-analytic tables (Python only, ~1 h on 12 cores)
#   make sherpa          free-nucleon Sherpa scans (needs Sherpa + the ChiEFT_UFO model)
#   make events          Achilles event samples (needs $ACHILLES_BUILD; many CPU hours)
#   make analysis        NUISANCE binning of the events (needs docker)
#   make figures         every paper figure from whatever of the above exists
#   make all-samples     every NUISANCE sample of analyses.yml, also those not in the paper

PY ?= python3

.PHONY: figures tables sherpa events analysis all-samples docker-image

figures:
	cd free_nucleon && $(PY) figures.py
	cd nuclear && $(PY) figures.py

tables:
	cd free_nucleon && $(PY) tables.py photo && $(PY) tables.py cc && $(PY) tables.py nc && $(PY) tables.py electron

sherpa:
	cd free_nucleon && $(PY) sherpa_scan.py photo && $(PY) sherpa_scan.py cc

events:
	cd nuclear && $(PY) run_nuclear.py generate

docker-image:
	cd nuclear && docker build -f docker/Dockerfile -t achilles-physval:chieft .

analysis: docker-image
	cd nuclear && $(PY) run_nuclear.py analyze --docker

all-samples:
	cd nuclear && $(PY) all_samples.py
