SHELL := /bin/bash

# Conda environment name and runner. Override: make ENV=myenv CONDA=mamba ...
ENV    ?= gapmind
CONDA  ?= conda
RUN    := $(CONDA) run --no-capture-output -n $(ENV)

# Arguments for `make analyze`
FAA     ?=
NAME    ?=
SET     ?=
THREADS ?=

.PHONY: help env setup setup-code setup-hmms databases all analyze check clean distclean

help:
	@echo "GapMind end-to-end pipeline"
	@echo
	@echo "  make env         create the '$(ENV)' conda environment from environment.yml"
	@echo "  make setup       clone code, wire conda tools, install usearch + Swissknife, fetch HMM DBs"
	@echo "  make databases   build curated + steps databases (Phases 1-2) for the configured SETS"
	@echo "  make all         setup + databases (full from-scratch build)"
	@echo "  make analyze FAA=proteins.faa NAME=myorg [SET='aa carbon'] [THREADS=8]"
	@echo "  make check       syntax-check the scripts (bash -n)"
	@echo "  make distclean   delete the cloned code base and all built databases"
	@echo
	@echo "Typical run:  make env && conda activate $(ENV) && make all"
	@echo "Then:         make analyze FAA=my_proteins.faa NAME=myorg"

env:
	@if $(CONDA) env list | awk '{print $$1}' | grep -qx '$(ENV)'; then \
	  echo ">> Updating env $(ENV)"; $(CONDA) env update -n $(ENV) -f environment.yml; \
	else \
	  echo ">> Creating env $(ENV)"; $(CONDA) env create -n $(ENV) -f environment.yml; \
	fi

setup-code:
	CONDA_ENV=$(ENV) $(RUN) bash scripts/setup_code.sh

setup-hmms:
	CONDA_ENV=$(ENV) $(RUN) bash scripts/setup_hmms.sh

setup: setup-code setup-hmms

databases:
	CONDA_ENV=$(ENV) $(RUN) bash scripts/build_databases.sh

all: setup databases

analyze:
	@if [ -z "$(FAA)" ] || [ -z "$(NAME)" ]; then \
	  echo "Usage: make analyze FAA=proteins.faa NAME=myorg [SET='aa carbon'] [THREADS=8]"; exit 1; fi
	CONDA_ENV=$(ENV) $(RUN) bash scripts/run_gapmind.sh -f "$(FAA)" -n "$(NAME)" \
	  $(if $(SET),-s "$(SET)") $(if $(THREADS),-t "$(THREADS)")

check:
	@for s in scripts/*.sh; do echo "bash -n $$s"; bash -n "$$s"; done; echo "OK"

clean:
	@echo "Nothing deleted. Use 'make distclean' to remove the cloned code + databases."

distclean:
	@bash -c 'source ./config.sh && echo "Removing $$CODE_DIR" && rm -rf "$$CODE_DIR"'
