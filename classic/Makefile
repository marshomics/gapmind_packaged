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
# Argument for `make batch-prepare` (or set MANIFEST in config.sh)
MANIFEST ?=
# Parallelism for `make batch-local`: how many batches to run at once
JOBS    ?= 1

.PHONY: help env setup setup-code setup-hmms databases all analyze check clean distclean \
        batch-prepare batch-calibrate batch-local batch-submit batch-status batch-merge plots presence

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
	@echo "Large-scale batch analysis (many genomes), on a cluster or locally:"
	@echo "  make batch-prepare MANIFEST=genomes.tsv   split a TSV of .faa paths into batches"
	@echo "  make batch-calibrate                      time one batch to size batch/mem/walltime"
	@echo "  make batch-local [JOBS=N]                 run every batch on THIS machine, no scheduler"
	@echo "  make batch-submit                         SGE only: qsub the array + a held merge job"
	@echo "  make batch-status                         progress (qstat + completed count)"
	@echo "  make batch-merge                          combine finished batches into merged tables"
	@echo "  make plots                                summary/QC/biological figures (PNG+SVG) from merged output"
	@echo "  make presence [PA_MODE=probably|strict]   genome x pathway 0/1 presence + 2/1/0 confidence matrices"
	@echo
	@echo "Typical run:       make env && conda activate $(ENV) && make all"
	@echo "One proteome:      make analyze FAA=my_proteins.faa NAME=myorg"
	@echo "Many, no cluster:  make batch-prepare MANIFEST=genomes.tsv && make batch-local JOBS=4 && make batch-merge && make plots && make presence"

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

batch-prepare:
	CONDA_ENV=$(ENV) $(RUN) bash scripts/batch_prepare.sh $(MANIFEST)

batch-calibrate:
	@echo ">> Calibrating on batch 1 (real output kept). Tune BATCH_SIZE in config.sh and rerun to retarget walltime."
	CONDA_ENV=$(ENV) $(RUN) bash -c 'if command -v /usr/bin/time >/dev/null 2>&1; then /usr/bin/time -v bash scripts/batch_run.sh 1; else time bash scripts/batch_run.sh 1; fi'
	@echo ">> Set SGE_H_RT from 'Elapsed (wall clock) time' and SGE_H_VMEM from 'Maximum resident set size' (KB) in config.sh."

batch-submit:
	bash scripts/batch_submit.sh

# Run all prepared batches on the local machine (no SGE). JOBS = how many batches
# at once; each still uses SGE_SLOTS cores (default 8), so keep JOBS x SGE_SLOTS
# near your core count, e.g. `SGE_SLOTS=4 make batch-local JOBS=4` on a 16-core box.
# Idempotent: finished batches (.done) are skipped, so re-run to resume after a failure.
batch-local:
	@CONDA_ENV=$(ENV) $(RUN) bash -c 'set -euo pipefail; source ./config.sh; nbf="$$BATCH_DIR/nbatches.txt"; [ -s "$$nbf" ] || { echo "No batches yet. Run: make batch-prepare MANIFEST=your.tsv" >&2; exit 1; }; n=$$(cat "$$nbf"); j=$(JOBS); echo ">> running $$n batches locally, $$j at a time, $${SGE_SLOTS:-8} threads each"; if [ "$$j" -gt 1 ]; then seq "$$n" | xargs -P "$$j" -I{} bash scripts/batch_run.sh {}; else for i in $$(seq 1 "$$n"); do bash scripts/batch_run.sh "$$i"; done; fi; echo ">> all batches done. Next: make batch-merge && make plots && make presence"'

batch-status:
	-@qstat -u "$$USER" 2>/dev/null || echo "(qstat unavailable here)"
	@CONDA_ENV=$(ENV) $(RUN) bash -c 'source ./config.sh; nb=$$(cat "$$BATCH_DIR/nbatches.txt" 2>/dev/null || echo 0); d=$$(find "$$BATCH_DIR/results" -maxdepth 2 -name .done 2>/dev/null | wc -l); echo "$$d / $$nb batches complete"'

batch-merge:
	CONDA_ENV=$(ENV) $(RUN) bash scripts/batch_merge.sh

plots:
	CONDA_ENV=$(ENV) $(RUN) bash scripts/plots.sh

presence:
	CONDA_ENV=$(ENV) $(RUN) bash scripts/presence.sh

check:
	@for s in scripts/*.sh; do echo "bash -n $$s"; bash -n "$$s"; done; echo "OK"

clean:
	@echo "Nothing deleted. Use 'make distclean' to remove the cloned code + databases."

distclean:
	@bash -c 'source ./config.sh && echo "Removing $$CODE_DIR" && rm -rf "$$CODE_DIR"'
