# GapMind — Nextflow + Singularity

A portable, containerized version of the GapMind batch pipeline. Same analysis as
the SGE scripts in the parent folder (amino-acid biosynthesis + carbon catabolism,
per-genome pathway calls, merged tables, presence/absence matrices, and figures),
but orchestrated by Nextflow and run inside a container — so it runs on **local,
SGE, and SLURM** with a profile switch, and the host environment no longer matters.

This is self-contained and separate from the SGE pipeline; nothing here touches it.

## Why this exists

Every environment problem the SGE version had to patch at setup time (system-perl
shebangs, Swissknife's path, dash vs bash, missing Perl modules) came from the host
differing from what the scripts assume. The container fixes that once: it carries a
pinned conda environment and the patched PaperBLAST code, so a task runs identically
on any machine. Nextflow handles scheduling, parallelism, retries, and resume, and
abstracts SGE/SLURM/local behind `-profile`.

It targets the **prebuilt-database** path only (the recommended one): the curated +
steps databases are downloaded ready-made, so there is no from-scratch build, no
Pfam/Swiss-Prot download, and no usearch/BLAST-legacy in the analysis itself.

## Requirements

- **Nextflow** ≥ 22.10 (needs Java; recent Nextflow needs Java 17+).
- A container runtime: **Singularity/Apptainer** (HPC) or **Docker** (workstation).

## 1. Build the container (once)

```bash
cd containers
./build.sh gapmind:latest <PAPERBLAST_COMMIT_SHA>   # pin the SHA for reproducibility
```

`build.sh` builds the Docker image and, if `apptainer`/`singularity` is present,
converts it to `gapmind.sif`. On a cluster without Docker, build the image on a
workstation and `apptainer build gapmind.sif docker://<registry>/gapmind:latest`,
or `apptainer build` directly from the Dockerfile-derived image. Point runs at it
with `--singularity_image /path/to/gapmind.sif`.

## 2. Samplesheet

A CSV with a header and a `faa` column (path to a protein FASTA); `sample` is
optional and derived from the filename if absent. Same inputs as the SGE manifest.

```
sample,faa
Dvulgaris_Hildenborough,/data/genomes/GCF_000195755.1.faa
Ecoli_K12_MG1655,/data/genomes/GCF_000005845.2.faa
```

## 3. Run

```bash
# SGE + Singularity
nextflow run . --input samplesheet.csv --outdir results \
         -profile sge,singularity --singularity_image containers/gapmind.sif

# SLURM + Singularity
nextflow run . --input samplesheet.csv --outdir results -profile slurm,singularity

# local + Docker (small runs / a workstation)
nextflow run . --input samplesheet.csv --outdir results -profile local,docker
```

Nextflow batches the genomes (`--batch_size`, default 500), runs one containerized
task per batch on the chosen executor, merges the results, and builds the
presence/absence matrices and figures. Failed tasks are retried and the run
resumes with `-resume` — no `.done` bookkeeping needed.

Key options (see `nextflow.config`): `--sets 'aa carbon'`, `--batch_size`,
`--db_dir` (reuse an existing database directory instead of downloading — e.g. the
`tmp/path.*` your SGE build already made), `--pa_mode probably|strict`,
`--keep_cand true|false`, `--knownsim true|false`, `--search_tool diamond|usearch`.

## Outputs (`--outdir`)

- `merged/` — `<set>.sum.rules`, `<set>.sum.steps` (and `.cand`, `.warn`,
  `.knownsim`), plus `orgs.tsv`, each row labelled by your sample name.
- `presence/` — `<set>.presence.tsv` (0/1), `<set>.confidence.tsv` (2/1/0),
  `<set>.pathways.tsv` legend.
- `plots/` — PNG + SVG figures and `summary_stats.tsv`.
- `pipeline_info/` — timeline, report, and trace (provenance for a shared resource).

## Validate the wiring

```bash
nextflow run . -stub -profile test
```

`-stub` runs the whole DAG with placeholder commands — no container or real data —
and publishes the expected `merged/`, `presence/`, and `plots/` structure. Good for
a first check and for CI. (This repo was validated this way.)

## Layout

```
main.nf                     workflow: PREPARE_DB -> GAPMIND_BATCH -> MERGE -> PRESENCE + PLOTS
nextflow.config             params + profiles (local/sge/slurm, singularity/apptainer/docker)
conf/test.config            tiny stub profile
bin/                        task scripts (auto on PATH inside processes)
  prepare_db.sh             download prebuilt DBs + extractHmms + diamond makedb
  run_batch.sh              analyze one batch (buildorgs -> gapsearch -> ... )
  merge_tables.sh           combine per-batch tables, join sample names
  presence_absence.py       genome x pathway 0/1 + 2/1/0 matrices
  make_plots.py             summary / QC / biological figures
containers/
  Dockerfile                miniforge + pinned env + patched PaperBLAST
  environment.yml           pinned conda environment
  setup_paperblast.sh       clone (pinned) + shebang/tool patches + usearch
  build.sh                  docker build (+ apptainer .sif)
  lock.sh                   pin environment.yml -> conda-lock.yml
assets/samplesheet.csv      example
nextflow_schema.json        typed, documented parameter schema
(repo root)                 LICENSE, NOTICE, CITATION.cff, .zenodo.json and
                            .github/workflows/ci.yml apply to the whole repository
```

## Notes

- **Resources.** `process_medium` (the batch tasks) defaults to 8 cpus / 24 GB /
  12 h, matching the SGE calibration. On SGE, `cpus` maps to `-pe parallel N` and
  `-q standard.q` via the `sge` profile; if your site doesn't honor `memory`/`time`
  automatically, uncomment the `clusterOptions` line in `nextflow.config` and set
  `h_vmem`/`h_rt` there.
- **Databases.** `PREPARE_DB` downloads and formats the prebuilt DBs once and caches
  them in the run's work dir; pass `--db_dir` to reuse a directory across runs.
- **ARM / license-clean image.** usearch (x86-only, non-free) is only pulled for the
  amino-acid known-gap step; drop that line in `setup_paperblast.sh` and run with
  `--knownsim false` to build a fully free, ARM-buildable image.
- **Pin PaperBLAST.** Build with a commit SHA (`build.sh <sha>`), not `master`, so
  the image is reproducible.

## Reproducibility

Three things pin a run: the PaperBLAST commit (`build.sh <sha>`), the conda
environment, and the database snapshot. To pin the environment exactly, run
`containers/lock.sh` on a machine that has conda (this repo's CI or a workstation —
the file `environment.yml` gives the intended versions; the lock gives the exact
build hashes). It writes `conda-lock.yml`; the header of `lock.sh` shows the
two-line Dockerfile change to install from the lock instead of solving fresh.

Parameters are described and type-checked in `nextflow_schema.json` (enums for
`search_tool`/`pa_mode`, `input`/`outdir` required). It doubles as documentation
and drives nf-core tooling if you later adopt it.

## Continuous integration

The repository's CI (`.github/workflows/ci.yml`, at the repo root) runs on every
push and PR: a stub run of this workflow (`nextflow run . -stub -profile test`, the
same validation described above, no container or data) plus linting (ShellCheck and
Ruff) across both pipelines. Green CI means the DAG still wires up and the scripts
still parse; it does not run the science.

## License and citation

GPL-3.0-or-later ([`../LICENSE`](../LICENSE)), matching PaperBLAST/GapMind, which
this workflow wraps rather than re-implements. [`../NOTICE`](../NOTICE) attributes
the wrapped software and the third-party tools in the image, and flags the two
redistribution items to confirm before you publish a public container (usearch, and
any archived database snapshot). If you use this, cite the GapMind papers (Price et
al., mSystems 2020 and PLOS Genetics 2022) as listed in
[`../CITATION.cff`](../CITATION.cff).

The license, citation, and CI files live at the repo root and cover both pipelines.
Before publishing, fill the placeholders in `../CITATION.cff`, `../.zenodo.json`,
`../NOTICE`, and this folder's `nextflow_schema.json` `$id`: author surname,
affiliation, ORCID, and the repository URL/owner.
