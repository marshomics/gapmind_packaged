# GapMind command-line pipeline

A conda-based, scripted version of the GapMind install from
[`gaps/SETUP`](https://github.com/morgannprice/PaperBLAST/blob/master/gaps/SETUP).
It takes you from nothing to "analyze a proteome for amino-acid biosynthesis and
carbon catabolism pathways" with one environment file and four scripts. By default
it downloads the prebuilt curated + steps databases; rebuilding those databases from scratch is an
opt-in (`DB_SOURCE=build`).

## Quickstart

```bash
# 1. create and activate the environment
make env
conda activate gapmind

# 2. set up code + tools, then get the databases (prebuilt download by default)
make all

# 3. analyze your proteome (a protein FASTA) for both pathway sets
make analyze FAA=/path/to/my_proteins.faa NAME=myorg
```

`make all` runs `setup` then `databases`; by default `databases` downloads the
prebuilt curated + steps databases. To rebuild them from scratch instead, use
`DB_SOURCE=build make databases` (see [How the databases are produced](#how-the-databases-are-produced)).

Outputs land in `PaperBLAST/tmp/myorg/`. Start with `aa.sum.rules` and
`carbon.sum.rules` (one row per pathway; `nHi`/`nMed`/`nLo` count high-, medium-,
and low-confidence steps), then `*.sum.steps` for the best candidate per step, and
`*.sum.db` for the combined sqlite3 database.

## How the databases are produced

In the default **prebuilt** mode, `make databases` downloads `curated.faa`,
`curated.db`, and `steps.db` for each set, extracts the HMM models that `steps.db`
references, and formats the curated DB for diamond. Phases 1 and 2 below are
skipped, and the databases match the official GapMind exactly.

In **build** mode (`DB_SOURCE=build`) those databases are rebuilt from scratch.
This is faithful to the sources but fragile: the pathway step files name curated
IDs from an older data snapshot, so a current-data rebuild can fail on UniProt or
BRENDA identifiers that have since changed. Mapped to the SETUP document:

Phase 1, the curated database, is `setupGaps.pl`. It pulls characterized proteins
out of the PaperBLAST database, adds curated Swiss-Prot sequences, clusters them
with usearch, runs Pfam hits, and writes `curated.faa` / `curated.db` /
`curated2.faa` into `tmp/path.<set>/`. This is the step that needs the PaperBLAST
data download, Swiss-Prot, usearch, and the Pfam HMMs.

Phase 2, the steps database, is `buildStepsDb.pl -doquery`. It runs `gapquery.pl`
for each pathway (matching curated proteins, Pfam/TIGRFam models, and some UniProt
sequences fetched live) and writes `steps.db`. `extractHmms.pl` then unpacks the
referenced models, and the curated DB is formatted for diamond (and BLAST).

Phase 3, the analysis, is `run_gapmind.sh`. For each proteome it runs
`buildorgs` → `gapsearch` → `gaprevsearch` → `gapsummary` →
`checkGapRequirements` → `orgsVsMarkers` → `buildGapsDb`, producing the summary
tables and sqlite database described above.

If you only ever want Phase 3 on new genomes after one build, just rerun
`make analyze` with different `FAA` and `NAME` values. The databases are built
once and reused.

## Configuration

Edit `config.sh`. The values you're most likely to change:

- `SETS` — which pathway groups to build and analyze (`aa`, `carbon`, or both).
- `SEARCH_TOOL` — `diamond` (default) or `usearch` for Phase 3.
- `THREADS` — CPU threads for searches.
- `CODE_DIR` — where the code base is cloned (default: `./PaperBLAST`).
- the `*_URL` / `*_BASE` variables — point these at a mirror or a pinned release if
  a download location changes.

Every variable also honors an environment override, for example
`SETS="aa" THREADS=16 make all`.

## Running steps individually

```bash
make setup-code     # clone + tool symlinks + usearch + Swissknife
make setup-hmms     # Pfam + TIGRFAMs download, index, metadata tables
make databases      # Phase 1 + Phase 2 for each set in $SETS
make analyze FAA=... NAME=... SET=aa THREADS=8
make check          # bash -n syntax check of every script
make distclean      # delete the clone and all built databases
```

## Large-scale analysis on an SGE cluster

The batch mode takes a tab-separated manifest, splits
the genomes into batches, runs one SGE array task per batch, and merges the
results. It's resumable and the batch size only affects speed, not the calls.

The manifest is a TSV with a header. One column must be `faa` (path to a protein
FASTA); an optional `name` column gives each genome a unique id (derived from the
filename if absent):

```
faa	name
/data/genomes/GCF_000195755.1.faa	Dvulgaris_H
/data/genomes/GCF_902167245.1.faa	Dvulgaris_K
...
```

Workflow:

```bash
# 1. split into batches (validates files, sanitizes names, stages messy paths)
make batch-prepare MANIFEST=genomes.tsv

# 2. time one batch to size things (keeps its real output)
make batch-calibrate

# 3. set SGE_* and BATCH_SIZE in config.sh from what you learned, then submit
make batch-submit          # array job + a merge job held until the array finishes

# 4. watch it; merge is automatic, or run it yourself for partial results
make batch-status
make batch-merge
```

Merged tables land in `batch/merged/<set>.sum.rules` (and `.steps`, `.warn`, and
`.cand` if `KEEP_CAND=1`), each row prefixed with the genome `name` from your
manifest. Start with `<set>.sum.rules` — one row per genome per pathway, with the
`nHi`/`nMed`/`nLo` counts.

### Sizing the job

Set these in `config.sh`:

- `SGE_SLOTS` — threads per task. 4-8 is the useful range; GapMind's per-batch HMM
  sweep and diamond search parallelize across these, with diminishing returns past
  ~8 for a single batch.
- `SGE_TC` — how many array tasks run at once. Keep `SGE_TC × SGE_SLOTS` at or below
  the slots you're allowed to use concurrently. If your budget is 400 slots and
  `SGE_SLOTS=4`, set `SGE_TC=100`.
- `BATCH_SIZE` — genomes per task. Use `make batch-calibrate` to time one batch, then
  pick a size that lands each task around 1-3 hours: comfortably inside your queue's
  runtime limit, short enough that a failure is cheap to retry, long enough that the
  per-task startup (loading the databases, forking the HMM searches) is amortized. If
  a batch of 300 takes 20 minutes, ~1000 gives roughly hourly tasks and fewer of them.
- `SGE_H_VMEM` / `SGE_H_RT` — from the calibration run's peak memory and elapsed time,
  with headroom (say 1.5-2x).

Rough throughput: total tasks = ceil(N_genomes / BATCH_SIZE), and wall-clock ≈
(total_tasks / SGE_TC) × per_task_time. For 340,000 genomes at BATCH_SIZE=1000
that's 340 tasks; at SGE_TC=100 and ~1 h/task, under ~4 hours of wall time. Raise
`SGE_TC` to go faster.

## Visualizing the results

The merge step also renders figures (or run them separately with `make plots`).
They're written to `batch/plots/` as **both PNG and SVG**, with the SVG text left
editable (`svg.fonttype=none`) so labels open as real text in Illustrator or
Inkscape. `scripts/make_plots.py` streams the big tables and holds only per-genome
scalars and a genome-by-pathway (never genome-by-step) matrix — about 20 MB at
340,000 proteomes — so the plots stay computable and meaningful at full scale.

Per pathway set (`aa`, `carbon`):

- `*_pathway_prevalence` — sorted stacked bars: fraction of proteomes where each
  pathway is present (all steps high-confidence), probably present (has a medium
  step), or a gap (has a low-confidence step). The biological headline.
- `*_pathways_per_proteome` — distribution of how many pathways are fully present
  per proteome.
- `*_qc_completeness` — per-proteome fraction of high-confidence steps (a quality
  signal; a low mode flags poor or unusual proteomes).
- `*_qc_size_vs_complete` — hexbin density of proteome size vs complete pathways
  (small/degraded proteomes and outliers stand out).
- `*_pathway_cooccurrence` — clustered pathway-by-pathway Jaccard heatmap
  (metabolic modules that tend to appear together).
- `*_top_gap_steps` — the steps most often missing across the dataset (which
  reactions are the common gaps).

Plus dataset-wide `qc_proteome_size` (proteins-per-proteome histogram),
`confidence_composition` (high/medium/low step fractions per set), and
`summary_stats.tsv` (N, median proteome size, median complete pathways, most/least
prevalent pathway per set).

The same script works on a single `make analyze` result:
`scripts/plots.sh PaperBLAST/tmp/<name> PaperBLAST/tmp/<name>/orgs.org PaperBLAST/tmp/<name>/plots`.

## Sources

- GapMind SETUP: https://github.com/morgannprice/PaperBLAST/blob/master/gaps/SETUP
- PaperBLAST README and data downloads: https://github.com/morgannprice/PaperBLAST
- "Conda installation" discussion: https://github.com/morgannprice/PaperBLAST/issues/16
- usearch 11 binaries (public domain): https://github.com/rcedgar/usearch_old_binaries
- TIGRFAMs 15.0: https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/
- Swissknife: https://swissknife.sourceforge.net/
