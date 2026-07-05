# GapMind command-line pipeline

A conda-based, scripted version of the GapMind install from
[`gaps/SETUP`](https://github.com/morgannprice/PaperBLAST/blob/master/gaps/SETUP).
It takes you from nothing to "analyze a proteome for amino-acid biosynthesis and
carbon catabolism pathways" with one environment file and four scripts. By default
it downloads the maintainer's prebuilt curated + steps databases (reliable, and
identical to the official GapMind); rebuilding those databases from scratch is an
opt-in (`DB_SOURCE=build`).

## What this replaces

The original SETUP asks you to compile HMMer from source, fetch the usearch
binary by hand, download legacy NCBI BLAST off an FTP server, install a pile of
Perl modules, and drop each executable into a specific `bin/` folder. Almost all
of that becomes one `environment.yml`.

The trick is that GapMind's Perl scripts call their tools by hardcoded path:
`bin/hmmfetch`, `bin/blast/formatdb`, `bin/usearch`, and so on. So after conda
installs the tools, `setup_code.sh` symlinks each one into the exact path the
scripts expect. That glue is what lets a conda env stand in for the manual
binary drops.

A few dependencies can't come from conda, so the pipeline fetches them: the
**usearch 11** binary (Phase 1 clustering and the `orgsVsMarkers` known-gap step
hardcode it), the prebuilt **PaperBLAST data** (`litsearch.db`, `uniq.faa` — only
for a from-scratch rebuild), and **Swissknife** (`SWISS::Entry`, the Swiss-Prot
parser — also only for a rebuild). In the default prebuilt mode, only usearch is
needed.

Phase 3 (the actual genome analysis) defaults to **diamond**, which is free,
installs from conda, and the GapMind author reports gives very similar results to
usearch.

## Requirements

A working `conda` (or `mamba`) and an internet connection. Designed for Linux
x86-64, the usual environment for this kind of work. The default prebuilt mode
needs a few GB of disk (the per-set databases are hundreds of MB each). The opt-in
from-scratch build (`DB_SOURCE=build`) needs ~30-50 GB and several hours, most of
it in Phase 1's hmmsearch over every Pfam model. See [Caveats](#caveats) for Apple
Silicon.

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

## Caveats

The build is large and network-bound. Pfam-A is several GB once unzipped and
indexed; the PaperBLAST data is about 1.5 GB; Swiss-Prot adds another ~0.6 GB
compressed. Phase 2 fetches some sequences from UniProt at run time, so it needs
network access, not just the local files.

The usearch binary is x86-64 only. On Apple Silicon it runs through Rosetta 2
(the pipeline picks the `osx-64` build); diamond, hmmer, and blast-legacy all have
native arm64 conda builds. On ARM Linux there is no usearch build at all, so the
from-scratch build (which needs usearch for clustering) won't run there. A Linux
x86-64 host avoids every one of these edge cases.

The conda environment covers the command-line build and analysis. The optional
`gapView.cgi` web viewer needs a few more Perl modules; they're listed, commented
out, at the bottom of `environment.yml`.

Downloads track "current" releases (Pfam, Swiss-Prot) and the live PaperBLAST
data, so two builds months apart can differ. Pin the `*_URL` values in `config.sh`
to a specific release if you need reproducibility. The pipeline verifies the
usearch checksum; it does not pin checksums for the large reference downloads.

This covers the GapMind setup only. It consumes the prebuilt PaperBLAST database
rather than rebuilding it; that separate build is a much larger, multi-day job
(see the PaperBLAST README).

## Large-scale analysis on an SGE cluster

For many genomes (thousands to hundreds of thousands), `make analyze` (one genome
at a time) is the wrong tool. The batch mode takes a tab-separated manifest, splits
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

### Finding your cluster's capacity (SGE)

You need three numbers: how many slots you can use at once, a parallel environment
that packs slots onto one node, and a per-task memory/time footprint.

```bash
qhost                       # every exec host: NCPU (cores) and MEMTOT (RAM)
qstat -g c                  # per queue: USED / AVAIL / TOTAL schedulable slots
qconf -spl                  # parallel environments; pick an SMP/threaded one
qquota -u "$USER"           # your resource-quota usage and limits (if any)
qconf -srqs                 # cluster-wide resource quotas (per-user slot caps)
```

`qhost` inventories the hardware (sum the NCPU column for total cores). `qstat -g c`
gives the schedulable slots — `TOTAL` is the ceiling, `AVAIL` is free right now.
GapMind threads within a single node (it is not MPI), so pick an SMP-style parallel
environment from `qconf -spl` (commonly `smp`, `thread`, `threaded`, or `shm`) and
put it in `SGE_PE`. If there's a per-user slot cap (from `qconf -srqs` / `qquota`,
or your admin), that cap — not the cluster total — is your real budget.

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
`SGE_TC` (within your slot budget) to go faster.

### Accuracy, robustness, and disk

Batch size does not change the confidence calls. The only batch-dependent knob is the
initial candidate E-value in `gapsearch.pl`, which the code deliberately scales by the
number of organisms in the batch to keep per-genome detection stable as the search
database grows; the Hi/Med/Lo scoring in `gapsummary.pl` uses absolute bit-score,
identity, and coverage thresholds and is unaffected. Use one consistent `BATCH_SIZE`
for a fully uniform run.

Every task is idempotent: a finished batch writes a `.done` marker and is skipped on
resubmit, so after failures you just fix the cause and run `make batch-submit` again —
only the incomplete batches rerun. Each task works in its own directory and its
search temp files are PID-scoped, so tasks sharing a node don't collide. `batch-merge`
reports how many batches are complete and lists any stragglers in
`batch/incomplete.txt`.

Two things to watch. `<set>.sum.cand` (one row per candidate per step per genome) is
by far the largest output — tens of GB at 340k genomes — so set `KEEP_CAND=0` if you
only need pathway and step calls. And GapMind keys file-based genomes by the MD5 of
their sequences, so two genomes with byte-identical proteomes collapse to one result
(the same call for both); dedupe upstream if you need a row per file regardless.

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

## Troubleshooting

`Can't locate DBI.pm in @INC` (or another missing module), with an `@INC` that
lists only system paths like `/usr/lib/x86_64-linux-gnu/perl/5.34`: the script
ran under the system perl instead of the conda env's perl. The repo ships
`#!/usr/bin/perl` shebangs, and `setupGaps.pl` / `buildStepsDb.pl` call the other
scripts by path, so the shebang chooses the interpreter. `setup_code.sh` rewrites
those shebangs to `#!/usr/bin/env perl`; re-run `make setup-code` to apply it to
an existing clone, then `make databases` again. The downloads and `sprot.curated_parsed`
are cached, so it resumes quickly. Confirm the env perl can load DBI with
`conda run -n gapmind perl -MDBI -MDBD::SQLite -e 'print "ok\n"'`. (The rewrite
covers `bin/*.pl`; the optional web viewer's `cgi/*.cgi` would need the same.)

`Can't locate SWISS/Entry.pm in @INC` during Phase 1: a Swiss-Prot helper that
`setupGaps.pl` runs couldn't find Swissknife. `config.sh` exports
`PERL5LIB=$CODE_DIR/SWISS/lib` so every perl process can import it; make sure
you're invoking the pipeline through the Makefile (or have `config.sh` sourced)
so that variable is set, then re-run `make databases`. Check the library is in
place with `ls "$CODE_DIR"/SWISS/lib/SWISS/Entry.pm`.

`Failed to fetch https://...uniprot.../<id>.txt` / `gapquery.pl failed` during
Phase 2: a query sequence couldn't be fetched. Two causes, both handled by the
`setup_code.sh` patch to `lib/Steps.pm`: the old `www.uniprot.org/uniprot/`
endpoint is retired (it now fetches via curl from `rest.uniprot.org`), and some
accessions in the step files have been DELETED from UniProtKB (redundant-proteome
cleanup), for which it recovers the archived sequence from UniParc via the
inactive-entry JSON. Re-run `make setup-code`, then `make databases`. Fetched
sequences are cached in `static/uniprotCache.tsv`. If it still fails, either the
node can't reach `rest.uniprot.org` (build where outbound HTTPS works) or a
deleted entry has no UniParc record; in that case download the prebuilt
`steps.db` from `papers.genomics.lbl.gov/tmp/path.<set>/` instead of rebuilding
Phase 2.

The Phase 1/2 build is resumable: `curated.db` marks Phase 1 done for a set and
`steps.db` marks Phase 2 done, so re-running `make databases` skips finished
sets and continues where it stopped.

## Sources

- GapMind SETUP: https://github.com/morgannprice/PaperBLAST/blob/master/gaps/SETUP
- PaperBLAST README and data downloads: https://github.com/morgannprice/PaperBLAST
- "Conda installation" discussion: https://github.com/morgannprice/PaperBLAST/issues/16
- usearch 11 binaries (public domain): https://github.com/rcedgar/usearch_old_binaries
- TIGRFAMs 15.0: https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/
- Swissknife: https://swissknife.sourceforge.net/
