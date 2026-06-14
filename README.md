# GapMind command-line pipeline

A conda-based, scripted version of the GapMind install from
[`gaps/SETUP`](https://github.com/morgannprice/PaperBLAST/blob/master/gaps/SETUP).
It takes you from nothing to "analyze a proteome for amino-acid biosynthesis and
carbon catabolism pathways" with one environment file and four scripts. The full
from-scratch build (all three phases) is wired end to end.

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

Three dependencies can't come from conda, and the scripts force all three, so the
pipeline fetches them automatically:

- the prebuilt **PaperBLAST data** (`litsearch.db`, `uniq.faa`); Phase 1 reads it
  to assemble the curated proteins, and GapMind never rebuilds it itself;
- the **usearch 11** binary; Phase 1 clustering (`-cluster_fast`) and `.udb`
  creation in `setupGaps.pl` are hardcoded to usearch, so a from-scratch build
  needs it even though Phase 3 defaults to diamond;
- **Swissknife** (`SWISS::Entry`), the Perl parser the Swiss-Prot steps use, which
  the code expects unpacked at `SWISS/lib/`.

Phase 3 (the actual genome analysis) defaults to **diamond**, which is free,
installs from conda, and the GapMind author reports gives very similar results to
usearch.

## Requirements

A working `conda` (or `mamba`) and an internet connection. Designed for Linux
x86-64, the usual environment for this kind of work. Plan for roughly 30-50 GB of
free disk and several hours for the from-scratch build, most of it in Phase 1's
hmmsearch over every Pfam model and in the large downloads. See
[Caveats](#caveats) for Apple Silicon.

## Quickstart

```bash
# 1. create and activate the environment
make env
conda activate gapmind

# 2. full from-scratch build: clone, fetch tools + HMMs, build all databases
make all

# 3. analyze your proteome (a protein FASTA) for both pathway sets
make analyze FAA=/path/to/my_proteins.faa NAME=myorg
```

Outputs land in `PaperBLAST/tmp/myorg/`. Start with `aa.sum.rules` and
`carbon.sum.rules` (one row per pathway; `nHi`/`nMed`/`nLo` count high-, medium-,
and low-confidence steps), then `*.sum.steps` for the best candidate per step, and
`*.sum.db` for the combined sqlite3 database.

## The three phases

`make all` runs setup and then the build. Mapped to the SETUP document:

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

## Sources

- GapMind SETUP: https://github.com/morgannprice/PaperBLAST/blob/master/gaps/SETUP
- PaperBLAST README and data downloads: https://github.com/morgannprice/PaperBLAST
- "Conda installation" discussion: https://github.com/morgannprice/PaperBLAST/issues/16
- usearch 11 binaries (public domain): https://github.com/rcedgar/usearch_old_binaries
- TIGRFAMs 15.0: https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/
- Swissknife: https://swissknife.sourceforge.net/
