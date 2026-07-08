# config.sh — single source of configuration for the GapMind pipeline.
# Sourced by the Makefile and every script. Edit values here, then run `make all`.
# Every value can also be overridden from the environment, e.g. SETS="aa" make all

# Resolve this file's own directory so paths work regardless of the current dir.
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Where the PaperBLAST/GapMind code base is cloned and built.
# All databases are built in-tree under $CODE_DIR/tmp.
CODE_DIR="${CODE_DIR:-$PIPELINE_DIR/PaperBLAST}"

# Swissknife (SWISS::Entry etc.) is installed under $CODE_DIR/SWISS/lib, which is
# the directory that contains the SWISS/ package. Putting it on PERL5LIB makes it
# importable by every perl process here -- including the helper scripts that
# setupGaps.pl spawns internally -- no matter what each script's own `use lib`
# line says. Swissknife is pure perl, so this is safe across perl versions.
export PERL5LIB="$CODE_DIR/SWISS/lib${PERL5LIB:+:$PERL5LIB}"

# Git source for the code base.
PAPERBLAST_REPO="${PAPERBLAST_REPO:-https://github.com/morgannprice/PaperBLAST.git}"
PAPERBLAST_REF="${PAPERBLAST_REF:-master}"

# Which GapMind sets (pathway groups) to build / analyze. Space-separated.
#   aa     = amino acid biosynthesis
#   carbon = carbon catabolism
SETS="${SETS:-aa carbon}"

# Search tool for the ANALYSIS phase (Phase 3): diamond (free, recommended) or usearch.
# usearch is ALWAYS downloaded because the from-scratch build (Phase 1 clustering)
# hardcodes it; this setting only controls which tool Phase 3 uses.
SEARCH_TOOL="${SEARCH_TOOL:-diamond}"

# How to obtain the per-set curated + steps databases:
#   prebuilt -- download the maintainer's consistent curated.faa/curated.db/steps.db
#               (reliable, matches the official/validated GapMind, fast). Default.
#   build    -- rebuild them from scratch (Phases 1-2). Faithful to the sources,
#               but the pathway step files reference curated IDs from an older data
#               snapshot, so a current-data rebuild can fail on drifted UniProt/
#               BRENDA identifiers. Only the curated DB is rebuilt; PaperBLAST
#               itself is never rebuilt.
DB_SOURCE="${DB_SOURCE:-prebuilt}"

# CPU threads for searches and HMM steps.
THREADS="${THREADS:-4}"

# Conda environment name (kept in sync with the Makefile ENV variable).
CONDA_ENV="${CONDA_ENV:-gapmind}"

# ---------------------------------------------------------------------------
# External data sources. Override these if a mirror or release changes.
# ---------------------------------------------------------------------------

# Prebuilt per-set databases (DB_SOURCE=prebuilt). Files fetched per set:
#   $PREBUILT_BASE/path.<set>/{curated.faa,curated.db,steps.db}
PREBUILT_BASE="${PREBUILT_BASE:-https://papers.genomics.lbl.gov/tmp}"

# PaperBLAST curated-literature database (DB_SOURCE=build, Phase 1 reads
# litsearch.db + uniq.faa). ~1.5 GB total. GapMind never rebuilds PaperBLAST.
PAPERBLAST_DATA_BASE="${PAPERBLAST_DATA_BASE:-http://papers.genomics.lbl.gov/data}"

# Swiss-Prot flat file (Phase 1: curated2 sequences + heteromer detection).
SPROT_URL="${SPROT_URL:-https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz}"

# Pfam-A HMM library (Phase 1 pfam hits + Phase 2 step queries).
PFAM_URL="${PFAM_URL:-https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz}"

# TIGRFAMs HMM library, release 15.0 — the last full release (Phase 2 queries).
TIGRFAM_URL="${TIGRFAM_URL:-https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/TIGRFAMs_15.0_HMM.LIB.gz}"

# Swissknife (SWISS::Entry) Perl library — required by the Swiss-Prot parsers.
SWISSKNIFE_URL="${SWISSKNIFE_URL:-https://sourceforge.net/projects/swissknife/files/latest/download}"

# usearch 11 public-domain mirror (the OS-specific file is chosen in setup_code.sh).
USEARCH_BASE="${USEARCH_BASE:-https://raw.githubusercontent.com/rcedgar/usearch_old_binaries/main/bin}"

# ---------------------------------------------------------------------------
# Large-scale batch analysis on an SGE cluster (scripts/batch_*.sh)
# ---------------------------------------------------------------------------

# Tab-separated manifest of genomes to analyze. Header line required; must have a
# column named "faa" (path to a protein FASTA). Optional column "name" (a unique
# organism id); if absent, a sanitized name is derived from the file's basename.
MANIFEST="${MANIFEST:-}"

# Working + output tree for batch runs (put this on scratch with plenty of space).
BATCH_DIR="${BATCH_DIR:-$PIPELINE_DIR/batch}"

# Genomes per batch = per SGE array task. Larger = fewer tasks and better amortized
# startup; smaller = shorter tasks, finer parallelism, cheaper retries. 500 suits
# this cluster's 24h walltime and 8-core tasks; confirm with `make batch-calibrate`
# (340,000 genomes / 500 = 680 tasks).
BATCH_SIZE="${BATCH_SIZE:-500}"

# Keep the per-candidate table (*.sum.cand) in the merged output. It is by far the
# largest output at scale; set to 0 to keep only rules + steps.
KEEP_CAND="${KEEP_CAND:-1}"

# Presence/absence thresholding for the genome x pathway 0/1 matrices:
#   probably (default) -- present iff nLo==0 (every step at least medium conf.);
#                         this is GapMind's own "probably present" call.
#   strict             -- present iff nMed==0 and nLo==0 (every step high conf.).
# A companion 2/1/0 confidence matrix is always written so you can re-threshold.
PA_MODE="${PA_MODE:-probably}"

# ---- SGE submission (tuned for this cluster from qhost / qstat -g c / qconf) ----
# Observed: standard.q has 865 schedulable slots, idle, 24h walltime cap (s_rt
# 23:55:00, h_rt 24:00:00); PE "parallel" has allocation_rule $pe_slots, so a
# task's slots all land on ONE node (what GapMind's threading needs); no per-user
# slot quota. Everything a task reads/writes (CODE_DIR, BATCH_DIR, the conda env,
# the genome FASTAs) must be on shared storage visible from every compute node.

# CPU cores per task, all on one node (passed to gapsearch/gaprevsearch -nCPU).
SGE_SLOTS="${SGE_SLOTS:-8}"
# Parallel environment (used because SGE_SLOTS > 1). "parallel" is $pe_slots here.
SGE_PE="${SGE_PE:-parallel}"
# Concurrent array tasks: 100 x 8 = 800 of standard.q's 865 slots. Raise toward 108
# to use all of standard.q; or widen SGE_QUEUE (below) to also use long.q.
SGE_TC="${SGE_TC:-100}"
# Queue(s). standard.q is idle and allows 24h. Use "standard.q,long.q" to also pull
# long.q's free slots (also 24h); leave empty to let SGE pick any queue you can use.
SGE_QUEUE="${SGE_QUEUE:-standard.q}"
# Per-task limits, set from a calibration run: a 500-genome/8-core batch peaked at
# ~15 GB and ~44 min. 24 GB gives headroom for heavier batches (RAM is abundant,
# 500 GB - 2 TB per node, and cores bind before memory); 12h sits far above the
# ~44-min actual and under the queue's 24h cap.
SGE_H_VMEM="${SGE_H_VMEM:-24G}"
SGE_H_RT="${SGE_H_RT:-12:00:00}"
