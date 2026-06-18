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
