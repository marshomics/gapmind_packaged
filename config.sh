# config.sh — single source of configuration for the GapMind pipeline.
# Sourced by the Makefile and every script. Edit values here, then run `make all`.
# Every value can also be overridden from the environment, e.g. SETS="aa" make all

# Resolve this file's own directory so paths work regardless of the current dir.
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Where the PaperBLAST/GapMind code base is cloned and built.
# All databases are built in-tree under $CODE_DIR/tmp.
CODE_DIR="${CODE_DIR:-$PIPELINE_DIR/PaperBLAST}"

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

# CPU threads for searches and HMM steps.
THREADS="${THREADS:-4}"

# Conda environment name (kept in sync with the Makefile ENV variable).
CONDA_ENV="${CONDA_ENV:-gapmind}"

# ---------------------------------------------------------------------------
# External data sources. Override these if a mirror or release changes.
# ---------------------------------------------------------------------------

# PaperBLAST curated-literature database (Phase 1 reads litsearch.db + uniq.faa).
# ~1.5 GB total. GapMind does not rebuild PaperBLAST itself; it consumes this.
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
