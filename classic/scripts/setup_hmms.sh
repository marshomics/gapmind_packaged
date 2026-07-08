#!/usr/bin/env bash
# Download and index the Pfam-A and TIGRFAMs HMM libraries, and build the
# metadata tables (pfam.tab, tigrinfo) that gapquery.pl reads. Needed for
# both Phase 1 (pfam hits) and Phase 2 (step queries).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$here/config.sh"

# In prebuilt mode the HMM databases aren't needed: the downloaded steps.db
# already carries its HMM models (extractHmms.pl pulls them out), and the
# Phase 1/2 steps that read hmm/Pfam-A.hmm and TIGRFAMs.hmm never run.
if [ "$DB_SOURCE" = "prebuilt" ]; then
  echo ">> DB_SOURCE=prebuilt: HMM databases not needed; skipping setup_hmms."
  exit 0
fi

cd "$CODE_DIR/hmm"

# --- Pfam-A ----------------------------------------------------------------
if [ ! -s Pfam-A.hmm ]; then
  echo ">> Downloading Pfam-A"
  curl -fSL "$PFAM_URL" -o Pfam-A.hmm.gz
  gunzip -f Pfam-A.hmm.gz
fi

# --- TIGRFAMs (canonical name expected by the code is TIGRFAMs.hmm) --------
if [ ! -s TIGRFAMs_15.0_HMM.LIB ]; then
  echo ">> Downloading TIGRFAMs 15.0"
  curl -fSL "$TIGRFAM_URL" -o TIGRFAMs_15.0_HMM.LIB.gz
  gunzip -f TIGRFAMs_15.0_HMM.LIB.gz
fi
ln -sf TIGRFAMs_15.0_HMM.LIB TIGRFAMs.hmm

# --- index both libraries --------------------------------------------------
# hmmpress builds the binary db for hmmscan; hmmfetch --index builds the SSI
# index so hmmfetch can pull a model by name/accession quickly (gapquery.pl
# and runPfamHits.pl fetch thousands of models, so the index is essential).
for f in Pfam-A.hmm TIGRFAMs.hmm; do
  [ -e "$f.h3i" ] || hmmpress -f "$f"
  [ -e "$f.ssi" ] || hmmfetch --index "$f"
done

# --- metadata tables gapquery.pl expects in hmm/ ---------------------------
cp -f ../gaps/tigrinfo ./tigrinfo
perl -e 'print "acc\tname\n";
  while(<STDIN>){ chomp;
    $name=$1 if m/^NAME\s+(.*)$/;
    $acc =$1 if m/^ACC\s+(.*)$/;
    next unless $_ eq "//";
    print "$acc\t$name\n" if $name && $acc;
    $name=undef; $acc=undef; }' < Pfam-A.hmm > pfam.tab

echo ">> HMM databases ready in $CODE_DIR/hmm"
