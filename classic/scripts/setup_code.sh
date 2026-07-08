#!/usr/bin/env bash
# Clone the code base and wire conda-provided tools into the paths the Perl
# scripts expect (bin/<tool> and bin/blast/<tool>). Also installs the two
# dependencies that are not on conda: the usearch 11 binary and Swissknife.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$here/config.sh"

# --- 1. clone (or update) PaperBLAST/GapMind ------------------------------
if [ ! -d "$CODE_DIR/.git" ]; then
  echo ">> Cloning $PAPERBLAST_REPO ($PAPERBLAST_REF) into $CODE_DIR"
  git clone --depth 1 -b "$PAPERBLAST_REF" "$PAPERBLAST_REPO" "$CODE_DIR"
else
  echo ">> Updating existing clone in $CODE_DIR"
  git -C "$CODE_DIR" pull --ff-only || echo "   (could not fast-forward; keeping current checkout)"
fi

cd "$CODE_DIR"

# --- 2. working directories the scripts assume ----------------------------
mkdir -p bin/blast fbrowse_data private tmp/downloaded hmm ind data

# --- 2a. replace the cluster-oriented runPfamHits.pl ----------------------
# Upstream runPfamHits.pl fans ~20k hmmfetch+hmmsearch jobs out through
# submitter.pl, which is built for ssh/cluster dispatch and load-throttles on
# uptime; in this environment the per-model .domtbl files never get produced.
# Drop in an equivalent that runs the same search -- every Pfam-A model at its
# trusted cutoff against curated.faa -- in one hmmsearch, emitting the identical
# pfam.hits.tab columns. (hmmsearch reads the multi-model Pfam-A.hmm directly.)
cat > bin/runPfamHits.pl <<'RUNPFAM'
#!/usr/bin/env bash
# runPfamHits.pl -- pipeline-managed drop-in replacement (see setup_code.sh).
# Same output as upstream (tmp/path.<set>/pfam.hits.tab), via one hmmsearch
# instead of submitter.pl job dispatch.
set -euo pipefail
set=${1:?Usage: runPfamHits.pl <set>}
work="tmp/path.$set"
[ -d "$work" ] || { echo "Directory $work does not exist!" >&2; exit 1; }
for f in hmm/Pfam-A.hmm bin/hmmsearch "$work/curated.faa"; do
  [ -e "$f" ] || { echo "file $f does not exist!" >&2; exit 1; }
done
echo "Computing $work/pfam.hits.tab"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/pfam.hits.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
bin/hmmsearch --cut_tc --cpu "${MC_CORES:-4}" \
  --domtblout "$tmp/all.domtbl" -o /dev/null \
  hmm/Pfam-A.hmm "$work/curated.faa"
perl -ane 'next if m/^[#-]/;
  ($ids,$dash1,$qlen,$hmmName,$hmmAcc,$hmmLen,$seqEval,$seqBits,$seqBias,$domI,$domN,$evalC,$evalE,$bits,$bias,$hmmFrom,$hmmTo,$seqFrom,$seqTo)=split / +/;
  next unless defined $seqTo && $dash1 eq "-" && $ids =~ m/:/;
  print join("\t",$ids,$hmmName,$hmmAcc,$evalC,$bits,$seqFrom,$seqTo,$qlen,$hmmFrom,$hmmTo,$hmmLen)."\n";' \
  "$tmp/all.domtbl" | sort > "$work/pfam.hits.tab"
echo "Wrote $work/pfam.hits.tab"
RUNPFAM
chmod +x bin/runPfamHits.pl
echo "   installed robust bin/runPfamHits.pl"

# --- 2b. make the scripts use the conda perl, not /usr/bin/perl -----------
# The repo scripts start with "#!/usr/bin/perl", which ignores PATH and runs
# the system perl (which lacks DBI, DBD::SQLite, etc.). setupGaps.pl and
# buildStepsDb.pl invoke the other scripts by path, so their shebang decides the
# interpreter. Rewrite every perl shebang to "#!/usr/bin/env perl" so they pick
# up the conda env's perl (and its modules) from PATH. Idempotent; leaves shell
# scripts such as runPfamHits.pl untouched.
n_sb=0
for f in bin/*.pl; do
  [ -f "$f" ] || continue
  if head -1 "$f" | grep -q '^#!.*perl'; then
    perl -i -pe 's{^#!.*\bperl\b.*$}{#!/usr/bin/env perl} if $. == 1' "$f"
    n_sb=$((n_sb + 1))
  fi
done
echo "   normalized perl shebang in $n_sb scripts"

# Some bin/ scripts are bash but ship with no shebang (first line is a bare "#"),
# so the kernel's ENOEXEC fallback runs them with /bin/sh. Where /bin/sh is dash
# (Debian/Ubuntu), their bash syntax ("[[ ... ]]", ">& file") breaks -- e.g.
# runPfamHits.pl, which setupGaps.pl runs in Phase 1. Give any shebang-less text
# script that uses such syntax an explicit bash shebang. The guards skip the conda
# tool symlinks and the usearch binary so neither is ever modified.
n_bash=0
for f in bin/*; do
  [ -f "$f" ] || continue                              # skip directories
  [ -L "$f" ] && continue                              # skip the tool symlinks
  grep -Iq . "$f" || continue                          # skip binary files (usearch)
  case "$(head -1 "$f")" in '#!'*) continue ;; esac    # already has a shebang
  if grep -Eq '\[\[|>&' "$f"; then
    { echo '#!/usr/bin/env bash'; cat "$f"; } > "$f.shebang.tmp" \
      && mv "$f.shebang.tmp" "$f" && chmod +x "$f"
    n_bash=$((n_bash + 1))
  fi
done
echo "   added bash shebang to $n_bash shebang-less script(s)"

# --- 2c. make UniProt fetches work (modern endpoint + UniParc fallback) ----
# gapquery.pl (Phase 2) fetches query sequences via Steps.pm. Two problems:
#  (1) it used the retired https://www.uniprot.org/uniprot/<id>.txt endpoint; and
#  (2) some accessions in the step files have since been DELETED from UniProtKB
#      ("redundant proteome" cleanup), so a live fetch returns nothing and
#      gapquery.pl aborts (e.g. G8ALI9 in carbon's D-alanine pathway).
# Modernize the URL (used in the error message), route the fetch through curl
# (reliable for HTTPS even where perl's TLS stack is not), and -- when an entry
# is gone -- recover its sequence from the UniParc archive named in the entry's
# inactive JSON. UniParc preserves the identical sequence, so the rebuild stays
# faithful.
if grep -q 'www\.uniprot\.org/uniprot/' lib/Steps.pm 2>/dev/null; then
  perl -i -pe 's{https://www\.uniprot\.org/uniprot/}{https://rest.uniprot.org/uniprotkb/}g' lib/Steps.pm
fi
if ! grep -q 'sub FetchUniProtFlat' lib/Steps.pm 2>/dev/null; then
  perl -0777 -i -pe 's{\$content = get\(\$URL\);}{\$content = Steps::FetchUniProtFlat(\$id);}' lib/Steps.pm
  cat >> lib/Steps.pm <<'STEPSPATCH'

# --- pipeline-added: robust UniProt fetch (curl + UniParc fallback) ---
# Returns the UniProtKB flatfile for $id, or -- if the entry has been deleted
# from UniProtKB -- a synthesized flatfile carrying the archived UniParc sequence
# so the existing parser in FetchUniProtSequence still works.
sub FetchUniProtFlat {
  my ($id) = @_;
  my $txt = `curl -fsSL "https://rest.uniprot.org/uniprotkb/${id}.txt" 2>/dev/null`;
  return $txt if defined $txt && $txt =~ /\nSQ /;
  my $json = `curl -fsSL "https://rest.uniprot.org/uniprotkb/${id}.json" 2>/dev/null`;
  if (defined $json && $json =~ /"uniParcId"\s*:\s*"(UPI[0-9A-Fa-f]+)"/) {
    my $upi = $1;
    my $fa = `curl -fsSL "https://rest.uniprot.org/uniparc/${upi}.fasta" 2>/dev/null`;
    if (defined $fa && $fa =~ /\n/) {
      my @fl = split /\n/, $fa;
      shift @fl;                       # drop the ">" header line
      my $seq = join("", @fl); $seq =~ s/\s//g;
      return "DE   UniParc $upi (UniProtKB $id deleted)\n"
           . "SQ   SEQUENCE " . length($seq) . " AA;\n     $seq\n//\n"
        if $seq =~ /^[A-Z]+\z/ && length($seq) > 0;
    }
  }
  return undef;
}
1;
STEPSPATCH
  echo "   patched lib/Steps.pm: curl + UniParc fallback for UniProt fetches"
fi

# --- 3. symlink conda tools into bin/ and bin/blast/ ----------------------
# The perl scripts call "$Bin/hmmfetch", "$Bin/blast/formatdb", "$Bin/usearch"
# etc. by hardcoded relative path, so the conda executables must appear there.
link_tool () {  # link_tool <dest-path> <executable-name>
  local dest="$1" name="$2" src
  src="$(command -v "$name" || true)"
  if [ -z "$src" ]; then
    echo "ERROR: '$name' is not on PATH. Activate the conda env first:" >&2
    echo "       conda activate $CONDA_ENV" >&2
    exit 1
  fi
  ln -sf "$src" "$dest"
  echo "   linked bin: $dest -> $src"
}
for t in hmmsearch hmmfetch hmmscan diamond; do link_tool "bin/$t" "$t"; done
for t in formatdb blastall fastacmd;        do link_tool "bin/blast/$t" "$t"; done

# --- 4. usearch 11 binary (required by the from-scratch build) ------------
os="$(uname -s)"; arch="$(uname -m)"
case "$os/$arch" in
  Linux/x86_64)        ub=usearch11.0.667_i86linux64; md5=fa050a3029d33b7b25036a2cc8c6da97 ;;
  Darwin/x86_64|Darwin/arm64) ub=usearch11.0.667_i86osx64; md5=9d3d2cf73deb1880976643e5abfe371d ;;
  *) echo "ERROR: no known usearch 11 binary for $os/$arch. Put one at bin/usearch manually." >&2; exit 1 ;;
esac
if [ ! -x bin/usearch ]; then
  echo ">> Downloading usearch ($ub)"
  curl -fSL "$USEARCH_BASE/$ub" -o bin/usearch
  chmod +x bin/usearch
  got="$( { md5sum bin/usearch 2>/dev/null || md5 -q bin/usearch; } | awk '{print $1}')"
  if [ "$got" != "$md5" ]; then
    echo "ERROR: usearch checksum mismatch (got $got, expected $md5)" >&2
    rm -f bin/usearch; exit 1
  fi
  echo "   usearch checksum OK"
fi

# --- 5. Swissknife (SWISS::Entry) -> SWISS/ so SWISS/lib/SWISS/*.pm exist --
if [ ! -e SWISS/lib/SWISS/Entry.pm ]; then
  echo ">> Installing Swissknife into SWISS/"
  tmp="$(mktemp -d)"
  curl -fSL "$SWISSKNIFE_URL" -o "$tmp/swissknife.tar.gz"
  mkdir -p "$tmp/x"
  tar -xzf "$tmp/swissknife.tar.gz" -C "$tmp/x"
  d="$(find "$tmp/x" -mindepth 1 -maxdepth 1 -type d | head -1)"
  if [ -z "$d" ]; then echo "ERROR: unexpected Swissknife archive layout" >&2; exit 1; fi
  rm -rf SWISS && mv "$d" SWISS
  rm -rf "$tmp"
fi
if [ ! -e SWISS/lib/SWISS/Entry.pm ]; then
  echo "ERROR: Swissknife not laid out at SWISS/lib/SWISS/Entry.pm" >&2
  exit 1
fi
echo "   Swissknife OK"

echo ">> Code base ready: $CODE_DIR"
