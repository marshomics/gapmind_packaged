#!/usr/bin/env bash
# Install the PaperBLAST/GapMind code base into the image and make it run under
# the container's perl + tools. Scoped to the PREBUILT analysis path, so this is
# far smaller than the from-scratch setup: no Swissknife, Pfam, legacy BLAST, or
# runPfamHits/UniProt patching (none of those run when the databases are
# downloaded ready-made). Assumes hmmsearch/hmmfetch/diamond are already on PATH
# (from the conda env).
#   setup_paperblast.sh <dest> <git-ref>
set -euo pipefail

dest="${1:?usage: setup_paperblast.sh <dest> <git-ref>}"
ref="${2:?pin a commit SHA (reproducible) or a branch}"
repo="${PAPERBLAST_REPO:-https://github.com/morgannprice/PaperBLAST.git}"

git clone "$repo" "$dest"
git -C "$dest" checkout "$ref"
cd "$dest"

# Directories buildorgs.pl expects to exist (it dies otherwise).
mkdir -p fbrowse_data private tmp/downloaded

# The scripts ship with "#!/usr/bin/perl", which would use a system perl that
# lacks DBI etc. Point every perl script at the conda perl via PATH. (Analysis
# scripts call each other by path, so their shebang decides the interpreter.)
for f in bin/*.pl; do
  [ -f "$f" ] || continue
  if head -1 "$f" | grep -q '^#!.*perl'; then
    sed -i '1 s|^#!.*perl.*$|#!/usr/bin/env perl|' "$f"
  fi
done

# The scripts call $Bin/<tool> and $Bin/usearch by hardcoded path -- link the
# conda tools into bin/ so those resolve.
for t in hmmsearch hmmfetch diamond; do
  ln -sf "$(command -v "$t")" "bin/$t"
done

# usearch 11 (public-domain mirror); only used by the amino-acid known-gap step
# (orgsVsMarkers). Drop this line to build an ARM-only, fully free image.
curl -fsSL "https://raw.githubusercontent.com/rcedgar/usearch_old_binaries/main/bin/usearch11.0.667_i86linux64" \
  -o bin/usearch && chmod +x bin/usearch
got="$(md5sum bin/usearch | awk '{print $1}')"
[ "$got" = "fa050a3029d33b7b25036a2cc8c6da97" ] || { echo "usearch checksum mismatch: $got" >&2; exit 1; }

echo "PaperBLAST ready at $dest ($(git -C "$dest" rev-parse --short HEAD))"
