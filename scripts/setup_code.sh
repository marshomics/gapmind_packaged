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
