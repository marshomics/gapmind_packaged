#!/usr/bin/env bash
# Produce the CANONICAL reference: run the original PaperBLAST/GapMind scripts one
# genome at a time (nOrgs=1, as the GapMind web server does) against a given
# database, and concatenate into merged reference tables to compare with the
# pipeline output. Point --code-dir at a PaperBLAST clone (the original code) and
# --db at the SAME databases the pipeline used. Needs perl/diamond/hmmer on PATH
# (run inside the conda env or the container).
#
#   run_reference.sh --manifest genomes.tsv --code-dir PaperBLAST --db PaperBLAST/tmp \
#                    --sets "aa carbon" --tool diamond --out reference [--ncpu 4]
set -euo pipefail

manifest=""; code="${GAPMIND_DIR:-}"; db=""; sets="aa carbon"; tool="diamond"; out="reference"; ncpu=4
while [ $# -gt 0 ]; do
  case "$1" in
    --manifest) manifest="$2"; shift 2;;
    --code-dir) code="$2";     shift 2;;
    --db)       db="$2";       shift 2;;
    --sets)     sets="$2";     shift 2;;
    --tool)     tool="$2";     shift 2;;
    --out)      out="$2";      shift 2;;
    --ncpu)     ncpu="$2";     shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[ -n "$manifest" ] && [ -s "$manifest" ] || { echo "ERROR: --manifest genomes.tsv required" >&2; exit 1; }
[ -n "$code" ] && [ -d "$code/bin" ]     || { echo "ERROR: --code-dir must be a PaperBLAST clone" >&2; exit 1; }
[ -n "$db" ] && [ -d "$db" ]             || { echo "ERROR: --db must be a dir with path.<set>/" >&2; exit 1; }

b="$code/bin"
export MC_CORES="$ncpu"
mkdir -p "$out" "$code/fbrowse_data" "$code/private" "$code/tmp/downloaded" 2>/dev/null || true
for s in $sets; do : > "$out/$s.sum.rules"; : > "$out/$s.sum.steps"; done
: > "$out/orgs.tsv"; owrote=0
work="$(mktemp -d)"

# manifest -> "name<TAB>faa" (derive+sanitize name from filename if absent)
awk -F'\t' '
  NR==1 { for (i=1;i<=NF;i++){ if($i=="faa")fa=i; if($i=="name")nm=i }; next }
  { p=$fa; n=(nm?$nm:"");
    if (n=="") { k=p; sub(/.*\//,"",k); sub(/\.(faa|fasta|fa|fna|pep)$/,"",k); n=k }
    gsub(/[^A-Za-z0-9_.-]/,"_",n); if (p!="") print n "\t" p }' "$manifest" > "$work/list"

while IFS=$'\t' read -r name faa; do
  [ -s "$faa" ] || { echo "  skip (missing): $faa" >&2; continue; }
  g="$work/$name"; mkdir -p "$g"
  printf 'file:%s:%s\n' "$faa" "$name" > "$g/orgfile"
  (
    cd "$g"
    "$b/buildorgs.pl" -out orgs -orgfile orgfile -cache cache >/dev/null
    [ "$tool" = usearch ] || diamond makedb --quiet --in orgs.faa -d orgs.faa.dmnd >/dev/null 2>&1
    for s in $sets; do ln -sfn "$db/path.$s" "path.$s"; done
    for s in $sets; do
      if [ "$tool" = usearch ]; then
        "$b/gapsearch.pl"    -orgs orgs -set "$s" -dir "$db/path.$s" -out "$s.hits" -nCPU "$ncpu" >/dev/null 2>&1
        "$b/gaprevsearch.pl" -orgs orgs -hits "$s.hits" -curated "$db/path.$s/curated.faa.udb" -out "$s.revhits" -nCPU "$ncpu" >/dev/null 2>&1
      else
        "$b/gapsearch.pl"    -diamond -orgs orgs -set "$s" -dir "$db/path.$s" -out "$s.hits" -nCPU "$ncpu" >/dev/null 2>&1
        "$b/gaprevsearch.pl" -diamond -orgs orgs -hits "$s.hits" -curated "$db/path.$s/curated.faa.dmnd" -out "$s.revhits" -nCPU "$ncpu" >/dev/null 2>&1
      fi
      "$b/gapsummary.pl" -orgs orgs -set "$s" -dbDir "$db/path.$s" -hits "$s.hits" -revhits "$s.revhits" -out "$s.sum" >/dev/null 2>&1
    done
  )
  for s in $sets; do
    for suf in rules steps; do
      f="$g/$s.sum.$suf"; [ -s "$f" ] || continue
      if [ ! -s "$out/$s.sum.$suf" ]; then cat "$f" >> "$out/$s.sum.$suf"
      else tail -n +2 "$f" >> "$out/$s.sum.$suf"; fi
    done
  done
  if [ -s "$g/orgs.org" ]; then
    if [ "$owrote" = 0 ]; then cat "$g/orgs.org" >> "$out/orgs.tsv"; owrote=1
    else tail -n +2 "$g/orgs.org" >> "$out/orgs.tsv"; fi
  fi
  rm -rf "$g"
  echo "  reference: $name" >&2
done < "$work/list"
rm -rf "$work"
echo ">> reference tables (nOrgs=1, tool=$tool) in $out"
