#!/usr/bin/env bash
# Prepare an SGE batch run from a TSV manifest of protein FASTAs.
#
#   scripts/batch_prepare.sh [manifest.tsv]
#
# The manifest needs a header line with a column "faa" (path to a .faa) and an
# optional column "name" (a unique organism id; derived from the filename if
# absent). One pass: validates each file, sanitizes/dedupes names, stages any
# path containing whitespace or ":" via a symlink, and splits the valid genomes
# into $BATCH_DIR/batches/NNNNN.orgfile of BATCH_SIZE each.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$here/config.sh"

manifest="${1:-$MANIFEST}"
[ -n "$manifest" ] || { echo "ERROR: no manifest. Pass one, or set MANIFEST in config.sh" >&2; exit 1; }
[ -s "$manifest" ] || { echo "ERROR: manifest not found or empty: $manifest" >&2; exit 1; }
[ -d "$CODE_DIR/tmp" ] || { echo "ERROR: databases not built yet. Run 'make databases' first." >&2; exit 1; }

mkdir -p "$BATCH_DIR/batches" "$BATCH_DIR/results" "$BATCH_DIR/logs" "$BATCH_DIR/staging"

# Symlinks so checkGapRequirements.pl (which reads <results>/path.<set>/steps.db)
# finds the shared, already-built databases without copying them.
for set in $SETS; do
  [ -d "$CODE_DIR/tmp/path.$set" ] || { echo "ERROR: missing $CODE_DIR/tmp/path.$set (build databases first)" >&2; exit 1; }
  ln -sfn "$CODE_DIR/tmp/path.$set" "$BATCH_DIR/results/path.$set"
done

echo ">> Preparing batches from $manifest (BATCH_SIZE=$BATCH_SIZE)"

BATCH_DIR="$BATCH_DIR" BATCH_SIZE="$BATCH_SIZE" perl -MCwd=abs_path -MFile::Basename -e '
  my $batchDir = $ENV{BATCH_DIR};
  my $batchSize = $ENV{BATCH_SIZE} + 0; $batchSize = 300 if $batchSize < 1;
  open(my $man, "<", $ARGV[0]) or die "Cannot read $ARGV[0]\n";
  my $header = <$man>; $header =~ s/[\r\n]+$//;
  my @cols = split /\t/, $header;
  my %idx; $idx{$cols[$_]} = $_ for 0..$#cols;
  die "Manifest has no \"faa\" column (header: $header)\n" unless exists $idx{faa};
  my $faaI = $idx{faa};
  my $nameI = exists $idx{name} ? $idx{name} : -1;

  open(my $skip, ">", "$batchDir/skipped.tsv") or die;
  print $skip "faa\treason\n";

  my %seen; my ($nValid,$nSkip,$b,$inBatch) = (0,0,0,0);
  my $fh;
  while (my $line = <$man>) {
    $line =~ s/[\r\n]+$//;
    next if $line eq "";
    my @F = split /\t/, $line, -1;
    my $path = defined $F[$faaI] ? $F[$faaI] : "";
    next if $path eq "";
    my $name = ($nameI >= 0 && defined $F[$nameI]) ? $F[$nameI] : "";

    if ($name eq "") { $name = basename($path); $name =~ s/\.(faa|fasta|fa|fna|pep)$//i; }
    $name =~ s/[^A-Za-z0-9_.-]/_/g;
    $name = "g" . $name if $name !~ /^[A-Za-z]/;   # start with a letter
    $name = "genome" if $name eq "";
    if (exists $seen{$name}) { $seen{$name}++; $name .= "__" . $seen{$name}; $seen{$name} = 1; }
    else { $seen{$name} = 1; }

    if ($path =~ /\.gz$/) { print $skip "$path\tgzipped (decompress first)\n"; $nSkip++; next; }
    unless (-s $path) { print $skip "$path\tmissing_or_empty\n"; $nSkip++; next; }
    my $first = "";
    if (open(my $ff, "<", $path)) { $first = <$ff>; close $ff; }
    unless (defined $first && $first =~ /^>/) { print $skip "$path\tnot_fasta\n"; $nSkip++; next; }

    my $abs = abs_path($path);
    my $use = $abs;
    if ($abs =~ /[\s:]/) {                          # buildorgs cannot take these
      my $link = "$batchDir/staging/$name.faa";
      unlink $link;
      symlink($abs, $link) or die "Cannot symlink $abs -> $link\n";
      $use = $link;
    }

    if ($inBatch == 0) { $b++; open($fh, ">", sprintf("%s/batches/%05d.orgfile", $batchDir, $b)) or die; }
    print $fh "file:$use:$name\n";
    $inBatch++; $nValid++;
    if ($inBatch >= $batchSize) { close $fh; $inBatch = 0; }
  }
  close $fh if $inBatch > 0;
  close $skip;
  open(my $nb, ">", "$batchDir/nbatches.txt") or die; print $nb "$b\n"; close $nb;
  print STDERR sprintf("   %d genomes -> %d batches of up to %d (%d skipped; see %s/skipped.tsv)\n",
                       $nValid, $b, $batchSize, $nSkip, $batchDir);
' "$manifest"

nb=$(cat "$BATCH_DIR/nbatches.txt")
echo ">> Prepared $nb batches in $BATCH_DIR/batches/"
echo ">> Next: calibrate one batch, then submit. See 'make batch-calibrate' and 'make batch-submit'."
