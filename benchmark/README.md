# Benchmarking against the original GapMind

This validates that the pipeline (SGE or Nextflow) reproduces GapMind and
quantifies the two deliberate departures from canonical GapMind — diamond as the
default search tool, and batched multi-genome runs. The pipeline calls the exact
PaperBLAST scripts (`gapsearch.pl`, `gapsummary.pl`, ...), so under matched
database, tool, and batch size the results should be **identical**; anything less
is a wrapper bug to find, not a scientific difference. That makes "faithful
re-implementation" a falsifiable claim.

Two tools:

- `run_reference.sh` — runs the original scripts one genome at a time (`nOrgs=1`,
  as the GapMind web server does) against a chosen database, giving the canonical
  reference tables.
- `concordance.py` — joins two output dirs on `orgId` (GapMind derives it from the
  proteome MD5, so identical sequences match exactly) and reports agreement.

Run everything inside the `gapmind` conda env (or the container): the scripts need
perl, diamond, hmmer, numpy on `PATH`.

## A defensible test panel

You do not need 340k genomes to prove equivalence — a few hundred is stronger and
tractable. Suggested panel (put paths in a `bench.tsv` with a `faa` column):

- a taxonomically diverse set (e.g. a few hundred GTDB species representatives
  across bacterial + archaeal phyla), plus
- the experimental validation organisms from the GapMind papers — the 35 bacteria
  grown in defined media (amino-acid auxotrophy is ground truth) and the 29
  fitness-profiled bacteria (carbon), and
- edge cases: a partial/tiny proteome, a curated known-gap organism, and a
  duplicated proteome.

Freeze one database release and one set of tool versions and use them for **both**
sides. Record the PaperBLAST commit SHA, the prebuilt-DB date/checksum, the
diamond/hmmer versions, the container image digest, and the Nextflow version.

## Tier 1 — exact concordance vs the original (the core claim)

Reference from the original scripts, then the pipeline on the same genomes with
the same DB and tool, then score:

```bash
# reference: original code, one genome at a time
benchmark/run_reference.sh --manifest bench.tsv \
    --code-dir PaperBLAST --db PaperBLAST/tmp --sets "aa carbon" \
    --tool diamond --out bench/reference

# test: the pipeline (SGE shown; or `nextflow run ... --input bench.csv`)
BATCH_DIR=$PWD/bench/run BATCH_SIZE=1 MANIFEST=bench.tsv make batch-prepare
BATCH_DIR=$PWD/bench/run make batch-submit        # wait, then it merges to bench/run/merged

# score
python benchmark/concordance.py --ref bench/reference --test bench/run/merged \
    --sets "aa carbon" --out bench/concordance_b1
```

Expect ~1.000 exact step and pathway agreement, κ≈1, zero discordances. Any
discordant rows are written to `bench/concordance_b1/<set>.discordant_*.tsv` —
investigate each; they should reduce to nondeterminism or a fixable bug.

## Tier 2 — quantify the deliberate departures

Batch size (only `gapsearch`'s initial E-value prefilter scales with `nOrgs`, by
design; the confidence calls should be stable):

```bash
for B in 1 100 500 2000; do
  BATCH_DIR=$PWD/bench/b$B BATCH_SIZE=$B MANIFEST=bench.tsv make batch-prepare
  BATCH_DIR=$PWD/bench/b$B make batch-submit
done
python benchmark/concordance.py --ref bench/b1/merged --test bench/b500/merged --sets "aa carbon" --out bench/batch_1_vs_500
```

Search tool (justifies the diamond default; the GapMind author reports "very
similar" results — put a number on it):

```bash
benchmark/run_reference.sh --manifest bench.tsv --code-dir PaperBLAST --db PaperBLAST/tmp --tool usearch  --out bench/ref_usearch
benchmark/run_reference.sh --manifest bench.tsv --code-dir PaperBLAST --db PaperBLAST/tmp --tool diamond  --out bench/ref_diamond
python benchmark/concordance.py --ref bench/ref_usearch --test bench/ref_diamond --sets "aa carbon" --out bench/diamond_vs_usearch
```

## Tier 3 — determinism and portability

- Re-run identical input; `md5sum bench/run/merged/*.sum.*` should be unchanged.
- Nextflow: run the panel under `-profile local`, `sge`, and `slurm`; concordance
  across them should be exact (the executor is scientifically inert).
- Kill a run mid-flight and `-resume` (Nextflow) or re-`make batch-submit` (SGE);
  concordance vs a clean run should be exact.

## Tier 4 — external biological ground truth

Reproduce the papers' validation: amino-acid auxotrophy predictions vs the 35
defined-media organisms, and carbon catabolism vs the 29 fitness-profiled ones;
report sensitivity/specificity and compare to the published numbers. Confirm the
curated known-gap cases are flagged (`*.sum.knownsim`). Reproducing a table/figure
from the papers with the pipeline is a strong, citable check.

## Tier 5 — performance and scaling

Record wall-clock, CPU-hours, and cost per 1,000 genomes; a strong-scaling curve
as concurrency rises; container overhead vs the bare scripts; and the full 343k
run as the throughput headline (vs GapMind-web's ~15 s per single genome).

## What `concordance.py` reports

Per set: step-level exact-match rate + 3×3 confusion (0/1/2) + Cohen's κ
(unweighted and quadratic-weighted) + best-candidate locus agreement; pathway-level
presence (`nLo==0`) accuracy + 2×2 confusion + McNemar's test + exact
`nHi`/`nMed`/`nLo` match; a per-pathway agreement table and plot (PNG+SVG); and
`*.discordant_steps.tsv` / `*.discordant_pathways.tsv` for inspection. Genomes or
rows present on only one side are counted (they should be zero for a matched run).

## Sources

- GapMind (amino acids): Price, Deutschbauer, Arkin, *mSystems* 2020,
  https://doi.org/10.1128/mSystems.00291-20
- GapMind for carbon sources: Price et al., https://pubmed.ncbi.nlm.nih.gov/35417463/
- GapMind 2024 update: https://www.biorxiv.org/content/10.1101/2024.10.14.618325v1
- PaperBLAST / GapMind source: https://github.com/morgannprice/PaperBLAST
