# GapMind pipelines

Two ways to run [GapMind](https://github.com/morgannprice/PaperBLAST) (amino-acid
biosynthesis and carbon catabolism) over your own proteomes. Both call the original
PaperBLAST/GapMind programs and produce the same analysis: per-genome pathway calls,
merged tables, presence/absence matrices, and figures. Pick whichever fits how you work.

## Which one

|                | `classic/`                                | `nextflow/`                                   |
|----------------|-------------------------------------------|-----------------------------------------------|
| Setup          | conda env + Make                          | build one container image                     |
| Runs on        | a workstation or an **SGE** cluster       | **local, SGE, or SLURM** (a profile switch)   |
| Scheduling     | local loop (`batch-local`) or SGE array jobs | Nextflow (retries, `-resume`)              |
| Reach for it   | you already have conda + SGE and want direct control | you want portability/reproducibility, or a non-SGE cluster |

Start at [`classic/README.md`](classic/README.md) for the conda + Make pipeline, or
[`nextflow/README.md`](nextflow/README.md) for the containerized workflow. They're
maintained together and cross-checked by `benchmark/`, which runs both on the same
genomes and reports concordance (Cohen's kappa, confusion matrices) —
see [`benchmark/README.md`](benchmark/README.md).

## License and citation

GPL-3.0-or-later ([`LICENSE`](LICENSE)), matching PaperBLAST/GapMind, which both
pipelines wrap rather than re-implement. [`NOTICE`](NOTICE) attributes the wrapped
software and the third-party tools. If you use either pipeline, cite the GapMind
papers (Price et al., mSystems 2020 and PLOS Genetics 2022), listed in
[`CITATION.cff`](CITATION.cff).
