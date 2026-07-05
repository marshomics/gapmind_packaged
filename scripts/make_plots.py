#!/usr/bin/env python3
"""Summary / QC / biological plots for GapMind output, aggregated over many
proteomes.

Reads the merged batch tables (<set>.sum.rules, <set>.sum.steps) plus an orgs
table (orgs.tsv or orgs.org), and writes PNG + SVG figures. SVG text is kept
editable (svg.fonttype=none) so labels can be tweaked in Illustrator/Inkscape.

Scales to hundreds of thousands of proteomes: it streams the big tables and
holds only per-genome scalars and a genome x PATHWAY (not x step) presence
matrix. For N genomes and P pathways that matrix is N*P bytes (~20 MB at
N=340,000, P=60), and the step table is consumed one row at a time.

Usage:
  make_plots.py --tables DIR --orgs orgs.tsv --sets "aa carbon" \
                --code-dir CODE_DIR --out PLOTS_DIR
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "svg.fonttype": "none",       # editable text in SVG output
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.constrained_layout.use": True,
})

GAP, PARTIAL, PRESENT = 0, 1, 2
STATUS_COLOR = {PRESENT: "#3182bd", PARTIAL: "#fdae61", GAP: "#d9d9d9"}
STATUS_LABEL = {PRESENT: "present (all steps high-conf.)",
                PARTIAL: "probably present (has medium)",
                GAP: "gap (has low-conf. step)"}
SETNAME = {"aa": "Amino-acid biosynthesis", "carbon": "Carbon catabolism"}


def setname(s):
    return SETNAME.get(s, s)


def savefig(fig, outdir, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, name + "." + ext), dpi=200)
    plt.close(fig)


def open_table(path):
    """Return (file handle positioned after header, {colname: index})."""
    f = open(path)
    hdr = f.readline().rstrip("\n").split("\t")
    return f, {h: i for i, h in enumerate(hdr)}


def load_orgs(path):
    """orgs table -> (orgId->row index, names array, nProteins array)."""
    f, idx = open_table(path)
    oi = idx["orgId"]
    ni = idx.get("genomeName", idx.get("name", oi))
    pi = idx.get("nProteins")
    o2i, names, nprot = {}, [], []
    for line in f:
        F = line.rstrip("\n").split("\t")
        oid = F[oi]
        if oid in o2i:
            continue
        o2i[oid] = len(names)
        names.append(F[ni] if ni < len(F) else oid)
        v = -1
        if pi is not None and pi < len(F) and F[pi].lstrip("-").isdigit():
            v = int(F[pi])
        nprot.append(v)
    f.close()
    return o2i, np.array(names, dtype=object), np.array(nprot, dtype=np.int64)


def load_pathway_descs(code_dir, s):
    """gaps/<set>/<set>.table -> {pathwayId: description}; empty if not found."""
    d = {}
    if not code_dir:
        return d
    p = os.path.join(code_dir, "gaps", s, s + ".table")
    if os.path.isfile(p):
        f, idx = open_table(p)
        pi = idx.get("pathwayId", 0)
        di = idx.get("desc", 1)
        for line in f:
            F = line.rstrip("\n").split("\t")
            if len(F) > max(pi, di) and F[pi] != "all":
                d[F[pi]] = F[di]
        f.close()
    return d


def aggregate(tables, s, o2i):
    """Stream <set>.sum.rules and <set>.sum.steps into small aggregates."""
    rules = os.path.join(tables, s + ".sum.rules")
    steps = os.path.join(tables, s + ".sum.steps")
    if not (os.path.isfile(rules) and os.path.isfile(steps)):
        return None
    N = len(o2i)

    # Pass 1 over rules(all): learn the pathway set (small).
    f, idx = open_table(rules)
    r = {k: idx[k] for k in ("orgId", "pathway", "rule", "nHi", "nMed", "nLo")}
    paths = {}
    for line in f:
        F = line.rstrip("\n").split("\t")
        if F[r["rule"]] != "all":
            continue
        p = F[r["pathway"]]
        if p not in paths:
            paths[p] = len(paths)
    f.close()

    # Pass 2 over rules(all): fill the genome x pathway status matrix.
    P = np.full((N, len(paths)), -1, dtype=np.int8)   # -1 = pathway not scored
    nmiss = 0
    f, idx = open_table(rules)
    for line in f:
        F = line.rstrip("\n").split("\t")
        if F[r["rule"]] != "all":
            continue
        g = o2i.get(F[r["orgId"]])
        if g is None:
            nmiss += 1
            continue
        nlo = int(F[r["nLo"]] or 0)
        nmed = int(F[r["nMed"]] or 0)
        st = PRESENT if (nlo == 0 and nmed == 0) else (PARTIAL if nlo == 0 else GAP)
        P[g, paths[F[r["pathway"]]]] = st
    f.close()

    # One pass over steps: overall confidence, per-step gap frequency,
    # per-genome high-confidence fraction (all as small arrays/dicts).
    f, idx = open_table(steps)
    st = {k: idx[k] for k in ("orgId", "pathway", "step", "score")}
    conf = np.zeros(3, dtype=np.int64)          # [low/absent, medium, high]
    stepgap = defaultdict(int)                  # (pathway, step) -> low count
    steptot = defaultdict(int)
    hi = np.zeros(N, dtype=np.int64)
    tot = np.zeros(N, dtype=np.int64)
    for line in f:
        F = line.rstrip("\n").split("\t")
        sc = F[st["score"]]
        sc = 0 if sc == "" else (int(sc) if sc in ("0", "1", "2") else None)
        if sc is None:
            continue
        conf[sc] += 1
        key = (F[st["pathway"]], F[st["step"]])
        steptot[key] += 1
        if sc == 0:
            stepgap[key] += 1
        g = o2i.get(F[st["orgId"]])
        if g is not None:
            tot[g] += 1
            if sc == 2:
                hi[g] += 1
    f.close()

    return dict(paths=paths, P=P, conf=conf, stepgap=stepgap, steptot=steptot,
                hi=hi, tot=tot, nmiss=nmiss)


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------
def plot_prevalence(agg, descs, s, N, outdir):
    paths, P = agg["paths"], agg["P"]
    inv = {v: k for k, v in paths.items()}
    npath = len(paths)
    frac = np.zeros((npath, 3))
    for j in range(npath):
        col = P[:, j]
        seen = col >= 0
        n = int(seen.sum())
        if n == 0:
            continue
        frac[j, PRESENT] = np.count_nonzero(col == PRESENT) / n
        frac[j, PARTIAL] = np.count_nonzero(col == PARTIAL) / n
        frac[j, GAP] = np.count_nonzero(col == GAP) / n
    order = np.argsort(frac[:, PRESENT])
    labels = [descs.get(inv[j], inv[j]) for j in order]
    fig, ax = plt.subplots(figsize=(8.5, max(3.0, 0.30 * npath)))
    y = np.arange(npath)
    left = np.zeros(npath)
    for stt in (PRESENT, PARTIAL, GAP):
        v = frac[order, stt]
        ax.barh(y, v, left=left, color=STATUS_COLOR[stt], label=STATUS_LABEL[stt],
                edgecolor="white", linewidth=0.3)
        left += v
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of proteomes")
    ax.set_title("%s: pathway prevalence (N=%s proteomes)" % (setname(s), format(N, ",")))
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    savefig(fig, outdir, "%s_pathway_prevalence" % s)


def plot_pathways_per_genome(agg, s, N, outdir):
    P = agg["P"]
    complete = np.count_nonzero(P == PRESENT, axis=1)
    npath = P.shape[1]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(complete, bins=np.arange(0, npath + 2) - 0.5, color="#3182bd")
    ax.set_xlabel("number of pathways fully present per proteome (of %d)" % npath)
    ax.set_ylabel("proteomes")
    ax.set_title("%s: complete pathways per proteome" % setname(s))
    med = int(np.median(complete))
    ax.axvline(med, color="k", ls="--", lw=1)
    ax.text(med, ax.get_ylim()[1] * 0.95, " median=%d" % med, fontsize=8, va="top")
    savefig(fig, outdir, "%s_pathways_per_proteome" % s)


def plot_confidence(aggs, outdir):
    sets = [s for s in aggs if aggs[s] is not None]
    if not sets:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    x = np.arange(len(sets))
    labels3 = ["high", "medium", "low / absent"]
    colors3 = [STATUS_COLOR[PRESENT], STATUS_COLOR[PARTIAL], STATUS_COLOR[GAP]]
    bottom = np.zeros(len(sets))
    for k, lab, col in [(2, labels3[0], colors3[0]), (1, labels3[1], colors3[1]),
                        (0, labels3[2], colors3[2])]:
        vals = np.array([aggs[s]["conf"][k] / max(1, aggs[s]["conf"].sum()) for s in sets])
        ax.bar(x, vals, bottom=bottom, color=col, label=lab, width=0.6)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([setname(s) for s in sets])
    ax.set_ylabel("fraction of pathway-steps")
    ax.set_ylim(0, 1)
    ax.set_title("Step-level confidence composition")
    ax.legend(fontsize=8, frameon=False, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 1.02))
    savefig(fig, outdir, "confidence_composition")


def plot_proteome_size(nprot, outdir):
    v = nprot[nprot > 0]
    if v.size == 0:
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.logspace(np.log10(max(1, v.min())), np.log10(v.max()), 60)
    ax.hist(v, bins=bins, color="#756bb1")
    ax.set_xscale("log")
    ax.set_xlabel("proteins per proteome")
    ax.set_ylabel("proteomes")
    ax.set_title("Proteome size distribution (N=%s)" % format(v.size, ","))
    med = int(np.median(v))
    ax.axvline(med, color="k", ls="--", lw=1)
    ax.text(med, ax.get_ylim()[1] * 0.95, " median=%s" % format(med, ","), fontsize=8, va="top")
    savefig(fig, outdir, "qc_proteome_size")


def plot_completeness(agg, s, outdir):
    hi, tot = agg["hi"], agg["tot"]
    m = tot > 0
    if not np.any(m):
        return
    frac = hi[m] / tot[m]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(frac, bins=np.linspace(0, 1, 51), color="#31a354")
    ax.set_xlabel("fraction of steps that are high-confidence")
    ax.set_ylabel("proteomes")
    ax.set_title("%s: per-proteome completeness (QC)" % setname(s))
    savefig(fig, outdir, "%s_qc_completeness" % s)


def plot_size_vs_complete(agg, nprot, s, outdir):
    P = agg["P"]
    complete = np.count_nonzero(P == PRESENT, axis=1)
    m = nprot > 0
    if np.count_nonzero(m) < 10:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    hb = ax.hexbin(nprot[m], complete[m], gridsize=45, xscale="log", bins="log",
                   cmap="viridis", mincnt=1)
    ax.set_xlabel("proteins per proteome (log)")
    ax.set_ylabel("complete pathways")
    ax.set_title("%s: proteome size vs complete pathways" % setname(s))
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label("log10(proteomes)")
    savefig(fig, outdir, "%s_qc_size_vs_complete" % s)


def plot_cooccurrence(agg, descs, s, outdir):
    paths, P = agg["paths"], agg["P"]
    inv = {v: k for k, v in paths.items()}
    npath = len(paths)
    if npath < 3:
        return
    B = (P == PRESENT).astype(np.float32)          # N x P
    inter = B.T @ B                                # P x P co-present counts
    cnt = np.diag(inter).copy()
    union = cnt[:, None] + cnt[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        jac = np.where(union > 0, inter / union, 0.0)
    order = np.arange(npath)
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform
        d = 1.0 - jac
        np.fill_diagonal(d, 0.0)
        order = leaves_list(linkage(squareform(d, checks=False), method="average"))
    except Exception:
        order = np.argsort(-cnt)
    J = jac[np.ix_(order, order)]
    labels = [descs.get(inv[j], inv[j]) for j in order]
    fig, ax = plt.subplots(figsize=(max(6, 0.28 * npath + 3), max(6, 0.28 * npath + 2)))
    im = ax.imshow(J, cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(range(npath))
    ax.set_yticks(range(npath))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title("%s: pathway co-occurrence (Jaccard)" % setname(s))
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Jaccard (co-present / either present)")
    savefig(fig, outdir, "%s_pathway_cooccurrence" % s)


def plot_top_gap_steps(agg, descs, s, N, outdir, topn=30):
    stepgap, steptot = agg["stepgap"], agg["steptot"]
    if not stepgap:
        return
    items = sorted(stepgap.items(), key=lambda kv: kv[1], reverse=True)[:topn]
    labels = ["%s: %s" % (descs.get(p, p), st) for (p, st), _ in items][::-1]
    frac = [c / N for (_, c) in items][::-1]
    fig, ax = plt.subplots(figsize=(8, max(3, 0.30 * len(items))))
    ax.barh(np.arange(len(items)), frac, color="#de2d26")
    ax.set_yticks(np.arange(len(items)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("fraction of proteomes where the step is a gap (low/absent)")
    ax.set_title("%s: most frequent gap steps" % setname(s))
    savefig(fig, outdir, "%s_top_gap_steps" % s)


def write_summary(aggs, nprot, N, outdir):
    with open(os.path.join(outdir, "summary_stats.tsv"), "w") as o:
        o.write("metric\tvalue\n")
        o.write("n_proteomes\t%d\n" % N)
        v = nprot[nprot > 0]
        if v.size:
            o.write("median_proteins_per_proteome\t%d\n" % int(np.median(v)))
            o.write("min_proteins\t%d\nmax_proteins\t%d\n" % (int(v.min()), int(v.max())))
        for s, agg in aggs.items():
            if agg is None:
                continue
            P = agg["P"]
            npath = P.shape[1]
            complete = np.count_nonzero(P == PRESENT, axis=1)
            prev = {inv: np.count_nonzero(P[:, j] == PRESENT) / max(1, np.count_nonzero(P[:, j] >= 0))
                    for inv, j in agg["paths"].items()}
            top = sorted(prev.items(), key=lambda kv: kv[1], reverse=True)
            o.write("%s_n_pathways\t%d\n" % (s, npath))
            o.write("%s_median_complete_pathways\t%d\n" % (s, int(np.median(complete))))
            if top:
                o.write("%s_most_prevalent\t%s (%.1f%%)\n" % (s, top[0][0], 100 * top[0][1]))
                o.write("%s_least_prevalent\t%s (%.1f%%)\n" % (s, top[-1][0], 100 * top[-1][1]))
            if agg["nmiss"]:
                o.write("%s_rows_with_unknown_orgId\t%d\n" % (s, agg["nmiss"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tables", required=True, help="dir with <set>.sum.rules/.steps")
    ap.add_argument("--orgs", required=True, help="orgs.tsv or orgs.org")
    ap.add_argument("--sets", default="aa carbon")
    ap.add_argument("--code-dir", default="", help="PaperBLAST dir (for pathway descriptions)")
    ap.add_argument("--out", required=True, help="output dir for figures")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    if not os.path.isfile(args.orgs):
        sys.exit("ERROR: orgs table not found: %s" % args.orgs)
    o2i, names, nprot = load_orgs(args.orgs)
    N = len(o2i)
    if N == 0:
        sys.exit("ERROR: no genomes in %s" % args.orgs)
    print(">> %d proteomes" % N)

    sets = args.sets.split()
    aggs = {}
    descs_by_set = {}
    for s in sets:
        descs = load_pathway_descs(args.code_dir, s)
        descs_by_set[s] = descs
        print(">> aggregating set %s ..." % s)
        agg = aggregate(args.tables, s, o2i)
        aggs[s] = agg
        if agg is None:
            print("   (no tables for %s; skipping)" % s)
            continue
        plot_prevalence(agg, descs, s, N, args.out)
        plot_pathways_per_genome(agg, s, N, args.out)
        plot_completeness(agg, s, args.out)
        plot_size_vs_complete(agg, nprot, s, args.out)
        plot_cooccurrence(agg, descs, s, args.out)
        plot_top_gap_steps(agg, descs, s, N, args.out)

    plot_confidence(aggs, args.out)
    plot_proteome_size(nprot, args.out)
    write_summary(aggs, nprot, N, args.out)
    print(">> figures (PNG+SVG) and summary_stats.tsv written to %s" % args.out)


if __name__ == "__main__":
    main()
