#!/usr/bin/env python3
"""Score the concordance of two GapMind runs -- typically a reference (the
original PaperBLAST scripts, one genome at a time) versus this pipeline's output.

Both inputs are directories holding per-set tables named <set>.sum.rules and
<set>.sum.steps (the merged output of either pipeline). Genomes are matched on
orgId, which GapMind derives from the MD5 of the proteome, so identical input
sequences get the same orgId in both runs and join exactly.

Reported per set:
  step level    -- exact-match rate over (genome,pathway,step) confidence calls,
                   the 3x3 confusion matrix (0/1/2), Cohen's kappa (unweighted and
                   quadratic-weighted), and best-candidate locus agreement.
  pathway level -- presence (nLo==0) 2x2 confusion, accuracy, kappa, and McNemar's
                   test for directional bias; plus exact nHi/nMed/nLo agreement.
  per pathway   -- presence agreement per pathway (to catch a systematically off one).
Discordant rows are written out for inspection. Under matched database, search
tool, and batch size, agreement should be ~100%; anything less is a wrapper bug
or a characterized, quantified departure (diamond vs usearch; batch size).
"""
import argparse
import os
import sys

import numpy as np

SCORE = {"": 0, "0": 0, "1": 1, "2": 2}


def open_table(path):
    f = open(path)
    idx = {h: i for i, h in enumerate(f.readline().rstrip("\n").split("\t"))}
    return f, idx


def load_steps(path):
    """(orgId,pathway,step) -> (score int, best locusId)."""
    if not os.path.isfile(path):
        return {}
    f, idx = open_table(path)
    need = ("orgId", "pathway", "step", "score")
    if not all(k in idx for k in need):
        sys.exit("ERROR: %s missing columns %s" % (path, [k for k in need if k not in idx]))
    oi, pi, si, sc = (idx[k] for k in need)
    li = idx.get("locusId")
    out = {}
    for line in f:
        F = line.rstrip("\n").split("\t")
        s = SCORE.get(F[sc])
        if s is None:
            continue
        loc = F[li] if (li is not None and li < len(F)) else ""
        out[(F[oi], F[pi], F[si])] = (s, loc)
    f.close()
    return out


def load_rules(path):
    """(orgId,pathway) -> (nHi,nMed,nLo) for rule==all."""
    if not os.path.isfile(path):
        return {}
    f, idx = open_table(path)
    need = ("orgId", "pathway", "rule", "nHi", "nMed", "nLo")
    if not all(k in idx for k in need):
        sys.exit("ERROR: %s missing columns %s" % (path, [k for k in need if k not in idx]))
    oi, pi, ri, hi, me, lo = (idx[k] for k in need)
    out = {}
    for line in f:
        F = line.rstrip("\n").split("\t")
        if F[ri] != "all":
            continue
        out[(F[oi], F[pi])] = (int(F[hi] or 0), int(F[me] or 0), int(F[lo] or 0))
    f.close()
    return out


def cohen_kappa(conf, weighted=False):
    conf = conf.astype(float)
    n = conf.sum()
    if n == 0:
        return float("nan")
    po_row = conf.sum(1) / n
    po_col = conf.sum(0) / n
    k = conf.shape[0]
    if weighted:
        w = np.array([[(i - j) ** 2 for j in range(k)] for i in range(k)], float)
        wmax = (k - 1) ** 2
        w = w / wmax
        po = 1 - (w * conf / n).sum()
        pe = 1 - (w * np.outer(po_row, po_col)).sum()
    else:
        po = np.trace(conf) / n
        pe = (po_row * po_col).sum()
    return (po - pe) / (1 - pe) if (1 - pe) else float("nan")


def mcnemar(b, c):
    """b, c = the two discordant cell counts. Returns (chi2 with continuity corr, note)."""
    if b + c == 0:
        return 0.0
    return (abs(b - c) - 1) ** 2 / (b + c)


def compare_set(refdir, testdir, s, outdir, o):
    o.write("\n=== set: %s ===\n" % s)
    # ---- steps ----
    R = load_steps(os.path.join(refdir, s + ".sum.steps"))
    T = load_steps(os.path.join(testdir, s + ".sum.steps"))
    common = R.keys() & T.keys()
    only_ref = len(R.keys() - T.keys())
    only_test = len(T.keys() - R.keys())
    conf3 = np.zeros((3, 3), dtype=np.int64)
    cand_common = cand_agree = 0
    disc = []
    for k in common:
        sr, lr = R[k]
        st, lt = T[k]
        conf3[sr, st] += 1
        if lr and lt:
            cand_common += 1
            cand_agree += (lr == lt)
        if sr != st:
            disc.append((k[0], k[1], k[2], sr, st, lr, lt))
    nstep = int(conf3.sum())
    exact = int(np.trace(conf3))
    o.write("steps compared: %d  (ref-only %d, test-only %d)\n" % (nstep, only_ref, only_test))
    if nstep:
        o.write("  exact confidence match: %d/%d = %.4f\n" % (exact, nstep, exact / nstep))
        o.write("  Cohen kappa: %.4f   quadratic-weighted kappa: %.4f\n"
                % (cohen_kappa(conf3), cohen_kappa(conf3, weighted=True)))
        o.write("  confusion (rows=ref 0/1/2, cols=test 0/1/2):\n")
        for i in range(3):
            o.write("    %d: %s\n" % (i, "  ".join("%8d" % conf3[i, j] for j in range(3))))
    if cand_common:
        o.write("  best-candidate locus agreement: %d/%d = %.4f\n"
                % (cand_agree, cand_common, cand_agree / cand_common))

    # ---- pathway presence ----
    Rr = load_rules(os.path.join(refdir, s + ".sum.rules"))
    Tr = load_rules(os.path.join(testdir, s + ".sum.rules"))
    pcommon = Rr.keys() & Tr.keys()
    conf2 = np.zeros((2, 2), dtype=np.int64)      # rows=ref absent/present, cols=test
    nhml_exact = 0
    per_path = {}                                  # pathway -> [n, agree]
    pdisc = []
    for k in pcommon:
        hr, mr, lr = Rr[k]
        ht, mt, lt = Tr[k]
        pr = 1 if lr == 0 else 0
        pt = 1 if lt == 0 else 0
        conf2[pr, pt] += 1
        nhml_exact += ((hr, mr, lr) == (ht, mt, lt))
        pa = per_path.setdefault(k[1], [0, 0])
        pa[0] += 1
        pa[1] += (pr == pt)
        if pr != pt or (hr, mr, lr) != (ht, mt, lt):
            pdisc.append((k[0], k[1], hr, mr, lr, ht, mt, lt))
    npath = int(conf2.sum())
    if npath:
        acc = (conf2[0, 0] + conf2[1, 1]) / npath
        o.write("pathway calls compared: %d\n" % npath)
        o.write("  presence (nLo==0) accuracy: %.4f   kappa: %.4f\n" % (acc, cohen_kappa(conf2)))
        o.write("  nHi/nMed/nLo exact match: %d/%d = %.4f\n" % (nhml_exact, npath, nhml_exact / npath))
        b, c = int(conf2[0, 1]), int(conf2[1, 0])
        o.write("  presence confusion  ref\\test: [absent,present]\n")
        o.write("    absent : %8d %8d\n" % (conf2[0, 0], conf2[0, 1]))
        o.write("    present: %8d %8d\n" % (conf2[1, 0], conf2[1, 1]))
        o.write("  McNemar chi2 (discordant %d vs %d): %.3f%s\n"
                % (b, c, mcnemar(b, c), "" if (b + c) else " (no discordance)"))

    # ---- per-pathway agreement table + optional plot ----
    with open(os.path.join(outdir, "%s.per_pathway_agreement.tsv" % s), "w") as pf:
        pf.write("pathway\tn\tpresence_agreement\n")
        for p in sorted(per_path):
            n, ag = per_path[p]
            pf.write("%s\t%d\t%.4f\n" % (p, n, ag / n if n else float("nan")))
    _plot(per_path, s, outdir)

    # ---- discordances ----
    if disc:
        with open(os.path.join(outdir, "%s.discordant_steps.tsv" % s), "w") as df:
            df.write("orgId\tpathway\tstep\tref_score\ttest_score\tref_locus\ttest_locus\n")
            for r in disc:
                df.write("\t".join(map(str, r)) + "\n")
    if pdisc:
        with open(os.path.join(outdir, "%s.discordant_pathways.tsv" % s), "w") as df:
            df.write("orgId\tpathway\tref_nHi\tref_nMed\tref_nLo\ttest_nHi\ttest_nMed\ttest_nLo\n")
            for r in pdisc:
                df.write("\t".join(map(str, r)) + "\n")
    return dict(nstep=nstep, step_exact=exact, npath=npath,
                pres_acc=(acc if npath else float("nan")),
                only_ref=only_ref, only_test=only_test)


def _plot(per_path, s, outdir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        matplotlib.rcParams["svg.fonttype"] = "none"
        import matplotlib.pyplot as plt
    except Exception:
        return
    if not per_path:
        return
    items = sorted(per_path.items(), key=lambda kv: kv[1][1] / max(1, kv[1][0]))
    labels = [p for p, _ in items]
    vals = [ag / max(1, n) for _, (n, ag) in items]
    fig, ax = plt.subplots(figsize=(8, max(3, 0.28 * len(items))))
    ax.barh(range(len(items)), vals, color="#3182bd")
    ax.set_yticks(range(len(items)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlim(0, 1)
    ax.set_xlabel("presence-call agreement (ref vs test)")
    ax.set_title("%s: per-pathway concordance" % s)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, "%s.per_pathway_agreement.%s" % (s, ext)),
                    dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="reference output dir (<set>.sum.rules/.steps)")
    ap.add_argument("--test", required=True, help="pipeline output dir")
    ap.add_argument("--sets", default="aa carbon")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    summary = os.path.join(args.out, "concordance_summary.txt")
    with open(summary, "w") as o:
        o.write("GapMind concordance: ref=%s  test=%s\n" % (args.ref, args.test))
        overall = []
        for s in args.sets.split():
            overall.append((s, compare_set(args.ref, args.test, s, args.out, o)))
        o.write("\n=== overall ===\n")
        for s, r in overall:
            o.write("%-8s steps %d exact=%.4f | pathways %d presence-acc=%.4f | unmatched ref=%d test=%d\n"
                    % (s, r["nstep"], (r["step_exact"] / r["nstep"] if r["nstep"] else float('nan')),
                       r["npath"], r["pres_acc"], r["only_ref"], r["only_test"]))
    print(open(summary).read())
    print(">> details + discordances in %s" % args.out)


if __name__ == "__main__":
    main()
