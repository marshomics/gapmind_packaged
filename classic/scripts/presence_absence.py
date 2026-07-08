#!/usr/bin/env python3
"""Coerce GapMind pathway calls into a genome x pathway presence/absence matrix.

GapMind does not emit a clean present/absent flag. For each pathway in each
genome it reports, on the best-scoring path, how many steps are high- (nHi),
medium- (nMed) and low-confidence (nLo). We threshold using GapMind's OWN
operational definition rather than an arbitrary cutoff:

  * step confidence (high/medium/low) is decided by gapsummary.pl from percent
    identity, coverage, and bit score vs the best reverse hit;
  * a pathway with nLo == 0 (every step at least medium confidence) is called
    "probably present" -- this is exactly how GapMind reports pathway presence;
  * a pathway with nMed == 0 and nLo == 0 (every step high confidence) is a
    "high-confidence" call;
  * nLo > 0 (at least one low-confidence / likely-missing step) is called absent.

Two thresholds are offered because no single cutoff is perfect -- medium steps
and curated "known gaps" are genuine grey zones:
  --mode probably  (default)  present = 1 iff nLo == 0        (GapMind's call)
  --mode strict               present = 1 iff nMed == 0 and nLo == 0

Outputs, per set:
  <set>.presence.tsv     genome x pathway 0/1 at the chosen threshold
  <set>.confidence.tsv   genome x pathway 2=high / 1=probably / 0=absent
  <set>.pathways.tsv     pathwayId  ->  description  (column legend)

The 2/1/0 confidence table preserves the full call so you can re-threshold (or
build your own rule) without recomputing. Memory-bounded: one int8 per
genome x pathway cell (~27 MB at 340,000 genomes x ~80 pathways).
"""
import argparse
import os
import sys

import numpy as np

HIGH, PROBABLY, ABSENT = 2, 1, 0


def open_table(path):
    f = open(path)
    hdr = f.readline().rstrip("\n").split("\t")
    return f, {h: i for i, h in enumerate(hdr)}


def load_orgs(path):
    f, idx = open_table(path)
    oi = idx["orgId"]
    ni = idx.get("genomeName", idx.get("name", oi))
    o2i, names = {}, []
    for line in f:
        F = line.rstrip("\n").split("\t")
        oid = F[oi]
        if oid not in o2i:
            o2i[oid] = len(names)
            names.append(F[ni] if ni < len(F) else oid)
    f.close()
    return o2i, names


def load_pathways(code_dir, s, rules):
    """Ordered pathway list + descriptions. Prefer gaps/<set>/<set>.table (the
    canonical, complete list); otherwise discover from the data."""
    order, desc = [], {}
    tbl = os.path.join(code_dir, "gaps", s, s + ".table") if code_dir else ""
    if tbl and os.path.isfile(tbl):
        f, idx = open_table(tbl)
        pi = idx.get("pathwayId", 0)
        di = idx.get("desc", 1)
        for line in f:
            F = line.rstrip("\n").split("\t")
            if len(F) > pi and F[pi] and F[pi] != "all":
                order.append(F[pi])
                desc[F[pi]] = F[di] if len(F) > di else F[pi]
        f.close()
        if order:
            return order, desc
    # fallback: discover from rule==all rows
    f, idx = open_table(rules)
    ri, rr = idx["pathway"], idx["rule"]
    seen = set()
    for line in f:
        F = line.rstrip("\n").split("\t")
        if F[rr] == "all" and F[ri] not in seen:
            seen.add(F[ri])
            order.append(F[ri])
            desc[F[ri]] = F[ri]
    f.close()
    return order, desc


def build_matrix(rules, o2i, porder):
    """N x P int8 of confidence codes (2/1/0); -1 = pathway not reported."""
    p2i = {p: j for j, p in enumerate(porder)}
    C = np.full((len(o2i), len(porder)), -1, dtype=np.int8)
    f, idx = open_table(rules)
    need = ("orgId", "pathway", "rule", "nHi", "nMed", "nLo")
    if not all(k in idx for k in need):
        sys.exit("ERROR: %s missing columns %s" % (rules, [k for k in need if k not in idx]))
    ci = {k: idx[k] for k in need}
    nmiss_g = nmiss_p = 0
    for line in f:
        F = line.rstrip("\n").split("\t")
        if F[ci["rule"]] != "all":
            continue
        g = o2i.get(F[ci["orgId"]])
        if g is None:
            nmiss_g += 1
            continue
        j = p2i.get(F[ci["pathway"]])
        if j is None:
            nmiss_p += 1
            continue
        nmed = int(F[ci["nMed"]] or 0)
        nlo = int(F[ci["nLo"]] or 0)
        C[g, j] = HIGH if (nmed == 0 and nlo == 0) else (PROBABLY if nlo == 0 else ABSENT)
    f.close()
    return C, nmiss_g, nmiss_p


def write_matrix(path, names, porder, rows_int, missing_fill):
    with open(path, "w") as o:
        o.write("genome\t" + "\t".join(porder) + "\n")
        for i, nm in enumerate(names):
            vals = rows_int[i]
            o.write(nm + "\t" + "\t".join(missing_fill if v < 0 else str(v) for v in vals) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", required=True, help="dir with <set>.sum.rules")
    ap.add_argument("--orgs", required=True, help="orgs.tsv or orgs.org")
    ap.add_argument("--sets", default="aa carbon")
    ap.add_argument("--code-dir", default="", help="PaperBLAST dir (for the canonical pathway list)")
    ap.add_argument("--mode", choices=["probably", "strict"], default="probably")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    o2i, names = load_orgs(args.orgs)
    N = len(o2i)
    if N == 0:
        sys.exit("ERROR: no genomes in %s" % args.orgs)

    thr = PROBABLY if args.mode == "probably" else HIGH
    print(">> %d genomes; mode=%s (present iff confidence >= %d)" % (N, args.mode, thr))

    for s in args.sets.split():
        rules = os.path.join(args.tables, s + ".sum.rules")
        if not os.path.isfile(rules):
            print("   (%s.sum.rules not found; skipping %s)" % (s, s))
            continue
        porder, desc = load_pathways(args.code_dir, s, rules)
        C, mg, mp = build_matrix(rules, o2i, porder)

        # presence (0/1 at threshold); missing pathway -> 0 (not called = absent)
        pres = np.where(C < 0, 0, (C >= thr).astype(np.int8))
        write_matrix(os.path.join(args.out, s + ".presence.tsv"), names, porder, pres, "0")
        # confidence (2/1/0); missing -> blank
        write_matrix(os.path.join(args.out, s + ".confidence.tsv"), names, porder, C, "")
        with open(os.path.join(args.out, s + ".pathways.tsv"), "w") as o:
            o.write("pathwayId\tdescription\n")
            for p in porder:
                o.write("%s\t%s\n" % (p, desc.get(p, p)))

        called = int(np.count_nonzero(pres.sum(axis=1) >= 0))
        per_genome = pres.sum(axis=1)
        print("   %s: %d pathways x %d genomes -> %s.presence.tsv "
              "(median present/genome = %d; %d unknown orgIds, %d off-list pathway rows)"
              % (s, len(porder), N, s, int(np.median(per_genome)), mg, mp))

    print(">> wrote presence/confidence matrices to %s" % args.out)


if __name__ == "__main__":
    main()
