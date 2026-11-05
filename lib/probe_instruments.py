#!/usr/bin/env python3
"""
The probe-level tests: differential probe usage, and the pre-mRNA ratio.

Everything runs on pseudobulk raw counts aggregated by (label x sample). Size factors are
deliberately NOT applied anywhere in the tests: both instruments are within-gene ratios, so
sequencing depth cancels by construction and normalising would only add noise. Size factors
are for display.

Instrument A -- differential probe usage (the headline)
    For each probe p on gene g, a quasi-binomial GLM over pseudobulk units:
        logit(n_p / N_g) ~ is_target + sample
    fitted as Binomial with var_weights = N_g and the scale estimated by Pearson chi-square,
    which is the overdispersion. A significant `is_target` coefficient means that probe takes
    a different share of its gene's counts in the target cells -- a change in *which*
    transcript is present, not in how much of the gene there is, because the response is a
    within-gene fraction. `tx_set` then names the isoforms that probe can come from.
    Testing is stage-wise: Bonferroni over a gene's probes gives a gene-level p (valid under
    the negative within-gene dependence that fractions induce, where Simes is not), BH across
    genes screens, then BH over probes within significant genes says which probes moved.

Instrument B -- pre-mRNA ratio
    A junction probe cannot hybridise to pre-mRNA; a within-exon probe can. For a gene with
    both, restricted to pairs whose isoform sets are *identical* so that any change in isoform
    usage cancels exactly, log((exon + 0.5) / (junction + 0.5)) is monotone in the pre-mRNA
    fraction. Probe efficiency cancels only across clusters for the same pair, so the statistic
    is always a cross-cluster difference, never an absolute level.

Instrument C -- annotation discordance (repurposed; see the plan doc)
    Probes that fall in an intron under GENCODE 2024-A were expected to be a pre-mRNA-only
    panel. Measured per-probe they carry ~the same UMIs as exonic probes on the same gene, so
    they are not reporting pre-mRNA -- they mark regions 10x's design annotation called exonic
    and GENCODE no longer does. Reported as such.

Usage:
  probe_instruments.py pseudobulk --probe-root DIR --label-root DIR --dataset DS --outdir DIR [--samples S ...]
  probe_instruments.py targets    --pb DIR --min-cells N --min-samples N --out TXT
  probe_instruments.py dpu        --pb DIR --target LABEL [--reference LABEL] --out TSV
  probe_instruments.py premrna    --pb DIR --out TSV
  probe_instruments.py discord    --pb DIR --out TSV
  probe_instruments.py genede     --pb DIR --target LABEL --dpu TSV --out TSV

Every input root is an argument.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

from numpy.linalg import LinAlgError
from statsmodels.tools.sm_exceptions import PerfectSeparationError
FIT_ERRORS = (LinAlgError, PerfectSeparationError)
DPU_COLUMNS = ["probe_id", "gene_id", "gene_name", "genomic_class", "n_tx_hit", "n_tx_gene",
               "tx_set", "n_probes_gene", "n_units", "df", "usage_target", "usage_ref",
               "delta_usage", "log_or", "se", "p", "gene_umi_target", "gene_umi_ref",
               "p_gene", "q_gene", "q_probe"]

MIN_UNIT_CELLS = 20          # a (label, sample) unit below this is too small to pseudobulk
MIN_GENE_UMI = 50            # per unit, for a gene to enter a within-gene fraction test
MIN_PROBE_UMI = 20           # per side of a ratio, for instrument B

# A probe must be *observed on both sides* to be tested for a usage shift. Without this the
# quasi-binomial fit hits perfect separation whenever a probe is absent from the smaller group:
# binomial variance at a fitted fraction of 0 collapses, the standard error goes to ~0, and the
# Wald test returns p ~ 1e-9 for what is really an expected zero. A germinal-centre cluster of
# ~110 cells per sample simply cannot observe a probe used at 0.2 % of its gene, so those zeros
# carry no evidence of a usage change and are not tested.
MIN_PROBE_UMI_PER_GROUP = 10


# ------------------------------------------------------------------ pseudobulk

def probe_paths(probe_root, dataset, samples=None):
    """-> sorted [(sample, h5ad path)]; a named sample without its h5ad is an error."""
    paths = sorted(glob.glob(os.path.join(probe_root, dataset, "*.probe.h5ad")))
    if not paths:
        raise SystemExit(f"no *.probe.h5ad under {probe_root}/{dataset}")
    have = {os.path.basename(p).replace(".probe.h5ad", ""): p for p in paths}
    if samples:
        miss = [s for s in samples if s not in have]
        if miss:
            raise SystemExit(f"samples without a probe matrix: {miss}")
        return [(s, have[s]) for s in samples]
    return sorted(have.items())


DOUBLET_LOG = []   # (sample, label, cells, doublets) of every read_labels(drop_doublets=True) call


def read_labels(label_root, dataset, sample, drop_doublets=False):
    """barcode -> label (the `label` column: CellTypist call refined by the marker programs).
    drop_doublets: leave out cells with scDblFinder_class == "doublet" (config P12_REMOVE_DOUBLETS);
    the counts removed per label are appended to DOUBLET_LOG."""
    lab_path = os.path.join(label_root, dataset, f"{sample}.labels.tsv.gz")
    if not os.path.exists(lab_path):
        raise SystemExit(f"{sample}: no labels at {lab_path}")
    t = pd.read_csv(lab_path, sep="\t", index_col=0)
    # the step-11 contract: index = 16 nt cell barcode named `barcode`, a `label` column
    if t.index.name != "barcode" or "label" not in t.columns:
        raise SystemExit(f"{lab_path}: expected index `barcode` and a `label` column, got "
                         f"index {t.index.name!r}, columns {list(t.columns)[:8]}")
    if not t.index.str.fullmatch(r"[ACGTN]{16}").all():
        raise SystemExit(f"{lab_path}: barcodes are not 16 nt cores, e.g. {t.index[0]!r}")
    if t.index.has_duplicates:
        raise SystemExit(f"{lab_path}: duplicated barcodes")
    if drop_doublets:
        if "scDblFinder_class" not in t.columns:
            raise SystemExit(f"{lab_path}: no scDblFinder_class column to remove doublets with")
        bad = set(t.scDblFinder_class.astype(str)) - {"singlet", "doublet"}
        if bad:
            raise SystemExit(f"{lab_path}: unexpected scDblFinder_class values {bad}")
        dbl = t.scDblFinder_class == "doublet"
        for lb, g in t.groupby("label"):
            DOUBLET_LOG.append((sample, str(lb), len(g), int(dbl[g.index].sum())))
        t = t[~dbl]
    return t["label"].astype(str)


def write_doublet_log(path):
    """The doublets removed per sample x label, as a TSV (empty header-only when none were read)."""
    pd.DataFrame(DOUBLET_LOG, columns=["sample", "label", "cells", "doublets_removed"]).to_csv(
        path, sep="\t", index=False)


def pseudobulk(probe_root, label_root, dataset, outdir, samples=None, drop_doublets=False):
    import anndata as ad

    os.makedirs(outdir, exist_ok=True)
    rows, mats, var = [], [], None
    for sample, path in probe_paths(probe_root, dataset, samples):
        lab = read_labels(label_root, dataset, sample, drop_doublets).to_frame("label")
        A = ad.read_h5ad(path)
        if var is None:
            var = A.var.copy()
        elif not A.var.index.equals(var.index):
            raise SystemExit(f"{sample}: probe axis differs from the first sample")
        common = [b for b in A.obs_names if b in lab.index]
        if not common:
            raise SystemExit(f"{sample}: no barcode shared between matrix and labels")
        A = A[common].copy()
        labels = lab.loc[common, "label"].astype(str).values
        X = A.X.tocsr()
        for lb in pd.unique(labels):
            m = labels == lb
            if m.sum() < MIN_UNIT_CELLS:
                continue
            mats.append(np.asarray(X[m].sum(axis=0)).ravel())
            rows.append({"label": lb, "sample": sample, "n_cells": int(m.sum())})
        print(f"[pb] {sample}: {A.n_obs} cells -> "
              f"{sum(1 for r in rows if r['sample'] == sample)} units", flush=True)

    if not rows:
        raise SystemExit(f"no pseudobulk units for {dataset}")
    meta = pd.DataFrame(rows)
    M = np.vstack(mats)
    np.save(os.path.join(outdir, "probe_counts.npy"), M)
    meta.to_csv(os.path.join(outdir, "units.tsv"), sep="\t", index=False)
    if drop_doublets:
        write_doublet_log(os.path.join(outdir, "doublets_removed.tsv"))
    var.to_csv(os.path.join(outdir, "probes.tsv"), sep="\t")
    print(f"[pb] {dataset}: {M.shape[0]} units x {M.shape[1]} probes -> {outdir}",
          flush=True)
    print(meta.pivot_table(index="label", columns="sample", values="n_cells",
                           fill_value=0).to_string(), flush=True)
    return M, meta, var


def load_pb(pb):
    M = np.load(os.path.join(pb, "probe_counts.npy"))
    meta = pd.read_csv(os.path.join(pb, "units.tsv"), sep="\t")
    var = pd.read_csv(os.path.join(pb, "probes.tsv"), sep="\t", index_col=0)
    var["tx_set"] = var["tx_set"].fillna("")
    return M, meta, var


def targets(pb, min_cells, min_samples):
    """Labels present in >= min_samples samples with >= min_cells cells in total."""
    u = pd.read_csv(os.path.join(pb, "units.tsv"), sep="\t")
    g = u.groupby("label").agg(cells=("n_cells", "sum"), samples=("sample", "nunique"))
    return sorted(g[(g.cells >= min_cells) & (g.samples >= min_samples)].index)


# ------------------------------------------------------------------ instrument A

def dpu(pb, target, reference=None, min_probes=2):
    import statsmodels.api as sm

    M, meta, var = load_pb(pb)
    if target not in set(meta.label):
        raise SystemExit(f"label {target!r} not among {sorted(set(meta.label))}")
    keep = np.ones(len(meta), bool) if reference is None else \
        meta.label.isin([target, reference]).values
    M, meta = M[keep], meta[keep].reset_index(drop=True)
    is_t = (meta.label == target).values.astype(float)
    if is_t.sum() == 0 or (1 - is_t).sum() == 0:
        raise SystemExit("target or comparison group is empty after filtering")

    # sample as a blocking factor; drop the first level for identifiability
    samples = pd.get_dummies(meta["sample"], drop_first=True).astype(float).values
    X = np.column_stack([np.ones(len(meta)), is_t, samples])

    gene = var["gene_id"].values
    order = np.argsort(gene, kind="stable")
    out, n_fail = [], 0
    for gi in np.split(order, np.unique(gene[order], return_index=True)[1][1:]):
        if len(gi) < min_probes:
            continue
        sub = M[:, gi]
        N = sub.sum(axis=1)
        ok = N >= MIN_GENE_UMI
        if ok.sum() < X.shape[1] + 1 or is_t[ok].sum() == 0 or (1 - is_t[ok]).sum() == 0:
            continue
        Xo, No = X[ok], N[ok]
        for k, pi in enumerate(gi):
            y = sub[ok, k] / No
            cnt = sub[ok, k]
            in_t = is_t[ok] == 1
            if cnt[in_t].sum() < MIN_PROBE_UMI_PER_GROUP or \
               cnt[~in_t].sum() < MIN_PROBE_UMI_PER_GROUP:
                continue
            if (cnt == No).all():
                continue
            # a degenerate design (singular, perfectly separated) is a legitimate non-test;
            # it is counted and reported, never swallowed silently
            try:
                fit = sm.GLM(y, Xo, family=sm.families.Binomial(),
                             var_weights=No).fit(scale="X2")
                coef, se = fit.params[1], fit.bse[1]
            except FIT_ERRORS:
                n_fail += 1
                continue
            if not np.isfinite(se) or se == 0:
                n_fail += 1
                continue
            from scipy import stats
            dfree = max(int(ok.sum() - Xo.shape[1]), 1)
            p = 2 * stats.t.sf(abs(coef / se), dfree)
            ft = float(np.average(y[is_t[ok] == 1], weights=No[is_t[ok] == 1]))
            fr = float(np.average(y[is_t[ok] == 0], weights=No[is_t[ok] == 0]))
            out.append({
                "probe_id": var.index[pi], "gene_id": var["gene_id"].iloc[pi],
                "gene_name": var["gene_name"].iloc[pi],
                "genomic_class": var["genomic_class"].iloc[pi],
                "n_tx_hit": var["n_tx_hit"].iloc[pi],
                "n_tx_gene": var["n_tx_gene"].iloc[pi],
                "tx_set": var["tx_set"].iloc[pi],
                "n_probes_gene": len(gi), "n_units": int(ok.sum()), "df": dfree,
                "usage_target": ft, "usage_ref": fr, "delta_usage": ft - fr,
                "log_or": float(coef), "se": float(se), "p": float(p),
                "gene_umi_target": float(No[is_t[ok] == 1].sum()),
                "gene_umi_ref": float(No[is_t[ok] == 0].sum())})
    print(f"[dpu] {target}: {len(out)} probes fitted, {n_fail} degenerate fits not tested",
          flush=True)
    if not out:
        # a legitimate outcome (no probe observed on both sides); the caller writes a
        # header-only table so a completed-but-empty target is distinguishable from a crash
        return pd.DataFrame(columns=DPU_COLUMNS)
    df = pd.DataFrame(out)

    # stage-wise: gene screen by Bonferroni over the gene's probes, BH across genes
    from statsmodels.stats.multitest import multipletests
    g = df.groupby("gene_id")["p"].agg(lambda s: min(1.0, s.min() * len(s)))
    g = g.to_frame("p_gene")
    g["q_gene"] = multipletests(g.p_gene.values, method="fdr_bh")[1]
    df = df.merge(g, left_on="gene_id", right_index=True)
    df["q_probe"] = np.nan
    sig = df.q_gene < 0.05
    if sig.any():
        df.loc[sig, "q_probe"] = multipletests(df.loc[sig, "p"].values,
                                               method="fdr_bh")[1]
    return df.sort_values(["q_gene", "p"])


# ------------------------------------------------------------------ instrument B

def premrna(pb):
    """log((exon+0.5)/(junction+0.5)) per (gene, unit) over tx_set-identical pairs."""
    M, meta, var = load_pb(pb)
    v = var.reset_index().rename(columns={"index": "probe_id"})
    if "probe_id" not in v.columns:
        v = v.rename(columns={v.columns[0]: "probe_id"})
    v["idx"] = np.arange(len(v))
    j = v[v.genomic_class == "junction"]
    x = v[v.genomic_class == "within_exon"]
    pairs = j.merge(x, on=["gene_id", "tx_set"], suffixes=("_j", "_x"))
    if pairs.empty:
        raise SystemExit("no tx_set-matched junction/exon pairs")

    rows = []
    for gid, gp in pairs.groupby("gene_id"):
        ji = np.unique(gp.idx_j.values)
        xi = np.unique(gp.idx_x.values)
        Jn = M[:, ji].sum(axis=1)
        Xn = M[:, xi].sum(axis=1)
        ok = (Jn >= MIN_PROBE_UMI) & (Xn >= MIN_PROBE_UMI)
        if ok.sum() < 4:
            continue
        # Per-probe normalised, so an unequal number of probes per side does not bias it.
        # No pseudocount: MIN_PROBE_UMI already guarantees both sides are non-zero, and a
        # pseudocount added after dividing by the probe count shrinks the ratio toward 0 by
        # an amount that depends on the unit's depth. That produced a spurious ordering of
        # cell types by group size -- the smallest group sat at a log-ratio of ~0 and the
        # largest at the top -- which is a property of the estimator, not of the cells.
        lr = np.log((Xn[ok] / len(xi)) / (Jn[ok] / len(ji)))
        for u, val in zip(np.flatnonzero(ok), lr):
            rows.append({"gene_id": gid, "gene_name": gp.gene_name_x.iloc[0],
                         "label": meta.label.iloc[u], "sample": meta["sample"].iloc[u],
                         "n_junction_probes": len(ji), "n_exon_probes": len(xi),
                         "junction_umi": float(Jn[u]), "exon_umi": float(Xn[u]),
                         "log_exon_over_junction": float(val)})
    if not rows:
        raise SystemExit("no gene/unit passed the depth filter")
    return pd.DataFrame(rows)


def premrna_table(premrna_tsv, pb):
    """Per label, the mean over genes measured in EVERY label of the per-gene mean
    log(exon/junction), with the label's cells and donors from the pseudobulk units. All labels
    are written; the manuscript table chooses its rows."""
    d = pd.read_csv(premrna_tsv, sep="\t")
    g = d.groupby(["gene_id", "label"]).log_exon_over_junction.mean().unstack().dropna()
    u = pd.read_csv(os.path.join(pb, "units.tsv"), sep="\t")
    cd = u.groupby("label").agg(cells=("n_cells", "sum"), donors=("sample", "nunique"))
    tb = g.mean().rename("log_ratio").to_frame().join(cd)
    tb = tb.sort_values("log_ratio", ascending=False).reset_index().rename(
        columns={"index": "label"})
    tb["n_genes"] = g.shape[0]
    return tb, g.shape[0]


# ------------------------------------------------------------------ instrument C

def discord(pb):
    """Per-probe UMI by genomic class, and the intronic-vs-exonic comparison per gene."""
    M, meta, var = load_pb(pb)
    tot = M.sum(axis=0)
    v = var.copy()
    v["umi"] = tot
    by_class = v.groupby("genomic_class").agg(
        probes=("umi", "size"), umi=("umi", "sum"))
    by_class["pct_umi"] = 100 * by_class.umi / by_class.umi.sum()
    by_class["umi_per_probe"] = by_class.umi / by_class.probes

    per_gene = v.pivot_table(index="gene_id", columns="genomic_class", values="umi",
                             aggfunc="sum").fillna(0.0)
    n_per = v.pivot_table(index="gene_id", columns="genomic_class", values="umi",
                          aggfunc="size").fillna(0.0)
    have = [c for c in ("intronic", "within_exon") if c in per_gene.columns]
    res = None
    if len(have) == 2:
        sel = per_gene[(per_gene.intronic > 0) & (per_gene.within_exon > 0)]
        ni = n_per.loc[sel.index, "intronic"].clip(lower=1)
        nx = n_per.loc[sel.index, "within_exon"].clip(lower=1)
        res = pd.DataFrame({
            "gene_id": sel.index,
            "intronic_umi": sel.intronic.values, "exon_umi": sel.within_exon.values,
            "n_intronic_probes": ni.values, "n_exon_probes": nx.values,
            "per_probe_ratio": ((sel.intronic.values / ni.values) /
                                (sel.within_exon.values / nx.values))})
    return by_class, res


# ------------------------------------------------------------------ cli

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pseudobulk"); p.add_argument("--dataset", required=True); p.add_argument("--outdir", required=True)
    p.add_argument("--probe-root", required=True); p.add_argument("--label-root", required=True)
    p.add_argument("--samples", nargs="+")
    p.add_argument("--drop-doublets", action="store_true")
    t = sub.add_parser("targets"); t.add_argument("--pb", required=True); t.add_argument("--out", required=True)
    t.add_argument("--min-cells", type=int, required=True); t.add_argument("--min-samples", type=int, required=True)
    bt = sub.add_parser("premrna_table"); bt.add_argument("--premrna", required=True)
    bt.add_argument("--pb", required=True); bt.add_argument("--out", required=True)
    d = sub.add_parser("dpu"); d.add_argument("--pb", required=True); d.add_argument("--target", required=True)
    d.add_argument("--reference"); d.add_argument("--out", required=True)
    b = sub.add_parser("premrna"); b.add_argument("--pb", required=True); b.add_argument("--out", required=True)
    c = sub.add_parser("discord"); c.add_argument("--pb", required=True); c.add_argument("--out", required=True)
    e = sub.add_parser("genede"); e.add_argument("--pb", required=True); e.add_argument("--target", required=True)
    e.add_argument("--reference"); e.add_argument("--dpu", required=True); e.add_argument("--out", required=True)
    a = ap.parse_args()

    if a.cmd == "targets":
        lbs = targets(a.pb, a.min_cells, a.min_samples)
        open(a.out, "w").write("".join(f"{lb}\n" for lb in lbs))
        print(f"[targets] {len(lbs)} labels in >= {a.min_samples} samples with >= {a.min_cells} cells")
        return
    if a.cmd == "premrna_table":
        tb, n_genes = premrna_table(a.premrna, a.pb)
        tb.to_csv(a.out, sep="\t", index=False)
        print(f"[premrna_table] {n_genes} genes measured in all {len(tb)} labels")
        print(tb.to_string(index=False))
        return
    if a.cmd == "genede":
        de = gene_de(a.pb, a.target, a.reference)
        dpu_df = pd.read_csv(a.dpu, sep="\t")
        cls = dpu_df[dpu_df.genomic_class.isin(["junction", "within_exon"])]
        # an empty dpu (no testable probe) has no class shift; every DE gene is then concordant
        gshift = (float(abs(cls.groupby("genomic_class").delta_usage.mean().diff().iloc[-1]))
                  if cls.genomic_class.nunique() == 2 else float("nan"))
        j = join_de_dpu(de, dpu_df, gshift)
        j.to_csv(a.out, sep="\t", index=False)
        sig = j[j.q < 0.05]
        print(f"[genede] {len(de)} genes tested; {len(sig)} DE at q<0.05")
        print(f"[genede] global junction-vs-exon usage shift for this cluster = {gshift:.5f} "
              f"(read per-probe effects against this)")
        print("\n--- how the probe level explains this cluster's DE genes ---")
        print(sig.probe_verdict.value_counts().to_string())
        print("\n--- all genes with a probe-level shift, DE or not ---")
        print(j[j.n_probes_shifted > 0].probe_verdict.value_counts().to_string())
        top = j[(j.q < 0.05) & (j.probe_verdict == "isoform_switch")].head(10)
        if len(top):
            print("\n--- DE genes whose isoform usage also switched ---")
            print(top[["gene_name", "log2fc", "q", "n_probes_shifted",
                       "max_abs_delta_usage", "shifted_classes"]].to_string(index=False))
        return

    if a.cmd == "pseudobulk":
        pseudobulk(a.probe_root, a.label_root, a.dataset, a.outdir, a.samples, a.drop_doublets)
    elif a.cmd == "dpu":
        df = dpu(a.pb, a.target, a.reference)
        df.to_csv(a.out, sep="\t", index=False)
        if df.empty:
            print(f"[dpu] {a.target}: no testable probe; header-only table written")
            return
        ng = df.gene_id.nunique()
        sg = df.loc[df.q_gene < 0.05, "gene_id"].nunique()
        sp_ = int((df.q_probe < 0.05).sum())
        print(f"[dpu] {len(df)} probes on {ng} genes tested; "
              f"{sg} genes with differential probe usage (q<0.05); "
              f"{sp_} individual probes significant")
        if sg:
            top = df[df.q_probe < 0.05].head(15)
            print(top[["gene_name", "probe_id", "genomic_class", "n_tx_hit", "n_tx_gene",
                       "usage_target", "usage_ref", "delta_usage", "q_gene",
                       "q_probe"]].to_string(index=False))
    elif a.cmd == "premrna":
        df = premrna(a.pb)
        df.to_csv(a.out, sep="\t", index=False)
        print(f"[premrna] {df.gene_id.nunique()} genes x {df.label.nunique()} labels, "
              f"{len(df)} gene/unit observations")
        print(df.groupby("label")["log_exon_over_junction"].describe()[
            ["count", "mean", "50%", "std"]].to_string())
    else:
        by_class, res = discord(a.pb)
        by_class.to_csv(a.out, sep="\t")
        print("[discord] UMI by genomic class"); print(by_class.to_string())
        if res is not None:
            res.to_csv(a.out.replace(".tsv", ".per_gene.tsv"), sep="\t", index=False)
            print(f"\n[discord] {len(res)} genes with both an intronic and an exonic probe; "
                  f"median per-probe intronic/exonic ratio = "
                  f"{np.median(res.per_probe_ratio):.3f}")




# ------------------------------------------------------------------ gene-level DE + join

def gene_de(pb, target, reference=None):
    """Pseudobulk gene-level DE from the probe matrix aggregated back to genes.

    Quasi-Poisson GLM per gene, `counts ~ is_target + sample` with offset log(unit total UMI),
    which is the standard pseudobulk contrast. The gene counts come from summing the probe
    matrix, which the validation gate already showed reproduces the campaign gene matrix, so
    gene DE and probe DE are guaranteed to be the same underlying counts.
    """
    import statsmodels.api as sm
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    M, meta, var = load_pb(pb)
    keep = np.ones(len(meta), bool) if reference is None else \
        meta.label.isin([target, reference]).values
    M, meta = M[keep], meta[keep].reset_index(drop=True)
    is_t = (meta.label == target).values.astype(float)

    genes, gmat = [], []
    gid = var["gene_id"].values
    order = np.argsort(gid, kind="stable")
    for gi in np.split(order, np.unique(gid[order], return_index=True)[1][1:]):
        genes.append(gid[gi[0]])
        gmat.append(M[:, gi].sum(axis=1))
    G = np.column_stack(gmat)
    name_of = dict(zip(var["gene_id"], var["gene_name"]))

    lib = G.sum(axis=1)
    samples = pd.get_dummies(meta["sample"], drop_first=True).astype(float).values
    X = np.column_stack([np.ones(len(meta)), is_t, samples])
    off = np.log(lib)

    rows, n_fail = [], 0
    for k, g in enumerate(genes):
        y = G[:, k]
        if y.sum() < MIN_GENE_UMI or (y[is_t == 1].sum() == 0 and y[is_t == 0].sum() == 0):
            continue
        try:
            fit = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit(scale="X2")
            coef, se = fit.params[1], fit.bse[1]
        except FIT_ERRORS:
            n_fail += 1
            continue
        if not np.isfinite(se) or se == 0:
            n_fail += 1
            continue
        dfree = max(len(meta) - X.shape[1], 1)
        rows.append({"gene_id": g, "gene_name": name_of.get(g, g),
                     "log2fc": float(coef / np.log(2)), "se": float(se / np.log(2)),
                     "p": float(2 * stats.t.sf(abs(coef / se), dfree)),
                     "cpm_target": float(1e6 * y[is_t == 1].sum() / lib[is_t == 1].sum()),
                     "cpm_ref": float(1e6 * y[is_t == 0].sum() / lib[is_t == 0].sum())})
    print(f"[genede] {target}: {len(rows)} genes fitted, {n_fail} degenerate fits not tested",
          flush=True)
    df = pd.DataFrame(rows)
    df["q"] = multipletests(df.p.values, method="fdr_bh")[1]
    return df.sort_values("q")


def join_de_dpu(de, dpu_df, global_class_shift=None):
    """Classify each cluster's DE genes by what the probe level says about them.

    concordant_expression  DE, and no probe deviates from its gene's overall change: the gene
                           simply changed level, every isoform in step.
    isoform_switch         a probe changed its share of the gene AND the gene's probes map to
                           different isoform sets, so the shift names transcripts.
    splicing_state         a probe changed its share but the gene's probes all map to the same
                           isoform set (or the gene has one isoform), so the shift is
                           junction-vs-exon, i.e. spliced/unspliced composition, not isoform
                           choice.
    """
    d = dpu_df.copy()
    d["tx_set"] = d["tx_set"].fillna("")
    ntx = d.groupby("gene_id").n_tx_gene.first()
    ndistinct = d.groupby("gene_id").tx_set.nunique()
    informative = set(ntx[ntx >= 2].index) & set(ndistinct[ndistinct >= 2].index)

    sig = d[d.q_probe < 0.05]
    per = {}
    for g, gg in sig.groupby("gene_id"):
        per[g] = {
            "n_probes_shifted": int(len(gg)),
            "max_abs_delta_usage": float(gg.delta_usage.abs().max()),
            "isoform_informative": g in informative,
            "shifted_classes": ",".join(sorted(set(gg.genomic_class))),
            "shifted_tx_sets": " | ".join(sorted(set(gg.tx_set))[:3]),
        }
    # A cluster with no significant probes is a legitimate outcome, not an error: build the
    # columns explicitly so the merge and the verdict still work on an empty result.
    cols = ["n_probes_shifted", "max_abs_delta_usage", "isoform_informative",
            "shifted_classes", "shifted_tx_sets"]
    p = pd.DataFrame(per).T if per else pd.DataFrame(columns=cols)
    p = p.reindex(columns=cols)
    out = de.merge(p, left_on="gene_id", right_index=True, how="left")
    out["n_probes_shifted"] = out["n_probes_shifted"].fillna(0).astype(int)
    out["probe_verdict"] = np.where(
        out.n_probes_shifted == 0, "concordant_expression",
        np.where(out.isoform_informative.fillna(False), "isoform_switch", "splicing_state"))
    if global_class_shift is not None:
        out["delta_vs_global_shift"] = out.max_abs_delta_usage / global_class_shift
    return out


if __name__ == "__main__":
    main()
