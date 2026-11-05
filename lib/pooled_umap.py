#!/usr/bin/env python3
"""Pooled multi-sample UMAP of one library for simpleaf, cyto and Cell Ranger.

The library's samples are pooled per tool on ONE
shared cell set and ONE gene axis, every tool goes through the same steps, and the cyto / simpleaf
clusters are colored by the Cell Ranger cluster they overlap most (Hungarian assignment on the shared
cells), so equal colors mean matched clusters.

variant basic   the plain scanpy path of flex-paper tools/umap_probe_matrix.py: normalize_total 1e4, log1p,
                2,000 HVGs, PCA 50, 15 neighbours, UMAP, Leiden (igraph, 2 iterations, res 1.0, seed 0).
                Cells: Cell Ranger called AND reported by all three tools. No batch correction.
variant book    the sc-best-practices pipeline, reusing step 09's per-sample objects (QC, SoupX, scDblFinder
                annotation, scran): cells = Cell Ranger called AND QC-passing in all three tools' step-09
                objects; genes = kept in every object; normalisation = each sample's own scran layer;
                deviance top 4000 on the pooled SoupX counts; PCA (arpack); 30-PC neighbours; UMAP; Leiden 1.0.
                Integration chapter: sample mixing is MEASURED first (mean entropy of sample labels among each
                cell's k nearest PCA neighbours, divided by the same after permuting sample labels within Leiden
                clusters, i.e. 1 = as mixed as the clusters allow; see mixing()). Batch correction (Harmony on the
                PCA, key = sample) is applied only if the mean of that ratio is below --mix-min in ANY tool,
                and then to every tool alike. The decision and both scores are written to the JSON.

Cluster counts are matched across tools exactly as step 09 does per sample (lib/kmatch_clustering.match_k:
bisection + dense-scan fallback, seed 0): the target k is Cell Ranger's k at --resolution (1.0), and every
tool's final `leiden` is re-clustered on its own graph to that k. A tool that cannot reach the target is an
error (exit 1). target_k, and per tool the k at --resolution and the matched k / resolution / search method,
are written to the JSON (`kmatch`).

Usage: pooled_umap.py --dataset DS --variant basic|book --out DIR [--objects DIR] [--resolution 1.0]
                      [--mix-k 30] [--mix-min 0.8] [--threads 8] [--samples "S1 S2"]
                      [--km-lo .02 --km-hi 6 --km-bisect-steps 18 --km-scan-step .005 --km-scan-min 1e-4
                       --km-fine-step 1e-4]
"""
import argparse, json, os, sys

if "--threads" in sys.argv:                       # cap BLAS before numpy loads (as scbp_pipeline.py does)
    for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
        os.environ[_v] = sys.argv[sys.argv.index("--threads") + 1]
# rpy2 must find the scbp env's R (the Bioconductor packages), also when run from the venv over it
os.environ.setdefault("R_HOME", os.path.join(sys.base_prefix, "lib", "R"))
os.environ.setdefault("R_LIBS", os.path.join(os.environ["R_HOME"], "library"))
os.environ.setdefault("R_LIBS_USER", os.environ["R_LIBS"])

import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad, scanpy as sc
from scipy.optimize import linear_sum_assignment
from sklearn import metrics as skm
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flexio

TOOLS = ("simpleaf", "cyto", "cellranger")
SEED = 0


# ---------------------------------------------------------------- inputs

def load_basic(ds, samples):
    """-> {tool: AnnData of raw counts} on identical (sample_barcode, gene_id) axes."""
    name2id = flexio.probe_name_to_id(flexio.dataset_row(ds)["probe_csv"])
    sa, cy = flexio.simpleaf_quant_dir(ds), flexio.cyto_root(ds)
    per = {t: [] for t in TOOLS}
    for s in samples:
        L = {"cellranger": flexio.load_cellranger(flexio.cellranger_path(ds, s)),
             "simpleaf": flexio.load_alevin(sa, s), "cyto": flexio.load_cyto(cy, s)}
        for t, m in L.items():
            if m is None:
                raise SystemExit(f"{ds}/{s}: no {t} matrix")
        L = {t: (M.tocsr(), list(bc), flexio.to_gene_id(list(g), name2id)) for t, (M, bc, g) in L.items()}
        others = [set(L[t][1]) for t in ("simpleaf", "cyto")]
        bcs = [b for b in L["cellranger"][1] if all(b in o for o in others)]
        for t, (M, bc, g) in L.items():
            bi = {b: i for i, b in enumerate(bc)}
            a = ad.AnnData(M[[bi[b] for b in bcs], :], obs=pd.DataFrame(index=[f"{s}_{b}" for b in bcs]))
            a.var_names = g
            a = a[:, ~a.var_names.duplicated()].copy()
            a.obs["sample"] = s
            per[t].append(a)
        print(f"{s}: {len(bcs)} shared CR-called cells", flush=True)
    return common_axes({t: ad.concat(v) for t, v in per.items()})


def load_book(ds, samples, objects):
    """-> {tool: AnnData} from step 09's per-sample objects: X = scran log layer, layers soupX_counts."""
    per = {t: [] for t in TOOLS}
    n = {}
    for s in samples:
        objs = {t: ad.read_h5ad(os.path.join(objects, f"{ds}__{s}__{t}__common.h5ad")) for t in TOOLS}
        cells = sorted(set.intersection(*(set(o.obs_names) for o in objs.values())))
        genes = sorted(set.intersection(*(set(o.var_names) for o in objs.values())))
        for t, o in objs.items():
            a = ad.AnnData(o[cells, genes].layers["scran_normalization"].copy(),
                           obs=pd.DataFrame(index=[f"{s}_{c}" for c in cells]))
            a.var_names = genes
            a.layers["soupX_counts"] = o[cells, genes].layers["soupX_counts"].copy()
            a.obs["sample"] = s
            a.obs["scDblFinder_class"] = o.obs.loc[cells, "scDblFinder_class"].astype(str).values
            per[t].append(a)
        n[s] = len(cells)
        print(f"{s}: {len(cells)} cells QC-passing in all three tools, {len(genes)} genes", flush=True)
    return common_axes({t: ad.concat(v, join="inner") for t, v in per.items()}), n


def common_axes(pooled):
    genes = sorted(set.intersection(*(set(a.var_names) for a in pooled.values())))
    out = {t: a[:, genes].copy() for t, a in pooled.items()}
    ref = out["cellranger"].obs_names
    assert all((a.obs_names == ref).all() for a in out.values()), "tools do not share one cell order"
    return out


# ---------------------------------------------------------------- pipelines

def basic(a, res):
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000)
    sc.pp.pca(a, n_comps=50, mask_var="highly_variable", random_state=SEED)
    sc.pp.neighbors(a, n_neighbors=15, random_state=SEED)
    return a


def book_reduce(a, ro):
    from scbp_pipeline import deviance_features          # the per-sample pipeline's own implementation
    b = ad.AnnData(a.layers["soupX_counts"].tocsr(), obs=a.obs[[]].copy(), var=pd.DataFrame(index=a.var_names))
    deviance_features(ro, b)                               # top 4000 on pooled SoupX counts
    a.var["highly_variable"] = b.var["highly_deviant"].values
    a.var["binomial_deviance"] = b.var["binomial_deviance"].values
    sc.pp.pca(a, svd_solver="arpack", mask_var="highly_variable", random_state=SEED)
    return a


def mixing(a, rep, k, cluster_key):
    """Sample mixing in the kNN graph, against a perfect-mixing null.

    Per cell: Shannon entropy of the sample labels among its k nearest neighbours (30 PCs). The null is the
    same statistic after permuting sample labels WITHIN each Leiden cluster (seed 0): it keeps cell-type
    composition and the finite-k sampling bias, and removes only sample structure. Score = mean observed /
    mean null entropy (1 = as mixed as the clusters allow). The kNN inverse Simpson (iLISI-like) is reported too.
    """
    X = a.obsm[rep][:, :30]                                # the 30 PCs the neighbour graph uses
    _, ind = NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X)
    smp = pd.Categorical(a.obs["sample"]).codes
    ns = smp.max() + 1
    cl = a.obs[cluster_key].astype(str).values

    def ent(lab):
        cnt = np.zeros((a.n_obs, ns))
        np.add.at(cnt, (np.repeat(np.arange(a.n_obs), k), lab[ind[:, 1:]].ravel()), 1)
        p = cnt / k
        return -(np.where(p > 0, p * np.log(np.where(p > 0, p, 1)), 0)).sum(1), 1 / (p ** 2).sum(1)

    rng = np.random.default_rng(SEED)
    perm = smp.copy()
    for c in np.unique(cl):
        idx = np.flatnonzero(cl == c)
        perm[idx] = rng.permutation(smp[idx])
    h, simpson = ent(smp)
    h0, simpson0 = ent(perm)
    return {"mixing_score": float(h.mean() / h0.mean()), "mean_entropy": float(h.mean()),
            "mean_entropy_null": float(h0.mean()), "mean_knn_inverse_simpson": float(simpson.mean()),
            "mean_knn_inverse_simpson_null": float(simpson0.mean()), "n_samples": int(ns), "k": int(k)}


def finish(a, res, use_rep):
    if use_rep != "X_pca":
        sc.pp.neighbors(a, n_pcs=30, use_rep=use_rep, random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    sc.tl.leiden(a, resolution=res, key_added="leiden_fixed_res", flavor="igraph", n_iterations=2, random_state=SEED)


# ---------------------------------------------------------------- comparison

def match_to(ref, lab):
    """Map each cluster of `lab` to the `ref` cluster maximizing overlap (Hungarian); unmatched -> u<i>."""
    ct = pd.crosstab(pd.Series(lab, name="t"), pd.Series(ref, name="r"))
    r, c = linear_sum_assignment(-ct.values)
    m = {ct.index[i]: ct.columns[j] for i, j in zip(r, c)}
    extra = {x: f"u{i}" for i, x in enumerate(sorted(set(ct.index) - set(m), key=str))}
    return np.array([m.get(x, extra.get(x)) for x in lab])


def concordance(x, y):
    return {"ARI": float(skm.adjusted_rand_score(x, y)),
            "NMI": float(skm.normalized_mutual_info_score(x, y, average_method="arithmetic")),
            "homogeneity": float(skm.homogeneity_score(x, y)),
            "completeness": float(skm.completeness_score(x, y)),
            "n_clusters": [int(len(set(x))), int(len(set(y)))], "n_cells": int(len(x))}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True); p.add_argument("--variant", required=True, choices=("basic", "book"))
    p.add_argument("--out", required=True); p.add_argument("--objects")
    p.add_argument("--resolution", type=float, default=1.0)
    p.add_argument("--mix-k", type=int, default=30); p.add_argument("--mix-min", type=float, default=0.8)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--samples", help="space-separated subset (default: every sample of the dataset)")
    p.add_argument("--km-lo", type=float, default=0.02); p.add_argument("--km-hi", type=float, default=6.0)
    p.add_argument("--km-bisect-steps", type=int, default=18); p.add_argument("--km-scan-step", type=float, default=0.005)
    p.add_argument("--km-scan-min", type=float, default=1e-4); p.add_argument("--km-fine-step", type=float, default=1e-4)
    p.add_argument("--tools", default=" ".join(TOOLS),
                   help="tools to embed (cells are still those that pass QC in all three); with only cellranger "
                        "there is no cross-tool matching or concordance (Supp. Figure 1 needs only its embedding)")
    a = p.parse_args(); os.makedirs(a.out, exist_ok=True)
    run = [t for t in TOOLS if t in a.tools.split()]
    if "cellranger" not in run:
        raise SystemExit("--tools must include cellranger (the reference embedding and cluster count)")
    sc.settings.n_jobs = a.threads
    samples = a.samples.split() if a.samples else [s for s, _ in flexio.samples_of(a.dataset)]
    out = {"dataset": a.dataset, "samples": samples, "variant": a.variant, "resolution": a.resolution, "seed": SEED, "pairs": {}}

    if a.variant == "basic":
        pooled = load_basic(a.dataset, samples)
        for t in run:
            basic(pooled[t], a.resolution); finish(pooled[t], a.resolution, "X_pca")
        out["cells"] = "Cell Ranger called AND reported by all three tools"
        out["batch_correction"] = "none (the basic scanpy path has none)"
        out["batch_corrected"] = False
    else:
        if not a.objects:
            raise SystemExit("--variant book needs --objects (step 09 objects dir)")
        pooled, nps = load_book(a.dataset, samples, a.objects)
        from scbp_pipeline import r_setup
        ro = r_setup(a.threads)
        for t in run:
            book_reduce(pooled[t], ro)
            sc.pp.neighbors(pooled[t], n_pcs=30, random_state=SEED)
            sc.tl.leiden(pooled[t], resolution=a.resolution, key_added="_mixcl", flavor="igraph",
                         n_iterations=2, random_state=SEED)
        mix = {t: mixing(pooled[t], "X_pca", a.mix_k, "_mixcl") for t in run}
        correct = any(m["mixing_score"] < a.mix_min for m in mix.values())
        out["cells"] = "Cell Ranger called AND QC-passing in all three step-09 objects"
        out["cells_per_sample"] = nps
        out["mixing_before"] = mix
        out["mixing_threshold"] = a.mix_min
        out["mixing_score_min"] = min(m["mixing_score"] for m in mix.values())
        out["batch_corrected"] = bool(correct)
        out["mixing_rule"] = (f"correct iff mixing_score = kNN sample entropy (k={a.mix_k}, 30 PCs) / the same "
                              f"after permuting samples within Leiden {a.resolution} clusters < {a.mix_min} in any tool")
        rep = "X_pca"
        if correct:
            import harmonypy
            for t in run:
                sc.external.pp.harmony_integrate(pooled[t], "sample", basis="X_pca",
                                                 adjusted_basis="X_pca_harmony", random_state=SEED)
            rep = "X_pca_harmony"
            out["mixing_after"] = {t: mixing(pooled[t], rep, a.mix_k, "_mixcl") for t in run}
            out["batch_correction"] = f"harmony (harmonypy {harmonypy.__version__ if hasattr(harmonypy, '__version__') else ''}) on X_pca, key=sample"
        else:
            out["batch_correction"] = "none (mixing above threshold in every tool)"
        for t in run:
            finish(pooled[t], a.resolution, rep)

    # equal k across tools, as in step 09: target = Cell Ranger's k at the fixed resolution
    from kmatch_clustering import match_k
    target = int(pooled["cellranger"].obs["leiden_fixed_res"].nunique())
    km = {"target_k": target, "target_rule": f"Cell Ranger's k at resolution {a.resolution}",
          "search": dict(lo=a.km_lo, hi=a.km_hi, tries=a.km_bisect_steps, scan_step=a.km_scan_step,
                         scan_min=a.km_scan_min, fine_step=a.km_fine_step), "tools": {}}
    out["tools"] = run
    for t in run:
        lab_t, r, k, method, n = match_k(pooled[t], target, **km["search"])
        pooled[t].obs["leiden"] = pd.Categorical([str(x) for x in lab_t])
        km["tools"][t] = {"k_at_fixed_resolution": int(pooled[t].obs["leiden_fixed_res"].nunique()),
                          "k": int(k), "resolution": float(r), "method": method, "leiden_runs": int(n)}
        print(f"kmatch {t}: k={k} (target {target}) at r={r:.4f} via {method}", flush=True)
    out["kmatch"] = km
    bad = [t for t, v in km["tools"].items() if v["k"] != target]
    if bad:
        raise SystemExit(f"FATAL: target k={target} not reached for {bad}: {km['tools']}")

    ref = pooled["cellranger"]
    out["n_cells"], out["n_genes"] = int(ref.n_obs), int(ref.n_vars)
    lab = {t: pooled[t].obs["leiden"].astype(str).values for t in run}
    for x, y in (("cellranger", "simpleaf"), ("cellranger", "cyto"), ("cyto", "simpleaf")):
        if x in run and y in run:
            out["pairs"][f"{x}_vs_{y}"] = concordance(lab[x], lab[y])
    for t in run:
        o = pooled[t].obs
        o["matched"] = lab["cellranger"] if t == "cellranger" else match_to(lab["cellranger"], lab[t])
        o[["sample", "leiden", "leiden_fixed_res", "matched"]].assign(
            umap1=pooled[t].obsm["X_umap"][:, 0], umap2=pooled[t].obsm["X_umap"][:, 1]
        ).to_csv(f"{a.out}/{a.dataset}_{t}_{a.variant}.obs.tsv.gz", sep="\t")
        out.setdefault("n_clusters", {})[t] = int(len(set(lab[t])))
        out.setdefault("n_unmatched_clusters", {})[t] = int(sum(str(m).startswith("u") for m in set(o["matched"])))
    with open(f"{a.out}/{a.dataset}_{a.variant}_concordance.json", "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
