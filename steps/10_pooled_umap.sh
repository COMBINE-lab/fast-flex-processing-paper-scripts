#!/usr/bin/env bash
# 10: pooled-sample UMAPs of $UMAP_DATASET for simpleaf, cyto and Cell Ranger, matched cluster colors, and pooled
#     clustering concordance (ARI, NMI, homogeneity, completeness; 3 pairs) on shared cells.
#   basic  lib/pooled_umap.py --variant basic: the scanpy path of flex-paper tools/umap_probe_matrix.py
#          (normalize_total 1e4, log1p, 2000 HVG, PCA 50, 15 neighbours, UMAP, Leiden 1.0); no batch correction.
#   book   lib/pooled_umap.py --variant book: step 09's per-sample best-practice objects pooled (QC, SoupX,
#          scDblFinder, scran reused), deviance top 4000 on the pool, PCA, 30-PC neighbours; sample mixing is
#          measured (kNN sample entropy vs the cluster's expectation, k=$UMAP_MIX_K) and Harmony is applied only
#          if the mean falls below $UMAP_MIX_MIN in any tool -- the decision is in the concordance JSON.
#   Both variants then match cluster counts across tools as step 09 does (target = Cell Ranger's k at res 1.0,
#   lib/kmatch_clustering.match_k per tool, $KMATCH_* knobs); a tool that misses the target fails the step.
#   plot   lib/pooled_umap_plot.py: one 1 x 3 figure per variant (a simpleaf, b cyto, c Cell Ranger), 400 dpi.
#   This is separate from step 09's per-sample metrics, which it does not replace.
# Inputs: step 06 simpleaf, $CYTO_OUT, $CR_OUT (basic); $DUR/09_clustering/objects (book).
# Outputs ($DUR/10_pooled_umap/): <ds>_<tool>_<variant>.obs.tsv.gz (leiden, matched, umap1/2),
#   <ds>_<variant>_concordance.json (pairs; book also mixing_before/after per tool, mixing_threshold,
#   mixing_score_min, batch_corrected, batch_correction), fig_pooled_umap_<variant>.png, DONE.
# Usage: 10_pooled_umap.sh [--step basic|book|plot|all]   ~30-60 min for 16plex_human; launch detached
#   (setsid nohup ... &) for the full library.
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/conc_common.sh"
source "$REPRO/lib/common.sh"

STEPS="basic book plot"; SEL="all"   # default (all) = the paper's: $P10_VARIANTS, plot only with all three tools
while [[ $# -gt 0 ]]; do case "$1" in
  --step) SEL="$2"; shift 2;;
  *) fail "unknown arg: $1";;
esac; done
if [[ "$SEL" == all ]]; then SEL="$P10_VARIANTS"; [[ "$UMAP_TOOLS" == *simpleaf* && "$UMAP_TOOLS" == *cyto* ]] && SEL+=" plot"; fi
for s in $SEL; do [[ " $STEPS " == *" $s "* ]] || fail "unknown --step $s"; done

DS="${OVR_UMAP_DATASET:-$UMAP_DATASET}"
DATASETS="$DS"
[[ -n "${OVR_SA_QUANT_ROOT:-}" || -f "$SA_QUANT_ROOT/$DS/.done" ]] || fail "step 06 not done for $DS"
mkdir -p "$SD"
build_bench_root "$SD/bench_root" || fail "bench root"
preflight || fail "preflight"
ensure_ext_env > "$LOGD/ext_env.log" 2>&1 || { cat "$LOGD/ext_env.log" >&2; fail "ext env"; }
OBJ="$CLUST_DIR/objects"

for v in basic book; do
  [[ " $SEL " == *" $v "* ]] || continue
  J="$DD/${DS}_${v}_concordance.json"
  # skip only a finished result made by the matched-k code; an older fixed-resolution result is recomputed
  if [[ -f "$J.done" ]] && grep -q '"kmatch"' "$J" && "$PY" -c "import json,sys; sys.exit(json.load(open(sys.argv[1])).get('tools') != sys.argv[2].split())" "$J" "$UMAP_TOOLS"; then
    say "$v: done -> skip"; continue; fi
  rm -f "$J.done"
  args=(--dataset "$DS" --variant "$v" --out "$DD" --resolution 1.0 --mix-k "$UMAP_MIX_K" --tools "$UMAP_TOOLS"
        --mix-min "$UMAP_MIX_MIN" --threads "$CLUST_THREADS" --samples "$(samples "$DS" | tr '\n' ' ')"
        --km-lo "$KMATCH_LO" --km-hi "$KMATCH_HI" --km-bisect-steps "$KMATCH_BISECT_STEPS"
        --km-scan-step "$KMATCH_SCAN_STEP" --km-scan-min "$KMATCH_SCAN_MIN" --km-fine-step "$KMATCH_FINE_STEP")
  if [[ "$v" == book ]]; then
    for s in $(samples "$DS"); do for t in cellranger simpleaf cyto; do
      [[ -f "$OBJ/${DS}__${s}__${t}__common.done" ]] || fail "book: no step-09 object $DS/$s/$t in $OBJ"
    done; done
    args+=(--objects "$OBJ")
  fi
  say "$v: pooling $DS"
  "$EXT_PY" "$REPRO/lib/pooled_umap.py" "${args[@]}" > "$LOGD/$v.log" 2>&1 || { tail -30 "$LOGD/$v.log" >&2; fail "$v"; }
  "$PY" - "$J" "$DD" "$DS" "$v" <<'PYEOF' || fail "$v: output check"
import json, os, sys
import pandas as pd
j, d, ds, v = sys.argv[1:]
c = json.load(open(j))
run = c["tools"]
want = {f"{x}_vs_{y}" for x, y in (("cellranger", "simpleaf"), ("cellranger", "cyto"), ("cyto", "simpleaf")) if x in run and y in run}
assert set(c["pairs"]) == want, c["pairs"]
assert "batch_correction" in c and isinstance(c["batch_corrected"], bool) and c["n_cells"] > 0
km = c["kmatch"]
assert all(v["k"] == km["target_k"] for v in km["tools"].values()), km
assert len(set(c["n_clusters"].values())) == 1 and c["n_clusters"]["cellranger"] == km["target_k"], c["n_clusters"]
if v == "book":
    assert {"mixing_before", "mixing_threshold", "mixing_score_min"} <= set(c), sorted(c)
ref = None
for t in run:
    o = pd.read_csv(f"{d}/{ds}_{t}_{v}.obs.tsv.gz", sep="\t", index_col=0)
    assert len(o) == c["n_cells"] and o[["umap1", "umap2"]].notna().all().all(), t
    ref = o.index if ref is None else ref
    assert (o.index == ref).all(), f"{t}: cell order differs"
print(f"{v}: {c['n_cells']} cells, tools {run}, batch correction: {c['batch_correction']}; " + ", ".join(
    f"{p} ARI={m['ARI']:.3f}" for p, m in c["pairs"].items()))
PYEOF
  date -u +%FT%TZ > "$J.done"
done

if [[ " $SEL " == *" plot "* ]]; then      # optional: the 1 x 3 per-tool figure (not in the paper)
  for v in $P10_VARIANTS; do
    [[ -f "$DD/${DS}_${v}_concordance.json.done" ]] || fail "plot: run --step $v first"
    "$PY" "$REPRO/lib/pooled_umap_plot.py" --dir "$DD" --dataset "$DS" --variant "$v" \
      --out "$DD/fig_pooled_umap_$v.png" > "$LOGD/plot_$v.log" 2>&1 || { cat "$LOGD/plot_$v.log" >&2; fail "plot $v"; }
  done
fi
# DONE once every requested variant is finished (the plot is optional)
LINES=()
for v in $P10_VARIANTS; do [[ -f "$DD/${DS}_${v}_concordance.json.done" ]] || exit 0; done
for f in "$DD"/*_concordance.json "$DD"/*.obs.tsv.gz "$DD"/fig_pooled_umap_*.png; do
  [[ -e "$f" ]] && LINES+=("$(basename "$f")	$(sha256sum "$f" | cut -c1-64)")
done
write_done "${LINES[@]}"
