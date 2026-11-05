#!/usr/bin/env bash
# 09: clustering + marker agreement between Cell Ranger, simpleaf and cyto (clustering and marker agreement, Supp. Table 4,
#     Methods "Data analysis").
#   cluster   lib/scbp_pipeline.py, once per (sample, tool) on Cell Ranger's called cells and
#             the common gene axis: MAD QC 5/3 + mito > 8% dropped, SoupX, genes in >= 20 cells, scDblFinder
#             annotate-only, scran pooling, scry deviance top 4000, PCA, 30-PC neighbours, Leiden 0.25/0.5/1.0,
#             seed 0 (parameters live in that file; see config.env). + cluster_concordance.py per resolution.
#   seedctl   the same pipeline on Cell Ranger's matrix under seed $SEEDCTL_SEED for the first
#             $SEEDCTL_MAX_SAMPLES samples of $SEEDCTL_DATASETS ("12 of 44"), + seed_control.py per resolution.
#             Runs before kmatch, because kmatch scores the control only if its object already exists.
#   kmatch    lib/kmatch_clustering.py (bisection + dense-scan fallback) to Cell
#             Ranger's k at resolution 1.0. A tool whose k cannot reach the target is listed in exceptions.tsv.
#   sweep     the same at k = $SWEEP_K.
#   deg       $REPRO/lib/deg_overlap.py: Wilcoxon top-$MARKER_TOP markers per matched cluster, Hungarian
#             cluster pairing, Jaccard + shared count.
#   assemble  lib/assemble_biol.py: Supp. Table 4, the 16plex Results numbers, marker overlap averaged BOTH over all
#             3 pairs (as Supp. Table 4) and over the 2 Cell Ranger pairs (as clustering and marker agreement), gates on counts.
# Inputs: step 06 simpleaf, $CYTO_OUT, $CR_OUT/<ds>/$CR_SUBPATH (lib/conc_common.sh). Objects (~360 MB per
#   16plex .h5ad, ~50 GB in all) go to $DUR/09_clustering/objects because steps 10b and 11 read them.
# Outputs ($DUR/09_clustering/): objects/, kmatch/, kmatch_k{4,8,12}/, deg/, tables/, DONE.
# Usage: 09_clustering.sh [--step cluster|seedctl|kmatch|sweep|deg|assemble|all]
#   cluster is 3-4 h for all 132 objects at CONC_JOBS=6: launch detached, e.g.
#   setsid nohup bash steps/09_clustering.sh --step all > $DUR/09_clustering/logs/run.log 2>&1 &
#   Every unit is resumable (.done written only after its output is verified).
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/conc_common.sh"
source "$REPRO/lib/common.sh"

STEPS="cluster seedctl kmatch sweep deg assemble"; SEL="all"
while [[ $# -gt 0 ]]; do case "$1" in
  --step) SEL="$2"; shift 2;;
  *) fail "unknown arg: $1";;
esac; done
[[ "$SEL" == all ]] && SEL="$STEPS"
for s in $SEL; do [[ " $STEPS " == *" $s "* ]] || fail "unknown --step $s"; done

for ds in $DATASETS; do
  [[ -n "${OVR_SA_QUANT_ROOT:-}" || -f "$SA_QUANT_ROOT/$ds/.done" ]] || fail "step 06 not done for $ds"
done
OBJ="$DD/objects"; mkdir -p "$SD" "$OBJ" "$DD/kmatch" "$DD/deg" "$DD/tables"
build_bench_root "$SD/bench_root" "$OBJ" || fail "bench root"
preflight || fail "preflight"
FAILS="$SD/failed.txt"; : > "$FAILS"
SCBP_VARIANT=(--n-deviant "$CLUST_N_DEVIANT"); (( CLUST_REMOVE_DOUBLETS )) && SCBP_VARIANT+=(--remove-doublets)
KM_ARGS=(--lo "$KMATCH_LO" --hi "$KMATCH_HI" --bisect-steps "$KMATCH_BISECT_STEPS"
         --scan-step "$KMATCH_SCAN_STEP" --scan-min "$KMATCH_SCAN_MIN" --fine-step "$KMATCH_FINE_STEP")

# sample list once: "ds sample" lines
PAIRS_DS_S=(); for ds in $DATASETS; do for s in $(samples "$ds"); do PAIRS_DS_S+=("$ds $s"); done; done
(( ${#PAIRS_DS_S[@]} > 0 )) || fail "no samples selected"
seedctl_samples(){   # ds -> its seed-control samples (first $SEEDCTL_MAX_SAMPLES in flexio order, if selected)
  [[ " $SEEDCTL_DATASETS " == *" $1 "* ]] || return 0
  "$PY" -c "import sys;sys.path.insert(0,'$REPRO/lib');import flexio
print('\n'.join(n for n,_ in flexio.samples_of('$1')[:$SEEDCTL_MAX_SAMPLES]))" | while read -r s; do
    if [[ -z "${OVR_SAMPLES:-}" || " $OVR_SAMPLES " == *" $1:$s "* ]]; then echo "$s"; fi
  done
}
scbp(){   # ds sample tool seed -> object + .done
  local ds="$1" s="$2" t="$3" seed="$4" tag="${1}__${2}__${3}__common"
  (( seed )) && tag+="__seed$seed"
  [[ -f "$OBJ/$tag.done" ]] && return 0
  bg_job "$tag" "$LOGD/clust_$tag.log" bash -c \
    "'$PY' '$REPRO/lib/scbp_pipeline.py' --dataset '$ds' --sample '$s' --tool '$t' --outdir '$OBJ' \
       --gene-axis common --threads '$CLUST_THREADS' --seed '$seed' ${SCBP_VARIANT[*]} \
     && test -s '$OBJ/$tag.h5ad' -a -s '$OBJ/$tag.obs.tsv.gz' -a -s '$OBJ/$tag.summary.json' \
     && date -u +%FT%TZ > '$OBJ/$tag.done'"
}

if [[ " $SEL " == *" cluster "* ]]; then
  say "cluster: ${#PAIRS_DS_S[@]} samples x 3 tools"
  for x in "${PAIRS_DS_S[@]}"; do read -r ds s <<< "$x"
    for t in cellranger simpleaf cyto; do scbp "$ds" "$s" "$t" 0; done
  done
  check_fails
  for x in "${PAIRS_DS_S[@]}"; do read -r ds s <<< "$x"
    o="$OBJ/cmp__${ds}__${s}__common.json"; [[ -f "$o.done" ]] && continue
    "$PY" "$REPRO/lib/cluster_concordance.py" --indir "$OBJ" --dataset "$ds" --sample "$s" \
      --gene-axis common --out "$o" > "$LOGD/cmp_${ds}__${s}.log" 2>&1 || fail "cluster_concordance $ds/$s"
    date -u +%FT%TZ > "$o.done"
  done
  say "cluster: done"
fi

if [[ " $SEL " == *" seedctl "* ]]; then
  for ds in $DATASETS; do for s in $(seedctl_samples "$ds"); do scbp "$ds" "$s" cellranger "$SEEDCTL_SEED"; done; done
  check_fails
  for ds in $DATASETS; do for s in $(seedctl_samples "$ds"); do
    o="$OBJ/seedctl__${ds}__${s}__common.json"; [[ -f "$o.done" ]] && continue
    "$PY" "$REPRO/lib/seed_control.py" --indir "$OBJ" --dataset "$ds" --sample "$s" --tool cellranger \
      --gene-axis common --seed "$SEEDCTL_SEED" --out "$o" > "$LOGD/seedctl_${ds}__${s}.log" 2>&1 \
      || fail "seed_control $ds/$s"
    date -u +%FT%TZ > "$o.done"
  done; done
  say "seedctl: done"
fi

kmatch(){   # outdir [target-k]
  local od="$1" k="${2:-}" x ds s
  mkdir -p "$od"
  for x in "${PAIRS_DS_S[@]}"; do read -r ds s <<< "$x"
    for t in cellranger simpleaf cyto; do
      [[ -f "$OBJ/${ds}__${s}__${t}__common.done" ]] || fail "no $t object for $ds/$s (run --step cluster)"
    done
    if [[ -z "$k" && " $(seedctl_samples "$ds" | tr '\n' ' ') " == *" $s "* ]]; then
      [[ -f "$OBJ/${ds}__${s}__cellranger__common__seed${SEEDCTL_SEED}.done" ]] \
        || fail "seed control object missing for $ds/$s (run --step seedctl before kmatch)"
    fi
    local j="$od/${ds}__${s}.json"; [[ -f "$j.done" ]] && continue
    bg_job "kmatch${k:+_k$k} $ds/$s" "$LOGD/kmatch${k:+_k$k}_${ds}__${s}.log" bash -c \
      "'$PY' '$REPRO/lib/kmatch_clustering.py' --indir '$OBJ' --dataset '$ds' --sample '$s' \
         --seed '$SEEDCTL_SEED' ${k:+--target-k $k} ${KM_ARGS[*]} --out '$j.part' \
       && mv '$j.part' '$j' && date -u +%FT%TZ > '$j.done'"
  done
  check_fails
}
if [[ " $SEL " == *" kmatch "* ]]; then kmatch "$DD/kmatch"; say "kmatch: done"; fi
if [[ " $SEL " == *" sweep "* ]]; then for K in $SWEEP_K; do kmatch "$DD/kmatch_k$K" "$K"; done; say "sweep: done"; fi

if [[ " $SEL " == *" deg "* ]]; then
  for x in "${PAIRS_DS_S[@]}"; do read -r ds s <<< "$x"
    [[ -f "$DD/kmatch/${ds}__${s}.json.done" ]] || fail "no kmatch for $ds/$s (run --step kmatch)"
    j="$DD/deg/${ds}__${s}.json"; [[ -f "$j.done" ]] && continue
    bg_job "deg $ds/$s" "$LOGD/deg_${ds}__${s}.log" bash -c \
      "'$PY' '$REPRO/lib/deg_overlap.py' --indir '$OBJ' --kmatch-dir '$DD/kmatch' --dataset '$ds' \
         --sample '$s' --top '$MARKER_TOP' --out '$j.part' \
       && mv '$j.part' '$j' && date -u +%FT%TZ > '$j.done'"
  done
  check_fails
  say "deg: done"
fi

if [[ " $SEL " == *" assemble "* ]]; then
  EXP="$SD/expected.tsv"; printf '%s\n' "${PAIRS_DS_S[@]}" | tr ' ' '\t' > "$EXP"
  SCX="$SD/expected_seedctl.tsv"; : > "$SCX"
  for ds in $DATASETS; do for s in $(seedctl_samples "$ds"); do printf '%s\t%s\n' "$ds" "$s" >> "$SCX"; done; done
  "$PY" "$REPRO/lib/assemble_biol.py" --root "$DD" --expected "$EXP" --expected-seedctl "$SCX" \
    --sweep-k "$SWEEP_K" --top "$MARKER_TOP" --out "$DD/tables" > "$LOGD/assemble.log" 2>&1 \
    || { cat "$LOGD/assemble.log" >&2; fail "assemble/gate"; }
  cat "$LOGD/assemble.log"
  if [[ -z "${OVR_SAMPLES:-}${OVR_DATASETS:-}" ]]; then
    (( ${#PAIRS_DS_S[@]} == 44 )) || fail "expected 44 samples, have ${#PAIRS_DS_S[@]}"
  fi
  LINES=()
  for f in "$DD"/tables/*; do LINES+=("tables/$(basename "$f")	$(sha256sum "$f" | cut -c1-64)"); done
  LINES+=("equal_k_exceptions	$(( $(wc -l < "$DD/tables/exceptions.tsv") - 1 ))")
  write_done "${LINES[@]}"
fi
