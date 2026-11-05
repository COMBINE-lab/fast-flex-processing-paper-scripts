#!/usr/bin/env bash
# 11: cell-type labels on the NEW simpleaf gene matrices (Methods "Cell-type annotation").
#   model     CellTypist $CELLTYPIST_VERSION (venv $EXT_ENV over scbp) downloads $CELLTYPIST_MODEL into
#             $DUR/11_celltypes/celltypist_home (CELLTYPIST_FOLDER); model_check.json records its sha256, its
#             models.json entry, its feature / cell-type counts and the Methods' mouse ortholog numbers.
#   annotate  $REPRO/lib/celltype_annotate.py, unchanged: majority voting over each sample's own leiden_res1_0
#             (step 09's simpleaf object), on the step-06 simpleaf gene matrix (normalize_total 1e4 + log1p);
#             mouse symbols -> human through MGI one-to-one homology classes ($MGI_HOMOLOGY).
#   labels    lib/cell_labels.py (copy, CellTypist dir as a parameter): erythroid / plasmablast / germinal-centre B
#             / proliferating programs scored per cluster, relabel at robust z >= 10 (+ erythroid needs 5x the
#             median cluster hemoglobin fraction). Both the CellTypist call and the refined label are kept.
#   toollabels  lib/supp1_tool_labels.py: the same labelling (CellTypist majority vote + marker programs) on EACH
#             tool's step-09 object (counts layer) of $UMAP_DATASET, for Supp. Figure 1 (step 14 draws it);
#             tool_labels/<ds>.tsv.gz = cell (sample_barcode), tool, label for simpleaf, cyto and Cell Ranger.
# Inputs: step 06 simpleaf; $DUR/09_clustering/objects/<ds>__<s>__simpleaf__common.obs.tsv.gz.
# Outputs ($DUR/11_celltypes/): celltypist_home/, model_check.json, celltypist/<ds>/<s>.celltypist.tsv.gz,
#   <ds>/<s>.labels.tsv.gz + <ds>/<s>.cluster_program_z.tsv (index barcode = 16-nt core,
#   columns ... celltypist (CellTypist majority vote), marker_program, label (refined)), DONE.
# Usage: 11_celltypes.sh [--step model|annotate|labels|toollabels|all]   ~30 min for all 44 samples. Resumable per dataset.
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/conc_common.sh"
source "$REPRO/lib/common.sh"

STEPS="model annotate labels toollabels"; SEL="all"
while [[ $# -gt 0 ]]; do case "$1" in
  --step) SEL="$2"; shift 2;;
  *) fail "unknown arg: $1";;
esac; done
[[ "$SEL" == all ]] && SEL="$STEPS"
for s in $SEL; do [[ " $STEPS " == *" $s "* ]] || fail "unknown --step $s"; done

for ds in $DATASETS; do
  [[ -n "${OVR_SA_QUANT_ROOT:-}" || -f "$SA_QUANT_ROOT/$ds/.done" ]] || fail "step 06 not done for $ds"
done
OBJ="$CLUST_DIR/objects"
mkdir -p "$SD"
build_bench_root "$SD/bench_root" "$OBJ" || fail "bench root"
preflight || fail "preflight"
ensure_ext_env > "$LOGD/ext_env.log" 2>&1 || { cat "$LOGD/ext_env.log" >&2; fail "ext env"; }
export CELLTYPIST_FOLDER="$DD/celltypist_home"
MODEL_DIR="$CELLTYPIST_FOLDER/data/models"

if [[ " $SEL " == *" model "* && ! -f "$DD/model_check.json.done" ]]; then
  MOUSE_CSV=$("$PY" -c "import sys;sys.path.insert(0,'$REPRO/lib');import flexio;print(flexio.dataset_row('16k_mouse_Spleno')['probe_csv'])")
  "$EXT_PY" "$REPRO/lib/celltypist_model_check.py" --model "$CELLTYPIST_MODEL" --mgi "$MGI_HOMOLOGY" \
    --mouse-probe-csv "$MOUSE_CSV" --out "$DD/model_check.json" > "$LOGD/model.log" 2>&1 \
    || { cat "$LOGD/model.log" >&2; fail "model download/check"; }
  grep -v Warning "$LOGD/model.log" | grep -v "scanpy import" | head -20
  date -u +%FT%TZ > "$DD/model_check.json.done"
fi
[[ -f "$DD/model_check.json.done" ]] || fail "model not checked (run --step model)"
# the server always serves its latest model file: a re-download that differs from the paper's model is an error
GOT=$("$PY" -c "import json,sys;print(json.load(open(sys.argv[1]))['sha256'])" "$DD/model_check.json")
[[ "$GOT" == "$CELLTYPIST_MODEL_SHA256" ]] || fail "$CELLTYPIST_MODEL sha256 $GOT != CELLTYPIST_MODEL_SHA256 $CELLTYPIST_MODEL_SHA256"

# per-sample jobs run up to $CONC_JOBS at a time (each is ~30 s of mostly I/O); a failed job leaves
# $SD/fail.<tag> with its log path, and the step fails after the batch
par_wait(){ wait; local f; for f in "$SD"/fail.*; do [[ -e "$f" ]] || continue; tail -20 "$(cat "$f")" >&2; fail "${f##*/fail.}"; done; }
par_slot(){ while (( $(jobs -rp | wc -l) >= CONC_JOBS )); do wait -n || true; done; }
rm -f "$SD"/fail.*
ann_one(){   # ann_one DS SAMPLE: CellTypist into a private dir (celltype_annotate.py writes a fixed annotate.json)
  local ds="$1" s="$2" o="$CT/.tmp_$2" L="$LOGD/annotate_${1}__${2}.log"
  rm -rf "$o"; mkdir -p "$o"
  ( cd "$SD" && "$EXT_PY" "$REPRO/lib/celltype_annotate.py" --dataset "$ds" --sample "$s" \
      --model "$MODEL_DIR/$CELLTYPIST_MODEL" --outdir "$o" ) > "$L" 2>&1 || { echo "$L" > "$SD/fail.celltypist_${ds}__$s"; return; }
  mv "$o/$s.celltypist.tsv.gz" "$CT/"; mv "$o/annotate.json" "$CT/$s.annotate.json"; rmdir "$o"
  date -u +%FT%TZ > "$CT/$s.celltypist.tsv.gz.done"
}
lab_one(){   # lab_one DS SAMPLE
  local ds="$1" s="$2" L="$LOGD/labels_${1}__${2}.log"
  "$EXT_PY" "$REPRO/lib/cell_labels.py" --dataset "$ds" --sample "$s" --celltypist-dir "$CT" \
    --outdir "$LB" > "$L" 2>&1 || { echo "$L" > "$SD/fail.labels_${ds}__$s"; return; }
  date -u +%FT%TZ > "$LB/$s.labels.tsv.gz.done"
}
for ds in $CELLTYPE_DATASETS; do   # default: the libraries the paper shows; --full: every library
  SL=$(samples "$ds" | tr '\n' ' ')
  for s in $SL; do
    [[ -f "$OBJ/${ds}__${s}__simpleaf__common.done" ]] || fail "no step-09 simpleaf object for $ds/$s in $OBJ"
  done
  CT="$DD/celltypist/$ds"; LB="$DD/$ds"; mkdir -p "$CT" "$LB"
  if [[ " $SEL " == *" annotate "* ]]; then
    for s in $SL; do
      [[ -f "$CT/$s.celltypist.tsv.gz.done" ]] && continue
      par_slot; ann_one "$ds" "$s" &
    done
    par_wait
    for s in $SL; do grep -h "^\[annot\]" "$LOGD/annotate_${ds}__${s}.log" 2>/dev/null || true; done
  fi
  if [[ " $SEL " == *" labels "* ]]; then
    for s in $SL; do
      [[ -f "$LB/$s.labels.tsv.gz.done" ]] && continue
      [[ -f "$CT/$s.celltypist.tsv.gz.done" ]] || fail "no CellTypist output for $ds/$s (run --step annotate)"
      par_slot; lab_one "$ds" "$s" &
    done
    par_wait
    for s in $SL; do grep -h "^\[labels\]" "$LOGD/labels_${ds}__${s}.log" 2>/dev/null || true; done
  fi
done

if [[ " $SEL " == *" toollabels "* ]]; then      # Supp. Figure 1: every tool labelled the same way
  TL="$DD/tool_labels"; TD="$TL/${UMAP_DATASET:?}"; mkdir -p "$TD"
  for s in $(samples "$UMAP_DATASET"); do
    for t in simpleaf cyto cellranger; do
      o="$TD/${s}__$t.tsv.gz"; [[ -f "$o.done" ]] && continue
      [[ -f "$OBJ/${UMAP_DATASET}__${s}__${t}__common.done" ]] || fail "no step-09 $t object for $UMAP_DATASET/$s"
      par_slot
      ( L="$LOGD/toollabels_${s}__$t.log"
        "$EXT_PY" "$REPRO/lib/supp1_tool_labels.py" --object "$OBJ/${UMAP_DATASET}__${s}__${t}__common.h5ad" \
          --model "$MODEL_DIR/$CELLTYPIST_MODEL" --out "$o" > "$L" 2>&1 || { echo "$L" > "$SD/fail.toollabels_${s}__$t"; exit 0; }
        date -u +%FT%TZ > "$o.done" ) &
    done
  done
  par_wait
  "$PY" - "$TD" "$TL/$UMAP_DATASET.tsv.gz" <<'PY'
import glob, os, sys
import pandas as pd
rows = []
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.tsv.gz"))):
    s, t = os.path.basename(f)[:-7].rsplit("__", 1)
    d = pd.read_csv(f, sep="\t", index_col=0)
    rows.append(pd.DataFrame({"cell": s + "_" + d.index, "tool": t, "label": d["label"].values}))
pd.concat(rows).to_csv(sys.argv[2], sep="\t", index=False)
print(f"[toollabels] {len(rows)} sample x tool tables -> {sys.argv[2]}")
PY
fi

if [[ " $SEL " == *" labels "* ]]; then
  LINES=("model	$(grep -o '"sha256": "[0-9a-f]*"' "$DD/model_check.json" | cut -d'"' -f4)")
  for f in "$DD"/*/*.labels.tsv.gz; do LINES+=("${f#"$DD"/}	$(sha256sum "$f" | cut -c1-64)"); done
  NEXP=0; for ds in $CELLTYPE_DATASETS; do NEXP=$(( NEXP + $(samples "$ds" | wc -l) )); done
  [[ -n "${OVR_SAMPLES:-}${OVR_DATASETS:-}" ]] || (( ${#LINES[@]} == NEXP + 1 )) || fail "expected $NEXP label files, have $(( ${#LINES[@]} - 1 ))"
  T="$DD/tool_labels/$UMAP_DATASET.tsv.gz"   # Supp. Figure 1 input is part of the step's result
  [[ -s "$T" ]] || fail "no $T (run --step toollabels)"
  LINES+=("tool_labels/$UMAP_DATASET.tsv.gz	$(sha256sum "$T" | cut -c1-64)")
  write_done "${LINES[@]}"
fi
