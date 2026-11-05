#!/usr/bin/env bash
# 07: probe-granularity matrices from step 06's collated RAD.
#   quant    : `$AF_BIN quant --umi-edit-dist $UMI_EDIT_DIST` on $DUR/06_gene_matrices/<ds>/af_quant
#              with a probe->probe t2g (column 1 of the index t2g, duplicated), -r $RESOLUTION --use-mtx,
#              into $SD/<ds>/probe; the MatrixMarket header nnz must equal the entry-line count (truncation check).
#   validate : lib/probe_matrix.py validate (probe sums back to the step-06 gene matrix). Gate = structural only:
#              |loss| < $PROBE_LOSS_MAX, identical non-zero gene set, gene log-corr > $PROBE_LOGCORR_MIN.
#              The per-cell check cannot hold at --umi-edit-dist 1; it is recorded
#              (PASS_strict, cells_losing_over_1pct), not gated.
#   export   : lib/probe_matrix.py export -> one <sample>.probe.h5ad per sample over Cell Ranger 10.1.0's called cells
#              ($CR_OUT/<ds>/$CR_SUBPATH/<sample>/sample_filtered_feature_bc_matrix.h5), .var = $PROBE_ANNOT_DIR/<probe-set>.tsv; every sample of the dataset must be written.
#              Only for $PROBE_EXPORT_DATASETS (the datasets step 12 analyses); for the others export writes
#              .export_done reading "skipped" and no h5ad (validate and ed0 still run for all of $DATASETS).
#   320k_GEMX: ~11 min quant, ~5 min validate, ~7 min export, so run it one --step at a time.
#   ed0      : (separate run_all node, after validate) re-quant the SAME collated RAD with the real alevin-fry
#              (no --umi-edit-dist) at --umi-edit-dist 0 with BOTH the gene t2g and the probe t2g, then
#              lib/probe_matrix.py validate on that pair -> validation_ed0.json. Its total_umi_loss_frac is the paper's
#              tied-probe share (UMIs whose reads hit >1 probe of one gene, collapsed at gene level, split at probe
#              level); it shows why ed1 probe sums exceed gene sums. Gate: the FULL validate must pass at ed0 (PASS_strict: no cell > 1 % off). The raw ed0
#              matrices are deleted after validation.
#   ed1gene  : (separate run_all node, after validate) a DIRECT gene-level quant of the SAME collated RAD with the
#              gene t2g, the same binary and options as the probe quant ($AF_BIN, --umi-edit-dist $UMI_EDIT_DIST,
#              -r $RESOLUTION --use-mtx), then lib/probe_matrix.py validate of the probe matrix against it ->
#              validation_ed1gene.json. Its total_umi_loss_frac (< 0: probe totals exceed gene totals) is the
#              paper's one-mismatch comparison. Step 06's gene matrix is the sum of the probes (per-probe UMI
#              resolution, 92ff2ef), so `validate` against it measures no difference (loss 0). Gate: the
#              structural validate. The raw gene matrix is deleted after validation.
# Outputs: $DD/<ds>/{quant.json,quant.stderr.log,validation.json,validation_ed0.json,validation_ed1gene.json,export.json,<sample>.probe.h5ad}, DONE.
# Usage: 07_probe_matrices.sh [--dataset DS] [--step quant|validate|export|ed0|ed1gene|all]
#   (all = quant validate export; ed0 and ed1gene run only when asked). Resumable.
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"
source "$REPRO/lib/quant_common.sh"
source "$REPRO/lib/marker_rules.sh"

DSL="$DATASETS"; STEPS="quant validate export"
while [[ $# -gt 0 ]]; do case "$1" in
  --dataset) DSL="$2"; shift 2;;
  --step) case "$2" in quant|validate|export|ed0|ed1gene) STEPS="$2";; all) ;; *) fail "bad --step $2";; esac; shift 2;;
  *) fail "unknown arg: $1";;
esac; done
T06="$DUR/06_gene_matrices"
# flexio (imported by lib/probe_matrix.py) reads the sample list of a dataset from $FLEX_DATASETS_TSV
export FLEX_DATASETS_TSV="$DATASETS_TSV" PYTHONWARNINGS=ignore \
       NUMBA_CACHE_DIR="$SCR/tmp/numba_cache" MPLCONFIGDIR="$SCR/tmp/mpl"
mkdir -p "$SD" "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR"

for ds in $DSL; do
  Q="$T06/$ds/af_quant"; P="$SD/$ds/probe"; D="$DD/$ds"; mkdir -p "$P" "$D"
  [[ -f "$T06/$ds/.done" ]] || fail "[$ds] step 06 not done ($T06/$ds/.done)"
  index_of_ds "$ds"   # -> DS_PROBE, DS_T2G (index t2g of this dataset's panel)
  for st in $STEPS; do case "$st" in
  quant)
    if [[ -f "$P/alevin/.quant_complete" ]]; then say "[$ds] quant already complete"; continue; fi
    [[ "$(matrix_digest "$Q/alevin" | cut -f5)" == "$(cut -f5 "$T06/$ds/digest.tsv")" ]] \
      || fail "[$ds] gene matrix under $Q no longer matches step 06's digest"
    awk -F'\t' '{print $1"\t"$1}' "$DS_T2G" > "$P/probe_t2g.tsv"
    rm -rf "$P/alevin"
    say "[$ds] probe quant from $Q (t=$PROBE_THREADS)"
    t0=$SECONDS
    "$AF_BIN" quant --umi-edit-dist "$UMI_EDIT_DIST" -i "$Q" -o "$P" -m "$P/probe_t2g.tsv" -t "$PROBE_THREADS" -r "$RESOLUTION" --use-mtx \
      2> "$D/quant.stderr.log" || { tail -20 "$D/quant.stderr.log"; fail "[$ds] probe quant failed"; }
    read -r _ _ NNZ < <(grep -v '^%' "$P/alevin/quants_mat.mtx" | head -1)
    ACTUAL=$(( $(wc -l < "$P/alevin/quants_mat.mtx") - $(head -100 "$P/alevin/quants_mat.mtx" | grep -c '^%') - 1 ))
    [[ "$NNZ" == "$ACTUAL" ]] || fail "[$ds] TRUNCATED probe matrix: header nnz $NNZ vs $ACTUAL entry lines"
    UED=$("$PY" -c "import json,sys;print(json.load(open(sys.argv[1]))['quant_options']['crlike_umi_edit'])" "$P/quant.json")
    [[ "$UED" == "$UMI_EDIT_DIST" ]] || fail "[$ds] probe quant crlike_umi_edit=$UED != $UMI_EDIT_DIST"
    [[ "$(head -1 "$P/alevin/quants_mat_cols.txt")" != "$(head -1 "$T06/$ds/af_quant/alevin/quants_mat_cols.txt")" ]] \
      || fail "[$ds] probe matrix columns look like genes (t2g not probe-level?)"
    cp -f "$P/quant.json" "$D/quant.json"
    touch "$P/alevin/.quant_complete"
    say "[$ds] probe quant: $NNZ entries, $(wc -l < "$P/alevin/quants_mat_cols.txt") probes, $((SECONDS-t0)) s"
    ;;
  validate)
    [[ -f "$P/alevin/.quant_complete" ]] || fail "[$ds] probe quant not complete (run --step quant)"
    if [[ -f "$D/.validate_done" ]]; then say "[$ds] validate already done"; continue; fi
    # probe_matrix.py exits 1 when its PASS (the structural gate) is false; the JSON is written first
    "$PY" "$REPRO/lib/probe_matrix.py" validate --probe-quant "$P" --gene-quant "$Q" \
      --loss-max "$PROBE_LOSS_MAX" --logcorr-min "$PROBE_LOGCORR_MIN" --json "$D/validation.json" \
      | sed "s/^/[$ds]   /" || fail "[$ds] probe->gene validation FAILED (see $D/validation.json)"
    "$PY" -c "import json,sys; m=json.load(open(sys.argv[1]))
print('[gate] PASS', {k: m[k] for k in ('total_umi_loss_frac','gene_logcorr','genes_nonzero_probe','genes_nonzero_gene','cells_losing_over_1pct','PASS_strict')})" "$D/validation.json"
    touch "$D/.validate_done"
    ;;
  export)
    [[ -f "$D/.validate_done" ]] || fail "[$ds] validation not passed (run --step validate)"
    if export_marker_current "$D/.export_done" "$ds"; then say "[$ds] export already done"; continue; fi
    if [[ -f "$D/.export_done" ]]; then   # stale "skipped" marker: $ds was added to PROBE_EXPORT_DATASETS since
      mv "$D/.export_done" "$D/.export_done.superseded_$(date -u +%Y%m%dT%H%M%SZ)"
      say "[$ds] stale skipped-export marker moved aside: $ds is now in PROBE_EXPORT_DATASETS"
    fi
    if [[ " $PROBE_EXPORT_DATASETS " != *" $ds "* ]]; then
      say "[$ds] export skipped: not in PROBE_EXPORT_DATASETS ($PROBE_EXPORT_DATASETS)"
      echo "skipped (not in PROBE_EXPORT_DATASETS)" > "$D/.export_done"; continue
    fi
    ANNOT="$PROBE_ANNOT_DIR/$(basename "${DS_PROBE%.csv}").tsv"
    [[ -s "$ANNOT" ]] || fail "[$ds] probe annotation $ANNOT missing"
    rm -f "$D"/*.probe.h5ad "$D/export.json"
    "$PY" "$REPRO/lib/probe_matrix.py" export --probe-quant "$P" --dataset "$ds" \
      --annot "$ANNOT" --outdir "$D" --cr-dir "$CR_OUT/$ds/$CR_SUBPATH" --json "$D/export.json" | sed "s/^/[$ds]   /"
    # every sample of the dataset must have been exported (probe_matrix.py skips a sample silently)
    "$PY" -c "
import json, sys; sys.path.insert(0, sys.argv[-1]); import flexio
ds, ej = sys.argv[1], json.load(open(sys.argv[2]))
want = [s for s, _ in flexio.samples_of(ds)]
miss = [s for s in want if s not in ej]
empty = [s for s, v in ej.items() if v['n_cells'] == 0 or v['nnz'] == 0]
assert not miss and not empty, f'{ds}: samples not exported {miss}, empty {empty}'
print(f'[{ds}] export OK: {len(ej)} samples, {sum(v[\"n_cells\"] for v in ej.values())} cells')" "$ds" "$D/export.json" "$REPRO/lib"
    touch "$D/.export_done"
    ;;
  ed0)
    [[ -f "$D/.validate_done" ]] || fail "[$ds] ed1 validation not passed (run --step validate)"
    if [[ -f "$D/.ed0_done" ]]; then say "[$ds] ed0 already done"; continue; fi
    AF_BIN="$SRC/alevin-fry/target/release/alevin-fry"; E="$SD/$ds/ed0"
    rm -rf "$E"; mkdir -p "$E"
    awk -F'\t' '{print $1"\t"$1}' "$DS_T2G" > "$E/probe_t2g.tsv"
    for arm in gene probe; do
      if [[ "$arm" == probe ]]; then t2g="$E/probe_t2g.tsv"; else t2g="$DS_T2G"; fi
      say "[$ds] ed0 $arm quant (real alevin-fry, --umi-edit-dist 0, t=$PROBE_THREADS)"
      "$AF_BIN" quant --umi-edit-dist 0 -i "$Q" -o "$E/$arm" -m "$t2g" -t "$PROBE_THREADS" -r "$RESOLUTION" --use-mtx \
        2> "$D/ed0_$arm.stderr.log" || { tail -20 "$D/ed0_$arm.stderr.log"; fail "[$ds] ed0 $arm quant failed"; }
      read -r _ _ NNZ < <(grep -v '^%' "$E/$arm/alevin/quants_mat.mtx" | head -1)
      ACTUAL=$(( $(wc -l < "$E/$arm/alevin/quants_mat.mtx") - $(head -100 "$E/$arm/alevin/quants_mat.mtx" | grep -c '^%') - 1 ))
      [[ "$NNZ" == "$ACTUAL" ]] || fail "[$ds] TRUNCATED ed0 $arm matrix: header nnz $NNZ vs $ACTUAL entry lines"
      UED=$("$PY" -c "import json,sys;print(json.load(open(sys.argv[1]))['quant_options']['crlike_umi_edit'])" "$E/$arm/quant.json")
      [[ "$UED" == 0 ]] || fail "[$ds] ed0 $arm quant ran with crlike_umi_edit=$UED"
      cp -f "$E/$arm/quant.json" "$D/ed0_$arm.quant.json"
    done
    "$PY" "$REPRO/lib/probe_matrix.py" validate --probe-quant "$E/probe" --gene-quant "$E/gene" \
      --loss-max "$PROBE_LOSS_MAX" --logcorr-min "$PROBE_LOGCORR_MIN" --json "$D/validation_ed0.json" \
      | sed "s/^/[$ds ed0]   /" || fail "[$ds] ed0 probe->gene validation FAILED (see $D/validation_ed0.json)"
    "$PY" -c "import json,sys; m=json.load(open(sys.argv[1]))
assert m['PASS_strict'], 'ed0 full validate failed: cells_losing_over_1pct=%s' % m['cells_losing_over_1pct']
print('[ed0 gate] PASS_strict; tied-probe share (probe-level UMI loss at ed0) = %.5f %%' % (100 * m['total_umi_loss_frac']))" \
      "$D/validation_ed0.json" || fail "[$ds] ed0 full validate (no cell > 1 % off) failed"
    rm -rf "$E"
    touch "$D/.ed0_done"
    ;;
  ed1gene)
    [[ -f "$D/.validate_done" ]] || fail "[$ds] ed1 validation not passed (run --step validate)"
    if [[ -f "$D/.ed1gene_done" ]]; then say "[$ds] ed1gene already done"; continue; fi
    E="$SD/$ds/ed1gene"; rm -rf "$E"; mkdir -p "$E"
    say "[$ds] ed1gene: direct gene quant (--umi-edit-dist $UMI_EDIT_DIST, t=$PROBE_THREADS)"
    "$AF_BIN" quant --umi-edit-dist "$UMI_EDIT_DIST" -i "$Q" -o "$E" -m "$DS_T2G" -t "$PROBE_THREADS" -r "$RESOLUTION" --use-mtx \
      2> "$D/ed1gene.stderr.log" || { tail -20 "$D/ed1gene.stderr.log"; fail "[$ds] ed1gene quant failed"; }
    read -r _ _ NNZ < <(grep -v '^%' "$E/alevin/quants_mat.mtx" | head -1)
    ACTUAL=$(( $(wc -l < "$E/alevin/quants_mat.mtx") - $(head -100 "$E/alevin/quants_mat.mtx" | grep -c '^%') - 1 ))
    [[ "$NNZ" == "$ACTUAL" ]] || fail "[$ds] TRUNCATED ed1gene matrix: header nnz $NNZ vs $ACTUAL entry lines"
    UED=$("$PY" -c "import json,sys;print(json.load(open(sys.argv[1]))['quant_options']['crlike_umi_edit'])" "$E/quant.json")
    [[ "$UED" == "$UMI_EDIT_DIST" ]] || fail "[$ds] ed1gene quant crlike_umi_edit=$UED != $UMI_EDIT_DIST"
    [[ "$(head -1 "$E/alevin/quants_mat_cols.txt")" == "$(head -1 "$Q/alevin/quants_mat_cols.txt")" ]] \
      || fail "[$ds] ed1gene matrix columns are not step 06's genes (t2g not gene-level?)"
    cp -f "$E/quant.json" "$D/ed1gene.quant.json"
    "$PY" "$REPRO/lib/probe_matrix.py" validate --probe-quant "$P" --gene-quant "$E" \
      --loss-max "$PROBE_LOSS_MAX" --logcorr-min "$PROBE_LOGCORR_MIN" --json "$D/validation_ed1gene.json" \
      | sed "s/^/[$ds ed1gene]   /" || fail "[$ds] ed1gene probe->gene validation FAILED (see $D/validation_ed1gene.json)"
    "$PY" -c "import json,sys; m=json.load(open(sys.argv[1]))
print('[ed1gene] probe-level minus direct gene-level total UMI = %.5f %% of the gene total' % (-100 * m['total_umi_loss_frac']))" \
      "$D/validation_ed1gene.json"
    rm -rf "$E"
    touch "$D/.ed1gene_done"
    ;;
  esac; done
done

ALLDONE=1; for ds in $DATASETS; do [[ -f "$DD/$ds/.export_done" && -f "$DD/$ds/.ed0_done" && -f "$DD/$ds/.ed1gene_done" ]] || ALLDONE=0; done
if (( ALLDONE )); then
  LINES=()
  for ds in $DATASETS; do
    LINES+=("$ds/validation.json	$(sha256sum "$DD/$ds/validation.json" | cut -c1-64)")
    LINES+=("$ds/validation_ed0.json	$(sha256sum "$DD/$ds/validation_ed0.json" | cut -c1-64)")
    LINES+=("$ds/validation_ed1gene.json	$(sha256sum "$DD/$ds/validation_ed1gene.json" | cut -c1-64)")
    if [[ " $PROBE_EXPORT_DATASETS " == *" $ds "* ]]; then
      for h in "$DD/$ds"/*.probe.h5ad; do LINES+=("$ds/$(basename "$h")	$(sha256sum "$h" | cut -c1-64)"); done
    fi
  done
  write_done "${LINES[@]}"
else
  say "DONE not written until every dataset is exported, ed0- and ed1gene-validated: $DATASETS"
fi
