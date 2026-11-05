#!/usr/bin/env bash
# 04: TOOLCHAIN GATE. One full simpleaf multiplex-quant on $GATE_DATASET with the gate arguments
#   ($GATE_THREADS threads under taskset $GATE_CPUSET -- untimed, digest is thread-invariant; cr-like, frequency/frequency bc correction, decoder auto,
#   umi-edit-dist $UMI_EDIT_DIST, forwarded by multiplex-quant) and the geometry $GATE_GEOM,
#   on the step-03 index. Plain host run (no timing).
#   Checks, in an order that localizes a failure:
#     1. map_info.json num_mapped/num_poisoned == $GATE_NUM_MAPPED/$GATE_NUM_POISONED  (index + mapper)
#     2. matrix nnz/total == $GATE_NNZ/$GATE_TOTAL                                   (af reading RAD)
#     3. matrix digest (lib matrix_digest.py) == the $GATE_DIGEST_REF row            (the gate)
#   PASS iff all three hold; writes $DUR/04_validate_toolchain/gate.tsv either way.
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"
SA="$SRC/simpleaf/target/release/simpleaf"
for s in 01_toolchain 03_index; do [[ -f "$DUR/$s/DONE" ]] || fail "step $s not done"; done
ds="$GATE_DATASET"
[[ "$GATE_GEOM" == *'x[0-4]hamming('* ]] || fail "GATE_GEOM is not the x[0-4] geometry: $GATE_GEOM"
row=$(awk -F'\t' -v d="$ds" '$1==d' "$DUR/02_stage_inputs/datasets_staged.tsv")
[[ -n "$row" && -f "$SCR/02_stage_inputs/$ds.done" ]] || fail "$ds not staged by step 02"
IFS=$'\t' read -r _ _ CBC SBC PROBE R1 R2 <<<"$row"
panel=$(basename "$PROBE" | sed -nE 's/.*_(Human|Mouse)_.*_(v[0-9]+\.[0-9]+\.[0-9]+)_.*/\L\1\E_\2/p')
IDX="$SCR/index/$panel"; [[ -f "$IDX/.done" ]] || fail "index $IDX not built by step 03"
REF=$(awk -F'\t' -v d="$ds" '$1==d{print $2}' "$GATE_DIGEST_REF"); [[ -n "$REF" ]] || fail "no $ds row in $GATE_DIGEST_REF"

OUT="$SD/$ds"; LOG="$LOGD/multiplex_quant_$ds.log"
if [[ ! -f "$OUT/.done" ]]; then
  rm -rf "$OUT"; mkdir -p "$OUT"
  say "[$ds] multiplex-quant t=$GATE_THREADS cpus=$GATE_CPUSET idx=$panel geom=$GATE_GEOM"
  ( cd "$OUT" && ALEVIN_FRY_HOME="$AF_HOME" taskset -c "$GATE_CPUSET" "$SA" multiplex-quant --threads "$GATE_THREADS" \
      --resolution "$RESOLUTION" --geometry "$GATE_GEOM" --cell-bc-list "$CBC" --sample-bc-list "$SBC" \
      --index "$IDX/index" \
      --cell-bc-correction "$CELL_BC_CORRECTION" --sample-bc-correction "$SAMPLE_BC_CORRECTION" \
      --decoder "$DECODER" --umi-edit-dist "$UMI_EDIT_DIST" --reads1 "$R1" --reads2 "$R2" --output "$OUT" ) > "$LOG" 2>&1 \
    || { tail -30 "$LOG"; fail "[$ds] multiplex-quant failed (see $LOG)"; }
  date -u +%FT%TZ > "$OUT/.done"
fi

# provenance of the run: --umi-edit-dist forwarded, lz4 collate, x[0-4] geometry
grep -qE -- "quant cmd: .* --umi-edit-dist $UMI_EDIT_DIST" "$LOG" || fail "[$ds] quant was not run with --umi-edit-dist $UMI_EDIT_DIST"
grep -qE -- "collate cmd: .*--compress[= ]$COMPRESS" "$LOG" || fail "[$ds] collate lacks --compress $COMPRESS"
grep -qF -- "x[0-4]hamming(f[TTGCTAGGACCG],1)" "$LOG" || fail "[$ds] log lacks the x[0-4] geometry"
MI="$OUT/af_map/map_info.json"; M="$OUT/af_quant/alevin"
[[ -s "$MI" && -s "$M/quants_mat_rows.txt" ]] || fail "[$ds] missing map_info.json or matrix"
read -r NMAP NPOI NREADS PLAN < <("$PY" -c "
import json,sys; m=json.load(open(sys.argv[1])); p=m.get('execution_plan') or {}
print(m['num_mapped'], m['num_poisoned'], m['num_reads'], '%s/%s/%s' % (p.get('allocation',{}).get('kind'), p.get('map_threads'), p.get('decode_slots')))" "$MI")
read -r NR NC NNZ TOT DIG < <("$PY" "$REPRO/lib/matrix_digest.py" "$M")
UE=$("$PY" -c "import json,sys;print(json.load(open(sys.argv[1]))['quant_options']['crlike_umi_edit'])" "$OUT/af_quant/quant.json")
[[ "$UE" == "$UMI_EDIT_DIST" ]] || fail "[$ds] quant.json crlike_umi_edit=$UE != $UMI_EDIT_DIST (wrapper not effective)"
mkdir -p "$DD/evidence"; cp -f "$MI" "$OUT"/af_quant/*.json "$DD/evidence/"
[[ -f "$OUT/thread_policy.json" ]] && cp -f "$OUT/thread_policy.json" "$DD/evidence/"

ok(){ [[ "$1" == "$2" ]] && echo PASS || echo FAIL; }
{
  printf 'check\tobserved\texpected\tresult\n'
  printf 'num_mapped\t%s\t%s\t%s\n'   "$NMAP" "$GATE_NUM_MAPPED"   "$(ok "$NMAP" "$GATE_NUM_MAPPED")"
  printf 'num_poisoned\t%s\t%s\t%s\n' "$NPOI" "$GATE_NUM_POISONED" "$(ok "$NPOI" "$GATE_NUM_POISONED")"
  printf 'nnz\t%s\t%s\t%s\n'          "$NNZ"  "$GATE_NNZ"          "$(ok "$NNZ" "$GATE_NNZ")"
  printf 'total_umis\t%s\t%s\t%s\n'   "$TOT"  "$GATE_TOTAL"        "$(ok "$TOT" "$GATE_TOTAL")"
  printf 'digest\t%s\t%s\t%s\n'       "$DIG"  "$REF"               "$(ok "$DIG" "$REF")"
  printf 'info:threads\t%s\t-\t-\n' "$GATE_THREADS"
  printf 'info:num_reads\t%s\t-\t-\ninfo:shape\t%sx%s\t-\t-\ninfo:execution_plan\t%s\t-\t-\n' "$NREADS" "$NR" "$NC" "$PLAN"
} > "$DD/gate.tsv"
column -t -s$'\t' "$DD/gate.tsv"
grep -qP '\tFAIL$' "$DD/gate.tsv" && fail "GATE FAILED for $ds (see $DD/gate.tsv); diagnose before any later step"
say "GATE PASS: $ds digest $DIG == per-probe reference $REF"
write_done "gate.tsv	$(sha256sum "$DD/gate.tsv" | cut -c1-64)" "$ds matrix_digest	$DIG"
