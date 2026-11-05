#!/usr/bin/env bash
# 08: per-sample gene-count concordance, 44 samples x 3 tool pairs (count concordance).
#   Metric set = the alevin-fry paper's (mean/median per-cell Spearman, MARD drop-NA and NA=0, rFP/rFN per cell)
#   + count concordance % (1 - L1/total), per-gene log1p Pearson / Spearman, total-UMI ratio test/ref.
#   Code: $REPRO/lib/compare_consistency.py, unchanged. Cell set = Cell Ranger's per-sample FILTERED barcodes,
#   applied to every tool; genes matched on gene_id; raw counts, no normalisation.
#   In simpleaf_vs_cyto, simpleaf is the reference (UMI ratio = cyto/simpleaf); see pair_definitions.tsv.
#
# Inputs: simpleaf = step 06 ($SA_QUANT_ROOT/<ds>/$SA_QUANT_SUBDIR); cyto = $CYTO_OUT/<ds>/cyto;
#   Cell Ranger = $CR_OUT/<ds>/$CR_SUBPATH/<sample>/ (step 05b rep1).
#   This step reads the
#   version from every h5 at run time and writes it to cr_provenance.tsv; with the config default it is
#   10.1.0 everywhere.
# Outputs ($DUR/08_concordance_counts/): per_sample/<ds>__<sample>.json (+ .done), per_sample.tsv,
#   table_s1_concordance_detail.tsv, pair_definitions.tsv, text_ranges.json, cr_provenance.tsv, DONE.
# Usage: 08_concordance_counts.sh        ~15 min at CONC_JOBS=6 (320k_GEMX dominates). Resumable per sample.
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/conc_common.sh"
source "$REPRO/lib/common.sh"

for ds in $DATASETS; do
  [[ -n "${OVR_SA_QUANT_ROOT:-}" || -f "$SA_QUANT_ROOT/$ds/.done" ]] || fail "step 06 not done for $ds ($SA_QUANT_ROOT/$ds/.done)"
done
mkdir -p "$SD" "$DD/per_sample"
build_bench_root "$SD/bench_root" || fail "bench root"
preflight || fail "preflight: an arm does not resolve to the configured inputs"
say "inputs: simpleaf=$SA_QUANT_ROOT cyto=$CYTO_OUT cellranger=$CR_OUT/<ds>/$CR_SUBPATH"

FAILS="$SD/failed.txt"; : > "$FAILS"; EXP="$SD/expected.tsv"; : > "$EXP"
for ds in $DATASETS; do
  for s in $(samples "$ds"); do
    printf '%s\t%s\n' "$ds" "$s" >> "$EXP"
    j="$DD/per_sample/${ds}__${s}.json"
    [[ -f "$j.done" ]] && continue
    rm -f "$j"
    bg_job "$ds/$s" "$LOGD/counts_${ds}__${s}.log" bash -c \
      "'$PY' '$REPRO/lib/compare_consistency.py' --dataset '$ds' --sample '$s' --out '$j.part' \
       && mv '$j.part' '$j' && date -u +%FT%TZ > '$j.done'"
  done
done
check_fails
(( $(wc -l < "$EXP") > 0 )) || fail "no samples selected"

"$PY" "$REPRO/lib/conc_tables.py" --jsons "$DD/per_sample" --expected "$EXP" --out "$DD" \
  > "$LOGD/tables.log" 2>&1 || { cat "$LOGD/tables.log" >&2; fail "tables/gate"; }
cat "$LOGD/tables.log"

if [[ -z "${OVR_SAMPLES:-}${OVR_DATASETS:-}" ]]; then
  n=$(wc -l < "$EXP"); (( n == 44 )) || fail "expected 44 samples, have $n"
fi
LINES=()
for f in per_sample.tsv table_s1_concordance_detail.tsv pair_definitions.tsv text_ranges.json cr_provenance.tsv; do
  LINES+=("$f	$(sha256sum "$DD/$f" | cut -c1-64)")
done
LINES+=("cellranger_versions	$(cut -f1,4 "$DD/cr_provenance.tsv" | tail -n +2 | sort -u | tr '\t\n' ': ' )")
write_done "${LINES[@]}"
