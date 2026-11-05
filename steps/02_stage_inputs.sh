#!/usr/bin/env bash
# 02: stage inputs from $DATA_ROOT (network storage) onto the local NVMe.
#   - GEX R1/R2 FASTQs ($NONGEX_RE excluded) from $FASTQ_SRC/<ds> -> $SCR/fastq/<ds>/<same relative path>
#   - cell-bc list, sample-bc list and probe CSV (columns 6/7/8 of $DATASETS_TSV) -> $SCR/refs/<ds>/
#   - every file byte-size verified against its source; per-dataset .done marker after verification
#   - writes $DUR/02_stage_inputs/datasets_staged.tsv (the $SCR paths later steps read)
# Usage: 02_stage_inputs.sh [--dataset DS]   (default: all of $DATASETS). Resumable.
#   $STAGE_STREAMS cp streams per dataset; several --dataset runs may go concurrently (flock on the TSV).
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"

DSL="$DATASETS"
while [[ $# -gt 0 ]]; do case "$1" in
  --dataset) DSL="$2"; shift 2;;
  *) fail "unknown arg: $1";;
esac; done
mkdir -p "$SD"
STAGED="$DD/datasets_staged.tsv"
[[ -s "$STAGED" ]] || printf 'dataset\tgeometry_tsv\tcell_bc_list\tsample_bc_list\tprobe_csv\treads1\treads2\n' > "$STAGED"
tsvcol(){ awk -F'\t' -v d="$1" -v c="$2" '!/^#/ && $1==d {print $c}' "$DATASETS_TSV"; }

# copy_verified SRC DST: copy via .part unless DST already has SRC's byte size; then verify
copy_verified(){
  local s="$1" t="$2"
  if [[ ! -f "$t" || "$(stat -c %s "$t")" != "$(stat -Lc %s "$s")" ]]; then
    mkdir -p "$(dirname "$t")"; cp -L "$s" "$t.part"; mv "$t.part" "$t"
  fi
  [[ "$(stat -c %s "$t")" == "$(stat -Lc %s "$s")" ]] || { echo "size mismatch: $t vs $s" >&2; return 1; }
}
export -f copy_verified

for ds in $DSL; do
  src="$FASTQ_SRC/$ds"; fq="$SCR/fastq/$ds"; rf="$SCR/refs/$ds"
  [[ -n "$(tsvcol "$ds" 1)" ]] || fail "$ds not in $DATASETS_TSV"
  [[ "$(tsvcol "$ds" 12)" == "$src" ]] || fail "$ds: TSV fastq_dir '$(tsvcol "$ds" 12)' != $src"
  mapfile -t FQS < <(find -L "$src" -name '*_R[12]_*.fastq.gz' | grep -ivE "$NONGEX_RE" | sort)
  (( ${#FQS[@]} > 0 )) || fail "$ds: no GEX R1/R2 under $src"
  if [[ -f "$SD/$ds.done" ]]; then
    say "[$ds] already staged -> re-verify sizes only"
  else
    say "[$ds] staging ${#FQS[@]} GEX FASTQs ($(du -chL "${FQS[@]}" | tail -1 | cut -f1)) -> $fq"
  fi
  for f in "${FQS[@]}"; do printf '%s\0%s\0' "$f" "$fq/${f#$src/}"; done \
    | xargs -0 -n2 -P"$STAGE_STREAMS" bash -c 'copy_verified "$1" "$2"' _ || fail "[$ds] FASTQ copy/verify failed"
  for c in 6 7 8; do
    s=$(tsvcol "$ds" "$c"); [[ -s "$s" ]] || fail "[$ds] TSV column $c file missing: $s"
    copy_verified "$s" "$rf/$(basename "$s")" || fail "[$ds] ref copy failed: $s"
  done
  # inventory: every staged GEX file has the source's byte size (re-checked on every run)
  for f in "${FQS[@]}"; do
    t="$fq/${f#$src/}"; [[ "$(stat -c %s "$t")" == "$(stat -Lc %s "$f")" ]] || fail "[$ds] $t size != $f"
  done
  R1=$(printf '%s\n' "${FQS[@]}" | grep '_R1_' | sed "s#^$src/#$fq/#" | paste -sd,)
  R2=$(printf '%s\n' "${FQS[@]}" | grep '_R2_' | sed "s#^$src/#$fq/#" | paste -sd,)
  [[ $(tr ',' '\n' <<<"$R1" | wc -l) == $(tr ',' '\n' <<<"$R2" | wc -l) ]] || fail "[$ds] R1/R2 count differ"
  exec 9>"$STAGED.lock"; flock 9   # concurrent --dataset runs share datasets_staged.tsv
  awk -F'\t' -v d="$ds" '$1!=d' "$STAGED" > "$STAGED.tmp"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$ds" "$(tsvcol "$ds" 5)" "$rf/$(basename "$(tsvcol "$ds" 6)")" \
    "$rf/$(basename "$(tsvcol "$ds" 7)")" "$rf/$(basename "$(tsvcol "$ds" 8)")" "$R1" "$R2" >> "$STAGED.tmp"
  mv "$STAGED.tmp" "$STAGED"; flock -u 9
  ( cd "$rf" && sha256sum * ) > "$DD/refs_$ds.sha256"
  printf '%s\n' "${FQS[@]}" | xargs stat -Lc '%s %n' > "$DD/fastq_$ds.sizes"
  date -u +%FT%TZ > "$SD/$ds.done"
  say "[$ds] staged + size-verified: ${#FQS[@]} FASTQs ($(du -sh "$fq" | cut -f1)), refs: $(ls "$rf" | paste -sd' ')"
done

ALLDONE=1; for ds in $DATASETS; do [[ -f "$SD/$ds.done" ]] || ALLDONE=0; done
if (( ALLDONE )); then
  LINES=("datasets_staged.tsv	$(sha256sum "$STAGED" | cut -c1-64)")
  for ds in $DATASETS; do LINES+=("refs_$ds.sha256	$(sha256sum "$DD/refs_$ds.sha256" | cut -c1-64)"); done
  write_done "${LINES[@]}"
else
  say "staged: $(cd "$SD" && ls *.done | sed 's/.done//' | paste -sd' ') ; DONE not written until all of: $DATASETS"
fi
