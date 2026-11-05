#!/usr/bin/env bash
# Writes data/perprobe_digest_ref.tsv: the reference matrix digest of every dataset under the per-probe default
# (simpleaf multiplex-quant -> alevin-fry quant-probe on the step-03 probe-set index, no --t2g-map).
# One UNTIMED run per dataset at $ANALYSIS_THREADS threads (steps 01-03 done); steps 04/05/06 then require the
# same digest at their own thread counts, so the reference also checks thread invariance.
# Usage: lib/make_perprobe_ref.sh [DATASET ...]   (default: $DATASETS). Output dirs: $SCR/perprobe_ref/<ds>.
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"; SD="$SCR/perprobe_ref"
source "$REPRO/lib/quant_common.sh"
OUT_TSV="$REPRO/data/perprobe_digest_ref.tsv"; DSL="${*:-$DATASETS}"
TMP="$SD/ref.tsv.part"; mkdir -p "$SD"; printf 'dataset\tdigest\tfirst_seen\n' > "$TMP"
for ds in $DSL; do
  resolve_ds "$ds"; O="$SD/$ds"; rm -rf "$O"; mkdir -p "$O"
  mq_args "$ds" "$O" "$ANALYSIS_THREADS"
  ( cd "$O" && ALEVIN_FRY_HOME="$AF_HOME" "${MQ_ARGS[@]}" ) > "$O.log" 2>&1 || { tail -20 "$O.log"; fail "[$ds] multiplex-quant failed"; }
  grep -qE "quant cmd: [^ ]+ quant-probe " "$O.log" || fail "[$ds] simpleaf did not run quant-probe (per-probe default)"
  check_mq_output "$ds" "$O" "$O.log"
  dg=$(matrix_digest "$O/af_quant/alevin" | cut -f5)
  printf '%s\t%s\t%s\n' "$ds" "$dg" "untimed/t$ANALYSIS_THREADS/simpleaf $(git -C "$SRC/simpleaf" rev-parse --short HEAD)" >> "$TMP"
  say "[$ds] per-probe digest $dg"
done
mv "$TMP" "$OUT_TSV"; column -t "$OUT_TSV"
