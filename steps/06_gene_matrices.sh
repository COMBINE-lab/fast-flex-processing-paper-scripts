#!/usr/bin/env bash
# 06: durable gene-level matrices (the paper's simpleaf matrices) + the collated RAD step 07 re-quantifies.
#   Source (default): step 05's kept rep1 output, $SCR/05_timing_simpleaf/<ds>/rep1; its matrix digest must
#   equal $DUR/05_timing_simpleaf/<ds>/rep1.done, i.e. it is the timed rep1, untouched.
#   --untimed (the PRIMARY path in run_all.sh: 06 runs before the timing): ONE untimed multiplex-quant with
#   step 05's argv (lib/quant_common.sh mq_args) except --threads $ANALYSIS_THREADS, into $SD/<ds>; if a
#   step-05 rep1 digest exists, the untimed matrix must equal it.
#   Either way the matrix digest must equal the dataset's per-probe reference digest in $DIGEST_REF_TSV (fails loudly).
#   The whole af_quant directory (alevin/ matrix, every json, map.collated.rad + .chunkidx,
#   collation_manifest.bin, permit/sample maps) plus af_map/map_info.json is copied to
#   $DD/<ds>/af_quant, and af_map/{unmapped_bc_count.bin,map_info.json} to $DD/<ds>/af_map (step 13 reads the
#   mapping counters from there instead of re-running the quant); the copy is verified file-by-file (byte sizes, sha256 of RAD + manifest) and by
#   re-digesting the copied matrix. 320k_GEMX's RAD is ~70 GiB, so its copy to durable storage can take a while.
# Outputs: $DD/<ds>/af_quant/, $DD/<ds>/af_map/, $DD/<ds>/sha256.txt, $DD/digest_ref.tsv (dataset, digest, source), DONE.
# Usage: 06_gene_matrices.sh [--dataset DS] [--untimed]      Resumable per dataset ($DD/<ds>/.done).
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"
source "$REPRO/lib/quant_common.sh"

DSL="$DATASETS"; UNTIMED=0
while [[ $# -gt 0 ]]; do case "$1" in
  --dataset) DSL="$2"; shift 2;;
  --untimed) UNTIMED=1; shift;;
  *) fail "unknown arg: $1";;
esac; done
T05="$DUR/05_timing_simpleaf"; DIGREF="$DD/digest_ref.tsv"
mkdir -p "$SD"

for ds in $DSL; do
  D="$DD/$ds"
  if [[ -f "$D/.done" ]]; then say "[$ds] done -> skip"; continue; fi
  REF05=""; if [[ -f "$T05/$ds/rep1.done" ]]; then REF05=$(cut -f5 "$T05/$ds/rep1.done"); fi
  if (( UNTIMED )); then
    resolve_ds "$ds"
    SRCD="$SD/$ds"; SOURCE="untimed"
    if [[ ! -f "$SRCD/.mq_done" ]]; then
      [[ -x "$SA_BIN" && -s "$AF_HOME/simpleaf_info.json" ]] || fail "simpleaf or $AF_HOME missing (re-run step 01)"
      rm -rf "$SRCD"; mkdir -p "$SRCD"
      mq_args "$ds" "$SRCD" "$ANALYSIS_THREADS"
      say "[$ds] untimed multiplex-quant (argv = step 05's except --threads $ANALYSIS_THREADS; digests are thread-invariant): ${MQ_ARGS[*]}"
      ALEVIN_FRY_HOME="$AF_HOME" "${MQ_ARGS[@]}" > "$LOGD/${ds}_multiplex_quant.log" 2>&1 \
        || { tail -30 "$LOGD/${ds}_multiplex_quant.log"; fail "[$ds] multiplex-quant failed"; }
      check_mq_output "$ds" "$SRCD" "$LOGD/${ds}_multiplex_quant.log"
      touch "$SRCD/.mq_done"
    fi
  else
    SRCD="$SCR/05_timing_simpleaf/$ds/rep1"; SOURCE="05_timing_simpleaf/rep1"
    [[ -n "$REF05" ]] || fail "[$ds] no step-05 rep1 ($T05/$ds/rep1.done); run step 05 or use --untimed"
    [[ -d "$SRCD/af_quant" ]] || fail "[$ds] step-05 rep1 output $SRCD is gone (scratch lost?); use --untimed"
  fi
  Q="$SRCD/af_quant"
  DIG=$(matrix_digest "$Q/alevin") || fail "[$ds] digest of $Q/alevin failed"
  DG=$(cut -f5 <<<"$DIG")
  if [[ -n "$REF05" && "$DG" != "$REF05" ]]; then fail "[$ds] source matrix digest $DG != step-05 rep1 $REF05 ($SOURCE)"; fi
  check_ref_digest "$ds" "$DG" "$SOURCE"
  say "[$ds] source $SOURCE digest $DG ($(du -sh "$Q" | cut -f1)) -> $D/af_quant"

  rm -rf "$D/af_quant.part" "$D/af_map.part"; mkdir -p "$D" "$D/af_map.part"
  cp -a "$Q" "$D/af_quant.part"
  cp -f "$SRCD/af_map/map_info.json" "$D/af_quant.part/map_info.json"
  cp -f "$SRCD/af_map/unmapped_bc_count.bin" "$SRCD/af_map/map_info.json" "$D/af_map.part/"
  # verify: same file set with the same byte sizes, sha256 of the RAD + manifest, re-digest of the copy
  diff <(cd "$Q" && find . -type f -printf '%P %s\n' | sort) \
       <(cd "$D/af_quant.part" && find . -type f ! -name map_info.json -printf '%P %s\n' | sort) \
    || fail "[$ds] copied af_quant differs from $Q in file set or sizes"
  : > "$D/sha256.part"
  for f in map.collated.rad collation_manifest.bin; do
    a=$(sha256sum "$Q/$f" | cut -c1-64); b=$(sha256sum "$D/af_quant.part/$f" | cut -c1-64)
    [[ "$a" == "$b" ]] || fail "[$ds] sha256 of copied $f differs"
    echo -e "$f\t$a" >> "$D/sha256.part"
  done
  a=$(sha256sum "$SRCD/af_map/unmapped_bc_count.bin" | cut -c1-64); b=$(sha256sum "$D/af_map.part/unmapped_bc_count.bin" | cut -c1-64)
  [[ "$a" == "$b" ]] || fail "[$ds] sha256 of copied unmapped_bc_count.bin differs"
  echo -e "af_map/unmapped_bc_count.bin\t$a" >> "$D/sha256.part"
  [[ "$(matrix_digest "$D/af_quant.part/alevin")" == "$DIG" ]] || fail "[$ds] copied matrix digest differs"
  rm -rf "$D/af_quant" "$D/af_map"; mv "$D/af_quant.part" "$D/af_quant"; mv "$D/af_map.part" "$D/af_map"; mv "$D/sha256.part" "$D/sha256.txt"
  echo -e "$DIG" > "$D/digest.tsv"
  [[ -s "$DIGREF" ]] || echo -e "dataset\tn_rows\tn_cols\tnnz\ttotal\tdigest\tsource" > "$DIGREF"
  awk -F'\t' -v d="$ds" '$1!=d' "$DIGREF" > "$DIGREF.tmp"
  echo -e "$ds\t$DIG\t$SOURCE" >> "$DIGREF.tmp"; mv "$DIGREF.tmp" "$DIGREF"
  date -u +%FT%TZ > "$D/.done"
  say "[$ds] DONE: $D/af_quant (digest $DG, RAD sha256 $(head -1 "$D/sha256.txt" | cut -f2 | cut -c1-12)...)"
done

ALLDONE=1; for ds in $DATASETS; do [[ -f "$DD/$ds/.done" ]] || ALLDONE=0; done
if (( ALLDONE )); then
  LINES=("digest_ref.tsv	$(sha256sum "$DIGREF" | cut -c1-64)")
  for ds in $DATASETS; do
    LINES+=("$ds/af_quant/alevin matrix digest	$(cut -f5 "$DD/$ds/digest.tsv")")
    while IFS=$'\t' read -r f h; do
      if [[ "$f" == af_map/* ]]; then LINES+=("$ds/$f	$h"); else LINES+=("$ds/af_quant/$f	$h"); fi
    done < "$DD/$ds/sha256.txt"
  done
  write_done "${LINES[@]}"
else
  say "DONE not written until all of: $DATASETS"
fi
