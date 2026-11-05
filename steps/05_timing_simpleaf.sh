#!/usr/bin/env bash
# 05: simpleaf timing on the config.env toolchain.
#   Per dataset: $WARMUP warmup + $REPS timed `simpleaf multiplex-quant` runs, each in a docker container
#   (--cpus $THREADS --cpuset $CPUSET --mem $MEM_LIMIT, lib/run_in_container.sh, peak = host-side
#   proc-tree VmHWM sum). Before EVERY run: idle guard (waits, never kills, while a foreign heavy job such
#   as A30_concordance_k19.sh runs), then every input (step-02 FASTQs + barcode lists, step-03 index) must
#   resolve to the $SCR filesystem (lib/input_mounts.py; table rep<r>.input_mounts.tsv, source in the rep
#   json meta), then the FASTQs of ALL datasets are evicted from the page cache (lib/evict_cache.py) and residency is
#   proven < $EVICT_MAX_RESIDENT_PCT % (lib/page_resident.py). The index is not primed or evicted.
#   Reported: MIN wall and MAX VmHWM (GiB) over the timed reps. Every run (warmup included) must give the
#   same order-independent matrix digest as rep1; a mismatch is written to digest_mismatch.tsv and fails.
#   rep1's digest must also equal $DIGEST_REF_TSV (per-probe reference) and, when step 06 ran first, step 06's digest.
#   Per-run gates: exit 0, VmHWM source, lz4 collate, geometry, crlike_umi_edit == $UMI_EDIT_DIST,
#   execution plan, non-empty matrix (lib/quant_common.sh check_mq_output).
#   rep1's output is KEPT at $SD/<ds>/rep1 (af_quant incl. map.collated.rad + collation_manifest.bin +
#   alevin/ + json) as the reference outputs for step 06; the other runs are deleted after their checks.
# Outputs ($DD = $DUR/05_timing_simpleaf):
#   <ds>/rep<r>.json, .series.csv, .evidence/ (container log + af json), rep<r>.done (digest line)
#   timing_reps.tsv (one row per run) and results.tsv (per dataset: MIN wall, MAX VmHWM, digest)
#   cgroup_enforcement.txt: whether docker enforced --cpuset-cpus / --memory here. On this host it does NOT
#   (cgroup v1 without cpuset/memory controllers: "Cpuset discarded"); only the --cpus quota applies.
#   This is recorded, not fatal.
# Runtime at 32t: ~2 min (Flexv2) to ~24 min (320k) per run, ~3.5 h for all 5 x 4 runs.
# Usage: 05_timing_simpleaf.sh [--dataset DS] [--detach]
#   --detach re-launches itself with setsid nohup (log: $DD/logs/run_<ts>.log, pid: $DD/run.pid); poll the log.
#   Resumable: a run with rep<r>.done is skipped; the warmup is re-run whenever a timed rep is pending.
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"
source "$REPRO/lib/quant_common.sh"

DSL="$DATASETS"; DETACH=0; ARGS=()
while [[ $# -gt 0 ]]; do case "$1" in
  --dataset) DSL="$2"; ARGS+=("$1" "$2"); shift 2;;
  --detach) DETACH=1; shift;;
  *) fail "unknown arg: $1";;
esac; done
if (( DETACH )); then
  L="$LOGD/run_$(date -u +%Y%m%dT%H%M%SZ).log"
  setsid nohup bash "$0" "${ARGS[@]}" > "$L" 2>&1 < /dev/null &
  echo $! > "$DD/run.pid"; say "detached pid $! -> $L"; exit 0
fi

[[ -f "$DUR/01_toolchain/DONE" ]] || fail "step 01 not done ($DUR/01_toolchain/DONE)"
[[ -x "$SA_BIN" && -s "$AF_HOME/simpleaf_info.json" ]] || fail "simpleaf or $AF_HOME/simpleaf_info.json missing (scratch lost? re-run step 01)"
mkdir -p "$SD"
cgroup_probe "$DD/cgroup_enforcement.txt" "$THREADS" "$CPUSET" "$MEM_LIMIT"
REPS_TSV="$DD/timing_reps.tsv"; MISM="$DD/digest_mismatch.tsv"
HDR="dataset\trep\twall_s\tvmhwm_gib\tmem_source\tcpu_seconds\tresident_after_evict_pct\tdecoder_engaged\tmap_threads\tdecode_slots\tnum_mapped\tn_rows\tn_cols\tnnz\ttotal\tdigest"

# container mounts: tools read-only, scratch + ALEVIN_FRY_HOME read-write, index read-only
MOUNTS=( --mount "$REPRO:$REPRO:ro" --mount "$SCR:$SCR" )
[[ "$AF_HOME" == "$SCR"/* ]] || MOUNTS+=( --mount "$AF_HOME:$AF_HOME" )
[[ "$INDEX_ROOT" == "$SCR"/* ]] || MOUNTS+=( --mount "$INDEX_ROOT:$INDEX_ROOT:ro" )
# every staged FASTQ of every dataset is evicted before each run: a run must not start with an earlier
# dataset's FASTQs in the page cache (a full cache made 320k reads wait on reclaim, +17-44% wall at equal CPU-s)
ALL_FQS=(); while IFS=$'\t' read -r _ _ _ _ _ r1 r2; do IFS=',' read -ra _f <<< "$r1,$r2"; ALL_FQS+=("${_f[@]}"); done \
  < <(tail -n +2 "$DUR/02_stage_inputs/datasets_staged.tsv")
(( ${#ALL_FQS[@]} > 0 )) || fail "no FASTQs in $DUR/02_stage_inputs/datasets_staged.tsv"

# run_rep DS REP -> appends one row to $D/reps.tsv and writes $D/rep$REP.done
run_rep(){
  local ds="$1" rep="$2" D="$DD/$1" OUT="$SD/$1/rep$2" M="$SD/$1/metrics/rep$2.json" FQS res ec pk msrc dig
  idle_guard noload
  rm -rf "$OUT"; mkdir -p "$OUT" "$(dirname "$M")" "$D"
  IFS=',' read -ra FQS <<< "$DS_R1,$DS_R2"
  scr_inputs_check "$ds rep$rep" "$D/rep$rep.input_mounts.tsv" "${FQS[@]}" "$DS_CBC" "$DS_SBC" "$DS_IDX" "$DS_T2G"
  evict_check "$ds rep$rep" "${ALL_FQS[@]}"; res="$RES"
  say "[$ds rep$rep] evicted all ${#ALL_FQS[@]} staged FASTQs (resident ${res}%; page cache $(awk '/^Cached:/{printf "%.1f", $2/1048576}' /proc/meminfo) GiB) -> multiplex-quant"
  mq_args "$ds" "$OUT"
  if ! bash "$REPRO/lib/run_in_container.sh" \
      --name "repro05_${ds}_rep${rep}" --image "$TIMING_IMAGE" --cpus "$THREADS" --cpuset "$CPUSET" --mem "$MEM_LIMIT" \
      "${MOUNTS[@]}" --workdir "$OUT" --env "ALEVIN_FRY_HOME=$AF_HOME" \
      --metrics-out "$M" --series-out "${M%.json}.series.csv" \
      --meta dataset="$ds" --meta tool=simpleaf --meta step=05_timing_simpleaf --meta rep="$rep" \
      --meta threads="$THREADS" --meta decoder="$DECODER" --meta kmer="$KMER" --meta excluded_probes="$EXCLUDED_PROBES" \
      --meta cell_bc_correction="$CELL_BC_CORRECTION" --meta sample_bc_correction="$SAMPLE_BC_CORRECTION" \
      --meta umi_edit_dist="$UMI_EDIT_DIST" --meta resident_after_evict_pct="$res" \
      --meta input_source="$INPUT_SOURCE" --meta input_mounts="$D/rep$rep.input_mounts.tsv" \
      -- "${MQ_ARGS[@]}"; then
    tail -30 "${M%.json}.container.log"; fail "[$ds rep$rep] multiplex-quant container failed"
  fi
  ec=$(jget "$M" exit_code); pk=$(jget "$M" peak_rss_gb); msrc=$(jget "$M" mem_source)
  [[ "$ec" == 0 ]] || fail "[$ds rep$rep] exit_code=$ec"
  [[ "$msrc" == "host:proc_tree_vmhwm_sum" ]] || fail "[$ds rep$rep] mem_source=$msrc, not the VmHWM proc-tree sum"
  awk -v p="$pk" 'BEGIN{exit !(p+0 > 0)}' || fail "[$ds rep$rep] null/zero peak '$pk'"
  check_mq_output "$ds" "$OUT" "${M%.json}.container.log"
  dig=$(matrix_digest "$OUT/af_quant/alevin") || fail "[$ds rep$rep] matrix digest failed"
  local nr nc nnz tot dg; IFS=$'\t' read -r nr nc nnz tot dg <<<"$dig"
  (( nnz > 0 )) || fail "[$ds rep$rep] empty matrix"
  # fail fast (before the next rep): every run must equal the per-probe reference and step 06's untimed matrix
  check_ref_digest "$ds" "$dg" "step-05 rep$rep"
  if [[ -f "$DUR/06_gene_matrices/$ds/digest.tsv" && "$(cut -f5 "$DUR/06_gene_matrices/$ds/digest.tsv")" != "$dg" ]]; then
    fail "[$ds rep$rep] digest $dg != step-06 digest $(cut -f5 "$DUR/06_gene_matrices/$ds/digest.tsv")"
  fi
  # evidence -> durable
  local EV="$D/rep$rep.evidence"; rm -rf "$EV"; mkdir -p "$EV"
  cp -f "$M" "${M%.json}.series.csv" "$D/"; cp -f "${M%.json}.container.log" "$EV/container.log"
  cp -f "$OUT"/af_quant/*.json "$OUT/af_map/map_info.json" "$EV/"
  ls -laR "$OUT/af_quant" > "$EV/af_quant.ls.txt"
  [[ -s "$D/reps.tsv" ]] || echo -e "$HDR" > "$D/reps.tsv"
  awk -F'\t' -v r="$rep" '$2!=r' "$D/reps.tsv" > "$D/reps.tsv.tmp"
  echo -e "$ds\t$rep\t$(jget "$M" wall_s)\t$pk\t$msrc\t$(jget "$M" cpu_seconds)\t$res\t$ENGAGED\t$MAP_THREADS\t$DECODE_SLOTS\t$NUM_MAPPED\t$nr\t$nc\t$nnz\t$tot\t$dg" >> "$D/reps.tsv.tmp"
  mv "$D/reps.tsv.tmp" "$D/reps.tsv"
  echo -e "$dig" > "$D/rep$rep.done"
  say "[$ds rep$rep] wall_s=$(jget "$M" wall_s) vmhwm_gib=$pk engaged=$ENGAGED map/decode=$MAP_THREADS/$DECODE_SLOTS digest=$dg"
  [[ "$rep" == 1 ]] || rm -rf "$OUT"     # rep1 = reference output for step 06
}

for ds in $DSL; do
  D="$DD/$ds"
  if [[ -f "$D/.done" ]]; then say "[$ds] done -> skip"; continue; fi
  resolve_ds "$ds"
  say "[$ds] geometry=$DS_GEOM index=$DS_IDX fastq_streams=$DS_NFQ"
  # resume safety: reps of one dataset must share argv + binaries (a config.env change mid-way would
  # otherwise surface only as a digest mismatch against an old-config rep1)
  mq_args "$ds" OUT
  FP=$( { printf '%s\n' "${MQ_ARGS[@]}"; sha256sum "$SA_BIN" "$SRC/alevin-fry/target/release/alevin-fry" "$SRC/piscem/target/release/piscem" "$DS_T2G"; } | sha256sum | cut -c1-16)
  mkdir -p "$D"
  if [[ -f "$D/fingerprint" && "$(cat "$D/fingerprint")" != "$FP" ]]; then
    fail "[$ds] multiplex-quant argv or binaries changed since this dataset's earlier reps ($(cat "$D/fingerprint") -> $FP); move $D away to re-time"
  fi
  echo "$FP" > "$D/fingerprint"
  TIMED=$(seq 1 "$REPS")
  PENDING=""; for r in $TIMED; do [[ -f "$D/rep$r.done" ]] || PENDING+="$r "; done
  if [[ -f "$D/rep1.done" && ! -d "$SD/$ds/rep1/af_quant" ]]; then
    say "[$ds] WARNING: rep1 is done but its output $SD/$ds/rep1 is gone (scratch lost); step 06 needs --untimed for $ds"
  fi
  idle_guard load
  for w in $(seq 1 "$WARMUP"); do run_rep "$ds" "warmup$w"; done
  for r in $PENDING; do run_rep "$ds" "$r"; done
  # digest invariance: every run of this dataset (warmups of this session included) == rep1
  REF=$(cut -f5 "$D/rep1.done")
  BAD=0
  for f in "$D"/rep*.done; do
    if [[ "$(cut -f5 "$f")" != "$REF" ]]; then
      [[ -s "$MISM" ]] || echo -e "dataset\trun\tdigest\trep1_digest" > "$MISM"
      echo -e "$ds\t$(basename "$f" .done)\t$(cut -f5 "$f")\t$REF" >> "$MISM"; BAD=1
    fi
  done
  (( BAD == 0 )) || fail "[$ds] COUNT MISMATCH across reps (see $MISM)"
  check_ref_digest "$ds" "$REF" "step-05 reps"
  if [[ -f "$DUR/06_gene_matrices/$ds/digest.tsv" ]]; then
    [[ "$(cut -f5 "$DUR/06_gene_matrices/$ds/digest.tsv")" == "$REF" ]] || fail "[$ds] step-05 digest $REF != step-06 digest $(cut -f5 "$DUR/06_gene_matrices/$ds/digest.tsv")"
  fi
  date -u +%FT%TZ > "$D/.done"
  say "[$ds] DONE: all $(ls "$D"/rep*.done | wc -l) runs digest $REF"
done

# ---------------- results (rebuilt from every completed dataset) ----------------
{ echo -e "$HDR"; for ds in $DATASETS; do if [[ -s "$DD/$ds/reps.tsv" ]]; then tail -n +2 "$DD/$ds/reps.tsv"; fi; done; } > "$REPS_TSV"
"$PY" - "$REPS_TSV" "$DD/results.tsv" "$DATASETS" "$DD" <<'PY'
import csv, os, sys
reps, out, order, dd = sys.argv[1], sys.argv[2], sys.argv[3].split(), sys.argv[4]
rows = [r for r in csv.DictReader(open(reps), delimiter="\t") if not r["rep"].startswith("warmup")]
with open(out, "w") as f:
    f.write("dataset\tn_reps\twall_min_s\tvmhwm_max_gib\twall_reps_s\tvmhwm_reps_gib\tnum_mapped\tdecoder_engaged\tdigest\n")
    for ds in order:
        if not os.path.exists(os.path.join(dd, ds, ".done")):
            continue
        rs = sorted((r for r in rows if r["dataset"] == ds), key=lambda r: r["rep"])
        assert len({r["digest"] for r in rs}) == 1 and len({r["num_mapped"] for r in rs}) == 1, ds
        f.write(f"{ds}\t{len(rs)}\t{min(float(r['wall_s']) for r in rs):.3f}\t"
                f"{max(float(r['vmhwm_gib']) for r in rs):.3f}\t"
                f"{','.join(r['wall_s'] for r in rs)}\t{','.join(r['vmhwm_gib'] for r in rs)}\t"
                f"{rs[0]['num_mapped']}\t{rs[0]['decoder_engaged']}\t{rs[0]['digest']}\n")
PY
column -t -s$'\t' "$DD/results.tsv"

ALLDONE=1; DONE_DS=""
for ds in $DATASETS; do if [[ -f "$DD/$ds/.done" ]]; then DONE_DS+="$ds "; else ALLDONE=0; fi; done
if (( ALLDONE )); then
  LINES=("results.tsv	$(sha256sum "$DD/results.tsv" | cut -c1-64)" "timing_reps.tsv	$(sha256sum "$REPS_TSV" | cut -c1-64)" \
         "cgroup_enforcement.txt	$(paste -sd' ' "$DD/cgroup_enforcement.txt")")
  for ds in $DATASETS; do LINES+=("$ds/rep1 matrix digest	$(cut -f5 "$DD/$ds/rep1.done")"); done
  write_done "${LINES[@]}"
else
  say "datasets done: ${DONE_DS:-none} ; DONE needs all of: $DATASETS"
fi
