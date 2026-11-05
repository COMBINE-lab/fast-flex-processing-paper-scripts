#!/usr/bin/env bash
# 05b: the paper's cyto $CYTO_VERSION and Cell Ranger $CELLRANGER_VERSION runs: wall time + peak VmHWM, and the
#   cyto / Cell Ranger matrices of every later step. It is the only producer of cyto and Cell Ranger outputs:
#   step 14 reads results.tsv (Figure 1), steps 07-13 read the rep1 matrices ($CYTO_OUT, $CR_OUT/<ds>/$CR_SUBPATH).
#
#   Protocol (the same as step 05 for simpleaf):
#   - every run in a docker container (lib/run_in_container.sh) with --cpus $THREADS --cpuset $CPUSET
#     --mem $MEM_LIMIT; cyto runs with -T $THREADS, Cell Ranger with `multi --localcores $THREADS
#     --localmem $CR_LOCALMEM --disable-ui` (create-bam false).
#   - datasets in $TIMING_ORDER (smallest GEX FASTQ first; must be a permutation of $DATASETS); all Cell Ranger
#     runs first, then all cyto runs.
#   - every input on the local NVMe ($SCR): FASTQs and barcode / probe lists from step 02
#     ($DUR/02_stage_inputs/datasets_staged.tsv); cyto's gene-probes TSV and the Cell Ranger reference are
#     staged to $SCR/refs by this step; the Cell Ranger multi CSV is written by this step from
#     data/cr_multi/<ds>.csv (lib/quant_common.sh cr_multi_csv) with every path on $SCR.
#   - before EVERY run (warmup and reps): every input (for Cell Ranger also every file the CSV names) must
#     resolve to the $SCR filesystem (lib/input_mounts.py; table rep<r>.input_mounts.tsv, source in the
#     rep json meta); then the FASTQs are evicted from the page cache and their residency is proven
#     < $EVICT_MAX_RESIDENT_PCT %. The reference / gene-probe files are not primed or evicted.
#   - cyto: $WARMUP warmup + $REPS timed reps. Cell Ranger: no warmup, $CR_REPS rep (~35 min to ~7.6 h each).
#   - reported: MIN wall and MAX VmHWM (GiB) over the timed reps. Cell Ranger's peak comes from
#     lib/cr_peak_sampler.py (Martian workers reparent away from the container init, so the container
#     sampler misses them).
#   - every run goes through the step-05 idle guard (waits, never kills); a run is valid only with exit 0
#     and a non-null VmHWM peak.
#   On this host docker discards --cpuset-cpus and --memory (only the --cpus quota applies); see
#   cgroup_enforcement.txt.
#
# Outputs ($DD = $DUR/05b_timing_cyto_cellranger):
#   <ds>/cyto/rep<r>.{json,series.csv,container.log,input_mounts.tsv}, <ds>/cyto/counts/ (rep1 matrices)
#   <ds>/cellranger/rep<r>.* (+ peak json), <ds>/cellranger/cr_multi_<ds>_gexonly.csv,
#   <ds>/cellranger/per_sample_outs_rep1/ (rep1 per-sample outputs: filtered matrices, metrics)
#   results.tsv, cgroup_enforcement.txt, DONE.
# Usage: 05b_timing_cyto_cellranger.sh [--tool cyto|cellranger|both] [--dataset DS] [--dry-run] [--detach]
#   --dry-run prints every container command; it writes only the log and the Cell Ranger CSVs (under $SD).
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"
source "$REPRO/lib/quant_common.sh"

TOOLS="cellranger cyto"; DSL="$TIMING_ORDER"; DRY=0; DETACH=0; ARGS=()   # Cell Ranger first, each tool smallest to largest
while [[ $# -gt 0 ]]; do case "$1" in
  --tool) case "$2" in cyto|cellranger) TOOLS="$2";; both) ;; *) fail "bad --tool $2";; esac; ARGS+=("$1" "$2"); shift 2;;
  --dataset) DSL="$2"; ARGS+=("$1" "$2"); shift 2;;
  --dry-run) DRY=1; ARGS+=("$1"); shift;;
  --detach) DETACH=1; shift;;
  *) fail "unknown arg: $1";;
esac; done
[[ "$(tr ' ' '\n' <<<"$TIMING_ORDER" | grep . | sort | paste -sd' ')" == "$(tr ' ' '\n' <<<"$DATASETS" | grep . | sort | paste -sd' ')" ]] \
  || fail "TIMING_ORDER '$TIMING_ORDER' is not a permutation of DATASETS '$DATASETS'"
for ds in $DSL; do [[ " $DATASETS " == *" $ds "* ]] || fail "--dataset $ds is not in DATASETS"; done

summarize(){ # OUT "ds tool dir"... : MIN wall / MAX VmHWM over the timed rep<N>.json of each dir -> one TSV row each
  "$PY" - "$@" <<'PY'
import glob, json, os, sys
out = sys.argv[1]
rows = [a.split() for a in sys.argv[2:]]
with open(out, "w") as f:
    f.write("dataset\ttool\tn_reps\twall_min_s\tvmhwm_max_gib\twall_reps_s\tvmhwm_reps_gib\tcpus\tcpuset\tmem_source\tsource_dir\n")
    for ds, tool, d in rows:
        js = sorted(j for j in glob.glob(os.path.join(d, "rep[0-9]*.json")) if os.path.basename(j).count(".") == 1)
        assert js, f"no rep*.json under {d}"
        m = [json.load(open(j)) for j in js]
        bad = [j for j, x in zip(js, m) if x.get("exit_code") != 0 or not x.get("peak_rss_gb")]
        assert not bad, f"invalid reps (exit != 0 or null peak): {bad}"
        w = [x["wall_s"] for x in m]; p = [x["peak_rss_gb"] for x in m]
        f.write(f"{ds}\t{tool}\t{len(m)}\t{min(w):.3f}\t{max(p):.3f}\t{','.join(map(str, w))}\t"
                f"{','.join(map(str, p))}\t{m[0]['cpus']}\t{m[0]['cpuset']}\t{'|'.join(sorted({x['mem_source'] for x in m}))}\t{d}\n")
PY
}

if (( DETACH )); then
  L="$LOGD/run_$(date -u +%Y%m%dT%H%M%SZ).log"
  setsid nohup bash "$0" "${ARGS[@]}" > "$L" 2>&1 < /dev/null &
  echo $! > "$DD/run.pid"; say "detached pid $! -> $L"; exit 0
fi

[[ -x "$CYTO_BIN" && -x "$CELLRANGER_BIN" ]] || fail "cyto ($CYTO_BIN) or cellranger ($CELLRANGER_BIN) missing (run step 01)"
[[ "$("$CYTO_BIN" --version 2>&1 | head -1)" == "cyto $CYTO_VERSION" ]] || fail "cyto is not $CYTO_VERSION"
# cyto's python helpers (cell-filter, pycyto) are uv tools that cyto installs itself (pinned versions);
# the login env sets the arch-suffixed XDG/UV dirs (e.g. ~/.local/share-x86_64-glibc2.39), defaults otherwise
UV_BIN_DIR="$HOME/.local/bin"; XDG_DATA_HOME="${XDG_DATA_HOME:-$HOME/.local/share}"; XDG_CACHE_HOME="${XDG_CACHE_HOME:-$HOME/.cache}"
UV_CACHE_DIR="${UV_CACHE_DIR:-$HOME/.cache/uv}"; UV_TOOL_DIR="${UV_TOOL_DIR:-$XDG_DATA_HOME/uv/tools}"
CR_ROOT="$(dirname "$(dirname "$CELLRANGER_BIN")")"
mkdir -p "$SD"
(( DRY )) || cgroup_probe "$DD/cgroup_enforcement.txt" "$THREADS" "$CPUSET" "$MEM_LIMIT"

run(){ # echo the command in dry-run mode, else run it
  if (( DRY )); then printf '[dry-run]'; printf ' %q' "$@"; echo; else "$@"; fi
}
pre_run(){ # LABEL MOUNTS_TSV FASTQS(comma list) INPUT... -> input mount check + FASTQ eviction; sets RES INPUT_SOURCE
  local label="$1" mt="$2" fqs; IFS=',' read -ra fqs <<<"$3"; shift 3
  if (( DRY )); then RES="dry"; INPUT_SOURCE="dry"; echo "[dry-run] input mount check of ${#fqs[@]} FASTQs + $# inputs; evict + residency check"; return 0; fi
  scr_inputs_check "$label" "$mt" "${fqs[@]}" "$@"
  evict_check "$label" "${fqs[@]}"
}
valid_rep(){ # metrics.json label
  local ec pk; ec=$(jget "$1" exit_code); pk=$(jget "$1" peak_rss_gb)
  [[ "$ec" == 0 ]] || fail "[$2] exit_code=$ec ($1)"
  awk -v p="$pk" 'BEGIN{exit !(p+0 > 0)}' || fail "[$2] null/zero peak '$pk' ($1)"
}

# ---------------- cyto ----------------
cyto_geo(){ case "$1" in
  v1*) echo "--preset gex-v1";;
  v2_ConfigA) echo "--geometry [barcode][umi:12][:13][probe]|[gex] --remap-window 5";;
  *) fail "no cyto geometry for chemistry $1";;
esac; }

run_cyto(){
  local ds="$1" D="$DD/$1/cyto" gp inputs=() r1 r2 i rep geo
  if [[ -f "$D/.done" ]]; then say "[$ds cyto] done -> skip"; return 0; fi
  resolve_ds "$ds"
  [[ -s "$REFS_SRC/$ds/cyto_gene_probes.tsv" ]] || fail "[$ds] cyto gene-probes $REFS_SRC/$ds/cyto_gene_probes.tsv missing"
  gp="$SCR/refs/$ds/cyto_gene_probes.tsv"
  (( DRY )) || stage_file "$REFS_SRC/$ds/cyto_gene_probes.tsv" "$gp"
  IFS=',' read -ra r1 <<<"$DS_R1"; IFS=',' read -ra r2 <<<"$DS_R2"
  (( ${#r1[@]} == ${#r2[@]} )) || fail "[$ds] R1/R2 counts differ"
  for i in "${!r1[@]}"; do inputs+=( "${r1[$i]}" "${r2[$i]}" ); done     # cyto takes R1 R2 R1 R2 ...
  read -ra geo <<<"$(cyto_geo "$DS_CHEM")"
  local tmpd="$SCR/tmp/05b_cyto/$ds"; mkdir -p "$D" "$tmpd"
  idle_guard load
  for rep in $(seq -f 'warmup%g' 1 "$WARMUP") $(seq 1 "$REPS"); do
    if [[ -f "$D/rep$rep.done" ]]; then continue; fi
    idle_guard noload
    local OUT="$SD/$ds/cyto/rep$rep" M="$SD/$ds/cyto/metrics/rep$rep.json"
    rm -rf "$OUT" "$tmpd"/*; mkdir -p "$(dirname "$M")"
    pre_run "$ds cyto rep$rep" "${M%.json}.input_mounts.tsv" "$DS_R1,$DS_R2" "$gp" "$DS_CBC" "$DS_SBC"
    say "[$ds cyto rep$rep] resident=${RES} -> cyto workflow gex"
    run bash "$REPRO/lib/run_in_container.sh" --name "repro05b_cyto_${ds}_${rep}" --image "$TIMING_IMAGE" \
      --cpus "$THREADS" --cpuset "$CPUSET" --mem "$MEM_LIMIT" \
      --mount "$(dirname "$CYTO_BIN"):$(dirname "$CYTO_BIN"):ro" --mount "$HOME/.local:$HOME/.local" --mount "$XDG_CACHE_HOME:$XDG_CACHE_HOME" \
      --mount "$SCR:$SCR" \
      --env "HOME=$HOME" --env "PATH=$UV_BIN_DIR:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
      --env "UV_CACHE_DIR=$UV_CACHE_DIR" --env "XDG_DATA_HOME=$XDG_DATA_HOME" --env "XDG_CACHE_HOME=$XDG_CACHE_HOME" \
      --env "UV_TOOL_DIR=$UV_TOOL_DIR" --env "TMPDIR=$tmpd" \
      --metrics-out "$M" --series-out "${M%.json}.series.csv" \
      --meta dataset="$ds" --meta tool=cyto --meta step=05b --meta rep="$rep" --meta resident_after_evict_pct="$RES" \
      --meta input_source="$INPUT_SOURCE" --meta input_mounts="$D/rep$rep.input_mounts.tsv" \
      -- "$CYTO_BIN" workflow gex -c "$gp" -w "$DS_CBC" -p "$DS_SBC" -o "$OUT" -T "$THREADS" --format mtx "${geo[@]}" "${inputs[@]}" \
      || { (( DRY )) || tail -20 "${M%.json}.container.log"; fail "[$ds cyto rep$rep] container failed"; }
    (( DRY )) && continue
    valid_rep "$M" "$ds cyto rep$rep"
    [[ -n "$(find "$OUT/counts" -name 'matrix.mtx*' -size +0 -print -quit)" ]] || fail "[$ds cyto rep$rep] no matrix under $OUT/counts"
    cp -f "$M" "${M%.json}.series.csv" "${M%.json}.container.log" "${M%.json}.input_mounts.tsv" "$D/"
    if [[ "$rep" == 1 ]]; then   # rep1 = the cyto matrices of steps 07-13 ($CYTO_OUT/<ds>/cyto/counts)
      rm -rf "$D/counts" "$D/counts.part"; cp -r "$OUT/counts" "$D/counts.part"; mv "$D/counts.part" "$D/counts"
    fi
    say "[$ds cyto rep$rep] wall_s=$(jget "$M" wall_s) vmhwm_gib=$(jget "$M" peak_rss_gb)"
    rm -rf "$OUT" "$tmpd"/*
    date -u +%FT%TZ > "$D/rep$rep.done"
  done
  (( DRY )) && return 0
  [[ -d "$D/counts" ]] || fail "[$ds cyto] rep1 matrices $D/counts missing"
  date -u +%FT%TZ > "$D/.done"
}

# ---------------- Cell Ranger ----------------
run_cr(){
  local ds="$1" D="$DD/$1/cellranger" rep stage=stage
  if [[ -f "$D/.done" ]]; then say "[$ds cellranger] done -> skip"; return 0; fi
  resolve_ds "$ds"; (( DRY )) && stage=nostage
  CSV="$SD/$ds/cellranger/cr_multi_${ds}_gexonly.csv"; cr_multi_csv "$ds" "$CSV" false "$stage"; mkdir -p "$D"
  local -a CSVP; mapfile -t CSVP <<<"$CSV_PATHS"
  idle_guard load
  for rep in $(seq 1 "$CR_REPS"); do
    if [[ -f "$D/rep$rep.done" ]]; then continue; fi
    idle_guard noload
    local id="repro05b_cr_${ds}_rep${rep}" OUTP="$SD/$ds/cellranger/out" M="$SD/$ds/cellranger/metrics/rep$rep.json"
    local cname="repro05b_cr_${ds}_${rep}"
    rm -rf "${OUTP:?}/$id"; mkdir -p "$OUTP" "$(dirname "$M")"
    pre_run "$ds cellranger rep$rep" "${M%.json}.input_mounts.tsv" "$DS_R1,$DS_R2" "$CSV" "${CSVP[@]}"
    run "$PY" "$REPRO/lib/cr_peak_sampler.py" --pipestance "$OUTP/$id" --run-id "$id" --container "$cname" \
      --wait-container-s 3600 --interval 1 --out "${M%.json}.peak.json" --series "${M%.json}.peak.series.csv" \
      --meta dataset="$ds" --meta tool=cellranger --meta rep="$rep" &
    local SAMP=$!
    run bash "$REPRO/lib/run_in_container.sh" --name "$cname" --image "$TIMING_IMAGE" \
      --cpus "$THREADS" --cpuset "$CPUSET" --mem "$MEM_LIMIT" \
      --mount "$CR_ROOT:$CR_ROOT:ro" --mount "$SCR:$SCR" --workdir "$OUTP" \
      --metrics-out "$M" --series-out "${M%.json}.series.csv" \
      --meta dataset="$ds" --meta tool=cellranger --meta step=05b --meta rep="$rep" --meta resident_after_evict_pct="$RES" \
      --meta input_source="$INPUT_SOURCE" --meta input_mounts="$D/rep$rep.input_mounts.tsv" \
      -- "$CELLRANGER_BIN" multi --id "$id" --csv "$CSV" --localcores "$THREADS" --localmem "$CR_LOCALMEM" --disable-ui \
      || { (( DRY )) || tail -20 "${M%.json}.container.log"; fail "[$ds cellranger rep$rep] container failed"; }
    wait "$SAMP" || fail "[$ds cellranger rep$rep] peak sampler failed"
    (( DRY )) && continue
    # the Martian-proof sampler's VmHWM replaces the container sampler's
    "$PY" - "$M" "${M%.json}.peak.json" <<'PY'
import json, sys
m, p = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
m["peak_rss_proctree_init_gb"] = m.get("peak_rss_gb")
for k in ("peak_rss_bytes", "peak_rss_gb", "peak_vmrss_gb", "mem_source"):
    m[k] = p[k]
json.dump(m, open(sys.argv[1], "w"), indent=2)
PY
    valid_rep "$M" "$ds cellranger rep$rep"
    cp -f "$M" "${M%.json}".* "$CSV" "$D/"
    # rep r's per-sample outputs (filtered matrices, barcodes, metrics) -> durable; rep1 = $CR_OUT/<ds>/$CR_SUBPATH
    rm -rf "$D/per_sample_outs_rep$rep" "$D/per_sample_outs_rep$rep.part"
    cp -r "$OUTP/$id/outs/per_sample_outs" "$D/per_sample_outs_rep$rep.part"; mv "$D/per_sample_outs_rep$rep.part" "$D/per_sample_outs_rep$rep"
    say "[$ds cellranger rep$rep] wall_s=$(jget "$M" wall_s) vmhwm_gib=$(jget "$M" peak_rss_gb) ($(jget "$M" mem_source))"
    rm -rf "${OUTP:?}/$id"
    date -u +%FT%TZ > "$D/rep$rep.done"
  done
  (( DRY )) || date -u +%FT%TZ > "$D/.done"
}

for t in $TOOLS; do for ds in $DSL; do
  case "$t" in cyto) run_cyto "$ds";; cellranger) run_cr "$ds";; esac
done; done
(( DRY )) && { say "dry run complete"; exit 0; }

ROWS=(); for ds in $DATASETS; do for t in cyto cellranger; do
  if [[ -f "$DD/$ds/$t/.done" ]]; then ROWS+=("$ds $t $DD/$ds/$t"); fi
done; done
summarize "$DD/results.tsv" "${ROWS[@]}"
column -t -s$'\t' "$DD/results.tsv" | cut -c1-160
ALLDONE=1; for ds in $DATASETS; do for t in cyto cellranger; do [[ -f "$DD/$ds/$t/.done" ]] || ALLDONE=0; done; done
if (( ALLDONE )); then write_done "results.tsv	$(sha256sum "$DD/results.tsv" | cut -c1-64)" \
  "cgroup_enforcement.txt	$(paste -sd' ' "$DD/cgroup_enforcement.txt")"
else say "DONE not written until cyto + cellranger are done for all of: $DATASETS"; fi
