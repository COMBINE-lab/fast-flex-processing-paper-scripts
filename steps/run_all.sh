#!/usr/bin/env bash
# run_all: the end-to-end driver, 00 -> 15 in README order, as a dependency graph.
#   00 -> 01 -> {02 per dataset (parallel), 03} -> 04 (needs 02 $GATE_DATASET + 03)
#   -> 05 timing ALONE, then 05b timing ALONE (cyto + Cell Ranger; each starts only on an idle driver and
#      nothing else starts while it runs; no analysis runs before them)
#   -> per dataset: 06 --untimed                  (06: $RUNALL_GENE_JOBS at once; 04/07/07ed0/07ed1: $RUNALL_QUANT_JOBS)
#   -> per dataset 07 -> {07ed0, 07ed1}; 08, 09, 13bam concurrently; 10 and 11 after 09; 12 after 07 + 11;
#      13a-d after 13bam; 13e-f after 08 + 13a-d          (07, 08, 09, 13 read 05b's cyto / Cell Ranger outputs)
#   -> 14 -> 15 (manuscript numerals + figures)
# 06 runs --untimed (06 still checks a step-05 rep1 digest whenever one exists).
# Resumable: a node whose completion marker exists is skipped; the steps themselves resume internally.
# A failed node stops its dependents; independent nodes continue; the driver exits non-zero at the end.
# Pre-flight (also shown by --dry-run): refuses to start if any 05-15 output dir / node marker under $DUR or
#   $SCR was not written by run_all (ledger $DUR/run_all/ledger.tsv; markers get a run_all_writer line).
# Logs: $DUR/run_all/run_all.log and $DUR/run_all/logs/<node>.log.
# --redo NODE[,NODE]: invalidate the named nodes and every node that depends on them (data dependencies;
#   the timing nodes' waits are order-only, so a redo of 06 never redoes 05 / 05b; --redo 05b redoes every
#   node that reads cyto / Cell Ranger outputs, then 14 and 15). Their
#   output dirs and markers are MOVED (never deleted) to $DUR/run_all/superseded/<ts>/ ($SCR ones to
#   $SCR/run_all_superseded/<ts>/, same filesystem), each move logged and ledgered. A dir shared with a node
#   outside the redo set is kept and only the marker is moved (the step resumes on its own sub-step markers).
#   Then the run proceeds as usual (pre-flight included).
# --hold NODE[,NODE]: never start these nodes (nor anything downstream of them); the run ends with exit 0
#   once only held work is left ("held: ..."), without writing run_all's DONE. E.g. --hold 15 stops before
#   the manuscript is written.
# --full: produce every optional analysis as well (config.env FULL block); default: only what the paper shows.
#   A node 09-12 marker written without --full is stale under --full, so a later --full run adds the extras.
# Usage: run_all.sh [--full] [--dry-run] [--redo NODE[,NODE]] [--hold NODE[,NODE]]     (launch detached: setsid nohup bash steps/run_all.sh > /dev/null 2>&1 &)
set -euo pipefail
for _a in "$@"; do [[ "$_a" == --full ]] && export FULL=1; done; unset _a   # before config.env: its FULL block
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"
source "$REPRO/lib/marker_rules.sh"
ARGS_ALL="$*"; DRY=0; REDO=""; HOLD=""
while [[ $# -gt 0 ]]; do case "$1" in
  --dry-run) DRY=1; shift;;
  --full) shift;;   # FULL=1 exported above
  --redo) REDO="${2//,/ }"; shift 2;;
  --hold) HOLD="${2//,/ }"; shift 2;;
  *) fail "unknown arg: $1";;
esac; done
S="$REPRO/steps"

# ---------------- the graph: name | deps | class | marker | command ----------------
declare -A DEPS CLASS MARK CMD VALID; ORDER=()   # VALID[n]: optional check that an existing marker is current
node(){ ORDER+=("$1"); DEPS[$1]="$2"; CLASS[$1]="$3"; MARK[$1]="$4"; CMD[$1]="$5"; }
all_of(){ local p="$1" ds o=""; for ds in ${2:-$DATASETS}; do o+="${p}_$ds "; done; echo "$o"; }

# 00, 01 and 03 keep their working state on $SCR (mount, $AF_HOME set-paths, indices), so their run_all
# markers are on $SCR too: after an instance stop they rerun (remount, re-register, restore the index from $DUR)
node 00 ""  light "$SCR/00_scratch_and_env/.done" "bash $S/00_scratch_and_env.sh"
node 01 "00" light "$AF_HOME/simpleaf_info.json" "bash $S/01_toolchain.sh"
for ds in $DATASETS; do
  node "02_$ds" "01" io "$SCR/02_stage_inputs/$ds.done" "bash $S/02_stage_inputs.sh --dataset $ds"
done
for p in $PANELS; do node "03_$p" "01" light "$SCR/index/$p/.done" "bash $S/03_index.sh --panel $p"; done
node 04 "02_$GATE_DATASET $(all_of 03 "$PANELS")" quant "$DUR/04_validate_toolchain/DONE" "bash $S/04_validate_toolchain.sh"
# 05 and 05b (timing) each run ALONE on an idle host: they wait (order-only, ODEPS: not followed by --redo)
# for every node that does not need their output, and nothing else starts while they run.
PRE=""; for n in "${ORDER[@]}"; do PRE+="$n "; done
declare -A ODEPS
node 05 "00 01 $(all_of 02) $(all_of 03 "$PANELS")" timing "$DUR/05_timing_simpleaf/DONE" "bash $S/05_timing_simpleaf.sh"
ODEPS[05]="$PRE"
# 05b is the only producer of the cyto and Cell Ranger outputs ($CYTO_OUT, $CR_OUT): every node that reads
# them depends on it, so --redo 05b re-runs them (and 14, 15)
node 05b "00 01 $(all_of 02)" timing "$DUR/05b_timing_cyto_cellranger/DONE" "bash $S/05b_timing_cyto_cellranger.sh"
ODEPS[05b]="$PRE 05"
# every analysis node starts only after both timing nodes: 06 waits (order-only) for 05 and 05b, and 07-13
# follow from 06 (13bam from 05b), so no analysis ever runs before or beside the benchmark
for ds in $DATASETS; do
  node "06_$ds" "04 02_$ds $(all_of 03 "$PANELS")" gene "$DUR/06_gene_matrices/$ds/.done" "bash $S/06_gene_matrices.sh --dataset $ds --untimed"
  ODEPS["06_$ds"]="05 05b"
done
for ds in $DATASETS; do
  node "07_$ds" "06_$ds 05b" quant "$DUR/07_probe_matrices/$ds/.export_done" "bash $S/07_probe_matrices.sh --dataset $ds"
  VALID["07_$ds"]="export_marker_current $DUR/07_probe_matrices/$ds/.export_done $ds"   # a skipped export is stale once $ds is in PROBE_EXPORT_DATASETS
  node "07ed0_$ds" "07_$ds" quant "$DUR/07_probe_matrices/$ds/.ed0_done" "bash $S/07_probe_matrices.sh --dataset $ds --step ed0"
  # its own node (marker .ed1gene_done): added after 07ed0 had run, so only this arm reruns on a done tree
  node "07ed1_$ds" "07_$ds" quant "$DUR/07_probe_matrices/$ds/.ed1gene_done" "bash $S/07_probe_matrices.sh --dataset $ds --step ed1gene"
done
node 08 "$(all_of 06) 05b" analysis "$DUR/08_concordance_counts/DONE" "bash $S/08_concordance_counts.sh"
node 09 "$(all_of 06) 05b" analysis "$DUR/09_clustering/DONE" "bash $S/09_clustering.sh --step all"
node 10 "09" analysis "$DUR/10_pooled_umap/DONE" "bash $S/10_pooled_umap.sh --step all"
node 11 "09" analysis "$DUR/11_celltypes/DONE" "bash $S/11_celltypes.sh --step all"
node 12 "$(all_of 07) 11" analysis "$DUR/12_probe_biology/DONE" "bash $S/12_probe_biology.sh"
# step 12's sub-step markers carry their settings (doublets, panel-a choice, search datasets): its read-only
# --check says whether DONE is still current; if not, run_all reruns 12 and the step redoes only what changed
step12_current(){
  grep -q -- '"--check"' "$S/12_probe_biology.sh" || fail "$S/12_probe_biology.sh has no --check mode (needed by run_all)"
  bash "$S/12_probe_biology.sh" --check > "$DD/check_12.txt" 2>&1
}
VALID[12]="step12_current"
# --full: a marker of 09-12 must record full 1 (written by write_done), or the node reruns for the extras
full_marker(){ (( ${FULL:-0} == 0 )) || grep -qx "full	1" "$1"; }
for _n in 09 10 11; do VALID[$_n]="full_marker ${MARK[$_n]}"; done; unset _n
step12_full_current(){ step12_current && full_marker "${MARK[12]}"; }
VALID[12]="step12_full_current"
# 13bam: untimed Cell Ranger run WITH BAMs (whole host) on the step-05b inputs; after 05b so it never overlaps a timing run
node 13bam "02_$READLEVEL_DS 05b" analysis "$(dirname "$CR_BAM_OUT")/.done" "bash $S/13_readlevel.sh --step bam"
node 13a-d "06_$READLEVEL_DS 13bam" analysis "$SCR/13_readlevel/d.done" \
  "bash $S/13_readlevel.sh --step a && bash $S/13_readlevel.sh --step b && bash $S/13_readlevel.sh --step c && bash $S/13_readlevel.sh --step d"
node 13e-f "13a-d 08 05b" analysis "$DUR/13_readlevel/DONE" \
  "bash $S/13_readlevel.sh --step e && bash $S/13_readlevel.sh --step f && bash $S/13_readlevel.sh --step numbers"
node 14 "05 05b 08 09 10 11 12 13e-f $(all_of 07ed0) $(all_of 07ed1)" light "$DUR/14_figures/DONE" "bash $S/14_figures_numbers.sh"

alldeps(){ echo "${DEPS[$1]} ${ODEPS[$1]:-}"; }
for n in "${ORDER[@]}"; do for d in $(alldeps "$n"); do [[ -n "${CMD[$d]+x}" ]] || fail "node $n depends on unknown $d"; done; done
for n in $HOLD; do [[ -n "${CMD[$n]+x}" ]] || fail "--hold: unknown node $n"; done

# class caps: quant = full-host multiplex-quant / af quant jobs; timing must run alone
# caps are re-read from config.env (RUNALL_* lines only) on every scheduling loop, so a live run can be retuned
declare -A CAP
load_caps(){
  eval "$(grep -E '^RUNALL_(IO|GENE|QUANT|ANALYSIS)_JOBS=[0-9]+' "$REPRO/config.env" | sed 's/[[:space:]]*#.*//')"
  local v; for v in RUNALL_IO_JOBS RUNALL_GENE_JOBS RUNALL_QUANT_JOBS RUNALL_ANALYSIS_JOBS; do
    [[ "${!v:-}" =~ ^[1-9][0-9]*$ ]] || fail "config.env $v='${!v:-}' is not a positive integer"; done
  local old="${CAP[gene]:-}/${CAP[quant]:-}/${CAP[analysis]:-}/${CAP[io]:-}"
  CAP=([io]=$RUNALL_IO_JOBS [gene]=$RUNALL_GENE_JOBS [quant]=$RUNALL_QUANT_JOBS [analysis]=$RUNALL_ANALYSIS_JOBS [light]=99 [timing]=1)
  local now="${CAP[gene]}/${CAP[quant]}/${CAP[analysis]}/${CAP[io]}"
  [[ -z "${old//\//}" || "$old" == "$now" ]] || say "caps changed (gene/quant/analysis/io): $old -> $now"
}
load_caps
declare -A STATE=() PIDN=(); declare -A RUNNING_BY_CLASS=([io]=0 [gene]=0 [quant]=0 [analysis]=0 [light]=0 [timing]=0)
STALE=""; mkdir -p "$DD"
for n in "${ORDER[@]}"; do
  if [[ -f "${MARK[$n]}" ]] && { [[ -z "${VALID[$n]:-}" ]] || ${VALID[$n]}; }; then STATE[$n]=done
  else STATE[$n]=pending; [[ -f "${MARK[$n]}" ]] && STALE+="$n "; fi
done

# ---------------- pre-flight: nothing under 05-15 may come from outside run_all ----------------
# Every node run_all starts gets its output dirs recorded in $LEDGER; every marker it sees written gets a
# "run_all_writer" line appended (only files named DONE / .done / *.done / .*_done; a TSV or JSON used as a
# marker is left untouched) and its sha256 recorded.
# A 05-15 dir under $DUR or $SCR that is not in the ledger, or a 05-15 node marker whose sha256 does not
# match its ledger entry, is FOREIGN (a test or a hand-run) and run_all refuses to start.
LEDGER="$DD/ledger.tsv"; [[ -s "$LEDGER" ]] || printf 'kind\tpath\tsha256\tnode\trun_id\tdate\n' > "$LEDGER"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)_$$"
protected(){ [[ "$1" =~ ^(0[5-9]|1[0-5]) ]]; }
node_dirs(){ local n="$1" sc m top
  for sc in $(grep -oE "$S/[0-9][0-9a-z_]*\.sh" <<<"${CMD[$n]}" | sort -u); do sc=$(basename "$sc" .sh); echo "$DUR/$sc"; echo "$SCR/$sc"; done
  m="${MARK[$n]}"; for top in "$DUR" "$SCR"; do [[ "$m" == "$top/"* ]] && { m="${m#$top/}"; echo "$top/${m%%/*}"; }; done; }
in_ledger(){ awk -F'\t' -v k="$1" -v p="$2" -v h="${3:-}" '$1==k && $2==p && (h=="" || $3==h){f=1} END{exit !f}' "$LEDGER"; }
preflight(){
  local foreign=() d n h
  say "pre-flight: markers under $DUR:"
  find "$DUR" -maxdepth 4 \( -name DONE -o -name .done -o -name '*.done' -o -name '.*_done' \) 2>/dev/null | sort | sed "s#^$DUR/#    #"
  for d in "$DUR"/* "$SCR"/*; do
    [[ -d "$d" ]] && protected "$(basename "$d")" || continue
    in_ledger dir "$d" || foreign+=("dir $d")
  done
  for n in "${ORDER[@]}"; do
    protected "$n" && [[ -f "${MARK[$n]}" ]] || continue
    h=$(sha256sum "${MARK[$n]}" | cut -c1-64); in_ledger marker "${MARK[$n]}" "$h" || foreign+=("marker ${MARK[$n]} (node $n)")
  done
  if (( ${#foreign[@]} )); then
    say "pre-flight: FOREIGN 05-15 outputs (not written by run_all):"; printf '    %s\n' "${foreign[@]}"; return 1
  fi
  say "pre-flight OK: no foreign 05-15 dir or marker"
}

# ---------------- --redo: closure over data dependents, then move (never delete) ----------------
redo_closure(){ local -A in=(); local n d changed=1
  for n in $REDO; do [[ -n "${CMD[$n]+x}" ]] || fail "--redo: unknown node $n"; in[$n]=1; done
  while (( changed )); do changed=0
    for n in "${ORDER[@]}"; do [[ -n "${in[$n]:-}" ]] && continue
      for d in ${DEPS[$n]}; do [[ -n "${in[$d]:-}" ]] && { in[$n]=1; changed=1; break; }; done; done; done
  for n in "${ORDER[@]}"; do if [[ -n "${in[$n]:-}" ]]; then echo "$n"; fi; done; }
# redo_plan: prints "kind<TAB>src<TAB>dst" for every existing dir / marker of the closure. A dir is moved only
# when every node owning it is in the redo set; a dir shared with a node outside it is KEPT, and that node's
# marker alone is moved, so the step reruns and resumes on its own sub-step markers (e.g. 07_<ds>: only the
# export marker of one dataset goes, the other datasets' outputs in 07_probe_matrices stay).
redo_plan(){ local ts="$1" n d o m sh; local -A owners=() inset=() seen=() moved=()
  for n in "${ORDER[@]}"; do for d in $(node_dirs "$n"); do owners[$d]+="$n "; done; done
  for n in $(redo_closure); do inset[$n]=1; done
  for n in "${!inset[@]}"; do for d in $(node_dirs "$n"); do
    [[ -n "${seen[$d]:-}" ]] && continue; seen[$d]=1
    sh=""; for o in ${owners[$d]}; do [[ -n "${inset[$o]:-}" ]] || { sh="$o"; break; }; done
    [[ -z "$sh" ]] || { echo "keep	$d	shared with $sh (outside the redo set): marker-only for its nodes" >&2; continue; }
    [[ -e "$d" ]] || continue; moved[$d]=1
    if [[ "$d" == "$SCR"/* ]]; then printf 'dir\t%s\t%s\n' "$d" "$SCR/run_all_superseded/$ts/${d#$SCR/}"
    else printf 'dir\t%s\t%s\n' "$d" "$DD/superseded/$ts/${d#$DUR/}"; fi
  done; done
  for n in "${!inset[@]}"; do m="${MARK[$n]}"; [[ -e "$m" ]] || continue
    for d in "${!moved[@]}"; do [[ "$m" == "$d"/* ]] && continue 2; done
    if [[ "$m" == "$SCR"/* ]]; then printf 'marker\t%s\t%s\n' "$m" "$SCR/run_all_superseded/$ts/${m#$SCR/}"
    else printf 'marker\t%s\t%s\n' "$m" "$DD/superseded/$ts/${m#$DUR/}"; fi
  done; }
REDO_TS="$(date -u +%Y%m%dT%H%M%SZ)"
if [[ -n "$REDO" ]]; then
  REDO_SET=$(redo_closure | paste -sd' '); REDO_MOVES=$(redo_plan "$REDO_TS" 2> "$DD/.redo_keep.$$") || exit 1
  REDO_KEEP=$(sort -u "$DD/.redo_keep.$$"); rm -f "$DD/.redo_keep.$$"
fi

# a stale node's done data dependents keep their outputs unless they are in --redo: say so loudly
stale_dependents(){ local n d o=""; local -A in=(); for n in $STALE; do in[$n]=1; done; local changed=1
  while (( changed )); do changed=0; for n in "${ORDER[@]}"; do [[ -n "${in[$n]:-}" ]] && continue
    for d in ${DEPS[$n]}; do [[ -n "${in[$d]:-}" ]] && { in[$n]=1; changed=1; break; }; done; done; done
  for n in "${ORDER[@]}"; do [[ -n "${in[$n]:-}" && " $STALE " != *" $n "* && "${STATE[$n]}" == done && " ${REDO_SET:-} " != *" $n "* ]] && o+="$n "; done
  echo "$o"; }
if (( DRY )); then
  [[ ! -s "$DD/check_12.txt" || "${STATE[12]}" == done ]] || { say "(dry run) step 12 --check:"; sed 's/^/    /' "$DD/check_12.txt"; }
  [[ -z "$STALE" ]] || say "(dry run) stale markers (exist but no longer valid under config.env; these nodes rerun): $STALE"
  [[ -z "$HOLD" ]] || say "(dry run) --hold $HOLD: these and their dependents will not start"
  SD_WARN=$(stale_dependents); [[ -z "$SD_WARN" ]] || say "(dry run) WARNING: done nodes downstream of a stale node keep their outputs (add them to --redo to rerun): $SD_WARN"
  if [[ -n "$REDO" ]]; then
    say "(dry run) --redo $REDO -> invalidates: $REDO_SET"
    [[ -n "$REDO_MOVES" ]] && printf '%s\n' "$REDO_MOVES" | while IFS=$'\t' read -r k a b; do say "(dry run)   would move $k $a -> $b"; done
    [[ -n "$REDO_KEEP" ]] && printf '%s\n' "$REDO_KEEP" | while IFS=$'\t' read -r k a b; do say "(dry run)   keep dir $a ($b)"; done
    for n in $REDO_SET; do STATE[$n]=pending; done
  fi
  preflight || say "(dry run) the real run would REFUSE to start until the FOREIGN entries above are resolved"
  # waves: a node's wave = 1 + max wave of its deps (done nodes are wave 0)
  declare -A WAVE
  for n in "${ORDER[@]}"; do
    w=0; [[ "${STATE[$n]}" == done ]] || { w=1; for d in $(alldeps "$n"); do (( WAVE[$d] + 1 > w )) && w=$(( WAVE[$d] + 1 )); done; }
    WAVE[$n]=$w
  done
  printf 'run_all plan (dry run; caps io=%s gene=%s quant=%s analysis=%s, timing alone; re-read every loop)\n' "${CAP[io]}" "${CAP[gene]}" "${CAP[quant]}" "${CAP[analysis]}"
  printf '%-5s %-22s %-8s %-9s %-40s %s\n' wave node class state deps command
  for n in "${ORDER[@]}"; do printf '%s\t%s\n' "${WAVE[$n]}" "$n"; done | sort -s -n -k1,1 | while IFS=$'\t' read -r w n; do
    dl="${DEPS[$n]}"; (( ${#dl} > 40 )) && dl="${dl:0:37}..."
    printf '%-5s %-22s %-8s %-9s %-40s %s\n' "$w" "$n" "${CLASS[$n]}" "${STATE[$n]}" "${dl:--}" "${CMD[$n]}"
  done
  exit 0
fi

# step processes writing to OUR $DUR (env DUR unset = config default = ours)
# (any checkout of the package, started by a relative or absolute path; a process without DUR in its
# environment writes to the config default)
FOREIGN=""; DEF_DUR=$(env -u DUR -u SCR bash -c 'source "$1/config.env"; echo "$DUR"' _ "$REPRO")
for p in $(pgrep -f "steps/[0-9]{2}[a-z]?_[a-z0-9_]+\.sh" || true); do
  pd=$(tr '\0' '\n' < "/proc/$p/environ" 2>/dev/null | sed -n 's/^DUR=//p')
  [[ "${pd:-$DEF_DUR}" == "$DUR" ]] && FOREIGN+="$p "
done
[[ -z "$FOREIGN" ]] || fail "step processes writing to $DUR already running outside run_all (pids $FOREIGN): wait for them or stop them first"
exec 8> "$DD/run_all.lock"; flock -n 8 || fail "another run_all holds $DD/run_all.lock (pid file: $(cat "$DD/run_all.pid" 2>/dev/null))"
REPRO_HEAD="$(git -C "$REPRO" rev-parse --short HEAD 2>/dev/null || echo unknown)"
printf 'pid\t%s\nstarted\t%s\nrepro_head\t%s\nargs\t%s\n' "$$" "$(date -u +%FT%TZ)" "$REPRO_HEAD" "$ARGS_ALL" > "$DD/run_all.pid"
trap 'rm -f "$DD/run_all.pid"' EXIT
exec > >(tee -a "$DD/run_all.log") 2>&1
[[ -z "$STALE" ]] || say "stale markers (exist but no longer valid under config.env; these nodes rerun): $STALE"
[[ ! -s "$DD/check_12.txt" || "${STATE[12]}" == done ]] || { say "step 12 --check:"; sed 's/^/    /' "$DD/check_12.txt"; }
SD_WARN=$(stale_dependents); [[ -z "$SD_WARN" ]] || say "WARNING: done nodes downstream of a stale node keep their outputs (add them to --redo to rerun): $SD_WARN"
if [[ -n "$REDO" ]]; then
  say "--redo $REDO -> invalidating: $REDO_SET (moved aside under superseded/$REDO_TS, nothing deleted)"
  [[ -n "$REDO_KEEP" ]] && printf '%s\n' "$REDO_KEEP" | while IFS=$'\t' read -r k a b; do say "  keep dir $a ($b)"; done
  if [[ -n "$REDO_MOVES" ]]; then
    while IFS=$'\t' read -r k a b; do
      mkdir -p "$(dirname "$b")"; mv "$a" "$b"
      printf 'superseded\t%s\t%s\t%s\t%s\t%s\n' "$a" "$b" "redo:$REDO" "$RUN_ID" "$(date -u +%FT%TZ)" >> "$LEDGER"
      say "  moved $k $a -> $b"
    done <<<"$REDO_MOVES"
  fi
  for n in $REDO_SET; do [[ ! -e "${MARK[$n]}" ]] || fail "--redo: marker ${MARK[$n]} of $n still present after the move"; STATE[$n]=pending; done
fi
preflight || fail "refusing to start: remove or account for the FOREIGN entries (they could be skipped into the paper)"
say "===== run_all START (pid $$): $(for n in "${ORDER[@]}"; do [[ ${STATE[$n]} == done ]] && echo -n "$n "; done | wc -w) of ${#ORDER[@]} nodes already done ====="
ready(){ local n="$1" d; [[ "${STATE[$n]}" == pending ]] || return 1
  [[ " $HOLD " != *" $n "* ]] || return 1
  for d in $(alldeps "$n"); do [[ "${STATE[$d]}" == done ]] || return 1; done; }
blocked(){ local n="$1" d; for d in $(alldeps "$n"); do [[ "${STATE[$d]}" == failed || "${STATE[$d]}" == blocked ]] && return 0; done; return 1; }
nrun(){ local c t=0; for c in "${!RUNNING_BY_CLASS[@]}"; do t=$(( t + RUNNING_BY_CLASS[$c] )); done; echo "$t"; }

while :; do
  changed=1
  while (( changed )); do changed=0
    for n in "${ORDER[@]}"; do
      if [[ "${STATE[$n]}" == pending ]] && blocked "$n"; then STATE[$n]=blocked; changed=1; say "[$n] BLOCKED by a failed dependency"; fi
    done
  done
  load_caps
  for n in "${ORDER[@]}"; do
    ready "$n" || continue
    c="${CLASS[$n]}"
    (( RUNNING_BY_CLASS[$c] < CAP[$c] )) || continue
    [[ "$c" == timing && "$(nrun)" != 0 ]] && continue     # timing starts on an otherwise idle driver
    (( RUNNING_BY_CLASS[timing] > 0 )) && continue          # and nothing starts while it runs
    say "[$n] START: ${CMD[$n]}"
    for d in $(node_dirs "$n"); do in_ledger dir "$d" || printf 'dir\t%s\t-\t%s\t%s\t%s\n' "$d" "$n" "$RUN_ID" "$(date -u +%FT%TZ)" >> "$LEDGER"; done
    ( set +e; bash -c "${CMD[$n]}" > "$LOGD/$n.log" 2>&1 < /dev/null; ec=$?
      if (( ec == 0 )) && [[ ! -f "${MARK[$n]}" ]]; then echo "run_all: exit 0 but marker ${MARK[$n]} missing" >> "$LOGD/$n.log"; ec=99; fi
      if (( ec == 0 )) && [[ -n "${VALID[$n]:-}" ]] && ! ${VALID[$n]}; then echo "run_all: exit 0 but marker ${MARK[$n]} is still stale (${VALID[$n]})" >> "$LOGD/$n.log"; ec=98; fi
      exit $ec ) &
    PIDN[$!]="$n"; STATE[$n]=running; RUNNING_BY_CLASS[$c]=$(( RUNNING_BY_CLASS[$c] + 1 ))
  done
  (( ${#PIDN[@]} > 0 )) || break
  set +e; wait -n -p FIN "${!PIDN[@]}"; ec=$?; set -e
  n="${PIDN[$FIN]}"; unset "PIDN[$FIN]"; c="${CLASS[$n]}"; RUNNING_BY_CLASS[$c]=$(( RUNNING_BY_CLASS[$c] - 1 ))
  if (( ec == 0 )); then
    m="${MARK[$n]}"
    [[ "$(basename "$m")" =~ ^(DONE|\.done|.*\.done|\..*_done)$ ]] && printf 'run_all_writer\t%s\t%s\t%s\n' "$RUN_ID" "$n" "$(date -u +%FT%TZ)" >> "$m"
    printf 'marker\t%s\t%s\t%s\t%s\t%s\n' "$m" "$(sha256sum "$m" | cut -c1-64)" "$n" "$RUN_ID" "$(date -u +%FT%TZ)" >> "$LEDGER"
    STATE[$n]=done; say "[$n] done"
  else STATE[$n]=failed; say "[$n] FAILED (exit $ec): tail of $LOGD/$n.log:"; tail -15 "$LOGD/$n.log" | sed 's/^/    /'; fi
done

# held: the --hold nodes and everything that only waits on them; anything else not done is a failure
held_closure(){ local -A in=(); local n d changed=1; for n in $HOLD; do in[$n]=1; done
  while (( changed )); do changed=0; for n in "${ORDER[@]}"; do [[ -n "${in[$n]:-}" ]] && continue
    for d in $(alldeps "$n"); do [[ -n "${in[$d]:-}" ]] && { in[$n]=1; changed=1; break; }; done; done; done
  for n in "${ORDER[@]}"; do if [[ -n "${in[$n]:-}" ]]; then echo "$n"; fi; done; }
HELD=" $(held_closure | paste -sd' ') "
bad=""; held=""; for n in "${ORDER[@]}"; do
  [[ "${STATE[$n]}" == done ]] && continue
  if [[ "${STATE[$n]}" == pending && "$HELD" == *" $n "* ]]; then held+="$n "; else bad+="$n(${STATE[$n]}) "; fi
done
[[ -z "$bad" ]] || fail "run_all incomplete: $bad${held:+ (held: $held)}"
if [[ -n "$held" ]]; then say "stopped at --hold: everything else is done; held (not started): $held"; exit 0; fi
LINES=(); for n in "${ORDER[@]}"; do LINES+=("$n	${MARK[$n]}"); done
write_done "${LINES[@]}"
