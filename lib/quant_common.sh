# Shared helpers for steps 05-07 (simpleaf timing, gene matrices, probe matrices).
# Source AFTER config.env and lib/common.sh. Defines:
#   jget JSON KEY                -> top-level key of a JSON file
#   index_of_ds DS               -> DS_PROBE DS_PANEL DS_IDX DS_T2G (no FASTQ checks; used by step 07)
#   resolve_ds DS                -> DS_CHEM DS_GEOM DS_CBC DS_SBC DS_PROBE DS_PANEL DS_IDX DS_T2G DS_R1 DS_R2 DS_NFQ
#   mq_args DS OUT [THREADS]     -> MQ_ARGS: the ONE multiplex-quant argv used by step 05 (timed) and 06 --untimed
#   check_ref_digest DS DIGEST L -> fail unless DIGEST == the dataset's row in $DIGEST_REF_TSV (per-probe reference)
#   check_mq_output DS OUT LOG   -> the per-run correctness gates (codec, geometry, umi-edit-dist, plan, matrix)
#   matrix_digest ALEVIN_DIR     -> "n_rows n_cols nnz total digest" (tab-separated, $REPRO/lib/matrix_digest.py)
#   idle_guard [load]            -> blocks (never kills) while a foreign heavy job runs; TIMING steps (05/05b) only
#   evict_check LABEL FILE...    -> evict FILEs from the page cache, fail unless residency < $EVICT_MAX_RESIDENT_PCT; sets RES
#   scr_inputs_check LABEL OUT PATH... -> fail unless every input (a dir: every file) resolves to $SCR's
#                                   filesystem; path -> device table in OUT; sets INPUT_SOURCE (e.g. /dev/md0)
#   stage_file SRC DST           -> byte-identical copy (resumable)
#   cr_multi_csv DS OUT BAM [nostage] -> the Cell Ranger multi CSV from data/cr_multi/<DS>.csv, every path on $SCR
SA_BIN="$SRC/simpleaf/target/release/simpleaf"      # = steps/01_toolchain.sh $SA
STAGED_TSV="$DUR/02_stage_inputs/datasets_staged.tsv"

jget(){ "$PY" -c "import json,sys;print(json.load(open(sys.argv[1]))[sys.argv[2]])" "$1" "$2"; }

_staged_row(){
  local row; row=$(awk -F'\t' -v d="$1" '$1==d' "$STAGED_TSV")
  [[ -n "$row" ]] || fail "[$1] not in $STAGED_TSV (run step 02)"
  echo "$row"
}

index_of_ds(){
  local ds="$1"
  DS_PROBE=$(_staged_row "$ds" | cut -f5)
  case "$(basename "$DS_PROBE")" in
    *Human*v1.1.0*) DS_PANEL="human_v1.1.0" ;;
    *Human*v2.0.0*) DS_PANEL="human_v2.0.0" ;;
    *Mouse*v2.0.0*) DS_PANEL="mouse_v2.0.0" ;;
    *) fail "[$ds] cannot map probe set $(basename "$DS_PROBE") to a panel of: $PANELS" ;;
  esac
  DS_IDX="$INDEX_ROOT/$DS_PANEL/index"; DS_T2G="$DS_IDX/t2g_2col.tsv"
  [[ -s "$DS_IDX/piscem_idx.ssi" && -s "$DS_T2G" ]] || fail "[$ds] index or t2g missing under $DS_IDX (run step 03)"
}

resolve_ds(){
  local ds="$1" row
  row=$(_staged_row "$ds")
  [[ -f "$SCR/02_stage_inputs/$ds.done" ]] || fail "[$ds] step 02 marker $SCR/02_stage_inputs/$ds.done missing (scratch lost? re-run step 02)"
  IFS=$'\t' read -r _ _ DS_CBC DS_SBC DS_PROBE DS_R1 DS_R2 <<<"$row"
  index_of_ds "$ds"
  DS_CHEM=$(awk -F'\t' -v d="$ds" '!/^#/ && $1==d {print $4}' "$DATASETS_TSV")
  case "$DS_CHEM" in
    v1) DS_GEOM="$GEOM_V1" ;;
    v2*) DS_GEOM="$GEOM_V2" ;;
    *) fail "[$ds] unknown chemistry '$DS_CHEM' in $DATASETS_TSV" ;;
  esac
  local f
  for f in ${DS_R1//,/ } ${DS_R2//,/ } "$DS_CBC" "$DS_SBC"; do
    [[ -s "$f" ]] || fail "[$ds] staged input missing: $f"
  done
  DS_NFQ=$(( $(tr ',' '\n' <<<"$DS_R1" | wc -l) + $(tr ',' '\n' <<<"$DS_R2" | wc -l) ))
}

mq_args(){
  local ds="$1" out="$2" t="${3:-$THREADS}"
  MQ_ARGS=( "$SA_BIN" multiplex-quant --threads "$t" --resolution "$RESOLUTION"
    --geometry "$DS_GEOM" --cell-bc-list "$DS_CBC" --sample-bc-list "$DS_SBC"
    --index "$DS_IDX"
    --cell-bc-correction "$CELL_BC_CORRECTION" --sample-bc-correction "$SAMPLE_BC_CORRECTION"
    --decoder "$DECODER" --compress "$COMPRESS" --umi-edit-dist "$UMI_EDIT_DIST"
    --reads1 "$DS_R1" --reads2 "$DS_R2" --output "$out" )
}

# check_mq_output DS OUT LOG -> sets ENGAGED MAP_THREADS DECODE_SLOTS NUM_MAPPED
check_mq_output(){
  local ds="$1" out="$2" log="$3" q="$2/af_quant"
  grep -qE -- "collate cmd: .*--compress[= ]$COMPRESS" "$log" || fail "[$ds] collate was not run with --compress $COMPRESS ($log)"
  [[ "$(jget "$q/collate.json" compressed_output)" == True && "$(jget "$q/collate.json" chunk_codec)" == "$COMPRESS" ]] \
    || fail "[$ds] $q/collate.json is not $COMPRESS-compressed"
  grep -qF -- "$DS_GEOM" "$log" || fail "[$ds] log lacks the configured geometry $DS_GEOM"
  local ued; ued=$("$PY" -c "import json,sys;print(json.load(open(sys.argv[1]))['quant_options']['crlike_umi_edit'])" "$q/quant.json")
  [[ "$ued" == "$UMI_EDIT_DIST" ]] || fail "[$ds] quant.json crlike_umi_edit=$ued != UMI_EDIT_DIST=$UMI_EDIT_DIST (--umi-edit-dist not forwarded?)"
  [[ -s "$q/map.collated.rad" && -s "$q/collation_manifest.bin" ]] || fail "[$ds] no collated RAD / manifest in $q"
  [[ -s "$q/alevin/quants_mat_rows.txt" ]] || fail "[$ds] no matrix under $q/alevin"
  read -r ENGAGED MAP_THREADS DECODE_SLOTS NUM_MAPPED < <("$PY" -c "
import json,sys; d=json.load(open(sys.argv[1])); p=d['execution_plan']
k=p['allocation']['kind']; s=p['decode_slots']
assert (k=='serial') == (s==0), (k, s)
print('serial' if k=='serial' else 'parallel_'+k, p['map_threads'], s, d['num_mapped'])" "$out/af_map/map_info.json") \
    || fail "[$ds] cannot parse execution_plan/num_mapped in $out/af_map/map_info.json"
}

matrix_digest(){ "$PY" "$REPRO/lib/matrix_digest.py" "$1"; }

# Timed runs (steps 05 and 05b), before EVERY run (warmup and reps) of every tool: the FASTQs are evicted
# and their residency proven below the limit; the reference / index files are not touched.
evict_check(){
  local label="$1"; shift
  "$PY" "$REPRO/lib/evict_cache.py" "$@" | tail -1
  RES=$("$PY" "$REPRO/lib/page_resident.py" "$@" | awk '/^TOTAL/{gsub(/[()%]/,"",$NF); print $NF}')
  awk -v r="$RES" -v m="$EVICT_MAX_RESIDENT_PCT" 'BEGIN{exit !(r != "" && r+0 < m)}' \
    || fail "[$label] eviction FAILED: ${RES}% of $# FASTQs still resident"
}
scr_inputs_check(){
  local label="$1" out="$2" o; shift 2
  mkdir -p "$(dirname "$out")"
  o=$("$PY" "$REPRO/lib/input_mounts.py" --scr "$SCR" --out "$out" "$@") \
    || fail "[$label] an input is not on the \$SCR filesystem (list above; table $out)"
  INPUT_SOURCE=$(awk '{print $5}' <<<"$o")
  say "[$label] inputs: $o"
}

check_ref_digest(){
  local ds="$1" dg="$2" label="$3" ref
  ref=$(awk -F'\t' -v d="$ds" 'NR>1 && $1==d {print $2}' "$DIGEST_REF_TSV")
  [[ -n "$ref" ]] || fail "[$ds] no reference digest in $DIGEST_REF_TSV"
  [[ "$dg" == "$ref" ]] || fail "[$ds] $label digest $dg != reference $ref ($DIGEST_REF_TSV)"
  say "[$ds] $label digest $dg == reference ($DIGEST_REF_TSV)"
}

# ---------------- idle guard: WAIT (never kill) while any foreign heavy job runs ----------------
# Excludes this shell's ancestor chain (a parent `bash -c '...'` whose text contains a pattern, e.g. the
# agent shell that launched the step) and its descendants (the step's own subshells), so only FOREIGN jobs block.
_ancestors(){ local p=$$; while [[ "$p" -gt 1 ]]; do echo "$p"; p=$(awk '/^PPid:/{print $2}' "/proc/$p/status"); done; }
_busy_reason(){
  local use_load="$1" cont procs ld anc dps pso
  # a failed probe counts as busy (reported), never as idle
  if ! dps=$(docker ps --format '{{.Names}}'); then echo "probe_error=docker_ps"; return 0; fi
  if ! pso=$(ps -eo pid=,ppid=,args=); then echo "probe_error=ps"; return 0; fi
  cont=$(RE="$IDLE_CONTAINER_RE" awk 'BEGIN{re=ENVIRON["RE"]} $0~re' <<<"$dps" | paste -sd,)
  anc=" $(_ancestors | paste -sd' ') "
  procs=$(RE="$IDLE_PROC_RE" ANC="$anc" ROOT="$$" awk '
    BEGIN { re = ENVIRON["RE"]; anc = ENVIRON["ANC"]; root = ENVIRON["ROOT"] }
    { pp[$1] = $2; a = $0; sub(/^ *[0-9]+ +[0-9]+ /, "", a); args[$1] = a }
    END { for (p in args) {
            if (args[p] !~ re || index(anc, " " p " ")) continue
            q = p; mine = 0; while ((q in pp) && q > 1) { if (q == root) { mine = 1; break }; q = pp[q] }
            if (!mine) printf "%s:%s; ", p, substr(args[p], 1, 80) } }' <<<"$pso")
  ld=$(awk '{print $1}' /proc/loadavg)
  if [[ -n "$cont" ]]; then echo "containers=$cont"
  elif [[ -n "$procs" ]]; then echo "procs=$procs"
  elif [[ "$use_load" == load ]] && awk -v l="$ld" -v m="$IDLE_LOAD_MAX" 'BEGIN{exit !(l >= m)}'; then echo "load1=$ld>=$IDLE_LOAD_MAX"
  fi
}
idle_guard(){
  local waited=0 r
  while r=$(_busy_reason "${1:-noload}"); [[ -n "$r" ]]; do
    (( waited % IDLE_LOG_EVERY_S == 0 )) && say "idle_guard: waiting ${waited}s ($r)"
    sleep "$IDLE_POLL_S"; waited=$((waited + IDLE_POLL_S))
  done
  (( waited > 0 )) && say "idle_guard: box idle after ${waited}s"
  return 0
}

# cgroup_probe OUT CPUS CPUSET MEM: records whether docker on this host ENFORCES the timing caps. On this
# host (cgroup v1, cgroupfs, no cpuset/memory controller mounted) docker prints "Cpuset discarded" /
# "Limitation discarded" and only --cpus (CFS quota) applies. Recorded, not fatal; the file is
# listed in the step's DONE.
cgroup_probe(){
  local out="$1" cpus="$2" cpuset="$3" mem="$4" o
  o=$(docker run --rm --cpus "$cpus" --cpuset-cpus "$cpuset" --memory "$mem" --memory-swap "$mem" "$TIMING_IMAGE" \
      bash -c 'grep Cpus_allowed_list /proc/self/status; nproc' 2>&1)
  {
    echo -e "requested_cpus\t$cpus"; echo -e "requested_cpuset\t$cpuset"; echo -e "requested_mem\t$mem"
    echo -e "container_cpus_allowed\t$(awk '/Cpus_allowed_list/{print $2}' <<<"$o")"
    echo -e "container_nproc\t$(tail -1 <<<"$o")"
    echo -e "cpuset_enforced\t$(grep -q 'Cpuset discarded' <<<"$o" && echo no || echo yes)"
    echo -e "memory_limit_enforced\t$(grep -q 'Limitation discarded' <<<"$o" && echo no || echo yes)"
    echo -e "probed_at\t$(date -u +%FT%TZ)"
  } > "$out"
  if grep -q $'\tno$' "$out"; then
    say "WARNING: docker does NOT enforce every timing cap on this host (only --cpus quota applies): $(grep $'\tno$' "$out" | cut -f1 | paste -sd,) -> $out"
  fi
  return 0
}

# ---------------- Cell Ranger inputs on $SCR (steps 05b and 13 bam) ----------------
# stage_file SRC DST: copy (via .part) unless DST already equals SRC byte for byte
stage_file(){
  if [[ ! -f "$2" ]] || ! cmp -s "$1" "$2"; then mkdir -p "$(dirname "$2")"; cp -L "$1" "$2.part"; mv "$2.part" "$2"; fi
  cmp -s "$1" "$2" || fail "staged $2 differs from $1"
}
# cr_multi_csv DS OUT CREATE_BAM(true|false) [nostage] -> writes OUT from data/cr_multi/<DS>.csv; sets REFD CSV_PATHS.
#   reference -> $SCR/refs/cr_ref/<name> (staged from $REFS_SRC/cr_ref unless "nostage"), probe-set -> the step-02
#   staged probe CSV, fastqs (@STAGED@) -> the step-02 staged dir of each fastq_id, create-bam -> CREATE_BAM.
#   Needs resolve_ds DS first. Fails unless every absolute path in OUT is under $SCR.
cr_multi_csv(){
  local ds="$1" out="$2" bam="$3" stage="${4:-stage}" tpl="$REPRO/data/cr_multi/$1.csv" line fid rest dir inlibs=0 refb p
  [[ -s "$tpl" ]] || fail "[$ds] no Cell Ranger template $tpl"
  [[ "$bam" == true || "$bam" == false ]] || fail "cr_multi_csv: create-bam must be true|false, not '$bam'"
  refb=$(awk -F, 'tolower($1)=="reference"{print $2}' "$tpl" | tr -d ' \r')
  [[ -n "$refb" && "$refb" != */* ]] || fail "[$ds] template reference '$refb' is not a bare reference name"
  [[ "$(awk -F, 'tolower($1)=="probe-set"{print $2}' "$tpl" | tr -d ' \r')" == "$(basename "$DS_PROBE")" ]] \
    || fail "[$ds] template probe-set != staged $(basename "$DS_PROBE")"
  grep -qx 'create-bam,false' "$tpl" || fail "[$ds] template has no 'create-bam,false' line"
  REFD="$SCR/refs/cr_ref/$refb"
  if [[ "$stage" == stage && ! -f "$REFD/.staged" ]]; then
    say "[$ds] staging Cell Ranger reference $refb -> $REFD"
    rm -rf "$REFD.part"; mkdir -p "$(dirname "$REFD")"; cp -rL "$REFS_SRC/cr_ref/$refb" "$REFD.part"
    rm -rf "$REFD"; mv "$REFD.part" "$REFD"; touch "$REFD/.staged"
  fi
  mkdir -p "$(dirname "$out")"; : > "$out"
  while IFS= read -r line; do
    line="${line%$'\r'}"
    [[ "$line" == "#"* ]] && continue
    case "$line" in
      "[libraries]"*) inlibs=1; echo "$line" >> "$out"; continue;;
      "["*) inlibs=0;;
    esac
    if [[ "$inlibs" == 1 && -n "${line// /}" && ! "$line" =~ ^fastq_id ]]; then
      fid="${line%%,*}"; rest=$(cut -d, -f3- <<<"$line")
      [[ "$(cut -d, -f2 <<<"$line")" == "@STAGED@" ]] || fail "[$ds] template fastqs of $fid is not @STAGED@"
      dir=$(tr ',' '\n' <<<"$DS_R1" | awk -v f="/${fid}_S" 'index($0, f) { sub(/\/[^\/]*$/, ""); print }' | sort -u)
      [[ $(wc -l <<<"$dir") == 1 && -n "$dir" ]] || fail "[$ds] fastq_id $fid: staged dir not unique/found ('$dir')"
      echo "$fid,$dir,$rest" >> "$out"; continue
    fi
    case "$line" in
      reference,*) line="reference,$REFD";;
      probe-set,*) line="probe-set,$DS_PROBE";;
      create-bam,*) line="create-bam,$bam";;
    esac
    echo "$line" >> "$out"
  done < "$tpl"
  # content check: every absolute path the CSV names is under $SCR (scr_inputs_check proves the filesystem)
  CSV_PATHS=$(tr ',' '\n' < "$out" | grep '^/' | sort -u)
  [[ -n "$CSV_PATHS" ]] || fail "[$ds] $out names no path"
  p=$(grep -v "^$SCR/" <<<"$CSV_PATHS" || true)
  [[ -z "$p" ]] || fail "[$ds] $out names paths outside $SCR: $p"
  grep -qx "reference,$REFD" "$out" && grep -qx "probe-set,$DS_PROBE" "$out" && grep -qx "create-bam,$bam" "$out" \
    || fail "[$ds] $out lacks the staged reference / probe set / create-bam,$bam"
}
