# Shared plumbing for steps 08-11. Source AFTER config.env and BEFORE lib/common.sh.
#
# The analysis code is lib/ (compare_consistency.py, scbp_pipeline.py, deg_overlap.py, ...).
# It finds every matrix through flexio.py's layout under $FLEX_LAYOUT_ROOT:
#   concordance/quant/<ds>/simpleaf/af_quant   concordance/quant/<ds>/cyto   concordance/cellranger/<ds>
#   concordance/clustering                      (per-sample scbp objects)
# so each step builds that layout as a tree of symlinks under its own scratch dir (bench_root) and
# points FLEX_LAYOUT_ROOT at it. af_quant/alevin is a REAL directory of per-file links: flexio writes its
# per-sample split cache (_by_sample/, for matrices > 2 GiB) next to the matrix, and that cache must land
# in this step's tree, not inside step 06's durable output.
#
# Test overrides (environment). Tests set DUR and SCR (config.env honours them) to a scratch dir; a step whose
# inputs are overridden REFUSES to run against the real (config default) $DUR, so a test can never leave
# outputs or .done markers in the paper tree.
#   OVR_SA_QUANT_ROOT, OVR_SA_QUANT_SUBDIR, OVR_CR_OUT, OVR_CR_SUBPATH, OVR_CYTO_OUT   input roots
#   OVR_CLUST_DIR=DIR      step 09's durable dir, as read by steps 10/11
#   OVR_SEEDCTL_DATASETS="a b"  datasets that get a seed control (step 09)
#   OVR_DATASETS="a b"     dataset list;  OVR_SAMPLES="ds:sample ds:sample"  restrict to these samples
SA_QUANT_ROOT="${OVR_SA_QUANT_ROOT:-$SA_QUANT_ROOT}"
SA_QUANT_SUBDIR="${OVR_SA_QUANT_SUBDIR-$SA_QUANT_SUBDIR}"
CR_OUT="${OVR_CR_OUT:-$CR_OUT}"; CR_SUBPATH="${OVR_CR_SUBPATH-$CR_SUBPATH}"
CYTO_OUT="${OVR_CYTO_OUT:-$CYTO_OUT}"
CLUST_DIR="${OVR_CLUST_DIR:-$DUR/09_clustering}"
DATASETS="${OVR_DATASETS:-$DATASETS}"
SEEDCTL_DATASETS="${OVR_SEEDCTL_DATASETS:-$SEEDCTL_DATASETS}"
# Same for a non-default clustering variant (CLUST_REMOVE_DOUBLETS / CLUST_N_DEVIANT from the environment).
_cfg(){ env -u DUR -u SCR -u CLUST_REMOVE_DOUBLETS -u CLUST_N_DEVIANT bash -c "source '$REPRO/config.env' && echo \"\$$1\""; }
if env | grep -q '^OVR_' || [[ "$CLUST_REMOVE_DOUBLETS" != "$(_cfg CLUST_REMOVE_DOUBLETS)" || "$CLUST_N_DEVIANT" != "$(_cfg CLUST_N_DEVIANT)" ]]; then
  if [[ "$(realpath -m "$DUR")" == "$(realpath -m "$(_cfg DUR)")"* ]]; then
    echo "FATAL: OVR_* input overrides or a non-default clustering variant are set but DUR=$DUR is the real deliverable tree; set DUR and SCR to a scratch dir" >&2
    exit 1
  fi
fi
export FLEX_DATASETS_TSV="$DATASETS_TSV"

# sa_dir DS -> step-06 af_quant dir of DS
sa_dir(){ echo "$SA_QUANT_ROOT/$1${SA_QUANT_SUBDIR:+/$SA_QUANT_SUBDIR}"; }
cr_dir(){ echo "$CR_OUT/$1${CR_SUBPATH:+/$CR_SUBPATH}"; }

# samples DS -> sample names (flexio order), restricted by OVR_SAMPLES when set
samples(){
  local s
  for s in $("$PY" -c "import sys;sys.path.insert(0,'$REPRO/lib');import flexio
print('\n'.join(n for n,_ in flexio.samples_of('$1')))"); do
    if [[ -z "${OVR_SAMPLES:-}" || " $OVR_SAMPLES " == *" $1:$s "* ]]; then echo "$s"; fi
  done
}

# build_bench_root ROOT [clustering-dir] -> symlink layout for $DATASETS; exports FLEX_LAYOUT_ROOT
build_bench_root(){
  local R="$1" CL="${2:-}" ds src f
  mkdir -p "$R/concordance/quant" "$R/concordance/cellranger"
  for ds in $DATASETS; do
    src="$(sa_dir "$ds")"
    [[ -d "$src/alevin" ]] || { echo "FATAL: no simpleaf quant $src/alevin" >&2; return 1; }
    [[ -d "$CYTO_OUT/$ds/cyto/counts" ]] || { echo "FATAL: no cyto counts $CYTO_OUT/$ds/cyto/counts" >&2; return 1; }
    [[ -d "$(cr_dir "$ds")" ]] || { echo "FATAL: no Cell Ranger dir $(cr_dir "$ds")" >&2; return 1; }
    local q="$R/concordance/quant/$ds/simpleaf/af_quant"
    mkdir -p "$q/alevin"
    for f in "$src"/*; do [[ "$(basename "$f")" == alevin ]] || ln -sfn "$f" "$q/$(basename "$f")"; done
    for f in "$src"/alevin/quants_mat*; do ln -sfn "$f" "$q/alevin/$(basename "$f")"; done
    clear_stale_split_lock "$q/alevin"
    # reuse a COMPLETE split cache if the source already has one (320k: rebuilding it takes long)
    if [[ -e "$src/alevin/_by_sample/_done" && ! -e "$q/alevin/_by_sample" ]]; then
      ln -sfn "$src/alevin/_by_sample" "$q/alevin/_by_sample"
    fi
    ln -sfn "$CYTO_OUT/$ds/cyto" "$R/concordance/quant/$ds/cyto"
    ln -sfn "$(cr_dir "$ds")" "$R/concordance/cellranger/$ds"
  done
  if [[ -n "$CL" ]]; then mkdir -p "$CL"; ln -sfn "$CL" "$R/concordance/clustering"; fi
  export FLEX_LAYOUT_ROOT="$R"
}

# clear_stale_split_lock ALEVIN_DIR: flexio._ensure_split takes an O_EXCL lock file (holding its pid) while it
# builds _by_sample/ and removes it on exit, but a SIGKILLed builder leaves it behind, and every later reader then
# waits 7200 s and dies. Remove the lock (and the builder's private _by_sample.tmp<pid>) only when the cache is
# incomplete and the pid in the lock is not alive; a live builder is left alone.
clear_stale_split_lock(){
  local a="$1" lk="$1/_by_sample.lock" pid
  [[ -f "$lk" && ! -e "$a/_by_sample/_done" ]] || return 0
  pid="$(tr -dc '0-9' < "$lk")"
  if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then return 0; fi
  echo "clearing stale split lock $lk (pid ${pid:-none} not running)" >&2
  rm -f "$lk"; rm -rf "$a"/_by_sample.tmp*
}

# preflight: every dataset's three arms resolve to the configured inputs (Cell Ranger: the per-sample dir
# $CR_OUT/<ds>/$CR_SUBPATH/<sample>, i.e. step 05b rep1), and every sample has all three
preflight(){
  "$PY" - "$DATASETS" "$SA_QUANT_ROOT" "${OVR_SAMPLES:-}" "$CR_OUT" "$CR_SUBPATH" <<'PYEOF'
import os, sys
sys.path.insert(0, os.path.join(os.environ["REPRO"], "lib"))
import flexio
dsl, sa_root, only = sys.argv[1].split(), os.path.realpath(sys.argv[2]), sys.argv[3].split()
cr_root, cr_sub = sys.argv[4], sys.argv[5]
bad = []
for ds in dsl:
    d = flexio.simpleaf_quant_dir(ds)
    mtx = flexio._find(os.path.join(d, "alevin"), "quants_mat.mtx*") if d else None
    if not mtx or not os.path.realpath(mtx).startswith(sa_root + os.sep):
        bad.append(f"{ds}: simpleaf arm {d} -> {mtx and os.path.realpath(mtx)} not under {sa_root}")
        continue
    for s, _ in flexio.samples_of(ds):
        if only and f"{ds}:{s}" not in only:
            continue
        miss = [t for t in ("cellranger", "simpleaf", "cyto") if not flexio.has_tool(ds, s, t)]
        if miss:
            bad.append(f"{ds}/{s}: missing {miss}")
            continue
        cr = os.path.realpath(flexio.cellranger_path(ds, s))
        if not cr.startswith(os.path.realpath(os.path.join(cr_root, ds, cr_sub)) + os.sep):
            bad.append(f"{ds}/{s}: Cell Ranger resolves to {cr}, not under $CR_OUT/{ds}/$CR_SUBPATH")
for b in bad:
    print("  PREFLIGHT !! " + b)
sys.exit(1 if bad else 0)
PYEOF
}
export REPRO

# Background jobs: bg_job TAG LOG cmd... runs cmd in the background, at most $CONC_JOBS at once. A failing
# job appends "TAG (LOG)" to $FAILS (set by the caller) instead of aborting its siblings; check_fails
# then exits non-zero listing every failure, so nothing is dropped silently.
bg_job(){
  local tag="$1" log="$2"; shift 2
  while (( $(jobs -rp | wc -l) >= CONC_JOBS )); do wait -n; done
  ( "$@" > "$log" 2>&1 || echo "$tag ($log)" >> "$FAILS" ) &
}
check_fails(){
  wait
  if [[ -s "$FAILS" ]]; then echo "FATAL: $(wc -l < "$FAILS") job(s) failed:" >&2; cat "$FAILS" >&2; exit 1; fi
}

# ensure_ext_env: $EXT_ENV (venv over scbp + celltypist $CELLTYPIST_VERSION / harmonypy $HARMONYPY_VERSION) is
# BUILT by step 00 (lib/ext_env.sh build_ext_env) from envs/scbp_ext/{freeze.txt,
# scbp_constraints.txt}; steps only check it here, so a run can never resolve different package versions.
ensure_ext_env(){
  local want="celltypist==$CELLTYPIST_VERSION harmonypy==$HARMONYPY_VERSION"
  [[ -x "$EXT_PY" && "$(cat "$EXT_ENV/.done" 2>/dev/null)" == "$want" ]] \
    || { echo "FATAL: $EXT_ENV is not built for '$want' (run steps/00_scratch_and_env.sh)" >&2; return 1; }
  diff -q <("$EXT_PY" -m pip freeze 2>/dev/null) "$EXT_ENV/freeze.txt" > /dev/null \
    || { echo "FATAL: $EXT_ENV pip freeze differs from $EXT_ENV/freeze.txt (re-run step 00)" >&2; return 1; }
}
