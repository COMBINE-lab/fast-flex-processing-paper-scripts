#!/usr/bin/env bash
# 03: build the three Flex probe-set indices (human v1.1.0, human v2.0.0, mouse v2.0.0) with the
#   step-01 simpleaf/piscem: simpleaf index --probe-csv CSV -k $KMER -m $MINIMIZER
#   --excluded-probes $EXCLUDED_PROBES -t $INDEX_THREADS (the build_idx19.log command; simpleaf adds
#   --seed 1 --dict auto --ram-limit-gib 8 --decoy-paths itself, asserted below).
#   Output $SCR/index/<panel>, durable copy $DUR/03_index/<panel>.
#   Verifies: piscem build cmd flags, poison table present, piscem_idx_ver.json = piscem-rs
#   $PISCEM_RS_VERSION, and t2g identical to the index in $PRIOR_INDEX_ROOT.
#   Generates index/t2g_2col.tsv (cut -f1,2 ref/t2g_3col.tsv, rows == probes in ref.fa): the
#   --t2g-map file every later step uses (sha256 in $DUR/03_index/<panel>.sha256).
#   The prior indices are compared on text files only: piscem-rs 0.11 cannot read them (PECTB03).
#   piscem-rs index binaries differ byte-wise between builds, so a lost $SCR copy is restored from
#   $DUR (sha256-checked), never rebuilt; delete $DUR/03_index/<panel> to force a rebuild (then re-gate).
# Usage: 03_index.sh [--panel P]   (default: all of $PANELS). Resumable (per-panel .done).
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"
SA="$SRC/simpleaf/target/release/simpleaf"
[[ -f "$DUR/01_toolchain/DONE" ]] || fail "step 01 not done"

PL="$PANELS"
while [[ $# -gt 0 ]]; do case "$1" in
  --panel) PL="$2"; shift 2;;
  *) fail "unknown arg: $1";;
esac; done

# panel name from a probe CSV basename: Chromium_Human_..._v2.0.0_... -> human_v2.0.0
panel_of(){ basename "$1" | sed -nE 's/.*_(Human|Mouse)_.*_(v[0-9]+\.[0-9]+\.[0-9]+)_.*/\L\1\E_\2/p'; }
probe_csv_for(){ # first probe CSV in $DATASETS_TSV (column 8) of panel $1
  local p; for p in $(awk -F'\t' '!/^#/ {print $8}' "$DATASETS_TSV"); do
    [[ "$(panel_of "$p")" == "$1" && -s "$p" ]] && { echo "$p"; return 0; }
  done; return 1
}

mkdir -p "$SCR/index"
for panel in $PL; do
  out="$SCR/index/$panel"; dst="$DD/$panel"
  if [[ -f "$out/.done" && -f "$dst/.done" ]]; then say "[$panel] already built -> skip"; continue; fi
  # piscem index binaries are not byte-deterministic across -t N builds: once a build is in $DUR,
  # restore THAT build (the one step 04 gated) instead of rebuilding after a scratch loss
  if [[ -f "$dst/.done" ]]; then
    rm -rf "$out"; cp -r "$dst" "$out"; rm -f "$out/.done"
    ( cd "$out" && sha256sum -c --quiet "$DD/$panel.sha256" ) || fail "[$panel] restored copy fails $DD/$panel.sha256"
    date -u +%FT%TZ > "$out/.done"; say "[$panel] restored from $dst (sha256 verified)"; continue
  fi
  csv=$(probe_csv_for "$panel") || fail "[$panel] no probe CSV in $DATASETS_TSV"
  rm -rf "$out"
  say "[$panel] simpleaf index -k $KMER -m $MINIMIZER --excluded-probes $EXCLUDED_PROBES -t $INDEX_THREADS --probe-csv $csv"
  ( cd "$SCR/index" && ALEVIN_FRY_HOME="$AF_HOME" "$SA" index --probe-csv "$csv" -o "$out" \
      -k "$KMER" -m "$MINIMIZER" --excluded-probes "$EXCLUDED_PROBES" -t "$INDEX_THREADS" ) > "$LOGD/index_$panel.log" 2>&1 \
    || { tail -20 "$LOGD/index_$panel.log"; fail "[$panel] simpleaf index failed"; }

  ix="$out/index"
  # simpleaf writes only ref/t2g_3col.tsv; --t2g-map takes the 2-column probe->gene map, generated here
  cut -f1,2 "$out/ref/t2g_3col.tsv" > "$ix/t2g_2col.tsv"
  nprobe=$(grep -c '^>' "$out/ref/ref.fa")
  [[ "$(wc -l < "$ix/t2g_2col.tsv")" == "$nprobe" ]] || fail "[$panel] t2g_2col rows $(wc -l < "$ix/t2g_2col.tsv") != $nprobe probes in ref.fa"
  cmd=$("$PY" -c "import json,sys;print(json.load(open(sys.argv[1]))['cmd'])" "$ix/simpleaf_index.json")
  for want in "-k $KMER " "-m $MINIMIZER " "--seed 1 " "--threads $INDEX_THREADS " "--decoy-paths $out/ref/probe_decoys.fa" "$SRC/piscem/target/release/piscem build"; do
    [[ "$cmd" == *"$want"* ]] || fail "[$panel] piscem build cmd lacks '$want': $cmd"
  done
  [[ -s "$out/ref/probe_decoys.fa" && -s "$ix/piscem_idx.poison" && -s "$ix/piscem_idx.poison.json" ]] || fail "[$panel] decoy FASTA / poison table missing"
  for fm in $INDEX_MAGICS; do
    mg=$(head -c 7 "$ix/${fm%%=*}"); [[ "$mg" == "${fm#*=}" ]] || fail "[$panel] ${fm%%=*} magic '$mg' != ${fm#*=}"
  done
  ver=$("$PY" -c "import json,sys;print(json.load(open(sys.argv[1])).get('piscem-rs',''))" "$ix/piscem_idx_ver.json")
  [[ "$ver" == "$PISCEM_RS_VERSION" ]] || fail "[$panel] piscem_idx_ver.json piscem-rs='$ver' != $PISCEM_RS_VERSION"
  cmp -s "$ix/t2g_2col.tsv" "$PRIOR_INDEX_ROOT/$panel/index/t2g_2col.tsv" || fail "[$panel] t2g_2col.tsv differs from $PRIOR_INDEX_ROOT/$panel"
  cmp -s "$out/ref/probe_decoys.fa" "$PRIOR_INDEX_ROOT/$panel/ref/probe_decoys.fa" || fail "[$panel] probe_decoys.fa differs from $PRIOR_INDEX_ROOT/$panel"
  say "[$panel] OK: $(du -sh "$out" | cut -f1); genes=$(cut -f2 "$ix/t2g_2col.tsv" | sort -u | wc -l); ver=$(tr -d ' \n' < "$ix/piscem_idx_ver.json"); poison=$(tr -d ' \n' < "$ix/piscem_idx.poison.json") (prior: $(tr -d ' \n' < "$PRIOR_INDEX_ROOT/$panel/index/piscem_idx.poison.json"))"
  date -u +%FT%TZ > "$out/.done"
  rm -rf "$dst"; cp -r "$out" "$dst"
  diff -rq "$out" "$dst" > /dev/null || fail "[$panel] durable copy differs"
  ( cd "$dst" && find . -type f ! -name .done -print0 | sort -z | xargs -0 sha256sum ) > "$DD/$panel.sha256"
  date -u +%FT%TZ > "$dst/.done"
done

ALLDONE=1; for p in $PANELS; do [[ -f "$DD/$p/.done" ]] || ALLDONE=0; done
if (( ALLDONE )); then
  LINES=(); for p in $PANELS; do LINES+=("$p.sha256	$(sha256sum "$DD/$p.sha256" | cut -c1-64)"); done
  write_done "${LINES[@]}"
fi
