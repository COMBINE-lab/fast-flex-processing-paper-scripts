#!/usr/bin/env bash
# 13: read-level accounting of the residual simpleaf-vs-Cell Ranger discordance on $READLEVEL_DS
#     (supplement "Sources of residual discordance").
#   bam  Cell Ranger multi with create-bam true on $READLEVEL_DS -> $CR_BAM_OUT (untimed, whole host, inputs on $SCR).
#   a  no re-quant: reads step 06's durable mouse run ($DUR/06_gene_matrices/<ds>: af_quant/{alevin,map_info.json,
#      sample_info.json,map.collated.rad} + af_map/unmapped_bc_count.bin); its matrix digest must equal step 06's
#      digest_ref.tsv and the $DIGEST_REF_TSV row (per-probe reference). Stage-wise read ledger from map_info.json +
#      unmapped_bc_count.bin + sample_info.json; every identity asserted, closes to 0.     (seconds)
#   b  lib/r1_geom_scan.c over every R1 read with the window/anchor/tolerance parsed from GEOM_V2;
#      its drop must equal a's  num_reads - num_mapped - num_poisoned - ledger  exactly.  (~5 min)
#   c  R2 rules, both directions: one pass per CR BAM ($CR_BAM_OUT) writes each CR-assigned read as a
#      synthetic pair (R1 cell barcode = read ordinal) that the step-01 piscem maps with the step-03
#      index; placements are joined back by ordinal (exact piscem, plus the k-mer proxy at K=$KMER).
#      Each BAM is split into $READLEVEL_CHUNKS coordinate chunks; $READLEVEL_JOBS workers run at once.
#      The same pass dumps CR's per-read (cell, raw UMI, gene) molecules for d.  (~1 h; run detached)
#   d  gene reassignment: lib/rad_dump mol on a's lz4 collated RAD vs c's CR molecules.   (~15 min)
#   e  UMI totals per tool (summed from step 08's per-sample JSONs) + step 08's mouse means (text_ranges.json).
#   f  Cell Ranger 'Valid probe barcodes', and proof that $CR_BAM_OUT is the same CR result as $CR_OUT (step 05b rep1).
#   numbers  $DD/readlevel_numbers.json (value + definition + basis per number) and DONE.
# RAD reading: piscem-rs 0.11.0 writes RAD v2 and the final collate is lz4 ($COMPRESS). Instead of
#   collating uncompressed, lib/rad_dump is built against
#   libradicl 0.21.0 (= alevin-fry edfd1aa) and reads the codec from the file tags.
# Usage: 13_readlevel.sh [--step bam|a|b|c|d|e|f|numbers|all] [--detach]     Resumable per sub-step.
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"
source "$REPRO/lib/quant_common.sh"
SUB="all"; DETACH=0
while [[ $# -gt 0 ]]; do case "$1" in
  --step) SUB="$2"; shift 2;;
  --detach) DETACH=1; shift;;
  *) fail "unknown arg: $1";;
esac; done
if (( DETACH )); then
  setsid nohup bash "$0" --step "$SUB" > "$LOGD/detached_${SUB}.log" 2>&1 < /dev/null &
  say "detached (pid $!): poll $LOGD/detached_${SUB}.log and $SD/<sub>.done"; exit 0
fi
DS="$READLEVEL_DS"
PISCEM="$SRC/piscem/target/release/piscem"
mkdir -p "$SD"
exec > >(tee -a "$LOGD/13_readlevel_${SUB}.log") 2>&1

resolve_ds "$DS"
REF_FA="$(dirname "$DS_IDX")/ref/ref.fa"
[[ -s "$REF_FA" ]] || fail "probe reference $REF_FA missing (step 03)"
mapfile -t SAMPLES < <(cut -f3 "$DS_SBC")
SBC0=$(head -1 "$DS_SBC" | cut -f1)
(( ${#SAMPLES[@]} > 0 )) || fail "no samples in $DS_SBC"
need_bams(){ local s
  for s in "${SAMPLES[@]}"; do
    [[ -s "$CR_BAM_OUT/per_sample_outs/$s/sample_alignments.bam" ]] || fail "no CR BAM for $s under $CR_BAM_OUT (run --step bam)"
  done; }
eval "$("$PY" "$REPRO/lib/readlevel_ledger.py" geom "$DS_GEOM")"   # CB_LEN UMI_LEN SP_LO SP_HI ANCHOR MAXD SBC_LEN
D6="$DUR/06_gene_matrices/$DS"; SA_Q="$D6/af_quant"; UNMAPPED="$D6/af_map/unmapped_bc_count.bin"   # step 06 = the paper's run

tools(){
  if [[ -f "$SD/tools.done" && -x "$SD/bin/r1_geom_scan" && -x "$SD/build/release/rad_dump" ]]; then return 0; fi
  mkdir -p "$SD/bin"
  cc -O2 -Wall -o "$SD/bin/r1_geom_scan" "$REPRO/lib/r1_geom_scan.c"
  ( cd "$REPRO/lib/rad_dump" && CARGO_TARGET_DIR="$SD/build" cargo build --release --locked ) > "$LOGD/build_rad_dump.log" 2>&1 \
    || { tail -20 "$LOGD/build_rad_dump.log"; fail "cargo build lib/rad_dump failed"; }
  grep -A1 '^name = "libradicl"$' "$REPRO/lib/rad_dump/Cargo.lock" | grep -q 'version = "0.21.0"' || fail "rad_dump not on libradicl 0.21.0"
  touch "$SD/tools.done"
}

# bam: Cell Ranger $CELLRANGER_VERSION multi WITH BAMs on $DS (untimed: --localcores $CR_BAM_THREADS --localmem
#   $CR_BAM_LOCALMEM, no container). Same inputs as step 05b's Cell Ranger run (FASTQs, staged probe set and
#   reference on $SCR; CSV from data/cr_multi/<ds>.csv with create-bam,true). Runs on $SD, then outs/ -> $CR_BAM_OUT.
step_bam(){
  local dst; dst="$(dirname "$CR_BAM_OUT")"
  if [[ -f "$dst/.done" ]]; then say "[bam] done -> skip"; return 0; fi
  local csv="$SD/cr_multi_${DS}_bam.csv" run="$SD/cr_mouse_bam" s
  cr_multi_csv "$DS" "$csv" true
  local -a P; mapfile -t P <<<"$CSV_PATHS"
  scr_inputs_check "$DS bam" "$SD/bam_input_mounts.tsv" ${DS_R1//,/ } ${DS_R2//,/ } "$csv" "${P[@]}"
  rm -rf "$run"
  say "[bam] cellranger multi --localcores $CR_BAM_THREADS --localmem $CR_BAM_LOCALMEM (create-bam true) -> $run"
  ( cd "$SD" && "$CELLRANGER_BIN" multi --id cr_mouse_bam --csv "$csv" --localcores "$CR_BAM_THREADS" \
      --localmem "$CR_BAM_LOCALMEM" --disable-ui ) > "$LOGD/cr_bam.log" 2>&1 || { tail -30 "$LOGD/cr_bam.log"; fail "[bam] cellranger failed"; }
  for s in "${SAMPLES[@]}"; do
    [[ -s "$run/outs/per_sample_outs/$s/sample_alignments.bam" ]] || fail "[bam] no BAM for $s"
  done
  rm -rf "$dst.part"; mkdir -p "$dst.part"
  cp -r "$run/outs" "$dst.part/outs"; cp "$csv" "$SD/bam_input_mounts.tsv" "$dst.part/"
  rm -rf "$dst"; mv "$dst.part" "$dst"; rm -rf "$run"
  date -u +%FT%TZ > "$dst/.done"; say "[bam] $(ls "$CR_BAM_OUT/per_sample_outs" | wc -l) samples with BAMs -> $CR_BAM_OUT"
}

step_a(){
  [[ -f "$SD/a.done" ]] && { say "a done -> skip"; return 0; }
  [[ -f "$D6/.done" ]] || fail "[a] step 06 not done for $DS ($D6/.done)"
  local ref6; ref6=$(awk -F'\t' -v d="$DS" '$1==d{print $6}' "$DUR/06_gene_matrices/digest_ref.tsv")
  [[ -n "$ref6" ]] || fail "step 06 digest for $DS missing ($DUR/06_gene_matrices/digest_ref.tsv)"
  local refa; refa=$(awk -F'\t' -v d="$DS" '$1==d{print $2}' "$DIGEST_REF_TSV")
  [[ -n "$refa" ]] || fail "no $DS row in $DIGEST_REF_TSV"
  local f
  for f in "$SA_Q/map_info.json" "$SA_Q/sample_info.json" "$SA_Q/quant.json" "$SA_Q/collate.json" "$SA_Q/map.collated.rad" "$UNMAPPED"; do
    [[ -s "$f" ]] || fail "[a] step 06 output $f missing"
  done
  local want; want=$(awk -F'\t' '$1=="af_map/unmapped_bc_count.bin"{print $2}' "$D6/sha256.txt")
  [[ -n "$want" && "$(sha256sum "$UNMAPPED" | cut -c1-64)" == "$want" ]] \
    || fail "[a] $UNMAPPED does not match its sha256 in $D6/sha256.txt"
  # the paper's run: step 06's matrix (= the per-probe reference), with the configured UMI collapse and collate codec
  local dig; dig=$(matrix_digest "$SA_Q/alevin" | cut -f5)
  [[ "$dig" == "$ref6" ]] || fail "[a] step 06 matrix digest $dig != its digest_ref.tsv $ref6"
  [[ "$dig" == "$refa" ]] || fail "[a] step 06 matrix digest $dig != $DIGEST_REF_TSV $refa"
  [[ "$(jget "$SA_Q/collate.json" chunk_codec)" == "$COMPRESS" ]] || fail "[a] collated RAD is not $COMPRESS"
  [[ "$("$PY" -c "import json,sys;print(json.load(open(sys.argv[1]))['quant_options']['crlike_umi_edit'])" "$SA_Q/quant.json")" == "$UMI_EDIT_DIST" ]] \
    || fail "[a] quant.json crlike_umi_edit != $UMI_EDIT_DIST"
  "$PY" "$REPRO/lib/readlevel_ledger.py" ledger "$SA_Q/map_info.json" "$UNMAPPED" "$SA_Q" "$DD/a_ledger.json.part"
  mv "$DD/a_ledger.json.part" "$DD/a_ledger.json"
  mkdir -p "$DD/a_run"
  echo -e "$dig\t$D6" > "$DD/a_run/matrix_digest.txt"
  touch "$SD/a.done"; say "[a] ledger -> $DD/a_ledger.json (step 06 output, matrix digest $dig = $refa)"
}

step_b(){
  [[ -f "$SD/b.done" ]] && { say "b done -> skip"; return 0; }
  [[ -f "$SD/a.done" ]] || fail "[b] needs sub-step a"
  tools
  local i=0 pids=() f
  mkdir -p "$SD/b"
  for f in ${DS_R1//,/ }; do
    i=$((i + 1))
    ( set -o pipefail; pigz -dc "$f" | "$SD/bin/r1_geom_scan" --cb-len "$CB_LEN" --umi-len "$UMI_LEN" \
        --win-lo "$SP_LO" --win-hi "$SP_HI" --anchor "$ANCHOR" --max-dist "$MAXD" --sbc-len "$SBC_LEN" \
        --sbc-list "$DS_SBC" > "$SD/b/part_$i.json" ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p" || fail "[b] an R1 scan failed"; done
  "$PY" "$REPRO/lib/readlevel_ledger.py" reconcile "$DD/a_ledger.json" "$SD"/b/part_*.json "$DD/b_r1_classes.json.part"
  mv "$DD/b_r1_classes.json.part" "$DD/b_r1_classes.json"
  touch "$SD/b.done"; say "[b] R1 classes -> $DD/b_r1_classes.json (reconciled exactly)"
}

step_c(){
  [[ -f "$SD/c.done" ]] && { say "c done -> skip"; return 0; }
  need_bams; tools
  mkdir -p "$SD/c"
  local s i n=0 pids=() names=() p r1 r2
  for s in "${SAMPLES[@]}"; do
    local d="$CR_BAM_OUT/per_sample_outs/$s"
    [[ -s "$d/sample_alignments.bam.bai" ]] || fail "[c] $s: BAM index missing (chunks are coordinate regions)"
    for (( i = 0; i < READLEVEL_CHUNKS; i++ )); do
      local c; c=$(printf '%s.c%02d' "$s" "$i")
      [[ -f "$SD/c/$c.pass.done" ]] && continue
      ( "$PY" "$REPRO/lib/mechanism_readlevel.py" --bam "$d/sample_alignments.bam" --ref "$REF_FA" --k "$KMER" \
          --geom "$DS_GEOM" --sbc "$SBC0" --barcodes "$d/sample_filtered_barcodes.csv" \
          --features "$d/sample_filtered_feature_bc_matrix/features.tsv.gz" --out "$SD/c/$c" \
          --chunk "$i" --nchunks "$READLEVEL_CHUNKS" --gzip-threads 1 --progress 0 \
          > "$LOGD/c_pass_$c.log" 2>&1 && touch "$SD/c/$c.pass.done" ) &
      pids+=($!); names+=("$c")
      if (( ${#pids[@]} >= READLEVEL_JOBS )); then
        wait "${pids[0]}" || fail "[c] BAM pass ${names[0]} failed ($LOGD/c_pass_${names[0]}.log)"
        pids=("${pids[@]:1}"); names=("${names[@]:1}")
      fi
    done
  done
  for i in "${!pids[@]}"; do wait "${pids[$i]}" || fail "[c] BAM pass ${names[$i]} failed ($LOGD/c_pass_${names[$i]}.log)"; done
  for s in "${SAMPLES[@]}"; do
    [[ -f "$SD/c/$s.map.done" ]] && continue
    # chunks with no universe read hold an empty gzip, which piscem cannot open: leave them out
    r1=$("$PY" -c "
import json, sys
d, s, n = sys.argv[1:4]
print(','.join(f'{d}/{s}.c{i:02d}_R1.fastq.gz' for i in range(int(n))
               if json.load(open(f'{d}/{s}.c{i:02d}.json'))['counters'].get('universe', 0)))" "$SD/c" "$s" "$READLEVEL_CHUNKS")
    [[ -n "$r1" ]] || fail "[c] $s: no universe read in any chunk"
    r2=${r1//_R1./_R2.}
    rm -rf "$SD/c/${s}_map"
    taskset -c "$READLEVEL_CPUSET" "$PISCEM" map-sc -i "$DS_IDX/piscem_idx" -g "$DS_GEOM" -o "$SD/c/${s}_map" \
      -t "$READLEVEL_THREADS" --skipping-strategy permissive --max-ec-card 4096 --dict auto --decoder "$DECODER" \
      -1 "$r1" -2 "$r2" > "$LOGD/c_map_$s.log" 2>&1 \
      || { tail -20 "$LOGD/c_map_$s.log"; fail "[c] piscem map-sc $s failed"; }
    "$SD/build/release/rad_dump" map "$SD/c/${s}_map/map.rad" "$SD/c/${s}_pairs" >> "$LOGD/c_map_$s.log" 2>&1 \
      || fail "[c] rad_dump map $s failed"
    touch "$SD/c/$s.map.done"
  done
  p=(); for s in "${SAMPLES[@]}"; do p+=("$SD/c/$s"); done
  "$PY" "$REPRO/lib/readlevel_r2_classify.py" "$DD/c_r2_mapping.json.part" "${p[@]}" > "$LOGD/c_classify.log" 2>&1 \
    || { tail -20 "$LOGD/c_classify.log"; fail "[c] classification failed"; }
  mv "$DD/c_r2_mapping.json.part" "$DD/c_r2_mapping.json"
  rm -f "$SD"/c/*_R[12].fastq.gz "$SD"/c/*_map/map.rad
  touch "$SD/c.done"; say "[c] R2 cross-tab -> $DD/c_r2_mapping.json"
}

step_d(){
  [[ -f "$SD/d.done" ]] && { say "d done -> skip"; return 0; }
  [[ -f "$SD/a.done" && -f "$SD/c.done" ]] || fail "[d] needs sub-steps a and c"
  tools
  mkdir -p "$SD/d"
  local s args=()
  for s in "${SAMPLES[@]}"; do
    cut -d, -f2 "$CR_BAM_OUT/per_sample_outs/$s/sample_filtered_barcodes.csv" | grep -v '^barcode$' | cut -c1-16
  done | sort -u > "$SD/d/cells16_union.txt"
  "$SD/build/release/rad_dump" mol "$SA_Q/map.collated.rad" "$SD/d/cells16_union.txt" "$SD/d/afmol" \
    > "$LOGD/d_rad_dump.log" 2>&1 || { tail -20 "$LOGD/d_rad_dump.log"; fail "[d] rad_dump mol failed"; }
  for s in "${SAMPLES[@]}"; do args+=("$s=$SD/c/$s,$CR_BAM_OUT/per_sample_outs/$s/sample_filtered_barcodes.csv"); done
  "$PY" "$REPRO/lib/readlevel_molclass.py" "$SD/d/afmol" "$SA_Q/sample_info.json" "$DD/d_reassignment.json.part" \
    "${args[@]}" > "$LOGD/d_molclass.log" 2>&1 || { tail -20 "$LOGD/d_molclass.log"; fail "[d] molecule classes failed"; }
  mv "$DD/d_reassignment.json.part" "$DD/d_reassignment.json"
  touch "$SD/d.done"; say "[d] reassignment -> $DD/d_reassignment.json"
}

step_e(){
  [[ -f "$SD/e.done" ]] && { say "e done -> skip"; return 0; }
  local s js=()
  for s in "${SAMPLES[@]}"; do
    local j="$DUR/08_concordance_counts/per_sample/${DS}__${s}.json"
    [[ -f "$j.done" ]] || fail "[e] step 08 output $j missing"
    js+=("$j")
  done
  local tr="$DUR/08_concordance_counts/text_ranges.json"
  [[ -s "$tr" ]] || fail "[e] step 08 $tr missing"
  "$PY" - "$DD/e_umi_totals.json" "$tr" "${js[@]}" <<'PY'
import json, sys
out, tr, files = sys.argv[1], sys.argv[2], sys.argv[3:]
m = json.load(open(tr))["mouse_cellranger_vs_simpleaf"]
if m["n_samples"] != len(files):
    raise SystemExit(f"step 08 mouse means cover {m['n_samples']} samples, expected {len(files)}")
tot = {"cellranger": 0, "simpleaf": 0, "cyto": 0}
cr_for = {"simpleaf": 0, "cyto": 0}
per = {}
for f in files:
    d = json.load(open(f))
    a, b = d["pairs"]["cellranger_vs_simpleaf"], d["pairs"]["cellranger_vs_cyto"]
    for p in (a, b):
        if "error" in p:
            raise SystemExit(f"{f}: {p['error']}")
    tot["cellranger"] += int(a["total_umi_ref"]); tot["simpleaf"] += int(a["total_umi_test"])
    tot["cyto"] += int(b["total_umi_test"])
    cr_for["simpleaf"] += int(a["total_umi_ref"]); cr_for["cyto"] += int(b["total_umi_ref"])
    per[d["sample"]] = {"cr_ref_simpleaf": a["total_umi_ref"], "simpleaf": a["total_umi_test"],
                        "cr_ref_cyto": b["total_umi_ref"], "cyto": b["total_umi_test"],
                        "count_concordance_pct": a["count_concordance_pct"], "gene_spearman": a["gene_spearman"],
                        "n_common_genes_simpleaf": a["n_common_genes"], "n_common_genes_cyto": b["n_common_genes"]}
res = {"totals": tot, "cellranger_total_on_genes_shared_with": cr_for,
       "diff_pct": {t: 100.0 * (tot[t] / cr_for[t] - 1) for t in ("simpleaf", "cyto")},
       "mean_per_sample": {"count_concordance_pct": m["mean_of_per_sample_count_concordance_pct"],
                           "gene_spearman": m["mean_of_per_sample_gene_spearman"], "source": tr},
       "per_sample": per, "sources": files}
json.dump(res, open(out, "w"), indent=2)
print(json.dumps({k: v for k, v in res.items() if k not in ("per_sample", "sources")}))
PY
  touch "$SD/e.done"; say "[e] UMI totals -> $DD/e_umi_totals.json"
}

step_f(){
  [[ -f "$SD/f.done" ]] && { say "f done -> skip"; return 0; }
  need_bams
  local s a b vals=()
  for s in "${SAMPLES[@]}"; do
    a="$CR_BAM_OUT/per_sample_outs/$s"; b="$CR_OUT/$DS/$CR_SUBPATH/$s"
    cmp -s "$a/sample_filtered_barcodes.csv" "$b/sample_filtered_barcodes.csv" \
      || fail "[f] $s: CR BAM run's called cells differ from $CR_OUT's run"
    # the two runs may differ only in CR's sampled gDNA estimate; every other metric and the matrix must match
    diff <(grep -v 'Estimated UMIs from genomic DNA' "$a/metrics_summary.csv") \
         <(grep -v 'Estimated UMIs from genomic DNA' "$b/metrics_summary.csv") > /dev/null \
      || fail "[f] $s: metrics_summary.csv differs between the two CR runs beyond the gDNA estimate"
    "$PY" -c "
import h5py, numpy as np, sys
a, b = (h5py.File(x + '/sample_filtered_feature_bc_matrix.h5', 'r')['matrix'] for x in sys.argv[1:3])
for k in ('data', 'indices', 'indptr', 'barcodes', 'shape', 'features/id'):
    if not np.array_equal(a[k][()], b[k][()]): raise SystemExit('filtered matrix differs in ' + k)
" "$a" "$b" || fail "[f] $s: CR BAM run's filtered matrix differs from $CR_OUT's"
    vals+=("$(awk -F, '$5=="Valid probe barcodes"{print $6}' "$a/metrics_summary.csv")")
  done
  [[ -n "${vals[0]}" ]] || fail "[f] no 'Valid probe barcodes' row"
  for v in "${vals[@]}"; do [[ "$v" == "${vals[0]}" ]] || fail "[f] Valid probe barcodes differs across samples: ${vals[*]}"; done
  printf '{"valid_probe_barcodes": %s, "samples": %s, "cr_bam_out": "%s", "cr_out": "%s"}\n' "${vals[0]}" \
    "$(printf '"%s",' "${SAMPLES[@]}" | sed 's/,$//; s/^/[/; s/$/]/')" "$CR_BAM_OUT" "$CR_OUT/$DS/$CR_SUBPATH" > "$DD/f_cr_probe_bc.json"
  touch "$SD/f.done"; say "[f] Valid probe barcodes ${vals[0]} -> $DD/f_cr_probe_bc.json"
}

step_numbers(){
  for x in a b c d e f; do [[ -f "$SD/$x.done" ]] || fail "[numbers] sub-step $x not done"; done
  "$PY" - "$DD/meta.json" <<PY
import json, sys
json.dump({"dataset": "$DS", "geometry": r"""$DS_GEOM""", "k": $KMER, "m": $MINIMIZER,
           "cell_bc_correction": "$CELL_BC_CORRECTION", "sample_bc_correction": "$SAMPLE_BC_CORRECTION",
           "umi_edit_dist": $UMI_EDIT_DIST, "compress": "$COMPRESS", "index": "$DS_IDX",
           "simpleaf": "$SIMPLEAF_COMMIT", "alevin_fry": "$ALEVINFRY_COMMIT", "piscem": "$PISCEM_COMMIT",
           "piscem_rs": "$PISCEM_RS_VERSION", "cr_bam_out": "$CR_BAM_OUT", "samples": "${SAMPLES[*]}".split()},
          open(sys.argv[1], "w"), indent=2)
PY
  "$PY" "$REPRO/lib/readlevel_numbers.py" "$DD" "$DD/meta.json" "$DD/readlevel_numbers.json" > "$LOGD/numbers.log" 2>&1 \
    || { cat "$LOGD/numbers.log"; fail "[numbers] assembly failed"; }
  cat "$LOGD/numbers.log"
  local LINES=() f
  for f in readlevel_numbers.json a_ledger.json b_r1_classes.json c_r2_mapping.json d_reassignment.json e_umi_totals.json f_cr_probe_bc.json; do
    LINES+=("$f	$(sha256sum "$DD/$f" | cut -c1-64)")
  done
  LINES+=("mouse matrix digest (step 06 = $DIGEST_REF_TSV)	$(cut -f1 "$DD/a_run/matrix_digest.txt")")
  write_done "${LINES[@]}"
}

case "$SUB" in
  bam|a|b|c|d|e|f|numbers) "step_$SUB" ;;
  all) step_bam; step_a; step_b; step_c; step_d; step_e; step_f; step_numbers ;;
  *) fail "unknown --step $SUB" ;;
esac
