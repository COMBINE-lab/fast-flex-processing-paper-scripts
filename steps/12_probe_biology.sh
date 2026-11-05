#!/usr/bin/env bash
# 12: probe annotation, differential probe usage, pre-mRNA ratio and differential expression.
# Produces Figure 2, Supp. Figure 2 and the probe numbers of the manuscript.
#
# Inputs (earlier steps; overridable for tests):
#   PROBE_ROOT=$DUR/07_probe_matrices   <ds>/<sample>.probe.h5ad  (var: gene_id gene_name region_csv
#                                        genomic_class n_tx_hit n_tx_gene constitutive tx_set)
#   LABEL_ROOT=$DUR/11_celltypes        <ds>/<sample>.labels.tsv.gz (barcode index, column `label`)
# Every output is under $DUR/12_probe_biology/ (the 4 GB-class embedding is built in $SCR first).
# Which annotation each analysis uses:
#   region_csv     the 10x CSV `region` column (spliced / unspliced), vendor-provided
#                  -> sub-step d: unspliced/spliced ratio, Figure 2 panel c, Figure 2 panel c
#   genomic_class  the 50-mer relocation against GRCh38/GRCm39-2024-A (lib/probe_annotate.py:
#                  within_exon / junction / intronic / not_found) and tx_set (isoform set)
#                  -> sub-step b: premrna instrument (tx_set-identical junction vs within-exon
#                  pairs), the junction/within-exon panels and table; dpu / genede verdicts;
#                  the isoform reading of the YBX3 / ZBTB16 / PHACTR1 probes
#   Both human libraries are produced (results/<ds>/premrna*.tsv, figures/<ds>/); config.env names
#   which one each manuscript item uses (PREMRNA_TABLE_DS, PREMRNA_FIG_HUMAN_DS).
#
# Sub-steps (each resumable through $DD/.done_<sub>_<ds>; run all with no argument):
#   --step a  probe annotation: relocate every probe of each panel (lib/probe_annotate.py),
#             summarise the Methods counts, gate that the step-07 matrices carry it
#   --step b  pseudobulk (label x sample) + premrna for every library; the junction/within-exon table for
#             PREMRNA_TABLE_DS; per-label dpu + genede for STORY_DS only
#   --step c  paired PyDESeq2 gene-vs-probe DE with the minimum-UMI floor for STORY_LABEL of STORY_DS
#             only (Figure 2 panels a and b), + story numbers
#   --step d  unspliced/spliced ratios from the 10x region column, RATIO_DS only
#   --step e  pooled embedding of ALL donors of $STORY_DS ($P12_EMBED_THREADS threads; ~30-60 min,
#             run beside a-d in `all` mode)
# Concurrency (untimed step): every sub-step builds a flat job list over (library, label) and runs
# it P12_JOBS at a time, each job with P12_THREADS BLAS / PyDESeq2 threads (config.env). Heavy
# intermediates (embedding build) are on $SCR; every deliverable is on $DUR.
#   --step f  junction/within-exon panels (RATIO_DS, PREMRNA_FIG_HUMAN_DS) and the DE-gene verdicts (STORY_DS)
#   --step g  Figure 2 (cache + layout) and Supp. Figure 2 on the step-e embedding
set -euo pipefail
source "$(dirname "$0")/../config.env"
PROBE_ROOT="${PROBE_ROOT:-$DUR/07_probe_matrices}"
LABEL_ROOT="${LABEL_ROOT:-$DUR/11_celltypes}"
# test inputs must never write into the real durable tree: overridden inputs need an overridden DUR
REAL_DUR="$(env -u DUR -u SCR bash -c 'source "$1/config.env"; echo "$DUR"' _ "$REPRO")"
if [[ "$PROBE_ROOT" != "$DUR/07_probe_matrices" || "$LABEL_ROOT" != "$DUR/11_celltypes" ]] && \
   [[ "$(readlink -m "$DUR")" == "$(readlink -m "$REAL_DUR")" ]]; then
  echo "[12_probe_biology] FATAL: PROBE_ROOT/LABEL_ROOT are overridden but DUR is the real $REAL_DUR; set DUR and SCR to a test dir" >&2; exit 1
fi
source "$REPRO/lib/common.sh"   # after the guard: it creates $DUR/<step>

DATASETS_P12="${P12_DATASETS:-$P12_DATASETS_ALL}"
read -r -a SAMPLES_ARG <<< "${P12_SAMPLES:-}"          # tests: restrict every dataset to these
SAMPLES_OPT=(); [[ ${#SAMPLES_ARG[@]} -gt 0 ]] && SAMPLES_OPT=(--samples "${SAMPLES_ARG[@]}")
# Figure 2 donors (config FIG_UMAP_DONORS = all | comma list): the SAME donors feed panel a and
# the panel-b embedding (sub-steps e and g); recorded in the e marker and numbers/umap_<ds>.json
FIG_DONORS="${FIG_UMAP_DONORS:?set FIG_UMAP_DONORS in config.env}"
# P12_REMOVE_DOUBLETS=1: scDblFinder doublets (step-11 label files) leave every pseudobulk and every
# cell-level analysis (b c d e g); the counts removed are written next to each output
case "${P12_REMOVE_DOUBLETS:?set P12_REMOVE_DOUBLETS in config.env}" in
  1) DBL_OPT=(--drop-doublets);; 0) DBL_OPT=();; *) fail "P12_REMOVE_DOUBLETS must be 0 or 1";; esac
if [[ "$FIG_DONORS" == "all" ]]; then FIG_OPT=("${SAMPLES_OPT[@]}")
else IFS=',' read -r -a _fd <<< "$FIG_DONORS"; FIG_OPT=(--samples "${_fd[@]}"); fi
read -r -a STORY_GENES_A <<< "$STORY_GENES"
read -r -a STORY_TX_A <<< "$STORY_TX"
read -r -a STORY_PANEL_A_GENES <<< "$STORY_PANEL_A"
SELF="$(readlink -f "$0")"
SUB="all"; JOB=""
[[ "${1:-}" == "--step" ]] && SUB="${2:?--step needs a|b|c|d|e|f|g|s}"
[[ "${1:-}" == "--job" ]] && { JOB="${2:?}"; shift 2; JOBARGS=("$@"); }
# --check (run_all's validity hook): read-only; one "<sub-step marker>\tdone|redo" line per marker the
# DONE condition needs, exit 0 iff all are done under the current settings and DONE exists
CHECK=0; [[ "${1:-}" == "--check" ]] && CHECK=1

export NUMBA_CACHE_DIR="$SD/numba" MPLCONFIGDIR="$SD/mpl" PYTHONWARNINGS=ignore
# per-job BLAS threads; P12_JOBS jobs run at once (P12_JOBS x P12_THREADS = the host)
export OMP_NUM_THREADS=$P12_THREADS OPENBLAS_NUM_THREADS=$P12_THREADS MKL_NUM_THREADS=$P12_THREADS
(( CHECK )) || mkdir -p "$SD" "$DD"/{annot,pb,results,numbers,figures,embedding}
LIBD="$REPRO/lib"

safe(){ printf '%s' "$1" | tr -c 'A-Za-z0-9+._-' '_'; }
# job NAME ARGS... -> one line that re-invokes this script for one unit of work (--job mode)
job(){ printf '%q ' bash "$SELF" --job "$@"; echo; }
# run_jobs N FILE: run FILE's lines N at a time; any failing line fails the step (xargs exit 123)
run_jobs(){
  local n="$1" f="$2"
  [[ -s "$f" ]] || return 0
  say "running $(wc -l < "$f") jobs, $n at a time"
  xargs -P "$n" -d '\n' -I{} bash -c '{}' < "$f" || fail "a job of $f failed (see $LOGD)"
}
done_mark(){ echo "$DD/.done_$1"; }
# markers of cell-level sub-steps record P12_REMOVE_DOUBLETS; a marker written under the other
# setting is not done, so a changed setting redoes b-g but not the
# doublet-independent relocation (a)
cellstep(){ [[ "$1" != a_* && "$1" != annot_* ]]; }
# the settings a marker must carry: P12_REMOVE_DOUBLETS for b-g; the panel-a examples for c and g
tags(){
  cellstep "$1" || return 0
  echo "P12_REMOVE_DOUBLETS=$P12_REMOVE_DOUBLETS"
  case "$1" in c_*|g_*) echo "FIG2_PANEL_A=$FIG2_PANEL_A_DILUTED;$FIG2_PANEL_A_CANCELLING";;
               s) echo "FIG2_SEARCH_DATASETS=$FIG2_SEARCH_DATASETS";; esac
}
is_done(){
  local f t; f=$(done_mark "$1"); [[ -f "$f" ]] || return 1
  while IFS= read -r t; do [[ -z "$t" ]] || grep -qxF "$t" "$f" || return 1; done < <(tags "$1")
}
# d writes outside results/<ds>, because b_prep wipes results/<ds> (a, d, b, ... order)
ratio_tsv(){ echo "${DD:?}/ratio/${1:?}/allgene_ratio.tsv"; }
mark(){ { date -u +%FT%TZ; tags "$1"; } > "$(done_mark "$1")"; }
probe_csv_of(){ local f; f=$(ls "${P12_REFS:-$REFS_SRC}/$1"/*Transcriptome_Probe_Set*.csv); [[ $(wc -l <<< "$f") -eq 1 ]] || fail "$1: not exactly one probe-set CSV"; echo "$f"; }
genome_of(){ case "$(basename "$1")" in *Mouse*) echo "$GENOME_MOUSE";; *Human*) echo "$GENOME_HUMAN";; *) fail "species of $1?";; esac; }
gtf_of(){ local g; g=$(ls "$1"/genes/genes.gtf*); [[ $(wc -l <<< "$g") -eq 1 ]] || fail "$1: not exactly one GTF"; echo "$g"; }
first_h5ad(){ local f; f=$(ls "$PROBE_ROOT/$1"/*.probe.h5ad | head -1); [[ -n "$f" ]] || fail "$1: no probe matrices under $PROBE_ROOT"; echo "$f"; }

# ---------------------------------------------------------------- a: annotation
a_panel(){   # one relocation per panel, shared by the libraries that use it
  local csv="$1" panel genome out
  panel=$(basename "$csv" .csv); out="$DD/annot/$panel.tsv"; genome=$(genome_of "$csv")
  say "a: relocating $panel against $genome"
  local tmp; tmp=$(mktemp "$out.XXXXXX.tmp")       # job-unique, same dir: the mv is atomic
  "$PY_PROBE" "$LIBD/probe_annotate.py" --probe-csv "$csv" --gtf "$(gtf_of "$genome")" \
    --fasta "$genome/fasta/genome.fa" --out "$tmp" > "$LOGD/a_annot_$panel.log" 2>&1
  mv "$tmp" "$out"; mark "annot_$panel"
}
a_summary(){
  local ds="$1" csv; csv=$(probe_csv_of "$ds")
  "$PY_PROBE" "$LIBD/annot_summary.py" --probe-csv "$csv" --annot "$DD/annot/$(basename "$csv" .csv).tsv" \
    --out "$DD/numbers/annot_$ds.json" --h5ad "$(first_h5ad "$ds")" > "$LOGD/a_summary_$ds.log" 2>&1
  mark "a_$ds"
}
step_a(){
  # one relocation job per UNIQUE panel: libraries sharing a panel (16plex_human and 320k_GEMX
  # both use human v1.1.0) each have their own copy of the CSV, so dedupe by panel name, not path
  local J="$SD/jobs_a1" J2="$SD/jobs_a2"; : > "$J"; : > "$J2"
  local -A seen=()
  for ds in $DATASETS_P12; do
    is_done "a_$ds" && { say "a $ds: done"; continue; }
    local csv panel; csv=$(probe_csv_of "$ds"); panel=$(basename "$csv" .csv)
    if ! is_done "annot_$panel" && [[ -z "${seen[$panel]:-}" ]]; then
      seen[$panel]="$ds"; job a_panel "$csv" >> "$J"
    fi
    job a_summary "$ds" >> "$J2"
  done
  run_jobs "$P12_JOBS" "$J"; run_jobs "$P12_JOBS" "$J2"
}

# ---------------------------------------------------------------- b: pseudobulk + instruments
b_prep(){
  local ds="$1" L="$LOGD/b_$1.log"
  local pb="${DD:?}/pb/${ds:?}" res="${DD:?}/results/${ds:?}"
  rm -rf "$pb" "$res"; mkdir -p "$res"            # a partial run is redone, never resumed
  say "b $ds: pseudobulk"
  { "$PY_PROBE" "$LIBD/probe_instruments.py" pseudobulk --probe-root "$PROBE_ROOT" \
      --label-root "$LABEL_ROOT" --dataset "$ds" --outdir "$pb" "${SAMPLES_OPT[@]}" "${DBL_OPT[@]}"
    "$PY_PROBE" "$LIBD/probe_instruments.py" premrna --pb "$pb" --out "$res/premrna.tsv"
    if [[ "$ds" == "$PREMRNA_TABLE_DS" ]]; then
      "$PY_PROBE" "$LIBD/probe_instruments.py" premrna_table --premrna "$res/premrna.tsv" --pb "$pb" \
        --out "$res/premrna_table.tsv"
    fi
    if [[ "$ds" == "$STORY_DS" ]]; then      # verdicts: every qualifying label of this library
      "$PY_PROBE" "$LIBD/probe_instruments.py" targets --pb "$pb" --min-cells "$DE_MIN_CELLS" \
        --min-samples "$DE_MIN_SAMPLES" --out "$res/targets.txt"
    fi; } > "$L" 2>&1
  [[ "$ds" == "$STORY_DS" ]] || return 0
  [[ -s "$res/targets.txt" ]] || fail "b $ds: no label qualifies for the per-label tests"
  if [[ -n "${P12_TARGETS:-}" ]]; then                # tests: restrict to these labels ('|'-separated)
    local want; IFS='|' read -r -a want <<< "$P12_TARGETS"
    printf '%s\n' "${want[@]}" | grep -Fxf "$res/targets.txt" > "$res/targets.sub" || fail "b $ds: none of P12_TARGETS qualifies"
    mv "$res/targets.sub" "$res/targets.txt"
  fi
}
b_label(){
  local ds="$1" lb="$2" pb="$DD/pb/$1" res="$DD/results/$1" s; s=$(safe "$2")
  say "b $ds: dpu + genede $lb"
  { "$PY_PROBE" "$LIBD/probe_instruments.py" dpu --pb "$pb" --target "$lb" --out "$res/dpu__${s}__vs_rest.tsv"
    "$PY_PROBE" "$LIBD/probe_instruments.py" genede --pb "$pb" --target "$lb" \
      --dpu "$res/dpu__${s}__vs_rest.tsv" --out "$res/genede_join__${s}.tsv"; } > "$LOGD/b_${ds}__$s.log" 2>&1
}
step_b(){
  local todo=() J="$SD/jobs_b1" J2="$SD/jobs_b2"; : > "$J"; : > "$J2"
  for ds in $DATASETS_P12; do
    if is_done "b_$ds"; then say "b $ds: done"; else todo+=("$ds"); job b_prep "$ds" >> "$J"; fi
  done
  run_jobs "$P12_JOBS" "$J"                         # every library's pseudobulk at once
  local ds="$STORY_DS"                              # dpu + genede feed only the verdict figure (STORY_DS)
  if [[ " ${todo[*]} " == *" $ds "* ]]; then
    while IFS= read -r lb; do job b_label "$ds" "$lb" >> "$J2"; done < "$DD/results/$ds/targets.txt"
    run_jobs "$P12_JOBS" "$J2"                      # every label's test at once
    [[ $(ls "$DD/results/$ds"/genede_join__*.tsv | wc -l) -eq $(wc -l < "$DD/results/$ds/targets.txt") ]] || fail "b $ds: genede tables missing"
  fi
  for ds in "${todo[@]}"; do mark "b_$ds"; done
}

# ---------------------------------------------------------------- c: masked DE + story numbers
c_label(){
  local ds="$1" lb="$2" pb="$DD/pb/$1" out="$DD/results/$1/masked_de" s n need cnt; s=$(safe "$2")
  cnt=$("$PY_PROBE" "$LIBD/pseudobulk_de.py" --pb "$pb" --target "$lb" --count-donors)
  read -r n need <<< "$cnt"
  if [[ "$n" -lt "$need" ]]; then
    # recorded, not silent: the paired design needs >= MIN_DONORS donors with enough cells
    # (one short line per append: concurrent O_APPEND writes do not interleave)
    printf '%s\t%s\tskipped: %s paired donors < %s\n' "$ds" "$lb" "$n" "$need" >> "$out/skipped.tsv"
    return 0
  fi
  say "c $ds: PyDESeq2 $lb ($n donors)"
  "$PY_PROBE" "$LIBD/pseudobulk_de.py" --pb "$pb" --target "$lb" --out "$out/${s}.tsv" \
    --min-probe-umi "$MIN_PROBE_UMI_FLOOR" --n-cpus "$P12_THREADS" > "$LOGD/c_${ds}__$s.log" 2>&1
}
step_c(){
  # the paired gene-vs-probe DESeq2 feeds Figure 2 panel b (STORY_LABEL of STORY_DS: PHACTR1) and
  # Figure 2 panel a (the FIG2_PANEL_A_DILUTED and _CANCELLING examples)
  # (the example libraries are part of P12_DATASETS_ALL by construction, config.env)
  local ex role ds lb gene
  for ds in $DATASETS_P12; do
    local jobs=() out="$DD/results/$ds/masked_de"
    [[ "$ds" == "$STORY_DS" ]] && jobs+=("story|$STORY_LABEL|$STORY_GENES")
    for ex in "diluted|$FIG2_PANEL_A_DILUTED" "cancelling|$FIG2_PANEL_A_CANCELLING"; do
      IFS='|' read -r role d2 lb gene <<< "$ex"; [[ "$d2" == "$ds" ]] && jobs+=("$role|$lb|$gene")
    done
    [[ ${#jobs[@]} -gt 0 ]] || continue
    is_done "c_$ds" && { say "c $ds: done"; continue; }
    is_done "b_$ds" || fail "c $ds: run --step b first"
    rm -rf "${DD:?}/results/${ds:?}/masked_de"; mkdir -p "$out"; : > "$out/skipped.tsv"
    local J="$SD/jobs_c_$ds"; : > "$J"
    printf '%s\n' "${jobs[@]}" | cut -d'|' -f2 | sort -u | while IFS= read -r lb; do job c_label "$ds" "$lb" >> "$J"; done
    run_jobs "$P12_JOBS" "$J"
    for j in "${jobs[@]}"; do
      IFS='|' read -r role lb gene <<< "$j"
      [[ -s "$out/$(safe "$lb").tsv" ]] || fail "c $ds: no DE table for $lb ($(cat "$out/skipped.tsv"))"
      local js="$DD/numbers/story_$ds.json"; [[ "$role" == story ]] || js="$DD/numbers/panelA_$role.json"
      "$PY_PROBE" "$LIBD/story_numbers.py" --masked "$out/$(safe "$lb").tsv" --pb "$DD/pb/$ds" \
        --label "$lb" --genes $gene --tx "${STORY_TX_A[@]}" \
        --min-probe-umi "$MIN_PROBE_UMI_FLOOR" --out "$js" > "$LOGD/c_${role}_$ds.log" 2>&1
    done
    mark "c_$ds"
  done
}

# ---------------------------------------------------------------- d: 10x region ratios
d_one(){
  local ds="$1" order=(); IFS='|' read -r -a order <<< "$RATIO_ORDER"
  mkdir -p "$(dirname "$(ratio_tsv "$ds")")"
  "$PY_PROBE" "$LIBD/allgene_ratio.py" --probe-root "$PROBE_ROOT" --label-root "$LABEL_ROOT" \
    --dataset "$ds" "${SAMPLES_OPT[@]}" "${DBL_OPT[@]}" --out "$(ratio_tsv "$ds")" \
    --json "$DD/numbers/ratio_$ds.json" --order "${order[@]}" > "$LOGD/d_$ds.log" 2>&1
  [[ -s "$(ratio_tsv "$ds")" ]] || fail "d $ds: $(ratio_tsv "$ds") not written"
  mark "d_$ds"
}
step_d(){
  # the 10x-region ratio feeds Figure 2 panel c and Figure 2 panel c: RATIO_DS only
  local J="$SD/jobs_d"; : > "$J"
  for ds in $RATIO_DS; do
    [[ " $DATASETS_P12 " == *" $ds "* ]] || continue
    if is_done "d_$ds" && [[ -s "$(ratio_tsv "$ds")" ]]; then say "d $ds: done"; else job d_one "$ds" >> "$J"; fi
  done
  run_jobs "$P12_JOBS" "$J"
}

# ---------------------------------------------------------------- e: embedding
step_e(){
  local ds="$STORY_DS"
  [[ " $DATASETS_P12 " == *" $ds "* ]] || return 0
  is_done "e_$ds" && { say "e $ds: done"; return 0; }
  mkdir -p "$SD/embedding"
  say "e $ds: pooled embedding of the labelled cells of FIG_UMAP_DONORS=$FIG_DONORS ${SAMPLES_ARG[*]:-}"
  OMP_NUM_THREADS=$P12_EMBED_THREADS NUMBA_NUM_THREADS=$P12_EMBED_THREADS OPENBLAS_NUM_THREADS=$P12_EMBED_THREADS \
  MKL_NUM_THREADS=$P12_EMBED_THREADS "$PY_PROBE" "$LIBD/umap_probe_matrix.py" --probe-root "$PROBE_ROOT" --label-root "$LABEL_ROOT" \
    --dataset "$ds" "${FIG_OPT[@]}" "${DBL_OPT[@]}" --out "$SD/embedding/umap_pooled_$ds.h5ad" \
    --n-neighbors "$UMAP_N_NEIGHBORS" --seed "$UMAP_SEED" 2>&1 | tee "$LOGD/e_$ds.log"
  cp "$SD/embedding/umap_pooled_$ds".* "$DD/embedding/"   # h5ad (+ doublets_removed.tsv)
  mark "e_$ds"; echo "FIG_UMAP_DONORS=$FIG_DONORS" >> "$(done_mark "e_$ds")"
}

# ---------------------------------------------------------------- s: panel-a candidate search
# every library of FIG2_SEARCH_DATASETS, every cell type with >= 3 paired donors: pseudobulk (same
# doublet setting) + gene-vs-probe DESeq2, then lib/fig2_simple_examples.py ranks CellTypist
# markers that are flat at gene level and have a DE probe (cancelling / diluted); the user picks
# FIG2_PANEL_A_* from search/flat_markers_with_DE_probes.tsv
s_pb(){
  local ds="$1" pb="${DD:?}/search/${1:?}/pb"
  rm -rf "$pb"
  "$PY_PROBE" "$LIBD/probe_instruments.py" pseudobulk --probe-root "$PROBE_ROOT" --label-root "$LABEL_ROOT" \
    --dataset "$ds" --outdir "$pb" "${SAMPLES_OPT[@]}" "${DBL_OPT[@]}" > "$LOGD/s_pb_$ds.log" 2>&1
  "$PY_PROBE" - "$pb" "$LIBD" > "$DD/search/$ds/labels.txt" <<'PYEOF'
import sys, pandas as pd
sys.path.insert(0, sys.argv[2])
import pseudobulk_de as p
u = pd.read_csv(sys.argv[1] + "/units.tsv", sep="\t")
for lb in sorted(u.label.unique()):
    if p.n_paired_donors(sys.argv[1], lb) >= p.MIN_DONORS: print(lb)
PYEOF
}
s_de(){
  local ds="$1" lb="$2" s; s=$(safe "$2")
  "$PY_PROBE" "$LIBD/pseudobulk_de.py" --pb "$DD/search/$ds/pb" --target "$lb" --out "$DD/search/$ds/de/$s.tsv" \
    --min-probe-umi "$MIN_PROBE_UMI_FLOOR" --n-cpus "$P12_THREADS" > "$LOGD/s_de_${ds}__$s.log" 2>&1
}
step_s(){
  is_done "s" && { say "s: done"; return 0; }
  local J="$SD/jobs_s1" J2="$SD/jobs_s2" ds; : > "$J"; : > "$J2"
  for ds in $FIG2_SEARCH_DATASETS; do mkdir -p "$DD/search/$ds/de"; job s_pb "$ds" >> "$J"; done
  run_jobs "$P12_JOBS" "$J"
  for ds in $FIG2_SEARCH_DATASETS; do
    [[ -s "$DD/search/$ds/labels.txt" ]] || fail "s $ds: no cell type with >= 3 paired donors"
    while IFS= read -r lb; do job s_de "$ds" "$lb" >> "$J2"; done < "$DD/search/$ds/labels.txt"
  done
  run_jobs "$P12_JOBS" "$J2"
  local DEARGS=()
  for ds in $FIG2_SEARCH_DATASETS; do
    [[ $(ls "$DD/search/$ds/de"/*.tsv | wc -l) -eq $(wc -l < "$DD/search/$ds/labels.txt") ]] || fail "s $ds: DE tables missing"
    DEARGS+=("$ds=$DD/search/$ds/de")
  done
  "$EXT_PY" "$LIBD/fig2_simple_examples.py" --model "$DUR/11_celltypes/celltypist_home/data/models/$CELLTYPIST_MODEL" \
    --out "$DD/search/flat_markers_with_DE_probes.tsv" --de "${DEARGS[@]}" --max-rank "$FIG2_SEARCH_MAX_RANK" \
    > "$LOGD/s_search.log" 2>&1
  [[ -s "$DD/search/flat_markers_with_DE_probes.tsv" ]] || fail "s: empty candidate table"
  mark "s"
}

# ---------------------------------------------------------------- f: probe figures
f_one(){   # f_one DS WHICH...  (premrna | verdicts)
  local ds="$1"; shift
  "$PY_PROBE" "$LIBD/probe_figures.py" --dataset "$ds" --results "$DD/results/$ds" \
    --outdir "$DD/figures/$ds" --which "$@" > "$LOGD/f_probe_figures_$ds.log" 2>&1
  mark "f_$ds"
}
step_f(){
  # junction/within-exon = RATIO_DS (mouse) + PREMRNA_FIG_HUMAN_DS panels; verdicts = STORY_DS
  local J="$SD/jobs_f"; : > "$J"
  for ds in $(f_datasets); do
    [[ " $DATASETS_P12 " == *" $ds "* ]] || continue
    is_done "f_$ds" && { say "f $ds: done"; continue; }
    is_done "b_$ds" || fail "f: run --step b for $ds first"
    local w=()
    [[ "$ds" == "$RATIO_DS" || "$ds" == "$PREMRNA_FIG_HUMAN_DS" ]] && w+=(premrna)
    [[ "$ds" == "$STORY_DS" ]] && w+=(verdicts)
    job f_one "$ds" "${w[@]}" >> "$J"
  done
  run_jobs "$P12_JOBS" "$J"
}
opt(){ [[ " $P12_OPTIONAL " == *" $1 "* ]]; }   # --full only: d, the extra f panels, s
f_datasets(){ { echo "$STORY_DS"; if opt f_extra; then printf '%s\n' "$RATIO_DS" "$PREMRNA_FIG_HUMAN_DS"; fi; } | sort -u; }

# ---------------------------------------------------------------- g: Figure 2, Supp. Figure 2
step_g(){
  local ds="$STORY_DS" F="$DD/figures/$STORY_DS"
  [[ " $DATASETS_P12 " == *" $ds "* ]] || return 0
  is_done "g_$ds" && { say "g $ds: done"; return 0; }
  is_done "e_$ds" || fail "g: run --step e first"
  for r in diluted cancelling; do [[ -s "$DD/numbers/panelA_$r.json" ]] || fail "g: run --step c first (no panelA_$r.json)"; done
  grep -qx "FIG_UMAP_DONORS=$FIG_DONORS" "$(done_mark "e_$ds")" || \
    fail "g: the embedding was built for another FIG_UMAP_DONORS ($(grep FIG_UMAP_DONORS "$(done_mark "e_$ds")")); redo --step e"
  is_done "b_$ds" && [[ -s "$DD/results/$ds/premrna.tsv" ]] || fail "g: run --step b for $ds first (panel c reads its premrna.tsv)"
  mkdir -p "$F"
  local emb="$DD/embedding/umap_pooled_$ds.h5ad"
  "$PY_PROBE" "$LIBD/umap_probe_plot.py" --probe-root "$PROBE_ROOT" --label-root "$LABEL_ROOT" \
    --dataset "$ds" "${FIG_OPT[@]}" "${DBL_OPT[@]}" --embedding "$emb" --out "$F/fig_umap_cache.pkl" \
    --label "$STORY_LABEL" --panel-a "$FIG2_PANEL_A_DILUTED" "$FIG2_PANEL_A_CANCELLING" --gene "$STORY_PANEL_B" \
    --panel-a-json "diluted=$DD/numbers/panelA_diluted.json" "cancelling=$DD/numbers/panelA_cancelling.json" \
    --tx "${STORY_TX_A[@]}" --k-pool "$UMAP_K_POOL" --n-jobs "$P12_EMBED_THREADS" 2>&1 | tee "$LOGD/g_umap_cache.log"
  cp "$F/fig_umap_cache.summary.json" "$DD/numbers/umap_$ds.json"
  "$PY_PROBE" "$LIBD/fig2_plot.py" --cache "$F/fig_umap_cache.pkl" \
    --premrna "$DD/results/$ds/premrna.tsv" --out "$F/fig_umap_probe.png" 2>&1 | tee "$LOGD/g_fig_umap.log"
  "$PY_PROBE" "$LIBD/fig_umap_labels.py" --embedding "$emb" --highlight "$STORY_LABEL" \
    --out "$F/fig_umap_labels.png" 2>&1 | tee "$LOGD/g_fig_umap_labels.log"
  mark "g_$ds"
}

if (( CHECK )); then
  ok=1
  need=(); for ds in $DATASETS_P12; do need+=("a_$ds" "b_$ds"); done
  for ds in $(f_datasets); do need+=("f_$ds"); done
  need+=("c_$STORY_DS" "e_$STORY_DS" "g_$STORY_DS"); if opt d; then need+=("d_$RATIO_DS"); fi; if opt s; then need+=("s"); fi
  for m in "${need[@]}"; do
    if is_done "$m"; then printf '%s\tdone\n' "$m"; else printf '%s\tredo\n' "$m"; ok=0; fi
  done
  if [[ -f "$DD/DONE" ]]; then printf 'DONE\tdone\n'; else printf 'DONE\tredo\n'; ok=0; fi
  (( ok )) && exit 0 || exit 1
fi
if [[ -n "$JOB" ]]; then                    # one unit of work, called by run_jobs
  case "$JOB" in a_panel|a_summary|b_prep|b_label|c_label|d_one|f_one|s_pb|s_de) "$JOB" "${JOBARGS[@]}";;
    *) fail "unknown job $JOB";; esac
  exit 0
fi
case "$SUB" in
  a) step_a;; b) step_b;; c) step_c;; d) step_d;; e) step_e;; f) step_f;; g) step_g;; s) step_s;;
  all)
    # the embedding (one $P12_EMBED_THREADS-thread process) runs beside the per-label jobs
    step_e > "$LOGD/all_e.log" 2>&1 & EPID=$!
    step_a; if opt d; then step_d; fi; step_b; step_c; step_f; if opt s; then step_s; fi
    wait "$EPID" || fail "sub-step e failed (see $LOGD/all_e.log)"
    step_g;;
  *) fail "unknown sub-step $SUB";;
esac

# DONE only when every sub-step of every dataset is complete
for ds in $DATASETS_P12; do for s in a b; do is_done "${s}_$ds" || exit 0; done; done
for ds in $(f_datasets); do is_done "f_$ds" || exit 0; done
is_done "c_$STORY_DS" && is_done "e_$STORY_DS" && is_done "g_$STORY_DS" || exit 0
if opt d; then is_done "d_$RATIO_DS" || exit 0; fi
if opt s; then is_done s || exit 0; fi
"$PY_PROBE" - "$DD" <<'PYEOF' > "$DD/numbers/probe_biology.json"
import json, glob, os, sys
dd = sys.argv[1]
print(json.dumps({os.path.basename(f)[:-5]: json.load(open(f))
                  for f in sorted(glob.glob(f"{dd}/numbers/*.json"))
                  if not f.endswith("probe_biology.json")}, indent=1))
PYEOF
mapfile -t DELIV < <(cd "$DD" && find numbers figures results ratio embedding pb -type f \( -name '*doublets_removed.tsv' -o -name '*.h5ad' -o -name '*.json' -o -name '*.png' -o -name 'premrna_table.tsv' -o -name 'allgene_ratio.tsv' \) | sort)
LINES=(); for f in "${DELIV[@]}"; do LINES+=("$f	$(sha256sum "$DD/$f" | cut -d' ' -f1)"); done
write_done "${LINES[@]}"
