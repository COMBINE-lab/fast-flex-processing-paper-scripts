#!/usr/bin/env bash
# 14: every manuscript figure, regenerated from the durable outputs.
#
# Outputs ($DUR/14_figures/; the dispatcher's name, not the slug), laid out as the manuscript
# references them:
#   plots/fig1_three_panel.png             make_fig1.py on step 05 results.tsv + step 05b results.tsv
#   plots/fig1_values.tsv, fig1_axis.json  the values behind the bars; the runtime axis break (caption)
#   figures/320k_GEMX/fig_umap_probe.png, fig_umap_labels.png   (step 12)
#   figures/<UMAP_DATASET>/fig_supp1_pooled_umap_labels.png  Supp. Figure 1: step 10's Cell Ranger pooled
#                                           embedding, cells colored by step 11's per-tool labels (tool_labels/)
#   supp1_values.json                      its samples, cells and per-tool agreement with Cell Ranger
# Parts (P14_PARTS, default all): fig1 probe pooled.
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"
# the dispatcher named the output $DUR/14_figures (not the slug-derived $DUR/14_figures_numbers)
find "$DD" -depth -type d -empty -delete     # the empty dirs common.sh made
DD="$DUR/14_figures"; LOGD="$DD/logs"; mkdir -p "$DD/plots" "$DD/figures" "$LOGD"
PARTS="${P14_PARTS:-fig1 probe pooled}"
LIBD="$REPRO/lib"; P12="$DUR/12_probe_biology/figures"
export MPLCONFIGDIR="$SD/mpl" PYTHONWARNINGS=ignore
mkdir -p "$SD"
want(){ [[ " $PARTS " == *" $1 "* ]]; }
need(){ [[ -s "$1" ]] || fail "missing input $1 (has its step run?)"; }
put(){ need "$1"; mkdir -p "$(dirname "$DD/$2")"; cp "$1" "$DD/$2"; say "  $2 <- $1"; }

if want fig1; then
  need "$DUR/05_timing_simpleaf/results.tsv"; need "$DUR/05b_timing_cyto_cellranger/results.tsv"
  "$PY_PROBE" "$LIBD/make_fig1.py" --simpleaf "$DUR/05_timing_simpleaf/results.tsv" \
    --others "$DUR/05b_timing_cyto_cellranger/results.tsv" \
    --out "$DD/plots/fig1_three_panel.png" --values-out "$DD/plots/fig1_values.tsv" \
    --axis-out "$DD/fig1_axis.json" 2>&1 | tee "$LOGD/fig1.log"
fi

if want probe; then
  for f in fig_umap_probe.png fig_umap_labels.png; do
    put "$P12/$STORY_DS/$f" "figures/$STORY_DS/$f"
  done
fi

if want pooled; then      # Supp. Figure 1: Cell Ranger's pooled embedding, colored by each tool's labels
  E="$DUR/10_pooled_umap/${UMAP_DATASET}_cellranger_book.obs.tsv.gz"; T="$DUR/11_celltypes/tool_labels/$UMAP_DATASET.tsv.gz"
  need "$E"; need "$T"; mkdir -p "$DD/figures/$UMAP_DATASET"
  "$PY_PROBE" "$LIBD/supp_fig1_pooled_umap.py" --embedding "$E" --labels "$T" \
    --out "$DD/figures/$UMAP_DATASET/fig_supp1_pooled_umap_labels.png" --values-out "$DD/supp1_values.json" \
    2>&1 | tee "$LOGD/supp1.log"
fi


# DONE only for a complete run
if [[ "$PARTS" == "fig1 probe pooled" ]]; then
  mapfile -t DELIV < <(cd "$DD" && find plots figures supp1_values.json fig1_axis.json -type f | sort)
  LINES=(); for f in "${DELIV[@]}"; do LINES+=("$f	$(sha256sum "$DD/$f" | cut -d' ' -f1)"); done
  write_done "${LINES[@]}"
fi
