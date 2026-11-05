# How the results in the manuscript were generated

This repository holds the code that produced every number, table and figure of the manuscript.
It is a record of what was run, not a turnkey package: the sequencing data, the 10x reference
files and the software environment are not included and are yours to provide.

`config.env` holds every path, version, flag and threshold. `REPRO` is the directory of
`config.env`.

## Layout

| path | content |
| --- | --- |
| `config.env` | every path, version, flag and threshold |
| `steps/` | the numbered steps and `run_all.sh` |
| `lib/` | the scripts the steps call |
| `data/` | the dataset table, the Cell Ranger multi templates and the cyto lockfile |

Host paths go in `site.env` (copy `site.env.example`): `DATA_ROOT` holds the inputs,
`SCR` is local scratch, `DUR` is where every output is kept.

## Input files you provide

Nothing below is hosted here. Place each file at the path shown, relative to `$DATA_ROOT`.

### Sequencing data

Five public 10x Genomics datasets. Download the Gene Expression FASTQs from each page into
`fastq/<dataset>/`.

| dataset | FASTQ path | 10x dataset page |
| --- | --- | --- |
| `Flexv2_16k_Human` | `fastq/Flexv2_16k_Human/` | https://www.10xgenomics.com/datasets/Flexv2_16k_Human_PBMCs_TotalseqC_4plex |
| `16k_mouse_Spleno` | `fastq/16k_mouse_Spleno/` | https://www.10xgenomics.com/datasets/Flexv2_16k_mouse_Splenocytes_TotalseqC_4plex |
| `40k_DTC` | `fastq/40k_DTC/` | https://www.10xgenomics.com/datasets/40k-mixture-of-dissociated-tumor-cells-from-3-donors-stained-with-totalseqc-antibodies |
| `16plex_human` | `fastq/16plex_human/` | https://www.10xgenomics.com/datasets/128k-human-pbmcs-stained-with-totalseqc-human-universal-cocktail |
| `320k_GEMX` | `fastq/320k_GEMX/` | https://www.10xgenomics.com/datasets/320k_Human_PBMCs_Sub_Pool_16-plex_GEM-X_FLEX |

### Probe set references

From the 10x Genomics download pages
(https://www.10xgenomics.com/support/software/cell-ranger/downloads). Exact filenames:

| file | used by |
| --- | --- |
| `refs/<dataset>/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv` | `40k_DTC`, `16plex_human`, `320k_GEMX` |
| `refs/<dataset>/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv` | `Flexv2_16k_Human` |
| `refs/<dataset>/Chromium_Mouse_Transcriptome_Probe_Set_v2.0.0_GRCm39-2024-A.csv` | `16k_mouse_Spleno` |

Step 12 `--step a` writes the per-probe annotation table for each panel into
`data/probe_annot/<probe-set basename>.tsv`, which steps 07 and 12 then read.

### Cell barcode whitelists

Shipped inside Cell Ranger under `lib/python/cellranger/barcodes/`. Copy them to:

| file | used by |
| --- | --- |
| `refs/<dataset>/737K-flex-v2.txt` | the two Flex v2 datasets |
| `refs/<dataset>/737K-fixed-rna-profiling.txt` | the three Flex v1 datasets |

### Sample barcode lists

`refs/<dataset>/sample_bc_list.tsv` (and `sample_bc_list_perprobe.tsv` for `320k_GEMX`):
the sample-to-probe-barcode assignment for that library, as published in the `[samples]`
section of its own Cell Ranger multi configuration on its 10x dataset page. The templates
in `data/cr_multi/<dataset>.csv` carry the same assignment.

### Cell Ranger transcriptome references

From https://www.10xgenomics.com/support/software/cell-ranger/downloads, unpacked under
`refs/cr_ref/`:

| file | used by |
| --- | --- |
| `refs/cr_ref/refdata-gex-GRCh38-2024-A` | the four human datasets |
| `refs/cr_ref/refdata-gex-GRCm39-2024-A` | `16k_mouse_Spleno` |

### Mouse-human homology

`HOM_MouseHumanSequence.rpt` from Mouse Genome Informatics
(https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt), placed at
`data/HOM_MouseHumanSequence.rpt` inside this folder, or anywhere `MGI_HOMOLOGY` points.
MGI data and annotations are licensed CC BY 4.0; copyright The Jackson Laboratory.

## Environment you build

No environment file is provided. Install the following yourself and point `site.env` at them.

### Command-line tools

| tool | version |
| --- | --- |
| simpleaf | 0.31.0 |
| alevin-fry | 0.19.0 |
| piscem | 0.24.0 (piscem-rs 0.11.1) |
| cyto | 0.4.6 (built with `data/cyto/Cargo.lock`) |
| Cell Ranger | 10.1.0 — download from https://www.10xgenomics.com/support/software/cell-ranger/downloads and set `CELLRANGER_BIN` |

Also `cargo`, `cc`, `pigz` and `uv` on `PATH`.

### Python

`PY` points at an environment with: scanpy, anndata, numpy, pandas, scipy, scikit-learn,
statsmodels, matplotlib, seaborn, h5py, pysam, rpy2, igraph, leidenalg.

`PY_PROBE` points at an environment that adds CellTypist 1.7.1, pydeseq2 0.5.4 and
harmonypy 0.0.10.

### R

Reachable through rpy2: SoupX, scDblFinder, scran, scry, SingleCellExperiment, Matrix,
BiocParallel.

## Configuration used for the published results

| item | value |
| --- | --- |
| quantification axis | per probe (`multiplex-quant --index <probe-set index>` without `--t2g-map`) |
| index | k = 19, m = 15, `included=FALSE` probes as poison decoys |
| geometry, Flex v1 | `1{b[16]u[12]x:}2{r[50]x[18]s[8]x:}` |
| geometry, Flex v2 | `1{b[16]u[12]x[0-4]hamming(f[TTGCTAGGACCG],1)s[10]x:}2{r:}` |
| barcode correction | frequency, for cell and sample barcodes |
| UMI collapse | one mismatch (`--umi-edit-dist 1`) |
| resolution / decoder / collate codec | cr-like / auto / lz4 |
| threads | 32 for every tool |

Timings were measured on an AWS i4i.16xlarge (Xeon Platinum 8375C, 64 vCPU, 495 GiB) with
its four instance-store NVMe drives in RAID0, each tool in a container limited to 32 CPUs
and 300 GB, reading its inputs from local NVMe.

## Steps

| step | what it does | outputs (`$DUR/<step>/`) |
| --- | --- | --- |
| 00_scratch_and_env | mounts the instance-store RAID0 and checks the environment | `host.txt` |
| 01_toolchain | builds simpleaf, alevin-fry, piscem and cyto at pinned commits; records the Cell Ranger version | `toolchain.tsv` |
| 02_stage_inputs | copies FASTQs, barcode lists and probe sets to scratch | `manifest` |
| 03_index | `simpleaf index` k19/m15 with decoys, for the three probe panels | indices |
| 04_validate_toolchain | checks the Flex v2 matrix digest | `verdict` |
| 05_timing_simpleaf | times simpleaf: 1 warmup, 3 reps, cache eviction; reports minimum runtime and maximum peak RSS | timing TSV, rep JSONs |
| 05b_timing_cyto_cellranger | times cyto and Cell Ranger, and keeps their matrices for steps 07-13 | `results.tsv`, matrices |
| 06_gene_matrices | gene matrices, collated RADs and digests | matrices, `digest_ref.tsv` |
| 07_probe_matrices | probe-level requantification from the collated RADs; per-sample h5ads over the called cells | h5ads, `validation*.json` |
| 08_concordance_counts | per-cell and per-gene count concordance over 44 samples x 3 tool pairs | Supp. Table S1 |
| 09_clustering | per-sample clustering, resolution matching and marker genes | `tbl-biol` |
| 10_pooled_umap | pooled 16-sample UMAPs with matched clusters | supplementary figure, concordance |
| 11_celltypes | CellTypist annotation and marker programs | labels |
| 12_probe_biology | probe annotation, differential probe usage, junction/within-exon ratios, differential expression | Figure 2, Supp. Figure 2 |
| 13_readlevel | Cell Ranger with BAMs on the mouse library, then read-level R1/R2 accounting against it | supplementary tables |
| 14_figures_numbers | every figure of the manuscript, regenerated from the step outputs | figures |

Run them in order with `bash steps/<step>.sh`, or `bash steps/run_all.sh`. Steps 07-13 need
06 and 05b (and 11 for 12), so they can run concurrently once 05b is done, but never during
step 05 or 05b, which need an idle host.

## Authors

This code was written at Altos Labs.

## License

BSD 3-Clause License; see `LICENSE`.

## Attribution

The five sequencing datasets are published by 10x Genomics under a Creative Commons Attribution
4.0 International License (CC BY 4.0); the dataset pages above are the source for each.

`data/cr_multi/<dataset>.csv` are derived from the Cell Ranger `multi` configurations 10x
publishes on those dataset pages.

`HOM_MouseHumanSequence.rpt` is from Mouse Genome Informatics, CC BY 4.0, copyright The Jackson
Laboratory.

`lib/rad_dump/src/main.rs` is derived from alevin-fry
(https://github.com/COMBINE-lab/alevin-fry), BSD 3-Clause License, Copyright (c) 2020-,
University of Maryland.

simpleaf, alevin-fry, piscem and cyto are each BSD 3-Clause licensed and are built from their own
repositories by step 01. Cell Ranger is proprietary 10x Genomics software, downloaded by the user
under 10x's own licence. No part of any of these tools, and no 10x probe-set or reference file, is
redistributed here. `data/cyto/Cargo.lock` is a generated lockfile recording the versions cyto was
built against.
