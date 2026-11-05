#!/usr/bin/env python3
"""Step 13: assemble readlevel_numbers.json (every number of the supplement section
"Sources of residual discordance").

usage: readlevel_numbers.py <durable step dir> <meta.json> <out.json>
Reads a_ledger.json, b_r1_classes.json, c_r2_mapping.json, d_reassignment.json, e_umi_totals.json,
f_cr_probe_bc.json from the step dir. Each number is {value, definition, basis, source}; basis is
"exact" (counted from the tools' own outputs) or "proxy" (a model of a tool's rule, not the tool).
"""
import json
import sys

D, META, OUT = sys.argv[1:4]
J = {k: json.load(open(f"{D}/{k}.json")) for k in
     ("a_ledger", "b_r1_classes", "c_r2_mapping", "d_reassignment", "e_umi_totals", "f_cr_probe_bc")}
a, b, c, d, e, f = (J[k] for k in J)
ct, dt = c["total"], d["total"]
N = a["num_reads"]
U = ct["universe"]
nums = {}


def put(name, value, definition, source, basis="exact"):
    if name in nums:
        raise SystemExit(f"duplicate number {name}")
    nums[name] = {"value": value, "definition": definition, "basis": basis, "source": source}


def pct(x, base):
    return 100.0 * x / base


LIB = "of all library read pairs (num_reads)"
A = "a_ledger.json (step 06's af_quant/map_info.json, af_map/unmapped_bc_count.bin, af_quant/sample_info.json)"

# ---------------- (a) stage-wise read ledger ----------------
put("library_read_pairs", N, "read pairs in the 16k mouse spleen GEX library (piscem num_reads)", A)
put("r1_parse_drop", a["r1_parse_drop"],
    "reads piscem drops while parsing R1 (anchor not found in the window, UMI not ACGT, barcode with >1 N): "
    "num_reads - num_mapped - num_poisoned - sum(unmapped_bc_count.bin)", A)
put("r1_parse_drop_pct", pct(a["r1_parse_drop"], N), "r1_parse_drop, % " + LIB, A)
put("r1_parse_pass_pct", 100.0 * a["r1_parse_pass_rate"],
    "reads that pass R1 geometry parsing (num_reads - r1_parse_drop), % " + LIB + " [L433, was 92.83%]", A)
put("poisoned_by_decoy", a["poisoned_by_decoy"],
    "R1-parsed reads discarded because R2 hits an included=FALSE probe decoy (piscem num_poisoned)", A)
put("poisoned_by_decoy_pct", pct(a["poisoned_by_decoy"], N), "poisoned_by_decoy, % " + LIB, A)
put("r2_no_placement", a["r2_no_placement"],
    "R1-parsed, not poisoned reads whose R2 piscem could not place on any probe (sum of unmapped_bc_count.bin)", A)
put("r2_no_placement_pct", pct(a["r2_no_placement"], N), "r2_no_placement, % " + LIB, A)
put("mapped", a["mapped"], "reads piscem placed on a probe (num_mapped)", A)
put("mapping_rate_pct", 100.0 * a["mapping_rate"], "num_mapped / num_reads x 100 [L439, was 92.24%]", A)
put("sample_bc_rejected", a["sample_bc_rejected"],
    "mapped reads whose probe (sample) barcode no sample-barcode correction could assign "
    "(sample_routing immediate_rejected + deferred_rejected)", A)
put("sample_bc_rejected_pct", pct(a["sample_bc_rejected"], N), "sample_bc_rejected, % " + LIB, A)
put("sample_bc_corrected", a["sample_bc_corrected"],
    "mapped reads assigned to a sample after 1-mismatch frequency correction of the probe barcode "
    "(structurally_unique + deferred_accepted)", A)
put("sample_bc_corrected_pct", pct(a["sample_bc_corrected"], N), "sample_bc_corrected, % " + LIB, A)
put("cell_bc_not_found", a["cell_not_found"],
    "sample-assigned reads whose cell barcode has no whitelisted permit within Hamming 1 (summed over samples)", A)
put("cell_bc_ambiguous", a["cell_ambiguous"],
    "sample-assigned reads whose cell barcode is within Hamming 1 of more than one permit", A)
put("cell_bc_uncorrectable_pct", pct(a["cell_not_found"] + a["cell_ambiguous"], N),
    "(cell_bc_not_found + cell_bc_ambiguous), % " + LIB + " [tbl:r1-parsing row 'Uncorrectable cell barcode']", A)
put("cell_bc_corrected", a["cell_corrected"], "collatable reads whose cell barcode was corrected (1 mismatch)", A)
put("cell_bc_corrected_pct", pct(a["cell_corrected"], N),
    "cell_bc_corrected, % " + LIB + " [tbl:r1-parsing row 'Correctable cell or sample barcode', cell part]", A)
put("collatable", a["collatable"], "reads written to the collated RAD (summed collatable_reads)", A)
put("collatable_pct", pct(a["collatable"], N), "collatable, % " + LIB, A)

# ---------------- (b) R1 read classes ----------------
B = "b_r1_classes.json (lib/r1_geom_scan.c over every R1 read; drop total reconciled exactly to a_ledger)"
g = b["geometry"]
lo, hi = g["window_abs_lo"], g["window_abs_hi"]
ok = b["parse_ok_by_offset_dist"]
by_off = {}
for k, v in ok.items():
    o, dist = (int(x) for x in k.split(":"))
    by_off.setdefault(o, [0, 0])[0 if dist == 0 else 1] += v
outside = b["anchor_absent_by_found_offset_dist"]
at = lambda o: sum(v for k, v in outside.items() if int(k.split(":")[0]) == o)  # noqa: E731
cl, sub = b["classes"], b["anchor_absent_subclasses"]
put("r1_window_anchor_positions", [lo, hi],
    f"0-based start offsets of the R1 anchor accepted by GEOM_V2 (spacer {lo - g['cb_len'] - g['umi_len']}"
    f"-{hi - g['cb_len'] - g['umi_len']} nt), max Hamming distance {g['max_dist']}", B)
put("r1_retained", cl["parse_ok"], f"reads with the anchor at {lo}-{hi} within Hamming {g['max_dist']} "
    "and valid UMI / barcode fields (piscem passes them to R2 mapping)", B)
put("r1_retained_pct", pct(cl["parse_ok"], N), "r1_retained, % " + LIB + " [tbl:r1-parsing row 1]", B)
put("r1_retained_exact_anchor", sum(v[0] for v in by_off.values()), "r1_retained with an exact anchor", B)
put("r1_retained_anchor_1mm", sum(v[1] for v in by_off.values()),
    "r1_retained whose anchor has exactly one mismatch (isolated anchor sequencing error)", B)
put("r1_retained_by_offset", {str(o): {"exact": v[0], "one_mismatch": v[1]} for o, v in sorted(by_off.items())},
    "r1_retained split by anchor start offset and anchor mismatches", B)
put("r1_anchor_before_window", sub["anchor_before_window"],
    f"no anchor at {lo}-{hi}; the first anchor occurrence in the read (exact, else within Hamming "
    f"{g['max_dist']}) starts before {lo} (deletion-shifted barcode/UMI) [tbl:r1-parsing 'Anchor before position {lo}']", B)
put("r1_anchor_before_window_pct", pct(sub["anchor_before_window"], N), "r1_anchor_before_window, % " + LIB, B)
put("r1_anchor_at_window_lo_minus_1", at(lo - 1),
    f"subset of r1_anchor_before_window with the anchor at {lo - 1}: the spacer -1 class (UMI interval holds "
    f"only {g['umi_len'] - 1} nt), EXCLUDED BY DESIGN", B)
put("r1_anchor_at_window_lo_minus_1_pct", pct(at(lo - 1), N), "r1_anchor_at_window_lo_minus_1, % " + LIB, B)
put("r1_anchor_after_window", sub["anchor_after_window"],
    f"no anchor at {lo}-{hi}; first anchor occurrence starts after {hi} (spacer longer than the chemistry range)", B)
put("r1_anchor_after_window_pct", pct(sub["anchor_after_window"], N), "r1_anchor_after_window, % " + LIB, B)
put("r1_anchor_at_window_hi_plus_1", at(hi + 1),
    f"subset of r1_anchor_after_window with the anchor at {hi + 1}: the spacer-{hi + 1 - g['cb_len'] - g['umi_len']} "
    "class, EXCLUDED BY DESIGN (outside GEOM_V2's spacer range; with a 54-nt R1 and x[0-4] its probe barcode would "
    "also run past the read end, so widening the window would not retain it)", B)
put("r1_anchor_at_window_hi_plus_1_pct", pct(at(hi + 1), N), "r1_anchor_at_window_hi_plus_1, % " + LIB, B)
put("r1_no_anchor_within_hamming", sub["no_anchor_within_max_dist"],
    f"reads with no anchor occurrence within Hamming {g['max_dist']} anywhere in R1 [tbl:r1-parsing]", B)
put("r1_no_anchor_within_hamming_pct", pct(sub["no_anchor_within_max_dist"], N),
    "r1_no_anchor_within_hamming, % " + LIB, B)
put("r1_umi_invalid", cl["umi_invalid"], "anchor found but the UMI holds a non-ACGT base", B)
put("r1_sample_bc_field_invalid", cl["sample_bc_invalid"],
    "anchor found but the probe-barcode field is truncated or holds >1 N / a non-ACGTN base", B)
put("r1_cell_bc_field_invalid", cl["cell_bc_invalid"], "anchor found but the cell-barcode field holds >1 N", B)
put("r1_reconciliation_difference", b["reconciliation"]["difference"],
    "scanner R1 drop minus a_ledger r1_parse_drop; must be 0", B)

# ---------------- (c) R2 bidirectional mapping ----------------
C = ("c_r2_mapping.json (each CR BAM read re-mapped by the step-01 piscem with the step-03 index; "
     "synthetic R1 carries the read ordinal)")
put("r2_universe", U,
    "primary, mapped CR BAM records with a pr probe tag, a 50-nt SEQ and the probe in the simpleaf reference "
    "(the reads Cell Ranger assigned to an included probe) [L480, was 1.06 billion]", C)
put("r2_simpleaf_only", ct["simpleaf_only"],
    "universe reads piscem places on CR's assigned gene although they fail CR's idealized first-stage rule "
    "(each 25-nt probe half within 1 substitution, ungapped). CR assigned these reads too (stage-2 rescue), "
    "so this is NOT 'discarded by Cell Ranger' [tbl:bidirectional-mapping row 1, was 13,816,526]", C)
put("r2_simpleaf_only_pct", pct(ct["simpleaf_only"], U), "r2_simpleaf_only, % of r2_universe [was 1.3%]", C)
put("r2_simpleaf_only_one_half_exact", ct["simpleaf_only_one_half_exact_other_gt1_sub"],
    "r2_simpleaf_only with one half exact and the other >1 substitution", C)
put("r2_simpleaf_only_both_halves_gt1", ct["simpleaf_only_both_halves_gt1_sub"],
    "r2_simpleaf_only with both halves >1 substitution", C)
put("r2_cr_only", ct["cr_only"],
    "universe reads that pass CR's idealized first-stage rule but piscem does not place on CR's gene "
    "[tbl:bidirectional-mapping row 2, was 1,153,283]", C)
put("r2_cr_only_pct", pct(ct["cr_only"], U), "r2_cr_only, % of r2_universe [was 0.11%]", C)
put("r2_cr_only_unplaced", ct["cr_only_unplaced"], "r2_cr_only that piscem placed nowhere", C)
put("r2_cr_only_poisoned", ct["cr_only_poisoned"], "r2_cr_only discarded by a decoy hit", C)
put("r2_cr_only_other_gene", ct["cr_only_other_gene_only"], "r2_cr_only placed only on other genes", C)
put("r2_cr_only_two_single_subs", ct["cr_only_halves_1_1_sub"], "r2_cr_only with exactly 1 substitution in each half", C)
put("r2_agree_keep", ct["agree_keep"], "universe reads kept by both rules", C)
put("r2_proxy_simpleaf_only", ct["proxy"]["simpleaf_only"],
    f"as r2_simpleaf_only but simpleaf modelled as 'an exact {c['k']}-mer of the assigned gene'",
    C, basis="proxy")
put("r2_proxy_cr_only", ct["proxy"]["cr_only"], "as r2_cr_only with the exact-k-mer model of simpleaf", C, basis="proxy")

# ---------------- (d) gene reassignment ----------------
DD = "d_reassignment.json (rad_dump mol on step 06's collated RAD vs CR BAM fx/UR molecules, CR called cells)"
put("mol_cr_molecules", dt["cr_molecules"],
    "distinct (sample, cell, raw UR UMI, gene) of reads CR confidently assigned (fx tag), CR called cells, fx in "
    "CR's FILTERED features.tsv.gz (19,070 genes = the simpleaf reference genes; the 333 excluded-only genes are "
    "not on this axis) [was 254,599,334]", DD)
put("mol_reassigned", dt["reassigned_in"],
    "simpleaf molecules whose (cell, UMI) is a CR molecule but carry a gene CR never gave that UMI "
    "(gene reassignment) [L480, was 28,870]", DD)
put("mol_reassigned_pct", pct(dt["reassigned_in"], dt["cr_molecules"]), "mol_reassigned, % of mol_cr_molecules [was 0.011%]", DD)
put("mol_reassigned_out_cr_side", dt["reassigned_out"],
    "CR (key, gene) whose simpleaf molecule has another gene; counts every gene of a multi-gene CR UMI "
    "(inflated by the pre-collapse triple convention, see r1-anchor memory)", DD)
put("mol_lost_absent", dt["lost_absent"], "CR molecules with no simpleaf molecule at the same key", DD)
put("mol_lost_tie", dt["lost_tie"], "CR molecules whose simpleaf UMI is a gene tie", DD)
put("mol_simpleaf_novel", dt["novel"], "simpleaf molecules at keys CR has no molecule for", DD)
put("mol_simpleaf_molecules", dt["simpleaf_molecules"], "non-tied simpleaf molecules in CR called cells", DD)

# ---------------- (e) UMI totals ----------------
E = "e_umi_totals.json (UMI totals: step 08 per-sample JSONs summed over the 4 samples; means: step 08 text_ranges.json)"
for tool in ("cellranger", "simpleaf", "cyto"):
    put(f"umi_total_{tool}", e["totals"][tool],
        f"{tool} UMIs in CR's filtered cells over the genes shared with " +
        ("simpleaf (cellranger_vs_simpleaf ref)" if tool == "cellranger" else "Cell Ranger") + " [L126]", E)
for tool in ("simpleaf", "cyto"):
    put(f"umi_diff_{tool}_vs_cellranger_pct", e["diff_pct"][tool],
        f"(umi_total_{tool} / CR total over the same genes - 1) x 100", E)
put("count_concordance_cr_vs_simpleaf_mean_pct", e["mean_per_sample"]["count_concordance_pct"],
    "mean over the 4 samples of compare_consistency count_concordance_pct (1 - L1/total) [L513]", E)
put("gene_spearman_cr_vs_simpleaf_mean", e["mean_per_sample"]["gene_spearman"],
    "mean over the 4 samples of the per-gene Spearman of summed counts [L513]", E)

# ---------------- (f) Cell Ranger valid probe barcodes ----------------
put("cr_valid_probe_barcodes_pct", 100.0 * f["valid_probe_barcodes"],
    "Cell Ranger metrics_summary 'Valid probe barcodes' (library level, identical in every per-sample file) [L439, was 95.87%]",
    "f_cr_probe_bc.json ($CR_BAM_OUT and $CR_OUT metrics_summary.csv)")

out = {"meta": json.load(open(META)), "numbers": nums}
json.dump(out, open(OUT, "w"), indent=2)
for k, v in nums.items():
    val = v["value"]
    print(f"{k:42s} {val if not isinstance(val, float) else round(val, 6)}{'  (proxy)' if v['basis'] == 'proxy' else ''}")
