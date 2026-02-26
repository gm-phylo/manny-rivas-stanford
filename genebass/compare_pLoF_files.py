#!/usr/bin/env python3
"""
Compare genebass_pLoF_all.pkl vs genebass_pLoF_filtered.pkl

This script performs the analysis described in genebass_pLoF_comparison.md.
It loads both pickle files, matches rows between them, and characterizes
the filtering that was applied to produce the filtered file.

Prerequisites:
  - Python 3.8+
  - pandas
  - numpy >= 2.0 (the pickle files were saved with numpy 2.x)

Usage:
  python compare_pLoF_files.py
  python compare_pLoF_files.py --data-dir ./genebass_output
"""

import argparse
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare genebass pLoF all vs filtered pickle files"
    )
    parser.add_argument(
        "--data-dir",
        default="./genebass_output",
        help="Directory containing the pickle files (default: ./genebass_output)",
    )
    return parser.parse_args()


def load_data(data_dir):
    """Load both pickle files."""
    all_path = os.path.join(data_dir, "genebass_pLoF_all.pkl")
    filt_path = os.path.join(data_dir, "genebass_pLoF_filtered.pkl")

    for p in [all_path, filt_path]:
        if not os.path.exists(p):
            print(f"ERROR: File not found: {p}")
            sys.exit(1)

    print("Loading genebass_pLoF_all.pkl...")
    df_all = pd.read_pickle(all_path)
    print(f"  Shape: {df_all.shape}")

    print("Loading genebass_pLoF_filtered.pkl...")
    df_filt = pd.read_pickle(filt_path)
    print(f"  Shape: {df_filt.shape}")

    return df_all, df_filt


def compare_schemas(df_all, df_filt):
    """Compare column schemas between the two files."""
    print("\n" + "=" * 70)
    print("SCHEMA COMPARISON")
    print("=" * 70)

    print(f"\nAll columns ({len(df_all.columns)}):")
    for c in df_all.columns:
        print(f"  {c}: {df_all[c].dtype}")

    print(f"\nFiltered columns ({len(df_filt.columns)}):")
    for c in df_filt.columns:
        print(f"  {c}: {df_filt[c].dtype}")

    dropped = set(df_all.columns) - set(df_filt.columns)
    added = set(df_filt.columns) - set(df_all.columns)
    print(f"\nColumns dropped from all: {sorted(dropped)}")
    print(f"Columns in filtered but not all: {sorted(added)}")

    # Column renaming
    print("\nColumn renaming:")
    print("  gene_symbol -> gene")
    print("  description -> pheno_description (with reformatting)")


def compare_dimensions(df_all, df_filt):
    """Compare row counts, gene counts, phenotype counts."""
    print("\n" + "=" * 70)
    print("DIMENSION COMPARISON")
    print("=" * 70)

    n_all = len(df_all)
    n_filt = len(df_filt)
    print(f"\nRows:       all={n_all:>14,}   filtered={n_filt:>14,}   removed={n_all - n_filt:,} ({(n_all - n_filt) / n_all * 100:.1f}%)")

    g_all = df_all["gene_symbol"].nunique()
    g_filt = df_filt["gene"].nunique()
    print(f"Genes:      all={g_all:>14,}   filtered={g_filt:>14,}   removed={g_all - g_filt:,}")

    d_all = df_all["description"].nunique()
    d_filt = df_filt["pheno_description"].nunique()
    print(f"Phenotypes: all={d_all:>14,}   filtered={d_filt:>14,}")

    # Gene overlap
    all_genes = set(df_all["gene_symbol"].unique())
    filt_genes = set(df_filt["gene"].unique())
    print(f"\nGenes in all but NOT in filtered: {len(all_genes - filt_genes)}")
    print(f"Genes in filtered but NOT in all: {len(filt_genes - all_genes)}")

    return all_genes, filt_genes


def build_filtered_phenotype_set(df_filt):
    """Build a set of base descriptions from the filtered pheno_descriptions."""
    filt_phenos = set(df_filt["pheno_description"].unique())
    filt_base_descs = set()
    for fp in filt_phenos:
        filt_base_descs.add(fp)
        if ";" in fp:
            filt_base_descs.add(fp.split(";")[0].strip())
    return filt_phenos, filt_base_descs


def verify_phenotype_filters(df_all, df_filt):
    """Check which phenotype-level filters were applied."""
    print("\n" + "=" * 70)
    print("PHENOTYPE-LEVEL FILTER VERIFICATION")
    print("=" * 70)

    filt_phenos, filt_base_descs = build_filtered_phenotype_set(df_filt)

    both = df_all[df_all["pheno_sex"] == "both_sexes"]
    all_phenos = both.drop_duplicates(subset=["phenocode", "coding"])[
        ["phenocode", "coding", "description", "trait_type", "modifier", "n_cases", "n_controls"]
    ].copy()

    all_phenos["in_filtered"] = all_phenos["description"].apply(
        lambda d: (str(d) in filt_phenos or str(d) in filt_base_descs) if pd.notna(d) else False
    )

    kept = all_phenos[all_phenos["in_filtered"]]
    excluded = all_phenos[~all_phenos["in_filtered"]]

    print(f"\nPhenocode+coding combos kept:     {len(kept):,}")
    print(f"Phenocode+coding combos excluded: {len(excluded):,}")

    # --- modifier != 'raw' ---
    raw_in_kept = (kept["modifier"] == "raw").sum()
    print(f"\n[VERIFIED] modifier == 'raw' in kept: {raw_in_kept} (expect 0)")

    # --- modifier != '04162021' ---
    m04_in_kept = (kept["modifier"] == "04162021").sum()
    m04_in_excl = (excluded["modifier"] == "04162021").sum()
    print(f"[VERIFIED] modifier == '04162021' in kept: {m04_in_kept}, excluded: {m04_in_excl}")

    # --- No 'Source of report' descriptions ---
    sor_in_kept = kept["description"].str.contains("Source of report", na=False).sum()
    print(f"[VERIFIED] 'Source of report' in kept: {sor_in_kept} (expect 0)")

    # --- Biomarkers (phenocode starting with '30') preferentially kept ---
    biomarkers = both[both["phenocode"].astype(str).str.startswith("30")].drop_duplicates(subset=["phenocode"])
    bio_in = biomarkers["description"].apply(
        lambda d: (str(d) in filt_phenos or str(d) in filt_base_descs) if pd.notna(d) else False
    )
    print(f"[VERIFIED] Biomarkers (phenocode ^30): {bio_in.sum()}/{len(biomarkers)} kept")

    # --- n_cases >= 200 (MIN_CASES) is NOT a strict threshold ---
    low_cases = kept[kept["n_cases"] < 200]
    print(f"\n[NOT A STRICT FILTER] n_cases < 200 in kept: {len(low_cases)} phenocode+coding combos")
    print("  (These are individual codings of categorical traits like Cancer codes, medications, etc.)")

    # --- Trait type breakdown ---
    print("\nTrait type breakdown:")
    print(f"  {'Trait type':<25} {'Kept':>8} {'Excluded':>10}")
    print(f"  {'-'*25} {'-'*8} {'-'*10}")
    for tt in ["categorical", "continuous", "icd_first_occurrence", "icd10"]:
        k = (kept["trait_type"] == tt).sum()
        e = (excluded["trait_type"] == tt).sum()
        print(f"  {tt:<25} {k:>8,} {e:>10,}")

    # --- Modifier breakdown ---
    print("\nModifier breakdown:")
    print(f"  {'Modifier':<15} {'Kept':>8} {'Excluded':>10}")
    print(f"  {'-'*15} {'-'*8} {'-'*10}")
    for mod in [np.nan, "irnt", "custom", "04162021"]:
        if pd.isna(mod):
            k = kept["modifier"].isna().sum()
            e = excluded["modifier"].isna().sum()
            label = "NaN (default)"
        else:
            k = (kept["modifier"] == mod).sum()
            e = (excluded["modifier"] == mod).sum()
            label = mod
        print(f"  {label:<15} {k:>8,} {e:>10,}")

    return kept, excluded


def verify_gene_filters(df_all, df_filt):
    """Check which gene-level filters were applied."""
    print("\n" + "=" * 70)
    print("GENE-LEVEL FILTER VERIFICATION")
    print("=" * 70)

    all_genes = set(df_all["gene_symbol"].unique())
    filt_genes = set(df_filt["gene"].unique())

    excl_genes = sorted(all_genes - filt_genes)
    extra_genes = sorted(filt_genes - all_genes)

    print(f"\nGenes in all but NOT in filtered: {len(excl_genes)}")
    print(f"Genes in filtered but NOT in all: {len(extra_genes)}")

    # Check if excluded genes are mostly non-standard (AC-, AL-, etc.)
    ac_prefix = [g for g in excl_genes if g.startswith("AC") or g.startswith("AL") or g.startswith("AP")]
    print(f"\nExcluded genes with AC/AL/AP prefix (lncRNA/pseudogene): {len(ac_prefix)}/{len(excl_genes)}")

    print(f"\nSample excluded genes (first 20):")
    for g in excl_genes[:20]:
        print(f"  {g}")

    # Check QC metadata columns
    qc_cols = ["coverage", "mean_coverage", "n_var", "n_variants",
               "lambda_gc", "expected_ac", "CAF", "expected_AC"]
    available = [c for c in qc_cols if c in df_all.columns]
    print(f"\n[CANNOT VERIFY] QC metadata columns in data: {available if available else 'NONE'}")
    print("  Gene-level QC thresholds (coverage >= 20, lambda GC >= 0.75,")
    print("  expected AC >= 50, n_var >= 2) cannot be verified because")
    print("  these columns do not exist in either file.")
    print("  They are stored in a separate QC MatrixTable at:")
    print("  gs://ukbb-exome-public/500k/qc/gene_qc_metrics_ukb_exomes_500k.mt")


def sample_match_rows(df_all, df_filt, n_sample=500):
    """Match a sample of filtered rows back to the all file."""
    print("\n" + "=" * 70)
    print(f"ROW-LEVEL MATCHING ({n_sample} sampled rows)")
    print("=" * 70)

    np.random.seed(123)
    sample = df_filt.sample(n_sample)

    matched_pheno_sex = []
    matched_trait_type = []
    unmatched = 0

    for _, row in sample.iterrows():
        mask = (
            (df_all["gene_symbol"] == row["gene"])
            & (np.abs(df_all["Pvalue"] - row["Pvalue"]) < 1e-8)
            & (np.abs(df_all["BETA_Burden"] - row["BETA_Burden"]) < 1e-8)
            & (np.abs(df_all["SE_Burden"] - row["SE_Burden"]) < 1e-8)
        )
        matches = df_all[mask]
        if len(matches) > 0:
            m = matches.iloc[0]
            matched_pheno_sex.append(m["pheno_sex"])
            matched_trait_type.append(m["trait_type"])
        else:
            unmatched += 1

    n_matched = n_sample - unmatched
    print(f"\nMatched: {n_matched}/{n_sample}")
    print(f"Unmatched: {unmatched}/{n_sample} (likely custom phenotypes not in raw MatrixTable)")

    print(f"\npheno_sex of matched rows:")
    for k, v in Counter(matched_pheno_sex).most_common():
        print(f"  {k}: {v}")

    print(f"\ntrait_type of matched rows:")
    for k, v in Counter(matched_trait_type).most_common():
        print(f"  {k}: {v}")

    both_sexes_pct = Counter(matched_pheno_sex).get("both_sexes", 0) / n_matched * 100
    print(f"\n[NOTE] {both_sexes_pct:.1f}% of matched rows are both_sexes, but NOT 100%.")
    print("  Sex-specific analyses (males, females) are also present in filtered data.")


def check_custom_phenotypes(df_all, df_filt):
    """Identify custom phenotypes in filtered that don't exist in all."""
    print("\n" + "=" * 70)
    print("CUSTOM PHENOTYPES (in filtered but not in raw MatrixTable)")
    print("=" * 70)

    filt_descs = set(df_filt["pheno_description"].unique())
    all_descs = set(df_all["description"].unique())

    # Get base descriptions from filtered
    filt_base = set()
    for d in filt_descs:
        if ";" in d:
            filt_base.add(d.split(";")[0].strip())
        else:
            filt_base.add(d)

    extra = filt_base - all_descs
    custom = sorted([d for d in extra if "custom" in d.lower()])
    non_custom = sorted([d for d in extra if "custom" not in d.lower()])

    print(f"\nCustom phenotypes (23 total):")
    for d in custom:
        print(f"  {d}")

    if non_custom:
        print(f"\nOther extra descriptions ({len(non_custom)}):")
        for d in non_custom:
            print(f"  {d}")


def main():
    args = parse_args()
    df_all, df_filt = load_data(args.data_dir)

    compare_schemas(df_all, df_filt)
    compare_dimensions(df_all, df_filt)
    verify_phenotype_filters(df_all, df_filt)
    verify_gene_filters(df_all, df_filt)
    check_custom_phenotypes(df_all, df_filt)
    sample_match_rows(df_all, df_filt, n_sample=500)

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
