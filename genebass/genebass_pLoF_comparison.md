# GeneBass pLoF File Comparison: `genebass_pLoF_all.pkl` vs `genebass_pLoF_filtered.pkl`

## Overview

This document compares the two GeneBass pLoF data files and describes the
filtering applied to produce the filtered version. An accompanying script
[`compare_pLoF_files.py`](compare_pLoF_files.py) reproduces the analysis below.

## File Origins

The two files come from **different sources**:

| File | Source | How produced |
|------|--------|-------------|
| `genebass_pLoF_all.pkl` | Raw Hail MatrixTable at `gs://ukbb-exome-public/500k/results/results.mt` | Downloaded locally via [`download_genebass.py`](download_genebass.py) |
| `genebass_pLoF_filtered.pkl` | QC-filtered results, likely from GeneBass's `get_qc_result_mt()` pipeline | Downloaded from a separate processed source; **no script in this repo produces it** |

There is no script in this repository that creates the filtered file from the
all file. The filtering was applied upstream, most likely by the GeneBass
pipeline code at [`Nealelab/ukb_exomes`](https://github.com/Nealelab/ukb_exomes)
(specifically the `get_qc_result_mt()` function in `utils/results.py`).

## Dimensions

| Metric | All | Filtered | Change |
|--------|-----|----------|--------|
| Rows | 78,342,883 | 24,383,843 | -68.9% |
| Columns | 16 | 8 | -8 |
| Genes | 18,356 | 18,092 | -264 |
| Phenotypes | 2,452 (phenocode+coding combos) | 1,414 (pheno_descriptions) | ~1,038 removed |

## Column Changes

### Columns in `_all` but not in `_filtered` (dropped)

```
gene_id, gene_symbol, trait_type, phenocode, pheno_sex, coding, modifier,
description, n_cases, n_controls
```

### Columns in `_filtered` but not in `_all` (renamed/added)

```
gene            (renamed from gene_symbol)
pheno_description  (renamed from description, with reformatting)
```

### Shared columns (unchanged)

```
annotation, Pvalue, Pvalue_Burden, Pvalue_SKAT, BETA_Burden, SE_Burden
```

### `pheno_description` format

The `pheno_description` column in the filtered file was constructed from
multiple fields in the raw data:

- **Continuous traits**: just the `description` value (e.g., `"Vitamin D"`)
- **Categorical traits**: `description; coding_meaning` (e.g., `"Alcohol drinker status; Never"`)
- **ICD codes**: the ICD description (e.g., `"Alcoholic cirrhosis"`)
- **Custom phenotypes**: a custom label (e.g., `"Alzheimers_custom1"`, `"FH_Diabetes_custom"`)

## Verified Filters

These filters were **confirmed by direct inspection** of the data:

### 1. `modifier != 'raw'`

Zero rows with `modifier == 'raw'` appear in the kept phenotypes. Confirmed.

### 2. `modifier != '04162021'`

Zero rows with `modifier == '04162021'` are in the kept set. All 7 such
phenocode+coding combos are excluded. Confirmed.

### 3. No "Source of report" descriptions

Zero phenotypes containing "Source of report" in their description appear
in the kept set. Confirmed.

### 4. Biomarkers (phenocode starting with `30`) preferentially included

64 out of 68 biomarker phenotypes are in the filtered file. Confirmed.

### 5. Phenotype selection is a curated list, NOT a simple numeric cutoff

This is the most important finding. The phenotype filtering is **not** based
on a single threshold (like `n_cases >= 200`). Evidence:

- "Smoking status" with 213,273 cases was **excluded**
- 290 phenocode+coding combos with `n_cases < 200` were **kept** (mostly
  individual codings of categorical traits like cancer codes, medication codes,
  and rare ICD first-occurrence dates)

The GeneBass source code (`Nealelab/ukb_exomes`, `hail/load_phenotype_data.py`)
reveals the actual phenotype selection uses a multi-criteria approach:

- An **external exclusion list** at `gs://ukbb-pharma-exome-analysis/500k/450K_phentoypes_removeflag_vF04132021.txt`
- A **priority score** assigned by pharma partners (score >= 1 to be included)
- Automatic inclusion of `icd_first_occurrence` trait types
- Automatic inclusion of biomarker phenocodes (starting with `30`)
- Automatic inclusion of custom pharma phenotypes (`biogen`, `abbvie`, `pfizer` modifiers)
- Exclusion of `modifier == 'raw'` and descriptions containing "Source of report"
- A minimum case count of `n_cases_defined >= 200` (applied as one of several
  criteria, not as the sole filter)
- Post-filtering removal of `n_cases_defined < 100`

### 6. Twenty-three custom phenotypes added

The filtered file contains 23 phenotypes that **do not exist** in the raw
MatrixTable. These are custom phenotype definitions from pharma partners:

```
Alzheimers_custom1
Corneal_hysteresis_left_custom
Corneal_hysteresis_right_custom
Corneal_resistance_factor_left_custom
Corneal_resistance_factor_right_custom
Depressive_symptoms_custom
FH_Alzheimer_disease_dementia_custom
FH_Bowel_cancer_custom
FH_Breast_cancer_custom
FH_Diabetes_custom
FH_Heart_disease_custom
FH_Lung_cancer_custom
FH_Parkinsons_disease_custom
FH_Prostate_cancer_custom
Glaucoma_custom
IOP_Corneal_left_custom
IOP_Corneal_right_custom
IOP_Goldmann_left_custom
IOP_Goldmann_right_custom
Ischemic_stroke_custom
Parkinsons_custom
Stroke_custom
Touchscreen_duration_custom
```

### 7. `pheno_sex` is mostly but not exclusively `both_sexes`

Row-level matching of 500 sampled filtered rows to the all file shows:
- 451/468 matched rows are `both_sexes` (96.4%)
- 10 are `males`, 7 are `females`
- 32 could not be matched (custom phenotypes)

This means sex-specific analyses are included for some phenotypes.

### 8. Trait type breakdown

| Trait type | Kept | Excluded |
|-----------|------|----------|
| categorical | 2,462 | 106 |
| continuous | 472 | 719 |
| icd_first_occurrence | 222 | 458 |
| icd10 | 29 | 16 |

Most categorical traits are kept. The majority of excluded phenotypes are
continuous traits and ICD first-occurrence codes.

### 9. Gene filtering

264 genes in `_all` are absent from `_filtered`. These are predominantly
non-standard gene symbols (AC-prefix lncRNAs, read-through transcripts, and
pseudogenes). 242 genes appear in `_filtered` but not in `_all`, suggesting
different gene symbol mappings between the two source MatrixTables.

## Unverifiable Filters

The following filters are described in the GeneBass source code
(`get_qc_result_mt()` in `Nealelab/ukb_exomes/utils/results.py`) but
**cannot be verified** from the two pickle files because the required QC
metadata columns are not present in either file:

| Filter | Threshold | Column needed |
|--------|-----------|---------------|
| Gene mean coverage | >= 20X | `mean_coverage` |
| Lambda GC (synonymous) | >= 0.75 | `synonymous_lambda_gc_skato` |
| Expected allele count | >= 50 | `expected_ac` |
| Minimum variants per gene | >= 2 | `total_variants` |
| Phenotype lambda GC | >= 0.75 | `lambda_gc_skato` |
| Correlated phenotype removal | r^2 >= 0.5 | `keep_pheno_unrelated` |
| Entry-level expected AC | >= 50 | `expected_AC` (per gene-phenotype) |

These QC annotations are pre-computed and stored in a separate MatrixTable at:
```
gs://ukbb-exome-public/500k/qc/gene_qc_metrics_ukb_exomes_500k.mt
```

To verify these filters, one would need to load that QC MatrixTable via Hail
and cross-reference the boolean QC flags against the gene/phenotype lists
in the filtered file.

## GeneBass Source Code References

The filtering pipeline lives in the [`Nealelab/ukb_exomes`](https://github.com/Nealelab/ukb_exomes) repository:

| File | Function | Purpose |
|------|----------|---------|
| `utils/results.py` | `get_qc_result_mt()` | Applies all QC filters to produce final results |
| `utils/results.py` | `annotate_qc_metric_mt()` | Defines QC thresholds and annotates boolean flags |
| `utils/results.py` | `modify_phenos_mt()` | Post-run phenotype removal (n_cases < 100, specific phenocodes) |
| `hail/load_phenotype_data.py` | (main logic) | Phenotype whitelist: priority scores, exclusion list, cancer codes |
| `resources/generic.py` | `MIN_CASES` | Set to 200 |

## How to Reproduce

Run the comparison script:

```bash
cd genebass
python compare_pLoF_files.py --data-dir ./genebass_output
```

This will print the full analysis to stdout. Requires `pandas` and `numpy >= 2.0`
(the pickle files were saved with numpy 2.x).
