# Disability-CKM-trends

Stata code for: **Trends in Cardiometabolic Risk Factor Prevalence, Treatment, and Control Among US Adults With Disability, NHANES 1999-2018.**

The analysis uses public-use data from the National Health and Nutrition Examination Survey (NHANES), 10 cycles from 1999-2000 through 2017-2018 (N = 53,540 adults; 18,158 with disability). Disability was defined from the Physical Functioning (PFQ, 1999-2012) and Disability (DLQ, 2013-2018) questionnaires.

## Data

The harmonized analysis dataset (`NHANES_Harmonized_data.dta`) is included in this repository. It contains only variables from public-use NHANES files (CDC/NCHS) and variables derived from them; definitions are in the manuscript Supplement. Use of the data is subject to the NCHS Data User Agreement: https://www.cdc.gov/nchs/policy/data-user-agreement.html. The original files are available from the CDC: https://wwwn.cdc.gov/nchs/nhanes/.

Two variables in the dataset carry older labels that do not match the analysis. `hyperlipid_controlled` and `joint_ctrl_all3` are labelled with a non-HDL cholesterol threshold of <160 mg/dL. The scripts do not use them. They use `lipid_controlled` (non-HDL cholesterol <130 mg/dL) and construct joint control of all three conditions from `dm_controlled`, `htn_controlled` and `lipid_controlled`.

## Requirements

Stata 17 or later. All estimates use the survey commands (`svy`) with examination weights, strata (SDMVSTRA) and PSU (SDMVPSU), and `singleunit(centered)`. The multiple imputation script uses Stata's `mi` commands.

## How to run

1. Download or clone this repository.
2. Open `00_config.do`, set `ROOT` to the repository folder, and run it.
3. Run any script below. Each script starts from the harmonized dataset and writes its log and result files to `results/`. The scripts do not depend on one another, so they can be run in any order after `00_config.do`.

Weights: the PFQ era (cycles 1-7) uses `wt_pfq` and the DLQ era (cycles 8-10) uses `wt_dlq`. Pooled and per-cycle age-standardized estimates use the examination weight `wtmec2yr`.

## Scripts

| Script | Manuscript item |
|---|---|
| `00_config.do` | Sets file paths (run first) |
| `02_table1_characteristics.do` | Table 1 |
| `03_table2_prevalence_trends.do` | Table 2, Supplemental Table S1 |
| `04_table3_treatment_control.do` | Table 3 (treatment and control), Supplemental Table S2, unweighted N |
| `s3_comparator.do` | Supplemental Table S3 (adults with and without disability) |
| `s4_comorbidity.do` | Supplemental Table S4 |
| `s5_fdr_correction.do` | Supplemental Table S5 (Benjamini-Hochberg false discovery rate) |
| `s6_sensitivity_analyses.do` | Supplemental Table S6 (S6a to S6d) |
| `s7_segmented_trends.do` | Supplemental Table S7 |
| `s8_evalues.do` | Supplemental Table S8 and Supplemental Figure S3 (E-values) |
| `s9_adjusted_trends.do` | Supplemental Table S9 (S9a to S9e) |
| `s10_multiple_imputation_alcohol.do` | Supplemental Table S10 |

## Citation

If you use this code, please cite the associated article (reference to be added on publication).
