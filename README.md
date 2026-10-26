# Disability-CKM-trends

Stata code for: **Trends in Cardiometabolic Risk Factor Prevalence, Treatment, and Control Among US Adults With Disability, NHANES 1999-2018.**

The analysis uses public-use data from the National Health and Nutrition Examination Survey (NHANES), 10 cycles from 1999-2000 through 2017-2018 (N = 53,540 adults; 18,158 with disability). Disability was defined from the Physical Functioning (PFQ, 1999-2012) and Disability (DLQ, 2013-2018) questionnaires.

## Data

The harmonized analysis dataset (`harmonized/nhanes_harmonized_primary_v2.dta`) is included in this repository. It contains only variables from public-use NHANES files (CDC/NCHS) and variables derived from them; definitions are in the manuscript Supplement. Use of the data is subject to the NCHS Data User Agreement: https://www.cdc.gov/nchs/policy/data-user-agreement.html. The original files are available from the CDC: https://wwwn.cdc.gov/nchs/nhanes/.

## Requirements

Stata 17 or later. All estimates use the survey commands (`svy`) with examination weights, strata (SDMVSTRA) and PSU (SDMVPSU), and `singleunit(centered)`.

## How to run

1. Download or clone this repository.
2. Open `00_config.do`, set `ROOT` to the repository folder, and run it.
3. Run the scripts in numerical order. Each script starts from the harmonized dataset and writes its log to `results/`.

## Scripts

| Script | Manuscript item |
|---|---|
| `02_table1_characteristics.do` | Table 1 |
| `03_table2_prevalence_trends.do` | Table 2, Supplemental Table S1 |
| `04_table3_treatment_control.do` | Table 3, Supplemental Table S2 |
| `05_table3_joint_control_trends.do` | Table 3 (joint control trends) |
| `06_table3_unweighted_n.do` | Unweighted N for Table 3 and Table S2 |
| `supplement/s3_comparator_disability.do` | Supplemental Table S3 (disability group) |
| `supplement/s3_comparator_nondisabled.do` | Supplemental Table S3 (non-disabled group) |
| `supplement/s4_comorbidity.do` | Supplemental Table S4 |
| `supplement/s5_fdr_correction.do` | Supplemental Table S5 (false discovery rate) |
| `supplement/s6a_age_group_sensitivity.do` | Supplemental Table S6a |
| `supplement/s6d_mobility_definition.do` | Supplemental Table S6d |
| `supplement/s7_segmented_trends.do` | Supplemental Table S7 |
| `supplement/s8_evalues.do` | E-values |
| `supplement/s9_adjusted_trends.do` | Supplemental Tables S9a, S9c, S9d |
| `supplement/s9_insurance_sensitivity.do` | Supplemental Tables S9b, S9e |
| `supplement/s10_multiple_imputation_alcohol.do` | Supplemental Table S10 |

## Citation

If you use this code, please cite the associated article (reference to be added on publication).
