*=============================================================================
* s10_multiple_imputation_alcohol.do
* Supplemental Table S10: multiple imputation for missing high-risk alcohol
* use, adults with disability, NHANES 1999-2018.
*
* Requires: global DATA (set in 00_config.do or in the Stata command window)
*
* Part A: disability x cycle interaction for high-risk alcohol use, complete
*         case, both groups pooled (PFQ era cycles 1-7, DLQ era cycles 8-10).
* Part B: missingness on high-risk alcohol use among adults with disability.
* Part C: complete-case trend models (adults with disability), as in
*         Table 2 / Supplemental Table S1.
* Part D: multiple imputation by chained equations (logistic model with survey
*         cycle, age group, sex and race/ethnicity; 20 imputations; seed
*         20240616). Same outcome model as Part C, pooled with Rubin's rules.
* Outcome model: survey-weighted Poisson, high-risk alcohol use on cycle,
*         age group, sex and race/ethnicity. PFQ era weight wt_pfq, DLQ era
*         weight wt_dlq.
* Design: strata SDMVSTRA, PSU SDMVPSU, singleunit(centered).
*=============================================================================

capture log close
version 17

use "$DATA", clear

log using "s10_multiple_imputation_alcohol_log.smcl", replace

*=============================================================================
* PART A: Disability x cycle interaction (complete case, both groups)
*=============================================================================
di "===== Interaction, PFQ era (cycles 1-7) ====="
svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)
quietly svy: poisson highrisk_alcohol i.disability##c.cycle_ord i.age_grp ///
    i.riagendr i.ridreth1 if inrange(cycle_ord,1,7)
lincom 1.disability#c.cycle_ord, eform

di "===== Interaction, DLQ era (cycles 8-10) ====="
svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)
quietly svy: poisson highrisk_alcohol i.disability##c.cycle_ord i.age_grp ///
    i.riagendr i.ridreth1 if inrange(cycle_ord,8,10)
lincom 1.disability#c.cycle_ord, eform

*=============================================================================
* Adults with disability from here on
*=============================================================================
keep if disability==1

*=============================================================================
* PART B: Missingness on high-risk alcohol use
*=============================================================================
di "===== Missingness, adults with disability ====="
misstable summarize highrisk_alcohol
count
count if highrisk_alcohol==.
di "Missing: " r(N)

*=============================================================================
* PART C: Complete-case trend models
*=============================================================================
di "===== Complete case: PFQ era trend ====="
svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)
quietly svy: poisson highrisk_alcohol c.cycle_ord i.age_grp i.riagendr ///
    i.ridreth1 if inrange(cycle_ord,1,7)
lincom cycle_ord, eform

di "===== Complete case: DLQ era trend ====="
svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)
quietly svy: poisson highrisk_alcohol c.cycle_ord i.age_grp i.riagendr ///
    i.ridreth1 if inrange(cycle_ord,8,10)
lincom cycle_ord, eform

*=============================================================================
* PART D: Multiple imputation
*=============================================================================
mi set wide
mi register imputed highrisk_alcohol
mi register regular cycle_ord age_grp riagendr ridreth1 sdmvpsu sdmvstra ///
    wt_pfq wt_dlq

set seed 20240616
mi impute chained (logit) highrisk_alcohol = c.cycle_ord i.age_grp ///
    i.riagendr i.ridreth1, add(20) rseed(20240616)

di "===== MI: PFQ era trend ====="
mi svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)
mi estimate, eform vartable: svy: poisson highrisk_alcohol c.cycle_ord ///
    i.age_grp i.riagendr i.ridreth1 if inrange(cycle_ord,1,7)

di "===== MI: DLQ era trend ====="
mi svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)
mi estimate, eform vartable: svy: poisson highrisk_alcohol c.cycle_ord ///
    i.age_grp i.riagendr i.ridreth1 if inrange(cycle_ord,8,10)

log close
