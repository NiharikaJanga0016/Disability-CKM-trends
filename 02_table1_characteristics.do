*=============================================================================
* 02_table1_characteristics.do
* Table 1: baseline characteristics of US adults with and without disability,
* NHANES 1999-2018 (pooled). Percentages and means with 95% CIs, P values.
*
* Requires: global DATA (set in 00_config.do or in the Stata command window)
* Weights: WTMEC2YR. Design: strata SDMVSTRA, PSU SDMVPSU, singleunit(centered).
* Insurance (HIQ011): cycles 4-10 (2005-2018); no/refused/dont know = uninsured.
*=============================================================================

capture log close
version 17

use "$DATA", clear

log using "02_table1_characteristics_log.smcl", replace

gen female = (riagendr==2) if riagendr!=.

svyset sdmvpsu [pweight=wtmec2yr], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

*-----------------------------------------------------------------------------
* A. Percentages / means with 95% CI, by disability status
*-----------------------------------------------------------------------------
di "===== Mean age ====="
svy: mean ridageyr, over(disability)

di "===== Female % ====="
svy: proportion female, over(disability)

di "===== Age group % ====="
svy: proportion age_grp, over(disability)

di "===== Race/ethnicity % ====="
svy: proportion ridreth1, over(disability)

di "===== Income-to-poverty tier % ====="
svy: proportion ipir_cat, over(disability)

di "===== Insurance coverage %, cycles 4-10 ====="
gen ins_c410 = (insured!=. & inrange(cycle_ord,4,10))
svy, subpop(ins_c410): proportion insured, over(disability)

foreach outcome in diabetes htn hyperlipid obese smoking highrisk_alcohol {
    quietly gen `outcome'_nonmiss = (`outcome'!=.)
    di "===== `outcome' %, pooled ====="
    svy, subpop(`outcome'_nonmiss): proportion `outcome', over(disability)
}

*-----------------------------------------------------------------------------
* B. P values (design-based tests)
*-----------------------------------------------------------------------------
di "===== P values ====="
di "--- Age (linear regression) ---"
svy: regress ridageyr disability
di "--- Female (logistic) ---"
svy: logistic female disability
di "--- Diabetes (logistic) ---"
svy: logistic diabetes disability
di "--- Hypertension (logistic) ---"
svy: logistic htn disability
di "--- Hyperlipidemia (logistic) ---"
svy: logistic hyperlipid disability
di "--- Obesity (logistic) ---"
svy: logistic obese disability
di "--- Current smoking (logistic) ---"
svy: logistic smoking disability
di "--- High-risk alcohol use (logistic) ---"
svy: logistic highrisk_alcohol disability
di "--- Race/ethnicity (Rao-Scott chi-square) ---"
svy: tab ridreth1 disability, pearson
di "--- Income-to-poverty tier (Rao-Scott chi-square) ---"
svy: tab ipir_cat disability, pearson
di "--- Age group (Rao-Scott chi-square) ---"
svy: tab age_grp disability, pearson
di "--- Insurance coverage (logistic, cycles 4-10) ---"
svy, subpop(ins_c410): logistic insured disability

*-----------------------------------------------------------------------------
* C. Unweighted N for the insurance row
*-----------------------------------------------------------------------------
count if ins_c410==1
di "N with insurance data, cycles 4-10: " r(N)

log close
