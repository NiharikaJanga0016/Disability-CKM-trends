*=============================================================================
* s5_fdr_correction.do
* Supplemental Table S5: Benjamini-Hochberg false discovery rate (FDR)
* correction across the 44 primary hypothesis tests.
*
* Requires: global DATA (set in 00_config.do or in the Stata command window)
*
* Each test is the trend (or interaction) term from a survey-weighted Poisson
* model, run separately for the PFQ era (cycles 1-7, weight wt_pfq) and the
* DLQ era (cycles 8-10, weight wt_dlq). Each family uses the same model as the
* table it corresponds to:
*   Domain 1  Prevalence trends, 12 tests          4-group age; Table 2, S1
*   Domain 2  Treatment and control, 16 tests      3-group age; Table 3, S2
*             (joint all-3 control uses the crude model, no covariates)
*   Domain 3  Comorbidity burden, 6 tests          4-group age; Table S4
*   Domain 4  Disability x cycle interaction, 10   4-group age; Table S3
* Covariates in every adjusted model: age group, sex, race/ethnicity.
* Correction: BH step-up across all 44 tests, adjusted p = raw p x 44 / rank,
* with monotonicity enforced. FDR threshold 5%.
* Design: strata SDMVSTRA, PSU SDMVPSU, singleunit(centered).
*=============================================================================

capture log close
version 17

use "$DATA", clear

log using "s5_fdr_correction_log.smcl", replace

* 3-group disability-specific age variable (Domain 2)
gen age_grp3 = 1 if inlist(age_grp,1,2)
replace age_grp3 = 2 if age_grp==3
replace age_grp3 = 3 if age_grp==4

* Joint control of all three conditions
gen joint_ctrl_all3_c = (dm_controlled==1 & htn_controlled==1 & ///
    lipid_controlled==1) if dm_htn_hyperlip==1

* Three or more concurrent risk factors
capture confirm variable rf3plus
if _rc {
    gen rf3plus = (n_rf>=3) if n_rf!=.
}

tempname results
postfile `results' str10 domain str30 outcome str3 era double raw_p ///
    double PR double ll double ul ///
    using s5_raw_results, replace

*=============================================================================
* DOMAIN 1: Prevalence trends (12 tests)
*=============================================================================
local prev_outcomes diabetes htn hyperlipid obese smoking highrisk_alcohol

foreach outcome of local prev_outcomes {
    svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
        singleunit(centered)
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp i.riagendr ///
        i.ridreth1 if disability==1 & inrange(cycle_ord,1,7)
    quietly lincom cycle_ord, eform
    post `results' ("Domain1") ("`outcome'") ("PFQ") (r(p)) (r(estimate)) (r(lb)) (r(ub))

    svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
        singleunit(centered)
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp i.riagendr ///
        i.ridreth1 if disability==1 & inrange(cycle_ord,8,10)
    quietly lincom cycle_ord, eform
    post `results' ("Domain1") ("`outcome'") ("DLQ") (r(p)) (r(estimate)) (r(lb)) (r(ub))

    di "Domain 1: `outcome' done"
}

*=============================================================================
* DOMAIN 2: Treatment and control (16 tests)
*=============================================================================
local tc_outcomes_adjusted htn_treated htn_controlled dm_treated ///
    dm_controlled lipid_treated lipid_controlled joint_ctrl_dm_htn

foreach outcome of local tc_outcomes_adjusted {
    svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
        singleunit(centered)
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp3 i.riagendr ///
        i.ridreth1 if disability==1 & inrange(cycle_ord,1,7)
    quietly lincom cycle_ord, eform
    post `results' ("Domain2") ("`outcome'") ("PFQ") (r(p)) (r(estimate)) (r(lb)) (r(ub))

    svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
        singleunit(centered)
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp3 i.riagendr ///
        i.ridreth1 if disability==1 & inrange(cycle_ord,8,10)
    quietly lincom cycle_ord, eform
    post `results' ("Domain2") ("`outcome'") ("DLQ") (r(p)) (r(estimate)) (r(lb)) (r(ub))

    di "Domain 2: `outcome' done"
}

* Joint all-3 control: crude model (no covariates)
svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)
quietly svy: poisson joint_ctrl_all3_c c.cycle_ord ///
    if disability==1 & dm_htn_hyperlip==1 & inrange(cycle_ord,1,7)
quietly lincom cycle_ord, eform
post `results' ("Domain2") ("joint_ctrl_all3") ("PFQ") (r(p)) (r(estimate)) (r(lb)) (r(ub))

svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)
quietly svy: poisson joint_ctrl_all3_c c.cycle_ord ///
    if disability==1 & dm_htn_hyperlip==1 & inrange(cycle_ord,8,10)
quietly lincom cycle_ord, eform
post `results' ("Domain2") ("joint_ctrl_all3") ("DLQ") (r(p)) (r(estimate)) (r(lb)) (r(ub))

di "Domain 2: joint_ctrl_all3 done"

*=============================================================================
* DOMAIN 3: Comorbidity burden (6 tests)
*=============================================================================
local comorbid_outcomes rf2plus rf3plus dm_htn_hyperlip

foreach outcome of local comorbid_outcomes {
    svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
        singleunit(centered)
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp i.riagendr ///
        i.ridreth1 if disability==1 & inrange(cycle_ord,1,7)
    quietly lincom cycle_ord, eform
    post `results' ("Domain3") ("`outcome'") ("PFQ") (r(p)) (r(estimate)) (r(lb)) (r(ub))

    svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
        singleunit(centered)
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp i.riagendr ///
        i.ridreth1 if disability==1 & inrange(cycle_ord,8,10)
    quietly lincom cycle_ord, eform
    post `results' ("Domain3") ("`outcome'") ("DLQ") (r(p)) (r(estimate)) (r(lb)) (r(ub))

    di "Domain 3: `outcome' done"
}

*=============================================================================
* DOMAIN 4: Comparator analysis, disability x cycle interaction (10 tests)
*=============================================================================
local comp_outcomes diabetes htn hyperlipid obese smoking

foreach outcome of local comp_outcomes {
    svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
        singleunit(centered)
    quietly svy: poisson `outcome' i.disability##c.cycle_ord ///
        i.age_grp i.riagendr i.ridreth1 if inrange(cycle_ord,1,7)
    quietly lincom 1.disability#c.cycle_ord, eform
    post `results' ("Domain4") ("`outcome'") ("PFQ") (r(p)) (r(estimate)) (r(lb)) (r(ub))

    svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
        singleunit(centered)
    quietly svy: poisson `outcome' i.disability##c.cycle_ord ///
        i.age_grp i.riagendr i.ridreth1 if inrange(cycle_ord,8,10)
    quietly lincom 1.disability#c.cycle_ord, eform
    post `results' ("Domain4") ("`outcome'") ("DLQ") (r(p)) (r(estimate)) (r(lb)) (r(ub))

    di "Domain 4: `outcome' done"
}

postclose `results'

*=============================================================================
* Benjamini-Hochberg step-up correction across all 44 tests
*=============================================================================
use s5_raw_results, clear
count
di "Total tests in primary family: " r(N)

gsort raw_p
gen rank = _n
local m = _N
gen bh_adj_p = raw_p * `m' / rank

* enforce monotonicity: running minimum from the largest rank downwards
gsort -rank
gen bh_adj_p_mono = bh_adj_p
replace bh_adj_p_mono = min(bh_adj_p_mono, bh_adj_p_mono[_n-1]) if _n>1
replace bh_adj_p_mono = min(bh_adj_p_mono, 1)
gsort rank

gen fdr_sig    = (bh_adj_p_mono < 0.05)
gen borderline = (raw_p < 0.05 & bh_adj_p_mono >= 0.05)

* PR and confidence limits are printed to 4 decimals; the tables report 3
format PR ll ul %6.4f
format raw_p bh_adj_p_mono %10.4g
di "===== All 44 tests, ordered by raw p (rank) ====="
list rank domain outcome era PR ll ul raw_p bh_adj_p_mono fdr_sig borderline, noobs

count if fdr_sig==1
di "Tests surviving FDR correction: " r(N)
count if borderline==1
di "Borderline tests (raw p<0.05, BH-adjusted p>=0.05): " r(N)
count if fdr_sig==0 & borderline==0
di "Not significant: " r(N)

save s5_fdr_results, replace

log close
