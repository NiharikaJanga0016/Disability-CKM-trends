*=============================================================================
* s9_adjusted_trends.do
* Supplemental Table S9: fully adjusted trend estimates among adults with
* disability, NHANES 1999-2018.
*
* Requires: global DATA (set in 00_config.do or in the Stata command window)
*
* Adjusted prevalence ratio (PR) per 2-year survey cycle from survey-weighted
* Poisson regression. Covariates: age group, sex, race/ethnicity, education,
* marital status, income-to-poverty category, current smoking and BMI;
* panels S9b and S9d also adjust for health insurance coverage (HIQ011).
*
*   S9a  PFQ era, cycles 1-7 (1999-2012), no insurance covariate
*   S9b  PFQ era, cycles 4-7 (2005-2012), with insurance covariate
*        (insurance data begin in cycle 4)
*   S9c  DLQ era, cycles 8-10 (2013-2018), no insurance covariate
*   S9d  DLQ era, cycles 8-10 (2013-2018), with insurance covariate
*   S9e  PFQ era, cycles 4-7, same sample as S9b, no insurance covariate
*
* Weights: wt_pfq (PFQ era) and wt_dlq (DLQ era). N is the number of
* participants with a positive examination weight and complete data on every
* covariate in the model. Stata's e(N) (printed as n_eN) is larger because it
* also counts interview-only participants with zero examination weight.
* Outcomes: n_rf is the number of risk factors (count outcome);
* lipid_treated is statin use.
* Design: strata SDMVSTRA, PSU SDMVPSU, singleunit(centered).
*=============================================================================

capture log close
version 17

use "$DATA", clear

log using "s9_adjusted_trends_log.smcl", replace

local outcomes diabetes htn hyperlipid obese highrisk_alcohol htn_treated ///
    htn_controlled dm_treated dm_controlled lipid_treated lipid_controlled ///
    rf2plus n_rf dm_htn_hyperlip

* Panel settings: first cycle, last cycle, weight, insurance covariate (1/0),
* same-sample restriction (1/0)
local lo_a 1
local hi_a 7
local wt_a wt_pfq
local ins_a 0
local same_a 0

local lo_b 4
local hi_b 7
local wt_b wt_pfq
local ins_b 1
local same_b 0

local lo_c 8
local hi_c 10
local wt_c wt_dlq
local ins_c 0
local same_c 0

local lo_d 8
local hi_d 10
local wt_d wt_dlq
local ins_d 1
local same_d 0

local lo_e 4
local hi_e 7
local wt_e wt_pfq
local ins_e 0
local same_e 1

foreach p in a b c d e {

    svyset sdmvpsu [pweight=`wt_`p''], strata(sdmvstra) vce(linearized) ///
        singleunit(centered)

    local insvar ""
    if `ins_`p'' == 1 local insvar "insured"

    local cond "disability==1 & inrange(cycle_ord,`lo_`p'',`hi_`p'')"
    if `same_`p'' == 1 local cond "`cond' & insured!=."

    tempname res
    postfile `res' str20 outcome double PR ll ul p n n_eN ///
        using s9`p'_results, replace

    foreach outcome of local outcomes {
        quietly svy: poisson `outcome' c.cycle_ord i.age_grp i.riagendr ///
            i.ridreth1 i.educ3 i.marital3 i.ipir_cat smoking bmxbmi ///
            `insvar' if `cond'
        * Stata's e(N) also counts interview-only participants with zero
        * examination weight, who contribute nothing to the estimates.
        * N reported here = participants in the model with positive weight.
        local nobs_eN = e(N)
        quietly count if e(sample) & `wt_`p'' > 0
        local nobs = r(N)
        quietly lincom cycle_ord, eform
        post `res' ("`outcome'") (r(estimate)) (r(lb)) (r(ub)) (r(p)) ///
            (`nobs') (`nobs_eN')
    }
    postclose `res'

    preserve
    use s9`p'_results, clear
    format PR ll ul %6.3f
    format p %7.4f
    format n n_eN %9.0fc
    di _n "===== S9`p': cycles `lo_`p''-`hi_`p'', insurance covariate = `ins_`p'', same sample = `same_`p'' ====="
    list, noobs
    restore
}

log close
