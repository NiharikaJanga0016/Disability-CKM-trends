*=============================================================================
* s8_evalues.do
* Supplemental Table S8 and Supplemental Figure S3: E-values (VanderWeele and
* Ding, 2017) for the 15 findings that survived FDR correction (Table S5).
*
* Requires: global DATA (set in 00_config.do or in the Stata command window)
*
* Step 1: re-estimate the PR and 95% CI for each of the 15 findings, using the
*         same model as in Supplemental Table S5:
*         - Domain 2 (treatment and control): 3-group age covariate (age_grp3),
*           sex, race/ethnicity; joint all-3 control is the crude model
*           (cycle only, no covariates).
*         - Domains 1, 3 and 4: 4-group age covariate (age_grp), sex,
*           race/ethnicity; Domain 4 is the disability x cycle interaction.
*         PFQ era = cycles 1-7 (weight wt_pfq); DLQ era = cycles 8-10 (wt_dlq).
* Step 2: E-value for the point estimate and for the confidence limit closest
*         to the null. For PR < 1 the PR is inverted first (RR = 1/PR).
*         E-value = RR + sqrt(RR x (RR - 1)); equals 1 if the CI includes 1.
* Design: strata SDMVSTRA, PSU SDMVPSU, singleunit(centered).
*=============================================================================

capture log close
version 17

use "$DATA", clear

log using "s8_evalues_log.smcl", replace

* 3-group disability-specific age variable (Domain 2)
gen age_grp3 = 1 if inlist(age_grp,1,2)
replace age_grp3 = 2 if age_grp==3
replace age_grp3 = 3 if age_grp==4

* Joint control of all three conditions
gen joint_ctrl_all3_c = (dm_controlled==1 & htn_controlled==1 & ///
    lipid_controlled==1) if dm_htn_hyperlip==1

capture confirm variable rf3plus
if _rc {
    gen rf3plus = (n_rf>=3) if n_rf!=.
}

tempname res
postfile `res' str40 finding double PR ll ul using s8_pr_ci_results, replace

*=============================================================================
* PFQ era (cycles 1-7)
*=============================================================================
svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

* Domain 2: 3-group age covariate
foreach v in lipid_treated lipid_controlled htn_controlled htn_treated ///
    joint_ctrl_dm_htn {
    quietly svy: poisson `v' c.cycle_ord i.age_grp3 i.riagendr i.ridreth1 ///
        if disability==1 & inrange(cycle_ord,1,7)
    quietly lincom cycle_ord, eform
    post `res' ("`v'_PFQ") (r(estimate)) (r(lb)) (r(ub))
}

* Joint all-3 control: crude model
quietly svy: poisson joint_ctrl_all3_c c.cycle_ord ///
    if disability==1 & dm_htn_hyperlip==1 & inrange(cycle_ord,1,7)
quietly lincom cycle_ord, eform
post `res' ("joint_ctrl_all3_PFQ") (r(estimate)) (r(lb)) (r(ub))

* Domains 1 and 3: 4-group age covariate
foreach v in dm_htn_hyperlip rf3plus diabetes obese hyperlipid {
    quietly svy: poisson `v' c.cycle_ord i.age_grp i.riagendr i.ridreth1 ///
        if disability==1 & inrange(cycle_ord,1,7)
    quietly lincom cycle_ord, eform
    post `res' ("`v'_PFQ") (r(estimate)) (r(lb)) (r(ub))
}

* Domain 4: smoking disparity (disability x cycle interaction)
quietly svy: poisson smoking i.disability##c.cycle_ord i.age_grp ///
    i.riagendr i.ridreth1 if inrange(cycle_ord,1,7)
quietly lincom 1.disability#c.cycle_ord, eform
post `res' ("smoking_disparity_PFQ") (r(estimate)) (r(lb)) (r(ub))

*=============================================================================
* DLQ era (cycles 8-10)
*=============================================================================
svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

quietly svy: poisson hyperlipid c.cycle_ord i.age_grp i.riagendr i.ridreth1 ///
    if disability==1 & inrange(cycle_ord,8,10)
quietly lincom cycle_ord, eform
post `res' ("hyperlipid_DLQ") (r(estimate)) (r(lb)) (r(ub))

quietly svy: poisson htn_controlled c.cycle_ord i.age_grp3 i.riagendr ///
    i.ridreth1 if disability==1 & inrange(cycle_ord,8,10)
quietly lincom cycle_ord, eform
post `res' ("htn_controlled_DLQ") (r(estimate)) (r(lb)) (r(ub))

quietly svy: poisson smoking i.disability##c.cycle_ord i.age_grp ///
    i.riagendr i.ridreth1 if inrange(cycle_ord,8,10)
quietly lincom 1.disability#c.cycle_ord, eform
post `res' ("smoking_disparity_DLQ") (r(estimate)) (r(lb)) (r(ub))

postclose `res'

*=============================================================================
* E-values
*=============================================================================
use s8_pr_ci_results, clear
count
di "Total findings: " r(N)

gen rr_point = PR
replace rr_point = 1/PR if PR<1
gen evalue_point = rr_point + sqrt(rr_point*(rr_point-1))

gen evalue_ci = 1
replace evalue_ci = ll + sqrt(ll*(ll-1)) if PR>=1 & ll>1
replace evalue_ci = (1/ul) + sqrt((1/ul)*((1/ul)-1)) if PR<1 & ul<1

gsort -evalue_point
* PR and confidence limits are printed to 4 decimals; the tables report 3
format PR ll ul %6.4f
format evalue_point evalue_ci %6.3f

di "===== E-values, ordered by point-estimate E-value ====="
list finding PR ll ul evalue_point evalue_ci, noobs

save s8_evalues_results, replace

log close
