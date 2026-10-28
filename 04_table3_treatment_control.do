*=============================================================================
* 04_table3_treatment_control.do
* Table 3 (treatment and control rows) and Supplemental Table S2: treatment and
* control of diabetes, hypertension and hyperlipidemia among US adults with
* disability, NHANES 1999-2018.
*
* Requires: global DATA (set in 00_config.do or in the Stata command window)
*
* Part A: age-standardized % by cycle for the six single-condition outcomes,
*         standardized to the disability-specific 2017-2018 reference
*         (20-59y 0.4956, 60-79y 0.3856, 80+y 0.1187), among adults with the
*         relevant condition.
* Part B: trend PR per 2-year cycle (95% CI, P) for the same six outcomes,
*         adjusted for 3-group age, sex and race/ethnicity.
*         PFQ era = cycles 1-7 (weight wt_pfq); DLQ era = cycles 8-10 (wt_dlq).
* Part C: joint control outcomes. Crude survey-weighted % by cycle; trend PR for
*         joint DM+HTN control (adjusted model) and joint all-3 control (crude
*         model, unadjusted).
* Part D: unweighted N by cycle for every denominator in Table 3 and Table S2.
* Design: strata SDMVSTRA, PSU SDMVPSU, singleunit(centered).
*=============================================================================

capture log close
version 17

use "$DATA", clear

log using "04_table3_treatment_control_log.smcl", replace

* Variables used throughout
gen joint_ctrl_all3_c = (dm_controlled==1 & htn_controlled==1 & ///
    lipid_controlled==1) if dm_htn_hyperlip==1

gen age_grp3 = 1 if inlist(age_grp,1,2)
replace age_grp3 = 2 if age_grp==3
replace age_grp3 = 3 if age_grp==4
label define age_grp3lbl 1 "20-59" 2 "60-79" 3 "80+"
label values age_grp3 age_grp3lbl

tempfile base
save `base'

* Disability-specific 2017-2018 reference weights (3 age groups)
local w1 = 0.4956
local w2 = 0.3856
local w3 = 0.1187

* outcome and its denominator condition
local spec1 "htn_treated htn==1"
local spec2 "htn_controlled htn==1"
local spec3 "dm_treated diabetes==1"
local spec4 "dm_controlled diabetes==1"
local spec5 "lipid_treated hyperlipid==1"
local spec6 "lipid_controlled hyperlipid==1"

*=============================================================================
* PART A: Age-standardized % by cycle, single-condition outcomes
*=============================================================================
svyset sdmvpsu [pweight=wtmec2yr], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

tempname resA
postfile `resA' str24 outcome cycle_ord double std_prev se ll ul ///
    using table3_by_cycle, replace

forvalues s = 1/6 {
    tokenize `spec`s''
    local outcome `1'
    local denomcond `2'

    forvalues c = 1/10 {
        quietly count if disability==1 & cycle_ord==`c' & `denomcond' & `outcome'!=.
        if r(N) > 0 {
            quietly svy: poisson `outcome' i.age_grp3 i.riagendr i.ridreth1 ///
                if disability==1 & `denomcond' & cycle_ord==`c'

            quietly margins, over(age_grp3) post

            quietly nlcom (std_prev: `w1'*_b[1.age_grp3] + `w2'*_b[2.age_grp3] ///
                + `w3'*_b[3.age_grp3])

            matrix R = r(b)
            matrix V = r(V)
            local est = R[1,1]
            local se  = sqrt(V[1,1])
            local lo  = `est' - 1.96*`se'
            local hi  = `est' + 1.96*`se'

            post `resA' ("`outcome'") (`c') (`est') (`se') (`lo') (`hi')
            di "`outcome', cycle `c': " %5.1f 100*`est' "%"
        }
    }
}
postclose `resA'

use table3_by_cycle, clear
gen pct    = 100*std_prev
gen pct_ll = 100*ll
gen pct_ul = 100*ul
format pct pct_ll pct_ul %5.1f
di "===== Age-standardized % by cycle (Table 3 shows cycles 1, 7, 8, 10) ====="
list outcome cycle_ord pct pct_ll pct_ul, sepby(outcome) noobs
save table3_by_cycle, replace

*=============================================================================
* PART B: Trend PR per 2-year cycle, single-condition outcomes
*=============================================================================
use `base', clear

tempname resB
postfile `resB' str24 outcome str3 era double PR ll ul p ///
    using table3_trend_results, replace

* PFQ era
svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

forvalues s = 1/6 {
    tokenize `spec`s''
    local outcome `1'
    local denomcond `2'
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp3 i.riagendr i.ridreth1 ///
        if disability==1 & `denomcond' & inrange(cycle_ord,1,7)
    quietly lincom cycle_ord, eform
    post `resB' ("`outcome'") ("PFQ") (r(estimate)) (r(lb)) (r(ub)) (r(p))
    di "PFQ `outcome': PR=" %6.3f r(estimate) "  p=" %6.3f r(p)
}

* DLQ era
svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

forvalues s = 1/6 {
    tokenize `spec`s''
    local outcome `1'
    local denomcond `2'
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp3 i.riagendr i.ridreth1 ///
        if disability==1 & `denomcond' & inrange(cycle_ord,8,10)
    quietly lincom cycle_ord, eform
    post `resB' ("`outcome'") ("DLQ") (r(estimate)) (r(lb)) (r(ub)) (r(p))
    di "DLQ `outcome': PR=" %6.3f r(estimate) "  p=" %6.3f r(p)
}
postclose `resB'

use table3_trend_results, clear
format PR ll ul %6.3f
format p %7.3f
di "===== Trend PR per 2-year cycle, single-condition outcomes ====="
list, sepby(outcome) noobs
save table3_trend_results, replace

*=============================================================================
* PART C: Joint control outcomes
*=============================================================================
use `base', clear

svyset sdmvpsu [pweight=wtmec2yr], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

di "===== Joint DM+HTN control: crude % by cycle ====="
forvalues c = 1/10 {
    quietly count if disability==1 & cycle_ord==`c' & dm_htn==1 & joint_ctrl_dm_htn!=.
    if r(N) > 0 {
        quietly svy: mean joint_ctrl_dm_htn if disability==1 & dm_htn==1 & cycle_ord==`c'
        matrix T = r(table)
        di "Cycle `c': " %5.1f 100*T[1,1] "% (" %5.1f 100*T[5,1] "-" %5.1f 100*T[6,1] ")"
    }
}

di "===== Joint all-3 control: crude % by cycle ====="
forvalues c = 1/10 {
    quietly count if disability==1 & cycle_ord==`c' & dm_htn_hyperlip==1 & joint_ctrl_all3_c!=.
    if r(N) > 0 {
        quietly svy: mean joint_ctrl_all3_c if disability==1 & dm_htn_hyperlip==1 & cycle_ord==`c'
        matrix T = r(table)
        di "Cycle `c': " %5.1f 100*T[1,1] "% (" %5.1f 100*T[5,1] "-" %5.1f 100*T[6,1] ")"
    }
}

di "===== Joint DM+HTN control: trend PR (adjusted) ====="
svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) singleunit(centered)
quietly svy: poisson joint_ctrl_dm_htn c.cycle_ord i.age_grp3 i.riagendr i.ridreth1 ///
    if disability==1 & dm_htn==1 & inrange(cycle_ord,1,7)
quietly lincom cycle_ord, eform
di "PFQ: PR=" %6.3f r(estimate) " (" %6.3f r(lb) "-" %6.3f r(ub) ")  p=" %6.3f r(p)

svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) singleunit(centered)
quietly svy: poisson joint_ctrl_dm_htn c.cycle_ord i.age_grp3 i.riagendr i.ridreth1 ///
    if disability==1 & dm_htn==1 & inrange(cycle_ord,8,10)
quietly lincom cycle_ord, eform
di "DLQ: PR=" %6.3f r(estimate) " (" %6.3f r(lb) "-" %6.3f r(ub) ")  p=" %6.3f r(p)

di "===== Joint all-3 control: trend PR (crude model, unadjusted) ====="
svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) singleunit(centered)
quietly svy: poisson joint_ctrl_all3_c c.cycle_ord ///
    if disability==1 & dm_htn_hyperlip==1 & inrange(cycle_ord,1,7)
quietly lincom cycle_ord, eform
di "PFQ: PR=" %6.3f r(estimate) " (" %6.3f r(lb) "-" %6.3f r(ub) ")  p=" %6.3f r(p)

svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) singleunit(centered)
quietly svy: poisson joint_ctrl_all3_c c.cycle_ord ///
    if disability==1 & dm_htn_hyperlip==1 & inrange(cycle_ord,8,10)
quietly lincom cycle_ord, eform
di "DLQ: PR=" %6.3f r(estimate) " (" %6.3f r(lb) "-" %6.3f r(ub) ")  p=" %6.3f r(p)

*=============================================================================
* PART D: Unweighted N by cycle (adults with disability), Table 3 / S2
* Columns are cycles 1 to 10.
*=============================================================================
use `base', clear

di "===== Unweighted N by cycle (cycles 1-10) ====="

local d1 "Joint DM+HTN control (dm_htn==1)"
local c1 "dm_htn==1"
local d2 "Joint all-3 control (dm_htn_hyperlip==1)"
local c2 "dm_htn_hyperlip==1"
local d3 "Glycemic control (diabetes==1 & lbxgh nonmissing)"
local c3 "diabetes==1 & lbxgh!=."
local d4 "Diabetes treatment (diabetes==1)"
local c4 "diabetes==1"
local d5 "HTN treatment and BP control (htn==1)"
local c5 "htn==1"
local d6 "Statin use and lipid control (hyperlipid==1)"
local c6 "hyperlipid==1"

forvalues k = 1/6 {
    local row ""
    forvalues c = 1/10 {
        quietly count if disability==1 & cycle_ord==`c' & `c`k''
        local row `row' `r(N)'
    }
    di "`d`k'': `row'"
}

log close
