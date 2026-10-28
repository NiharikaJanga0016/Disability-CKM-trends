*=============================================================================
* s6_sensitivity_analyses.do
* Supplemental Table S6: sensitivity analyses, NHANES 1999-2018.
*
* Requires: global DATA (set in 00_config.do or in the Stata command window)
*
* S6a  BP control, DLQ era, excluding the 2017-2018 cycle (new BP device).
*      3-group age covariate, as in Table 3 / Supplemental Table S2.
* S6b  Hypertension defined as SBP>=130 or DBP>=80 (ACC/AHA) instead of
*      the primary JNC 7 definition (variable htn_130).
* S6c  Hyperlipidemia defined as total cholesterol >=200 mg/dL or diagnosis
*      instead of the primary non-HDL definition (variable hyperlipid_tc).
* S6d  Mobility-restricted disability (dlq_mobility==1), DLQ era only.
*
* S6b and S6c, adults with disability:
*   survey-weighted prevalence (%) by cycle (crude, not age-standardized), and
*   trend PR per 2-year cycle from survey-weighted Poisson regression
*   adjusted for age group, sex and race/ethnicity.
*   PFQ era = cycles 1-7 (weight wt_pfq); DLQ era = cycles 8-10 (wt_dlq).
* Design: strata SDMVSTRA, PSU SDMVPSU, singleunit(centered).
*=============================================================================

capture log close
version 17

use "$DATA", clear

log using "s6_sensitivity_analyses_log.smcl", replace

* 3-group disability-specific age variable (S6a)
gen age_grp3 = 1 if inlist(age_grp,1,2)
replace age_grp3 = 2 if age_grp==3
replace age_grp3 = 3 if age_grp==4

tempfile base
save `base'

*=============================================================================
* S6a: BP control, excluding the 2017-2018 device-change cycle
*=============================================================================
svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

di "===== S6a primary: all 3 DLQ cycles (2013-2018) ====="
quietly svy: poisson htn_controlled c.cycle_ord i.age_grp3 i.riagendr ///
    i.ridreth1 if disability==1 & htn==1 & inrange(cycle_ord,8,10)
lincom cycle_ord, eform

di "===== S6a sensitivity: excluding 2017-2018 (cycles 8-9) ====="
quietly svy: poisson htn_controlled c.cycle_ord i.age_grp3 i.riagendr ///
    i.ridreth1 if disability==1 & htn==1 & inrange(cycle_ord,8,9)
lincom cycle_ord, eform

*=============================================================================
* S6b and S6c: alternative definitions, adults with disability
*=============================================================================
local alt_outcomes htn_130 hyperlipid_tc

* Part 1: survey-weighted prevalence (%) by cycle (crude)
svyset sdmvpsu [pweight=wtmec2yr], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

tempname resA
postfile `resA' str20 outcome cycle_ord double prev se ll ul ///
    using s6_by_cycle, replace

foreach outcome of local alt_outcomes {
    forvalues c = 1/10 {
        quietly count if disability==1 & cycle_ord==`c' & `outcome'!=.
        if r(N) > 0 {
            quietly svy: mean `outcome' if disability==1 & cycle_ord==`c'
            matrix R = r(table)
            post `resA' ("`outcome'") (`c') (R[1,1]) (R[2,1]) (R[5,1]) (R[6,1])
        }
    }
}
postclose `resA'

use s6_by_cycle, clear
gen pct    = 100*prev
gen pct_ll = 100*ll
gen pct_ul = 100*ul
format pct pct_ll pct_ul %5.1f
sort outcome cycle_ord
di "===== S6b (htn_130) and S6c (hyperlipid_tc): prevalence (%) by cycle ====="
list outcome cycle_ord pct pct_ll pct_ul, sepby(outcome) noobs
save s6_by_cycle, replace

* Part 2: trend PR per 2-year cycle
use `base', clear

tempname resB
postfile `resB' str20 outcome str3 era double PR ll ul p ///
    using s6_trend_results, replace

svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)
foreach outcome of local alt_outcomes {
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp i.riagendr i.ridreth1 ///
        if disability==1 & inrange(cycle_ord,1,7)
    quietly lincom cycle_ord, eform
    post `resB' ("`outcome'") ("PFQ") (r(estimate)) (r(lb)) (r(ub)) (r(p))
}

svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)
foreach outcome of local alt_outcomes {
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp i.riagendr i.ridreth1 ///
        if disability==1 & inrange(cycle_ord,8,10)
    quietly lincom cycle_ord, eform
    post `resB' ("`outcome'") ("DLQ") (r(estimate)) (r(lb)) (r(ub)) (r(p))
}
postclose `resB'

use s6_trend_results, clear
format PR ll ul %6.3f
format p %8.4f
di "===== S6b and S6c: trend PR per 2-year cycle ====="
list, sepby(outcome) noobs
save s6_trend_results, replace

*=============================================================================
* S6d: mobility-restricted disability, DLQ era only
*=============================================================================
use `base', clear

svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

foreach outcome in diabetes htn hyperlipid obese smoking highrisk_alcohol {
    di "===== S6d: `outcome', mobility-restricted, DLQ era ====="
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp i.riagendr ///
        i.ridreth1 if dlq_mobility==1 & inrange(cycle_ord,8,10)
    lincom cycle_ord, eform
}

log close
