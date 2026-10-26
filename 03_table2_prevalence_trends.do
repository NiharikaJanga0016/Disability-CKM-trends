*=============================================================================
* 03_table2_prevalence_trends.do
* Table 2 and Supplemental Table S1: age-standardized prevalence by cycle and
* prevalence-ratio (PR) trends among US adults with disability, NHANES 1999-2018.
*
* Requires: global DATA (set in 00_config.do or in the Stata command window)
*
* Part A: age-standardized prevalence (%) by cycle, standardized to the
*         disability-specific 2017-2018 reference population
*         (20-39y 0.2305, 40-59y 0.2651, 60-79y 0.3856, 80+y 0.1187).
* Part B: trend PR per 2-year cycle, 95% CI and P value from survey-weighted
*         Poisson regression adjusted for age group, sex and race/ethnicity.
*         PFQ era = cycles 1-7 (weight wt_pfq); DLQ era = cycles 8-10 (wt_dlq).
* Design: strata SDMVSTRA, PSU SDMVPSU, singleunit(centered).
*=============================================================================

capture log close
version 17

use "$DATA", clear

log using "03_table2_prevalence_trends_log.smcl", replace

* Disability-specific 2017-2018 reference population weights
local w1 = 0.2305   // age_grp 1: 20-39
local w2 = 0.2651   // age_grp 2: 40-59
local w3 = 0.3856   // age_grp 3: 60-79
local w4 = 0.1187   // age_grp 4: 80+

local outcomes diabetes htn hyperlipid obese smoking highrisk_alcohol

*=============================================================================
* PART A: Age-standardized prevalence (%) by cycle
*=============================================================================
svyset sdmvpsu [pweight=wtmec2yr], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

tempname resA
postfile `resA' str20 outcome cycle_ord double std_prev se ll ul ///
    using table2_prevalence_by_cycle, replace

foreach outcome of local outcomes {
    forvalues c = 1/10 {

        quietly count if disability==1 & cycle_ord==`c' & `outcome'!=.
        if r(N) > 0 {

            quietly svy: poisson `outcome' i.age_grp i.riagendr i.ridreth1 ///
                if disability==1 & cycle_ord==`c'

            quietly margins, over(age_grp) post

            quietly nlcom (std_prev: `w1'*_b[1.age_grp] + `w2'*_b[2.age_grp] ///
                + `w3'*_b[3.age_grp] + `w4'*_b[4.age_grp])

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

use table2_prevalence_by_cycle, clear
foreach v in std_prev se ll ul {
    replace `v' = 100*`v'
}
format std_prev se ll ul %5.1f
sort outcome cycle_ord
di "===== Age-standardized prevalence (%), all cycles ====="
list, sepby(outcome) noobs
di "===== Cycles shown in Table 2: 1 (1999-00), 7 (2011-12), 8 (2013-14), 10 (2017-18) ====="
list if inlist(cycle_ord,1,7,8,10), sepby(outcome) noobs
save table2_prevalence_by_cycle, replace

*=============================================================================
* PART B: Trend PR per cycle, 95% CI and P value
*=============================================================================
use "$DATA", clear

tempname resB
postfile `resB' str20 outcome str3 era double PR ll ul p ///
    using table2_trend_results, replace

* PFQ era (cycles 1-7)
svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

foreach outcome of local outcomes {
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp i.riagendr ///
        i.ridreth1 if disability==1 & inrange(cycle_ord,1,7)
    quietly lincom cycle_ord, eform
    post `resB' ("`outcome'") ("PFQ") (r(estimate)) (r(lb)) (r(ub)) (r(p))
    di "PFQ `outcome': PR=" %6.3f r(estimate) "  p=" %6.3f r(p)
}

* DLQ era (cycles 8-10)
svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

foreach outcome of local outcomes {
    quietly svy: poisson `outcome' c.cycle_ord i.age_grp i.riagendr ///
        i.ridreth1 if disability==1 & inrange(cycle_ord,8,10)
    quietly lincom cycle_ord, eform
    post `resB' ("`outcome'") ("DLQ") (r(estimate)) (r(lb)) (r(ub)) (r(p))
    di "DLQ `outcome': PR=" %6.3f r(estimate) "  p=" %6.3f r(p)
}

postclose `resB'

use table2_trend_results, clear
format PR ll ul %6.3f
format p %7.3f
di "===== Trend PR per 2-year cycle ====="
list, sepby(outcome) noobs
save table2_trend_results, replace

log close
