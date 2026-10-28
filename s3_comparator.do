*=============================================================================
* s3_comparator.do
* Supplemental Table S3 and the comparator rows of Table 3: adults with versus
* without disability, NHANES 1999-2018.
*
* Requires: global DATA (set in 00_config.do or in the Stata command window)
*
* Part A: age-standardized prevalence (%) by cycle for adults with disability
*         and adults without disability, both standardized to the 2017-2018
*         US all-adults reference population
*         (20-39y 0.3555, 40-59y 0.3497, 60-79y 0.2479, 80+y 0.0468).
*         Outcomes: diabetes, hypertension, hyperlipidemia, obesity, smoking.
* Part B: disability x cycle interaction PR (95% CI, P) from a pooled
*         Poisson model including both groups, adjusted for age group, sex
*         and race/ethnicity. PFQ era = cycles 1-7 (weight wt_pfq);
*         DLQ era = cycles 8-10 (weight wt_dlq).
* Design: strata SDMVSTRA, PSU SDMVPSU, singleunit(centered).
*=============================================================================

capture log close
version 17

use "$DATA", clear

log using "s3_comparator_log.smcl", replace

* 2017-2018 US all-adults reference population weights
local w1 = 0.3555   // age_grp 1: 20-39
local w2 = 0.3497   // age_grp 2: 40-59
local w3 = 0.2479   // age_grp 3: 60-79
local w4 = 0.0468   // age_grp 4: 80+

local outcomes diabetes htn hyperlipid obese smoking

*=============================================================================
* PART A: Age-standardized prevalence (%) by cycle and disability group
*=============================================================================
svyset sdmvpsu [pweight=wtmec2yr], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

tempname resA
postfile `resA' str20 outcome byte disability cycle_ord double std_prev se ll ul ///
    using s3_by_cycle, replace

foreach outcome of local outcomes {
    foreach g in 1 0 {
        forvalues c = 1/10 {

            quietly count if disability==`g' & cycle_ord==`c'
            if r(N) > 0 {

                quietly svy: poisson `outcome' i.age_grp i.riagendr i.ridreth1 ///
                    if disability==`g' & cycle_ord==`c'

                quietly margins, over(age_grp) post

                quietly nlcom (std_prev: `w1'*_b[1.age_grp] + `w2'*_b[2.age_grp] ///
                    + `w3'*_b[3.age_grp] + `w4'*_b[4.age_grp])

                matrix R = r(b)
                matrix V = r(V)
                local est = R[1,1]
                local se  = sqrt(V[1,1])
                local lo  = `est' - 1.96*`se'
                local hi  = `est' + 1.96*`se'

                post `resA' ("`outcome'") (`g') (`c') (`est') (`se') (`lo') (`hi')
                di "`outcome', disability=`g', cycle `c': " %5.1f 100*`est' "%"
            }
        }
    }
}
postclose `resA'

use s3_by_cycle, clear
gen pct    = 100*std_prev
gen pct_ll = 100*ll
gen pct_ul = 100*ul
format pct pct_ll pct_ul %5.1f
sort outcome disability cycle_ord
di "===== Age-standardized prevalence (%), disability=1 with, disability=0 without ====="
list outcome disability cycle_ord pct pct_ll pct_ul, sepby(outcome disability) noobs
save s3_by_cycle, replace

*=============================================================================
* PART B: Disability x cycle interaction PR (95% CI, P)
*=============================================================================
use "$DATA", clear

tempname resB
postfile `resB' str20 outcome str3 era double int_PR ll ul p ///
    using s3_interaction_results, replace

* PFQ era (cycles 1-7)
svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

foreach outcome of local outcomes {
    quietly svy: poisson `outcome' i.disability##c.cycle_ord ///
        i.age_grp i.riagendr i.ridreth1 if inrange(cycle_ord,1,7)
    quietly lincom 1.disability#c.cycle_ord, eform
    post `resB' ("`outcome'") ("PFQ") (r(estimate)) (r(lb)) (r(ub)) (r(p))
    di "PFQ `outcome': interaction PR=" %6.3f r(estimate) "  p=" %7.4f r(p)
}

* DLQ era (cycles 8-10)
svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) vce(linearized) ///
    singleunit(centered)

foreach outcome of local outcomes {
    quietly svy: poisson `outcome' i.disability##c.cycle_ord ///
        i.age_grp i.riagendr i.ridreth1 if inrange(cycle_ord,8,10)
    quietly lincom 1.disability#c.cycle_ord, eform
    post `resB' ("`outcome'") ("DLQ") (r(estimate)) (r(lb)) (r(ub)) (r(p))
    di "DLQ `outcome': interaction PR=" %6.3f r(estimate) "  p=" %7.4f r(p)
}
postclose `resB'

use s3_interaction_results, clear
format int_PR ll ul %6.3f
format p %8.4f
di "===== Disability x cycle interaction PR ====="
list, sepby(outcome) noobs
save s3_interaction_results, replace

log close
