*=============================================================================
* s7_segmented_trends.do
* Supplemental Table S7: segmented trend analysis, adults with disability,
* NHANES 1999-2018.
*
* Requires: global DATA (set in 00_config.do or in the Stata command window)
*
* Step 1: age-standardized prevalence (%) and its standard error for each
*         outcome and cycle, from a survey-weighted Poisson model with age
*         group as the only covariate, standardized to the disability-specific
*         2017-2018 reference population
*         (20-39y 0.2305, 40-59y 0.2651, 60-79y 0.3856, 80+y 0.1187).
*         Cycles 1-7 use weight wt_pfq; cycles 8-10 use wt_dlq.
* Step 2: inverse-variance-weighted linear regression (vwls) of the cycle-level
*         prevalence on cycle number, fit three ways per outcome:
*         overall (cycles 1-10), PFQ era (cycles 1-7), DLQ era (cycles 8-10).
*         The slope is the average percentage-point change (APC) per
*         2-year cycle.
* Design: strata SDMVSTRA, PSU SDMVPSU, singleunit(centered).
*=============================================================================

capture log close
version 17

use "$DATA", clear

log using "s7_segmented_trends_log.smcl", replace

local outcomes diabetes htn hyperlipid obese smoking

*=============================================================================
* Step 1: cycle-level age-standardized prevalence and SE
*=============================================================================
tempname results
matrix `results' = J(50,4,.)   // columns: outcome number, cycle, pct, se
local row = 0
local oi = 0

foreach outcome of local outcomes {
    local ++oi
    forvalues c = 1/10 {
        local ++row
        matrix `results'[`row',1] = `oi'
        matrix `results'[`row',2] = `c'

        quietly count if disability==1 & cycle_ord==`c' & `outcome'!=.
        if r(N) > 0 {
            if `c' <= 7 {
                quietly svyset sdmvpsu [pweight=wt_pfq], strata(sdmvstra) ///
                    vce(linearized) singleunit(centered)
            }
            else {
                quietly svyset sdmvpsu [pweight=wt_dlq], strata(sdmvstra) ///
                    vce(linearized) singleunit(centered)
            }
            quietly svy: poisson `outcome' i.age_grp ///
                if disability==1 & cycle_ord==`c'
            quietly margins age_grp, post
            quietly nlcom (pct: 100*(0.2305*_b[1.age_grp] + 0.2651*_b[2.age_grp] ///
                + 0.3856*_b[3.age_grp] + 0.1187*_b[4.age_grp])), post
            matrix `results'[`row',3] = _b[pct]
            matrix `results'[`row',4] = _se[pct]
            di "`outcome' cycle `c': " %5.2f _b[pct] "% (SE " %5.3f _se[pct] ")"
        }
    }
}

*=============================================================================
* Step 2: inverse-variance-weighted regression, APC per 2-year cycle
*=============================================================================
clear
svmat `results', names(col)
rename c1 outcome_n
rename c2 cyc
rename c3 pct
rename c4 se
drop if pct==.

local oi = 0
foreach outcome of local outcomes {
    local ++oi
    di _n "===== `outcome': overall (cycles 1-10) ====="
    vwls pct cyc if outcome_n==`oi', sd(se)

    di _n "===== `outcome': PFQ era (cycles 1-7) ====="
    vwls pct cyc if outcome_n==`oi' & cyc<=7, sd(se)

    di _n "===== `outcome': DLQ era (cycles 8-10) ====="
    vwls pct cyc if outcome_n==`oi' & cyc>=8, sd(se)
}

log close
