/*
	Education and Political Alignment
	Eduardo Fe
	University of Manchester.
	
	This version: Version 2.0. 
	Date: 15 December 2025
	
	Part 1: Functions
	Part 2: Data Preparation
	Part 3: Analysis. 
	
		
*/
qui {
	
capture program drop rdrandinf_bandgraph
program define rdrandinf_bandgraph
    version 15.1
    syntax varlist(min=2 max=2) ,  C(real) ///
        MINBW(real) MAXBW(real) STEP(real) CI(real) [FUZZY(varname)]

    // Parse variables
    gettoken yvar rest : varlist
    gettoken xvar : rest
    local fuzzyvar `fuzzy'

    // Create temporary dataset to store results
    matrix aux =J(1,5,.)
    matrix list aux

    // Loop over bandwidths
    forvalues bw = `minbw'(`step')`maxbw' {
		
		local lowerbw = `c'-`bw'
		local upperbw = `c'+`bw'
        //rdrandinf `yvar' `xvar', fuzzy(`fuzzyvar' tsls)  c(`c') wl(`lowerbw') wr(`upperbw') ci(`ci')

		if "`fuzzyvar'" != "" {
            rdrandinf `yvar' `xvar', fuzzy(`fuzzyvar' tsls) ///
                c(`c') wl(`lowerbw') wr(`upperbw') ci(`ci')
        }
        else {
            rdrandinf `yvar' `xvar', ///
                c(`c') wl(`lowerbw') wr(`upperbw') ci(`ci')
        }

		
		
        local est   = r(obs_stat)
        local pval  = r(asy_pval)
        matrix CI = r(CI)
        local cilow = CI[1,1]
        local cihigh = CI[1,2]
		matrix newrow = (`est', `pval', `cilow',`cihigh', `bw')
		matrix list newrow
        // Add row
        matrix aux = aux \ newrow
    }

    // Save for later if needed
   svmat aux

    // Graph: estimate with confidence intervals
   
		
	gen est_label = string(aux1, "%9.3f")  // 3 decimals

* Add stars based on p-value thresholds (say aux2 = pval)
replace est_label = est_label + "***" if aux2 < 0.01
replace est_label = est_label + "**"  if aux2 >= 0.01 & aux2 < 0.05
replace est_label = est_label + "*"   if aux2 >= 0.05 & aux2 < 0.1

* Now plot with labels
twoway ///
    (rarea aux3 aux4 aux5, color(gs14%50)) ///
    (line aux1 aux5, lcolor(blue) lwidth(medthick)) ///
    (scatter aux1 aux5, mlabel(est_label) msymbol(i) mcolor(blue) mlabcolor(black) mlabsize(small) mlabangle(90)), ///
    ytitle("Estimate") xtitle("Bandwidth") ///
    legend(off) ///
    name(mainplot, replace)
	
	drop aux* est_label
		
end
	
	
capture program drop my_rdrandinf_function
capture program drop standardize_var
capture program drop trimMean

program define my_rdrandinf_function, rclass
    args wl wr c fuzzy_var running outcomes
    
	local num_outcomes : word count `outcomes'
    
    // Create the matrix with the correct number of rows
    matrix mRes = J(`num_outcomes', 12, .)
	
	
    
    local i=1
    foreach j of local outcomes {
        rdrandinf  `j' `running' if `j'!=., fuzzy(`fuzzy_var') ///
			c(`c') wl(`wl') wr(`wr') ci(0.05) interfci(0.05)

        local itt =   r(obs_stat) 
        local pval =  r(randpval) 
        matrix mRes[`i',1] = `itt'
        matrix mRes[`i',2] = `pval'
		
		matrix aux = r(CI)
		matrix mRes[`i',3] = aux[1,1]
		matrix mRes[`i',4] = aux[1,2]
		
		local interference_lb = r(int_lb)
		local interference_ub = r(int_ub)
		matrix mRes[`i',5] =  `interference_lb'
		matrix mRes[`i',6] =  `interference_ub'
		
		
		
		
		rdrandinf  `j' `running' if `j'!=., fuzzy(`fuzzy_var' tsls) ///
			c(`c') wl(`wl') wr(`wr') ci(0.05)
		 
		 local tsls = r(obs_stat)
		 local pval = r(asy_pval)
		 matrix aux = r(CI)
		 matrix mRes[`i',9] =`tsls'
		 matrix mRes[`i',10] =`pval'
		 matrix mRes[`i',11] =aux[1,1]
		 matrix mRes[`i',12] =aux[1,2]
		 
		 
		rdrandinf  `fuzzy_var'  `running' if `j'!=.,  ///
			c(`c') wl(`wl') wr(`wr') 
		  
        local itt =   r(obs_stat) 
        local pval =  r(randpval) 
		matrix mRes[`i',7 ]=`itt'
		matrix mRes[`i',8 ]=`pval'
		
		
        local i = `i' + 1
    }
    
	 matrix list mRes
	 matrix colnames mRes = itt pval ci_l ci_r ci_lInt ci_rInt first_stage pval tsls pval ci_l ci_r

	 return matrix mRes = mRes
end




program define standardize_var
    // Declare the syntax for the program
    syntax varlist(min=1 max=1) [if] [in], BYVAR(varname) NEWVAR(string)

    // Parse the input
    
	sum `varlist'
	sum `byvar'
	tempvar y
	
	gen `y'= .
	sum `y'
	
	levelsof `byvar', local(levels)
	foreach level of local levels {
		summarize `varlist' if `byvar' == `level'
		local mu = r(mean)
        local sd = r(sd)

        quietly replace `y' = (`varlist' - `mu') / `sd' if `byvar' == `level'
        
	} 
	
	gen `newvar' = `y'
	
	
end




program define trimMean, rclass
	
	* Input is "prop" in proportion, but centile uses percentage.
    syntax varlist(min=1 max=1)  [if] [in], prop(real) [below(int 0)]
    
	tempvar touse
    marksample touse
	
 	
	local cent = `prop'* 100
	
	
    qui centile `varlist' if `touse', centile(`cent') 
	
	local cutoff = r(c_1)
	 
	qui{
	if(`below' == 0){
		sum `varlist' if `varlist'  >= `cutoff' & `touse'
	}
	else{
		sum `varlist' if `varlist' <= `cutoff' & `touse'
	}
	}
	local sol = r(mean)
	return scalar trimmedMean = `sol'
	
end



capture program drop llrplot
program define llrplot, rclass
    version 15.1
    syntax varlist(min=2 max=2) , ///
        Cutoff(real) ///
        BW(real) ///
        Xrange(real) ///
        [ Yrange(numlist min=2 max=2) ]

    // Parse variables
    gettoken yvar xvar : varlist

    // Extract y-range if specified
    if ("`yrange'" != "") {
        local ymin = word("`yrange'", 1)
        local ymax = word("`yrange'", 2)
        local step = (`ymax' - `ymin')/10
        local yopts yscale(range(`ymin' `ymax')) ylabel(`= `ymin''(`step')`= `ymax'' , format(%4.2f))
    }
    else {
        local yopts
    }

    // Clean up previous gen variables if they exist
    cap drop x_left y_left se_left ciL_left ciU_left
    cap drop x_right y_right se_right ciL_right ciU_right

    // Local variables
    local xlow = `cutoff' - `xrange'
    local xhigh = `cutoff' + `xrange'

    // Left of cutoff
    lpoly `yvar' `xvar' if `xvar' < `cutoff', ///
        gen(x_left y_left) at(`xvar') bw(`bw') epanechnikov se(se_left)
    replace x_left = . if `xvar' >= `cutoff'
    replace y_left = . if `xvar' >= `cutoff'
    gen ciL_left = y_left - 1.96*se_left
    gen ciU_left = y_left + 1.96*se_left

    // Right of cutoff
    lpoly `yvar' `xvar' if `xvar' >= `cutoff', ///
        gen(x_right y_right) at(`xvar') bw(`bw') epanechnikov se(se_right)
    replace x_right = . if `xvar' < `cutoff'
    replace y_right = . if `xvar' < `cutoff'
    gen ciL_right = y_right - 1.96*se_right
    gen ciU_right = y_right + 1.96*se_right

    // Compute means of two observations just below and above cutoff
    preserve
        // Below cutoff
        keep if `xvar' < `cutoff'
        gsort -`xvar'
        summarize `yvar' in 1/3, meanonly
        local mean_below = r(mean)
    restore

    preserve
        // Above cutoff
        keep if `xvar' >= `cutoff'
        gsort `xvar'
        summarize `yvar' in 1/3, meanonly
        local mean_above = r(mean)
    restore
	
	 
	local cutM2 = `cutoff'-2
	local cutP2 = `cutoff'+2

    // Plot
    twoway ///
      (rarea ciU_left ciL_left x_left if x_left >= `xlow' & x_left <= `xhigh', color(gs12%40)) ///
      (rarea ciU_right ciL_right x_right if x_right >= `xlow' & x_right <= `xhigh', color(gs12%40)) ///
      (line y_left x_left if x_left >= `xlow' & x_left <= `xhigh', lcolor(black) lwidth(0.5) ) ///
      (line y_right x_right if x_right >= `xlow' & x_right <= `xhigh', lcolor(black) lwidth(0.5)) ///
      (scatter `yvar' `xvar' if `xvar' >= `xlow' & `xvar' <= `xhigh', mcolor(gs10%30) msize(2)) ///
      (function y = `mean_below', range(`xlow' `cutoff')  lpattern("-")  lcolor(black) lwidth(0.2)) ///
      (function y = `mean_above', range(`cutoff' `xhigh')  lpattern("-")  lcolor(black)  lwidth(0.2)) ///
      , xlabel(`xlow'(5)`xhigh') ///
        xline(`cutoff' `cutM2' `cutP2', lpattern(dash)  ) ///
        legend(off) ///
        `yopts'
end





 

capture program drop make_modal
program define make_modal
    syntax varlist(min=1 max=1) , newname(string) codes(numlist)

    * Turn numlist into a comma-separated string
    local clist : subinstr local codes " " ",", all

    * Indicator = 1 if value is in list, 0 otherwise
    gen byte __tmp = inlist(`varlist', `clist')
    replace __tmp = . if `varlist' < 0 | missing(`varlist')

    * Person-level modal behaviour
    by pidp: egen `newname' = mode(__tmp)
    replace `newname' = . if missing(`newname')

    drop __tmp
end



capture program drop make_mean
program define make_mean
    syntax varlist(min=1 max=1) , newname(string) codes(numlist)

    * Turn numlist into a comma-separated string
    local clist : subinstr local codes " " ",", all

    * Indicator = 1 if value is in list, 0 otherwise
    gen byte __tmp = inlist(`varlist', `clist')
    replace __tmp = . if `varlist' < 0 | missing(`varlist')

    * Person-level mean behaviour
    by pidp: egen `newname' = mean(__tmp)
    replace `newname' = . if missing(`newname')

    drop __tmp
end



capture program drop plot_target_vs_others
program define plot_target_vs_others
    
  syntax varlist(min=2 max=2), years(numlist) target(real) [groupvar(varname) bw(real 10)]

    tokenize `varlist'
    local yvar `1'
    local xvar `2'

    // default group variable
    if "`groupvar'" == "" local gvar "doby_fix" 
    else local gvar `groupvar'

    // default bandwidth
    if "`bw'" == "" local bw 10

    // Ensure target is among the supplied years
    local yearlist `years'
    local found = 0
    foreach yy of numlist `yearlist' {
        if "`yy'" == "`target'" local found = 1
    }
    if `found' == 0 {
        di as err "Error: target `target' is not in years(`years')."
        exit 198
    }

    // Build list of others by looping (do NOT try `: list years - target')
    local others ""
    foreach yy of numlist `yearlist' {
        if "`yy'" != "`target'" local others `others' `yy'
    }

    // If no others remain, warn and stop
    if "`others'" == "" {
        di as err "Error: no non-target years (only target provided)."
        exit 198
    }

    // Build faint grey background lines for each non-target year
    local plots ""
    foreach yr of numlist `others' {
        local plots `plots' (lpoly `yvar' `xvar' if `gvar' == `yr', ///
            lcolor(gs12) lwidth(thin) bw(`bw'))
    }

    // Add target line (solid black)
    local plots `plots' (lpoly `yvar' `xvar' if `gvar' == `target', ///
        lcolor(black) lwidth(medthick) bw(`bw'))

		* Compute mean for non-target years by wave_t
	
	 * Build the "if" condition for all non-target years
	local cond ""
	foreach yr of local others {
		if "`cond'" == "" {
			local cond "doby_fix == `yr'"
		}
		else {
			local cond "`cond' | doby_fix == `yr'"
		}
	}

	* Now you can call lpoly with that condition
	local plots `plots'  (lpoly `yvar' `xvar' if `cond', lcolor(black) lpattern(dash) lwidth(medthick) bw(`bw'))


	
    // Draw the graph
    twoway `plots', ///
        legend(off) ///
        ytitle("Proportion") xtitle("`xvar'") ///
        graphregion(color(white)) plotregion(color(white)) bgcolor(white) ///
        scheme(s1color)
		
	 
end

}



capture program drop do_rdd
program define do_rdd, rclass
    // No varlist. Outcomes supplied explicitly.
    args  running fuzzy_var outcomes
    
	local num_outcomes : word count `outcomes'
    
    // Create the matrix with the correct number of rows
    matrix mRES = J(`num_outcomes', 5, .)

    
	// Loop over outcomes
    local i = 0
    foreach y of local outcomes {
        local ++i

        rdrobust `y' `running', c(0) fuzzy(`fuzzy_var') ///
			bwselect(msetwo)

        matrix mRES[`i',1] = e(tau_bc)
        matrix mRES[`i',2] = e(se_tau_rb)
        matrix mRES[`i',3] = e(pv_rb)
        matrix mRES[`i',4] = e(ci_l_rb)
        matrix mRES[`i',5] = e(ci_r_rb)
    }

    // Label columns
    matrix colnames mRES = tau_bc se_rb p_rb ci_l ci_r

    // Label rows with outcomes
    matrix rownames mRES = `outcomes'

    // Return matrix
    return matrix results = mRES
end



**# Compilation of the Panels.

// Define the directory where BHPS data files are stored
local data_dir "/Users/user/Documents/datasets/UKDA-6931-stata_SL/stata/stata13_se/bhps"

// Define a list of waves and the variables you want to extract from each wave
local waves "a b c d e f g h i j k l m n o p q r "
local vars "age* pidp hidp sex doby vote* qfedhi hiqualb_dv istrtdaty scend istrtdaty gor_dv age_dv fimnb jbstat f116 f139 maju dobm*   paju plbornc sctype  "

// Loop over each wave

local i =1
foreach wave of local waves {
	
	if(`i' >= 8){
		local vars "fimnlabgrs_dv fimngrs_dv  istrtdathh istrtdatd pidp hidp sex doby vote* qfedhi  istrtdaty scend istrtdaty gor_dv dobm* age* fimnb jbstat f116 f139 maju   paju plbornc sctype tuin1 "

	}
	else if(`i' >= 11){
	local vars "fimnlabgrs_dv fimngrs_dv j1soc90 istrtdathh istrtdatd pidp hidp sex doby vote* qfedhi  istrtdaty scend istrtdaty gor_dv dobm* age* fimnb jbstat f116 f139 maju   paju plbornc sctype tuin1 masoc00_cc pasoc00_cc"

	}
	else if(`i' >= 12){
	local vars "fimnlabgrs_dv fimngrs_dv j1soc90 istrtdathh istrtdatd pidp hidp sex doby vote* qfedhi  istrtdaty scend istrtdaty gor_dv dobm* age* fimnb jbstat f116 f139 maju   paju plbornc sctype fednt1 tuin1  "

	}
	else{
		local vars "fimnlabgrs_dv fimngrs_dv race istrtdathh istrtdatd pidp hidp sex doby vote* qfedhi  istrtdaty scend istrtdaty gor_dv dobm* age* fimnb jbstat f116 f139 plbornc sctype tuin1 "
	}
	
	if (`i' == 2 | `i' == 4 | `i' == 6 | `i' == 8 | `i' == 11 | `i' == 13 |`i' == 16   ){
		local vars "`vars' oppolc"
	}
	
	if (`i' == 1 | `i' == 3 | `i' == 5 | `i' == 7 | `i' == 10 | `i' == 14 | `i' == 17  ){
		local vars "`vars'   opsoc*  " 
	}
    // Load the dataset for the current wave
	
    use "`data_dir'/b`wave'_indresp_protect.dta", clear

	rename b`wave'_* *
	
	rename doby doby_dv
	
	
    // Keep only the desired variables
    keep `vars'

    // Add a wave identifier variable
    gen wave = "`wave'"
	
	gen wave_t = `i'

    // Save the modified dataset with a new name
    save "`data_dir'/bhps_wave_`wave'_extracted_protect.dta", replace
	
	local i = `i' +1
}

/*
Note on the "vote" variables. 
These range from vote1 to vote9 in the bhps.

 - They first ask if a person supports a party (vote1).
 - If the person says no, then they ask if that person fells closer to one party
than another (vote2)

 - vote4 (which party) is then asked if vote1 =1 or vote1=2 AND vote2 = 1
 - vote5 (how strongly) if vote1 =1 or vote1=2 AND vote2 = 1
*/

rename qfedhi qfhigh
save "`data_dir'/bhps_wave_`wave'_extracted_protect.dta", replace


// to merge all waves into a single dataset
clear
foreach wave of local waves {
    use "`data_dir'/bhps_wave_`wave'_extracted_protect.dta", clear
    if "`wave'" == "a" {
        save "`data_dir'/bhps_all_waves_protect.dta", replace
    }
    else {
        append using "`data_dir'/bhps_all_waves_protect.dta"
        save "`data_dir'/bhps_all_waves_protect.dta", replace
    }
}


* Repeat now for Understanding society
* ------------------------------------------------------------------------------
* ------------------------------------------------------------------------------


// Define the directory where BHPS data files are stored
local data_dir "/Users/user/Documents/datasets/UKDA-6931-stata_SL/stata/stata13_se/ukhls"


// Define a list of waves and the variables you want to extract from each wave
local waves "a b c d e f g h i j k l "

// Loop over each wave

 
foreach wave of local waves {
    // Load the dataset for the current wave
	
    use "`data_dir'/`wave'_indresp_protect.dta", clear

	rename `wave'_* *
	
    // Keep only the desired variables
    if("`wave'" == "c"){
	
	local vars "big5* cg*_dv fimnlabgrs_dv fimngrs_dv j1soc90 age_dv ethn_dv istrtdathh istrtdatd  pidp hidp sex doby_dv vote* jbstat   istrtdaty scend gor_dv age_dv fimnprben_dv fimnsben_dv bendis3 bendis5 maju dobm*  paju plbornc fimngrs_dv qfhigh "
	
	}
	else if ("`wave'" == "a"){

		local vars "fimnlabgrs_dv fimngrs_dv j1soc90  age_dv ethn_dv istrtdathh istrtdatd  pidp hidp sex doby_dv vote* jbstat   istrtdaty scend gor_dv age_dv fimnprben_dv fimnsben_dv bendis3 bendis5 maju   paju plbornc dobm* fimngrs_dv qfhigh "
		
	}
	else {

		local vars "fimnlabgrs_dv fimngrs_dv j1soc90  age_dv ethn_dv istrtdathh istrtdatd  pidp hidp sex doby_dv vote* jbstat dobm*  istrtdaty scend gor_dv age_dv fimnprben_dv fimnsben_dv bendis3 bendis5 maju   paju plbornc fimngrs_dv ftqual qfhigh masoc00_cc pasoc00_cc "
		
	}
	
	if("`wave'" == "b" | "`wave'" == "d"| "`wave'" == "f"| "`wave'" == "h" | "`wave'" == "j"){
		local vars "`vars' chargv volun tuin1"
	}
	if("`wave'" == "d" | "`wave'" == "h"| "`wave'" == "l"){
		local vars "`vars'  oprlg2"
	}
	
	if("`wave'" == "d" | "`wave'" == "j"){
		local vars "`vars'  scenv*"
	}
	
	if("`wave'" == "h" | "`wave'" == "j"| "`wave'" == "k" | "`wave'" == "l" ){
		local vars "`vars'  eumem"
	}
	
	
	if( "`wave'" == "l" ){
		local vars "`vars'  opsoc* tuin1 " 
	}
   	
	
	keep `vars'

    // Add a wave identifier variable
    gen wave = "us_`wave'"
	
	
	gen wave_t = `i'
    // Save the modified dataset with a new name
    save "`data_dir'/us_wave_`wave'_extracted_protect.dta", replace
	
	local i = `i' +1
}

// to merge all waves into a single dataset
clear
foreach wave of local waves {
    use "`data_dir'/us_wave_`wave'_extracted_protect.dta", clear
    if "`wave'" == "a" {
        save "`data_dir'/us_all_waves_protect.dta", replace
    }
    else {
        append using "`data_dir'/us_all_waves_protect.dta"
        save "`data_dir'/us_all_waves_protect.dta", replace
    }
}
/*

	Note on the "vote" variables in US. 
	These range from vote1 to vote10 in US.

	 - They first ask if a person supports a party (vote1).
	 - If the person says no, then they ask if that person fells closer to one party
	than another (vote2)

	 - vote4 (which party) is then asked if vote1 =1 or vote1=2 AND vote2 = 1
	 - vote5 (how strongly) if vote1 =1 or vote1=2 AND vote2 = 1
	 - Then vote6 is asked to everyone "interestest in politics".
 
*/

local data_dir "/Users/user/Documents/datasets/UKDA-6931-stata_SL/stata/stata13_se/bhps"

append using "`data_dir'/bhps_all_waves_protect.dta"

save "/Users/user/Dropbox/Econometrics/educationLiberalism/code/all_waves_protect.dta", replace




**# Prepare the outcome variables.
*=============================================================================== 

use "/Users/user/Dropbox/Econometrics/educationLiberalism/code/all_waves_protect.dta", clear
   
xtset pidp wave_t

* Modal behaviour

make_modal vote6, newname(noInterestedPolitics) codes(3 4)
make_modal vote4, newname(mostlyLabour) codes(2)
make_modal vote4, newname(mostlyConservative) codes(1)
make_modal vote4, newname(mostlyLeft) codes(2 4 5 6 8 11)
make_modal vote4, newname(mostlyRight) codes(1 7 10 12 13 14)
make_modal vote4, newname(mostlyCentrist) codes(3 9)

make_modal vote8, newname(voted_mostlyLabour) codes(2)
make_modal vote8, newname(voted_mostlyConservative) codes(1)
make_modal vote8, newname(voted_mostlyLeft) codes(2 4 5 6 8 11)
make_modal vote8, newname(voted_mostlyRight) codes(1 7 10 12 13 14)
make_modal vote8, newname(voted_mostlyCentrist) codes(3 9)


make_mean vote4, newname(onAverageLeft) codes(2 4 5 6 8 11)
make_mean vote4, newname(onAverageRight) codes(1 7 10 12 13 14)

 


* Vote dynamics.
gen aparty = vote4 
replace aparty =. if aparty <0 | aparty >90

label var aparty "Which party supports "

by pidp: gen aparty_lag = aparty[_n-1]
gen change_affiliation = 0
replace change_affiliation = 1 if aparty != aparty_lag & !missing(aparty) & !missing(aparty_lag)
by pidp: egen ever_chaffiliation = max(change_affiliation)

* Voty dynamics descriptives.
tab  aparty_lag aparty if  (aparty <=3 | aparty ==12 | aparty ==6) ///
	& (aparty_lag <=3 | aparty_lag ==12 | aparty_lag ==6) & gor <10 & gor >0 , row matcell(count) 
	
tab  aparty_lag aparty if change_aff==1& (aparty <=3 | aparty ==12 | aparty ==6) ///
	& (aparty_lag <=3 | aparty_lag ==12 | aparty_lag ==6) & gor <10 & gor >0 , row matcell(count) 
	
tab  aparty_lag aparty if ever_chaff==1 & (aparty <=3 | aparty ==12 | aparty ==5| aparty ==6) & (aparty_lag <=3 | aparty_lag ==12 | aparty_lag ==5| aparty_lag ==6) & gor ==11 & gor >0, row

* Fix doby (sometimes missing)
by pidp: egen doby_fix = max(doby_)

* Fix sex
by pidp: egen sex_fix = max(sex)

gen female = sex_fix == 2
replace female = . if sex_fix <0

* Leaving age (no instrument)

by pidp: egen schleavingage = max(scend)
replace schleavingage = . if schleavingage <0

gen dropAt14 = sch == 14
replace dropAt14 =. if sch ==.

gen dropAt15 = sch == 15
replace dropAt15 =. if sch ==.

gen dropAt16 = sch == 16
replace dropAt16 =. if sch ==.

gen dropAfter16 = sch >16
replace dropAfter16 =. if sch ==.
* Instruments

gen left15=sch >=15 
replace left15 =. if sch ==.

gen left16=sch >=16
replace left16 =. if sch ==.


* Unemployed

gen unemployed = (jbstat == 3) 
replace unemployed = . if jbstat <0 | jbstat == 97
by pidp: egen everunemployed = max(unemployed)

* Ethnicity
 
by pidp: egen ethg = max(race)
by pidp: egen ethg2 = max(ethn_dv)
gen ethg3 = max(ethg2, ethg)
by pidp: egen aux = max(ethg3)
gen white = aux == 1
drop aux

* Place born
gen aux = plbornc>0
by pidp: egen nouk = max(aux)
drop aux

* Parental background 
gen aux = maju == 1
replace aux =. if maju <0 | maju ==.
by pidp: egen mumworked14 = max(aux)
drop aux

gen aux = paju == 1
replace aux =. if paju <0 | paju ==.
by pidp: egen dadworked14 = max(aux)
drop aux

gen aux = maju == 3
replace aux =. if maju <0 | maju ==.
by pidp: egen mumdead14 = max(aux)
drop aux

gen aux = paju == 3
replace aux =. if paju <0 | paju ==.
by pidp: egen daddead14 = max(aux)
drop aux

gen aux = maju == 4
replace aux =. if maju <0 | maju ==.
by pidp: egen mumaway14 = max(aux)
drop aux

gen aux = paju == 4
replace aux =. if paju <0 | paju ==.
by pidp: egen dadaway14 = max(aux)
drop aux

 
* Cognition and personality
egen agecat =cut(age_dv), group(10)

gen cgnumability = cgna_dv 
replace cgnumability = . if cgnumability<0

gen cgverbfluency = cgvfc_dv
replace cgverbfluency =. if cgverbf<0

egen cgmemory = rowmean(cgwri_dv cgwrd_dv)
replace cgmemory = . if cgmemory <0

standardize_var cgnumability, byvar(agecat) newvar(numability_sd)
standardize_var cgmemory, byvar(agecat) newvar(memory_sd)
standardize_var cgverbfluency, byvar(agecat) newvar(verbfluency_sd)

by pid: egen cg_num = mean(numability_sd)
by pid: egen cg_vrb = mean(memory_sd)
by pid: egen cg_mem = mean(verbfluency_sd)

gen cg_num_dis = cg_num >0 & cg_num !=.
replace cg_num_dis=. if cg_num ==.

gen cg_vrb_dis = cg_vrb >0 & cg_vrb !=.
replace cg_vrb_dis = . if cg_vrb==.

gen cg_mem_dis = cg_mem >0 & cg_mem !=.
replace cg_mem_dis =. if cg_mem ==.

* Compute historical income, conditional on age, sex

gen year = istrtdaty
replace year =. if istrtdaty<=0
merge m:1 year using "/Users/user/Documents/datasets/UKDA-6614-stata/stata/stata13_se/deflator2022prices.dta"
gen income2024 = fimngrs_dv/cpi_index
replace income2024=. if income2024 <0

gen loginc = log(income2024)

drop _merge

gen age2 =age_dv^2
gen age3 = age_dv^3
gen ageSex = age_dv*sex
gen age4= age_dv^4
gen age5= age_dv^5
rlasso loginc age_dv age2 age3 age4 age5 ageSex sex if fimngrs_dv !=. 
predict xb
gen resY = fimngrs_dv- xb
xtset pidp wave_t

by pidp: egen aveinc = mean(resY)


/*
	* Key UK elections (including Bridgwater by-election)
	* Format: name date year

	* General Election        	18 Jun 1970  => 1970.47
	* General Election        	28 Feb 1974  => 1974.16
	* General Election         	10 Oct 1974  => 1974.78
	* General Election           3 May 1979  => 1979.34
	* General Election           9 Jun 1983  => 1983.44
	* General Election          11 Jun 1987  => 1987.45
	* General Election           9 Apr 1992  => 1992.27
	* General Election           1 May 1997  => 1997.33


	gen first_election = ""

	* Approximate age 18 as birthyear + 18.5 (assume July 1 birthday)
	gen voting_year = doby_fix + 18.5

	* Assume you already have year of birth in variable `yob`
	* Step 1: Calculate the year the person turns 18 (assuming July 1 birthday)

	* Step 2: Create a numeric variable for first election (coded by year)
	gen first_election_num = .

	* Election codes (you can change to use codes like 1970.13 for Bridgwater if desired)
	* We'll use the **year of the election** as the code

	replace first_election_num = 1970.13 if voting_year <= 1970.20   // Bridgwater by-election
	replace first_election_num = 1970    if voting_year > 1970.20  & voting_year <= 1974.16
	replace first_election_num = 1974.1  if voting_year > 1974.16  & voting_year <= 1974.78   // Feb 1974
	replace first_election_num = 1974.2  if voting_year > 1974.78  & voting_year <= 1979.34   // Oct 1974
	replace first_election_num = 1979    if voting_year > 1979.34  & voting_year <= 1983.44
	replace first_election_num = 1983    if voting_year > 1983.44  & voting_year <= 1987.45
	replace first_election_num = 1987    if voting_year > 1987.45  & voting_year <= 1992.27
	replace first_election_num = 1992    if voting_year > 1992.27  & voting_year <= 1997.33
	replace first_election_num = 1997    if voting_year > 1997.33
*/


* Step 1: Estimate year of birth
* Step 1: Create variable to store first eligible election
gen first_election = .

* Step 2: Assign first general election based on year of birth and voting age rules
replace first_election = 1900 if missing(first_election) & doby_fix <= 1879  // age 21 by 1900
replace first_election = 1906 if missing(first_election) & doby_fix <= 1885
replace first_election = 1910 if missing(first_election) & doby_fix <= 1889
replace first_election = 1918 if missing(first_election) & doby_fix <= 1897  // men 21, women 30+
replace first_election = 1922 if missing(first_election) & doby_fix <= 1901
replace first_election = 1923 if missing(first_election) & doby_fix <= 1902
replace first_election = 1924 if missing(first_election) & doby_fix <= 1903
replace first_election = 1929 if missing(first_election) & doby_fix <= 1908  // men & women 21+
replace first_election = 1931 if missing(first_election) & doby_fix <= 1910
replace first_election = 1935 if missing(first_election) & doby_fix <= 1914
replace first_election = 1945 if missing(first_election) & doby_fix <= 1924
replace first_election = 1950 if missing(first_election) & doby_fix <= 1929
replace first_election = 1951 if missing(first_election) & doby_fix <= 1930
replace first_election = 1955 if missing(first_election) & doby_fix <= 1934
replace first_election = 1959 if missing(first_election) & doby_fix <= 1938
replace first_election = 1964 if missing(first_election) & doby_fix <= 1943
replace first_election = 1966 if missing(first_election) & doby_fix <= 1945
replace first_election = 1970 if missing(first_election) & doby_fix <= 1952  // voting age lowered to 18
replace first_election = 1974 if missing(first_election) & doby_fix <= 1956
replace first_election = 1979 if missing(first_election) & doby_fix <= 1961
replace first_election = 1983 if missing(first_election) & doby_fix <= 1965
replace first_election = 1987 if missing(first_election) & doby_fix <= 1969
replace first_election = 1992 if missing(first_election) & doby_fix <= 1974
replace first_election = 1997 if missing(first_election) & doby_fix <= 1979

save "/Users/user/Dropbox/Econometrics/educationLiberalism/code/all_waves_protect.dta", replace




**# Start Analysis Here.

use "/Users/user/Dropbox/Econometrics/educationLiberalism/code/all_waves_protect.dta", clear

* ========================= Part 3: Analysis =========================

encode wave, gen(wave_num)
xtset pidp wave_num

replace dobm_dv = dobm if dobm_dv ==.

by pidp: gen doby_fwd = doby_fix[1] if !missing(doby_fix[1])
by pidp: replace doby_fwd = doby_fix if missing(doby_fwd) & !missing(doby_fix)

by pidp: gen dobm_fwd = dobm_dv[1] if !missing(dobm_dv[1])
by pidp: replace dobm_fwd = dobm_dv if missing(dobm_fwd) & !missing(dobm_dv)

by pidp (wave_num): gen doby_bwd = doby_fix[_N] if !missing(doby_fix[_N])
by pidp (wave_num): replace doby_bwd = doby_fix if missing(doby_bwd) & !missing(doby_fix)

by pidp (wave_num): gen dobm_bwd = dobm_dv[_N] if !missing(dobm_dv[_N])
by pidp (wave_num): replace dobm_bwd = dobm_dv if missing(dobm_bwd) & !missing(dobm_dv)

gen doby_clean = doby_fix
replace doby_clean = doby_fwd if missing(doby_clean) & !missing(doby_fwd)
replace doby_clean = doby_bwd if missing(doby_clean) & !missing(doby_bwd)
gen dobm_clean = dobm_dv
replace dobm_clean = dobm_fwd if missing(dobm_clean) & !missing(dobm_fwd)
replace dobm_clean = dobm_bwd if missing(dobm_clean) & !missing(dobm_bwd)
drop doby_fwd doby_bwd dobm_fwd dobm_bwd



 
gen vulnerable = dadworked14 ==0 | mumdead14==1 | daddead14==1| mumaway14 ==1 | dadaway14==1  | mumworked14 ==0


* Ordinary people share nations wealth
gen aux = opsoca == 1 | opsoca == 2
replace aux = . if opsoca < 0 | opsoca == .
bysort pidp: egen shareWealth = max(aux)
drop aux

* One law for rich and one for poor
gen aux = opsocb == 1 | opsocb == 2
replace aux = . if opsocb < 0 | opsocb == .
bysort pidp: egen oneLaw = max(aux)
drop aux

* Private enterprise solves economic probs
gen aux = opsocc == 1 | opsocc == 2
replace aux = . if opsocc < 0 | opsocc == .
bysort pidp: egen privateEnterprise = max(aux)
drop aux

* Public services ought to be state owned
gen aux = opsocd == 1 | opsocd == 2
replace aux = . if opsocd < 0 | opsocd == .
bysort pidp: egen stateOwned = max(aux)
drop aux

* Govt. has obligation to provide jobs
gen aux = opsoce == 1 | opsoce == 2
replace aux = . if opsoce < 0 | opsoce == .
bysort pidp: egen provideJobs = max(aux)
drop aux

* Strong trade unions protect employees
gen aux = opsocf == 1 | opsocf == 2
replace aux = . if opsocf < 0 | opsocf == .
bysort pidp: egen strongUnions = max(aux)
drop aux

* Young people don't respect British values
gen aux = opsock == 1 | opsock == 2
replace aux = . if opsock < 0 | opsock == .
bysort pidp: egen britishValues = max(aux)
drop aux

* Censorship necessary to uphold moral standards
gen aux = opsocl == 1 | opsocl == 2
replace aux = . if opsocl < 0 | opsocl == .
bysort pidp: egen censorship = max(aux)
drop aux

* Public meetings protest against government allowed
gen aux = opsocm == 1 | opsocm == 2
replace aux = . if opsocm < 0 | opsocm == .
bysort pidp: egen protestAllowed = max(aux)
drop aux

* More tolerant of unconventional lives
gen aux = opsocn == 1 | opsocn == 2
replace aux = . if opsocn < 0 | opsocn == .
bysort pidp: egen tolerantLives = max(aux)
drop aux

* Death penalty appropriate some crimes
gen aux = opsoco == 1 | opsoco == 2
replace aux = . if opsoco < 0 | opsoco == .
bysort pidp: egen deathPenalty = max(aux)
drop aux

* Stiffer sentences
gen aux = opsocp == 1 | opsocp == 2
replace aux = . if opsocp < 0 | opsocp == .
bysort pidp: egen stiffSentences = max(aux)
drop aux


* Qualifications in the US panel

gen aux =.
replace aux = qfhigh if qfhigh>0 &qfhigh!=.
bysort pidp: egen qualx = max(aux)
drop aux

* qualifications in the BHPS
bysort pidp:  egen qualReceived = max(qfedhi)
replace qualReceived=. if qualReceived<0


gen noqualifications = qualReceived==12 | qualx == 96
gen gce25 = qualReceived ==9 | qualx == 13
gen gceO = qualReceived ==7 | qualx == 12
gen degree= qualReceived== 1 | qualReceived==2 | qualx == 1 |qualx == 2


bysort pidp: egen noQual = max(noqualifications)
bysort pidp: egen gceGrade25 = max(gce25)
bysort pidp: egen gceOlevel = max(gceO)
bysort pidp: egen uniDegree = max(degree)



collapse cg_*_dis  aveinc female dropAt14 dropAt16 dropAt15 dropAfter noInterestedPolitics mostly* left15 left16 doby_fix nouk white mumworked14 dadworked14 mumdead14 daddead14 mumaway14 dadaway14 everunemployed loginc first_election   dob*_clean voted_* vulnerable shareWealth oneLaw privateEnterprise stateOwned provideJobs strongUnions britishValues censorship protestAllowed tolerantLives deathPenalty stiffSentences gce* noQual uniDegree onAverageLeft onAverageRight , by(pidp)

replace dobm_clean = floor(dobm_clean)

drop if dobm_clean ==. | doby_clean ==.  | dobm_clean <0 | doby_clean <0
 

* Alternative quarters 

* Step 1: Define custom quarter from month
gen custom_qtr = .
replace custom_qtr = 0 if inlist(dobm_clean,12,1,2)   // Q1 = Dec, Jan, Feb
replace custom_qtr = 1 if inlist(dobm_clean,3,4,5)    // Q2 = Mar, Apr, May
replace custom_qtr = 2 if inlist(dobm_clean,6,7,8)    // Q3 = Jun, Jul, Aug
replace custom_qtr = 3 if inlist(dobm_clean,9,10,11)  // Q4 = Sep, Oct, Nov

local base_year = 1957
local base_qtr = 3   // your custom Q4 = 3
gen aux = dobm_clean == 12

* Running variable (manual calculation)
gen q_birth_run = (doby_clean - `base_year')*4 + (custom_qtr - `base_qtr') + 4*aux



//gen q_birth_run = (doby_clean - `base_year')*4 + (custom_qtr - `base_qtr') + ///
                //  (custom_qtr < `base_qtr')
				  
* Step 1: Define bimesters starting in January
gen custom_bim = .
replace custom_bim = 0 if inlist(dobm_clean,1,2)      // Bim1 = Jan, Feb
replace custom_bim = 1 if inlist(dobm_clean,3,4)      // Bim2 = Mar, Apr
replace custom_bim = 2 if inlist(dobm_clean,5,6)      // Bim3 = May, Jun
replace custom_bim = 3 if inlist(dobm_clean,7,8)      // Bim4 = Jul, Aug
replace custom_bim = 4 if inlist(dobm_clean,9,10)     // Bim5 = Sep, Oct
replace custom_bim = 5 if inlist(dobm_clean,11,12)    // Bim6 = Nov, Dec
local base_year = 1957
local base_bim = 4   // Sep-Oct 1957 = 4

* Running variable for bimesters
gen bim_birth_run = (doby_clean - `base_year')*6 + (custom_bim - `base_bim')



* Step 1: Define custom 12-month periods starting in September
* Months: Sep=1, Oct=2, ..., Aug=12

gen custom_12m = .
replace custom_12m = 1  if inlist(dobm_clean, 1)
replace custom_12m = 2  if inlist(dobm_clean,2)
replace custom_12m = 3  if inlist(dobm_clean,3)
replace custom_12m = 4  if inlist(dobm_clean,4)
replace custom_12m = 5  if inlist(dobm_clean, 5)
replace custom_12m = 6  if inlist(dobm_clean, 6)
replace custom_12m = 7  if inlist(dobm_clean, 7)
replace custom_12m = 8  if inlist(dobm_clean, 8)
replace custom_12m = 9  if inlist(dobm_clean, 9)
replace custom_12m = 10 if inlist(dobm_clean, 10)
replace custom_12m = 11 if inlist(dobm_clean, 11)
replace custom_12m = 12 if inlist(dobm_clean, 12)

* Base period: September 1957 -> custom_12m = 1
local base_year = 1957
local base_period = 9

* Running variable: each period is 1 year long
gen month_birth_run = (doby_clean - `base_year')*12 + (custom_12m - `base_period')



gen year_birth_run = doby_clean - 1957
replace year_birth_run = year_birth_run - 1 if dobm_clean < 9





**# The outcomes.

gen wrong_run = doby_fix - 1958

local outcomes loginc mostlyLabour mostlyConservative ///
    mostlyLeft mostlyRight mostlyCentrist voted_mostlyLabour ///
    voted_mostlyConservative voted_mostlyLeft voted_mostlyRight ///
    voted_mostlyCentrist shareWealth oneLaw privateEnterprise ///
    stateOwned provideJobs strongUnions britishValues ///
    censorship protestAllowed tolerantLives deathPenalty stiffSentences ///
	gce25 gceGrade25 gceOlevel noQual noInterestedPolitics onAverageLeft onAverageRight

* List of running variables and corresponding matrix names
local runs month_birth_run bim_birth_run q_birth_run year_birth_run wrong_run
local num_outcomes : word count `outcomes'

foreach running of local runs {
	
	matrix mRES = J(`num_outcomes', 10, .)
	local i = 0
	foreach y of local outcomes{
		local ++i
		rdrobust `y' `running',  bwselect(msetwo)
		matrix mRES[`i',1] = e(tau_bc)
        matrix mRES[`i',2] = e(se_tau_rb)
        matrix mRES[`i',3] = e(pv_rb)
        //matrix mRES[`i',4] = e(ci_l_rb)
       // matrix mRES[`i',5] = e(ci_r_rb)
		rdrobust `y' `running', fuzzy(left16) bwselect(msetwo)
		matrix mRES[`i',6] = e(tau_bc)
        matrix mRES[`i',7] = e(se_tau_rb)
        matrix mRES[`i',8] = e(pv_rb)
       // matrix mRES[`i',9] = e(ci_l_rb)
       // matrix mRES[`i',10] = e(ci_r_rb)
		
	}
	matrix colnames mRES = tau_bc se_rb p_rb ci_l ci_r tau_bc se_rb p_rb ci_l ci_r
    matrix rownames mRES = `outcomes'
	outtable using "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/`running'", mat(mRES) replace nobox
}


/*
. rdwinselect month_birth_run white nouk mumworked14 dadworked14 mumdead14 daddead14 mumaway14 dadaway14  first , wmass level(0.05)  wasym  nw(20)
Missing values detected in covariates
Consider dropmissing option to exclude missing values
Mass points detected in running variable
You may use wmasspoints option for constructing windows at each mass point


Window selection for RD under local randomization


Cutoff c = 0.00   | Left of c   Right of c        Number of obs  =        108606
------------------+-----------------------        Order of poly  =             0
    Number of obs |     34765        73841        Kernel type    =       uniform
   1st percentile |       364          857        Reps           =          1000
   5th percentile |      1744         3761        Testing method =     rdrandinf
  10th percentile |      3560         7527        Balance test   =     diffmeans
  20th percentile |      6980        14855


                  |   Bal. test         Var. name    Bin. test 
      Window      |    p-value        (min p-value)   p-value     Obs<c   Obs>=c
------------------+-------------------------------------------------------------
  -1.000|   0.000 |      0.252             nouk         0.353        64       76
  -2.000|   1.000 |      0.300        mumaway14         0.390       144      160
  -3.000|   2.000 |      0.162        mumdead14         0.745       228      236
  -4.000|   3.000 |      0.174        mumdead14         1.000       303      303
  -5.000|   4.000 |      0.222        mumdead14         0.642       384      398
  -6.000|   5.000 |      0.214        mumdead14         0.672       463      477
  -7.000|   6.000 |      0.300        mumdead14         0.382       535      565
  -8.000|   7.000 |      0.214        mumdead14         0.325       614      650
  -9.000|   8.000 |      0.000       first_election     0.232       686      732
 -10.000|   9.000 |      0.000       first_election     0.122       746      808
 
Variable used in binomial test (running variable): month_birth_run
Covariates used in balance test: white nouk mumworked14 dadworked14 mumdead14 da
> ddead14 mumaway14 dadaway14 first_election

Recommended window is [-8.000;  7.000] with 1264 observations (614 below, 650 above).


. rdwinselect q_birth_run white nouk mumworked14 dadworked14 mumdead14 daddead14 mumaway14 dadaway14  first , wmass level(0.05)  wasym  nw(20)
Missing values detected in covariates
Consider dropmissing option to exclude missing values
Mass points detected in running variable
You may use wmasspoints option for constructing windows at each mass point


Window selection for RD under local randomization


Cutoff c = 0.00   | Left of c   Right of c        Number of obs  =        108606
------------------+-----------------------        Order of poly  =             0
    Number of obs |     34765        73841        Kernel type    =       uniform
   1st percentile |       364         1118        Reps           =          1000
   5th percentile |      1744         3761        Testing method =     rdrandinf
  10th percentile |      3790         7527        Balance test   =     diffmeans
  20th percentile |      6980        15014


                  |   Bal. test         Var. name    Bin. test 
      Window      |    p-value        (min p-value)   p-value     Obs<c   Obs>=c
------------------+-------------------------------------------------------------
  -1.000|   0.000 |      0.222        mumdead14         0.745       228      236
  -2.000|   1.000 |      0.232        mumdead14         0.672       463      477
  -3.000|   2.000 |      0.000       first_election     0.232       686      732
  -4.000|   3.000 |      0.000       first_election     0.136       889      954
  -5.000|   4.000 |      0.000       first_election     0.083      1105     1189
  -6.000|   5.000 |      0.000       first_election     0.290      1374     1431
  -7.000|   6.000 |      0.000       first_election     0.125      1598     1687
  -8.000|   7.000 |      0.000       first_election     0.012      1797     1951
  -9.000|   8.000 |      0.000       first_election     0.007      1999     2175
 -10.000|   9.000 |      0.000       first_election     0.002      2237     2452


Variable used in binomial test (running variable): q_birth_run
Covariates used in balance test: white nouk mumworked14 dadworked14 mumdead14 da
> ddead14 mumaway14 dadaway14 first_election

Recommended window is [-2.000;  1.000] with 940 observations (463 below, 477 above).




*/
local outcomes loginc mostlyLabour mostlyConservative ///
    mostlyLeft mostlyRight voted_mostlyLabour ///
    voted_mostlyConservative voted_mostlyLeft voted_mostlyRight ///
    voted_mostlyCentrist shareWealth oneLaw privateEnterprise ///
    stateOwned provideJobs strongUnions britishValues ///
    censorship protestAllowed tolerantLives deathPenalty stiffSentences ///
    gceGrade25 gceOlevel noQual noInterestedPolitics

* List of running variables and corresponding bandwidths
local runs month_birth_run q_birth_run year_birth_run wrong_run
local wl_list -8 -2 -1 -1
local wr_list 7 1 0.5 0.5

local num_outcomes : word count `outcomes'

* Loop over running variables
forvalues j = 1/`=wordcount("`runs'")' {
    local running : word `j' of `runs'
    local wl : word `j' of `wl_list'
    local wr : word `j' of `wr_list'
    
    matrix mRES = J(`num_outcomes', 10, .)
    local i = 0
    
    foreach y of local outcomes {
        local ++i
        
        * ITT
        rdrandinf `y' `running', wl(`wl') wr(`wr') //ci(0.05)
        matrix V = .
        matrix CI = r(CI)
        matrix mRES[`i',1] = r(obs_stat)
        matrix mRES[`i',2] = .
        matrix mRES[`i',3] = r(asy_pval)
        //matrix mRES[`i',4] = CI[1,1]
        //matrix mRES[`i',5] = CI[1,2]
        
        * LATE (fuzzy)
        rdrandinf `y' `running', fuzzy(left16 tsls) wl(`wl') wr(`wr')  //ci(0.05)
        matrix V = e(V)
        matrix CI = r(CI)
        matrix mRES[`i',6] = r(obs_stat)
        matrix mRES[`i',7] = V[1,1]
        matrix mRES[`i',8] = r(asy_pval)
        //matrix mRES[`i',9] = CI[1,1]
        //matrix mRES[`i',10] = CI[1,2]
    }
    
    matrix colnames mRES = tau_bc se_rb p_rb ci_l ci_r tau_bc se_rb p_rb ci_l ci_r
    matrix rownames mRES = `outcomes'
    
    outtable using "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/ri_`running'", mat(mRES) replace nobox
}




**# Descriptive Statistics

* descriptive statistics.

local explanatoryVariables female dropAt14 dropAt16 dropAt15 doby_fix nouk white mumworked14 ///
	dadworked14 mumdead14 daddead14 mumaway14 dadaway14 everunemployed ///
	 first_election aveinc loginc  mostlyLabour mostlyConservative mostlyLeft mostlyRight mostlyCentrist voted_mostlyLabour voted_mostlyConservative voted_mostlyLeft voted_mostlyRight voted_mostlyCentrist noInterestedPolitics

foreach j of local explanatoryVariables {

	replace `j' = .  if `j' <0
}
estpost sum `explanatoryVariables' 
est store des1

estpost sum `explanatoryVariables'  if inrange(year_birth_run, -1,0)
est store des2

estpost sum `explanatoryVariables'  if year_birth_run ==-1
est store des3

estpost sum `explanatoryVariables'  if year_birth_run==0
est store des4

esttab des1 des2 des3 des4 using "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/descriptives.tex", ///
	cells(mean(fmt(2)) sd(par fmt(2))  count() ) label booktabs  collabels(none) gaps f noobs ///
	replace se star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) se(%9.3f) 
	
gen below_bw = .
replace  below_bw = 1 if  year_birth_run==-1
replace below_bw = 0 if year_birth_run==0

local explanatoryVariables female dropAt14 dropAt16 dropAt15   white mumworked14 ///
	dadworked14 mumdead14 daddead14 mumaway14 dadaway14 everunemployed  first_election ///
	mostlyLabour mostlyConservative mostlyLeft mostlyRight mostlyCentrist voted_mostlyLabour voted_mostlyConservative voted_mostlyLeft voted_mostlyRight voted_mostlyCentrist noInterestedPolitics
	
matrix pvalues = J(24,1,.)
local i = 1
foreach j of local explanatoryVariables {

	ttest  `j' , by(below)
	matrix pvalues[`i',1] = r(p)
	local i = `i'+1
}
 
matrix rownames pvalues = `explanatoryVariables'
 
**# Plots
 
 
local agevar year_birth_run

graph bar if inrange(`agevar', -24, 24),  ylabel(, nogrid) over(`agevar', label(angle(vertical) labsize(small))  ) bar(1, fcolor(black) lcolor(none)) ytitle(Percentage)  graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white) lcolor(white) ) 
 

 


 preserve
 local agevar wrong_run

collapse left15 left16 dropAt15 dropAt16 dropAfter16, by(`agevar')

if "`agevar'" == "year_birth_run" {
    local xtag "Years from cut-off"
}
else if "`agevar'" == "q_birth_run" {
    local xtag "Quarters from cut-off"
}
else if "`agevar'" == "wrong_run"{
	local xtag "Year of birth"
}
else {
    local xtag "Months from cut-off"
}
local lb = -24
local ub = 24
local cutoff = 0

twoway ///
    (lpoly dropAt15 `agevar' if `agevar' < `cutoff' & `agevar' >= `lb' & `agevar' <= `ub', ///
        lcolor(black) lpattern(.#.) lwidth(0.5)) ///
    (lpoly dropAt15 `agevar' if `agevar' >= `cutoff' & `agevar' >= `lb' & `agevar' <= `ub', ///
        lcolor(black) lpattern(.#.) lwidth(0.5)) ///
    (scatter dropAt15 `agevar' if `agevar' >= `lb' & `agevar' <= `ub', ///
        mcolor(gray) msymbol(Oh) msize(small) ///
        xtitle(`xtag') ytitle(Percentage) ///
        graphregion(fcolor(white) lcolor(white)) ///
        plotregion(fcolor(white) lcolor(white)) ///
        legend(off) xlabel(, angle(45)) ///
        xscale(range(`lb' `ub'))) ///
    (lpoly dropAt16 `agevar' if `agevar' < `cutoff' & `agevar' >= `lb' & `agevar' <= `ub', ///
        lcolor(black)) ///
    (lpoly dropAt16 `agevar' if `agevar' >= `cutoff' & `agevar' >= `lb' & `agevar' <= `ub', ///
        lcolor(black)) ///
    (scatter dropAt16 `agevar' if `agevar' >= `lb' & `agevar' <= `ub', ///
        xline(`cutoff', lcolor(black)) ///
        mcolor(gray) msymbol(O) msize(small)) ///
    (lpoly dropAfter16 `agevar' if `agevar' < `cutoff' & `agevar' >= `lb' & `agevar' <= `ub', ///
        lcolor(black) lpattern(shortdash)) ///
    (lpoly dropAfter16 `agevar' if `agevar' >= `cutoff' & `agevar' >= `lb' & `agevar' <= `ub', ///
        lcolor(black) lpattern(shortdash)) ///
    (scatter dropAfter16 `agevar' if `agevar' >= `lb' & `agevar' <= `ub', ///
        xline(`cutoff', lcolor(black)) ///
        mcolor(black) msymbol(diamond_hollow) msize(small))

graph save "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/thegap_`agevar'.gph", replace
graph export "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/thegap_`agevar'.pdf", replace




restore
 

 
 preserve

local pretreatments first_election  white nouk mumworked14 dadworked14 mumdead14 daddead14 mumaway14 dadaway14 
 local agevar year_birth_run

collapse `pretreatments', by(`agevar')

foreach pretreat of local pretreatments{
	
	sum `agevar'
	local sd = r(sd)
	local n = r(N)^-0.2
	
	local h= 1.5*`n'*`sd'
	display `h'
	
	qui sum `pretreat' if abs(`agevar')<=10
	local ylow = r(mean) - 1.96*r(sd)
	local yhigh = r(max) + 1.96*r(sd)
	
	 
	llrplot  `pretreat' `agevar', cutoff(0) bw(2) xrange(10) yrange(`ylow' `yhigh')

graph save Graph "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/graph_`pretreat'_`agevar'.gph", replace

graph export "/Users/user/Dropbox/Econometrics/educationLiberalism//outcomes/graph_`pretreat'_`agevar'.pdf", replace

	
}

restore

preserve

collapse dropAt15 custom_qtr , by(q_birth_run)

twoway ///
    (lpoly d q if (q) < 0 & q>=-50, lcolor(black) lpattern(solid)) ///
	(lpoly d q if (q) >= 0 & q<=50, lcolor(black) lpattern(solid)) ///
    (scatter d q if abs(q) <= 50 & custom == 2, ///
        mcolor(black) msymbol(Oh) msize(big)) ///
    (scatter d q if abs(q) <= 50 & custom != 2, ///
        mcolor(gs8) msymbol(X) msize(big)), ///
    legend(order(2 "June-August" 3 "September to May") ///
           ring(0) pos(1) col(1)) ///
    graphregion(color(white)) ///
    plotregion(color(white)) ///
    xtitle("Quarters from 9/1957") ytitle("Proportion dropped at age 15")
	
graph save Graph "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/invalidIv.gph", replace

graph export "/Users/user/Dropbox/Econometrics/educationLiberalism//outcomes/invalidIv.pdf", replace


restore 




histogram month_birth_run if inrange(month_birth_run, -10, 10), ///
    discrete width(0.5)  ///
    ylabel(, angle(horizontal)) fcolor(black) lcolor(black) ///
    xlabel(-10(1)10) ///
	xtitle("Months to/from September 1957") ///
    graphregion(color(white)) bgcolor(white)

graph save Graph "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/distribution_age_month.gph", replace

graph export "/Users/user/Dropbox/Econometrics/educationLiberalism//outcomes/distribution_age_month.pdf", replace


histogram q_birth_run if inrange(q_birth_run, -10, 10), ///
    discrete width(.5)  ///
     ylabel(, angle(horizontal)) fcolor(black) lcolor(black) ///
    xlabel(-10(1)10) ///
	xtitle("Quarters to/from September 1957") ///
    graphregion(color(white)) bgcolor(white)
graph save Graph "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/distribution_age_quarter.gph", replace

graph export "/Users/user/Dropbox/Econometrics/educationLiberalism//outcomes/distribution_age_quarter.pdf", replace

	
histogram year_birth_run if inrange(year_birth_run, -10, 10), ///
    discrete width(.5)  ///
     ylabel(, angle(horizontal)) fcolor(black) lcolor(black) ///
    xlabel(-10(1)10) ///
	xtitle("Years to/from September 1957") ///
    graphregion(color(white)) bgcolor(white)

graph save Graph "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/distribution_age_year.gph", replace

graph export "/Users/user/Dropbox/Econometrics/educationLiberalism//outcomes/distribution_age_year.pdf", replace

	
	 

**# Q3 guys are different.

* Define the young dummy
gen young = (custom_qtr==2)

* Variables to process
local vars cg_num_dis cg_vrb_dis cg_mem_dis loginc everunemployed

* Create an empty matrix: rows = number of vars, cols = 4
matrix results = J(`: word count `vars'', 4, .)

local i = 1

foreach v of local vars {
    
    quietly ttest `v' if dropAt15==1, by(young)

    * Extract results
    local mean0 = r(mu_1)
    local se0   = r(se_1)

    local mean1 = r(mu_2)
    local se1   = r(se_2)

    local diff  = r(mu_2) - r(mu_1)
    local pval  = r(p)

    * Fill the matrix
    matrix results[`i',1] = `mean0'
    matrix results[`i',2] = `mean1'
    matrix results[`i',3] = `diff'
    matrix results[`i',4] = `pval'

    * Display formatted line with SEs in parentheses
    di as txt "`v' : " ///
        %6.3f `mean0' " (" %6.3f `se0' ")   " ///
        %6.3f `mean1' " (" %6.3f `se1' ")   " ///
        %6.3f `diff' "   p=" %6.4f `pval'

    local ++i
}

* Name matrix columns
matrix colnames results = mean0 mean1 diff pvalue

* Display matrix
matlist results, format(%8.3f)




* Variables to process (Unconditional)
local vars cg_num_dis cg_vrb_dis cg_mem_dis loginc everunemployed

* Create an empty matrix: rows = number of vars, cols = 4
matrix results = J(`: word count `vars'', 4, .)

local i = 1

foreach v of local vars {
    
    quietly ttest `v' , by(young)

    * Extract results
    local mean0 = r(mu_1)
    local se0   = r(se_1)

    local mean1 = r(mu_2)
    local se1   = r(se_2)

    local diff  = r(mu_2) - r(mu_1)
    local pval  = r(p)

    * Fill the matrix
    matrix results[`i',1] = `mean0'
    matrix results[`i',2] = `mean1'
    matrix results[`i',3] = `diff'
    matrix results[`i',4] = `pval'

    * Display formatted line with SEs in parentheses
    di as txt "`v' : " ///
        %6.3f `mean0' " (" %6.3f `se0' ")   " ///
        %6.3f `mean1' " (" %6.3f `se1' ")   " ///
        %6.3f `diff' "   p=" %6.4f `pval'

    local ++i
}

* Name matrix columns
matrix colnames results = mean0 mean1 diff pvalue

* Display matrix
matlist results, format(%8.3f)


* Make sure estout is installed
ssc install estout, replace

* Store regression results
eststo clear

reg cg_num_dis i.young##(dropAt15 dropAt16 dropAfter16), vce(robust)
eststo num

reg cg_vrb_dis i.young##(dropAt15 dropAt16 dropAfter16), vce(robust)
eststo vrb

reg cg_mem_dis i.young##(dropAt15 dropAt16 dropAfter16), vce(robust)
eststo mem

reg loginc i.young##(dropAt15 dropAt16 dropAfter16), vce(robust)
eststo inc

reg everunemployed i.young##(dropAt15 dropAt16 dropAfter16), vce(robust)
eststo unemploy

reg gceGrade25 i.young##(dropAt15 dropAt16 dropAfter16), vce(robust)
eststo gce25

reg gceO i.young##(dropAt15 dropAt16 dropAfter16), vce(robust)
eststo gceO

reg noQual i.young##(dropAt15 dropAt16 dropAfter16), vce(robust)
eststo noQual



esttab num vrb mem inc unemploy gce25 gceO noQual using "/Users/user/Dropbox/Econometrics/educationLiberalism//outcomes/results.tex", replace ///
    label se star(* 0.10 ** 0.05 *** 0.01) ///
    title("Effects of Being Young by Dropout Group") ///
    mtitles("Numerical" "Verbal" "Memory" "Log-Income" "Ever Unemployed") ///
	///
    compress

**# Bunching tests (Cattaneo, Jansson and Ma)

rddensity month_birth_run, c(0)  all noreg plot
rddensity q_birth_run, c(0)  all noreg
rddensity year_birth_run , c(0)  all plot  noreg 
	

**# NCDS


* Fomr the NCDS58
use "/Users/user/Documents/datasets/NCDS1958/UKDA-5565-stata/stata/stata13/ncds0123.dta", clear




* ==============================================================
* NCDS Age-16 (1974) – Socioeconomic background composite
* Recommended 8-component version (the "gold standard")
* ==============================================================


* --------------------------------------------------------------
* 1. Father's social class (Registrar General 1970) – primary
* --------------------------------------------------------------
gen father_sc = n2384
recode father_sc (-9/-1 = .)   // missing codes in NCDS
label var father_sc "Father's RG social class 1970 (1=I, 2=II, 3=IIINM, 4=IIIM, 5=IV, 6=V)"

* If father absent/missing → use mother's class if she works
gen father_sc_full = father_sc
replace father_sc_full = n2393 if n2393 > 0 & n2393 < . & (father_sc >= . | father_sc <= 0)
label var father_sc_full "Father/male head SC, or mother if father missing"

* --------------------------------------------------------------
* 2. Parental education (age left full-time education)
* --------------------------------------------------------------
recode n2396 n2397 (-9/-1 = .), gen(age_dad_left_educ age_mum_left_educ)
gen par_ed_max = max(age_dad_left_educ, age_mum_left_educ)
recode par_ed_max (min/15=1) (16=2) (17/18=3) (19/max=4) (.=.), ///
       gen(par_ed_cat)
label define pared 1 "≤15" 2 "16" 3 "17-18" 4 "≥19"
label values par_ed_cat pared
label var par_ed_cat "Highest parental age leaving full-time education (categories)"

* --------------------------------------------------------------
* 3. Housing tenure
* --------------------------------------------------------------
gen tenure = n2470
recode tenure (-9/-1=.) (1/2=1) (3=2) (4/5=3), gen(tenure3)
label define tenure3 1 "Owner-occupied" 2 "Council rented" 3 "Private rented/other"
label values tenure3 tenure3
label var tenure3 "Housing tenure (3 cats)"

* --------------------------------------------------------------
* 4. Crowding (persons per room)
* --------------------------------------------------------------
gen crowding = n1734   // version A, most commonly used
replace crowding = n1735 if crowding >= . & n1735 < .
recode crowding (min/1=1) (1.01/1.5=2) (1.51/max=3) (.=.), gen(crowd3)
label define crowd3 1 "≤1 ppr" 2 "1.01-1.5" 3 ">1.5"
label values crowd3 crowd3
label var crowd3 "Persons per room (3 categories)"

* --------------------------------------------------------------
* 5. Access to household amenities (derived variable – best one)
* --------------------------------------------------------------
gen amenities = n1736
recode amenities (-9/-1=.) (0/2=1) (3=2) (4=3), gen(amenities3)
label define amen3 1 "None/shared" 2 "Some" 3 "All three (bath, indoor WC, hot water)"
label values amenities3 amen3
label var amenities3 "Access to basic amenities (0-3)"

* --------------------------------------------------------------
* 6. Free school meals (strong poverty marker)
* --------------------------------------------------------------
gen fsm = (n2440 == 1) if n2440 >= 0 & n2440 < .
label var fsm "Any child in household receives free school meals (1974)"

* --------------------------------------------------------------
* 7. Financial difficulties in past year
* --------------------------------------------------------------
gen finhard = (n2441 == 1) if n2441 >= 0 & n2441 < .
label var finhard "Serious financial difficulties in past 12 months"

* --------------------------------------------------------------
* 8. Car ownership (very strong class marker in 1970s Britain)
*    (not directly listed by you but almost always used – usually derived)
* --------------------------------------------------------------
capture confirm variable n2476   // sometimes in derived files
if _rc == 0 {
    gen car = (n2476 >= 1 & n2476 < .) if n2476 >= 0
}
else {
    gen car = .
}
label var car "Household has use of a car/van (1974)"

* ==============================================================
* Final composite scores
* ==============================================================

* A. Simple additive deprivation index (0–7, higher = more disadvantaged)
gen ses_deprivation = 0
replace ses_deprivation = ses_deprivation + 1 if inlist(father_sc_full,4,5) | father_sc_full==6
replace ses_deprivation = ses_deprivation + 1 if par_ed_cat == 1
replace ses_deprivation = ses_deprivation + 1 if tenure3 == 2 | tenure3 == 3
replace ses_deprivation = ses_deprivation + 1 if crowd3 >= 2
replace ses_deprivation = ses_deprivation + 1 if amenities3 < 3
replace ses_deprivation = ses_deprivation + 1 if fsm == 1
replace ses_deprivation = ses_deprivation + 1 if finhard == 1
* replace ses_deprivation = ses_deprivation + 1 if car == 0   // uncomment if you have car

label var ses_deprivation "NCDS age-16 SES deprivation index (0-7 or 0-8)"

* B. Quintiles of the deprivation index (most common in publications)
xtile ses_quintile = ses_deprivation, nq(5)
label var ses_quintile "SES deprivation quintiles (1=least deprived, 5=most deprived)"

* C. Traditional 3-class scheme (very widely used)
gen social_class_3 = .
replace social_class_3 = 1 if inlist(father_sc_full,1,2) | (father_sc_full==3 & n2385<=5)  // I, II, IIINM professional/managerial
replace social_class_3 = 2 if inlist(father_sc_full,3) & !inlist(n2385,1,2,3,4,5)          // IIINM skilled non-manual (intermediate)
replace social_class_3 = 3 if inlist(father_sc_full,4,5,6)                                   // IIIM, IV, V
replace social_class_3 = 3 if social_class_3==. & tenure3==2 & fsm==1   // push remaining council + FSM into working class
label define sc3 1 "Middle class" 2 "Intermediate" 3 "Working class"
label values social_class_3 sc3
label var social_class_3 "NCDS standard 3-category social class at 16"

* --------------------------------------------------------------
* Final labels and notes
* --------------------------------------------------------------
notes ses_deprivation: NCDS 1974 (age 16) 7-item deprivation index (father SC manual, low par ed, rented tenure, crowding>1, missing amenities, FSM, financial hardship)
notes social_class_3: Standard 3-class scheme used in most NCDS social mobility papers






merge 1:1 ncdsid using "/Users/user/Documents/datasets/NCDS1958/UKDA-5566-stata/stata/stata13/ncds4.dta", keepusing(n5957 n5959 n5960 n5962 n5964 agelsch)

drop _merge

merge 1:1 ncdsid using "/Users/user/Documents/datasets/NCDS1958/UKDA-5567-stata/stata/stata13_se/ncds5cmi.dta", keepusing(n509514 n509516 n509523 n509527 n509528 n509556 n509560 n509562 n509563 n509565 n509567 n509569 n509663 n509667 n509670 n509714 n509715 n509752 n509753 n509755 n509756 n504635  n504636 n504638 n504640 n504641 n504642 n504644)

drop _merge

merge 1:1 ncdsid using "/Users/user/Documents/datasets/NCDS1958/UKDA-5578-stata/stata/stata13_se/ncds6_v2.dta", keepusing(vote97 votewho) norepor


gen complier = (n2741 == 1 & n2729 == 1)
replace complier =. if n2741 <0 | n2729 <0

gen complierB = (n2741 == 1 & agelsch == 1)
replace complierB =. if n2741 <0 | agelsch<0


gen votedIn79 =  n5959 == 1
replace votedIn79 =. if  n5959 >=8

gen votedLabour79 = n5960 == 2
replace votedLabour79 =. if n5960>=11

gen votedConservative79 = n5960 == 1
replace votedConservative79 =. if n5960>=11

gen pastUnionistWave4 = 1 if  n5964 == 1.
replace pastUnionistWave4 =. if n5964 >2


gen votedIn87 =  n504635 == 1
replace votedIn87 =. if  n504635==. 

gen votedLabour87 = n504636 == 2
replace votedLabour87 =. if n504636 >8

gen votedConservative87 = n504636 == 1
replace votedConservative87 =. if n504636 >8


gen votedIn97 =  vote97 == 1
replace votedIn97 =. if  vote97 >=8

gen votedLabour97 = votewho == 2
replace votedLabour97 =. if votewho >8

gen votedConservative97 = votewho == 1
replace votedConservative97 =. if votewho >8



* Views in 1991
local opinions n509514 n509516 n509523 n509527 n509528 n509556 n509560 n509562 n509563 n509565 n509567 n509569 n509663 n509667 n509670 n509714 n509715 votedIn79  votedLabour79 votedConservative79 votedIn87 votedLabour87  votedConservative87 votedIn97 votedLabour97 votedConservative97 ses_deprivation social_class_3

local n : word count `opinions'
display `n'
matrix mRES = J(2*`n', 4, .)

local i = 1
foreach j of local opinions {
	
	ttest `j' if `j' !=. , by(complierB)
	
	matrix mRES[`i', 1 ] = r(mu_1)
	matrix mRES[`i', 2 ] = r(mu_2)
	matrix mRES[`i', 3 ] = r(mu_2) -r(mu_1)
	matrix mRES[`i', 4 ] = r(p)
	
	matrix mRES[`i'+1, 1 ] = r(se_1)
	matrix mRES[`i'+1, 2 ] = r(se_2)


	local i =`i' +2
	
}

 

matrix rownames mRES = "at_expense_of_workers" "." "no_say_in_what_Govt_does" "." "Env_problems_notas_serious" "." "Law_should_be_obeyed_if_wrong" "." "Environment_vs_econ_growth" "." "Ordinary_people_fair_share" "." "Death_penalty_for_some_crimes" "." "Party_in_power_no_difference" "." "Women_right_to_choose_abortion" "." "Everyone_needs_private_health" "." "Environment_is_top_issue" "." "Private_schools_abolished" "." "Law_rich_vs_poor" "." "Politicians_self_interest" "." "Youth_no_respect_trad_values" "." "Govt_redistribute_income" "." "No_party_would_help_me" "." "votedIn79" "." "votedLabour79" "." "votedConservative79" "." "votedIn87" "." "votedLabour87" "." "votedConservative87" "." "votedIn97" "." "votedLabour97" "." "votedConservative97" "." "ses_deprivation" "."   "social_class_3" "."

matrix list mRES, format(%9.3f)

outtable using "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/ncds_views_B", mat(mRES) replace nobox



 
**# 1933 Reform




use "/Users/user/Dropbox/Econometrics/educationLiberalism/code/all_waves_protect.dta", clear

* ========================= Part 3: Analysis =========================

encode wave, gen(wave_num)
xtset pidp wave_num

replace dobm_dv = dobm if dobm_dv ==.

by pidp: gen doby_fwd = doby_fix[1] if !missing(doby_fix[1])
by pidp: replace doby_fwd = doby_fix if missing(doby_fwd) & !missing(doby_fix)

by pidp: gen dobm_fwd = dobm_dv[1] if !missing(dobm_dv[1])
by pidp: replace dobm_fwd = dobm_dv if missing(dobm_fwd) & !missing(dobm_dv)

by pidp (wave_num): gen doby_bwd = doby_fix[_N] if !missing(doby_fix[_N])
by pidp (wave_num): replace doby_bwd = doby_fix if missing(doby_bwd) & !missing(doby_fix)

by pidp (wave_num): gen dobm_bwd = dobm_dv[_N] if !missing(dobm_dv[_N])
by pidp (wave_num): replace dobm_bwd = dobm_dv if missing(dobm_bwd) & !missing(dobm_dv)

gen doby_clean = doby_fix
replace doby_clean = doby_fwd if missing(doby_clean) & !missing(doby_fwd)
replace doby_clean = doby_bwd if missing(doby_clean) & !missing(doby_bwd)
gen dobm_clean = dobm_dv
replace dobm_clean = dobm_fwd if missing(dobm_clean) & !missing(dobm_fwd)
replace dobm_clean = dobm_bwd if missing(dobm_clean) & !missing(dobm_bwd)
drop doby_fwd doby_bwd dobm_fwd dobm_bwd



 
gen vulnerable = dadworked14 ==0 | mumdead14==1 | daddead14==1| mumaway14 ==1 | dadaway14==1  | mumworked14 ==0


* Ordinary people share nations wealth
gen aux = opsoca == 1 | opsoca == 2
replace aux = . if opsoca < 0 | opsoca == .
bysort pidp: egen shareWealth = max(aux)
drop aux

* One law for rich and one for poor
gen aux = opsocb == 1 | opsocb == 2
replace aux = . if opsocb < 0 | opsocb == .
bysort pidp: egen oneLaw = max(aux)
drop aux

* Private enterprise solves economic probs
gen aux = opsocc == 1 | opsocc == 2
replace aux = . if opsocc < 0 | opsocc == .
bysort pidp: egen privateEnterprise = max(aux)
drop aux

* Public services ought to be state owned
gen aux = opsocd == 1 | opsocd == 2
replace aux = . if opsocd < 0 | opsocd == .
bysort pidp: egen stateOwned = max(aux)
drop aux

* Govt. has obligation to provide jobs
gen aux = opsoce == 1 | opsoce == 2
replace aux = . if opsoce < 0 | opsoce == .
bysort pidp: egen provideJobs = max(aux)
drop aux

* Strong trade unions protect employees
gen aux = opsocf == 1 | opsocf == 2
replace aux = . if opsocf < 0 | opsocf == .
bysort pidp: egen strongUnions = max(aux)
drop aux

* Young people don't respect British values
gen aux = opsock == 1 | opsock == 2
replace aux = . if opsock < 0 | opsock == .
bysort pidp: egen britishValues = max(aux)
drop aux

* Censorship necessary to uphold moral standards
gen aux = opsocl == 1 | opsocl == 2
replace aux = . if opsocl < 0 | opsocl == .
bysort pidp: egen censorship = max(aux)
drop aux

* Public meetings protest against government allowed
gen aux = opsocm == 1 | opsocm == 2
replace aux = . if opsocm < 0 | opsocm == .
bysort pidp: egen protestAllowed = max(aux)
drop aux

* More tolerant of unconventional lives
gen aux = opsocn == 1 | opsocn == 2
replace aux = . if opsocn < 0 | opsocn == .
bysort pidp: egen tolerantLives = max(aux)
drop aux

* Death penalty appropriate some crimes
gen aux = opsoco == 1 | opsoco == 2
replace aux = . if opsoco < 0 | opsoco == .
bysort pidp: egen deathPenalty = max(aux)
drop aux

* Stiffer sentences
gen aux = opsocp == 1 | opsocp == 2
replace aux = . if opsocp < 0 | opsocp == .
bysort pidp: egen stiffSentences = max(aux)
drop aux


* Qualifications in the US panel

gen aux =.
replace aux = qfhigh if qfhigh>0 &qfhigh!=.
bysort pidp: egen qualx = max(aux)
drop aux

* qualifications in the BHPS
bysort pidp:  egen qualReceived = max(qfedhi)
replace qualReceived=. if qualReceived<0


gen noqualifications = qualReceived==12 | qualx == 96
gen gce25 = qualReceived ==9 | qualx == 13
gen gceO = qualReceived ==7 | qualx == 12
gen degree= qualReceived== 1 | qualReceived==2 | qualx == 1 |qualx == 2


bysort pidp: egen noQual = max(noqualifications)
bysort pidp: egen gceGrade25 = max(gce25)
bysort pidp: egen gceOlevel = max(gceO)
bysort pidp: egen uniDegree = max(degree)



collapse cg_*_dis  aveinc female dropAt14 dropAt16 dropAt15 dropAfter noInterestedPolitics mostly* left15 left16 doby_fix nouk white mumworked14 dadworked14 mumdead14 daddead14 mumaway14 dadaway14 everunemployed loginc first_election   dob*_clean voted_* vulnerable shareWealth oneLaw privateEnterprise stateOwned provideJobs strongUnions britishValues censorship protestAllowed tolerantLives deathPenalty stiffSentences gce* noQual uniDegree onAverageLeft onAverageRight , by(pidp)

replace dobm_clean = floor(dobm_clean)

drop if dobm_clean ==. | doby_clean ==.  | dobm_clean <0 | doby_clean <0
 
/*

*** Custom quarters
gen custom_qtr = .
replace custom_qtr = 0 if inlist(dobm_clean,12,1,2)   // Q1 = Dec–Feb
replace custom_qtr = 1 if inlist(dobm_clean,3,4,5)    // Q2 = Mar–May
replace custom_qtr = 2 if inlist(dobm_clean,6,7,8)    // Q3 = Jun–Aug
replace custom_qtr = 3 if inlist(dobm_clean,9,10,11)  // Q4 = Sep–Nov

local base_year = 1933
local base_qtr  = 3        // September 1933 falls in Q4=3

gen aux = dobm_clean == 12

gen q_birth_run = (doby_clean - `base_year')*4 + (custom_qtr - `base_qtr') + 4*aux


*** Bimesters
gen custom_bim = .
replace custom_bim = 0 if inlist(dobm_clean,1,2)
replace custom_bim = 1 if inlist(dobm_clean,3,4)
replace custom_bim = 2 if inlist(dobm_clean,5,6)
replace custom_bim = 3 if inlist(dobm_clean,7,8)
replace custom_bim = 4 if inlist(dobm_clean,9,10)     // Sep–Oct
replace custom_bim = 5 if inlist(dobm_clean,11,12)    // Nov–Dec

local base_year = 1933
local base_bim  = 4

gen bim_birth_run = (doby_clean - `base_year')*6 + (custom_bim - `base_bim')

 
 
 *** Custom 12-month periods starting in September
gen custom_12m = .
replace custom_12m = 1  if dobm_clean == 9
replace custom_12m = 2  if dobm_clean == 10
replace custom_12m = 3  if dobm_clean == 11
replace custom_12m = 4  if dobm_clean == 12
replace custom_12m = 5  if dobm_clean == 1
replace custom_12m = 6  if dobm_clean == 2
replace custom_12m = 7  if dobm_clean == 3
replace custom_12m = 8  if dobm_clean == 4
replace custom_12m = 9  if dobm_clean == 5
replace custom_12m = 10 if dobm_clean == 6
replace custom_12m = 11 if dobm_clean == 7
replace custom_12m = 12 if dobm_clean == 8

local base_year   = 1933
local base_period = 1       // September = period 1

gen month_birth_run = (doby_clean - `base_year')*12 + (custom_12m - `base_period')

gen year_birth_run = doby_clean - 1933
replace year_birth_run = year_birth_run - 1 if dobm_clean < 9
gen wrong_run = doby_fix - 1933

*/




*** Custom quarters
gen custom_qtr = .
replace custom_qtr = 0 if inlist(dobm_clean,1,2,3)   // Q1 = Dec–Feb
replace custom_qtr = 1 if inlist(dobm_clean,4,5, 6)    // Q2 = Mar–May
replace custom_qtr = 2 if inlist(dobm_clean,7,8,9)    // Q3 = Jun–Aug
replace custom_qtr = 3 if inlist(dobm_clean,10,11,12)  // Q4 = Sep–Nov

local base_year = 1933
local base_qtr  = 1        // April 1933 is in Q2 (=1)

gen aux = dobm_clean == 12

gen q_birth_run = (doby_clean - `base_year')*4 + (custom_qtr - `base_qtr') + 4*aux

local base_bim = 1
 
*** Custom 12-month periods starting in September
* Mapping months to 12 periods, Sep = 1, Oct = 2, ..., Aug = 12
gen custom_12m = .
replace custom_12m = 1  if dobm_clean == 9
replace custom_12m = 2  if dobm_clean == 10
replace custom_12m = 3  if dobm_clean == 11
replace custom_12m = 4  if dobm_clean == 12
replace custom_12m = 5  if dobm_clean == 1
replace custom_12m = 6  if dobm_clean == 2
replace custom_12m = 7  if dobm_clean == 3
replace custom_12m = 8  if dobm_clean == 4
replace custom_12m = 9  if dobm_clean == 5
replace custom_12m = 10 if dobm_clean == 6
replace custom_12m = 11 if dobm_clean == 7
replace custom_12m = 12 if dobm_clean == 8

local base_year = 1933
local base_period = 8     // April 1933 = period 8

* Handle December wrap-around as in quarter code
 
* Running variable
gen month_birth_run = (doby_clean - `base_year')*12 + (custom_12m - `base_period') + 12*aux
gen year_birth_run = doby_clean - 1933
replace year_birth_run = year_birth_run - 1 if dobm_clean < 4




local outcomes  mostlyLeft mostlyRight mostlyCentrist noInterestedPolitics onAverageLeft onAverageRight

* List of running variables and corresponding matrix names
local runs month_birth_run q_birth_run year_birth_run 
local num_outcomes : word count `outcomes'

foreach running of local runs {
	
	matrix mRES = J(`num_outcomes', 10, .)
	local i = 0
	foreach y of local outcomes{
		local ++i
		rdrobust `y' `running',  bwselect(msetwo)
		matrix mRES[`i',1] = e(tau_bc)
        matrix mRES[`i',2] = e(se_tau_rb)
        matrix mRES[`i',3] = e(pv_rb)
        //matrix mRES[`i',4] = e(ci_l_rb)
       // matrix mRES[`i',5] = e(ci_r_rb)
		rdrobust `y' `running', fuzzy(left15) bwselect(msetwo)
		matrix mRES[`i',6] = e(tau_bc)
        matrix mRES[`i',7] = e(se_tau_rb)
        matrix mRES[`i',8] = e(pv_rb)
       // matrix mRES[`i',9] = e(ci_l_rb)
       // matrix mRES[`i',10] = e(ci_r_rb)
		
	}
	matrix colnames mRES = tau_bc se_rb p_rb ci_l ci_r tau_bc se_rb p_rb ci_l ci_r
    matrix rownames mRES = `outcomes'
	outtable using "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/`running'_1933", mat(mRES) replace nobox
}






/*
rdwinselect month_birth_run white nouk mumworked14 dadworked14 mumdead14 daddead14 mumaway14 dadaway14  first , wmass level(0.05)  wasym  nw(10)
Missing values detected in covariates
Consider dropmissing option to exclude missing values
Mass points detected in running variable
You may use wmasspoints option for constructing windows at each mass point


Window selection for RD under local randomization


Cutoff c = 0.00   | Left of c   Right of c        Number of obs  =        108606
------------------+-----------------------        Order of poly  =             0
    Number of obs |     34765        73841        Kernel type    =       uniform
   1st percentile |       364          857        Reps           =          1000
   5th percentile |      1744         3761        Testing method =     rdrandinf
  10th percentile |      3560         7527        Balance test   =     diffmeans
  20th percentile |      6980        14855


                  |   Bal. test         Var. name    Bin. test 
      Window      |    p-value        (min p-value)   p-value     Obs<c   Obs>=c
------------------+-------------------------------------------------------------
  -1.000|   0.000 |      0.252             nouk         0.353        64       76
  -2.000|   1.000 |      0.300        mumaway14         0.390       144      160
  -3.000|   2.000 |      0.162        mumdead14         0.745       228      236
  -4.000|   3.000 |      0.174        mumdead14         1.000       303      303
  -5.000|   4.000 |      0.222        mumdead14         0.642       384      398
  -6.000|   5.000 |      0.214        mumdead14         0.672       463      477
  -7.000|   6.000 |      0.300        mumdead14         0.382       535      565
  -8.000|   7.000 |      0.214        mumdead14         0.325       614      650
  -9.000|   8.000 |      0.000       first_election     0.232       686      732
 -10.000|   9.000 |      0.000       first_election     0.122       746      808
 
Variable used in binomial test (running variable): month_birth_run
Covariates used in balance test: white nouk mumworked14 dadworked14 mumdead14 da
> ddead14 mumaway14 dadaway14 first_election

Recommended window is [-8.000;  7.000] with 1264 observations (614 below, 650 above).


. rdwinselect q_birth_run white nouk mumworked14 dadworked14 mumdead14 daddead14 mumaway14 dadaway14  first , wmass level(0.05)  wasym  nw(20)
Missing values detected in covariates
Consider dropmissing option to exclude missing values
Mass points detected in running variable
You may use wmasspoints option for constructing windows at each mass point


Window selection for RD under local randomization


Cutoff c = 0.00   | Left of c   Right of c        Number of obs  =        108606
------------------+-----------------------        Order of poly  =             0
    Number of obs |     34765        73841        Kernel type    =       uniform
   1st percentile |       364         1118        Reps           =          1000
   5th percentile |      1744         3761        Testing method =     rdrandinf
  10th percentile |      3790         7527        Balance test   =     diffmeans
  20th percentile |      6980        15014


                  |   Bal. test         Var. name    Bin. test 
      Window      |    p-value        (min p-value)   p-value     Obs<c   Obs>=c
------------------+-------------------------------------------------------------
  -1.000|   0.000 |      0.222        mumdead14         0.745       228      236
  -2.000|   1.000 |      0.232        mumdead14         0.672       463      477
  -3.000|   2.000 |      0.000       first_election     0.232       686      732
  -4.000|   3.000 |      0.000       first_election     0.136       889      954
  -5.000|   4.000 |      0.000       first_election     0.083      1105     1189
  -6.000|   5.000 |      0.000       first_election     0.290      1374     1431
  -7.000|   6.000 |      0.000       first_election     0.125      1598     1687
  -8.000|   7.000 |      0.000       first_election     0.012      1797     1951
  -9.000|   8.000 |      0.000       first_election     0.007      1999     2175
 -10.000|   9.000 |      0.000       first_election     0.002      2237     2452


Variable used in binomial test (running variable): q_birth_run
Covariates used in balance test: white nouk mumworked14 dadworked14 mumdead14 da
> ddead14 mumaway14 dadaway14 first_election

Recommended window is [-2.000;  1.000] with 940 observations (463 below, 477 above).




*/
local outcomes  mostlyLeft mostlyRight noInterestedPolitics

* List of running variables and corresponding bandwidths
local runs month_birth_run q_birth_run year_birth_run 
local wl_list -8 -2 -1 -1
local wr_list 7 1 0.5 0.5

local num_outcomes : word count `outcomes'

* Loop over running variables
forvalues j = 1/`=wordcount("`runs'")' {
    local running : word `j' of `runs'
    local wl : word `j' of `wl_list'
    local wr : word `j' of `wr_list'
    
    matrix mRES = J(`num_outcomes', 10, .)
    local i = 0
    
    foreach y of local outcomes {
        local ++i
        
        * ITT
        rdrandinf `y' `running', wl(`wl') wr(`wr') //ci(0.05)
        matrix V = .
        matrix CI = r(CI)
        matrix mRES[`i',1] = r(obs_stat)
        matrix mRES[`i',2] = .
        matrix mRES[`i',3] = r(asy_pval)
        //matrix mRES[`i',4] = CI[1,1]
        //matrix mRES[`i',5] = CI[1,2]
        
        * LATE (fuzzy)
        rdrandinf `y' `running', fuzzy(left15 tsls) wl(`wl') wr(`wr')  //ci(0.05)
        matrix V = e(V)
        matrix CI = r(CI)
        matrix mRES[`i',6] = r(obs_stat)
        matrix mRES[`i',7] = V[1,1]
        matrix mRES[`i',8] = r(asy_pval)
        //matrix mRES[`i',9] = CI[1,1]
        //matrix mRES[`i',10] = CI[1,2]
    }
    
    matrix colnames mRES = tau_bc se_rb p_rb ci_l ci_r tau_bc se_rb p_rb ci_l ci_r
    matrix rownames mRES = `outcomes'
    
    outtable using "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/ri_`running'_1933", mat(mRES) replace nobox
}





 preserve
 local agevar q_birth_run

collapse  dropAt14 dropAt15 dropAt16, by(`agevar')

if "`agevar'" == "year_birth_run" {
    local xtag "Years from cut-off"
}
else if "`agevar'" == "q_birth_run" {
    local xtag "Quarters from cut-off"
}
else if "`agevar'" == "wrong_run"{
	local xtag "Year of birth"
}
else {
    local xtag "Months from cut-off"
}
local lb = -24
local ub = 24
local cutoff = 0

twoway ///
    (lpoly dropAt14 `agevar' if `agevar' < `cutoff' & `agevar' >= `lb' & `agevar' <= `ub', ///
        lcolor(black) lpattern(.#.) lwidth(0.5)) ///
    (lpoly dropAt14 `agevar' if `agevar' >= `cutoff' & `agevar' >= `lb' & `agevar' <= `ub', ///
        lcolor(black) lpattern(.#.) lwidth(0.5)) ///
    (scatter dropAt14 `agevar' if `agevar' >= `lb' & `agevar' <= `ub', ///
        mcolor(gray) msymbol(Oh) msize(small) ///
        xtitle(`xtag') ytitle(Percentage) ///
        graphregion(fcolor(white) lcolor(white)) ///
        plotregion(fcolor(white) lcolor(white)) ///
        legend(off) xlabel(, angle(45)) ///
        xscale(range(`lb' `ub'))) ///
    (lpoly dropAt15 `agevar' if `agevar' < `cutoff' & `agevar' >= `lb' & `agevar' <= `ub', ///
        lcolor(black)) ///
    (lpoly dropAt15 `agevar' if `agevar' >= `cutoff' & `agevar' >= `lb' & `agevar' <= `ub', ///
        lcolor(black)) ///
    (scatter dropAt15 `agevar' if `agevar' >= `lb' & `agevar' <= `ub', ///
        xline(`cutoff', lcolor(black)) ///
        mcolor(gray) msymbol(O) msize(small)) ///
    (lpoly dropAt16 `agevar' if `agevar' < `cutoff' & `agevar' >= `lb' & `agevar' <= `ub', ///
        lcolor(black) lpattern(shortdash)) ///
    (lpoly dropAt16 `agevar' if `agevar' >= `cutoff' & `agevar' >= `lb' & `agevar' <= `ub', ///
        lcolor(black) lpattern(shortdash)) ///
    (scatter dropAt16 `agevar' if `agevar' >= `lb' & `agevar' <= `ub', ///
        xline(`cutoff', lcolor(black)) ///
        mcolor(black) msymbol(diamond_hollow) msize(small))

graph save "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/thegap_`agevar'_1933.gph", replace
graph export "/Users/user/Dropbox/Econometrics/educationLiberalism/outcomes/thegap_`agevar'_1933.pdf", replace
