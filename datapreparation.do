/*
	Education and Political Alignment
	Eduardo Fe
	University of Manchester.
	
	This version: Version 1.0. 
	Date: 21 June 2024
	
	Part 1: Functions
	Part 2: Data Preparation
	Part 3: Analysis. 
	
*/

program drop my_rdrandinf_function
program drop standardize_var
program drop trimMean

program define my_rdrandinf_function, rclass
    args wl wr c fuzzy_var running outcomes
    
	local num_outcomes : word count `outcomes'
    
    // Create the matrix with the correct number of rows
    matrix mRes = J(`num_outcomes', 6, .)
	
    
    local i=1
    foreach j of local outcomes {
        rdrandinf  `j' `running' if `j'!=., fuzzy(`fuzzy_var') c(`c') wl(`wl') wr(`wr') ci(0.05) interfci(0.05)

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
		
        local i = `i' + 1
    }
    
    matrix list mRes
	
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

/*
program define ittbounds
    // Declare the syntax for the program
    syntax varlist(min=3 max= 3) 

	tokenize `varlist'
	tempvar y
	tempvar d
	tempvar z
	
	gen `y'   = `1'
	gen `d'   = `2'
	gen `z'   = `3'
	
	  
	 qui{
		count if `z'==1
		local numZ1 = r(N)
		
		count if `z'==0
		local numZ0 = r(N)
		
		count if `z' == 1 & `d'==1
		local numerator = r(N)
		local pD1Z1 = `numerator' / `numZ1'
		
		count if `z' == 1 & `d'==0
		local numerator = r(N)
		local pD0Z1 = `numerator' / `numZ1'
		
		count if `z' == 0 & `d'==1
		local numerator = r(N)
		local pD1Z0 = `numerator' / `numZ0'
		
		count if `z' == 0 & `d'==0
		local numerator = r(N)
		local pD0Z0 = `numerator' / `numZ0'
		
		local thetavar = (`pD1Z1' - `pD1Z0' / `pD1Z1')
		local thetavarM1 = 1-`thetavar'
		
		local den = 1-`pD1Z0'
		local rhovar =   (`pD1Z1' - `pD1Z0' ) / `den'
		local rhovarM1 = 1-`rhovar'
		
    }
	display "thatavarm1" `thetavarM1'
	display "thetavar " `thetavar'
	
	
	if `pD0Z1' / `pD0Z0' > 1 | `pD1Z0' / `pD1Z1' > 1 {
        display as error "Invalid probabilities: p01/p00 > 1 or p10/p11 > 1"
        exit
    }
	
	trimMean `y'  if `z'==1 & `d'==1, prop(`thetavar')
	local UbYc1 = r(trimmedMean)
	display `UbYc1'
	
	trimMean `y'  if `z'==1 & `d'==1, prop(`thetavarM1') below(1)
	local LbYc1 = r(trimmedMean)
	display "LbYc1 "`LbYc1'
	
	trimMean `y'  if `z'==0 & `d'==0, prop(`rhovar')
	local UbYc0 = r(trimmedMean)
	display `UbYc0'
	
	trimMean `y'  if `z'==0 & `d'==0, prop(`rhovarM1') below(1)
	local LbYc0 = r(trimmedMean)
	display "LbYc0 " `LbYc0'
 
    display "Upper bound: " `UbYc1' - `LbYc0'
	display "Lower bound: " `LbYc1' - `UbYc0'
	
	sum `y' if `z' == 1 & `d' == 1
	local yZ1D1 = r(mean)
	
	sum `y' if `z' == 0 & `d' == 0
	local yZ0D0 = r(mean)
	
	display "Lower bound: " `yZ1D1' - `yZ0D0'
	
	display "Upper bound new ass: " min(`UbYc1',`UbYc0') - max(`LbYc0',`LbYc1')
	
end
*/


*Part 2: Data preparation =========================

// Define the directory where BHPS data files are stored
local data_dir "/Users/user/Documents/datasets/UKDA-6614-stata/stata/stata13_se/bhps/"

// Define a list of waves and the variables you want to extract from each wave
local waves "a b c d e f g h i j k l m n o p q r "
local vars "age_dv pidp hidp sex doby vote* qfedhi hiqualb_dv istrtdaty scend istrtdaty gor_dv age_dv fimnb jbstat f116 f139 maju    paju plbornc   "

// Loop over each wave

local i =1
foreach wave of local waves {
	
	if(`i' >= 8){
		local vars "fimnlabgrs_dv fimngrs_dv  istrtdathh istrtdatd pidp hidp sex doby vote* qfedhi  istrtdaty scend istrtdaty gor_dv age_dv fimnb jbstat f116 f139 maju   paju plbornc  "

	}
	else if(`i' >= 11){
	local vars "fimnlabgrs_dv fimngrs_dv j1soc90 istrtdathh istrtdatd pidp hidp sex doby vote* qfedhi  istrtdaty scend istrtdaty gor_dv age_dv fimnb jbstat f116 f139 maju   paju plbornc  "

	}
	else{
		local vars "fimnlabgrs_dv fimngrs_dv race istrtdathh istrtdatd pidp hidp sex doby vote* qfedhi  istrtdaty scend istrtdaty gor_dv age_dv fimnb jbstat f116 f139 plbornc  "
	}
    // Load the dataset for the current wave
	
    use "`data_dir'/b`wave'_indresp.dta", clear

	rename b`wave'_* *
	
	rename doby doby_dv
	
	
    // Keep only the desired variables
    keep `vars'

    // Add a wave identifier variable
    gen wave = "`wave'"
	
	gen wave_t = `i'

    // Save the modified dataset with a new name
    save "`data_dir'/bhps_wave_`wave'_extracted.dta", replace
	
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
save "`data_dir'/bhps_wave_`wave'_extracted.dta", replace


// to merge all waves into a single dataset
clear
foreach wave of local waves {
    use "`data_dir'/bhps_wave_`wave'_extracted.dta", clear
    if "`wave'" == "a" {
        save "`data_dir'/bhps_all_waves.dta", replace
    }
    else {
        append using "`data_dir'/bhps_all_waves.dta"
        save "`data_dir'/bhps_all_waves.dta", replace
    }
}


* Repeat now for Understanding society

* ------------------------------------------------------------------------------
* ------------------------------------------------------------------------------

// Define the directory where BHPS data files are stored
local data_dir "/Users/user/Documents/datasets/UKDA-6614-stata/stata/stata13_se/ukhls/"


// Define a list of waves and the variables you want to extract from each wave
local waves "a b c d e f g h i j k l "

// Loop over each wave

foreach wave of local waves {
display "`wave'"

	
}

foreach wave of local waves {
    // Load the dataset for the current wave
	
    use "`data_dir'/`wave'_indresp.dta", clear

	rename `wave'_* *
	
    // Keep only the desired variables
    if("`wave'" == "c"){
	display "Eureka!"
	local vars "big5* cg*_dv fimnlabgrs_dv fimngrs_dv j1soc90 age_dv ethn_dv istrtdathh istrtdatd  pidp hidp sex doby_dv vote* jbstat   istrtdaty scend gor_dv age_dv fimnprben_dv fimnsben_dv bendis3 bendis5 maju   paju plbornc   "
	display "`wave'"
	display "`vars'"
	}
	else{

		local vars "fimnlabgrs_dv fimngrs_dv j1soc90  age_dv ethn_dv istrtdathh istrtdatd  pidp hidp sex doby_dv vote* jbstat   istrtdaty scend gor_dv age_dv fimnprben_dv fimnsben_dv bendis3 bendis5 maju   paju plbornc   "
		display "`wave'"
		display "`vars'"
	}
	keep `vars'

    // Add a wave identifier variable
    gen wave = "us_`wave'"
	
	
	gen wave_t = `i'
    // Save the modified dataset with a new name
    save "`data_dir'/us_wave_`wave'_extracted.dta", replace
	
	local i = `i' +1
}

// to merge all waves into a single dataset
clear
foreach wave of local waves {
    use "`data_dir'/us_wave_`wave'_extracted.dta", clear
    if "`wave'" == "a" {
        save "`data_dir'/us_all_waves.dta", replace
    }
    else {
        append using "`data_dir'/us_all_waves.dta"
        save "`data_dir'/us_all_waves.dta", replace
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

local data_dir "/Users/user/Documents/datasets/UKDA-6614-stata/stata/stata13_se/bhps/"

append using "`data_dir'/bhps_all_waves.dta"

save "/Users/user/Dropbox/Econometrics/educationLiberalism/code/all_waves.dta", replace
**# Bookmark #1
 
use "/Users/user/Dropbox/Econometrics/educationLiberalism/code/all_waves.dta", clear
  
 
* School leaving age allocated to all rows within individual
xtset pidp wave_t

by pidp: egen schleavingage = max(scend)


 

gen nopolitical_sr = (vote6 ==4 | vote6 == 3)
replace nopolitical_sr = . if vote6 <0 | vote6==.
label var nopolitical_sr "Not very/at all interested in politics "
by pid: egen noInterestedPolitics = mode(nopolitical_sr)

 
gen labour = vote4== 2 
replace labour =. if vote4 <0 | vote4 ==.
by pid: egen mostlyLabour = mode(labour)

gen conservative = vote4== 1
replace conservative =. if vote4<0| vote4 ==.
by pid: egen mostlyConservative = mode(conser)

gen left = vote4==2| vote4==4 |vote4==5|vote4== 6|vote4== 8|vote4==11
replace left =. if vote4<0| vote4 ==.
by pid: egen mostlyLeft = mode(left)

gen right = vote4==1 |vote4==7|vote4==10|vote4==12|vote4==13|vote4==14
replace right=. if vote4<0| vote4 ==.
by pid: egen mostlyRight = mode(right)

gen extreme =  vote4==10|vote4==11|vote4==12|vote4==13|vote4==14
replace extreme=. if vote4<0| vote4 ==.
by pid: egen mostlyExtreme = mode(extreme)

gen centre = vote4== 3 | vote4== 9
replace centre=. if vote4<0 | vote4 ==.
by pid: egen mostlyCentrist = mode(centre)

* Actual vote

by pid: egen whovotesfor = mode(vote8), minmode
by pid: egen whovotesfor2 = mode(vote8), maxmode


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

gen dropAt14 = sch == 14
gen dropAt15 = sch == 15
gen dropAt16 = sch == 16
gen dropAfter16 = sch >16
* Instruments

gen left15=sch >=15 
gen left16=sch >=16

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
by pidp: egen mumworked14 = max(aux)
drop aux

gen aux = paju == 1
by pidp: egen dadworked14 = max(aux)
drop aux

gen aux = maju == 3
by pidp: egen mumdead14 = max(aux)
drop aux

gen aux = paju == 3
by pidp: egen daddead14 = max(aux)
drop aux

gen aux = maju == 4
by pidp: egen mumaway14 = max(aux)
drop aux

gen aux = paju == 4
by pidp: egen dadaway14 = max(aux)
drop aux

* Compute historical income, conditional on age, sex
gen age2 =age_dv^2
gen age3 = age_dv^3
gen ageSex = age_dv*sex
gen age4= age_dv^4
gen age5= age_dv^5
rlasso fimngrs_dv age_dv age2 age3 age4 age5 ageSex sex if fimngrs_dv !=. 
predict xb
gen resY = fimngrs_dv- xb
by pidp: egen aveinc = mean(resY)

* First job class
gen aux =j1soc90_cc
forvalues k =1/9{
	local lb = `k'*10
	local ub = (`k'+1)*10
	replace aux = `k' if aux >=`lb' & aux<`ub'
}
by pid: egen firstJobClass = max(aux)
replace firstJobClass =. if firstJobClass <0
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


save "/Users/user/Dropbox/Econometrics/educationLiberalism/code/all_waves.dta", replace
use "/Users/user/Dropbox/Econometrics/educationLiberalism/code/all_waves.dta", clear

* ========================= Part 3: Analysis =========================

* collapse data: 1 entry per person.

collapse cg_*_dis firstJobClass aveinc female dropAt14 dropAt16 dropAt15 dropAfter noInterestedPolitics mostly* left15 left16 doby_fix nouk white mumworked14 dadworked14 mumdead14 daddead14 mumaway14 dadaway14 everunemployed , by(pidp)

drop if doby_fix <0

tab firstJob, gen(firstJob_)

* descriptive statistics.

local explanatoryVariables female dropAt14 dropAt16 dropAt15 doby_fix nouk white mumworked14 ///
	dadworked14 mumdead14 daddead14 mumaway14 dadaway14 everunemployed ///
	firstJob_1 firstJob_2 firstJob_3 firstJob_4 firstJob_5 firstJob_6 ///
	firstJob_7 firstJob_8 firstJob_9 

foreach j of local explanatoryVariables {

	replace `j' = .  if `j' <0
}
estpost sum `explanatoryVariables' 
est store des1

estpost sum `explanatoryVariables'  if doby_fix >=1956 & doby_fix<=1959
est store des2

estpost sum `explanatoryVariables'  if doby_fix >=1956 & doby_fix<=1957
est store des3

estpost sum `explanatoryVariables'  if doby_fix >=1958 & doby_fix<=1959
est store des4

esttab des1 des2 des3 des4 using "/Users/user/Dropbox/Econometrics/educationLiberalism/descriptives.tex", ///
	cells(mean(fmt(2)) sd(par fmt(2))  count() ) label booktabs  collabels(none) gaps f noobs ///
	replace se star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) se(%9.3f) 
	

graph bar if dob >=1940 & doby <=1980,  ylabel(, nogrid) over(doby_fix, label(angle(vertical) labsize(small))  ) bar(1, fcolor(black) lcolor(none)) ytitle(Percentage)  graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white) lcolor(white) ) 
 
graph save Graph "/Users/user/Dropbox/Econometrics/educationLiberalism/code/graph_age_prop.gph", replace
	
* The jump.

preserve
collapse left15 left16 dropAt15 dropAt16, by(doby_fix)
 

twoway(lpoly dropAt15 doby_fix if doby_fix<1958 & doby_fix>=1940 & ///
	doby_fix<=1980 , lcolor(black) lpattern(.#.)  lwidth(0.5)) ///
(lpoly dropAt15 doby_fix if doby_fix>=1958 & doby_fix>=1940 & doby_fix<=1980, lcolor(black) lpattern(.#.)  lwidth(0.5)) ///
(scatter dropAt15 doby_fix if doby_fix>=1940 & doby_fix<=1980, ///
	mcolor(gray) msymbol(Oh) msize(small) xtitle(Year of birth) ytitle(Percentage) graphregion(fcolor(white) ///
	lcolor(white)) plotregion(fcolor(white) lcolor(white)) legend(off) xlabel(, angle(45)) ///
	xscale(range(1940 1980))) ///
	(lpoly dropAt16 doby_fix if doby_fix<1958 & doby_fix>=1940 & ///
	doby_fix<=1980 , lcolor(black))(lpoly dropAt16 doby_fix if ///
	doby_fix>=1958 & doby_fix>=1940 & doby_fix<=1980, lcolor(black)) ///
	(scatter dropAt16 doby_fix if doby_fix>=1940 & doby_fix<=1980, xline(1958, lcolor(black)) ///
	mcolor(gray) msymbol(O) msize(small) )
restore


* Testing jumps at the cutoof 
preserve

contract doby
rename _freq freq
gen running_var = doby - 1958
egen sumtot = sum(freq)
gen phat = freq/sumtot

/*rdrandinf  phat running_var , ci(0.05) interfci(0.05) wl(-1) wr(1) statistic(ranksum) 
rdrandinf  phat running_var , ci(0.05) interfci(0.05) wl(-2) wr(2) statistic(ranksum) 
rdrandinf  phat running_var , ci(0.05) interfci(0.05) wl(-3) wr(3) statistic(ranksum) 
rdrandinf  phat running_var , ci(0.05) interfci(0.05) wl(-4) wr(4) statistic(ranksum) 
rdrandinf  phat running_var , ci(0.05) interfci(0.05) wl(-5) wr(5) statistic(ranksum) 
rdrandinf  phat running_var , ci(0.05) interfci(0.05) wl(-6) wr(6) statistic(ranksum) 
*/
rdrobust phat running_var

rdplot phat running_var, graph_options( legend(off))

local resDir "/Users/user/Dropbox/Econometrics/educationLiberalism/"
restore 












local outcomes noInterestedPolitics mostlyLabour mostlyConservative mostlyLeft mostlyRight mostlyExtreme mostlyCentrist 

foreach j of local explanatoryVariables {

	replace `j' = .  if `j' <0
}
estpost sum `outcomes' 
est store des1


esttab des1 using "/Users/user/Dropbox/Econometrics/educationLiberalism/descriptives_outcomes.tex", ///
	cells(mean(fmt(2)) sd(par fmt(2)) min() max()  ) label booktabs  collabels(none) gaps f ///
	replace se star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) se(%9.3f) 

	
	
* Exclusion 

gen earlydropper = (dropAt14 & doby_f <1958) | (dropAt15 & doby_f >=1958)


/*esttab margins_effects_1 margins_effects_2 margins_effects_3   using "/Users/user/Dropbox/Econometrics/educationLiberalism/probit_exclusion.tex", /// 
replace se star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) se(%9.3f) 
*/
gen deadparent = mumdead14==1 | daddead14==1

biprobit mostlyLeft earlydropper female nouk white doby_fix mumworked14 dadworked14 deadparent mumaway14 dadaway14 , vce(robust)
margins, dydx(*) predict(pmarg1) predict(pmarg2) post
eststo mybiprobit

esttab mybiprobit using "/Users/user/Dropbox/Econometrics/educationLiberalism/mybiprobit.tex", replace se ///
title("Biprobit Regression Results") ///
label note("Robust standard errors in parentheses.")

// Create the contingency table and get chi-squared statistic
tabulate mostlyLeft earlydropper, chi2

// Assuming chi2_stat = 10 and n = 100 for demonstration purposes
local  n = e(N)
local chi2_stat =r(chi2) 
// Compute the phi coefficient
scalar phi = sqrt(`chi2_stat' / `n')

// Display the result
display "Phi coefficient: " phi


* 1972 Law

/*

Uncomment to verify bandwidth choice..

local pretreatments  white nouk mumworked14 dadworked14 mumdead14 daddead14 mumaway14 dadaway14 

 rdwinselect doby_fix `pretreatments' , c(1958) wmass level(0.05)  wasym 
*/
local bw58l = 1956 //r(w_left) 
local bw58r = 1959 //r(w_right) 

local resDir "/Users/user/Dropbox/Econometrics/educationLiberalism/"

local outcomes  noInterestedPolitics mostlyLabour mostlyConservative mostlyLeft mostlyRight mostlyExtreme mostlyCentrist 

my_rdrandinf_function `bw58l' `bw58r' 1958 left16 doby_fix "`outcomes'"

matrix mRes58 = r(mRes)
matrix rown mRes58=  noInterestedPolitics mostlyLabour mostlyConservative mostlyLeft mostlyRight mostlyExtreme mostlyCentrist

outtable using "/Users/user/Dropbox/Econometrics/educationLiberalism/mRes1958.tex", mat(mRes58) replace


local resDir "/Users/user/Dropbox/Econometrics/educationLiberalism/"

rdsensitivity  noInterestedPolitics doby_fix, fuzzy(left16) c(1958)  saving(`resDir'sens1958_nopol)  

rdsensitivity  mostlyLeft doby_fix, fuzzy(left16) c(1958) saving(`resDir'sens1958_mostlyLeft)

rdsensitivity  mostlyRight doby_fix, fuzzy(left16) c(1958) saving(`resDir'sens1958_mostlyRight)

preserve
local resDir "/Users/user/Dropbox/Econometrics/educationLiberalism/"
use "/Users/user/Dropbox/Econometrics/educationLiberalism/sens1958_nopol.dta", clear
twoway (contour pvalue t w, levels(20) format(%8.0g) interp(shepard) crule(linear) scolor(white) ecolor(black)), ytitle(Estimate) ylabel(minmax) xtitle(Bandwidth) xtitle(, margin(zero)) xlabel(1(1)10) ztitle(P-value) zlabel(#10) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white))
graph export "`resDir'sensitivity_nonpol.pdf", as(pdf) replace

use "/Users/user/Dropbox/Econometrics/educationLiberalism/sens1958_mostlyLeft.dta", clear
twoway (contour pvalue t w, levels(20) format(%8.0g) interp(shepard) crule(linear) scolor(white) ecolor(black)), ytitle(Estimate) ylabel(minmax) xtitle(Bandwidth) xtitle(, margin(zero)) xlabel(1(1)10) ztitle(P-value) zlabel(#10) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white))
graph export "`resDir'sensitivity_left.pdf", as(pdf) replace

use "/Users/user/Dropbox/Econometrics/educationLiberalism/sens1958_mostlyRight.dta", clear
twoway (contour pvalue t w, levels(20) format(%8.0g) interp(shepard) crule(linear) scolor(white) ecolor(black)), ytitle(Estimate) ylabel(minmax) xtitle(Bandwidth) xtitle(, margin(zero)) xlabel(1(1)10) ztitle(P-value) zlabel(#10) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white))
graph export "`resDir'sensitivity_right.pdf", as(pdf) replace

restore



 
* vulnerable

gen vulnerable =  dadworked14 ==0 |deadparent ==1| mumaway14 ==1| dadaway14==1

local pretreatments  white nouk mumworked14 

rdwinselect doby_fix `pretreatments' if vulnerable==1 , c(1958) wmass level(0.05)  wasym 

preserve 
keep if vulnerable ==1 
local bw58l = 1951 //r(w_left) 
local bw58r = 1964 //r(w_right) 

local resDir "/Users/user/Dropbox/Econometrics/educationLiberalism/vulnerable/"

local outcomes  noInterestedPolitics mostlyLabour mostlyConservative mostlyLeft mostlyRight mostlyExtreme mostlyCentrist 

my_rdrandinf_function `bw58l' `bw58r' 1958 left16 doby_fix "`outcomes'"
matrix mRes58 = r(mRes)

matrix rown mRes58=  noInterestedPolitics mostlyLabour mostlyConservative mostlyLeft mostlyRight mostlyExtreme mostlyCentrist

outtable using "/Users/user/Dropbox/Econometrics/educationLiberalism/vulnerable/mRes1958.tex", mat(mRes58) replace


local resDir "/Users/user/Dropbox/Econometrics/educationLiberalism/vulnerable/"

rdsensitivity  noInterestedPolitics doby_fix, fuzzy(left16) c(1958)  saving(`resDir'sens1958_nopol)  

rdsensitivity  mostlyLeft doby_fix, fuzzy(left16) c(1958) saving(`resDir'sens1958_mostlyLeft)

rdsensitivity  mostlyRight doby_fix, fuzzy(left16) c(1958) saving(`resDir'sens1958_mostlyRight)

restore

preserve
local resDir "/Users/user/Dropbox/Econometrics/educationLiberalism/vulnerable/"
use "/Users/user/Dropbox/Econometrics/educationLiberalism/vulnerable/sens1958_nopol.dta", clear
twoway (contour pvalue t w, levels(20) format(%8.0g) interp(shepard) crule(linear) scolor(white) ecolor(black)), ytitle(Estimate) ylabel(minmax) xtitle(Bandwidth) xtitle(, margin(zero)) xlabel(1(1)10) ztitle(P-value) zlabel(#10) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white))
graph export "`resDir'sensitivity_nonpol.pdf", as(pdf) replace

use "/Users/user/Dropbox/Econometrics/educationLiberalism/vulnerable/sens1958_mostlyLeft.dta", clear
twoway (contour pvalue t w, levels(20) format(%8.0g) interp(shepard) crule(linear) scolor(white) ecolor(black)), ytitle(Estimate) ylabel(minmax) xtitle(Bandwidth) xtitle(, margin(zero)) xlabel(1(1)10) ztitle(P-value) zlabel(#10) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white))
graph export "`resDir'sensitivity_left.pdf", as(pdf) replace

use "/Users/user/Dropbox/Econometrics/educationLiberalism/vulnerable/sens1958_mostlyRight.dta", clear
twoway (contour pvalue t w, levels(20) format(%8.0g) interp(shepard) crule(linear) scolor(white) ecolor(black)), ytitle(Estimate) ylabel(minmax) xtitle(Bandwidth) xtitle(, margin(zero)) xlabel(1(1)10) ztitle(P-value) zlabel(#10) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white))
graph export "`resDir'sensitivity_right.pdf", as(pdf) replace

restore





* non vulnerable


local pretreatments  white nouk mumworked14 

rdwinselect doby_fix `pretreatments' if vulnerable==0 , c(1958) wmass level(0.05)  wasym 

preserve 
keep if vulnerable ==0 
local bw58l = 1956 //r(w_left) 
local bw58r = 1959 //r(w_right) 

local resDir "/Users/user/Dropbox/Econometrics/educationLiberalism/novulnerable/"

local outcomes  noInterestedPolitics mostlyLabour mostlyConservative mostlyLeft mostlyRight mostlyExtreme mostlyCentrist 

my_rdrandinf_function `bw58l' `bw58r' 1958 left16 doby_fix "`outcomes'"
matrix mRes58 = r(mRes)

matrix rown mRes58=  noInterestedPolitics mostlyLabour mostlyConservative mostlyLeft mostlyRight mostlyExtreme mostlyCentrist

outtable using "/Users/user/Dropbox/Econometrics/educationLiberalism/novulnerable/mRes1958.tex", mat(mRes58) replace


local resDir "/Users/user/Dropbox/Econometrics/educationLiberalism/novulnerable/"

rdsensitivity  noInterestedPolitics doby_fix, fuzzy(left16) c(1958)  saving(`resDir'sens1958_nopol)  

rdsensitivity  mostlyLeft doby_fix, fuzzy(left16) c(1958) saving(`resDir'sens1958_mostlyLeft)

rdsensitivity  mostlyRight doby_fix, fuzzy(left16) c(1958) saving(`resDir'sens1958_mostlyRight)

restore

preserve
local resDir "/Users/user/Dropbox/Econometrics/educationLiberalism/novulnerable/"
use "/Users/user/Dropbox/Econometrics/educationLiberalism/novulnerable/sens1958_nopol.dta", clear
twoway (contour pvalue t w, levels(20) format(%8.0g) interp(shepard) crule(linear) scolor(white) ecolor(black)), ytitle(Estimate) ylabel(minmax) xtitle(Bandwidth) xtitle(, margin(zero)) xlabel(1(1)10) ztitle(P-value) zlabel(#10) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white))
graph export "`resDir'sensitivity_nonpol.pdf", as(pdf) replace

use "/Users/user/Dropbox/Econometrics/educationLiberalism/novulnerable/sens1958_mostlyLeft.dta", clear
twoway (contour pvalue t w, levels(20) format(%8.0g) interp(shepard) crule(linear) scolor(white) ecolor(black)), ytitle(Estimate) ylabel(minmax) xtitle(Bandwidth) xtitle(, margin(zero)) xlabel(1(1)10) ztitle(P-value) zlabel(#10) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white))
graph export "`resDir'sensitivity_left.pdf", as(pdf) replace

use "/Users/user/Dropbox/Econometrics/educationLiberalism/novulnerable/sens1958_mostlyRight.dta", clear
twoway (contour pvalue t w, levels(20) format(%8.0g) interp(shepard) crule(linear) scolor(white) ecolor(black)), ytitle(Estimate) ylabel(minmax) xtitle(Bandwidth) xtitle(, margin(zero)) xlabel(1(1)10) ztitle(P-value) zlabel(#10) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white))
graph export "`resDir'sensitivity_right.pdf", as(pdf) replace

restore

