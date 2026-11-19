*! specbunch 0.1.0
program define specbunch, eclass
    version 14.1
    gettoken subcmd 0 : 0, parse(" ,")
    if inlist("`subcmd'","nonbuncher_mean","eti") {
        _specbunch_post `subcmd' `0'
        exit
    }
    if "`subcmd'" == "npmle" {
        _specbunch_npmle `0'
        exit
    }
    if "`subcmd'" == "windows" {
        if "`e(cmd)'" != "specbunch" error 301
        matrix list e(windows)
        exit
    }
    if !inlist("`subcmd'","mass","buncher") {
        display as error "syntax: specbunch mass|buncher varlist, options"
        exit 198
    }
    syntax varlist(numeric min=1 max=2) [if] [in] [pweight fweight], [Kink(real 0) Method(string) Friction(string) Policy(string) Threads(integer 0) Seed(string) Checkallwindows LPiterations(string) LPretryiterations(string) Draws(integer 499) NOPLOT PLOT(string) NAME(string asis) SAVing(string asis)]
    if "`plot'"=="" local plot scaled
    local plot = lower("`plot'")
    if !inlist("`plot'","scaled","raw","none") {
        display as error "plot() must be scaled, raw or none"
        exit 198
    }
    if "`noplot'"!="" local plot none
    // Keep graph options apart from the result-name loop locals below.
    local graph_name `"`name'"'
    local graph_saving `"`saving'"'
    local inputs `varlist'
    local nv : word count `inputs'
    if ("`subcmd'"=="mass" & `nv'!=1) | ("`subcmd'"=="buncher" & `nv'!=2) error 198
    if "`method'" == "" local method cosine
    if "`friction'" == "" local friction gaussian
    if !inlist("`method'","cosine","joint") error 198
    if !inlist("`friction'","gaussian","voigt") error 198
    local entry `friction'
    if "`method'"=="joint" local entry gaussian_jump
    if "`method'"=="joint" & "`friction'"=="voigt" local entry joint_voigt
    if "`subcmd'" == "buncher" local entry bmean
    if "`seed'" == "" {
        local seed 1
        if "`entry'" == "bmean" local seed 1
    }
    if "`lpiterations'" == "" local lpiterations 5000
    if "`lpretryiterations'" == "" local lpretryiterations 31250
    if "`policy'" == "" local policy tiered
    if !inlist("`policy'","tiered","first") error 198
    local pol = cond("`policy'"=="first",2,0)
    local c2 1
    local candidates = cond("`checkallwindows'"!="",2,0)
    local bandrule 2
    local p0rule 1
    local center : display %24.17g `kink'

    {
        local plat linux64
        if strpos(c(machine_type), "Mac") {
            local plat macintel64
            if strpos(c(machine_type), "Apple Silicon") local plat macarm64
        }
        else if c(os) == "Windows" local plat win64
        if c(bit) != 64 error 198
        capture program _specbunch_plugin, plugin using("specbunch_`plat'.plugin")
    }
    if _rc != 0 & _rc != 110 exit _rc
    // Pass if/in directly: do not mark out missing values or zero weights.
    tempvar wt
    if "`weight'" != "" {
        quietly generate double `wt' `exp'
        local wvar `wt'
    }
    tempvar sample_use
    marksample sample_use, novarlist
    markout `sample_use' `inputs' `wvar'
    if "`wvar'" != "" quietly replace `sample_use'=0 if `wvar'<=0
    local reuse = 0 & ("`e(cmd)'"=="specbunch" & "`e(method)'"=="`entry'" & e(core_status)==0)
    capture noisily plugin call _specbunch_plugin `inputs' `wvar' `if' `in', run "`entry'" "`center'" "`threads'" "`seed'" "`pol'" `c2' "`candidates'" "`lpiterations'" "`lpretryiterations'" "`bandrule'" `p0rule' "`draws'"
    if _rc {
        local rc = _rc
        ereturn clear
        exit `rc'
    }
    local reuse = `reuse' & `core_status'==0
    if `reuse' {
        mata: _specbunch_prepare_reuse()
    }
    else ereturn clear
    tempname S V B
    local ns : word count `scalar_names'
    matrix `S' = J(1, `ns', .)
    plugin call _specbunch_plugin, fetch scalars `S'
    local scalar_bits `bits'
    local estimate B
    if "`entry'" == "bmean" local estimate mu
    if "`entry'" == "rd" local estimate tau
    if `core_status' == 0 {
        local idx : list posof "`estimate'" in scalar_names
        matrix `B' = (`S'[1,`idx'])
        matrix colnames `B' = `estimate'
        if "`entry'" == "bmean" | "`entry'" == "rd" {
            local idx : list posof "se" in scalar_names
            matrix `V' = (`S'[1,`idx']^2)
            matrix colnames `V' = `estimate'
            matrix rownames `V' = `estimate'
            if `reuse' ereturn repost b=`B' V=`V'
            else ereturn post `B' `V', esample(`sample_use')
        }
        else {
            if `reuse' ereturn repost b=`B'
            else ereturn post `B', esample(`sample_use')
        }
    }
    local i 0
    foreach key of local scalar_names {
        local ++i
        local name = subinstr("`key'", ".", "_", .)
        ereturn scalar `name' = `S'[1,`i']
    }
    ereturn scalar core_status = `core_status'
    ereturn local scalar_names "`scalar_names'"
    ereturn local scalar_bits "`scalar_bits'"
    ereturn local vector_names "`vector_names'"
    ereturn local string_names "`string_names'"
    // One allocation and one plugin transfer; no ado per-column matrix work.
    local nvectors : word count `vector_names'
    if `nvectors' > 0 {
        tempname vectors
        mata: _specbunch_vector_allocate(st_local("vectors"))
        plugin call _specbunch_plugin, fetch_vectors `vectors'
        mata: _specbunch_vector_post(st_local("vectors"))
    }
    local i 0
    foreach key of local string_names {
        local ++i
        local name = subinstr("`key'", ".", "_", .)
        ereturn local str_`name' `"`str`i''"'
    }
    plugin call _specbunch_plugin, release
    ereturn local cmd "specbunch"
    ereturn local method "`method'"
    ereturn local friction "`friction'"
    if "`entry'"=="bmean" ereturn local method "buncher_mean"
    local running : word 1 of `inputs'
    ereturn local running_variable "`running'"
    ereturn local input_weight_exp "`exp'"
    ereturn scalar kink = `kink' 
    ereturn local policy "`policy'"
    ereturn local flags_text `"`e(str_flags_str)'"'
    capture local tier = string(e(tier))
    ereturn local tier_label "`tier'"
    display as text "specbunch (`method', `friction')" _col(40) "N = " as result %10.0g e(n)
    display as text "{hline 64}"
    if e(core_status) == 0 {
        if inlist("`entry'","gaussian","voigt","gaussian_jump","joint_voigt") {
            display as text "      B             SE(log B)         sigma"
            display as result %12.7g e(B) "    " %12.7g e(se_log_B) "    " %12.7g e(sigma)
            display as text "Policy: `policy'; tier " as result e(tier) as text `"; flags: `e(flags_text)'"'
        }
        else {
            display as text "  Estimate           Model SE"
            display as result %12.7g e(`estimate') "    " %12.7g e(se)
        }
    }
    else display as text "Refused (status " as result e(core_status) as text `"): `e(str_reason)'"'
    display as text "{hline 64}"
    display as text "Standard errors are not confidence intervals."
    if "`subcmd'"=="mass" & "`plot'"!="none" _specbunch_transform_plot, frequency(`plot') name(`graph_name') saving(`graph_saving')
end

// Transform fit of a mass estimate: empirical transform with a +/-2 SE band, the fitted
// curve over the whole grid and the shaded selected window. Joint fits add t*S(t).
program define _specbunch_transform_plot
    version 14.1
    syntax, frequency(string) [name(string asis) saving(string asis)]
    capture local nrows = rowsof(e(transform))
    if _rc exit
    if `nrows'<2 exit
    tempname T
    matrix `T'=e(transform)
    local joint = colnumb(`T',"sine_scaled")<.
    local sfx = cond("`frequency'"=="scaled","_scaled","")
    local wlo = e(window_lo`sfx')
    local whi = e(window_hi`sfx')
    local haswin = `wlo'<. & `whi'<.
    local fmt %9.4g
    local head = strproper("`e(method)'")+" fit, `e(friction)' friction"
    if e(core_status)!=0 local note1 `"`head'; refused: `e(str_reason)'"'
    else if "`e(friction)'"=="voigt" local note1 "`head':  {it:B} = `=strtrim(string(e(B),"`fmt'"))',  {it:a} = `=strtrim(string(e(a),"`fmt'"))',  {it:s}{sup:2} = `=strtrim(string(e(s2),"`fmt'"))'"
    else local note1 "`head':  {it:B} = `=strtrim(string(e(B),"`fmt'"))',  {it:σ} = `=strtrim(string(e(sigma),"`fmt'"))'"
    local note2
    if e(core_status)==0 & `joint' local note2 "jump = `=strtrim(string(e(jump),"`fmt'"))'"
    local xtitle = cond("`frequency'"=="scaled","Scaled frequency {it:t}·{it:σ}{sub:pilot}","Frequency {it:t}")
    local fit `"31 59 115"'
    local edge `"184 134 11"'
    local win `"251 236 200"'
    local wlabel
    if `haswin' local wlabel `"1 "Selected window [`=strtrim(string(`wlo',"`fmt'"))', `=strtrim(string(`whi',"`fmt'"))']""'
    preserve
    quietly {
        drop _all
        svmat double `T', names(col)
        local x = cond("`frequency'"=="scaled","t_scaled","t")
        summarize `x', meanonly
        local xmin=r(min)
        local xmax=r(max)
        // Ticks inside the grid, so the axis spans exactly the display grid.
        _natscale `xmin' `xmax' 7
        local xdelta=r(delta)
        local xt1=ceil(`xmin'/`xdelta')*`xdelta'
        local xt2=floor(`xmax'/`xdelta')*`xdelta'
    }
    local graphs
    local panels = cond(`joint',"cosine sine","cosine")
    foreach p of local panels {
        if "`p'"=="cosine" {
            local y cosine
            local se cosine_se
            local f cosine_fit
            local ytitle "Cosine transform Ĉ(t)"
            local emp "Empirical Ĉ(t)"
            local flab "Fit {it:B}·φ(t)"
            local note `"`note1'"'
        }
        else {
            local y sine_scaled
            local se sine_scaled_se
            local f sine_fit
            local ytitle "Scaled sine transform {it:t}·Ŝ(t)"
            local emp "Empirical {it:t}·Ŝ(t)"
            local flab "Fit jump·φ(t)"
            local note `"`note2'"'
        }
        tempvar lo hi wl wh
        quietly {
            generate double `lo'=`y'-2*`se'
            generate double `hi'=`y'+2*`se'
            summarize `lo', meanonly
            local ymin=min(r(min),0)
            summarize `hi', meanonly
            local ymax=max(r(max),0)
            summarize `f', meanonly
            local nfit=r(N)
            if `nfit'>0 {
                local ymin=min(`ymin',r(min))
                local ymax=max(`ymax',r(max))
            }
            local ymax=`ymax'+.18*(`ymax'-`ymin')
            // Nice tick range, so the window shading spans the whole plot region.
            _natscale `ymin' `ymax' 8
            local ydelta=r(delta)
            local ymin=floor(`ymin'/`ydelta')*`ydelta'
            local ymax=ceil(`ymax'/`ydelta')*`ydelta'
            generate double `wl'=`ymin' if inrange(`x',`wlo',`whi')
            generate double `wh'=`ymax' if inrange(`x',`wlo',`whi')
        }
        local fitlayer
        if `nfit'>0 local fitlayer (line `f' `x', lcolor(`"`fit'"') lwidth(medthin))
        local winlayer
        local winline
        local k 0
        if `haswin' {
            local k 1
            local winlayer (rarea `wh' `wl' `x', color(`"`win'"') lwidth(none))
            local winline xline(`wlo' `whi', lpattern(dash) lcolor(`"`edge'"') lwidth(vthin))
        }
        local order `=`k'+2' "`emp'" `=`k'+1' "± 2 SE"
        if `nfit'>0 local order `order' `=`k'+3' "`flab'"
        local g
        if `joint' {
            tempname g
            local nodraw nodraw
        }
        twoway `winlayer' (rarea `hi' `lo' `x', color(gs13) lwidth(none)) ///
               (line `y' `x', lcolor(gs6) lwidth(vthin)) `fitlayer', ///
               graphregion(color(white) margin(small)) plotregion(color(white) lstyle(none) margin(zero)) ///
               ylabel(`ymin'(`ydelta')`ymax', angle(0) nogrid labsize(small)) xlabel(`xt1'(`xdelta')`xt2', nogrid labsize(small)) ///
               yline(0, lcolor(gs10) lwidth(vthin) lpattern(solid)) `winline' ///
               xscale(range(`xmin' `xmax')) yscale(range(`ymin' `ymax')) ///
               xtitle(`"`xtitle'"', size(small)) ytitle(`"`ytitle'"', size(small)) ///
               title(`"`note'"', size(small) position(11) justification(left) color(black)) ///
               legend(order(`order' `wlabel') ring(0) position(1) cols(1) size(vsmall) symxsize(*.6) region(lstyle(none) fcolor(white))) ///
               `nodraw' `=cond(`joint',"name(`g')",cond(`"`name'"'!="","name(`name')",""))' ///
               `=cond(!`joint' & `"`saving'"'!="","saving(`saving')","")'
        local graphs `graphs' `g'
    }
    restore
    if `joint' {
        graph combine `graphs', cols(2) xsize(10) ysize(3.9) iscale(*1.25) graphregion(color(white) margin(zero)) ///
            `=cond(`"`name'"'!="","name(`name')","")' `=cond(`"`saving'"'!="","saving(`saving')","")'
        graph drop `graphs'
    }
end

// These helpers only transport C results. st_matrix supplies the legacy
// r1...rN/c1 labels; compatible existing matrices avoid repeated allocation.
program define _specbunch_post, eclass
    version 14.1
    gettoken operation 0 : 0
    syntax varlist(numeric min=1 max=2) [if] [in] [pweight fweight], [MASSfit(name) BUNCHERfit(name) MUb(string) SEMUb(real 0) BANDwidth(real .05) DEGree(integer 1) BINWidth(real .0025) TAXBelow(string) TAXAbove(string) NETBelow(string) NETAbove(string) INTegrationintervals(integer 200) Threads(integer 0)]
    local inputs `varlist'
    local nv : word count `inputs'
    if ("`operation'"=="nonbuncher_mean" & `nv'!=2) | ("`operation'"=="eti" & `nv'!=1) error 198
    if "`massfit'"!="" quietly estimates restore `massfit'
    if "`e(cmd)'"!="specbunch" | !inlist("`e(method)'","cosine","joint") error 301
    if e(core_status)!=0 error 498
    foreach key in B sigma kink window_lo window_hi {
        local `key' : display %24.17g e(`key')
    }
    local a 0
    if "`e(friction)'"=="voigt" local a : display %24.17g e(a)
    if "`operation'"=="nonbuncher_mean" {
        if "`buncherfit'"!="" {
            if "`mub'"!="" error 198
            quietly estimates restore `buncherfit'
            if "`e(cmd)'"!="specbunch" | "`e(method)'"!="buncher_mean" error 301
            if e(core_status)!=0 | e(kink)!=real("`kink'") error 498
            local mub : display %24.17g e(mu)
            local semub : display %24.17g e(se)
        }
        if "`mub'"=="" {
            display as error "buncherfit() or mub() is required"
            exit 198
        }
        local tail `""`mub'" "`semub'" "`bandwidth'" "`degree'" "`binwidth'""'
    }
    else {
        if "`netbelow'`netabove'"=="" {
            if "`taxbelow'"=="" | "`taxabove'"=="" error 198
            local netbelow : display %24.17g (1-real("`taxbelow'"))
            local netabove : display %24.17g (1-real("`taxabove'"))
        }
        else if "`taxbelow'`taxabove'"!="" | "`netbelow'"=="" | "`netabove'"=="" error 198
        local tail `""`window_lo'" "`window_hi'" "`netbelow'" "`netabove'" "`integrationintervals'""'
    }
    local plat linux64
    if strpos(c(machine_type),"Mac") {
        local plat macintel64
        if strpos(c(machine_type),"Apple Silicon") local plat macarm64
    }
    else if c(os)=="Windows" local plat win64
    capture program _specbunch_plugin, plugin using("specbunch_`plat'.plugin")
    if _rc!=0 & _rc!=110 exit _rc
    tempvar wt
    if "`weight'"!="" {
        quietly generate double `wt' `exp'
        local wvar `wt'
    }
    capture noisily plugin call _specbunch_plugin `inputs' `wvar' `if' `in', `operation' "`B'" "`sigma'" "`a'" "`kink'" "`threads'" `tail'
    if _rc {
        local rc=_rc
        ereturn clear
        exit `rc'
    }
    ereturn clear
    tempname S
    local ns : word count `scalar_names'
    matrix `S'=J(1,`ns',.)
    plugin call _specbunch_plugin, fetch scalars `S'
    local scalar_bits `bits'
    if `core_status'==0 {
        tempname coefficient
        local target m_NB
        if "`operation'"=="eti" local target saez.e
        local idx : list posof "`target'" in scalar_names
        matrix `coefficient'=(`S'[1,`idx'])
        local label=subinstr("`target'",".","_",.)
        matrix colnames `coefficient'=`label'
        ereturn post `coefficient'
    }
    local i 0
    foreach key of local scalar_names {
        local ++i
        local name = subinstr("`key'", ".", "_", .)
        ereturn scalar `name' = `S'[1,`i']
    }
    ereturn scalar core_status = `core_status'
    ereturn local scalar_names "`scalar_names'"
    ereturn local scalar_bits "`scalar_bits'"
    ereturn local vector_names "`vector_names'"
    ereturn local string_names "`string_names'"
    // One allocation and one plugin transfer; no ado per-column matrix work.
    local nvectors : word count `vector_names'
    if `nvectors' > 0 {
        tempname vectors
        mata: _specbunch_vector_allocate(st_local("vectors"))
        plugin call _specbunch_plugin, fetch_vectors `vectors'
        mata: _specbunch_vector_post(st_local("vectors"))
    }
    local i 0
    foreach key of local string_names {
        local ++i
        local name = subinstr("`key'", ".", "_", .)
        ereturn local str_`name' `"`str`i''"'
    }
    plugin call _specbunch_plugin, release
    ereturn local cmd "specbunch"
    ereturn local method "`operation'"
    ereturn local flags_text `"`e(str_flags)'"'
    if "`operation'"=="nonbuncher_mean" {
        display as text "Nonbuncher mean: " as result e(m_NB) as text "; SE: " as result e(se)
        display as text "Conditional standard error, not a confidence interval. Report the Gaussian/Voigt pair."
    }
    else display as text "ETI: Saez = " as result e(saez_e) as text "; continuous completion = " as result e(continuous_e)
    display as text `"Flags: `e(flags_text)'"'
end

program define _specbunch_npmle, eclass
    version 14.1
    if "`e(cmd)'"!="specbunch" | !inlist("`e(method)'","cosine","joint") error 301
    if e(core_status)!=0 {
        display as error "npmle requires a coherent mass result"
        exit 498
    }
    syntax [, Bins(name) Threads(integer 0) COREintervals(integer 12) GRIDpoints(integer 1001) EXACTthreshold(integer 1000) MAXiterations(integer 200) LAMBDAmin(real .01) LAMBDAmax(real 10000) LAMBDApoints(integer 15) REFINEpoints(integer 9) HISTogramwidth(real 0) QUADraturetolerance(real 1e-9) OPTimalitytolerance(real 1e-7) NOPLOT PLOT(string) NAME(string asis) SAVing(string asis)]
    if "`plot'"=="" local plot both
    local plot = lower("`plot'")
    if !inlist("`plot'","both","pdf","cdf","none") {
        display as error "plot() must be both, pdf, cdf or none"
        exit 198
    }
    if "`noplot'"!="" local plot none
    local z `e(running_variable)'
    local wexp `"`e(input_weight_exp)'"' 
    tempvar weight sample_use
    quietly generate byte `sample_use'=e(sample)
    if "`wexp'"!="" {
        quietly generate double `weight' `wexp'
        local wvar `weight'
    }
    foreach key in B sigma kink {
        local `key' : display %24.17g e(`key')
    }
    local a 0
    if "`e(friction)'"=="voigt" local a : display %24.17g e(a)
    if "`bins'"=="" local bins _none
    plugin call _specbunch_plugin `z' `wvar' if `sample_use', npmle "`B'" "`sigma'" "`a'" "`kink'" "`threads'" "`coreintervals'" "`gridpoints'" "`exactthreshold'" "`maxiterations'" "`lambdamin'" "`lambdamax'" "`lambdapoints'" "`refinepoints'" "`histogramwidth'" "`quadraturetolerance'" "`optimalitytolerance'" "`bins'"
    tempname scalars vectors
    local ns : word count `scalar_names'
    matrix `scalars'=J(1,`ns',.)
    plugin call _specbunch_plugin, fetch scalars `scalars'
    ereturn local npmle_scalar_bits "`bits'"
    local i 0
    foreach key of local scalar_names {
        local ++i
        ereturn scalar npmle_`key'=`scalars'[1,`i']
    }
    ereturn local npmle_scalar_names "`scalar_names'"
    ereturn local npmle_vector_names "`vector_names'"
    ereturn local npmle_string_names "`string_names'"
    mata: _specbunch_vector_allocate(st_local("vectors"))
    plugin call _specbunch_plugin, fetch_vectors `vectors'
    mata: _specbunch_npmle_post(st_local("vectors"))
    local i 0
    foreach key of local string_names {
        local ++i
        ereturn local npmle_`key' `"`str`i''"'
    }
    plugin call _specbunch_plugin, release
    display as text "NPMLE: `e(npmle_status)'; lambda = " as result e(npmle_lambda) as text "; EDF = " as result e(npmle_edf)
    display as text "The latent density jump is weakly identified."
    if "`plot'"!="none" _specbunch_npmle_plot, kind(`plot') name(`name') saving(`saving')
end

// Density view: data histogram, fitted observed density (bunchers included) and fitted
// latent density without the atom. CDF view: empirical, observed and latent CDF (jump B).
program define _specbunch_npmle_plot
    version 14.1
    syntax, kind(string) [name(string asis) saving(string asis)]
    tempname C H
    matrix `C'=e(npmle_curves)
    matrix `H'=e(npmle_histogram)
    local kink=e(npmle_kink)
    local xtitle "`e(running_variable)'"
    if "`xtitle'"=="" local xtitle z
    local note "Estimated {it:B} = `=strtrim(string(e(npmle_B),"%9.4g"))',  {it:σ} = `=strtrim(string(e(npmle_sigma),"%9.4g"))'"
    if e(npmle_a)>0 & e(npmle_a)<. local note "`note',  {it:a} = `=strtrim(string(e(npmle_a),"%9.4g"))'"
    local status "`e(npmle_status)'"
    local obs `"31 59 115"'
    local lat `"217 130 43"'
    local common graphregion(color(white) margin(small)) plotregion(color(white) lstyle(none)) ylabel(, angle(0) nogrid labsize(small)) xlabel(, nogrid labsize(small)) xtitle(`"`xtitle'"', size(small)) xline(`kink', lpattern(dot) lcolor(gs6) lwidth(vthin))
    preserve
    quietly {
        drop _all
        local nc=rowsof(`C')
        local nh=rowsof(`H')
        set obs `=max(`nc',4*`nh')'
        svmat double `C', names(col)
        // Bars as one stepped area: (lo,0) (lo,d) (hi,d) (hi,0) for each bin.
        generate long _bin=ceil(_n/4) if _n<=4*`nh'
        local jlo=colnumb(`H',"bin_lo")
        local jhi=colnumb(`H',"bin_hi")
        local jsh=colnumb(`H',"observed_share")
        generate double _hx=cond(inlist(mod(_n-1,4),0,1),`H'[_bin,`jlo'],`H'[_bin,`jhi']) if _bin<.
        generate double _hy=cond(inlist(mod(_n-1,4),1,2),`H'[_bin,`jsh']/(`H'[_bin,`jhi']-`H'[_bin,`jlo']),0) if _bin<.
        // View: .1% to 99.9% of the empirical CDF (display only; the fit uses every observation).
        summarize grid if empirical_cdf>=.001, meanonly
        local lo=min(r(min),`kink')
        summarize grid if empirical_cdf>=.999, meanonly
        local hi=max(r(min),`kink')
        summarize grid, meanonly
        local pad=.04*cond(`hi'>`lo',`hi'-`lo',1)
        local lo=max(`lo'-`pad',r(min))
        local hi=min(`hi'+`pad',r(max))
        replace _hx=min(max(_hx,`lo'),`hi') if _hx<.
        summarize _hy, meanonly
        local top=r(max)
        summarize observed_pdf if inrange(grid,`lo',`hi'), meanonly
        local top=1.12*max(`top',r(max))
    }
    local inx inrange(grid,`lo',`hi')
    local graphs
    if inlist("`kind'","both","pdf") {
        local g pdf
        if "`kind'"=="both" {
            tempname g
            local nodraw nodraw
        }
        else local g `name'
        local stat
        if "`status'"!="converged" local stat `""status: `status'""'
        twoway (area _hy _hx, color(gs12) lwidth(none)) ///
               (line observed_pdf grid if `inx', lcolor(`"`obs'"') lwidth(medthin)) ///
               (line latent_pdf grid if `inx', lcolor(`"`lat'"') lwidth(medthin) lpattern(dash)), ///
               `common' xscale(range(`lo' `hi')) yscale(range(0 `top')) ytitle("Density", size(small)) ///
               text(`top' `lo' `"`note'"' "Dotted line: kink" `stat', placement(se) size(small) justification(left)) ///
               legend(order(1 "Data (histogram)" 2 "Observed fit" 3 "Latent fit, excl. bunchers") ///
               ring(0) position(1) cols(1) size(vsmall) symxsize(*.6) region(lstyle(none) fcolor(white))) ///
               `nodraw' `=cond("`kind'"=="both","name(`g')",cond(`"`name'"'!="","name(`name')",""))' ///
               `=cond("`kind'"!="both" & `"`saving'"'!="","saving(`saving')","")'
        local graphs `graphs' `g'
    }
    if inlist("`kind'","both","cdf") {
        local g cdf
        if "`kind'"=="both" {
            tempname g
            local nodraw nodraw
        }
        twoway (line empirical_cdf grid if `inx', connect(J) lcolor(gs12) lwidth(thick)) ///
               (line observed_cdf grid if `inx', lcolor(`"`obs'"') lwidth(medthin)) ///
               (line latent_cdf grid if `inx', lcolor(`"`lat'"') lwidth(medthin) lpattern(dash)), ///
               `common' xscale(range(`lo' `hi')) yscale(range(0 1)) ytitle("Cumulative probability", size(small)) ///
               legend(order(1 "Empirical CDF" 2 "Observed fit" 3 "Latent fit (jump {it:B} at kink)") ///
               ring(0) position(11) cols(1) size(vsmall) symxsize(*.6) region(lstyle(none) fcolor(white))) ///
               `nodraw' `=cond("`kind'"=="both","name(`g')",cond(`"`name'"'!="","name(`name')",""))' ///
               `=cond("`kind'"!="both" & `"`saving'"'!="","saving(`saving')","")'
        local graphs `graphs' `g'
    }
    restore
    if "`kind'"=="both" {
        graph combine `graphs', cols(2) xsize(10) ysize(3.9) iscale(*1.25) graphregion(color(white) margin(zero)) ///
            `=cond(`"`name'"'!="","name(`name')","")' `=cond(`"`saving'"'!="","saving(`saving')","")'
        graph drop `graphs'
    }
end

mata:
void _specbunch_vector_allocate(string scalar target)
{
    real scalar j, n
    string rowvector keys
    keys=tokens(st_local("vector_names")); n=0
    for(j=1;j<=cols(keys);j++) n=max((n,strtoreal(st_local("len"+strofreal(j)))))
    st_matrix(target,J(max((1,n)),cols(keys),.))
}
void _specbunch_vector_post(string scalar source)
{
    real matrix values
    real scalar j, n, nc
    real rowvector candidates
    string rowvector keys, candidate_names
    string scalar key, target

    values=st_matrix(source); keys=tokens(st_local("vector_names"))
    candidates=J(1,cols(keys),0); candidate_names=J(1,cols(keys),""); nc=0
    for(j=1;j<=cols(keys);j++) {
        n=strtoreal(st_local("len"+strofreal(j)))
        if(n==0) continue
        key="vector"+strtrim(strofreal(j)); target="e("+key+")"
        _specbunch_put_matrix(target,values[|1,j\n,j|])
        st_global("e(bits_"+key+")",st_local("vbits"+strofreal(j)))
        if(substr(keys[j],1,8)=="windows.") {
            nc++; candidates[nc]=j; candidate_names[nc]=substr(keys[j],9,.)
        }
    }
    if(nc>0) {
        candidates=candidates[|1\nc|]; candidate_names=candidate_names[|1\nc|]
        n=strtoreal(st_local("len"+strofreal(candidates[1])))
        _specbunch_put_matrix("e(windows)",values[|1,1\n,cols(values)|][,candidates])
        st_matrixcolstripe("e(windows)",(J(nc,1,""),candidate_names'))
    }
    nc=0; candidates=J(1,cols(keys),0); candidate_names=J(1,cols(keys),"")
    for(j=1;j<=cols(keys);j++) {
        if(substr(keys[j],1,10)=="transform." & strtoreal(st_local("len"+strofreal(j)))>0) {
            nc++; candidates[nc]=j; candidate_names[nc]=substr(keys[j],11,.)
        }
    }
    if(nc>0) {
        candidates=candidates[|1\nc|]; candidate_names=candidate_names[|1\nc|]
        n=strtoreal(st_local("len"+strofreal(candidates[1])))
        _specbunch_put_matrix("e(transform)",values[|1,1\n,cols(values)|][,candidates])
        st_matrixcolstripe("e(transform)",(J(nc,1,""),candidate_names'))
    }
    nc=0; candidates=J(1,cols(keys),0); candidate_names=J(1,cols(keys),"")
    for(j=1;j<=cols(keys);j++) {
        if(substr(keys[j],1,10)=="local_fit." & strtoreal(st_local("len"+strofreal(j)))>0) {
            nc++; candidates[nc]=j; candidate_names[nc]=substr(keys[j],11,.)
        }
    }
    if(nc>0) {
        candidates=candidates[|1\nc|]; candidate_names=candidate_names[|1\nc|]
        n=strtoreal(st_local("len"+strofreal(candidates[1])))
        _specbunch_put_matrix("e(local_fit)",values[|1,1\n,cols(values)|][,candidates])
        st_matrixcolstripe("e(local_fit)",(J(nc,1,""),candidate_names'))
    }

}
void _specbunch_put_matrix(string scalar target, real matrix value)
{
    real matrix old
    string matrix stripe
    old=st_matrix(target)
    if(rows(old)==rows(value) & cols(old)==cols(value)) {
        st_replacematrix(target,value)
        stripe=(J(rows(value),1,""),"r":+strtrim(strofreal(1::rows(value),"%18.0g")))
        if(any(st_matrixrowstripe(target):!=stripe)) st_matrixrowstripe(target,stripe)
        if(cols(value)==1) st_matrixcolstripe(target,("","c1"))
    }
    else st_matrix(target,value)
}
void _specbunch_prepare_reuse()
{
    string colvector names
    string rowvector keep, keys
    real scalar i, nc
    keys=tokens(st_local("vector_names"))
    keep=("b","V","candidates",J(1,cols(keys),""))
    for(i=1;i<=cols(keys);i++) keep[i+3]=subinstr(keys[i],".","_")
    names=st_dir("e()","matrix","*")
    for(i=1;i<=rows(names);i++)
        if(!anyof(keep,names[i])) st_matrix("e("+names[i]+")",J(0,0,.))
    // Empty vectors must not leave a matrix from a previous successful call.
    for(i=1;i<=cols(keys);i++)
        if(strtoreal(st_local("len"+strofreal(i)))==0)
            st_matrix("e("+subinstr(keys[i],".","_")+")",J(0,0,.))
    nc=0
    for(i=1;i<=cols(keys);i++)
        if(substr(keys[i],1,5)=="cand." & strtoreal(st_local("len"+strofreal(i)))>0) nc++
    if(nc==0) st_matrix("e(windows)",J(0,0,.))
    names=st_dir("e()","numscalar","*")
    for(i=1;i<=rows(names);i++) st_numscalar("e("+names[i]+")",J(0,0,.))
    names=st_dir("e()","macro","*")
    for(i=1;i<=rows(names);i++)
        if(names[i]!="properties" & names[i]!="cmd") st_global("e("+names[i]+")","")
}

void _specbunch_npmle_post(string scalar source)
{
    real matrix values
    string rowvector keys, selected_names
    string scalar key, prefix
    real rowvector selected
    real scalar j,n,k,count
    values=st_matrix(source);keys=tokens(st_local("vector_names"))
    for(j=1;j<=cols(keys);j++) {
        n=strtoreal(st_local("len"+strofreal(j)))
        key="npmle_vector"+strtrim(strofreal(j))
        if(n>0) _specbunch_put_matrix("e("+key+")",values[|1,j\n,j|])
        st_global("e(npmle_bits"+strtrim(strofreal(j))+")",st_local("vbits"+strofreal(j)))
    }
    for(k=1;k<=3;k++) {
        prefix=("curves","histogram","lambda_path")[k]
        selected=J(1,0,.);selected_names=J(1,0,"")
        for(j=1;j<=cols(keys);j++) if(substr(keys[j],1,strlen(prefix)+1)==prefix+".") {
            selected=selected,j;selected_names=selected_names,substr(keys[j],strlen(prefix)+2,.)
        }
        n=strtoreal(st_local("len"+strofreal(selected[1])))
        key="e(npmle_"+prefix+")"
        _specbunch_put_matrix(key,values[|1,1\n,cols(values)|][,selected])
        st_matrixcolstripe(key,(J(cols(selected),1,""),selected_names'))
    }
}

end
