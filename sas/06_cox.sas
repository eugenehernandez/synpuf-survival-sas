/*==============================================================================
  Program   : 06_cox.sas
  Project   : synpuf-survival-sas
  Purpose   : Cox proportional hazards models for CHF and all-cause death.
              Unadjusted and fully adjusted models, the effect of the ties
              method, adjusted survival curves, CHF by age and CHF by sex
              interactions, and a CHF by log(time) term as a first look at a
              time-varying effect.
  Inputs    : SYNWORK.ANALYSIS (row level, one record per person)
  Outputs   : SYNWORK.COX_RESULTS,  results/tables/cox_models.csv
              SYNWORK.COX_INT_HR,   results/tables/cox_interaction_hr.csv
              results/figures/adjusted_survival_chf.png
  Seminar   : section 5 (5.1 Cox model, 5.2 survival curves after Cox,
              5.3 interactions, 5.5 time-varying covariates)
  Notes     : Every model uses TIES=EFRON, the method lifelines and the
              Python cross-check use. Death dates resolve only to the month,
              so roughly 85 deaths share each event time and the ties method
              is not a technicality. The Breslow fit of the unadjusted model is
              kept for comparison, and the plasmode analysis returns to it.
              GCONV=1E-12 tightens convergence so that estimates agree with
              the Python reference to six decimal places.
              In this synthetic release death is unrelated to CHF and to the
              other covariates, so hazard ratios near 1 are the expected
              result. The program demonstrates the modeling workflow and makes
              no clinical claim.
==============================================================================*/

/*------------------------------------------------------------------------------
  Load 00_config.sas from this program's folder every time, so the program
  can run alone or from run_all.sas and always sees the current settings.
------------------------------------------------------------------------------*/
%macro _boot;
  %local p;
  %if %symexist(_SASPROGRAMFILE) %then %let p = %qsysfunc(dequote(&_SASPROGRAMFILE));
  %if %length(&p) = 0 and %sysfunc(envlen(SAS_EXECFILEPATH)) > 0 %then
    %let p = %qsysfunc(sysget(SAS_EXECFILEPATH));
  %if %length(&p) = 0 %then %do;
    %put ERROR: [boot] Program location unknown. Open this file from disk in SAS.;
    %return;
  %end;
  %include "%unquote(%qsubstr(&p, 1, %sysfunc(findc(&p, %str(\/), b))))00_config.sas";
%mend _boot;
%_boot

/* Covariate list for the adjusted models. Race enters as a CLASS variable
   with White as the reference level. */
%let adjvars = chf age female race esrd alzh ckd cancer copd depress diab ihd
               osteo raoa stroke logreimb;

/* Age centered at 75, so that in the interaction model the CHF coefficient is
   the log hazard ratio for a 75-year-old man. */
data synwork._cox / view=synwork._cox;
  set synwork.analysis;
  agec = age - 75;
  label agec = 'Age minus 75 (years)';
run;

/*------------------------------------------------------------------------------
  Helper: pull one parameter from an ODS ParameterEstimates table and compare
  it with the Python reference. HR=1 also checks the hazard ratio and its 95
  percent Wald limits (from RISKLIMITS). CLS selects a CLASS level.
------------------------------------------------------------------------------*/
%macro pe_check(ds=, tag=, param=, name=, cls=, hr=1, tol=0.00001);
  %local b s h l u;
  proc sql noprint;
    select Estimate format=best16., StdErr format=best16.
           %if &hr %then %do; , HazardRatio format=best16., HRLowerCL format=best16.,
                                HRUpperCL format=best16. %end;
      into :b trimmed, :s trimmed
           %if &hr %then %do; , :h trimmed, :l trimmed, :u trimmed %end;
    from &ds
    where upcase(Parameter) = upcase("&param")
          %if %length(&cls) %then %do; and ClassVal0 = "&cls" %end;
    ;
  quit;
  %check(metric=&tag._&name._beta, actual=&b, tol=&tol)
  %check(metric=&tag._&name._se,   actual=&s, tol=&tol)
  %if &hr %then %do;
    %check(metric=&tag._&name._hr,    actual=&h, tol=&tol)
    %check(metric=&tag._&name._hrlcl, actual=&l, tol=&tol)
    %check(metric=&tag._&name._hrucl, actual=&u, tol=&tol)
  %end;
%mend pe_check;

/* -2 log partial likelihood with covariates. */
%macro m2ll_check(ds=, tag=);
  %local v;
  proc sql noprint;
    select WithCovariates format=best16. into :v trimmed
    from &ds where compress(upcase(Criterion)) = '-2LOGL';
  quit;
  %check(metric=&tag._m2ll, actual=&v, tol=0.001)
%mend m2ll_check;

/*------------------------------------------------------------------------------
  5.1a Unadjusted model, Efron ties.
------------------------------------------------------------------------------*/
title 'Unadjusted Cox model for CHF, Efron ties';
ods output ParameterEstimates = synwork._pe_unadj
           FitStatistics      = synwork._fit_unadj;
proc phreg data=synwork.analysis;
  model time * event(0) = chf / ties=efron risklimits gconv=1e-12;
  hazardratio 'CHF vs no CHF' chf / cl=wald;
run;
title;
%pe_check(ds=synwork._pe_unadj, tag=cox_unadj, param=chf, name=chf)
%m2ll_check(ds=synwork._fit_unadj, tag=cox_unadj)

/*------------------------------------------------------------------------------
  5.1b Unadjusted model, Breslow ties. With about 85 deaths per event time
  the two approximations differ in the third decimal of the log hazard
  ratio here. The plasmode analysis shows a case where the gap is larger.
------------------------------------------------------------------------------*/
title 'Unadjusted Cox model for CHF, Breslow ties';
ods output ParameterEstimates = synwork._pe_unadj_b
           FitStatistics      = synwork._fit_unadj_b;
proc phreg data=synwork.analysis;
  model time * event(0) = chf / ties=breslow risklimits gconv=1e-12;
run;
title;
%pe_check(ds=synwork._pe_unadj_b, tag=cox_unadj_breslow, param=chf, name=chf, hr=0)
%m2ll_check(ds=synwork._fit_unadj_b, tag=cox_unadj_breslow)

/*------------------------------------------------------------------------------
  5.1c Fully adjusted model.
------------------------------------------------------------------------------*/
title 'Adjusted Cox model, Efron ties';
ods output ParameterEstimates = synwork._pe_adj
           FitStatistics      = synwork._fit_adj
           HazardRatios       = synwork._hr_adj;
proc phreg data=synwork.analysis;
  class race(ref='White') / param=ref;
  model time * event(0) = &adjvars / ties=efron risklimits gconv=1e-12;
  hazardratio 'CHF vs no CHF'      chf    / cl=wald;
  hazardratio 'Age, per 10 years'  age    / units=10 cl=wald;
  hazardratio 'Female vs male'     female / cl=wald;
run;
title;

%pe_check(ds=synwork._pe_adj, tag=cox_adj, param=chf, name=chf)
%macro adj_checks;
  %local i v;
  %let i = 1;
  %let v = %scan(age female esrd alzh ckd cancer copd depress diab ihd osteo raoa stroke logreimb, &i);
  %do %while (%length(&v) > 0);
    %pe_check(ds=synwork._pe_adj, tag=cox_adj, param=&v, name=&v, hr=0)
    %let i = %eval(&i + 1);
    %let v = %scan(age female esrd alzh ckd cancer copd depress diab ihd osteo raoa stroke logreimb, &i);
  %end;
%mend adj_checks;
%adj_checks
%pe_check(ds=synwork._pe_adj, tag=cox_adj, param=race, cls=Black,    name=race_black, hr=0)
%pe_check(ds=synwork._pe_adj, tag=cox_adj, param=race, cls=Other,    name=race_other, hr=0)
%pe_check(ds=synwork._pe_adj, tag=cox_adj, param=race, cls=Hispanic, name=race_hisp,  hr=0)
%m2ll_check(ds=synwork._fit_adj, tag=cox_adj)

/*------------------------------------------------------------------------------
  5.3 Interactions. CHF by age (centered at 75) and CHF by sex, adjusted for
  the same covariates. The HAZARDRATIO statement reports the CHF hazard ratio
  at ages 70, 80, and 90 for men and for women.
------------------------------------------------------------------------------*/
title 'Adjusted Cox model with CHF by age and CHF by sex interactions';
ods output ParameterEstimates = synwork._pe_int
           FitStatistics      = synwork._fit_int
           HazardRatios       = synwork._hr_int;
proc phreg data=synwork._cox;
  class race(ref='White') / param=ref;
  model time * event(0) = chf agec female race esrd alzh ckd cancer copd depress
                          diab ihd osteo raoa stroke logreimb chf*agec chf*female
                          / ties=efron gconv=1e-12;
  hazardratio 'CHF vs no CHF at ages 70, 80, 90 by sex' chf
              / at(agec=-5 5 15 female=0 1) cl=wald;
run;
title;

%pe_check(ds=synwork._pe_int, tag=cox_int, param=chf,        name=chf,        hr=0)
%pe_check(ds=synwork._pe_int, tag=cox_int, param=chf*agec,   name=chf_agec,   hr=0)
%pe_check(ds=synwork._pe_int, tag=cox_int, param=chf*female, name=chf_female, hr=0)
%m2ll_check(ds=synwork._fit_int, tag=cox_int)

/* The same six hazard ratios computed from the coefficients, which gives a
   check that does not depend on the layout of the HazardRatios table. */
proc sql noprint;
  select Estimate format=best16. into :b_chf trimmed from synwork._pe_int where upcase(Parameter) = 'CHF';
  select Estimate format=best16. into :b_ca  trimmed from synwork._pe_int where upcase(Parameter) = 'CHF*AGEC';
  select Estimate format=best16. into :b_cf  trimmed from synwork._pe_int where upcase(Parameter) = 'CHF*FEMALE';
quit;

data synwork.cox_int_hr;
  length sex $6;
  do age = 70, 80, 90;
    do female = 0, 1;
      sex = ifc(female = 1, 'Female', 'Male');
      hr_chf = exp(&b_chf + &b_ca * (age - 75) + &b_cf * female);
      output;
    end;
  end;
  label hr_chf = 'Hazard ratio for CHF vs no CHF';
run;

%macro int_hr_checks;
  %local a f v;
  %do a = 70 %to 90 %by 10;
    %do f = 0 %to 1;
      proc sql noprint;
        select hr_chf format=best16. into :v trimmed from synwork.cox_int_hr
        where age = &a and female = &f;
      quit;
      %check(metric=cox_int_hr_chf_age&a._female&f, actual=&v, tol=0.00001)
    %end;
  %end;
%mend int_hr_checks;
%int_hr_checks

title 'CHF hazard ratio by age and sex (from the HAZARDRATIO statement)';
proc print data=synwork._hr_int noobs;
run;
title;
%export_csv(data=synwork._hr_int, name=cox_interaction_hr)

/*------------------------------------------------------------------------------
  5.5 A time-varying CHF effect. A programming statement creates
  chf_lt = chf * log(time), which PROC PHREG evaluates at every event time.
  A coefficient on chf_lt different from 0 means the CHF hazard ratio changes
  with time (a proportional hazards violation). 08_ph_assumption.sas tests
  this formally. The plasmode data build in a known effect of this form.
------------------------------------------------------------------------------*/
title 'Cox model with a CHF by log(time) term';
ods output ParameterEstimates = synwork._pe_tv
           FitStatistics      = synwork._fit_tv;
proc phreg data=synwork.analysis;
  model time * event(0) = chf chf_lt / ties=efron risklimits gconv=1e-12;
  chf_lt = chf * log(time);
run;
title;
%pe_check(ds=synwork._pe_tv, tag=cox_tv, param=chf,    name=chf,    hr=0)
%pe_check(ds=synwork._pe_tv, tag=cox_tv, param=chf_lt, name=chf_lt, hr=0)
%m2ll_check(ds=synwork._fit_tv, tag=cox_tv)

/*------------------------------------------------------------------------------
  Collected results table.
------------------------------------------------------------------------------*/
data synwork.cox_results;
  /* PHREG sizes Label to fit each model's longest label. Declaring one length
     here before SET prevents the multiple-lengths warning. Label is not kept. */
  length model $40 Parameter $20 ClassVal0 $16 Label $256;
  set synwork._pe_unadj  (in=a)
      synwork._pe_unadj_b(in=b)
      synwork._pe_adj    (in=c)
      synwork._pe_int    (in=d)
      synwork._pe_tv     (in=e);
  if a then model = 'Unadjusted, Efron';
  else if b then model = 'Unadjusted, Breslow';
  else if c then model = 'Adjusted, Efron';
  else if d then model = 'Adjusted with interactions, Efron';
  else model = 'CHF by log(time), Efron';
  /* Wald limits computed here for every row, including models fitted
     without RISKLIMITS. */
  hr     = exp(Estimate);
  hr_lcl = exp(Estimate - probit(0.975) * StdErr);
  hr_ucl = exp(Estimate + probit(0.975) * StdErr);
  keep model Parameter ClassVal0 Estimate StdErr ProbChiSq hr hr_lcl hr_ucl;
  format Estimate StdErr 10.6 hr hr_lcl hr_ucl 8.4 ProbChiSq pvalue6.4;
run;
%export_csv(data=synwork.cox_results, name=cox_models)

/* Hazard ratios and 95 percent Wald limits for every coefficient that has a
   Python reference value but was checked above only on the log scale. They
   are exp() of values already checked, and are checked here so that the
   numbers in results/tables/cox_models.csv are covered too. */
data _null_;
  set synwork.cox_results;
  length tag $24 pname $24;
  select (model);
    when ('Unadjusted, Breslow')               tag = 'cox_unadj_breslow';
    when ('Adjusted, Efron')                   tag = 'cox_adj';
    when ('Adjusted with interactions, Efron') tag = 'cox_int';
    when ('CHF by log(time), Efron')           tag = 'cox_tv';
    otherwise tag = '';
  end;
  pname = lowcase(translate(strip(Parameter), '_', '*'));
  if pname = 'race' then pname = cats('race_', ifc(ClassVal0 = 'Black', 'black',
                                         ifc(ClassVal0 = 'Other', 'other', 'hisp')));
  /* Rows checked with hazard ratios already, or with no reference value. */
  if tag = 'cox_adj' and pname = 'chf' then tag = '';
  if tag = 'cox_int' and pname not in ('chf' 'chf_agec' 'chf_female') then tag = '';
  if tag ne '' then do;
    call execute(cats('%nrstr(%check)(metric=', tag, '_', pname, '_hr, actual=',
                      put(hr, best16.), ', tol=0.00001)'));
    call execute(cats('%nrstr(%check)(metric=', tag, '_', pname, '_hrlcl, actual=',
                      put(hr_lcl, best16.), ', tol=0.00001)'));
    call execute(cats('%nrstr(%check)(metric=', tag, '_', pname, '_hrucl, actual=',
                      put(hr_ucl, best16.), ', tol=0.00001)'));
  end;
run;

/*------------------------------------------------------------------------------
  5.2 Adjusted survival curves. Two covariate profiles that differ only in
  CHF: a woman of median age with median 2008 reimbursement, White, and none
  of the other conditions. The BASELINE statement gives the model-based
  survival curve and 95 percent limits for each profile.
------------------------------------------------------------------------------*/
proc means data=synwork.analysis noprint;
  var age logreimb;
  output out=synwork._med median=age logreimb;
run;

data synwork._covs;
  set synwork._med(keep=age logreimb);
  length profile $40;
  female = 1; race = 1; esrd = 0; alzh = 0; ckd = 0; cancer = 0; copd = 0;
  depress = 0; diab = 0; ihd = 0; osteo = 0; raoa = 0; stroke = 0;
  do chf = 0, 1;
    profile = ifc(chf = 1, 'CHF', 'No CHF');
    output;
  end;
  format race racef. chf chff.;
run;

title 'Covariate profiles for the adjusted survival curves';
proc print data=synwork._covs noobs;
run;
title;

proc phreg data=synwork.analysis plots=none;
  class race(ref='White') / param=ref;
  model time * event(0) = &adjvars / ties=efron gconv=1e-12;
  baseline covariates=synwork._covs out=synwork._adjsurv
           survival=surv lower=surv_lcl upper=surv_ucl / rowid=profile;
run;

data synwork._adjsurv;
  set synwork._adjsurv;
  months = time * 24 / &admin_day;
  label months = 'Months since January 1, 2009';
run;

%fig(adjusted_survival_chf)
title 'Adjusted survival by CHF status (woman of median age, median reimbursement)';
proc sgplot data=synwork._adjsurv;
  band x=months lower=surv_lcl upper=surv_ucl / type=step group=chf transparency=0.75;
  step x=months y=surv / group=chf name='s';
  xaxis values=(0 to 24 by 3) label='Months since January 1, 2009';
  yaxis min=0.95 max=1 label='Model-based survival probability';
  keylegend 's' / title='';
run;
title;
%fig_end
