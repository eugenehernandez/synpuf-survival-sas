/*==============================================================================
  Program   : 10_sensitivity.sas
  Project   : synpuf-survival-sas
  Purpose   : Design sensitivity analyses. The adjusted Cox model of
              06_cox.sas is fitted in five study designs, and the CHF hazard
              ratio from each is shown in one table and one figure next to
              the primary result.
                Primary        Landmark January 1, 2009, 65 and older,
                               fee-for-service in 2008 (61,184).
                Naive          Time zero January 1, 2008, with the 2008 flags
                               and follow-up through December 31, 2010. Fitted
                               three ways (unadjusted, adjusted without 2008
                               reimbursement, fully adjusted).
                Coverage       Primary cohort with coverage-based censoring,
                               excluding the people censored at day 0 (55,075).
                Full 65+       Landmark design without the fee-for-service
                               restriction (95,730).
                All ages       Landmark design, every age, with an indicator
                               for younger than 65 (114,538).
  Inputs    : SYNWORK.ANALYSIS, SYNWORK.COHORT_COV, SYNWORK.COHORT_65,
              SYNWORK.COHORT_ALLAGES, SYNWORK.LANDMARK_ALL (naive design)
  Outputs   : SYNWORK.SENS_RESULTS  results/tables/sensitivity_chf.csv
              SYNWORK.SENS_NAIVE_MECHANISM
                        results/tables/sensitivity_naive_mechanism.csv
              results/figures/sensitivity_chf.png
  Seminar   : extends section 5 (no seminar counterpart)
  Notes     : The naive design measures CHF and the other conditions over
              2008 and starts follow-up on January 1, 2008, so the condition
              flags partly come from claims made during follow-up (look-ahead).
              In claims data with real clinical signal, a person who is
              diagnosed shortly before death in 2008 is counted as exposed
              from day 1, which biases the CHF hazard ratio. The landmark
              design avoids this. The naive cohort uses the same eligibility
              as the primary cohort, measured at January 1, 2008 (65 or older
              on that date and fee-for-service in all 12 months of 2008).
              Profiling found no loss to follow-up, and the program checks
              that this holds for the naive cohort as well.
              In this synthetic release death is unrelated to CHF. Four designs
              give a hazard ratio near 1. The fully adjusted naive model gives
              a modest spurious effect, and section 4 shows why.
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

/* Naive design time scale. Day 1 is January 1, 2008 and administrative
   censoring on December 31, 2010 is day 1,096. */
%let naive_origin = '31DEC2007'd;
%let naive_start  = '01JAN2008'd;
%let naive_admin  = 1096;

/*------------------------------------------------------------------------------
  1. Naive cohort. Everyone in the 2008 file was alive on January 1, 2008.
------------------------------------------------------------------------------*/
data synwork.cohort_naive;
  set synwork.landmark_all;
  age_int08 = intck('year', BENE_BIRTH_DT, &naive_start, 'c');
  if age_int08 >= 65 and ffs08 = 1;
  age = (&naive_start - BENE_BIRTH_DT) / 365.25;
  event = (&naive_start <= death_dt <= '31DEC2010'd);
  if event then time = death_dt - &naive_origin;
  else time = &naive_admin;
  /* Alive with no death date but missing from a later file. */
  ltfu_naive = (missing(death_dt) and not (in2009 and in2010));
  keep DESYNPUF_ID time event ltfu_naive chf age female race esrd alzh ckd cancer
       copd depress diab ihd osteo raoa stroke logreimb;
  label age = 'Age on January 1, 2008 (years)'
        time = 'Days from December 31, 2007 to death or censoring';
run;

proc sql noprint;
  select sum(ltfu_naive) into :_nv_ltfu trimmed from synwork.cohort_naive;
quit;
%check(metric=sens_naive_ltfu, actual=&_nv_ltfu)

/*------------------------------------------------------------------------------
  2. Fit the adjusted model in one design and keep one result row.
     Race enters as indicators, which gives the same model as the CLASS
     statement with White as the reference in 06.
------------------------------------------------------------------------------*/
%macro sens_fit(tag=, label=, data=, time=time, event=event, extra=, ord=,
                covs=&dxvars);
  data synwork._sv_&tag / view=synwork._sv_&tag;
    set &data;
    race_black = (race = 2);
    race_other = (race = 3);
    race_hisp  = (race = 5);
    %if %length(&extra) %then %do;
      under65 = (age_int < 65);
    %end;
  run;

  title "Adjusted Cox model, &label";
  ods output ParameterEstimates = synwork._spe_&tag FitStatistics = synwork._sfit_&tag;
  proc phreg data=synwork._sv_&tag;
    model &time * &event(0) = &covs &extra / ties=efron gconv=1e-12 risklimits;
  run;
  title;

  proc sql;
    create table synwork._srow_&tag as
    select &ord as ord, "&label" as design length=48,
           (select count(*) from synwork._sv_&tag) as n,
           (select sum(&event) from synwork._sv_&tag) as deaths,
           Estimate as beta, StdErr as se,
           exp(Estimate) as hr,
           exp(Estimate - probit(0.975) * StdErr) as hr_lcl,
           exp(Estimate + probit(0.975) * StdErr) as hr_ucl,
           ProbChiSq as p
    from synwork._spe_&tag where upcase(Parameter) = 'CHF';
  quit;

  /* Reference checks. */
  %local n ev b s h l u;
  proc sql noprint;
    select n, deaths, beta format=best16., se format=best16., hr format=best16.,
           hr_lcl format=best16., hr_ucl format=best16.
      into :n trimmed, :ev trimmed, :b trimmed, :s trimmed, :h trimmed, :l trimmed, :u trimmed
    from synwork._srow_&tag;
  quit;
  %check(metric=sens_&tag._n,         actual=&n)
  %check(metric=sens_&tag._events,    actual=&ev)
  %check(metric=sens_&tag._chf_beta,  actual=&b, tol=0.00001)
  %check(metric=sens_&tag._chf_se,    actual=&s, tol=0.00001)
  %check(metric=sens_&tag._chf_hr,    actual=&h, tol=0.00001)
  %check(metric=sens_&tag._chf_hrlcl, actual=&l, tol=0.00001)
  %check(metric=sens_&tag._chf_hrucl, actual=&u, tol=0.00001)
  %m2ll_check(ds=synwork._sfit_&tag, tag=sens_&tag)
  %if %length(&extra) %then %do;
    %pe_check(ds=synwork._spe_&tag, tag=sens_&tag, param=under65, name=under65, hr=0)
  %end;
%mend sens_fit;

/* -2 log partial likelihood and single-parameter checks, as in 06. */
%macro m2ll_check(ds=, tag=);
  %local v;
  proc sql noprint;
    select WithCovariates format=best16. into :v trimmed
    from &ds where compress(upcase(Criterion)) = '-2LOGL';
  quit;
  %check(metric=&tag._m2ll, actual=&v, tol=0.001)
%mend m2ll_check;

%macro pe_check(ds=, tag=, param=, name=, hr=0, tol=0.00001);
  %local b s;
  proc sql noprint;
    select Estimate format=best16., StdErr format=best16. into :b trimmed, :s trimmed
    from &ds where upcase(Parameter) = upcase("&param");
  quit;
  %check(metric=&tag._&name._beta, actual=&b, tol=&tol)
  %check(metric=&tag._&name._se,   actual=&s, tol=&tol)
%mend pe_check;

/* The adjusted covariates without 2008 reimbursement. */
%let dx_noreimb = chf age female race_black race_other race_hisp esrd alzh ckd
                  cancer copd depress diab ihd osteo raoa stroke;

/* The naive design appears three times, so the table shows how the CHF
   hazard ratio moves as covariates are added: unadjusted, adjusted without
   2008 reimbursement, and fully adjusted. */
%sens_fit(tag=primary,       ord=1, data=synwork.analysis,
          label=%str(Primary: landmark 2009, 65+, fee-for-service))
%sens_fit(tag=naive_unadj,   ord=2, data=synwork.cohort_naive, covs=chf,
          label=%str(Naive 2008: unadjusted))
%sens_fit(tag=naive_noreimb, ord=3, data=synwork.cohort_naive, covs=&dx_noreimb,
          label=%str(Naive 2008: adjusted, without reimbursement))
%sens_fit(tag=naive,         ord=4, data=synwork.cohort_naive,
          label=%str(Naive 2008: adjusted, with reimbursement))
%sens_fit(tag=cov,           ord=5, data=synwork.cohort_cov, time=time_cov, event=event_cov,
          label=%str(Coverage-based censoring))
%sens_fit(tag=full65,        ord=6, data=synwork.cohort_65,
          label=%str(All 65+, no fee-for-service restriction))
%sens_fit(tag=allages,       ord=7, data=synwork.cohort_allages, extra=under65,
          label=%str(All ages, with under-65 indicator))

/*------------------------------------------------------------------------------
  3. Comparison table and figure.
------------------------------------------------------------------------------*/
data synwork.sens_results;
  set synwork._srow_primary synwork._srow_naive_unadj synwork._srow_naive_noreimb
      synwork._srow_naive synwork._srow_cov synwork._srow_full65 synwork._srow_allages;
  format beta se 10.6 hr hr_lcl hr_ucl 8.4 p pvalue6.4 n deaths comma9.;
  label design = 'Design' n = 'People' deaths = 'Deaths'
        hr = 'CHF hazard ratio' hr_lcl = 'Lower 95% limit' hr_ucl = 'Upper 95% limit';
run;

title 'CHF hazard ratio (adjusted) across study designs';
proc print data=synwork.sens_results noobs label;
  var design n deaths hr hr_lcl hr_ucl p;
run;
title;
%export_csv(data=synwork.sens_results, name=sensitivity_chf)

%fig(sensitivity_chf)
title 'Adjusted CHF hazard ratio by study design';
proc sgplot data=synwork.sens_results noautolegend;
  refline 1 / axis=x lineattrs=(pattern=dot);
  scatter y=design x=hr / xerrorlower=hr_lcl xerrorupper=hr_ucl
          markerattrs=(symbol=diamondfilled size=9);
  xaxis type=log logbase=10 min=0.7 max=1.5 label='CHF hazard ratio (95% CI, log scale)';
  yaxis discreteorder=data reverse display=(nolabel);
run;
title;
%fig_end

/*------------------------------------------------------------------------------
  4. Why the naive adjusted model moves. In the naive design 2008
     reimbursement is measured during the first year of follow-up. People who
     die early in 2008 have less time to accumulate claims, so the outcome
     partly causes the covariate. Adjusting for it biases the coefficients of
     covariates correlated with it, and CHF is strongly correlated with
     reimbursement. The table gives the numbers behind that explanation.
------------------------------------------------------------------------------*/
data synwork._nvm / view=synwork._nvm;
  set synwork.cohort_naive;
  died08 = (event = 1 and time <= 366);
run;

proc sql;
  create table synwork.sens_naive_mechanism as
  select 'All' as chf_group length=8, died08, count(*) as n,
         mean(logreimb) as mean_logreimb
    from synwork._nvm group by died08
  union all
  select ifc(chf = 1, 'CHF', 'No CHF'), died08, count(*), mean(logreimb)
    from synwork._nvm group by chf, died08
  order by chf_group, died08;

  create table synwork._nvrate as
  select chf, mean(died08) as died08_rate from synwork._nvm group by chf;
quit;

proc corr data=synwork._nvm noprint outp=synwork._nvcorr;
  var chf logreimb;
run;

data synwork.sens_naive_mechanism;
  set synwork.sens_naive_mechanism;
  format mean_logreimb 8.3;
  label died08 = 'Died in 2008' mean_logreimb = 'Mean log(1 + 2008 reimbursement)';
run;
title 'Naive design: 2008 reimbursement by death in 2008 and CHF status';
proc print data=synwork.sens_naive_mechanism noobs label;
run;
title;
%export_csv(data=synwork.sens_naive_mechanism, name=sensitivity_naive_mechanism)

%macro mech_checks;
  %local g d v n r;
  %do g = 1 %to 3;
    %do d = 0 %to 1;
      proc sql noprint;
        select mean_logreimb format=best16. into :v trimmed
        from synwork.sens_naive_mechanism
        where chf_group = "%scan(All|CHF|No CHF, &g, |)" and died08 = &d;
      quit;
      %let n = %scan(all chf1 chf0, &g);
      %check(metric=sens_naive_mlogreimb_&n._died&d, actual=&v, tol=0.000001)
    %end;
  %end;
  %do d = 0 %to 1;
    proc sql noprint;
      select died08_rate format=best16. into :r trimmed from synwork._nvrate where chf = &d;
    quit;
    %check(metric=sens_naive_died08_rate_chf&d, actual=&r, tol=0.000001)
  %end;
  proc sql noprint;
    select logreimb format=best16. into :r trimmed from synwork._nvcorr
    where _TYPE_ = 'CORR' and upcase(_NAME_) = 'CHF';
  quit;
  %check(metric=sens_naive_corr_chf_logreimb, actual=&r, tol=0.000001)
%mend mech_checks;
%mech_checks
