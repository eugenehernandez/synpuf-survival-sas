/*==============================================================================
  Program   : 08_ph_assumption.sas
  Project   : synpuf-survival-sas
  Purpose   : Checks of the proportional hazards assumption for the adjusted
              Cox model, on the real cohort and the plasmode data side by side.
              The plasmode CHF hazard ratio falls over time by design, so the
              checks should flag CHF there and not on the real data.
  Inputs    : SYNWORK.ANALYSIS, SYNWORK.PLASMODE (through views from %dx_views)
  Outputs   : SYNWORK.PH_TV_TESTS  results/tables/ph_time_interactions.csv
              SYNWORK.PH_PIECE     results/tables/ph_piecewise_chf.csv
              SYNWORK.PH_SUPTEST   results/tables/ph_assess_suptest.csv
              results/figures/ph_loglog_chf.png
              results/figures/ph_schoenfeld_chf.png
              results/figures/ph_piecewise_chf.png
              results/figures/ph_score_chf_real.png
              results/figures/ph_score_chf_plasmode.png
  Seminar   : section 7 (proportional hazards: log-log survival plots,
              Schoenfeld residuals, covariate by time interactions, ASSESS PH)
  Notes     : Every part uses the full cohort except ASSESS PH, which uses the
              same fixed subsample of &diag_n people as 07, for the same
              reason (1,000 resamples over every person and a score process
              plot with one point per person). ASSESS is run with Breslow ties,
              and its p-values are reported without a Python check because
              they come from resampling inside PROC PHREG.
              Death dates fall on the first of the month, so there are 24
              event times. The Schoenfeld residuals are therefore summarized
              as a mean for each month rather than smoothed person by person.
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

%dx_views

/* Covariates other than CHF, for models that replace the CHF term. */
%let dx_nochf = age female race_black race_other race_hisp esrd alzh ckd cancer
                copd depress diab ihd osteo raoa stroke logreimb;

/*------------------------------------------------------------------------------
  1. Log of the negative log of Kaplan-Meier survival by CHF, against log
     time. Under proportional hazards the two curves are parallel.
------------------------------------------------------------------------------*/
%macro loglog(src);
  proc lifetest data=synwork._dx_&src outsurv=synwork._ll_&src noprint plots=none;
    time time * event(0);
    strata chf;
  run;
  data synwork._ll_&src;
    length source $8;
    set synwork._ll_&src;
    source = ifc("&src" = 'real', 'Real', 'Plasmode');
    if time > 0 and 0 < survival < 1;
    lls  = log(-log(survival));
    lmon = log(time * 24 / &admin_day);
    keep source chf time lls lmon;
  run;
%mend loglog;
%loglog(real)
%loglog(plas)

data synwork._ll;
  set synwork._ll_real synwork._ll_plas;
  label lls = 'log(-log(survival))' lmon = 'log(months since January 1, 2009)';
  format chf chff.;
run;

%fig(ph_loglog_chf)
title 'Log-log survival by CHF status';
proc sgpanel data=synwork._ll;
  panelby source / novarname columns=2 sort=descending;
  step x=lmon y=lls / group=chf lineattrs=(thickness=2) name='s';
  colaxis label='log(months since January 1, 2009)';
  rowaxis label='log(-log(survival))';
  keylegend 's' / title='';
run;
title;
%fig_end

/*------------------------------------------------------------------------------
  2. Scaled Schoenfeld residuals for CHF. The PH estimate plus the scaled
     residual approximates the CHF log hazard ratio at each event time
     (Grambsch and Therneau, 1994). A flat pattern supports proportional
     hazards. The plasmode panel also shows the true value at the middle of
     each month.
------------------------------------------------------------------------------*/
%macro schoen(src);
  ods output ParameterEstimates = synwork._phpe_&src;
  proc phreg data=synwork._dx_&src;
    model time * event(0) = &dxvars / ties=efron gconv=1e-12;
    output out=synwork._sch_&src wtressch=%prefix(ws_, &dxvars);
  run;
  proc sql noprint;
    select Estimate format=best16. into :_bchf trimmed
    from synwork._phpe_&src where upcase(Parameter) = 'CHF';
  quit;
  proc sql;
    create table synwork._schm_&src as
    select source, time,
           1 + intck('month', &landmark_dt, &origin_dt + time) as month,
           count(*) as deaths,
           &_bchf + mean(ws_chf) as beta_t,
           &_bchf + mean(ws_chf) - probit(0.975) * std(ws_chf) / sqrt(count(*)) as lcl,
           &_bchf + mean(ws_chf) + probit(0.975) * std(ws_chf) / sqrt(count(*)) as ucl,
           &_bchf as beta_ph,
           case when source = 'Plasmode'
                then &plas_chf_a + &plas_chf_b * log(time + 15) else . end as truth
    from synwork._sch_&src
    where event = 1
    group by source, time;
  quit;
%mend schoen;
%schoen(real)
%schoen(plas)

data synwork._schm;
  set synwork._schm_real synwork._schm_plas;
  label beta_t = 'PH estimate plus mean scaled Schoenfeld residual'
        beta_ph = 'PH estimate' truth = 'True value (middle of month)';
run;

%fig(ph_schoenfeld_chf)
title 'CHF log hazard ratio by month of death, from scaled Schoenfeld residuals';
proc sgpanel data=synwork._schm;
  panelby source / novarname columns=2 sort=descending;
  scatter x=month y=beta_t / yerrorlower=lcl yerrorupper=ucl
          markerattrs=(symbol=circlefilled size=5) legendlabel='Monthly mean (95% CI)' name='m';
  series x=month y=beta_ph / lineattrs=(pattern=dot thickness=2)
          legendlabel='Proportional hazards estimate' name='p';
  series x=month y=truth / lineattrs=(pattern=dash color=black thickness=2)
          legendlabel='True value (plasmode)' name='t';
  refline 0 / axis=y lineattrs=(pattern=shortdash color=gray);
  colaxis label='Month of follow-up (1 = January 2009)' values=(1 to 24 by 1)
          valueattrs=(size=7);
  rowaxis label='CHF log hazard ratio';
  keylegend 'm' 'p' 't' / title='';
run;
title;
%fig_end

/*------------------------------------------------------------------------------
  3. Covariate by log(time) interactions, one covariate at a time, each added
     to a model with a programming statement. A coefficient different from 0
     means that covariate's hazard ratio changes with time. Two versions.
       Naive                Both data sets. The adjusted working model
                            (linear age, constant CHF effect), with
                            covariate x log(time). This is what an analyst
                            without knowledge of the truth would fit.
       Correctly specified  Plasmode only. The true model form (age minus 75
                            and its square, and chf_lt = chf x log(time + 15)
                            for the decaying CHF effect), with covariate x
                            log(time + 15). Only here can a flag be judged,
                            because only here is the truth known.
     In the naive version the plasmode data flag diabetes even though its
     true hazard ratio is constant. Part of the misfit of the omitted age
     curvature and CHF decay is absorbed by any covariate that is correlated
     with them and given a time-varying term, and log(time) on dates set to
     the first of the month distorts early event times. The correctly
     specified version removes both sources.
------------------------------------------------------------------------------*/
%let dx_true = chf chf_lt agec agec2 female race_black race_other race_hisp esrd
               alzh ckd cancer copd depress diab ihd osteo raoa stroke logreimb;

%macro tv_test(src, v, version=Naive);
  %local lab;
  %if &src = real %then %let lab = Real;
  %else %let lab = Plasmode;
  ods output ParameterEstimates = synwork._tvpe;
  proc phreg data=synwork._dx_&src;
    %if &version = Naive %then %do;
      model time * event(0) = &dxvars tv_lt / ties=efron gconv=1e-12;
      tv_lt = &v * log(time);
    %end;
    %else %do;
      model time * event(0) = &dx_true tv_lt / ties=efron gconv=1e-12;
      chf_lt = chf * log(time + 15);
      tv_lt  = &v * log(time + 15);
    %end;
  run;
  data synwork._tvrow;
    length source $8 version $20 covariate $12 time_term $16;
    set synwork._tvpe(where=(upcase(Parameter) = 'TV_LT'));
    source    = "&lab";
    version   = ifc("&version" = 'Naive', 'Naive', 'Correctly specified');
    covariate = "&v";
    time_term = ifc("&version" = 'Naive', 'log(time)', 'log(time + 15)');
    z = Estimate / StdErr;
    keep source version covariate time_term Estimate StdErr z ChiSq ProbChiSq;
  run;
  proc append base=synwork.ph_tv_tests data=synwork._tvrow;
  run;
%mend tv_test;

title 'Covariate by log(time) interactions';
%macro tv_all;
  %local s i v;
  /* Rows are appended, so start from an empty table. */
  %if %sysfunc(exist(synwork.ph_tv_tests)) %then %do;
    proc datasets library=synwork nolist;
      delete ph_tv_tests;
    quit;
  %end;
  %do s = 1 %to 2;
    %do i = 1 %to 4;
      %let v = %scan(chf age female diab, &i);
      %tv_test(%scan(real plas, &s), &v)
    %end;
  %end;
  /* CHF is not tested in the correct version, because its time-varying
     term is already part of the true model. */
  %do i = 1 %to 3;
    %let v = %scan(agec female diab, &i);
    %tv_test(plas, &v, version=Correct)
  %end;
%mend tv_all;
%tv_all
title;

data synwork.ph_tv_tests;
  set synwork.ph_tv_tests;
  label Estimate = 'Coefficient of covariate x time term' z = 'Estimate / SE';
  format z 8.2;
run;
title 'Covariate by time interaction terms (one model per row)';
proc print data=synwork.ph_tv_tests noobs label;
run;
title;
%export_csv(data=synwork.ph_tv_tests, name=ph_time_interactions)

%macro tv_checks(source, pre, version=Naive, vars=chf age female diab, tag=ph_tv);
  %local i v b s;
  %let i = 1;
  %let v = %scan(&vars, &i);
  %do %while (%length(&v) > 0);
    proc sql noprint;
      select Estimate format=best16., StdErr format=best16. into :b trimmed, :s trimmed
      from synwork.ph_tv_tests
      where source = "&source" and covariate = "&v" and version = "&version";
    quit;
    %check(metric=&pre.&tag._&v._beta, actual=&b, tol=0.00001)
    %check(metric=&pre.&tag._&v._se,   actual=&s, tol=0.00001)
    %let i = %eval(&i + 1);
    %let v = %scan(&vars, &i);
  %end;
%mend tv_checks;
%tv_checks(Real, )
%tv_checks(Plasmode, plas_)
%tv_checks(Plasmode, plas_, version=Correctly specified, vars=agec female diab, tag=ph_tvc)

/* The covariates simulated with a constant hazard ratio should not be
   flagged once the model has the true form (|z| below 1.96). */
%macro tvc_null_checks;
  %local i v z;
  %do i = 1 %to 3;
    %let v = %scan(agec female diab, &i);
    proc sql noprint;
      select z format=best16. into :z trimmed from synwork.ph_tv_tests
      where source = 'Plasmode' and covariate = "&v" and version = 'Correctly specified';
    quit;
    %check_range(metric=plas_ph_tvc_&v._z, actual=&z, lo=-1.96, hi=1.96)
  %end;
%mend tvc_null_checks;
%tvc_null_checks

/*------------------------------------------------------------------------------
  4. Piecewise CHF hazard ratio. Separate CHF coefficients for months 1 to 6
     (days 1 to 181), 7 to 12 (days 182 to 365), and 13 to 24 (days 366 to
     730), created by programming statements. For the plasmode data the true
     log hazard ratio across each period runs from its value in the middle of
     the first month to its value in the middle of the last month.
------------------------------------------------------------------------------*/
%macro piece(src);
  title "Piecewise CHF effect (&src)";
  ods output ParameterEstimates = synwork._pc_&src;
  proc phreg data=synwork._dx_&src;
    model time * event(0) = chf_p1 chf_p2 chf_p3 &dx_nochf / ties=efron gconv=1e-12;
    chf_p1 = chf * (time <= 181);
    chf_p2 = chf * (181 < time <= 365);
    chf_p3 = chf * (time > 365);
  run;
  title;
  data synwork._pc_&src;
    length source $8 period $12;
    set synwork._pc_&src(where=(upcase(Parameter) in ('CHF_P1' 'CHF_P2' 'CHF_P3')));
    source = ifc("&src" = 'real', 'Real', 'Plasmode');
    p = input(substr(Parameter, 6, 1), 1.);
    period = choosec(p, 'Months 1-6', 'Months 7-12', 'Months 13-24');
    hr     = exp(Estimate);
    hr_lcl = exp(Estimate - probit(0.975) * StdErr);
    hr_ucl = exp(Estimate + probit(0.975) * StdErr);
    if source = 'Plasmode' then do;
      /* Middle of the first and last month of each period (day of the 16th). */
      t_first = choosen(p, 16, 197, 381);
      t_last  = choosen(p, 166, 349, 715);
      true_hr_first = exp(&plas_chf_a + &plas_chf_b * log(t_first));
      true_hr_last  = exp(&plas_chf_a + &plas_chf_b * log(t_last));
    end;
    keep source p period Estimate StdErr hr hr_lcl hr_ucl true_hr_first true_hr_last;
  run;
%mend piece;
%piece(real)
%piece(plas)

data synwork.ph_piece;
  set synwork._pc_real synwork._pc_plas;
  format hr hr_lcl hr_ucl true_hr_first true_hr_last 8.3;
run;
title 'CHF hazard ratio by period of follow-up';
proc print data=synwork.ph_piece noobs;
run;
title;
%export_csv(data=synwork.ph_piece, name=ph_piecewise_chf)

%macro piece_checks(source, pre);
  %local p b s;
  %do p = 1 %to 3;
    proc sql noprint;
      select Estimate format=best16., StdErr format=best16. into :b trimmed, :s trimmed
      from synwork.ph_piece where source = "&source" and p = &p;
    quit;
    %check(metric=&pre.ph_piece_p&p._beta, actual=&b, tol=0.00001)
    %check(metric=&pre.ph_piece_p&p._se,   actual=&s, tol=0.00001)
  %end;
%mend piece_checks;
%piece_checks(Real, )
%piece_checks(Plasmode, plas_)

%fig(ph_piecewise_chf)
title 'CHF hazard ratio by period of follow-up';
proc sgpanel data=synwork.ph_piece;
  panelby source / novarname columns=2 sort=descending;
  highlow x=period low=true_hr_last high=true_hr_first / type=bar barwidth=0.35
          fillattrs=(color=lightgray) legendlabel='True range (plasmode)' name='t';
  scatter x=period y=hr / yerrorlower=hr_lcl yerrorupper=hr_ucl
          markerattrs=(symbol=diamondfilled size=9) legendlabel='Estimate (95% CI)' name='e';
  refline 1 / axis=y lineattrs=(pattern=dot);
  colaxis label='Period of follow-up' discreteorder=data;
  rowaxis type=log logbase=10 label='CHF hazard ratio (log scale)';
  keylegend 'e' 't' / title='';
run;
title;
%fig_end

/*------------------------------------------------------------------------------
  5. ASSESS PH on the fixed subsample. The table covers every covariate
     of the adjusted model. ODS Graphics is off for that step, because it
     would draw one score process plot per covariate. The figure comes from a
     second fit with CHF alone, which draws a single plot.
------------------------------------------------------------------------------*/
%macro ph_assess(src, name);
  %hash_sample(data=synwork._dx_&src, out=synwork._sub_&src, n=&diag_n)

  ods graphics off;
  title "Supremum tests for proportional hazards (&name data, n = &diag_n)";
  ods output ProportionalHazardsSupTest = synwork._phsup_&src;
  proc phreg data=synwork._sub_&src;
    model time * event(0) = &dxvars / ties=breslow;
    assess ph / resample seed=%eval(&seed + 5);
  run;
  title;
  %graphics_on

  %fig(ph_score_chf_&name)
  title "Standardized score process for CHF (&name data, n = &diag_n)";
  proc phreg data=synwork._sub_&src;
    model time * event(0) = chf / ties=breslow;
    assess ph / resample seed=%eval(&seed + 6);
  run;
  title;
  %fig_end
%mend ph_assess;
/* Skipped when run_assess = 0 in 00_config.sas (see the note there). */
%macro ph_assess_all;
  %if &run_assess = 0 %then %do;
    %put NOTE: [assess] run_assess = 0, proportional hazards supremum tests skipped.;
    %return;
  %end;
  %ph_assess(real, real)
  %ph_assess(plas, plasmode)

  data synwork.ph_suptest;
    length source $8;
    set synwork._phsup_real(in=a) synwork._phsup_plas;
    source = ifc(a, 'Real', 'Plasmode');
  run;
  title "Supremum tests for proportional hazards, adjusted model (subsample n = &diag_n)";
  proc print data=synwork.ph_suptest noobs;
  run;
  title;
  %export_csv(data=synwork.ph_suptest, name=ph_assess_suptest)
%mend ph_assess_all;
%ph_assess_all
