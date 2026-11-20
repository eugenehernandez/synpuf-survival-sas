/*==============================================================================
  Program   : 05_nonparametric.sas
  Project   : synpuf-survival-sas
  Purpose   : Nonparametric survival analysis of the primary cohort by CHF
              status. Kaplan-Meier survival, Nelson-Aalen cumulative hazard,
              median and restricted mean survival time, tests that compare the
              two curves, and a kernel-smoothed hazard.
  Inputs    : SYNWORK.ANALYSIS
  Outputs   : results/tables/km_estimates.csv     survival and cumulative hazard at 365 and 730 days
              results/tables/rmst.csv             restricted mean survival time to day 730
              results/tables/group_tests.csv      log-rank and weighted tests
              results/figures/km_chf_lifetest.png  PROC LIFETEST plot with at-risk table
              results/figures/km_chf.png          README version with a zoomed y axis
              results/figures/cumhaz_chf.png      Nelson-Aalen cumulative hazard
              results/figures/hazard_chf.png      kernel-smoothed hazard
  Seminar   : section 3 (3.1 Kaplan-Meier, 3.2 Nelson-Aalen, 3.3 survival
              times, 3.4 comparing groups, 3.5 hazard function)
  Notes     : Death dates resolve only to the month, so the 61,184 records hold
              at most 25 distinct times per group (24 death months plus day
              730). The data are collapsed to one row per CHF group, time, and
              status with a FREQ count. Every estimate and test is identical to
              the row-level analysis, the output tables stay short, and the
              procedures run in seconds. The collapsed data are used only in
              this program. From 06_cox.sas onward every model uses the
              row-level SYNWORK.ANALYSIS, because residuals and influence
              statistics in 07 to 09 are computed per person.
              Estimates are reported on the day scale (day 365 is the end of
              2009 and day 730 the end of 2010). Plots use a month scale,
              MONTHS = TIME x 24 / 730, so axis ticks fall on whole months and
              day 365 is month 12.
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

/*------------------------------------------------------------------------------
  0. Collapse to one row per CHF group, time, and status.
------------------------------------------------------------------------------*/
proc sql;
  create table synwork._km as
  select chf, time, event, count(*) as count
  from synwork.analysis
  group by chf, time, event
  order by chf, time, event;
quit;

data synwork._km;
  set synwork._km;
  months = time * 24 / &admin_day;
  format chf chff.;
  label chf    = 'CHF status'
        time   = 'Days since December 31, 2008'
        months = 'Months since January 1, 2009'
        count  = 'Number of beneficiaries';
run;

/*------------------------------------------------------------------------------
  3.1 and 3.2. Kaplan-Meier and Nelson-Aalen estimates at 365 and 730 days,
  plus 3.4, tests comparing CHF groups.
    Log-rank weights every event time equally and is most powerful when the
    hazards are proportional. Wilcoxon (Gehan) weights by the number at risk,
    so early differences count more. Tarone-Ware (square root of the number
    at risk) sits between the two. Peto and modified Peto weight by a survival
    estimate and are less sensitive to differences in censoring between
    groups. Fleming-Harrington (0,1) weights late differences more heavily.
------------------------------------------------------------------------------*/
ods output ProductLimitEstimates = synwork._km_tl
           Quartiles             = synwork._quart
           HomTests              = synwork._tests;
title 'Kaplan-Meier and Nelson-Aalen estimates at 365 and 730 days, by CHF';
proc lifetest data=synwork._km timelist=365 730 nelson plots=none;
  time time * event(0);
  freq count;
  strata chf / test=(logrank wilcoxon tarone peto modpeto fleming(0,1));
run;

/* Same estimates for the whole cohort, for the README and the cross-check. */
ods output ProductLimitEstimates = synwork._km_tl_all;
title 'Kaplan-Meier and Nelson-Aalen estimates at 365 and 730 days, whole cohort';
proc lifetest data=synwork._km timelist=365 730 nelson plots=none;
  time time * event(0);
  freq count;
run;
title;

data synwork.km_estimates;
  length group $8;
  set synwork._km_tl(in=a) synwork._km_tl_all;
  if a then group = ifc(chf = 1, 'chf1', 'chf0');
  else group = 'all';
  keep group Timelist Survival StdErr CumHaz StdErrCumHaz Left;
run;

title 'Survival and cumulative hazard at 365 and 730 days';
proc print data=synwork.km_estimates noobs;
run;
title;
%export_csv(data=synwork.km_estimates, name=km_estimates)

%macro km_checks;
  %local g t grp tt s se h;
  %do g = 1 %to 3;
    %let grp = %scan(chf1 chf0 all, &g);
    %do t = 1 %to 2;
      %let tt = %scan(365 730, &t);
      proc sql noprint;
        select Survival format=best16., StdErr format=best16., CumHaz format=best16.
          into :s trimmed, :se trimmed, :h trimmed
        from synwork.km_estimates where group = "&grp" and Timelist = &tt;
      quit;
      %check(metric=km_s&tt._&grp,  actual=&s,  tol=0.000001)
      %check(metric=km_se&tt._&grp, actual=&se, tol=0.000001)
      %check(metric=na_h&tt._&grp,  actual=&h,  tol=0.000001)
    %end;
  %end;
%mend km_checks;
%km_checks

/* Group comparison tests. */
data synwork.group_tests;
  set synwork._tests;
  keep Test ChiSq DF ProbChiSq;
run;
%export_csv(data=synwork.group_tests, name=group_tests)

proc sql noprint;
  select ChiSq format=best16. into :lr trimmed from synwork._tests
    where compress(upcase(Test), '- ') = 'LOGRANK';
  select ChiSq format=best16. into :wx trimmed from synwork._tests
    where upcase(strip(Test)) = 'WILCOXON';
  select ChiSq format=best16. into :tw trimmed from synwork._tests
    where index(upcase(Test), 'TARONE') > 0;
  select ChiSq format=best16. into :pp trimmed from synwork._tests
    where upcase(strip(Test)) = 'PETO';
quit;
%check(metric=test_logrank_chisq,  actual=&lr, tol=0.0001)
%check(metric=test_wilcoxon_chisq, actual=&wx, tol=0.0001)
%check(metric=test_tarone_chisq,   actual=&tw, tol=0.0001)
%check(metric=test_peto_chisq,     actual=&pp, tol=0.0001)

/* P-values of the same four tests. */
proc sql noprint;
  select ProbChiSq format=best16. into :lrp trimmed from synwork._tests
    where compress(upcase(Test), '- ') = 'LOGRANK';
  select ProbChiSq format=best16. into :wxp trimmed from synwork._tests
    where upcase(strip(Test)) = 'WILCOXON';
  select ProbChiSq format=best16. into :twp trimmed from synwork._tests
    where index(upcase(Test), 'TARONE') > 0;
  select ProbChiSq format=best16. into :ppp trimmed from synwork._tests
    where upcase(strip(Test)) = 'PETO';
  /* Distinct event times in the primary cohort (24 months of deaths). */
  select count(distinct time) into :ndt trimmed from synwork.analysis where event = 1;
quit;
%check(metric=test_logrank_p,  actual=&lrp, tol=0.0001)
%check(metric=test_wilcoxon_p, actual=&wxp, tol=0.0001)
%check(metric=test_tarone_p,   actual=&twp, tol=0.0001)
%check(metric=test_peto_p,     actual=&ppp, tol=0.0001)
%check(metric=distinct_event_times_primary, actual=&ndt)

/*------------------------------------------------------------------------------
  3.3 Survival times.
  The median survival time is not reached. About 3 percent of the cohort dies
  within two years, so the survival curve never falls to 0.5 and PROC LIFETEST
  reports a missing median. The restricted mean survival time (RMST) to day 730
  is the useful summary instead. It is the area under the Kaplan-Meier curve
  from day 0 to day 730, read as "expected days alive out of the first 730".
  On SAS/STAT 15.1 (SAS 9.4M6) or later, the RMST comes from the RMST option
  of PROC LIFETEST. On older releases, a DATA step fallback computes the same
  quantity as the area under the Kaplan-Meier step function from OUTSURV=.
  The fallback is also used if the PROC LIFETEST table is not produced. The
  same reference checks run either way, and the SOURCE column of
  results/tables/rmst.csv records which method produced the numbers. The MEANS
  table that PROC LIFETEST prints by default is restricted to the largest
  event time (day 700), not day 730, so it is not the same quantity.
------------------------------------------------------------------------------*/
title 'Quartiles of survival time (the median is not reached)';
proc print data=synwork._quart noobs;
run;
title;

/* DATA step RMST: area under a step function. On [t(k), t(k+1)) the curve
   equals S(t(k)). */
%macro rmst_datastep(in, out, by=);
  proc sort data=&in out=synwork._rm;
    by &by time;
  run;
  data &out;
    set synwork._rm end=eof;
    %if %length(&by) %then by &by;;
    retain prev_t prev_s;
    %if %length(&by) %then %do;
      if first.&by then do; prev_t = 0; prev_s = 1; rmst730 = 0; end;
    %end;
    %else %do;
      if _n_ = 1 then do; prev_t = 0; prev_s = 1; rmst730 = 0; end;
    %end;
    if time > prev_t then do;
      rmst730 + prev_s * (min(time, &admin_day) - prev_t);
      prev_t = min(time, &admin_day);
    end;
    if not missing(survival) then prev_s = survival;
    %if %length(&by) %then %do;
      if last.&by then do;
    %end;
    %else %do;
      if eof then do;
    %end;
      rmst730 + prev_s * (&admin_day - prev_t);
      output;
    end;
    keep &by rmst730;
  run;
%mend rmst_datastep;

%macro rmst_all;
  %local mpos mlev use_proc;
  %let mpos = %index(&sysvlong, M);
  %if &mpos > 0 %then %let mlev = %substr(&sysvlong, %eval(&mpos + 1), 1);
  %else %let mlev = 0;
  %let use_proc = 0;
  %if %sysevalf(&sysver >= 9.4) and &mlev >= 6 %then %let use_proc = 1;

  %if &use_proc %then %do;
    /* PROC LIFETEST RMST, by CHF group and for the whole cohort. */
    %if %sysfunc(exist(synwork._rmst_p))  %then %do; proc delete data=synwork._rmst_p;  run; %end;
    %if %sysfunc(exist(synwork._rmst_pa)) %then %do; proc delete data=synwork._rmst_pa; run; %end;
    ods output RMST = synwork._rmst_p;
    proc lifetest data=synwork._km rmst(tau=&admin_day) plots=none;
      time time * event(0);
      freq count;
      strata chf;
    run;
    ods output RMST = synwork._rmst_pa;
    proc lifetest data=synwork._km rmst(tau=&admin_day) plots=none;
      time time * event(0);
      freq count;
    run;
    /* Fall back if either table is missing or lacks an Estimate column. */
    %if not %sysfunc(exist(synwork._rmst_p)) or not %sysfunc(exist(synwork._rmst_pa)) %then
      %let use_proc = 0;
    %else %do;
      %local dsid v1 v2 rc;
      %let dsid = %sysfunc(open(synwork._rmst_p));
      %let v1   = %sysfunc(varnum(&dsid, Estimate));
      %let rc   = %sysfunc(close(&dsid));
      %let dsid = %sysfunc(open(synwork._rmst_pa));
      %let v2   = %sysfunc(varnum(&dsid, Estimate));
      %let rc   = %sysfunc(close(&dsid));
      %if &v1 = 0 or &v2 = 0 %then %let use_proc = 0;
    %end;
    %if &use_proc = 0 %then
      %put WARNING: [05] PROC LIFETEST RMST table not found. Using the DATA step fallback.;
  %end;
  %else %put NOTE: [05] SAS release older than 9.4M6. Using the DATA step RMST fallback.;

  %if &use_proc %then %do;
    data synwork.rmst;
      length group $8 source $24;
      set synwork._rmst_p(in=a) synwork._rmst_pa;
      if a then group = ifc(chf = 1, 'chf1', 'chf0');
      else group = 'all';
      rmst730 = Estimate;
      source = 'PROC LIFETEST RMST';
      keep group rmst730 source;
    run;
  %end;
  %else %do;
    /* PLOTS=NONE as well as NOPRINT, because NOPRINT alone still draws the
       default survival plot when ODS Graphics is on. */
    proc lifetest data=synwork._km outsurv=synwork._outs noprint plots=none;
      time time * event(0);
      freq count;
      strata chf;
    run;
    proc lifetest data=synwork._km outsurv=synwork._outs_all noprint plots=none;
      time time * event(0);
      freq count;
    run;
    %rmst_datastep(synwork._outs, synwork._rmst_chf, by=chf)
    %rmst_datastep(synwork._outs_all, synwork._rmst_all)
    data synwork.rmst;
      length group $8 source $24;
      set synwork._rmst_chf(in=a) synwork._rmst_all;
      if a then group = ifc(chf = 1, 'chf1', 'chf0');
      else group = 'all';
      source = 'DATA step (KM area)';
      keep group rmst730 source;
    run;
  %end;
%mend rmst_all;
%rmst_all

proc sql noprint;
  select rmst730 format=best16. into :rm1 trimmed from synwork.rmst where group = 'chf1';
  select rmst730 format=best16. into :rm0 trimmed from synwork.rmst where group = 'chf0';
  select rmst730 format=best16. into :rma trimmed from synwork.rmst where group = 'all';
quit;
%let rmd = %sysevalf(&rm1 - &rm0);

data synwork.rmst;
  set synwork.rmst end=eof;
  label rmst730 = 'Restricted mean survival time to day 730 (days)';
  output;
  if eof then do;
    group = 'diff'; rmst730 = &rmd; output;
  end;
run;

title 'Restricted mean survival time to day 730 (diff = CHF minus no CHF)';
proc print data=synwork.rmst noobs label;
  format rmst730 9.3;
run;
title;
%export_csv(data=synwork.rmst, name=rmst)

%check(metric=rmst730_chf1, actual=&rm1, tol=0.0001)
%check(metric=rmst730_chf0, actual=&rm0, tol=0.0001)
%check(metric=rmst730_all,  actual=&rma, tol=0.0001)
%check(metric=rmst730_diff_chf1_minus_chf0, actual=&rmd, tol=0.0001)

/*------------------------------------------------------------------------------
  Figures. Month scale throughout.
------------------------------------------------------------------------------*/

/* 3.1 The standard PROC LIFETEST plot with confidence limits and an at-risk
   table every three months. Its y axis runs from 0 to 1, so with 3 percent
   mortality the curves sit near the top. */
%fig(km_chf_lifetest)
title 'Kaplan-Meier survival by CHF status';
proc lifetest data=synwork._km plots(only)=survival(cl atrisk=0 to 24 by 3);
  time months * event(0);
  freq count;
  strata chf / notest;
run;
title;
%fig_end

/* Nelson-Aalen estimates on the month scale, for the cumulative hazard plot. */
ods output ProductLimitEstimates = synwork._ple_m;
proc lifetest data=synwork._km nelson plots=none;
  time months * event(0);
  freq count;
  strata chf / notest;
run;

/* README version. Same estimates, with the y axis zoomed so the curves can be
   compared, step-shaped 95 percent bands, and an at-risk table. */
proc lifetest data=synwork._km outsurv=synwork._outs_m noprint plots=none;
  time months * event(0);
  freq count;
  strata chf;
run;

data synwork._grid;
  do tick = 0 to 24 by 3;
    output;
  end;
run;

proc sql;
  create table synwork._atrisk as
  select g.tick, a.chf, sum(a.time * 24 / &admin_day >= g.tick) as atrisk
  from synwork._grid as g, synwork.analysis as a
  group by g.tick, a.chf;
quit;

data synwork._kmplot;
  set synwork._outs_m synwork._atrisk;
  format chf chff.;
run;

%fig(km_chf)
title 'Kaplan-Meier survival by CHF status, primary cohort';
proc sgplot data=synwork._kmplot;
  band x=months lower=sdf_lcl upper=sdf_ucl / type=step group=chf transparency=0.75 name='band';
  step x=months y=survival / group=chf name='km';
  xaxistable atrisk / x=tick class=chf colorgroup=chf title='Number at risk' valueattrs=(size=8);
  xaxis values=(0 to 24 by 3) label='Months since January 1, 2009';
  yaxis min=0.95 max=1 label='Survival probability';
  keylegend 'km' / title='';
run;
title;
%fig_end

/* 3.2 Nelson-Aalen cumulative hazard. Censored rows carry no estimate. */
%fig(cumhaz_chf)
title 'Nelson-Aalen cumulative hazard by CHF status';
proc sgplot data=synwork._ple_m(where=(not missing(CumHaz)));
  step x=months y=CumHaz / group=chf;
  xaxis values=(0 to 24 by 3) label='Months since January 1, 2009';
  yaxis label='Cumulative hazard';
run;
title;
%fig_end

/* 3.5 Kernel-smoothed hazard. BW= is in the units of the TIME variable,
   which here is MONTHS, so BW=3 is three months (about 91 days). On the day
   scale the same smoothing would need BW=91.
   The automatic bandwidth (about 8.4 months)
   spans a third of the 24-month window, so most of the curve was shaped by
   the boundary correction at 0 and 24 months and showed an artificial U
   shape. A fixed 3-month bandwidth keeps the smoothing local. Estimates
   within one bandwidth of either end are still less reliable than those in
   the middle. With deaths only on the first of each month, the smoothing
   spreads each monthly spike over neighboring months. */
%fig(hazard_chf)
title 'Kernel-smoothed hazard by CHF status';
proc lifetest data=synwork._km plots(only)=hazard(bw=3);
  time months * event(0);
  freq count;
  strata chf / notest;
run;
title;
%fig_end
