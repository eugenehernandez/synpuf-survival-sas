/*==============================================================================
  Program   : 04_explore.sas
  Project   : synpuf-survival-sas
  Purpose   : Describe the primary cohort before any survival modeling.
              Table 1 by CHF status, the distribution of follow-up time among
              people who died, and correlations among the covariates.
  Inputs    : SYNWORK.ANALYSIS
  Outputs   : SYNWORK.TABLE1,        results/tables/table1.csv
              SYNWORK.CORRELATIONS,  results/tables/correlations.csv
              results/figures/followup_deaths_hist.png
  Seminar   : section 2, data exploration
  Notes     : People with CHF in 2008 differ from people without it (more
              kidney disease, ischemic heart disease, and higher reimbursement)
              even though death in this release is unrelated to these
              covariates. That pattern is what makes the adjusted models in
              06_cox.sas worth fitting despite a null crude association.
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
  1. Table 1 by CHF status.
     One PROC MEANS pass computes every statistic for the overall cohort
     (_TYPE_ = 0) and for each CHF group (_TYPE_ = 1). For a 0/1 variable the
     SUM is the count and the MEAN is the proportion. Percentiles use the SAS
     default definition (QNTLDEF=5), which the Python check reproduces.
------------------------------------------------------------------------------*/
data synwork._t1in / view=synwork._t1in;
  set synwork.analysis;
  race_white = (race = 1);
  race_black = (race = 2);
  race_other = (race = 3);
  race_hisp  = (race = 5);
run;

%let t1bin = female race_white race_black race_other race_hisp esrd alzh ckd cancer
             copd depress diab ihd osteo raoa stroke event;

proc means data=synwork._t1in noprint;
  class chf;
  types () chf;
  var age logreimb reimb08 time &t1bin;
  output out=synwork._t1 n(age)=n_people mean= std= median= p25= p75= sum= / autoname;
run;

/* One output row per characteristic and group. BYTE(37) is the percent sign,
   written that way so the macro processor never sees a bare % in the text. */
%macro t1bin(var, label, ord);
  ord = &ord; characteristic = "&label"; statistic = 'n (' || byte(37) || ')';
  value = strip(put(&var._Sum, comma9.)) || ' (' || strip(put(100 * &var._Mean, 5.1)) || ')';
  output;
%mend t1bin;

data synwork._t1long;
  set synwork._t1;
  length grp $8 characteristic $60 statistic $16 value $40;
  if _type_ = 0 then grp = 'Overall';
  else if chf = 1 then grp = 'CHF';
  else grp = 'NoCHF';

  ord = 1; characteristic = 'Beneficiaries'; statistic = 'N';
    value = strip(put(n_people, comma9.)); output;
  ord = 2; characteristic = 'Age on January 1, 2009 (years)'; statistic = 'Mean (SD)';
    value = strip(put(age_Mean, 6.1)) || ' (' || strip(put(age_StdDev, 5.1)) || ')'; output;
  %t1bin(female,     Female,                                   3)
  %t1bin(race_white, Race White,                               4)
  %t1bin(race_black, Race Black,                               5)
  %t1bin(race_other, Race Other,                               6)
  %t1bin(race_hisp,  Race Hispanic,                            7)
  %t1bin(esrd,       End-stage renal disease,                  8)
  %t1bin(alzh,       Alzheimer or related dementia,            9)
  %t1bin(ckd,        Chronic kidney disease,                  10)
  %t1bin(cancer,     Cancer,                                  11)
  %t1bin(copd,       COPD,                                    12)
  %t1bin(depress,    Depression,                              13)
  %t1bin(diab,       Diabetes,                                14)
  %t1bin(ihd,        Ischemic heart disease,                  15)
  %t1bin(osteo,      Osteoporosis,                            16)
  %t1bin(raoa,       Rheumatoid arthritis or osteoarthritis,  17)
  %t1bin(stroke,     Stroke or TIA,                           18)
  ord = 19; characteristic = 'Total 2008 Medicare reimbursement (USD)'; statistic = 'Median [Q1, Q3]';
    value = strip(put(reimb08_Median, comma9.)) || ' [' || strip(put(reimb08_P25, comma9.))
            || ', ' || strip(put(reimb08_P75, comma9.)) || ']'; output;
  ord = 20; characteristic = 'log(1 + 2008 reimbursement)'; statistic = 'Mean (SD)';
    value = strip(put(logreimb_Mean, 6.2)) || ' (' || strip(put(logreimb_StdDev, 5.2)) || ')'; output;
  %t1bin(event,      Died during follow-up,                   21)
  ord = 22; characteristic = 'Follow-up (days)'; statistic = 'Mean (SD)';
    value = strip(put(time_Mean, 6.1)) || ' (' || strip(put(time_StdDev, 5.1)) || ')'; output;
  keep ord grp characteristic statistic value;
run;

proc sort data=synwork._t1long;
  by ord characteristic statistic;
run;

proc transpose data=synwork._t1long out=synwork._t1wide(drop=_name_);
  by ord characteristic statistic;
  id grp;
  var value;
run;

data synwork.table1;
  retain characteristic statistic NoCHF CHF Overall;
  set synwork._t1wide;
  label characteristic = 'Characteristic' statistic = 'Statistic'
        NoCHF = 'No CHF' CHF = 'CHF' Overall = 'Overall';
  drop ord;
run;

title 'Table 1. Primary cohort by congestive heart failure (CHF) status in 2008';
proc print data=synwork.table1 noobs label;
run;
title;
%export_csv(data=synwork.table1, name=table1)

/* Reference checks on the numeric values behind Table 1. */
%macro t1check(where, tag);
  proc sql noprint;
    select n_people format=best16., age_Mean format=best16., age_StdDev format=best16.,
           female_Sum format=best16., race_white_Sum format=best16., esrd_Sum format=best16.,
           diab_Sum format=best16., ckd_Sum format=best16., reimb08_Median format=best16.,
           reimb08_P25 format=best16., reimb08_P75 format=best16., logreimb_Mean format=best16.,
           event_Sum format=best16.
      into :c_n trimmed, :c_am trimmed, :c_asd trimmed, :c_fem trimmed, :c_wh trimmed,
           :c_esrd trimmed, :c_diab trimmed, :c_ckd trimmed, :c_med trimmed, :c_p25 trimmed,
           :c_p75 trimmed, :c_lr trimmed, :c_ev trimmed
    from synwork._t1 where &where;
  quit;
  %check(metric=t1_n_&tag,             actual=&c_n)
  %check(metric=t1_age_mean_&tag,      actual=&c_am,  tol=0.0001)
  %check(metric=t1_age_sd_&tag,        actual=&c_asd, tol=0.0001)
  %check(metric=t1_female_n_&tag,      actual=&c_fem)
  %check(metric=t1_white_n_&tag,       actual=&c_wh)
  %check(metric=t1_esrd_n_&tag,        actual=&c_esrd)
  %check(metric=t1_diab_n_&tag,        actual=&c_diab)
  %check(metric=t1_ckd_n_&tag,         actual=&c_ckd)
  %check(metric=t1_reimb_median_&tag,  actual=&c_med)
  %check(metric=t1_reimb_p25_&tag,     actual=&c_p25)
  %check(metric=t1_reimb_p75_&tag,     actual=&c_p75)
  %check(metric=t1_logreimb_mean_&tag, actual=&c_lr,  tol=0.0001)
  %check(metric=t1_deaths_&tag,        actual=&c_ev)
%mend t1check;

%t1check(%str(_type_ = 1 and chf = 1), chf1)
%t1check(%str(_type_ = 1 and chf = 0), chf0)
%t1check(%str(_type_ = 0),             all)

/*------------------------------------------------------------------------------
  2. Follow-up time among people who died.
     Deaths are recorded on the first of a month, so the histogram uses
     one-month bins on a month scale (MONTHS = TIME x 24 / 730, so day 365 is
     month 12 and day 730 is month 24). Summary statistics stay in days.
------------------------------------------------------------------------------*/
data synwork._deaths / view=synwork._deaths;
  set synwork.analysis(where=(event = 1));
  months = time * 24 / &admin_day;
  label time = 'Days from December 31, 2008 to death'
        months = 'Months since January 1, 2009';
run;

/* Distribution summaries in days (shown in the results tab). */
title 'Follow-up time among beneficiaries who died (primary cohort)';
proc univariate data=synwork._deaths;
  var time;
run;

/* Histogram on the month scale (written to results/figures). HISTOGRAM can
   only plot variables named in the VAR statement. */
%fig(followup_deaths_hist)
proc univariate data=synwork._deaths;
  var months;
  histogram months / endpoints=0 to 24 by 1 odstitle='Deaths by month of follow-up';
run;
%fig_end
title;

proc means data=synwork._deaths noprint;
  var time;
  output out=synwork._evt mean=m std=s median=med;
run;
proc sql noprint;
  select m format=best16., s format=best16., med format=best16.
    into :ev_m trimmed, :ev_s trimmed, :ev_med trimmed
  from synwork._evt;
quit;
%check(metric=evtime_mean,   actual=&ev_m, tol=0.0001)
%check(metric=evtime_sd,     actual=&ev_s, tol=0.0001)
%check(metric=evtime_median, actual=&ev_med)

/*------------------------------------------------------------------------------
  3. Correlations among the covariates.
     Pearson correlations on 0/1 flags are phi coefficients. They show which
     conditions travel together, which matters for interpreting adjusted
     hazard ratios. Graphics are off because a scatter matrix of 0/1 data
     carries no information.
------------------------------------------------------------------------------*/
ods output PearsonCorr = synwork._corr;
title 'Pearson correlations among covariates (primary cohort)';
proc corr data=synwork.analysis pearson nosimple noprob plots=none;
  var age logreimb female esrd chf alzh ckd cancer copd depress diab ihd osteo raoa stroke;
run;
title;

data synwork.correlations;
  set synwork._corr;
  keep Variable age logreimb female esrd chf alzh ckd cancer copd depress diab ihd osteo raoa stroke;
  format age logreimb female esrd chf alzh ckd cancer copd depress diab ihd osteo raoa stroke 6.3;
run;
%export_csv(data=synwork.correlations, name=correlations)

proc sql noprint;
  select logreimb format=best16. into :r_al trimmed from synwork._corr where upcase(Variable) = 'AGE';
  select ckd format=best16., ihd format=best16. into :r_ck trimmed, :r_ih trimmed
    from synwork._corr where upcase(Variable) = 'CHF';
quit;
%check(metric=corr_age_logreimb, actual=&r_al, tol=0.0001)
%check(metric=corr_chf_ckd,      actual=&r_ck, tol=0.0001)
%check(metric=corr_chf_ihd,      actual=&r_ih, tol=0.0001)
