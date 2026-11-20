/*==============================================================================
  Program   : 03_build_cohort.sas
  Project   : synpuf-survival-sas
  Purpose   : Build the landmark cohort, the attrition table, and the analysis
              datasets for the primary analysis and the cohort-based
              sensitivity analyses.
  Inputs    : SYNWORK.BENE2008, SYNWORK.BENE2009, SYNWORK.BENE2010
  Outputs   : SYNWORK.LANDMARK_ALL   every 2008 beneficiary with derived fields
              SYNWORK.COHORT_ALLAGES alive and enrolled on January 1, 2009, any age
              SYNWORK.COHORT_65      as above, age 65 or older
              SYNWORK.ANALYSIS       primary cohort (65 or older, fee-for-service 2008)
              SYNWORK.COHORT_COV     primary cohort minus people who switch coverage
                                     in 2009 (coverage-censoring sensitivity analysis)
              SYNWORK.ATTRITION, results/tables/attrition.csv
              SYNWORK.ATTRITION_COV, results/tables/attrition_coverage.csv
              SYNWORK.COHORT_SUMMARY, results/tables/cohort_summary.csv
  Seminar   : data preparation (precedes section 2)

  Design (see docs/methods.md for the reasoning)
    Baseline     The 2008 file. Condition flags summarize the whole calendar
                 year, so they are only "known" at its end. Follow-up therefore
                 starts at a landmark on January 1, 2009, which avoids using
                 information from after time zero.
    Eligibility  Alive on January 1, 2009, present in the 2009 file, age 65 or
                 older on January 1, 2009, and in 2008 enrolled for 12 months of
                 Part A and 12 months of Part B with no Medicare Advantage months.
                 The last rule makes the claims-derived 2008 flags complete.
    Exposure     Congestive heart failure (SP_CHF = 1 in 2008).
    Outcome      All-cause death, earliest valid death date across 2008 to 2010.
    Time         Days since December 31, 2008. Day 1 is January 1, 2009 and
                 administrative censoring falls on day 730 (December 31, 2010).
                 Death dates are recorded to the month (always the first of the
                 month), so event times take at most 24 distinct values.
    Censoring    Administrative at day 730. Loss to follow-up (in the 2009 file,
                 absent from 2010, no death) is censored at day 365. That rule
                 never fires in Sample 1 but stays in the code and the table.
    Coverage     Sensitivity variables TIME_COV and EVENT_COV. A follow-up year
                 "switches" when it has any Medicare Advantage months or fewer
                 than 12 Part A months. Follow-up is censored at the start of the
                 first switching year (day 0 for 2009, day 365 for 2010), unless a
                 death is recorded in that year, in which case the death counts.
                 People who would be censored at day 0 contribute no follow-up,
                 so they are excluded from that analysis (SYNWORK.COHORT_COV)
                 rather than kept as zero-length records. This keeps the result
                 independent of how each procedure treats a time of zero.
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
  1. One row per 2008 beneficiary with the follow-up years attached.
     LEFT JOINs on a unique key, so there are no repeated BY values and no row
     multiplication. IN2009 and IN2010 record presence in each later file.
------------------------------------------------------------------------------*/
proc sql;
  create table synwork._base as
  select a.*,
         not missing(b.DESYNPUF_ID)    as in2009,
         not missing(c.DESYNPUF_ID)    as in2010,
         b.BENE_DEATH_DT               as death09 format=date9.,
         c.BENE_DEATH_DT               as death10 format=date9.,
         b.BENE_HI_CVRAGE_TOT_MONS     as hi09,
         b.BENE_HMO_CVRAGE_TOT_MONS    as hmo09,
         c.BENE_HI_CVRAGE_TOT_MONS     as hi10,
         c.BENE_HMO_CVRAGE_TOT_MONS    as hmo10
  from synwork.bene2008 as a
       left join synwork.bene2009 as b on a.DESYNPUF_ID = b.DESYNPUF_ID
       left join synwork.bene2010 as c on a.DESYNPUF_ID = c.DESYNPUF_ID
  order by a.DESYNPUF_ID;
quit;

/*------------------------------------------------------------------------------
  2. Derived analysis variables.
     Every IF that touches a possibly missing date is guarded, so the log has
     no "missing values were generated" notes.
------------------------------------------------------------------------------*/
data synwork.landmark_all;
  set synwork._base;

  /* Earliest valid death date across the three files. MIN ignores missing
     arguments, but returns missing with a log note when all are missing, so
     it runs only when at least one death date exists. */
  if n(BENE_DEATH_DT, death09, death10) > 0 then
    death_dt = min(BENE_DEATH_DT, death09, death10);
  else death_dt = .;
  if not missing(death_dt) then death_yr = year(death_dt);
  else death_yr = 0;

  /* Age. AGE_INT (completed years) defines eligibility. AGE (exact years) is
     the continuous covariate. Birth dates are always the first of a month. */
  age_int = intck('year', BENE_BIRTH_DT, &landmark_dt, 'c');
  age     = (&landmark_dt - BENE_BIRTH_DT) / 365.25;
  under65 = (age_int < 65);

  /* Demographics and ESRD, recoded to 0/1. */
  female = (BENE_SEX_IDENT_CD = 2);
  race   = BENE_RACE_CD;
  esrd   = (BENE_ESRD_IND = 'Y');

  /* Chronic condition flags, recoded from 1 = yes, 2 = no to 1/0. */
  array src{11} SP_ALZHDMTA SP_CHF SP_CHRNKIDN SP_CNCR SP_COPD SP_DEPRESSN
                SP_DIABETES SP_ISCHMCHT SP_OSTEOPRS SP_RA_OA SP_STRKETIA;
  array cc{11}  alzh chf ckd cancer copd depress diab ihd osteo raoa stroke;
  do _i = 1 to 11;
    cc{_i} = (src{_i} = 1);
  end;
  drop _i;

  /* Severity proxy. A few reimbursement amounts are negative (claim
     adjustments), so each component is floored at 0 before summing. */
  reimb08  = max(0, MEDREIMB_IP) + max(0, MEDREIMB_OP) + max(0, MEDREIMB_CAR);
  logreimb = log(1 + reimb08);

  /* Fee-for-service in 2008: complete Part A and B, no Medicare Advantage. */
  ffs08 = (BENE_HI_CVRAGE_TOT_MONS = 12 and BENE_SMI_CVRAGE_TOT_MONS = 12
           and BENE_HMO_CVRAGE_TOT_MONS = 0);

  /* Eligibility flags, in attrition order. */
  alive_landmark = (missing(death_dt) or death_dt > &origin_dt);
  elig_step3 = alive_landmark and in2009;
  elig_step4 = elig_step3 and not under65;
  elig_step5 = elig_step4 and ffs08;

  /* Primary outcome and follow-up time. */
  event = (&landmark_dt <= death_dt <= '31DEC2010'd);
  ltfu  = (in2009 and not in2010 and missing(death_dt));
  if event then time = death_dt - &origin_dt;
  else if ltfu then time = 365;
  else time = &admin_day;

  /* Coverage-based censoring (sensitivity analysis 2). */
  sw09 = (in2009 and (hmo09 > 0 or hi09 < 12));
  sw10 = (in2010 and (hmo10 > 0 or hi10 < 12));
  time_cov  = time;
  event_cov = event;
  if sw09 and death_yr ne 2009 then do;
    time_cov = 0;  event_cov = 0;
  end;
  else if not sw09 and sw10 and death_yr not in (2009, 2010) then do;
    time_cov = 365; event_cov = 0;
  end;

  format death_dt date9. female femalef. race racef. esrd chf alzh ckd cancer copd
         depress diab ihd osteo raoa stroke under65 ffs08 ynf. event event_cov eventf.;
  label death_dt  = 'Earliest recorded death date'
        age       = 'Age on January 1, 2009 (years)'
        age_int   = 'Age on January 1, 2009 (completed years)'
        under65   = 'Younger than 65 on January 1, 2009'
        female    = 'Female'
        race      = 'Race'
        esrd      = 'End-stage renal disease (2008)'
        alzh      = 'Alzheimer or related dementia (2008)'
        chf       = 'Congestive heart failure (2008)'
        ckd       = 'Chronic kidney disease (2008)'
        cancer    = 'Cancer (2008)'
        copd      = 'COPD (2008)'
        depress   = 'Depression (2008)'
        diab      = 'Diabetes (2008)'
        ihd       = 'Ischemic heart disease (2008)'
        osteo     = 'Osteoporosis (2008)'
        raoa      = 'Rheumatoid arthritis or osteoarthritis (2008)'
        stroke    = 'Stroke or TIA (2008)'
        reimb08   = 'Total 2008 Medicare reimbursement, IP + OP + carrier'
        logreimb  = 'log(1 + total 2008 Medicare reimbursement)'
        ffs08     = 'Fee-for-service all 12 months of 2008'
        event     = 'Died during follow-up'
        time      = 'Days from December 31, 2008 to death or censoring'
        ltfu      = 'Censored for loss to follow-up'
        time_cov  = 'Days to death or censoring, coverage-based censoring'
        event_cov = 'Died during follow-up, coverage-based censoring';
run;

/*------------------------------------------------------------------------------
  3. Analysis datasets. Only analysis variables are kept, so later programs
     cannot accidentally use the raw 2009 or 2010 fields.
------------------------------------------------------------------------------*/
%let keepvars = DESYNPUF_ID time event time_cov event_cov ltfu chf age age_int under65
                female race esrd alzh ckd cancer copd depress diab ihd osteo raoa
                stroke reimb08 logreimb ffs08 SP_STATE_CODE;

data synwork.cohort_allages(keep=&keepvars)
     synwork.cohort_65(keep=&keepvars)
     synwork.analysis(keep=&keepvars)
     synwork.cohort_cov(keep=&keepvars);
  set synwork.landmark_all;
  if elig_step3 then output synwork.cohort_allages;
  if elig_step4 then output synwork.cohort_65;
  if elig_step5 then output synwork.analysis;
  /* Coverage-censoring cohort. Drop anyone whose follow-up would end at day 0. */
  if elig_step5 and not (time_cov = 0 and event_cov = 0) then output synwork.cohort_cov;
run;

/*------------------------------------------------------------------------------
  4. Attrition table (CONSORT-style flow).
------------------------------------------------------------------------------*/
proc sql noprint;
  select count(*) format=best16., sum(alive_landmark) format=best16., sum(elig_step3) format=best16.,
         sum(elig_step4) format=best16., sum(elig_step5) format=best16.
    into :n1 trimmed, :n2 trimmed, :n3 trimmed, :n4 trimmed, :n5 trimmed
  from synwork.landmark_all;
quit;

data synwork.attrition;
  length step 8 description $90 n_remaining n_excluded 8;
  step = 1; description = 'Beneficiaries in the 2008 Beneficiary Summary File';
    n_remaining = &n1; n_excluded = .; output;
  step = 2; description = 'Excluded: death recorded on or before December 31, 2008';
    n_remaining = &n2; n_excluded = &n1 - &n2; output;
  step = 3; description = 'Excluded: alive but absent from the 2009 file';
    n_remaining = &n3; n_excluded = &n2 - &n3; output;
  step = 4; description = 'Excluded: younger than 65 on January 1, 2009';
    n_remaining = &n4; n_excluded = &n3 - &n4; output;
  step = 5; description = 'Excluded: not 12 months Part A and Part B with no Medicare Advantage in 2008';
    n_remaining = &n5; n_excluded = &n4 - &n5; output;
run;

%check(metric=attrition_step1_n, actual=&n1)
%check(metric=attrition_step2_n, actual=&n2)
%check(metric=attrition_step3_n, actual=&n3)
%check(metric=attrition_step4_n, actual=&n4)
%check(metric=attrition_step5_n, actual=&n5)

title 'Attrition from the 2008 file to the primary analysis cohort';
proc print data=synwork.attrition noobs label;
  format n_remaining n_excluded comma9.;
run;
title;
%export_csv(data=synwork.attrition, name=attrition)

/* Sensitivity attrition for the coverage-censoring analysis. Steps 1 to 5 are
   the primary flow. Step 6 removes people whose 2009 record shows Medicare
   Advantage months or fewer than 12 Part A months and who did not die in 2009,
   because their follow-up would be censored at day 0. */
%nobs(synwork.cohort_cov, n6)
data synwork.attrition_cov;
  set synwork.attrition end=last;
  output;
  if last then do;
    step = 6;
    description = 'Excluded: switched to Medicare Advantage or partial Part A in 2009, no 2009 death';
    n_excluded  = n_remaining - &n6;
    n_remaining = &n6;
    output;
  end;
run;
%check(metric=attrition_cov_step6_n, actual=&n6)

title 'Attrition for the coverage-censoring sensitivity analysis';
proc print data=synwork.attrition_cov noobs label;
  format n_remaining n_excluded comma9.;
run;
title;
%export_csv(data=synwork.attrition_cov, name=attrition_coverage)

/*------------------------------------------------------------------------------
  5. Cohort summary for each analysis population: size, deaths, CHF, loss to
     follow-up, and person-years. One row per cohort, aggregate only.
------------------------------------------------------------------------------*/
%macro summarize(ds, label, tag, time=time, event=event);
  proc sql noprint;
    create table synwork._sum_&tag as
    select "&label" as cohort length=40,
           count(*)                 as n,
           sum(&event)              as deaths,
           sum(chf)                 as chf_n,
           sum(chf * &event)        as chf_deaths,
           sum(ltfu)                as ltfu_censored,
           sum(&time)               as sum_days,
           sum(&time) / 365.25      as person_years format=comma12.1
    from &ds;
    select n format=best16., deaths format=best16., chf_n format=best16.,
           chf_deaths format=best16., sum_days format=best16.
      into :&tag._n trimmed, :&tag._ev trimmed, :&tag._chfn trimmed,
           :&tag._chfev trimmed, :&tag._days trimmed
    from synwork._sum_&tag;
  quit;
  %check(metric=n_&tag,          actual=&&&tag._n)
  %check(metric=events_&tag,     actual=&&&tag._ev)
  %check(metric=chf_n_&tag,      actual=&&&tag._chfn)
  %check(metric=chf_events_&tag, actual=&&&tag._chfev)
  %check(metric=sum_time_&tag,   actual=&&&tag._days)
%mend summarize;

%summarize(synwork.analysis,       %str(Primary: 65+, fee-for-service 2008), primary)
%summarize(synwork.cohort_65,      %str(Sensitivity: all 65+),               full65)
%summarize(synwork.cohort_allages, %str(Sensitivity: all ages),              all_ages)
%summarize(synwork.cohort_cov,     %str(Sensitivity: coverage censoring),    cov,
           time=time_cov, event=event_cov)

data synwork.cohort_summary;
  set synwork._sum_primary synwork._sum_full65 synwork._sum_all_ages synwork._sum_cov;
run;

title 'Analysis populations';
proc print data=synwork.cohort_summary noobs;
run;
title;
%export_csv(data=synwork.cohort_summary, name=cohort_summary)

/*------------------------------------------------------------------------------
  6. Remaining reference checks on the primary cohort and the coverage-based
     censoring variables.
------------------------------------------------------------------------------*/
proc sql noprint;
  /* FORMAT=BEST16. keeps full precision in the macro variables (the INTO
     default is BEST8.). */
  select sum(ltfu) format=best16., mean(age) format=best16., mean(logreimb) format=best16.,
         sum(female) format=best16., sum(esrd) format=best16.,
         sum(event = 1 and time <= 365) format=best16.,
         sum(event_cov = 0 and time_cov = 0) format=best16.,
         sum(event_cov = 0 and time_cov = 365) format=best16.,
         sum(event_cov) format=best16., sum(time_cov) format=best16.
    into :p_ltfu trimmed, :p_age trimmed, :p_logr trimmed, :p_fem trimmed, :p_esrd trimmed,
         :p_ev365 trimmed, :p_cov0 trimmed, :p_cov365 trimmed, :p_covev trimmed, :p_covdays trimmed
  from synwork.analysis;
  select sum(under65), sum(under65 * event)
    into :u65_n trimmed, :u65_ev trimmed
  from synwork.cohort_allages;
quit;

%check(metric=ltfu_censored_2009_primary,         actual=&p_ltfu)
%check(metric=mean_age_primary,                   actual=&p_age,  tol=0.0001)
%check(metric=mean_logreimb_primary,              actual=&p_logr, tol=0.0001)
%check(metric=female_n_primary,                   actual=&p_fem)
%check(metric=esrd_n_primary,                     actual=&p_esrd)
%check(metric=events_at_or_before_day365_primary, actual=&p_ev365)
%check(metric=cov_censored_day0_primary,          actual=&p_cov0)
%check(metric=cov_censored_day365_primary,        actual=&p_cov365)
%check(metric=cov_events_primary,                 actual=&p_covev)
%check(metric=cov_sum_time_primary,               actual=&p_covdays)

/* The deaths that leave with the people excluded at day 0 (expected 107),
   and the deaths among the remaining 55,075 that coverage censoring turns
   into censored records (expected 0). */
proc sql noprint;
  select sum(event) format=best16. into :p_covexd trimmed
    from synwork.analysis where time_cov = 0 and event_cov = 0;
  select sum(event - event_cov) format=best16. into :p_covconv trimmed
    from synwork.cohort_cov;
quit;
%check(metric=cov_excluded_primary_deaths,           actual=&p_covexd)
%check(metric=cov_kept_deaths_converted_to_censored, actual=&p_covconv)
%check(metric=n_under65_alive09,                  actual=&u65_n)
%check(metric=events_under65_alive09,             actual=&u65_ev)

/* Distinct event times in the primary cohort. Expect 24 (one per month). */
title 'Deaths by follow-up day in the primary cohort';
proc freq data=synwork.analysis(where=(event = 1));
  tables time / nocum;
run;
title;
