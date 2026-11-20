/*==============================================================================
  Program   : 02_profile.sas
  Project   : synpuf-survival-sas
  Purpose   : Data quality profiling and cross-year consistency checks that the
              study design depends on. Each count is compared with the Python
              reference value, and the full set is saved as a small CSV.
  Inputs    : SYNWORK.BENE2008, SYNWORK.BENE2009, SYNWORK.BENE2010
  Outputs   : SYNWORK.PROFILE_CHECKS
              results/tables/profile_checks.csv (aggregate counts only)
  Seminar   : data preparation (precedes section 2)
  Notes     : What these checks found in Sample 1, and how the design responds,
                - Each yearly file is the prior year minus that year's deaths.
                  No one enters late and no one is lost, so the loss to
                  follow-up rule is kept in 03 but never fires.
                - Every birth and death date is the first of a month, so time
                  is effectively measured in months and ties are heavy.
                - Death dates never conflict across files and no one appears
                  after death, so "earliest valid death date" is a safeguard
                  rather than a correction.
                - Chronic condition and ESRD flags switch on and off between
                  years, which supports fixing exposure at baseline (2008).
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

/* Stacked view of the three years for frequency tables. A view avoids a
   third copy of the data. */
data synwork._all / view=synwork._all;
  set synwork.bene2008 synwork.bene2009 synwork.bene2010;
run;

/*------------------------------------------------------------------------------
  1. Code values by year. Confirms the coding assumed in the design
     (sex 1/2, race 1/2/3/5, ESRD Y/0, condition flags 1/2).
------------------------------------------------------------------------------*/
title 'Code values by file year';
proc freq data=synwork._all;
  tables (BENE_SEX_IDENT_CD BENE_RACE_CD BENE_ESRD_IND
          SP_ALZHDMTA SP_CHF SP_CHRNKIDN SP_CNCR SP_COPD SP_DEPRESSN
          SP_DIABETES SP_ISCHMCHT SP_OSTEOPRS SP_RA_OA SP_STRKETIA) * YEAR
         / norow nopercent missing;
run;

title 'Coverage months by file year';
proc freq data=synwork._all;
  tables (BENE_HI_CVRAGE_TOT_MONS BENE_SMI_CVRAGE_TOT_MONS BENE_HMO_CVRAGE_TOT_MONS) * YEAR
         / norow nopercent missing;
run;

title 'Annual payment fields by file year';
proc means data=synwork._all n nmiss min p25 median p75 max maxdec=0;
  class YEAR;
  var MEDREIMB_IP BENRES_IP PPPYMT_IP MEDREIMB_OP BENRES_OP PPPYMT_OP
      MEDREIMB_CAR BENRES_CAR PPPYMT_CAR;
run;

/* Day of month for every date. All values should be 1. */
title 'Day of month of birth and death dates';
proc sql;
  select YEAR,
         sum(day(BENE_BIRTH_DT) ne 1) as birth_not_day1,
         sum(day(coalesce(BENE_DEATH_DT, 0)) ne 1) as death_not_day1,
         min(BENE_DEATH_DT) as first_death format=date9.,
         max(BENE_DEATH_DT) as last_death  format=date9.
  from synwork._all
  group by YEAR;
quit;

/* CHF flag transitions from 2008 to 2009 among people present in both years.
   If the flag were cumulative ("ever had CHF"), the 1 to 2 cell would be 0. */
title 'CHF flag in 2008 by CHF flag in 2009';
proc sql;
  create table synwork._chf89 as
  select a.SP_CHF as chf_2008, b.SP_CHF as chf_2009
  from synwork.bene2008 as a inner join synwork.bene2009 as b
    on a.DESYNPUF_ID = b.DESYNPUF_ID;
quit;
proc freq data=synwork._chf89;
  tables chf_2008 * chf_2009 / norow nocol nopercent;
run;
title;

/*------------------------------------------------------------------------------
  2. Consistency checks, one row per check. Every metric name matches a row in
     validation/expected_values.csv.
------------------------------------------------------------------------------*/
/* COALESCE gives a blank death date a harmless stand-in (0 is January 1, 1960,
   day 1 of a month), so date functions never receive a missing argument and
   the log stays free of missing-value notes. */
%macro yr_checks(yr);
  select "rows_&yr" as metric length=40, count(*) as value from synwork.bene&yr
  union all
  select "deaths_recorded_&yr", n(BENE_DEATH_DT) from synwork.bene&yr
  union all
  select "death_outside_file_year_&yr", sum(year(coalesce(BENE_DEATH_DT, mdy(1, 1, &yr))) ne &yr)
    from synwork.bene&yr
  union all
  select "death_before_birth_&yr", sum(not missing(BENE_DEATH_DT) and BENE_DEATH_DT < BENE_BIRTH_DT)
    from synwork.bene&yr
  union all
  select "death_not_first_of_month_&yr", sum(day(coalesce(BENE_DEATH_DT, 0)) ne 1)
    from synwork.bene&yr
  union all
  select "neg_medreimb_rows_&yr", sum(MEDREIMB_IP < 0 or MEDREIMB_OP < 0 or MEDREIMB_CAR < 0)
    from synwork.bene&yr
%mend yr_checks;

proc sql;
  create table synwork.profile_checks as
  %yr_checks(2008) union all
  %yr_checks(2009) union all
  %yr_checks(2010)
  union all
  /* Late entrants: IDs in a later file with no record in the prior file. */
  select 'ids_2009_not_in_2008', sum(missing(a.DESYNPUF_ID))
    from synwork.bene2009 b left join synwork.bene2008 a on b.DESYNPUF_ID = a.DESYNPUF_ID
  union all
  select 'ids_2010_not_in_2009', sum(missing(b.DESYNPUF_ID))
    from synwork.bene2010 c left join synwork.bene2009 b on c.DESYNPUF_ID = b.DESYNPUF_ID
  union all
  /* A person with a 2008 death date who still has a 2009 or 2010 record. */
  select 'appear_after_death_2008',
         (select count(*) from synwork.bene2008 a inner join synwork.bene2009 b
            on a.DESYNPUF_ID = b.DESYNPUF_ID where not missing(a.BENE_DEATH_DT))
       + (select count(*) from synwork.bene2008 a inner join synwork.bene2010 c
            on a.DESYNPUF_ID = c.DESYNPUF_ID where not missing(a.BENE_DEATH_DT))
    from synwork.bene2008(obs=1)
  union all
  select 'appear_after_death_2009', count(*)
    from synwork.bene2009 b inner join synwork.bene2010 c on b.DESYNPUF_ID = c.DESYNPUF_ID
    where not missing(b.BENE_DEATH_DT)
  union all
  select 'birth_dt_differs_2008_2009', sum(a.BENE_BIRTH_DT ne b.BENE_BIRTH_DT)
    from synwork.bene2008 a inner join synwork.bene2009 b on a.DESYNPUF_ID = b.DESYNPUF_ID
  union all
  select 'sex_differs_2008_2009', sum(a.BENE_SEX_IDENT_CD ne b.BENE_SEX_IDENT_CD)
    from synwork.bene2008 a inner join synwork.bene2009 b on a.DESYNPUF_ID = b.DESYNPUF_ID
  union all
  select 'race_differs_2008_2009', sum(a.BENE_RACE_CD ne b.BENE_RACE_CD)
    from synwork.bene2008 a inner join synwork.bene2009 b on a.DESYNPUF_ID = b.DESYNPUF_ID
  union all
  select 'esrd_differs_2008_2009', sum(a.BENE_ESRD_IND ne b.BENE_ESRD_IND)
    from synwork.bene2008 a inner join synwork.bene2009 b on a.DESYNPUF_ID = b.DESYNPUF_ID
  union all
  select 'chf_1_to_2_2008_2009', sum(a.SP_CHF = 1 and b.SP_CHF = 2)
    from synwork.bene2008 a inner join synwork.bene2009 b on a.DESYNPUF_ID = b.DESYNPUF_ID
  union all
  /* People with a death date in two or more yearly files (possible conflicts). */
  select 'ids_with_death_in_2plus_files', count(*) from
    (select DESYNPUF_ID from synwork._all where not missing(BENE_DEATH_DT)
     group by DESYNPUF_ID having count(*) >= 2)
  union all
  /* Implausible ages. Age in completed years on the landmark date. */
  select 'age_min_at_landmark_2008file', min(intck('year', BENE_BIRTH_DT, &landmark_dt, 'c'))
    from synwork.bene2008
  union all
  select 'age_max_at_landmark_2008file', max(intck('year', BENE_BIRTH_DT, &landmark_dt, 'c'))
    from synwork.bene2008
  union all
  select 'birth_not_first_of_month_2008', sum(day(BENE_BIRTH_DT) ne 1) from synwork.bene2008
  ;
quit;

/* Compare every check with the Python reference. CALL EXECUTE with %NRSTR
   queues each %check call to run after this step, one at a time. */
data _null_;
  set synwork.profile_checks;
  call execute(cats('%nrstr(%check)(metric=', metric, ', actual=', value, ')'));
run;

title 'Data quality and consistency checks';
proc print data=synwork.profile_checks noobs;
run;
title;

%export_csv(data=synwork.profile_checks, name=profile_checks)
