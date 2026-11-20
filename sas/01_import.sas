/*==============================================================================
  Program   : 01_import.sas
  Project   : synpuf-survival-sas
  Purpose   : Read the 2008, 2009, and 2010 Beneficiary Summary files with an
              explicit DATA step (fixed lengths, informats, and labels) instead
              of PROC IMPORT, so variable types never depend on a guess from the
              first rows of the file.
  Inputs    : &file_2008, &file_2009, &file_2010 (defined in 00_config.sas)
  Outputs   : SYNWORK.BENE2008, SYNWORK.BENE2009, SYNWORK.BENE2010
              (one row per beneficiary, plus YEAR)
  Seminar   : data preparation (precedes section 2)
  Notes     : The files use Windows line endings (CRLF), a quoted header row,
              and unquoted data rows. TERMSTR=CRLF is stated explicitly so the
              program reads the same way on Windows and Unix SAS servers.
              Dates are YYYYMMDD and every date in this release falls on the
              first of the month (see 02_profile.sas).
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

%macro import_year(yr);

  data synwork.bene&yr;
    infile "&&file_&yr" dsd dlm=',' firstobs=2 truncover termstr=crlf lrecl=1024;

    /* Lengths come from the profile of all three files (longest observed value). */
    length DESYNPUF_ID $16
           BENE_BIRTH_DT BENE_DEATH_DT 8
           BENE_SEX_IDENT_CD BENE_RACE_CD 8
           BENE_ESRD_IND $1
           SP_STATE_CODE $2 BENE_COUNTY_CD $3
           BENE_HI_CVRAGE_TOT_MONS BENE_SMI_CVRAGE_TOT_MONS
           BENE_HMO_CVRAGE_TOT_MONS PLAN_CVRG_MOS_NUM 8
           SP_ALZHDMTA SP_CHF SP_CHRNKIDN SP_CNCR SP_COPD SP_DEPRESSN
           SP_DIABETES SP_ISCHMCHT SP_OSTEOPRS SP_RA_OA SP_STRKETIA 8
           MEDREIMB_IP BENRES_IP PPPYMT_IP
           MEDREIMB_OP BENRES_OP PPPYMT_OP
           MEDREIMB_CAR BENRES_CAR PPPYMT_CAR 8
           YEAR 8;

    /* A blank death date under DSD is read as a missing value without an
       invalid-data note. */
    informat BENE_BIRTH_DT BENE_DEATH_DT yymmdd8.;
    format   BENE_BIRTH_DT BENE_DEATH_DT date9.;

    input DESYNPUF_ID $
          BENE_BIRTH_DT BENE_DEATH_DT
          BENE_SEX_IDENT_CD BENE_RACE_CD
          BENE_ESRD_IND $
          SP_STATE_CODE $ BENE_COUNTY_CD $
          BENE_HI_CVRAGE_TOT_MONS BENE_SMI_CVRAGE_TOT_MONS
          BENE_HMO_CVRAGE_TOT_MONS PLAN_CVRG_MOS_NUM
          SP_ALZHDMTA SP_CHF SP_CHRNKIDN SP_CNCR SP_COPD SP_DEPRESSN
          SP_DIABETES SP_ISCHMCHT SP_OSTEOPRS SP_RA_OA SP_STRKETIA
          MEDREIMB_IP BENRES_IP PPPYMT_IP
          MEDREIMB_OP BENRES_OP PPPYMT_OP
          MEDREIMB_CAR BENRES_CAR PPPYMT_CAR;

    YEAR = &yr;

    label DESYNPUF_ID              = 'Synthetic beneficiary ID'
          BENE_BIRTH_DT            = 'Date of birth'
          BENE_DEATH_DT            = 'Date of death'
          BENE_SEX_IDENT_CD        = 'Sex (1=male, 2=female)'
          BENE_RACE_CD             = 'Race (1=White, 2=Black, 3=Other, 5=Hispanic)'
          BENE_ESRD_IND            = 'End-stage renal disease (Y/0)'
          SP_STATE_CODE            = 'State code'
          BENE_COUNTY_CD           = 'County code'
          BENE_HI_CVRAGE_TOT_MONS  = 'Part A coverage months'
          BENE_SMI_CVRAGE_TOT_MONS = 'Part B coverage months'
          BENE_HMO_CVRAGE_TOT_MONS = 'Medicare Advantage (HMO) coverage months'
          PLAN_CVRG_MOS_NUM        = 'Part D coverage months'
          SP_ALZHDMTA = 'Alzheimer or related dementia (1=yes, 2=no)'
          SP_CHF      = 'Congestive heart failure (1=yes, 2=no)'
          SP_CHRNKIDN = 'Chronic kidney disease (1=yes, 2=no)'
          SP_CNCR     = 'Cancer (1=yes, 2=no)'
          SP_COPD     = 'COPD (1=yes, 2=no)'
          SP_DEPRESSN = 'Depression (1=yes, 2=no)'
          SP_DIABETES = 'Diabetes (1=yes, 2=no)'
          SP_ISCHMCHT = 'Ischemic heart disease (1=yes, 2=no)'
          SP_OSTEOPRS = 'Osteoporosis (1=yes, 2=no)'
          SP_RA_OA    = 'Rheumatoid arthritis or osteoarthritis (1=yes, 2=no)'
          SP_STRKETIA = 'Stroke or TIA (1=yes, 2=no)'
          MEDREIMB_IP  = 'Inpatient Medicare reimbursement'
          BENRES_IP    = 'Inpatient beneficiary responsibility'
          PPPYMT_IP    = 'Inpatient primary payer amount'
          MEDREIMB_OP  = 'Outpatient Medicare reimbursement'
          BENRES_OP    = 'Outpatient beneficiary responsibility'
          PPPYMT_OP    = 'Outpatient primary payer amount'
          MEDREIMB_CAR = 'Carrier Medicare reimbursement'
          BENRES_CAR   = 'Carrier beneficiary responsibility'
          PPPYMT_CAR   = 'Carrier primary payer amount'
          YEAR         = 'File year';
  run;

  /* Sorted by ID so later steps can rely on the order. NODUPKEY with a DUPOUT
     data set makes any duplicate ID visible instead of silently dropped. */
  proc sort data=synwork.bene&yr nodupkey dupout=synwork._dup&yr;
    by DESYNPUF_ID;
  run;

  %nobs(synwork.bene&yr, n_&yr)
  %nobs(synwork._dup&yr, dup_&yr)
  %check(metric=rows_&yr, actual=&&n_&yr)
  %if &&dup_&yr > 0 %then %put WARNING: [import] &&dup_&yr duplicate IDs in &yr;
  %else %put NOTE: [import] No duplicate IDs in &yr;

%mend import_year;

%import_year(2008)
%import_year(2009)
%import_year(2010)

/* Structure check that lands in the results window for review. */
title 'Imported Beneficiary Summary files';
proc sql;
  select 2008 as year, count(*) as rows, count(distinct DESYNPUF_ID) as ids,
         nmiss(BENE_BIRTH_DT) as missing_birth, n(BENE_DEATH_DT) as deaths
    from synwork.bene2008
  union all
  select 2009, count(*), count(distinct DESYNPUF_ID), nmiss(BENE_BIRTH_DT), n(BENE_DEATH_DT)
    from synwork.bene2009
  union all
  select 2010, count(*), count(distinct DESYNPUF_ID), nmiss(BENE_BIRTH_DT), n(BENE_DEATH_DT)
    from synwork.bene2010;
quit;
title;
