/*==============================================================================
  Program   : 00_config.sas
  Project   : synpuf-survival-sas
  Purpose   : Central configuration. Defines the data location, the three input
              files, the repository folders, the WORK-independent library, global
              options, formats, ODS Graphics settings, the random seed, and small
              utility macros used by every other program.
  Inputs    : none
  Outputs   : macro variables, library SYNWORK, format catalog SYNWORK.FORMATS
  Seminar   : setup (no seminar section)
  Notes     : Every program in sas/ includes this file automatically, so any one
              of them can be run alone or through run_all.sas.
==============================================================================*/

/*------------------------------------------------------------------------------
  0. Declare every configuration macro variable global first. Programs load
     this file from inside a small macro (%_boot), and without this line each
     %LET below would create a variable local to that macro, which disappears
     as soon as the macro ends.
------------------------------------------------------------------------------*/
%global data_root file_2008 file_2009 file_2010 sas_dir tab_dir fig_dir
        seed origin_dt landmark_dt admin_day dxvars diag_n
        plas_shape plas_chf_a plas_chf_b plas_age1 plas_age2 plas_hr_diab
        plas_hr_female plas_cens_2yr plas_target plas_reps plas_subn run_assess;

/*------------------------------------------------------------------------------
  1. DATA LOCATION
     No folder on your computer is written in this file, so nothing personal
     is committed to the repository. The folder that holds your unzipped
     CMS DE-SynPUF Sample 1 files is set in one of two ways, and is resolved
     in section 3b once the repository root is known.
       a. A file named sas/local_config.sas containing one line,
            %let data_root = C:\path\to\sample_01;
          Copy sas/local_config_template.sas to that name and edit it.
          .gitignore excludes sas/local_config.sas, so it never leaves your
          computer.
       b. The SYNPUF_DATA_ROOT environment variable, used when no
          local_config.sas exists (useful on a shared SAS server).
     Keep the data outside the repository.
------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  3. REPOSITORY ROOT
     Derived from the location of the program being run, so no second path is
     needed. SAS Enterprise Guide and SAS Studio expose that location as
     &_SASPROGRAMFILE. The SAS 9.4 windowing environment exposes it as the
     SAS_EXECFILEPATH environment variable. All programs live in sas/, so the
     repository root is one folder up. It is recomputed every time this file
     runs, so moving or re-pointing the repository takes effect at once. Only
     when no program location is available (for example a batch job) is a
     repo_root set beforehand with %let used instead.
------------------------------------------------------------------------------*/
%macro _resolve_repo_root;
  %global repo_root;
  %local prog pos sasdir;

  %if %symexist(_SASPROGRAMFILE) %then
    %let prog = %qsysfunc(dequote(&_SASPROGRAMFILE));
  %if %length(&prog) = 0 %then %do;
    %if %sysfunc(envlen(SAS_EXECFILEPATH)) > 0 %then
      %let prog = %qsysfunc(sysget(SAS_EXECFILEPATH));
  %end;
  %if %length(&prog) > 0 %then %do;
    %let pos    = %sysfunc(findc(&prog, %str(\/), b));
    %let sasdir = %qsubstr(&prog, 1, %eval(&pos - 1));
    %let pos    = %sysfunc(findc(&sasdir, %str(\/), b));
    %let repo_root = %unquote(%qsubstr(&sasdir, 1, %eval(&pos - 1)));
  %end;
  %else %if %length(&repo_root) = 0 %then %do;
    %put ERROR: [config] Could not determine the program location. Set repo_root before including 00_config.sas.;
    %abort cancel;
  %end;

  %if not %sysfunc(fileexist(&repo_root/sas/00_config.sas)) %then %do;
    %put ERROR: [config] &repo_root does not look like the repository root (sas/00_config.sas not found).;
    %abort cancel;
  %end;
  %put NOTE: [config] repo_root = &repo_root;
%mend _resolve_repo_root;
%_resolve_repo_root

/*------------------------------------------------------------------------------
  3b. DATA LOCATION AND INPUT FILES
     data_root comes from sas/local_config.sas if that file exists, otherwise
     from the SYNPUF_DATA_ROOT environment variable (see section 1). The
     include runs in its own macro so that the checks in _resolve_data_root
     run after the included %LET has taken effect.
------------------------------------------------------------------------------*/
%let data_root = ;
%macro _include_local_config;
  %if %sysfunc(fileexist(&repo_root/sas/local_config.sas)) %then %do;
    %include "&repo_root/sas/local_config.sas";
  %end;
%mend _include_local_config;
%_include_local_config

%macro _resolve_data_root;
  %local src;
  %if %length(&data_root) > 0 %then %let src = sas/local_config.sas;
  %else %if %sysfunc(envlen(SYNPUF_DATA_ROOT)) > 0 %then %do;
    %let data_root = %sysget(SYNPUF_DATA_ROOT);
    %let src = the SYNPUF_DATA_ROOT environment variable;
  %end;
  %if %length(&data_root) = 0 %then %do;
    %put ERROR: [config] The data location is not set. Copy sas/local_config_template.sas to sas/local_config.sas and set data_root, or set the SYNPUF_DATA_ROOT environment variable.;
    %abort cancel;
  %end;
  %if not %sysfunc(fileexist(&data_root)) %then %do;
    %put ERROR: [config] The data folder &data_root (from &src) does not exist.;
    %abort cancel;
  %end;
  %put NOTE: [config] data_root set from &src;
%mend _resolve_data_root;
%_resolve_data_root

/* The three input files, relative to &data_root. CMS ships each CSV inside its
   own folder, and the folder names are not consistent across years (two carry
   a numeric prefix, one does not). If your copy is laid out differently, edit
   only the part after &data_root. */
%let file_2008 = &data_root/176541_DE1_0_2008_Beneficiary_Summary_File_Sample_1/DE1_0_2008_Beneficiary_Summary_File_Sample_1.csv;
%let file_2009 = &data_root/176600_DE1_0_2009_Beneficiary_Summary_File_Sample_1/DE1_0_2009_Beneficiary_Summary_File_Sample_1.csv;
%let file_2010 = &data_root/DE1_0_2010_Beneficiary_Summary_File_Sample_1/DE1_0_2010_Beneficiary_Summary_File_Sample_1.csv;

%let sas_dir = &repo_root/sas;
%let tab_dir = &repo_root/results/tables;
%let fig_dir = &repo_root/results/figures;

/* Git does not track empty folders, so a fresh clone may lack results/,
   results/tables, or results/figures. PROC EXPORT and ODS LISTING fail if the
   folder is missing, so create any that do not exist. */
%macro _make_dir(parent, name);
  %local rc;
  %if not %sysfunc(fileexist(&parent/&name)) %then %do;
    %let rc = %sysfunc(dcreate(&name, &parent));
    %if %length(&rc) = 0 %then %put ERROR: [config] Could not create &parent/&name;
    %else %put NOTE: [config] Created &parent/&name;
  %end;
%mend _make_dir;
%_make_dir(&repo_root, results)
%_make_dir(&repo_root/results, tables)
%_make_dir(&repo_root/results, figures)

/*------------------------------------------------------------------------------
  4. LIBRARY AND OPTIONS
     SYNWORK holds intermediate datasets so that stages can be run in separate
     sessions. It lives in _work/ under the repository root, which .gitignore
     excludes. DLCREATEDIR creates the folder on first use.
------------------------------------------------------------------------------*/
options dlcreatedir;
libname synwork "&repo_root/_work";
options nodlcreatedir;

options validvarname=v7 msglevel=i mprint nomlogic nosymbolgen
        fmtsearch=(synwork work) nofmterr mergenoby=error;

/* Fixed seed for every random step (subsamples, plasmode simulation). */
%let seed = 20081231;

/* Study design constants. Day 1 is January 1, 2009, so a death recorded as
   January 1, 2009 is day 1, the first of each later month falls on a whole day
   count, and administrative censoring on December 31, 2010 is day 730. */
%let origin_dt = '31DEC2008'd;
%let landmark_dt = '01JAN2009'd;
%let admin_day = 730;

/* Covariates of the adjusted working model used by the diagnostics in 07 to
   09. Race enters as three indicators (White is the reference) rather than a
   CLASS variable, because the ASSESS statement and the per-covariate residual
   names of the OUTPUT statement are simpler with plain numeric variables. The
   model is the same as the adjusted model in 06_cox.sas. */
%let dxvars = chf age female race_black race_other race_hisp esrd alzh ckd cancer
              copd depress diab ihd osteo raoa stroke logreimb;

/* Size of the fixed-seed subsample used only by the resampling-based ASSESS
   statement in 07 and 08. Every other diagnostic uses the full cohort. */
%let diag_n = 10000;

/*------------------------------------------------------------------------------
  4b. PLASMODE TRUE PARAMETERS (used by 06b_plasmode.sas, and by 07 and 08 to
      draw the true curves). Time t is in days since December 31, 2008.
        Baseline hazard   Weibull, h0(t) = lambda * k * t**(k - 1), k = plas_shape.
                          lambda is calibrated in 06b so that the expected
                          two-year mortality equals plas_target.
        CHF               log HR(t) = plas_chf_a + plas_chf_b * log(t).
                          The CHF hazard ratio is about 2.5 at one month,
                          1.5 at one year, and 1.3 at two years.
        Age               log HR = plas_age1 * (age - 75) + plas_age2 * (age - 75)**2.
        Diabetes, female  hazard ratios plas_hr_diab and plas_hr_female.
        Everything else   hazard ratio 1 (no effect).
        Random censoring  exponential, with probability plas_cens_2yr of a
                          random censoring time before day 730.
------------------------------------------------------------------------------*/
%let plas_shape     = 1.2;
%let plas_chf_a     = 1.6;
%let plas_chf_b     = -0.2;
%let plas_age1      = 0.03;
%let plas_age2      = 0.0015;
%let plas_hr_diab   = 1.3;
%let plas_hr_female = 0.8;
%let plas_cens_2yr  = 0.05;
%let plas_target    = 0.175;
/* Replicate study: number of simulated data sets and the size of the fixed
   subsample they share. Set plas_reps to 0 to skip it. */
%let plas_reps      = 100;
%let plas_subn      = 20000;
/* The resampling tests of the ASSESS statement (07 and 08) are the slowest
   steps of the pipeline. The proportional hazards test on the plasmode
   subsample alone took about 35 minutes on SAS 9.4M7. They have no Python reference values, so setting run_assess to
   0 skips them without changing the check tally. Their tables and figures are
   then left as they were from the last full run. */
%let run_assess     = 1;

/*------------------------------------------------------------------------------
  5. FORMATS (defined once, used everywhere)
------------------------------------------------------------------------------*/
proc format library=synwork;
  value sexf    1 = 'Male'  2 = 'Female';
  value femalef 0 = 'Male'  1 = 'Female';
  value racef   1 = 'White' 2 = 'Black' 3 = 'Other' 5 = 'Hispanic';
  value ynf     0 = 'No'    1 = 'Yes';
  value chff    0 = 'No CHF' 1 = 'CHF';
  value eventf  0 = 'Censored' 1 = 'Died';
run;

/*------------------------------------------------------------------------------
  6. ODS GRAPHICS
     Figures are written as PNG files to results/figures with fixed names.
     Wrap each figure-producing procedure in %fig(name) and %fig_end.
       - When another destination (HTML in the SAS windowing environment,
         SASREPORT or HTML in Enterprise Guide, and so on) renders the same
         image, ODS gives the LISTING copy a numeric suffix (name1.png). To
         keep the name independent of the session, %fig reads the open
         destinations from SASHELP.VDEST and turns every one except LISTING
         and OUTPUT off with EXCLUDE ALL. Only LISTING renders the figure, so
         the file is exactly name.png. %fig_end turns them back on with
         SELECT ALL. The figure therefore appears only as a file, not in the
         results viewer.
       - %fig first deletes any earlier name.png and name1.png to name9.png.
       - %fig_end confirms with FILEEXIST that exactly name.png exists. If it
         does not, or if a numbered copy appears, it writes a WARNING and
         counts a FAIL in the [check] tally.
------------------------------------------------------------------------------*/
ods graphics on / reset=all imagefmt=png width=7in height=5in;
ods listing gpath="&fig_dir" image_dpi=150;

%macro _delfile(path);
  %local fr rc;
  %if %sysfunc(fileexist(&path)) %then %do;
    %let rc = %sysfunc(filename(fr, &path));
    %let rc = %sysfunc(fdelete(&fr));
    %let rc = %sysfunc(filename(fr));
  %end;
%mend _delfile;

%macro fig(name);
  %global _fig_name _fig_dests;
  %local i d;
  %let _fig_name = &name;
  %_delfile(&fig_dir/&name..png)
  %do i = 1 %to 9;
    %_delfile(&fig_dir/&name.&i..png)
  %end;

  /* Open destinations other than LISTING and OUTPUT, for example
     HTML or TAGSETS.SASREPORT13(EGSR). */
  %let _fig_dests = ;
  proc sql noprint;
    select destination into :_fig_dests separated by '|'
    from sashelp.vdest
    where upcase(destination) not in ('LISTING', 'OUTPUT');
  quit;
  %let i = 1;
  %let d = %qscan(%superq(_fig_dests), &i, |);
  %do %while (%length(&d) > 0);
    ods %unquote(&d) exclude all;
    %let i = %eval(&i + 1);
    %let d = %qscan(%superq(_fig_dests), &i, |);
  %end;

  ods listing select all;
  ods graphics / reset=index imagename="&name";
%mend fig;

%macro fig_end;
  %local i d extra;
  /* Restore every destination that %fig turned off. */
  %let i = 1;
  %let d = %qscan(%superq(_fig_dests), &i, |);
  %do %while (%length(&d) > 0);
    ods %unquote(&d) select all;
    %let i = %eval(&i + 1);
    %let d = %qscan(%superq(_fig_dests), &i, |);
  %end;

  /* Any numbered copy means a second destination still rendered the image. */
  %let extra = ;
  %do i = 1 %to 9;
    %if %sysfunc(fileexist(&fig_dir/&_fig_name.&i..png)) %then
      %let extra = &extra &_fig_name.&i..png;
  %end;

  %if %sysfunc(fileexist(&fig_dir/&_fig_name..png)) and %length(&extra) = 0 %then
    %put NOTE: [fig] Wrote &_fig_name..png;
  %else %do;
    %let _chk_fail = %eval(&_chk_fail + 1);
    %if not %sysfunc(fileexist(&fig_dir/&_fig_name..png)) %then
      %put WARNING: [fig] FAIL &_fig_name..png was not written to &fig_dir;
    %if %length(&extra) > 0 %then
      %put WARNING: [fig] FAIL unexpected numbered copies:&extra;
  %end;
%mend fig_end;

/*------------------------------------------------------------------------------
  7. UTILITY MACROS
------------------------------------------------------------------------------*/

/* Write a small aggregated table to results/tables as CSV. */
%macro export_csv(data=, name=);
  proc export data=&data outfile="&tab_dir/&name..csv" dbms=csv replace;
  run;
%mend export_csv;

/* Python reference values (validation/expected_values.csv, written by
   validation/cross_check.py). Loaded once so programs can check themselves. */
%macro _load_expected;
  %if %sysfunc(fileexist(&repo_root/validation/expected_values.csv)) %then %do;
    data synwork._expected;
      infile "&repo_root/validation/expected_values.csv" dsd firstobs=2 truncover termstr=lf;
      length stage $24 metric $64 vchar $32 value 8;
      input stage $ metric $ vchar $;
      /* Strip a carriage return in case the file was saved with Windows line endings. */
      value = input(compress(vchar, '0D'x), best32.);
      drop vchar;
    run;
  %end;
  %else %put WARNING: [config] validation/expected_values.csv not found. Checks will be skipped.;
%mend _load_expected;
%_load_expected

/* Running tally of check results. Created once per SAS session and reset by
   run_all.sas at the start of a full run, so it survives the reload of this
   file at the top of every program. */
%macro _init_tally;
  %global _chk_pass _chk_fail _chk_skip;
  %if %length(&_chk_pass) = 0 %then %let _chk_pass = 0;
  %if %length(&_chk_fail) = 0 %then %let _chk_fail = 0;
  %if %length(&_chk_skip) = 0 %then %let _chk_skip = 0;
%mend _init_tally;
%_init_tally

/* Record of the metrics compared in this session, used by %check_summary.
   %_reset_checked starts an empty table (run_all.sas calls it at the start
   of a full run). %_record_check adds one name. */
%macro _reset_checked;
  /* PROC SQL rather than a DATA step with STOP, which would write an
     "uninitialized variable" note. */
  proc sql;
    create table synwork._checked (metric char(64));
  quit;
%mend _reset_checked;

%macro _record_check(metric);
  %local notes;
  %let notes = %sysfunc(getoption(notes));
  options nonotes;
  %if not %sysfunc(exist(synwork._checked)) %then %_reset_checked;
  proc sql;
    insert into synwork._checked set metric = "&metric";
  quit;
  options &notes;
%mend _record_check;

/* Compare a value computed in SAS with the Python reference value of the same
   metric and write the result to the log. Call between steps, not inside a
   DATA step.
     PASS  |SAS - Python| <= tol                                      (NOTE)
     FAIL  the values differ, or the SAS value is blank or missing    (WARNING)
     SKIP  no Python value exists for the metric name, which usually
           means a typo in the name                                   (WARNING)
   A missing SAS value is tested explicitly, because a missing number sorts
   below every number and abs(. - x) <= tol would not catch it reliably. */
%macro check(metric=, actual=, tol=0);
  %local expected act_ok;
  %let expected = ;
  %if %sysfunc(exist(synwork._expected)) %then %do;
    proc sql noprint;
      select value format=best16. into :expected trimmed
      from synwork._expected where metric = "&metric" and not missing(value);
    quit;
  %end;
  /* Blank, or a SAS missing value such as . or .A, counts as missing. */
  %let act_ok = 1;
  %if %length(%superq(actual)) = 0 %then %let act_ok = 0;
  %else %if %qsubstr(%superq(actual), 1, 1) = . %then %do;
    %if %length(%superq(actual)) = 1 %then %let act_ok = 0;
    %else %if %sysfunc(anydigit(%qsubstr(%superq(actual), 2, 1))) = 0 %then %let act_ok = 0;
  %end;

  /* Record the metric name in SYNWORK._CHECKED, one row per call, so
     %check_summary can report the reference values that no program checked.
     A dataset is used instead of a growing macro variable, which would have
     to be rewritten (and re-quoted) on every call. NOTES are switched off for
     this one step to keep the log short. */
  %_record_check(&metric)

  %if %length(&expected) = 0 %then %do;
    %let _chk_skip = %eval(&_chk_skip + 1);
    %put WARNING: [check] SKIP &metric: no Python reference value with this name (SAS=&actual);
  %end;
  %else %if &act_ok = 0 %then %do;
    %let _chk_fail = %eval(&_chk_fail + 1);
    %put WARNING: [check] FAIL &metric: SAS value is blank or missing (Python=&expected);
  %end;
  %else %if %sysevalf(%sysfunc(abs(%sysevalf(&actual - &expected))) <= &tol) %then %do;
    %let _chk_pass = %eval(&_chk_pass + 1);
    %put NOTE: [check] PASS &metric: SAS=&actual Python=&expected;
  %end;
  %else %do;
    %let _chk_fail = %eval(&_chk_fail + 1);
    %put WARNING: [check] FAIL &metric: SAS=&actual Python=&expected;
  %end;
%mend check;

/* Range check for quantities that have no single reference value, such as
   the realized death rate of the plasmode data or a true parameter falling
   within a stated distance of its estimate. Counts in the same tally.
     PASS  lo <= actual <= hi    (NOTE)
     FAIL  otherwise or missing  (WARNING) */
%macro check_range(metric=, actual=, lo=, hi=);
  %if %length(%superq(actual)) = 0 %then %do;
    %let _chk_fail = %eval(&_chk_fail + 1);
    %put WARNING: [check] FAIL &metric: SAS value is blank (range &lo to &hi);
  %end;
  %else %if %qsubstr(%superq(actual), 1, 1) = . and %length(%superq(actual)) = 1 %then %do;
    %let _chk_fail = %eval(&_chk_fail + 1);
    %put WARNING: [check] FAIL &metric: SAS value is missing (range &lo to &hi);
  %end;
  %else %if %sysevalf(&actual >= &lo) and %sysevalf(&actual <= &hi) %then %do;
    %let _chk_pass = %eval(&_chk_pass + 1);
    %put NOTE: [check] PASS &metric: SAS=&actual within &lo to &hi;
  %end;
  %else %do;
    %let _chk_fail = %eval(&_chk_fail + 1);
    %put WARNING: [check] FAIL &metric: SAS=&actual outside &lo to &hi;
  %end;
%mend check_range;

/* Deterministic uniform draw from the beneficiary ID, used for every random
   draw of the plasmode simulation and for the fixed subsamples, so that the
   Python cross-check can generate exactly the same values.
     u = (H + 0.5) / 2**52, where H is the integer given by the first 13 hex
     digits (52 bits) of MD5(DESYNPUF_ID || seed || k).
   The 13 digits are the leading 52 bits of the first 8 bytes of the digest.
   Keeping 52 bits makes H + 0.5 exactly representable as a double, so the
   uniform is identical in SAS and Python and lies strictly inside (0, 1).
   K is a number or a numeric expression that names the draw (see 06b). */
%macro hash_u(k, id=DESYNPUF_ID);
  ((input(substr(put(md5(cats(&id, &seed, &k)), $hex32.), 1, 13), hex13.) + 0.5) / 2**52)
%mend hash_u;

/* Fixed subsample of N people: those with the smallest draw-0 uniform
   (SEL_U = %hash_u(0)). Equivalent to a simple random sample, and the same
   people in SAS and Python and in the real and plasmode data. */
%macro hash_sample(data=, out=, n=);
  proc sort data=&data out=&out;
    by sel_u;
  run;
  data &out;
    set &out(obs=&n);
  run;
%mend hash_sample;

/* Analysis views for the diagnostics in 07 to 09. Both views have the same
   variables, so every diagnostic runs unchanged on the real cohort and on the
   plasmode data. SOURCE labels the panels of the side-by-side figures. */
%macro dx_views;
  %local covs;
  %let covs = DESYNPUF_ID chf age age_int female race esrd alzh ckd cancer copd
              depress diab ihd osteo raoa stroke logreimb;
  data synwork._dx_real / view=synwork._dx_real;
    length source $8;
    set synwork.analysis(keep=time event &covs);
    source = 'Real';
    %_dx_derive
  run;
  %if %sysfunc(exist(synwork.plasmode)) %then %do;
    data synwork._dx_plas / view=synwork._dx_plas;
      length source $8;
      set synwork.plasmode(keep=time event &covs);
      source = 'Plasmode';
      %_dx_derive
    run;
  %end;
  %else %put ERROR: [dx] SYNWORK.PLASMODE not found. Run 06b_plasmode.sas first.;
%mend dx_views;

%macro _dx_derive;
  _row = _n_;
  /* Draw 0 orders people for the fixed subsamples (see 06b_plasmode.sas). */
  sel_u = %hash_u(0);
  race_black = (race = 2);
  race_other = (race = 3);
  race_hisp  = (race = 5);
  agec   = age - 75;
  agec2  = agec**2;
  months = time * 24 / &admin_day;
  label race_black = 'Black (vs White)' race_other = 'Other race (vs White)'
        race_hisp  = 'Hispanic (vs White)' agec = 'Age minus 75 (years)'
        agec2      = 'Age minus 75, squared'
        months     = 'Months since January 1, 2009'
        _row       = 'Row number (identifies a record within this view only)';
%mend _dx_derive;

/* ODS Graphics back on with the project settings, after a step that turned
   it off to avoid drawing one plot per covariate. */
%macro graphics_on;
  ods graphics on / imagefmt=png width=7in height=5in;
%mend graphics_on;

/* Prefix every word of a list, for example the per-covariate names of the
   DFBETA= and WTRESSCH= options: %prefix(db_, chf age) is db_chf db_age. */
%macro prefix(pre, list);
  %local i w out;
  %let i = 1;
  %let w = %scan(&list, &i, %str( ));
  %do %while (%length(&w) > 0);
    %let out = &out &pre.&w;
    %let i = %eval(&i + 1);
    %let w = %scan(&list, &i, %str( ));
  %end;
  &out
%mend prefix;

/* One-line summary of every check since the tally was last reset. */
/* Also counts the reference values in validation/expected_values.csv that
   no %check compared in this session, and lists up to 20 of them. After a
   full run_all.sas this should be 0. */
%macro check_summary(errors=0);
  %local unchecked nref;
  %let unchecked = 0;
  %let nref = 0;
  %if %sysfunc(exist(synwork._expected)) %then %do;
    %if not %sysfunc(exist(synwork._checked)) %then %_reset_checked;
    /* Anti-join: reference values with no row in SYNWORK._CHECKED. */
    proc sql noprint;
      create table synwork._unchecked as
      select e.metric from synwork._expected as e
      where not exists (select 1 from synwork._checked as c
                        where upcase(c.metric) = upcase(e.metric))
      order by e.metric;
      select count(*) into :unchecked trimmed from synwork._unchecked;
      select count(*) into :nref trimmed from synwork._expected;
    quit;
    data _null_;
      set synwork._unchecked(obs=20);
      put 'NOTE: [check] Not checked: ' metric;
    run;
  %end;
  /* The line is also kept in _CHK_SUMMARY, so run_all.sas can repeat it in
     the Log window after the log file is closed. A run with errors never
     gets a clean NOTE. SYSCC above 4 means at least one step ended with an
     error, and ERRORS= passes in a count of ERROR lines in the log file
     when the caller has one (run_all.sas does). */
  %global _chk_summary _chk_summary_type;
  %let _chk_summary = [check] SUMMARY PASS=&_chk_pass FAIL=&_chk_fail SKIP=&_chk_skip UNCHECKED=&unchecked of &nref reference values;
  %if &syscc > 4 or &errors > 0 %then %do;
    %let _chk_summary = &_chk_summary. ERRORS OCCURRED (SYSCC=&syscc, ERROR lines=&errors);
    %let _chk_summary_type = WARNING;
  %end;
  %else %if &_chk_fail > 0 or &_chk_skip > 0 or &unchecked > 0 %then %let _chk_summary_type = WARNING;
  %else %let _chk_summary_type = NOTE;
  %put &_chk_summary_type: &_chk_summary;
%mend check_summary;

/* Row count of a dataset into a macro variable. */
%macro nobs(data, mvar);
  %global &mvar;
  proc sql noprint;
    select count(*) into :&mvar trimmed from &data;
  quit;
%mend nobs;

%put NOTE: [config] 00_config.sas loaded.;
