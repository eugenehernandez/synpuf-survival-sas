/*==============================================================================
  Program   : run_all.sas
  Project   : synpuf-survival-sas
  Purpose   : Run the full pipeline end to end, in order.
  Inputs    : the three Beneficiary Summary CSVs named in 00_config.sas
  Outputs   : SYNWORK datasets, CSV tables in results/tables, PNG figures in results/figures
  Seminar   : all sections
  Notes     : Set data_root in sas/local_config.sas (copy it from
              sas/local_config_template.sas), open this file from disk in SAS
              Enterprise Guide, SAS Studio, or the SAS 9.4 windowing
              environment, and submit it. Each program can also be run alone.
              The full log is written to a time-stamped file in _work
              (run_all_YYYYMMDDTHHMMSS.log), because it is too long for the
              Log window of the SAS 9.4 windowing environment. _work is in
              .gitignore, so log files stay out of the repository. When the
              run ends, the Log window shows the log file name, the highest
              SAS condition code (SYSCC), and the tally of every comparison
              with the Python reference in validation/expected_values.csv,
              "[check] SUMMARY PASS=n FAIL=n SKIP=n UNCHECKED=n of N".
              The Log window also shows how many ERROR and WARNING lines the
              log file holds. FAIL, SKIP, UNCHECKED, SYSCC, and the ERROR count
              should all be 0. If any error occurred, the summary line is a
              WARNING that says ERRORS OCCURRED.
==============================================================================*/

/* Start a fresh tally of [check] results for this run. */
%global _chk_pass _chk_fail _chk_skip;
%let _chk_pass = 0;
%let _chk_fail = 0;
%let _chk_skip = 0;

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
  Send the log to a time-stamped file. The name is shown in the Log window
  first, so it can be found even if the run stops early. If a run is
  stopped before the end, submit  proc printto; run;  to bring the log back
  to the Log window. SYSCC is reset so that its value at the end reflects
  this run only.
------------------------------------------------------------------------------*/
%global _run_log;
%let _run_log = &repo_root/_work/run_all_%sysfunc(datetime(), B8601DT15.).log;
%put NOTE: [run_all] Writing the full log to &_run_log;
%let syscc = 0;
proc printto log="&_run_log" new;
run;
/* Empty record of checked metrics for this run (see %check_summary). */
%_reset_checked
%put NOTE: [run_all] Log started %sysfunc(datetime(), datetime20.);

%include "&sas_dir/01_import.sas";
%include "&sas_dir/02_profile.sas";
%include "&sas_dir/03_build_cohort.sas";
%include "&sas_dir/04_explore.sas";
%include "&sas_dir/05_nonparametric.sas";
%include "&sas_dir/06_cox.sas";
%include "&sas_dir/06b_plasmode.sas";
%include "&sas_dir/07_functional_form.sas";
%include "&sas_dir/08_ph_assumption.sas";
%include "&sas_dir/09_influence.sas";
%include "&sas_dir/10_sensitivity.sas";

/* One line that says whether every value matched the Python reference. */
%check_summary
%put NOTE: [run_all] Highest condition code SYSCC=&syscc (0 means no warnings or errors);
%put NOTE: [run_all] Log ended %sysfunc(datetime(), datetime20.);

/* Back to the Log window. */
proc printto;
run;

/* Count the ERROR and WARNING lines in the log file. Real messages start in
   column 1. Echoed source lines start with a line number and MPRINT lines
   with MPRINT, so neither is counted. */
%global _n_err _n_warn;
%let _n_err = 0;
%let _n_warn = 0;
data _null_;
  infile "&_run_log" lrecl=32767 truncover end=eof;
  input;
  if _infile_ =: 'ERROR' then nerr + 1;
  else if _infile_ =: 'WARNING' then nwarn + 1;
  if eof then do;
    call symputx('_n_err', nerr);
    call symputx('_n_warn', nwarn);
  end;
run;

/* Repeat the result in the Log window, now including the error count. */
%put NOTE: [run_all] Full log: &_run_log;
%put NOTE: [run_all] Log file has &_n_err ERROR lines and &_n_warn WARNING lines;
%put NOTE: [run_all] Highest condition code SYSCC=&syscc (0 means no warnings or errors);
%check_summary(errors=&_n_err)
