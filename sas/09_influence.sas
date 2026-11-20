/*==============================================================================
  Program   : 09_influence.sas
  Project   : synpuf-survival-sas
  Purpose   : Influence diagnostics for the adjusted Cox model, on the real
              cohort and the plasmode data side by side. DFBETA measures how
              much one person moves each coefficient. Likelihood displacement
              (LD) measures how much one person moves the whole fit.
  Inputs    : SYNWORK.ANALYSIS, SYNWORK.PLASMODE (through views from %dx_views)
  Outputs   : SYNWORK.INF_SUMMARY  results/tables/influence_summary.csv
              SYNWORK.INF_REFIT    results/tables/influence_refit.csv
              results/figures/inf_dfbeta_chf.png
              results/figures/inf_ld.png
  Seminar   : section 8 (influential observations: DFBETA and LD)
  Notes     : Everything here uses the full cohort. The OUTPUT statement
              computes the one-step DFBETA and LD approximations for all 61,184
              people quickly, and the scatter plots draw every person. No
              person-level values are printed or exported. The tables report
              only the distribution of the influence statistics and the effect
              of removing the most influential people as a group.
              In a cohort this size no single person should matter much, and
              the results show that for both data sets.
              The flagging threshold is |DFBETA / SE| > 0.1, that is, one
              person moving a coefficient by a tenth of its standard error.
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

/* Number of the most influential people (by LD) removed together in the
   refit check. */
%let n_top = 20;

/*------------------------------------------------------------------------------
  1. DFBETA for every coefficient and LD for every person. SAS reports DFBETA
     as the approximate change b - b(i) when person i is removed. Dividing by
     the standard error puts every coefficient on the same scale.
------------------------------------------------------------------------------*/
%macro infl(src);
  title "Adjusted model for influence diagnostics (&src)";
  ods output ParameterEstimates = synwork._ipe_&src;
  proc phreg data=synwork._dx_&src;
    model time * event(0) = &dxvars / ties=efron gconv=1e-12;
    output out=synwork._inf_&src dfbeta=%prefix(db_, &dxvars) ld=ld xbeta=xb;
  run;
  title;

  proc sql noprint;
    select Estimate format=best16., StdErr format=best16. into :_b_chf trimmed, :_s_chf trimmed
      from synwork._ipe_&src where upcase(Parameter) = 'CHF';
    select StdErr format=best16. into :_s_age trimmed
      from synwork._ipe_&src where upcase(Parameter) = 'AGE';
    select StdErr format=best16. into :_s_fem trimmed
      from synwork._ipe_&src where upcase(Parameter) = 'FEMALE';
    select StdErr format=best16. into :_s_dia trimmed
      from synwork._ipe_&src where upcase(Parameter) = 'DIAB';
  quit;

  data synwork._inf_&src;
    set synwork._inf_&src;
    sdb_chf    = db_chf    / &_s_chf;
    sdb_age    = db_age    / &_s_age;
    sdb_female = db_female / &_s_fem;
    sdb_diab   = db_diab   / &_s_dia;
    label sdb_chf = 'DFBETA for CHF / SE' ld = 'Likelihood displacement'
          xb = 'Linear predictor';
    format event eventf.;
  run;

  /* Distribution of the scaled DFBETA values and of LD. */
  proc sql;
    create table synwork._isum_&src as
    select source, 'CHF' as covariate length=8,
           max(abs(sdb_chf)) as max_abs_sdfbeta, sum(abs(sdb_chf) > 0.1) as n_over_0_1,
           max(ld) as max_ld, sum(ld) as sum_ld, count(*) as n
      from synwork._inf_&src group by source
    union all
    select source, 'Age', max(abs(sdb_age)), sum(abs(sdb_age) > 0.1), max(ld), sum(ld), count(*)
      from synwork._inf_&src group by source
    union all
    select source, 'Female', max(abs(sdb_female)), sum(abs(sdb_female) > 0.1), max(ld), sum(ld), count(*)
      from synwork._inf_&src group by source
    union all
    select source, 'Diabetes', max(abs(sdb_diab)), sum(abs(sdb_diab) > 0.1), max(ld), sum(ld), count(*)
      from synwork._inf_&src group by source;
  quit;

  /* Refit without the &n_top people with the largest LD and compare the
     actual change in the CHF coefficient with the sum of their DFBETA
     values, which is the one-step approximation of that change. */
  proc sort data=synwork._inf_&src(keep=_row ld db_chf) out=synwork._itop_&src;
    by descending ld;
  run;
  data synwork._itop_&src;
    set synwork._itop_&src(obs=&n_top);
  run;
  proc sql noprint;
    select sum(db_chf) format=best16. into :_sumdb trimmed from synwork._itop_&src;
    create view synwork._irefit_&src as
      select * from synwork._dx_&src
      where _row not in (select _row from synwork._itop_&src);
  quit;

  title "Adjusted model without the &n_top people with the largest LD (&src)";
  ods output ParameterEstimates = synwork._ipe2_&src;
  proc phreg data=synwork._irefit_&src;
    model time * event(0) = &dxvars / ties=efron gconv=1e-12;
  run;
  title;

  %local lab;
  %if &src = real %then %let lab = Real;
  %else %let lab = Plasmode;
  proc sql;
    create table synwork._iref_&src as
    select "&lab" as source length=8,
           &n_top as n_removed,
           &_b_chf as beta_chf_all,
           Estimate as beta_chf_without,
           Estimate - &_b_chf as actual_change,
           -&_sumdb as approx_change,
           (Estimate - &_b_chf) / &_s_chf as actual_change_in_se
    from synwork._ipe2_&src where upcase(Parameter) = 'CHF';
  quit;
%mend infl;
%infl(real)
%infl(plas)

data synwork.inf_summary;
  set synwork._isum_real synwork._isum_plas;
  format max_abs_sdfbeta 8.4 max_ld sum_ld 10.4;
  label max_abs_sdfbeta = 'Largest |DFBETA / SE|'
        n_over_0_1 = 'People with |DFBETA / SE| > 0.1';
run;
title 'Influence summary (no single person shown)';
proc print data=synwork.inf_summary noobs label;
run;
title;
%export_csv(data=synwork.inf_summary, name=influence_summary)

data synwork.inf_refit;
  set synwork._iref_real synwork._iref_plas;
  format beta_chf_all beta_chf_without actual_change approx_change 10.6
         actual_change_in_se 8.3;
run;
title "CHF coefficient after removing the &n_top people with the largest LD";
proc print data=synwork.inf_refit noobs;
run;
title;
%export_csv(data=synwork.inf_refit, name=influence_refit)

/* Reference checks for both data sets. The Python cross-check computes the
   same one-step approximations (score residuals under Efron ties, DFBETA as
   the inverse information times the score residual, and LD as the score
   residual weighted by the inverse information), so the largest LD, the
   CHF coefficient before and after removing the &n_top people with the
   largest LD, and the number of people with |DFBETA / SE| > 0.1 are all
   compared. Metric names for the plasmode data contain plas_. */
%macro inf_checks(source, pre);
  %local ld b0 b1 i c lab n;
  proc sql noprint;
    select max(max_ld) format=best16. into :ld trimmed
      from synwork.inf_summary where source = "&source";
    select beta_chf_all format=best16., beta_chf_without format=best16.
      into :b0 trimmed, :b1 trimmed
      from synwork.inf_refit where source = "&source";
  quit;
  %check(metric=inf_&pre.max_ld,                 actual=&ld, tol=1e-6)
  %check(metric=inf_&pre.beta_chf_all,           actual=&b0, tol=0.00001)
  %check(metric=inf_&pre.beta_chf_without_top20, actual=&b1, tol=0.00001)
  %do i = 1 %to 4;
    %let c   = %scan(chf age female diab, &i);
    %let lab = %scan(CHF Age Female Diabetes, &i);
    proc sql noprint;
      select n_over_0_1 into :n trimmed
      from synwork.inf_summary where source = "&source" and covariate = "&lab";
    quit;
    %check(metric=inf_&pre.n_over_0_1_&c, actual=&n)
  %end;
%mend inf_checks;
%inf_checks(Real, )
%inf_checks(Plasmode, plas_)

/*------------------------------------------------------------------------------
  2. Figures. Every person is drawn. ANTIALIASMAX is raised so the scatter
     plots are drawn with antialiasing at this size.
------------------------------------------------------------------------------*/
data synwork._inf;
  set synwork._inf_real(keep=source age event sdb_chf ld xb)
      synwork._inf_plas(keep=source age event sdb_chf ld xb);
run;

%fig(inf_dfbeta_chf)
ods graphics / antialiasmax=130000;
title 'Scaled DFBETA for CHF by age';
proc sgpanel data=synwork._inf;
  panelby source / novarname columns=2 sort=descending;
  scatter x=age y=sdb_chf / group=event markerattrs=(symbol=circle size=3)
          transparency=0.6 name='s';
  refline 0 / axis=y lineattrs=(pattern=dot);
  colaxis label='Age on January 1, 2009';
  rowaxis label='DFBETA for CHF / standard error';
  keylegend 's' / title='';
run;
title;
%fig_end

%fig(inf_ld)
ods graphics / antialiasmax=130000;
title 'Likelihood displacement by linear predictor';
proc sgpanel data=synwork._inf;
  panelby source / novarname columns=2 sort=descending;
  scatter x=xb y=ld / group=event markerattrs=(symbol=circle size=3)
          transparency=0.6 name='s';
  colaxis label='Linear predictor (log relative hazard)';
  rowaxis label='Likelihood displacement';
  keylegend 's' / title='';
run;
title;
%fig_end
