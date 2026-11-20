/*==============================================================================
  Program   : 07_functional_form.sas
  Project   : synpuf-survival-sas
  Purpose   : Functional form of age in the adjusted Cox model, on the real
              cohort and on the plasmode data side by side. The real outcome
              has no age effect to find. The plasmode outcome has a known
              quadratic age effect, so the same diagnostics should find it.
  Inputs    : SYNWORK.ANALYSIS, SYNWORK.PLASMODE (through views from %dx_views)
  Outputs   : SYNWORK.FF_SPLINE   results/tables/ff_spline_tests.csv
              SYNWORK.FF_SUPTEST  results/tables/ff_assess_suptest.csv
              results/figures/ff_martingale_age.png
              results/figures/ff_spline_age.png
              results/figures/ff_assess_age_real.png
              results/figures/ff_assess_age_plasmode.png
  Seminar   : section 6 (functional form: martingale residuals, ASSESS VAR=)
  Notes     : Three views of the same question, each on the full cohort
              unless stated.
                1. Martingale residuals from the model without age, averaged
                   by single year of age. Full cohort.
                2. A restricted cubic spline for age (5 knots at the 5th,
                   27.5th, 50th, 72.5th, and 95th percentiles) and a likelihood
                   ratio test of nonlinearity. Full cohort.
                3. The cumulative martingale residual check of Lin, Wei, and
                   Ying (ASSESS VAR=(age) with RESAMPLE). This is the only part
                   that uses a subsample. Each of the 1,000 resamples and each
                   of the 20 plotted paths runs over every person, so on 61,184
                   people the step is slow and the figure holds over a million
                   points. A fixed subsample of &diag_n people is used instead
                   (%hash_sample, the people with the smallest draw-0 hashed
                   uniforms), and it is the same people in both data sets.
                   ASSESS is run with Breslow ties, the setting under which the
                   method was developed. Its p-values come from resampling
                   inside PROC PHREG, so they are reported without a Python
                   check.
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

/* The adjusted model without age, for the martingale residuals. */
%let dx_noage = chf female race_black race_other race_hisp esrd alzh ckd cancer
                copd depress diab ihd osteo raoa stroke logreimb;

/*------------------------------------------------------------------------------
  1. Martingale residuals by age. If age enters the hazard as f(age), the
     average martingale residual from a model that omits age traces the shape
     of f. A flat line means no age effect.
------------------------------------------------------------------------------*/
%macro mart(src);
  proc phreg data=synwork._dx_&src noprint;
    model time * event(0) = &dx_noage / ties=efron;
    output out=synwork._mres_&src resmart=mres;
  run;

  proc sql;
    create table synwork._mbin_&src as
    select source, min(age_int, 95) as agebin, count(*) as n,
           mean(mres) as mean_mres,
           mean(mres) - probit(0.975) * std(mres) / sqrt(count(*)) as lcl,
           mean(mres) + probit(0.975) * std(mres) / sqrt(count(*)) as ucl
    from synwork._mres_&src
    group by source, calculated agebin;
  quit;
%mend mart;
%mart(real)
%mart(plas)

data synwork._mbin;
  set synwork._mbin_real synwork._mbin_plas;
  label agebin = 'Age on January 1, 2009 (95 = 95 and older)'
        mean_mres = 'Mean martingale residual';
run;

%fig(ff_martingale_age)
title 'Mean martingale residual by age, from the adjusted model without age';
proc sgpanel data=synwork._mbin;
  panelby source / novarname columns=2 sort=descending;
  refline 0 / axis=y lineattrs=(pattern=dot);
  scatter x=agebin y=mean_mres / yerrorlower=lcl yerrorupper=ucl markerattrs=(symbol=circlefilled size=5);
  pbspline x=agebin y=mean_mres / nomarkers lineattrs=(thickness=2);
  colaxis label='Age on January 1, 2009 (95 = 95 and older)' values=(65 to 95 by 5);
  rowaxis label='Mean martingale residual (95% CI)';
run;
title;
%fig_end

/*------------------------------------------------------------------------------
  2. Restricted cubic spline for age. Knots come from the age distribution of
     the real cohort, which the plasmode data share. The basis follows
     Harrell (Regression Modeling Strategies, section 2.4.5), scaled by the
     squared distance between the outer knots.
------------------------------------------------------------------------------*/
proc univariate data=synwork.analysis noprint;
  var age;
  output out=synwork._knots pctlpts=5 27.5 50 72.5 95 pctlpre=k;
run;

%global kn1 kn2 kn3 kn4 kn5;
data _null_;
  set synwork._knots;
  array k[5] k5 k27_5 k50 k72_5 k95;
  do i = 1 to 5;
    call symputx(cats('kn', i), put(k[i], best32.), 'G');
  end;
run;
%check(metric=ff_knot1, actual=&kn1, tol=1e-8)
%check(metric=ff_knot2, actual=&kn2, tol=1e-8)
%check(metric=ff_knot3, actual=&kn3, tol=1e-8)
%check(metric=ff_knot4, actual=&kn4, tol=1e-8)
%check(metric=ff_knot5, actual=&kn5, tol=1e-8)

/* Spline basis for the variable X, written into variables PRE1-PRE3.
   %rcs_knots declares the knots once per DATA step, before any %rcs call. */
%macro rcs_knots;
  array _kn[5] _temporary_ (&kn1 &kn2 &kn3 &kn4 &kn5);
%mend rcs_knots;
%macro rcs(x, pre);
  array &pre.[3] &pre.1 - &pre.3;
  do _j = 1 to 3;
    &pre.[_j] = ( max(&x - _kn[_j], 0)**3
                - max(&x - _kn[4], 0)**3 * (_kn[5] - _kn[_j]) / (_kn[5] - _kn[4])
                + max(&x - _kn[5], 0)**3 * (_kn[4] - _kn[_j]) / (_kn[5] - _kn[4]) )
                / (_kn[5] - _kn[1])**2;
  end;
  drop _j;
%mend rcs;

%macro spline(src);
  data synwork._spl_&src / view=synwork._spl_&src;
    set synwork._dx_&src;
    %rcs_knots
    %rcs(age, s)
  run;

  title "Adjusted model with linear age (&src)";
  ods output FitStatistics = synwork._ffit_lin_&src;
  proc phreg data=synwork._spl_&src;
    model time * event(0) = &dxvars / ties=efron gconv=1e-12;
  run;

  title "Adjusted model with a restricted cubic spline for age (&src)";
  ods output FitStatistics = synwork._ffit_spl_&src
             ParameterEstimates = synwork._fpe_spl_&src
             TestStmts = synwork._ftest_&src;
  proc phreg data=synwork._spl_&src outest=synwork._fest_&src covout;
    model time * event(0) = &dxvars s1 s2 s3 / ties=efron gconv=1e-12;
    nonlinear: test s1, s2, s3;
  run;
  title;

  /* Likelihood ratio test of the three nonlinear terms. The Wald test of
     the same terms is printed by the TEST statement above. */
  proc sql;
    create table synwork._fflr_&src as
    select "&src" as src length=4,
           a.WithCovariates as m2ll_linear,
           b.WithCovariates as m2ll_spline,
           a.WithCovariates - b.WithCovariates as lr_chisq,
           3 as df,
           sdf('chisq', calculated lr_chisq, 3) as lr_p format=pvalue6.4
    from synwork._ffit_lin_&src a, synwork._ffit_spl_&src b
    where compress(upcase(a.Criterion)) = '-2LOGL'
      and compress(upcase(b.Criterion)) = '-2LOGL';
  quit;

  /* Log hazard ratio for age relative to 75, with a 95 percent confidence
     band from the covariance of the age and spline coefficients. */
  data _null_;
    set synwork._fpe_spl_&src;
    _k = whichc(upcase(Parameter), 'AGE', 'S1', 'S2', 'S3');
    if _k then call symputx(cats('_fb', _k), put(Estimate, best32.));
  run;
  /* OUTEST= with COVOUT holds the covariance matrix with one row per
     parameter (_TYPE_ = 'COV', _NAME_ = parameter) and one column per
     parameter. */
  data _null_;
    set synwork._fest_&src(where=(upcase(_TYPE_) = 'COV'));
    _k = whichc(upcase(_NAME_), 'AGE', 'S1', 'S2', 'S3');
    if _k then do;
      call symputx(cats('_fv', _k, '_1'), put(age, best32.));
      call symputx(cats('_fv', _k, '_2'), put(s1,  best32.));
      call symputx(cats('_fv', _k, '_3'), put(s2,  best32.));
      call symputx(cats('_fv', _k, '_4'), put(s3,  best32.));
    end;
  run;

  data synwork._fcurve_&src;
    length source $8;
    source = ifc("&src" = 'real', 'Real', 'Plasmode');
    array b[4] _temporary_ (&_fb1 &_fb2 &_fb3 &_fb4);
    array v[4, 4] _temporary_ (&_fv1_1 &_fv1_2 &_fv1_3 &_fv1_4
                               &_fv2_1 &_fv2_2 &_fv2_3 &_fv2_4
                               &_fv3_1 &_fv3_2 &_fv3_3 &_fv3_4
                               &_fv4_1 &_fv4_2 &_fv4_3 &_fv4_4);
    array d[4] _temporary_;
    %rcs_knots
    /* Basis at the reference age 75. */
    age = 75;
    %rcs(age, r)
    do age = 65 to 100 by 0.5;
      %rcs(age, s)
      d[1] = age - 75;
      d[2] = s1 - r1;
      d[3] = s2 - r2;
      d[4] = s3 - r3;
      est = 0;
      var = 0;
      do i = 1 to 4;
        est = est + d[i] * b[i];
        do j = 1 to 4;
          var = var + d[i] * v[i, j] * d[j];
        end;
      end;
      lcl = est - probit(0.975) * sqrt(var);
      ucl = est + probit(0.975) * sqrt(var);
      if source = 'Plasmode' then
        truth = &plas_age1 * (age - 75) + &plas_age2 * (age - 75)**2;
      else truth = .;
      output;
    end;
    keep source age est lcl ucl truth;
  run;
%mend spline;
%spline(real)
%spline(plas)

data synwork.ff_spline;
  set synwork._fflr_real synwork._fflr_plas;
run;
title 'Test of nonlinearity in age (restricted cubic spline, 3 df)';
proc print data=synwork.ff_spline noobs;
run;
title;
%export_csv(data=synwork.ff_spline, name=ff_spline_tests)

/* Reference checks. Metric names for the plasmode data start with plas_. */
%macro spl_checks(src, pre);
  %local i v b s m2l m2s lr;
  proc sql noprint;
    select m2ll_linear format=best16., m2ll_spline format=best16., lr_chisq format=best16.
      into :m2l trimmed, :m2s trimmed, :lr trimmed
    from synwork.ff_spline where src = "&src";
  quit;
  %check(metric=&pre.ff_linear_m2ll,     actual=&m2l, tol=0.001)
  %check(metric=&pre.ff_spline_m2ll,     actual=&m2s, tol=0.001)
  %check(metric=&pre.ff_lr_nonlin_chisq, actual=&lr,  tol=0.001)
  %do i = 1 %to 4;
    %let v = %scan(age s1 s2 s3, &i);
    proc sql noprint;
      select Estimate format=best16., StdErr format=best16. into :b trimmed, :s trimmed
      from synwork._fpe_spl_&src where upcase(Parameter) = upcase("&v");
    quit;
    %check(metric=&pre.ff_spline_&v._beta, actual=&b, tol=0.00001)
    %check(metric=&pre.ff_spline_&v._se,   actual=&s, tol=0.00001)
  %end;
%mend spl_checks;
%spl_checks(real, )
%spl_checks(plas, plas_)

data synwork._fcurve;
  set synwork._fcurve_real synwork._fcurve_plas;
  label est = 'Spline estimate' truth = 'True effect';
run;

%fig(ff_spline_age)
title 'Log hazard ratio for age relative to age 75 (restricted cubic spline)';
proc sgpanel data=synwork._fcurve;
  panelby source / novarname columns=2 sort=descending;
  band x=age lower=lcl upper=ucl / transparency=0.7 legendlabel='95% CI' name='b';
  series x=age y=est / lineattrs=(thickness=2) legendlabel='Spline estimate' name='e';
  series x=age y=truth / lineattrs=(thickness=2 pattern=dash color=black)
                         legendlabel='True effect (plasmode)' name='t';
  refline 0 / axis=y lineattrs=(pattern=dot);
  colaxis label='Age on January 1, 2009' values=(65 to 100 by 5);
  rowaxis label='Log hazard ratio vs age 75';
  keylegend 'e' 'b' 't' / title='';
run;
title;
%fig_end

/*------------------------------------------------------------------------------
  3. ASSESS VAR=(age) on the fixed subsample. The supremum test asks
     whether the observed cumulative martingale residual path is unusual
     compared with paths simulated under a correct linear age term.
------------------------------------------------------------------------------*/
%macro ff_assess(src, name);
  %hash_sample(data=synwork._dx_&src, out=synwork._sub_&src, n=&diag_n)

  %fig(ff_assess_age_&name)
  title "Cumulative martingale residuals for age (&name data, n = &diag_n)";
  ods output FunctionalFormSupTest = synwork._ffsup_&src;
  proc phreg data=synwork._sub_&src;
    model time * event(0) = &dxvars / ties=breslow;
    assess var=(age) / resample seed=%eval(&seed + 4);
  run;
  title;
  %fig_end
%mend ff_assess;
/* Skipped when run_assess = 0 in 00_config.sas (see the note there). */
%macro ff_assess_all;
  %if &run_assess = 0 %then %do;
    %put NOTE: [assess] run_assess = 0, functional form supremum tests skipped.;
    %return;
  %end;
  %ff_assess(real, real)
  %ff_assess(plas, plasmode)

  data synwork.ff_suptest;
    length source $8;
    set synwork._ffsup_real(in=a) synwork._ffsup_plas;
    source = ifc(a, 'Real', 'Plasmode');
  run;
  title "Supremum test for the functional form of age (subsample n = &diag_n)";
  proc print data=synwork.ff_suptest noobs;
  run;
  title;
  %export_csv(data=synwork.ff_suptest, name=ff_assess_suptest)
%mend ff_assess_all;
%ff_assess_all
