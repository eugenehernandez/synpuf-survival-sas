# Methods and Design Decisions

This document records every design decision in the analysis, the reason for it, and every place where the analysis departs from the original plan. The README gives the results. The sections follow the order of the programs in `sas/`, and the last section lists the deviations from the plan in one place.

## Data source and profiling

The analysis reads the 2008, 2009, and 2010 Beneficiary Summary Files from the Centers for Medicare and Medicaid Services (CMS) Data Entrepreneurs' Synthetic Public Use File (DE-SynPUF), Sample 1. The claims files are not used. Profiling in `sas/02_profile.sas` (mirrored in `validation/cross_check.py`) established the following facts, and each one shaped a design decision.

- Each yearly file contains exactly the people from the prior year who did not die in that year. No one enters after 2008 and no one disappears without a death date. The loss to follow-up rule is kept in the code but censors no one.
- Every birth date and death date falls on the first of a month. Time to death is therefore resolved to the month, and the primary cohort has at most 24 distinct event times.
- Death dates appear only in the file for the year the death occurred. No beneficiary has conflicting death dates, appears in a later file after death, or has a death date before the birth date.
- The chronic condition flags and the end-stage renal disease (ESRD) indicator are not cumulative. Among people flagged with congestive heart failure (CHF) in 2008, 11,595 are not flagged in 2009. Exposure is therefore fixed at baseline.
- In this release the death indicator is essentially unrelated to every covariate, including age. The primary analysis is expected to return a null result.

## Time scale

Time is counted in days from December 31, 2008, so a death recorded on January 1, 2009 is day 1 and administrative censoring on December 31, 2010 is day 730. The origin is the day before the landmark so that a death on the landmark date itself has a positive time. Because every death date is the first of a month, the 2,040 deaths in the primary cohort fall on 24 distinct days, and roughly 85 deaths share each event time. Ties are therefore heavy, and every Cox model uses Efron ties unless stated otherwise, with Breslow ties fitted for comparison. Plots label the time axis in months (day times 24 / 730), and survival at day 365 means survival through the end of 2009. The kernel-smoothed hazard in `05_nonparametric.sas` uses a bandwidth of 3 months on this months scale, which corresponds to about 91 days.

The median survival time is not reached, because only about 3 percent of the primary cohort dies within two years. The restricted mean survival time (RMST) to day 730 is reported instead. PROC LIFETEST computes it on SAS 9.4M6 and later, and `05_nonparametric.sas` falls back to an equivalent DATA step calculation on earlier releases.

## Cohort

### Landmark design

The condition flags for a year are built from that whole year of claims. A design that starts follow-up on January 1, 2008 and uses the 2008 flags would classify exposure with information gathered during follow-up, which is a look-ahead problem. The analysis therefore uses a landmark. Every covariate comes from the 2008 file, the landmark is January 1, 2009, and follow-up runs to December 31, 2010. A person must be alive on the landmark date and present in the 2009 file. The naive design without a landmark is kept as a sensitivity analysis, described below, and it shows what the landmark prevents.

### Eligibility and the fee-for-service restriction

Eligibility is applied in the order of the attrition table (`results/tables/attrition.csv`). A person must have no death recorded on or before December 31, 2008, appear in the 2009 file, be 65 or older on January 1, 2009 in completed years, and be in fee-for-service Medicare for all of 2008. Fee-for-service means 12 months of Part A coverage, 12 months of Part B coverage, and no Medicare Advantage months. People in Medicare Advantage have claims that do not reach traditional Medicare, so their 2008 condition flags and reimbursement would be incomplete. The restriction removes 34,546 of the 95,730 people aged 65 and older and leaves 61,184. The full cohort of 95,730 without the restriction is kept as a sensitivity analysis.

## Covariates

Every covariate is measured in 2008, before follow-up begins. The exposure is the 2008 CHF flag. The adjusted model also includes age on January 1, 2009 in exact years, female sex, race (Black, Other, and Hispanic, with White as the reference), ESRD, ten other chronic conditions (Alzheimer disease or related dementia, chronic kidney disease, cancer, chronic obstructive pulmonary disease, depression, diabetes, ischemic heart disease, osteoporosis, rheumatoid arthritis or osteoarthritis, and stroke or transient ischemic attack), and total 2008 Medicare reimbursement. Condition flags are recoded from 1 = yes and 2 = no to 1 and 0. Programs `06` and `10` enter race through a CLASS statement or three indicator variables, which give the same model. Programs `07` to `09` use the indicators, because the ASSESS statement and the per-covariate names of the OUTPUT statement are simpler with plain numeric variables.

Total 2008 Medicare reimbursement is the sum of the inpatient, outpatient, and carrier reimbursement fields. A small number of these amounts are negative, which reflects claim adjustments, so each component is floored at 0 before summing and the covariate is log(1 + total).

## Censoring rules

In the primary analysis a death between January 1, 2009 and December 31, 2010 is an event, and everyone else is censored at day 730. A rule censors at day 365 anyone alive at the end of 2009 who is missing from the 2010 file. Profiling showed that no such person exists, so the rule censors no one, and the code keeps it so that the analysis stays correct on data where it would apply.

### Coverage-based censoring

For the coverage-based sensitivity analysis, coverage months are counts rather than dates, so the analysis assumes they run continuously from January. A follow-up year counts as a switch when it has any Medicare Advantage months or fewer than 12 Part A months. Follow-up is censored at the start of the first switching year unless a death is recorded in that year, in which case the death is counted as an event.

People whose first switch is in 2009 and who did not die in 2009 would be censored at day 0 and contribute no follow-up. Rather than keep them as zero-length records, whose handling differs between procedures, the analysis excludes them. This removes 6,109 people from the 61,184 in the primary cohort and leaves 55,075 (see `results/tables/attrition_coverage.csv`). The excluded group includes 107 people who die in 2010, so the coverage rule reduces the number of deaths from 2,040 in the primary analysis to 1,933. Those 107 deaths leave the analysis together with the people who were excluded. Among the 55,075 who remain, the rule never turns a death into a censored record, because a death in the year of a switch always counts as an event. Its effect on them is to end follow-up at day 365 for 6,563 survivors who switch in 2010. The tradeoff the sensitivity analysis measures is therefore a smaller and more selected cohort with shorter follow-up, in exchange for analyzing only time in which fee-for-service coverage, and therefore the claims-derived information, is complete.

## Design sensitivity analyses

`10_sensitivity.sas` fits the adjusted Cox model of `06_cox.sas` in five designs and reports the CHF hazard ratio from each in `results/tables/sensitivity_chf.csv`. The naive design is fitted three ways (unadjusted, adjusted without 2008 reimbursement, and fully adjusted), so the table has seven rows.

| Design | Definition | People |
|---|---|---:|
| Primary | Landmark on January 1, 2009, age 65 or older, fee-for-service in all of 2008 | 61,184 |
| Naive | Time zero January 1, 2008, the 2008 condition flags and 2008 reimbursement, follow-up through December 31, 2010 (day 1,096). Same eligibility as the primary cohort, measured on January 1, 2008 | 59,380 |
| Coverage-based censoring | Primary cohort, censored at the first follow-up year with Medicare Advantage or partial Part A, without the people censored at day 0 | 55,075 |
| All 65 and older | Landmark design without the fee-for-service restriction | 95,730 |
| All ages | Landmark design at every age, with an indicator for younger than 65 | 114,538 |

In the naive design the 2008 condition flags and the 2008 reimbursement are measured during the first year of follow-up, which is the look-ahead problem the landmark design avoids. In this release the share of 2008 beneficiaries with complete fee-for-service coverage in 2008 is the same, about 63 percent, among people who die in 2008 and people who do not, so the fee-for-service rule does not select on 2008 deaths.

| Design | People | Deaths | CHF hazard ratio (95% confidence interval) |
|---|---:|---:|---|
| Primary | 61,184 | 2,040 | 1.024 (0.914 to 1.146) |
| Naive, unadjusted | 59,380 | 2,904 | 1.021 (0.945 to 1.103) |
| Naive, adjusted without 2008 reimbursement | 59,380 | 2,904 | 0.998 (0.909 to 1.096) |
| Naive, fully adjusted | 59,380 | 2,904 | 1.121 (1.017 to 1.237) |
| Coverage-based censoring | 55,075 | 1,933 | 1.006 (0.896 to 1.129) |
| All 65 and older | 95,730 | 3,139 | 0.977 (0.886 to 1.076) |
| All ages | 114,538 | 3,647 | 0.974 (0.890 to 1.066) |

### Why the naive design gives a spurious CHF effect

Four of the five designs give a CHF hazard ratio near 1, as expected in data where death is unrelated to CHF. The fully adjusted naive model gives 1.12. The progression across the three naive rows shows where this comes from. Unadjusted, the hazard ratio is 1.02. Adjusted for every covariate except 2008 reimbursement, it is 0.998. Adding 2008 reimbursement moves it to 1.12.

In the naive design, 2008 reimbursement is measured during follow-up and is partly caused by the outcome. People who die in 2008 have less time to accumulate claims before they die. Their mean log(1 + reimbursement) is 4.31, against 6.23 for people who do not die in 2008. The gap appears in both CHF groups (5.87 against 8.39 with CHF, and 3.60 against 5.18 without CHF). CHF itself is not related to dying in 2008. The share who die in 2008 is 1.53 percent with CHF and 1.63 percent without it.

Adjusting for a covariate that is affected by the outcome biases the coefficients of the other covariates. Low reimbursement marks people who died early, and the model attributes part of that to the covariates that predict reimbursement. CHF is strongly correlated with reimbursement (r = 0.47), so the bias lands on CHF. Among people with the same reimbursement, those with CHF tend to have been the ones whose reimbursement was cut short by an early death, which makes CHF look harmful.

The landmark design prevents this by measuring every covariate in 2008, before follow-up begins on January 1, 2009. No covariate can then be affected by a death that happens during follow-up. The effect in the naive design is modest and its 95 percent confidence interval only narrowly excludes 1, so it should be read as a demonstration of a known mechanism rather than as a precise estimate of the bias. The numbers behind this explanation are in `results/tables/sensitivity_naive_mechanism.csv`, and each is checked against the Python reference values.

## Plasmode simulation

A plasmode data set keeps the real covariates of the primary cohort and replaces the outcome with one simulated from a model whose parameters are known. In this release the real outcome carries almost no signal, so the plasmode data are what show that the models and diagnostics find an effect when one exists. `sas/06b_plasmode.sas` generates the data in a DATA step. Every random number is derived from the beneficiary ID rather than from a random number generator, as described under "Random draws" below, so the simulated outcome is identical on every run, on every SAS release, and in Python. The true parameters are set once, in section 4b of `sas/00_config.sas`, and are listed below. Time t is counted in days from December 31, 2008, as in the primary analysis.

| Component | True value |
|---|---|
| Baseline hazard | Weibull with shape k = 1.2, so the baseline hazard rises slowly over time |
| Baseline scale | Calibrated so that the expected two-year mortality of the cohort, before random censoring, is 17.5 percent. The calibrated value is log(lambda) = -9.947540 |
| Congestive heart failure (CHF) | Log hazard ratio 1.6 - 0.2 log(t). The hazard ratio is about 2.5 at one month, 1.5 at one year, and 1.3 at two years |
| Age | Log hazard ratio 0.03 (age - 75) + 0.0015 (age - 75)^2. Relative to age 75, the hazard ratio is about 0.86 at 65, 1.57 at 85, and 3.3 at 95 |
| Diabetes | Hazard ratio 1.3 |
| Female sex | Hazard ratio 0.8 |
| Every other covariate | Hazard ratio 1 (race, end-stage renal disease, the other chronic conditions, and 2008 reimbursement) |
| Random censoring | Exponential, with a 5 percent probability of a censoring time before day 730 |
| Administrative censoring | Day 730 (December 31, 2010) |
| Date resolution | Every simulated date is moved to the first of its month, as in DE-SynPUF. The daily time is kept for comparison |

Under these values the expected two-year mortality is 24.8 percent for people with CHF and 14.0 percent for people without it. Because exp(b log t) equals t to the power b, the CHF term changes only the Weibull shape for people with CHF (from 1.2 to 1.0). The cumulative hazard therefore has a closed form, H(t) = lambda k exp(eta) t^k_i / k_i, where eta is the linear predictor without the time-varying part and k_i is the person's shape. Each survival time is generated by solving H(t) = -log(U) for a uniform draw U.

The recovery model has the same form as the truth. It includes CHF, CHF by log(time), age minus 75 and its square, and the same covariates as the adjusted model in `06_cox.sas`. It is fitted four ways on the full plasmode cohort.

1. Daily dates with Efron ties. This is the reference fit.
2. Monthly dates with Efron ties.
3. Monthly dates with Breslow ties. Comparing it with the second fit shows the effect of the ties method alone.
4. Monthly dates with Efron ties, with the CHF term evaluated at log(time + 15). A death recorded on the first of a month happened on average about 15 days later, and this fit uses that midpoint.

The time-varying CHF terms are the ones that monthly dates can distort. Every death in January is recorded as day 1, where log(time) is 0, even though those deaths happened between day 1 and day 31. Fits 2 and 3 therefore read the CHF effect at the wrong times, and the size of that distortion is one of the things the recovery table shows. The program checks that fits 1 and 4 place each built-in effect (the two CHF terms, the two age terms, diabetes, and female sex) within 3.29 standard errors of its true value, which is a two-sided 99.9 percent band. It also checks that the realized two-year mortality falls between 15 and 20 percent.

A replicate study redraws the outcome 100 times on one fixed subsample of 20,000 people. It reports bias, the empirical standard deviation of the estimates, the average model standard error, and 95 percent confidence interval coverage for the three monthly-date fits. The daily-date fit is left out of the replicates because, with 730 distinct event times and a programming statement that PROC PHREG evaluates at every event time, it takes far longer than the monthly fits. The number of replicates and the subsample size are set in `sas/00_config.sas`, and setting the number of replicates to 0 skips the study.

### Main result: monthly dates and the time-varying CHF effect

The plasmode results show that monthly date resolution, not the ties method, decides whether the time-varying CHF effect is recovered. When the CHF term is written as CHF times log(time), every death in a month is evaluated at the first day of that month. For January that is day 1, where log(time) is 0, although the deaths happened on days 1 to 31. The fit then misreads how fast the CHF effect decays. Evaluating the term at the middle of the month, log(time + 15), removes most of that error.

On the full plasmode cohort, the daily-date fit gives a = 1.381 and b = -0.163 (true values 1.6 and -0.2), each within 1.8 standard errors of the truth. The monthly-date fit with log(time) gives a = 1.086 and b = -0.114, 5.9 and 5.6 standard errors from the truth. The mid-month fit gives a = 1.436 and b = -0.173, within 1.3 standard errors. Efron and Breslow ties differ by about 0.005 in each coefficient, which is small next to the effect of the date resolution.

The replicate study, 100 simulated outcomes on 20,000 people, shows the same pattern in terms of coverage.

| Fit | Mean of a (true 1.6) | Mean of b (true -0.2) | 95% confidence interval coverage, a and b |
|---|---|---|---|
| Monthly dates, log(time), Efron ties | 1.246 | -0.141 | 0.35 and 0.42 |
| Monthly dates, log(time), Breslow ties | 1.240 | -0.140 | 0.34 and 0.40 |
| Monthly dates, log(time + 15), Efron ties | 1.664 | -0.211 | 0.98 and 0.97 |

In every fit the average model standard error is close to the spread of the estimates across replicates (0.154 and 0.154 for a with log(time), 0.224 and 0.212 with log(time + 15)), so the low coverage comes from bias rather than from underestimated standard errors. Age, age squared, diabetes, and female sex are recovered in every fit, with coverage between 0.91 and 0.97. The results come from `results/tables/plasmode_recovery.csv` and `results/tables/plasmode_replicates.csv`.

### Random draws

Each uniform random number is computed from a hash of the beneficiary ID. For draw k of a person, the program takes the MD5 digest of the text formed by joining the ID, the seed 20081231 (the `seed` setting in `sas/00_config.sas`), and k. The leading 13 hexadecimal digits of the digest, which are the top 52 bits of its first 8 bytes, form an integer H, and the uniform is (H + 0.5) / 2^52. Keeping 52 bits rather than all 64 bits of the first 8 bytes makes the conversion exact in double precision, so SAS and Python produce the same number to the last bit, and the value is always strictly between 0 and 1. The macro `%hash_u` in `sas/00_config.sas` and the function `hash_h` in `validation/cross_check.py` implement the same rule. The draws are numbered as follows.

| Draw k | Use |
|---|---|
| 0 | Orders people for the fixed subsamples. A subsample of n people is the n people with the smallest draw-0 uniforms, which is equivalent to a simple random sample |
| 1 and 2 | Event time and random censoring time in the main plasmode data |
| 2r + 1 and 2r + 2 | Event time and random censoring time in replicate r of the replicate study |

Because Python builds exactly the same data set, the cross-check compares every deterministic plasmode result with SAS. These include the counts of deaths and censored records, every coefficient and standard error of the four recovery fits and the -2 log partial likelihood of each, the bias, empirical standard deviation, average standard error, and coverage of the replicate study, and the functional form and proportional hazards models of `07` and `08` on the plasmode data. A self-test in `06b_plasmode.sas` first compares the 52-bit integers and the uniforms for draws 0, 1, and 2 of the five smallest beneficiary IDs. The continuous survival times can differ between the two languages in the last binary digit, because each language has its own logarithm and power functions, but the recorded times are whole days and are identical.

### Results without a Python check

The supremum tests from the ASSESS statements in `07` and `08` are not compared with Python. Their p-values come from 1,000 residual paths simulated inside PROC PHREG with its own random number generator, which cannot be reproduced outside SAS, and the cross-check does not implement the method. The program reports them in tables and figures only. The influence statistics in `09` are checked. The cross-check computes the score residuals under Efron ties as PROC PHREG defines them, DFBETA as the inverse information times the score residual, and likelihood displacement as the score residual weighted by the inverse information. It compares the largest likelihood displacement, the CHF coefficient before and after removing the 20 people with the largest displacement, and the number of people whose DFBETA moves a coefficient by more than a tenth of its standard error.

## Proportional hazards tests on the plasmode data

`08_ph_assumption.sas` tests each covariate for a hazard ratio that changes over time by adding covariate times a time term to the model, one covariate at a time. The plasmode data are tested in two versions, and both are reported in `results/tables/ph_time_interactions.csv`.

The naive version adds covariate times log(time) to the adjusted working model, which has a linear age term and a constant CHF effect. This is the model an analyst would fit without knowing the truth, and it is also the version applied to the real cohort. The correctly specified version uses the true model form, with age minus 75 and its square and with CHF times log(time + 15) for the decaying CHF effect, and adds covariate times log(time + 15). CHF is not tested in the correctly specified version because its time-varying term is already in the model.

| Covariate | Naive coefficient (SE) | Correctly specified coefficient (SE) |
|---|---|---|
| Diabetes | -0.064 (0.015) | -0.029 (0.024) |
| Female sex | -0.022 (0.015) | -0.028 (0.022) |
| Age | -0.0011 (0.0008) | 0.00003 (0.0012) |

SE is the standard error. Diabetes and female sex were simulated with constant hazard ratios. The naive version flags diabetes (coefficient 4.1 standard errors from 0), which is a false positive. It has two sources. First, the naive model omits the curvature in age and the decay in the CHF effect, and a time-varying term for a covariate that is correlated with age or CHF absorbs part of that misfit. Second, log(time) on dates set to the first of the month distorts the earliest event times, as described above. In the correctly specified version neither diabetes nor female sex is flagged (both within 1.3 standard errors of 0), and the program checks that each lies within 1.96 standard errors. The lesson the example demonstrates is that a proportional hazards test can only be read correctly for a model whose other terms are specified well.

## Influence diagnostics and the flagging threshold

`09_influence.sas` computes DFBETA (the approximate change in each coefficient when one person is removed) and likelihood displacement (LD, the approximate change in the whole fit) for all 61,184 people in each data set. A person is flagged when the absolute DFBETA divided by the coefficient's standard error exceeds 0.1, that is, when one person moves a coefficient by more than a tenth of its standard error. The common size-adjusted cutoff of 2 divided by the square root of n would be 0.008 here, which flags thousands of people whose individual effect is negligible, so an absolute threshold on the standard error scale is used instead. No one is flagged for CHF, age, female sex, or diabetes in either data set, and the largest value is 0.064. As a second check, the model is refitted without the 20 people with the largest LD. The CHF coefficient moves by 0.13 standard errors on the real data and 0.03 standard errors on the plasmode data, and the one-step approximation (the sum of their DFBETA values) predicts the actual change closely. The results are in `results/tables/influence_summary.csv` and `results/tables/influence_refit.csv`.

## Data used by each program

Because death dates resolve only to the month, `05_nonparametric.sas` collapses the cohort to one row per CHF group, follow-up time, and status, with a frequency count. Kaplan-Meier estimates, Nelson-Aalen estimates, and the group tests are identical to the row-level results. From `06_cox.sas` onward every model uses the row-level analysis dataset, because the residual and influence diagnostics in `07` to `09` are computed for each person. The diagnostics in `07` to `09` run the same code on the real cohort and on the plasmode data, and each figure shows the two side by side.

## Diagnostics on a large cohort

Every diagnostic uses the full cohort of 61,184 people except the cumulative residual checks of Lin, Wei, and Ying (1993). Those are the ASSESS statements in `07_functional_form.sas` (functional form of age) and `08_ph_assumption.sas` (proportional hazards). For each one, PROC PHREG simulates 1,000 residual paths over every person to compute the supremum test and draws 20 of them, so on the full cohort the step is slow and each figure holds over a million points. These two checks therefore use a fixed subsample of 10,000 people, the people with the smallest draw-0 uniforms described above, which is the same set of people in the real and plasmode data. The seeds of the resampling inside PROC PHREG are also fixed, so the p-values are stable from run to run on the same SAS release. The ASSESS steps use Breslow ties, the setting under which the method was developed, while every other model uses Efron ties. The sample size is the `diag_n` setting in `sas/00_config.sas`. Even on the subsample, the proportional hazards check on the plasmode data takes about 35 minutes, roughly half of a full run.

The supremum tests are reported in `results/tables/ff_assess_suptest.csv` and `results/tables/ph_assess_suptest.csv`. For the functional form of age, the p-value is 0.059 on the real data and 0.001 on the plasmode data, which has a built-in quadratic age effect. For proportional hazards on the real data, ischemic heart disease has p = 0.010 and every other covariate has p above 0.05. With 18 covariates tested, about one p-value below 0.05 is expected by chance. On the plasmode data, CHF is flagged (p = 0.002) as it should be, but race (Black), ESRD, cancer, and diabetes are also flagged with p-values between 0.011 and 0.047, although none of them has a time-varying effect. The adjusted working model omits the curvature in age and the decay in the CHF effect, and this is the same misspecification that produces the false positive in the time-interaction test described above.

The influence diagnostics in `09_influence.sas` use the full cohort. PROC PHREG computes the one-step approximations to DFBETA (the change in each coefficient when one person is removed) and likelihood displacement for every person quickly, and the scatter plots draw every person. No person-level values are printed or exported. The tables report the distribution of the influence statistics and the change in the CHF coefficient when the 20 most influential people are removed together.

Because there are only 24 event times, the scaled Schoenfeld residuals in `08_ph_assumption.sas` are summarized as a mean for each month with a 95 percent confidence interval, rather than smoothed person by person. The martingale residuals in `07_functional_form.sas` are averaged by single year of age for the same reason of scale. The restricted cubic spline for age has five knots at the 5th, 27.5th, 50th, 72.5th, and 95th percentiles of age in the real cohort (Harrell, *Regression Modeling Strategies*), and the plasmode data use the same knots.

## Deviations from the original plan

The following decisions differ from the original plan or add to it. Each is described in the section named.

1. The primary cohort is restricted to fee-for-service enrollees, and the full 65 and older cohort became a sensitivity analysis (Cohort).
2. Time is counted from December 31, 2008, so that January 1, 2009 is day 1, and plots use a months axis because dates resolve only to the month (Time scale).
3. In the coverage-based censoring analysis, the 6,109 people who would be censored at day 0 are excluded rather than kept with zero follow-up (Censoring rules).
4. The plasmode simulation, the replicate study, and the side-by-side diagnostics were added, because the real outcome carries no signal (Plasmode simulation).
5. The plasmode random numbers are hashed from the beneficiary ID instead of drawn with a seeded random number generator, so that Python can rebuild the same data. The hash keeps 52 of the 64 bits of the first 8 bytes of the digest, so that the conversion to a uniform is exact (Random draws).
6. The naive design is fitted three ways to show the source of its spurious effect (Design sensitivity analyses).
7. The proportional hazards tests on the plasmode data are reported in a naive and a correctly specified version (Proportional hazards tests on the plasmode data).
8. The cumulative residual checks use a fixed subsample of 10,000 people and Breslow ties (Diagnostics on a large cohort).
9. The functional form check with ASSESS covers age only. The plan also named the reimbursement covariate, which was not assessed this way. Age was chosen because the plasmode data have a known nonlinear age effect to recover.
10. No MAXOBS option was set for ODS Graphics. The large scatter plots are drawn with PROC SGPLOT and PROC SGPANEL from output data sets, the ASSESS plots use the subsample, and the full run wrote every figure without a warning.
11. The real data show no violation of proportional hazards for CHF, so no stratification or time interaction was added to the primary model. The time interaction is demonstrated on the plasmode data instead.
12. The influence diagnostics use an absolute threshold of 0.1 standard errors and report LD and DFBETA but not LMAX (Influence diagnostics and the flagging threshold).
13. The Python cross-check uses its own survival code in numpy and scipy rather than the lifelines or statsmodels packages, which were not available. The code is checked against published results for the Freireich leukemia data.
14. Negative reimbursement components are floored at 0 before summing (Covariates).
