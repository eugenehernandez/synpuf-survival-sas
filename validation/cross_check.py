#!/usr/bin/env python3
"""
cross_check.py
Independent Python reproduction of the key numbers produced by the SAS programs.

The SAS programs are the analysis. This script is a second, independent
implementation of the same cohort rules and survival estimators, so that every
headline number in the README can be checked across two languages.

Usage
    python validation/cross_check.py              # rebuild validation/expected_values.csv
    python validation/cross_check.py --selftest   # check the survival code against published answers

The data location is read the same way SAS reads it. The SYNPUF_DATA_ROOT
environment variable is used if it is set. Otherwise the data_root line of
sas/local_config.sas is used (copy it from sas/local_config_template.sas).
The three file names come from the file block of sas/00_config.sas.

Survival engine
    Kaplan-Meier, Nelson-Aalen, weighted log-rank tests, restricted mean
    survival time, and the Cox model (Efron and Breslow ties, time-varying
    terms, Efron score residuals, DFBETA, and likelihood displacement) are
    implemented below in numpy and scipy. No survival package is used. The
    estimators are checked by --selftest against the published Freireich (1963)
    leukemia results.

Plasmode data
    Every random draw is hashed from the beneficiary ID (hash_h), the same rule
    as %hash_u in 00_config.sas, so this script builds exactly the same
    plasmode data as 06b_plasmode.sas and checks every model fitted to it.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CONFIG = REPO / "sas" / "00_config.sas"
OUT = REPO / "validation" / "expected_values.csv"

LANDMARK = pd.Timestamp("2009-01-01")
ORIGIN = pd.Timestamp("2008-12-31")      # day 1 = January 1, 2009
ADMIN_DAYS = 730                         # December 31, 2010
NAIVE_ADMIN = 1096                       # naive design: day 1 = January 1, 2008, day 1096 = December 31, 2010
CC = {  # SAS name : source column (1 = yes, 2 = no in the source)
    "alzh": "SP_ALZHDMTA", "chf": "SP_CHF", "ckd": "SP_CHRNKIDN", "cancer": "SP_CNCR",
    "copd": "SP_COPD", "depress": "SP_DEPRESSN", "diab": "SP_DIABETES", "ihd": "SP_ISCHMCHT",
    "osteo": "SP_OSTEOPRS", "raoa": "SP_RA_OA", "stroke": "SP_STRKETIA",
}


# ----------------------------------------------------------------------------
# Configuration and import
# ----------------------------------------------------------------------------
def data_root() -> Path:
    """The data folder. The SYNPUF_DATA_ROOT environment variable if it is
    set (for example when Python runs on a different machine from SAS),
    otherwise the data_root line of sas/local_config.sas."""
    import os
    env = os.environ.get("SYNPUF_DATA_ROOT")
    if env:
        return Path(env.replace("\\", "/"))
    local = REPO / "sas" / "local_config.sas"
    if not local.exists():
        sys.exit("The data location is not set. Copy sas/local_config_template.sas to "
                 "sas/local_config.sas and set data_root, or set SYNPUF_DATA_ROOT.")
    m = re.search(r"^\s*%let\s+data_root\s*=\s*(.+?);", local.read_text(), re.I | re.M)
    if not m:
        sys.exit(f"No %let data_root line found in {local}")
    return Path(m.group(1).strip().replace("\\", "/"))


def read_config() -> dict[int, Path]:
    """Read data_root and the three file macro variables from 00_config.sas."""
    txt = CONFIG.read_text()
    root = str(data_root())
    files = {}
    for yr in (2008, 2009, 2010):
        rel = re.search(rf"%let\s+file_{yr}\s*=\s*&data_root\.?(.+?);", txt, re.I).group(1).strip()
        files[yr] = Path(root.replace("\\", "/") + rel.replace("\\", "/"))
    return files


def load(files: dict[int, Path]) -> dict[int, pd.DataFrame]:
    out = {}
    for yr, f in files.items():
        if not f.exists():
            # Allows a folder layout that differs from the CMS download, by
            # searching the data root for the file name.
            alt = sorted(data_root().glob(f"**/{f.name}"))
            if not alt:
                sys.exit(f"Cannot find {f}")
            f = alt[0]
        df = pd.read_csv(f, dtype=str, keep_default_na=False)
        df["birth_dt"] = pd.to_datetime(df.BENE_BIRTH_DT, format="%Y%m%d")
        df["death_dt"] = pd.to_datetime(df.BENE_DEATH_DT.replace("", None), format="%Y%m%d")
        for c in ["BENE_HI_CVRAGE_TOT_MONS", "BENE_SMI_CVRAGE_TOT_MONS",
                  "BENE_HMO_CVRAGE_TOT_MONS", "PLAN_CVRG_MOS_NUM"]:
            df[c] = df[c].astype(int)
        for c in [c for c in df.columns if c.startswith(("MEDREIMB", "BENRES", "PPPYMT"))]:
            df[c] = df[c].astype(float)
        out[yr] = df.set_index("DESYNPUF_ID")
    return out


# ----------------------------------------------------------------------------
# Profiling checks (mirror of 02_profile.sas)
# ----------------------------------------------------------------------------
def profile(d) -> dict[str, float]:
    a, b, c = d[2008], d[2009], d[2010]
    v = {}
    for yr, df in d.items():
        v[f"rows_{yr}"] = len(df)
        v[f"deaths_recorded_{yr}"] = int(df.death_dt.notna().sum())
        v[f"death_outside_file_year_{yr}"] = int((df.death_dt.notna() & (df.death_dt.dt.year != yr)).sum())
        v[f"death_before_birth_{yr}"] = int((df.death_dt < df.birth_dt).sum())
        v[f"neg_medreimb_rows_{yr}"] = int((df[["MEDREIMB_IP", "MEDREIMB_OP", "MEDREIMB_CAR"]] < 0).any(axis=1).sum())
        v[f"death_not_first_of_month_{yr}"] = int((df.death_dt.notna() & (df.death_dt.dt.day != 1)).sum())
    v["ids_2009_not_in_2008"] = len(b.index.difference(a.index))
    v["ids_2010_not_in_2009"] = len(c.index.difference(b.index))
    dead08 = a.index[a.death_dt.notna()]
    dead09 = b.index[b.death_dt.notna()]
    v["appear_after_death_2008"] = int(dead08.isin(b.index).sum() + dead08.isin(c.index).sum())
    v["appear_after_death_2009"] = int(dead09.isin(c.index).sum())
    both = a.join(b, how="inner", lsuffix="_8", rsuffix="_9")
    v["birth_dt_differs_2008_2009"] = int((both.birth_dt_8 != both.birth_dt_9).sum())
    v["sex_differs_2008_2009"] = int((both.BENE_SEX_IDENT_CD_8 != both.BENE_SEX_IDENT_CD_9).sum())
    v["race_differs_2008_2009"] = int((both.BENE_RACE_CD_8 != both.BENE_RACE_CD_9).sum())
    v["chf_1_to_2_2008_2009"] = int(((both.SP_CHF_8 == "1") & (both.SP_CHF_9 == "2")).sum())
    v["esrd_differs_2008_2009"] = int((both.BENE_ESRD_IND_8 != both.BENE_ESRD_IND_9).sum())
    # Conflicting death dates: a person with a death date in more than one file.
    dd = pd.concat([df.death_dt.dropna() for df in d.values()])
    v["ids_with_death_in_2plus_files"] = int((dd.index.value_counts() >= 2).sum())
    ag = age_int(a.birth_dt, LANDMARK)
    v["age_min_at_landmark_2008file"] = int(ag.min())
    v["age_max_at_landmark_2008file"] = int(ag.max())
    v["birth_not_first_of_month_2008"] = int((a.birth_dt.dt.day != 1).sum())
    return v


# ----------------------------------------------------------------------------
# Cohort (mirror of 03_build_cohort.sas)
# ----------------------------------------------------------------------------
def age_int(birth: pd.Series, at: pd.Timestamp) -> pd.Series:
    """Completed years of age, the same as SAS intck('year', birth, at, 'c')."""
    yrs = at.year - birth.dt.year
    before = (birth.dt.month > at.month) | ((birth.dt.month == at.month) & (birth.dt.day > at.day))
    return yrs - before.astype(int)


def build(d):
    a, b, c = d[2008], d[2009], d[2010]
    x = a.copy()
    x["in2009"] = x.index.isin(b.index)
    x["in2010"] = x.index.isin(c.index)
    death = pd.concat([df.death_dt for df in (a, b, c)]).dropna().groupby(level=0).min()
    x["death_dt"] = death.reindex(x.index)

    x["age_int"] = age_int(x.birth_dt, LANDMARK)
    x["age"] = (LANDMARK - x.birth_dt).dt.days / 365.25
    x["female"] = (x.BENE_SEX_IDENT_CD == "2").astype(int)
    x["race"] = x.BENE_RACE_CD.astype(int)
    x["esrd"] = (x.BENE_ESRD_IND == "Y").astype(int)
    for k, col in CC.items():
        x[k] = (x[col] == "1").astype(int)
    reimb = x[["MEDREIMB_IP", "MEDREIMB_OP", "MEDREIMB_CAR"]].clip(lower=0).sum(axis=1)
    x["reimb08"] = reimb
    x["logreimb"] = np.log1p(reimb)
    x["ffs08"] = ((x.BENE_HI_CVRAGE_TOT_MONS == 12) & (x.BENE_SMI_CVRAGE_TOT_MONS == 12)
                  & (x.BENE_HMO_CVRAGE_TOT_MONS == 0)).astype(int)

    # Primary outcome. Day 1 = January 1, 2009. Day 730 = December 31, 2010.
    in_window = x.death_dt.between(pd.Timestamp("2009-01-01"), pd.Timestamp("2010-12-31"))
    x["event"] = in_window.astype(int)
    x["ltfu"] = (x.in2009 & ~x.in2010 & x.death_dt.isna()).astype(int)
    t = np.where(x.event == 1, (x.death_dt - ORIGIN).dt.days, ADMIN_DAYS)
    t = np.where(x.ltfu == 1, 365, t)
    x["time"] = t.astype(int)

    # Coverage-based censoring (sensitivity 2). A follow-up year "switches" when it
    # has any Medicare Advantage months or fewer than 12 Part A months. Censor at the
    # start of the first switching year, unless a death is recorded in that year.
    def sw(df):
        s = (df.BENE_HMO_CVRAGE_TOT_MONS > 0) | (df.BENE_HI_CVRAGE_TOT_MONS < 12)
        return s.reindex(x.index).fillna(False).astype(bool)
    sw09, sw10 = sw(b), sw(c)
    dy = x.death_dt.dt.year
    tc, ec = x.time.copy(), x.event.copy()
    c09 = sw09 & ~(dy == 2009)
    c10 = ~sw09 & sw10 & ~(dy == 2010) & ~(dy == 2009)
    tc[c09], ec[c09] = 0, 0
    tc[c10], ec[c10] = 365, 0
    x["sw09"], x["sw10"] = sw09.astype(int), sw10.astype(int)
    x["time_cov"], x["event_cov"] = tc, ec
    return x


def attrition(x) -> pd.DataFrame:
    rows = []
    n0 = len(x)
    rows.append((1, "Beneficiaries in the 2008 Beneficiary Summary File", n0, None))
    s = x[x.death_dt.isna() | (x.death_dt > pd.Timestamp("2008-12-31"))]
    rows.append((2, "Excluded: death recorded on or before December 31, 2008", len(s), n0 - len(s)))
    s2 = s[s.in2009]
    rows.append((3, "Excluded: alive but absent from the 2009 file", len(s2), len(s) - len(s2)))
    s3 = s2[s2.age_int >= 65]
    rows.append((4, "Excluded: younger than 65 on January 1, 2009", len(s3), len(s2) - len(s3)))
    s4 = s3[s3.ffs08 == 1]
    rows.append((5, "Excluded: not 12 months Part A and Part B with no Medicare Advantage in 2008",
                 len(s4), len(s3) - len(s4)))
    df = pd.DataFrame(rows, columns=["step", "description", "n_remaining", "n_excluded"])
    return df, s2, s3, s4


# ----------------------------------------------------------------------------
# Survival estimators
# ----------------------------------------------------------------------------
def km(time, event, at=None):
    """Kaplan-Meier with Greenwood standard errors. Returns a table by event time."""
    time, event = np.asarray(time, float), np.asarray(event, int)
    ut = np.unique(time[event == 1])
    n = np.array([(time >= t).sum() for t in ut])
    dth = np.array([((time == t) & (event == 1)).sum() for t in ut])
    s = np.cumprod(1 - dth / n)
    gw = np.cumsum(dth / (n * (n - dth)))
    tab = pd.DataFrame({"time": ut, "at_risk": n, "events": dth, "surv": s, "se": s * np.sqrt(gw)})
    if at is None:
        return tab
    res = {}
    for t in at:
        k = tab[tab.time <= t]
        res[t] = (1.0, 0.0) if k.empty else (k.surv.iloc[-1], k.se.iloc[-1])
    return res


def nelson_aalen(time, event, at):
    time, event = np.asarray(time, float), np.asarray(event, int)
    ut = np.unique(time[event == 1])
    h = np.cumsum([((time == t) & (event == 1)).sum() / (time >= t).sum() for t in ut])
    return {t: (h[ut <= t][-1] if (ut <= t).any() else 0.0) for t in at}


def rmst(time, event, tau):
    tab = km(time, event)
    t = np.concatenate([[0], tab.time[tab.time <= tau], [tau]])
    s = np.concatenate([[1], tab.surv[tab.time <= tau]])
    return float(np.sum(np.diff(t) * s))


def logrank(time, event, group):
    """Two-or-more-group log-rank chi-square (the same statistic as PROC LIFETEST)."""
    time, event, group = np.asarray(time, float), np.asarray(event, int), np.asarray(group)
    gs = np.unique(group)
    ut = np.unique(time[event == 1])
    U = np.zeros(len(gs))
    V = np.zeros((len(gs), len(gs)))
    for t in ut:
        r = time >= t
        n = r.sum()
        dd = ((time == t) & (event == 1)).sum()
        ng = np.array([(r & (group == g)).sum() for g in gs])
        dg = np.array([((time == t) & (event == 1) & (group == g)).sum() for g in gs])
        U += dg - dd * ng / n
        if n > 1:
            V += dd * (n - dd) / (n - 1) * (np.diag(ng / n) - np.outer(ng, ng) / n**2)
    chi = float(U[:-1] @ np.linalg.solve(V[:-1, :-1], U[:-1]))
    from scipy.stats import chi2
    return chi, len(gs) - 1, float(chi2.sf(chi, len(gs) - 1))


def weighted_logrank(time, event, group, weight="logrank"):
    """Two-group weighted log-rank chi-square, matching PROC LIFETEST STRATA tests.

    logrank  w = 1
    wilcoxon w = n (Gehan-Breslow generalized Wilcoxon, the LIFETEST "Wilcoxon")
    tarone   w = sqrt(n) (Tarone-Ware)
    peto     w = S~(t), where S~ is the product of (1 - d/(n + 1)) over earlier times
    """
    time, event, group = np.asarray(time, float), np.asarray(event, int), np.asarray(group)
    g1 = group == np.unique(group)[-1]
    ut = np.unique(time[event == 1])
    U = V = 0.0
    stilde = 1.0
    for t in ut:
        r = time >= t
        n = r.sum(); n1 = (r & g1).sum()
        dd = ((time == t) & (event == 1)).sum(); d1 = ((time == t) & (event == 1) & g1).sum()
        if weight == "logrank":
            w = 1.0
        elif weight == "wilcoxon":
            w = float(n)
        elif weight == "tarone":
            w = np.sqrt(n)
        elif weight == "peto":
            stilde *= 1 - dd / (n + 1)
            w = stilde
        U += w * (d1 - dd * n1 / n)
        if n > 1:
            V += w**2 * dd * (n - dd) / (n - 1) * (n1 / n) * (1 - n1 / n)
    from scipy.stats import chi2
    chi = U**2 / V
    return float(chi), float(chi2.sf(chi, 1))


def sas_pctl(x, q):
    """Percentile by the SAS default definition (PCTLDEF=5)."""
    x = np.sort(np.asarray(x, float)); n = len(x)
    np_ = n * q / 100; j = int(np.floor(np_)); g = np_ - j
    if g > 0:
        return x[j]
    return (x[j - 1] + x[j]) / 2 if j > 0 else x[0]


def cox(X, time, event, ties="efron", tol=1e-9, maxit=50):
    """Cox partial likelihood by Newton-Raphson with exact gradient and Hessian.

    Computation is grouped by distinct time, so it is fast when times are coarse
    (this data has at most 25 distinct times). Returns (beta, se, loglik).
    """
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    time, event = np.asarray(time, float), np.asarray(event, int)
    n, p = X.shape
    ut, inv = np.unique(time, return_inverse=True)
    U = len(ut)
    ev_t = np.bincount(inv, weights=event, minlength=U)
    xd_sum = np.zeros((U, p))
    np.add.at(xd_sum, inv, X * event[:, None])
    beta = np.zeros(p)
    for it in range(maxit):
        eta = X @ beta
        r = np.exp(eta)
        # Sums by distinct time over everyone, then reverse cumulative = risk set.
        s0 = np.bincount(inv, weights=r, minlength=U)
        s1 = np.zeros((U, p)); np.add.at(s1, inv, X * r[:, None])
        s2 = np.zeros((U, p, p)); np.add.at(s2, inv, np.einsum("i,ij,ik->ijk", r, X, X))
        re = r * event
        d0 = np.bincount(inv, weights=re, minlength=U)
        d1 = np.zeros((U, p)); np.add.at(d1, inv, X * re[:, None])
        d2 = np.zeros((U, p, p)); np.add.at(d2, inv, np.einsum("i,ij,ik->ijk", re, X, X))
        R0 = np.cumsum(s0[::-1])[::-1]
        R1 = np.cumsum(s1[::-1], axis=0)[::-1]
        R2 = np.cumsum(s2[::-1], axis=0)[::-1]
        ll = float((eta * event).sum())
        g = xd_sum.sum(0).copy()
        H = np.zeros((p, p))
        for k in np.nonzero(ev_t)[0]:
            dk = int(ev_t[k])
            fs = np.arange(dk) / dk if ties == "efron" else np.zeros(dk)
            for f in fs:
                den = R0[k] - f * d0[k]
                num1 = R1[k] - f * d1[k]
                num2 = R2[k] - f * d2[k]
                ll -= np.log(den)
                g -= num1 / den
                H -= num2 / den - np.outer(num1, num1) / den**2
        if it > 0 and abs(ll - ll_old) < tol:
            break
        ll_old = ll
        step = np.linalg.solve(-H, g)
        beta = beta + step
    cov = np.linalg.inv(-H)
    return beta, np.sqrt(np.diag(cov)), ll


def efron_score_residuals(X, time, event, beta):
    """Score residuals of a Cox model with fixed covariates under Efron ties.

    Follows the SAS/STAT PHREG definition: at an event time with d deaths the
    k-th of d steps uses risk weights 1 - (k - 1) / d for the people who die
    there, and each person's score process is averaged over the d steps.
    """
    X = np.asarray(X, float); time = np.asarray(time, float); event = np.asarray(event, int)
    n, p = X.shape
    r = np.exp(X @ beta)
    ut = np.unique(time[event == 1])
    A = np.zeros(len(ut)); B = np.zeros((len(ut), p))      # every at-risk person
    Ad = np.zeros(len(ut)); Bd = np.zeros((len(ut), p))    # people dying at t
    Cbar = np.zeros((len(ut), p))                          # mean of Zbar_k over k
    order = np.argsort(-time, kind="stable")
    for j, t in enumerate(ut):
        R = time >= t
        D = R & (time == t) & (event == 1)
        d = int(D.sum())
        s0 = r[R].sum(); s1 = r[R] @ X[R]
        d0 = r[D].sum(); d1 = r[D] @ X[D]
        for k in range(d):
            f = k / d
            S0 = s0 - f * d0; S1 = s1 - f * d1; zb = S1 / S0
            w = 1 - f
            A[j] += 1 / S0; B[j] += zb / S0
            Ad[j] += w / S0; Bd[j] += w * zb / S0
            Cbar[j] += zb / d
    CA = np.cumsum(A); CB = np.cumsum(B, axis=0)
    # Index of the last event time at or before each person's time.
    idx = np.searchsorted(ut, time, side="right") - 1
    L = np.zeros((n, p))
    has = idx >= 0
    L[has] = -r[has, None] * (X[has] * CA[idx[has], None] - CB[idx[has]])
    dead = event == 1
    j = idx[dead]
    # For deaths, replace the step at their own time by the Efron-weighted step.
    L[dead] += r[dead, None] * (X[dead] * A[j, None] - B[j]) - r[dead, None] * (X[dead] * Ad[j, None] - Bd[j])
    L[dead] += X[dead] - Cbar[j]
    return L


def cox_influence(X, time, event, ties="efron"):
    """DFBETA (I^-1 L_i) and likelihood displacement (L_i' I^-1 L_i), as in PROC PHREG."""
    beta, se, _ = cox(X, time, event, ties)
    L = efron_score_residuals(X, time, event, beta)
    cov = cox_cov(X, time, event, beta)
    dfb = L @ cov
    ld = np.einsum("ij,ij->i", dfb, L)
    return beta, se, dfb, ld


def cox_cov(X, time, event, beta):
    """Inverse observed information (Efron) at beta."""
    X = np.asarray(X, float); time = np.asarray(time, float); event = np.asarray(event, int)
    p = X.shape[1]; r = np.exp(X @ beta); H = np.zeros((p, p))
    for t in np.unique(time[event == 1]):
        R = time >= t; D = R & (time == t) & (event == 1); d = int(D.sum())
        s0 = r[R].sum(); s1 = r[R] @ X[R]; s2 = (X[R] * r[R, None]).T @ X[R]
        d0 = r[D].sum(); d1 = r[D] @ X[D]; d2 = (X[D] * r[D, None]).T @ X[D]
        for k in range(d):
            f = k / d
            S0 = s0 - f * d0; S1 = s1 - f * d1; S2 = s2 - f * d2
            H -= S2 / S0 - np.outer(S1, S1) / S0**2
    return np.linalg.inv(-H)


def cox_tv(X, time, event, tv, ties="efron", tol=1e-9, maxit=50):
    """Cox model with time-varying columns, by direct risk-set sums.

    X holds the fixed covariates. tv is a list of (column_index, function) pairs:
    at event time t the column becomes column * function(t), which matches a
    PROC PHREG programming statement such as chf_lt = chf * log(time).
    Practical because event times are few (24 here).
    """
    X = np.asarray(X, float)
    time, event = np.asarray(time, float), np.asarray(event, int)
    p = X.shape[1]
    ut = np.unique(time[event == 1])
    beta = np.zeros(p)
    for it in range(maxit):
        ll = 0.0; g = np.zeros(p); H = np.zeros((p, p))
        for t in ut:
            Xt = X.copy()
            for j, f in tv:
                Xt[:, j] = X[:, j] * f(t)
            R = time >= t
            D = R & (time == t) & (event == 1)
            d = int(D.sum())
            eta = Xt[R] @ beta; r = np.exp(eta)
            XR = Xt[R]
            S0 = r.sum(); S1 = r @ XR; S2 = (XR * r[:, None]).T @ XR
            etaD = Xt[D] @ beta; rD = np.exp(etaD); XD = Xt[D]
            D0 = rD.sum(); D1 = rD @ XD; D2 = (XD * rD[:, None]).T @ XD
            ll += etaD.sum(); g += XD.sum(0)
            fs = np.arange(d) / d if ties == "efron" else np.zeros(d)
            for f in fs:
                den = S0 - f * D0; n1 = S1 - f * D1; n2 = S2 - f * D2
                ll -= np.log(den); g -= n1 / den
                H -= n2 / den - np.outer(n1, n1) / den**2
        if it > 0 and abs(ll - ll_old) < tol:
            break
        ll_old = ll
        beta = beta + np.linalg.solve(-H, g)
    return beta, np.sqrt(np.diag(np.linalg.inv(-H))), ll


def plasmode_config() -> dict[str, float]:
    """Read the plasmode true parameters from section 4b of 00_config.sas."""
    txt = CONFIG.read_text()
    return {m.group(1): float(m.group(2))
            for m in re.finditer(r"%let\s+(plas_\w+)\s*=\s*([-0-9.eE]+)\s*;", txt)}


def read_seed() -> int:
    return int(re.search(r"%let\s+seed\s*=\s*(\d+)\s*;", CONFIG.read_text()).group(1))


def hash_h(ids, seed, k) -> np.ndarray:
    """52-bit integers from MD5(ID || seed || k), the same as %hash_u in 00_config.sas.

    The leading 13 hex digits of the digest (the top 52 bits of its first 8
    bytes) are used, so (H + 0.5) / 2**52 is exact in double precision.
    """
    import hashlib
    return np.array([int(hashlib.md5(f"{i}{seed}{k}".encode("ascii")).hexdigest()[:13], 16)
                     for i in ids], dtype=np.int64)


MONTH_START = np.array([(pd.Timestamp(2009, 1, 1) + pd.DateOffset(months=m) - ORIGIN).days
                        for m in range(24)])


def plasmode_generate(ids, eta0, shape, loglam, pc, seed, ku, kc):
    """Mirror of %_plas_generate in 06b_plasmode.sas.

    Returns daily time, monthly (first of month) time, event, and censoring
    type (0 death, 1 random censoring, 2 administrative).
    """
    u = (hash_h(ids, seed, ku) + 0.5) / 2**52
    v = (hash_h(ids, seed, kc) + 0.5) / 2**52
    crate = -np.log(1 - pc["plas_cens_2yr"]) / ADMIN_DAYS
    t = (-np.log(u) * shape / (np.exp(loglam) * pc["plas_shape"] * np.exp(eta0))) ** (1 / shape)
    c = -np.log(v) / crate
    x = np.minimum(t, c)
    event = ((t <= c) & (t <= ADMIN_DAYS)).astype(int)
    admin = (event == 0) & (x >= ADMIN_DAYS)
    tday = np.where(admin, ADMIN_DAYS, np.maximum(1, np.ceil(np.minimum(x, ADMIN_DAYS)))).astype(int)
    tmon = np.where(admin, ADMIN_DAYS, MONTH_START[np.searchsorted(MONTH_START, tday, side="right") - 1])
    ctype = np.where(event == 1, 0, np.where(admin, 2, 1))
    return tday, tmon.astype(int), event, ctype


def rcs_basis(x, knots):
    """Restricted cubic spline basis (Harrell), scaled by (t_k - t_1)**2.

    Matches the %rcs macro in 07_functional_form.sas.
    """
    x = np.asarray(x, float); t = np.asarray(knots, float); k = len(t)
    p3 = lambda u: np.maximum(u, 0) ** 3
    cols = []
    for j in range(k - 2):
        cols.append((p3(x - t[j]) - p3(x - t[k - 2]) * (t[k - 1] - t[j]) / (t[k - 1] - t[k - 2])
                     + p3(x - t[k - 1]) * (t[k - 2] - t[j]) / (t[k - 1] - t[k - 2])) / (t[k - 1] - t[0]) ** 2)
    return np.column_stack(cols)


# ----------------------------------------------------------------------------
# Self-test against published answers
# ----------------------------------------------------------------------------
def selftest():
    """Freireich et al. (1963) leukemia remission data, 42 children.

    Published reference values (Kleinbaum and Klein, Survival Analysis, 3rd ed.,
    chapters 2 and 3, and R survival::coxph on the same data)
      6-MP group Kaplan-Meier S(t) at 6, 7, 10, 13, 16, 22, 23 =
          .857 .807 .753 .690 .627 .538 .448
      Log-rank chi-square = 16.79 (1 df). Wilcoxon (Gehan) 13.46, Tarone-Ware
      15.12, Peto 14.08 (Kleinbaum and Klein, chapter 2)
      Cox coefficient for placebo vs 6-MP, Breslow ties = 1.509 (SE 0.410, HR 4.523)
      Cox coefficient for placebo vs 6-MP, Efron ties   = 1.572 (SE 0.412, HR 4.817)
    """
    t_mp = [6, 6, 6, 7, 10, 13, 16, 22, 23, 6, 9, 10, 11, 17, 19, 20, 25, 32, 32, 34, 35]
    e_mp = [1] * 9 + [0] * 12
    t_pl = [1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12, 15, 17, 22, 23]
    e_pl = [1] * 21
    t = np.array(t_mp + t_pl); e = np.array(e_mp + e_pl); rx = np.array([0] * 21 + [1] * 21)

    pub_km = {6: .857, 7: .807, 10: .753, 13: .690, 16: .627, 22: .538, 23: .448}
    est_km = km(t[rx == 0], e[rx == 0], at=list(pub_km))
    chi, _, _ = logrank(t, e, rx)
    bE, seE, _ = cox(rx, t, e, "efron")
    bB, seB, _ = cox(rx, t, e, "breslow")

    # Each published value is compared at its own printed precision (half a unit
    # in the last reported decimal place).
    rows = [(f"KM 6-MP S({k})", v, 3, est_km[k][0]) for k, v in pub_km.items()]
    rows += [("Log-rank chi-square", 16.79, 2, chi),
             ("Log-rank chi-square, weighted code", 16.79, 2, weighted_logrank(t, e, rx, "logrank")[0]),
             ("Wilcoxon chi-square", 13.46, 2, weighted_logrank(t, e, rx, "wilcoxon")[0]),
             ("Tarone-Ware chi-square", 15.12, 2, weighted_logrank(t, e, rx, "tarone")[0]),
             ("Peto chi-square", 14.08, 2, weighted_logrank(t, e, rx, "peto")[0]),
             ("RMST to 23, 6-MP (hand calculation, not published)", 17.909, 3, rmst(t[rx == 0], e[rx == 0], 23)),
             ("Cox beta, Breslow", 1.509, 3, bB[0]), ("Cox SE, Breslow", 0.410, 3, seB[0]),
             ("Cox HR, Breslow", 4.523, 3, np.exp(bB[0])),
             ("Cox beta, Efron", 1.572, 3, bE[0]), ("Cox SE, Efron", 0.412, 3, seE[0]),
             ("Cox HR, Efron", 4.817, 3, np.exp(bE[0]))]
    df = pd.DataFrame(rows, columns=["quantity", "published", "decimals", "this_code"])
    df["abs_diff"] = (df.published - df.this_code).abs()
    df["match"] = df.abs_diff <= 0.5 * 10.0 ** (-df.decimals) + 1e-12
    pd.set_option("display.width", 120)
    print(df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    ok = bool(df.match.all())
    print("\nSELFTEST", "PASSED" if ok else "FAILED")
    return ok


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        sys.exit(0 if selftest() else 1)

    d = load(read_config())
    vals = []
    add = lambda stage, name, v: vals.append({"stage": stage, "metric": name, "value": v})

    for k, v in profile(d).items():
        add("02_profile", k, v)

    x = build(d)
    att, alive09, full65, prim = attrition(x)
    for _, r in att.iterrows():
        add("03_cohort", f"attrition_step{r.step}_n", r.n_remaining)
    add("03_cohort", "ltfu_censored_2009_primary", int(prim.ltfu.sum()))
    for name, df in [("primary", prim), ("full65", full65), ("all_ages", alive09)]:
        add("03_cohort", f"n_{name}", len(df))
        add("03_cohort", f"events_{name}", int(df.event.sum()))
        add("03_cohort", f"chf_n_{name}", int(df.chf.sum()))
        add("03_cohort", f"chf_events_{name}", int(df.event[df.chf == 1].sum()))
        add("03_cohort", f"sum_time_{name}", int(df.time.sum()))
    add("03_cohort", "n_under65_alive09", int((alive09.age_int < 65).sum()))
    add("03_cohort", "events_under65_alive09", int(alive09.event[alive09.age_int < 65].sum()))
    add("03_cohort", "mean_age_primary", round(prim.age.mean(), 4))
    add("03_cohort", "mean_logreimb_primary", round(prim.logreimb.mean(), 4))
    add("03_cohort", "female_n_primary", int(prim.female.sum()))
    add("03_cohort", "esrd_n_primary", int(prim.esrd.sum()))
    add("03_cohort", "events_at_or_before_day365_primary", int(((prim.event == 1) & (prim.time <= 365)).sum()))
    add("03_cohort", "cov_censored_day0_primary", int(((prim.event_cov == 0) & (prim.time_cov == 0)).sum()))
    add("03_cohort", "cov_censored_day365_primary", int(((prim.event_cov == 0) & (prim.time_cov == 365)).sum()))
    add("03_cohort", "cov_events_primary", int(prim.event_cov.sum()))
    add("03_cohort", "cov_sum_time_primary", int(prim.time_cov.sum()))

    # Coverage-censoring cohort: people censored at day 0 are excluded.
    cov = prim[~((prim.time_cov == 0) & (prim.event_cov == 0))]
    add("03_cohort", "attrition_cov_step6_n", len(cov))
    add("03_cohort", "n_cov", len(cov))
    add("03_cohort", "events_cov", int(cov.event_cov.sum()))
    add("03_cohort", "chf_n_cov", int(cov.chf.sum()))
    add("03_cohort", "chf_events_cov", int(cov.event_cov[cov.chf == 1].sum()))
    add("03_cohort", "sum_time_cov", int(cov.time_cov.sum()))
    excl = prim[(prim.time_cov == 0) & (prim.event_cov == 0)]
    add("03_cohort", "cov_excluded_primary_deaths", int(excl.event.sum()))
    add("03_cohort", "cov_kept_deaths_converted_to_censored", int((cov.event - cov.event_cov).sum()))

    # ---------------- Stage 3: 04_explore (seminar section 2) ----------------
    st = "04_explore"
    for lab, sub in [("chf1", prim[prim.chf == 1]), ("chf0", prim[prim.chf == 0]), ("all", prim)]:
        add(st, f"t1_n_{lab}", len(sub))
        add(st, f"t1_age_mean_{lab}", round(sub.age.mean(), 4))
        add(st, f"t1_age_sd_{lab}", round(sub.age.std(ddof=1), 4))
        add(st, f"t1_female_n_{lab}", int(sub.female.sum()))
        add(st, f"t1_white_n_{lab}", int((sub.race == 1).sum()))
        add(st, f"t1_esrd_n_{lab}", int(sub.esrd.sum()))
        add(st, f"t1_diab_n_{lab}", int(sub.diab.sum()))
        add(st, f"t1_ckd_n_{lab}", int(sub.ckd.sum()))
        add(st, f"t1_reimb_median_{lab}", sas_pctl(sub.reimb08, 50))
        add(st, f"t1_reimb_p25_{lab}", sas_pctl(sub.reimb08, 25))
        add(st, f"t1_reimb_p75_{lab}", sas_pctl(sub.reimb08, 75))
        add(st, f"t1_logreimb_mean_{lab}", round(sub.logreimb.mean(), 4))
        add(st, f"t1_deaths_{lab}", int(sub.event.sum()))
    ev = prim[prim.event == 1].time
    add(st, "evtime_mean", round(ev.mean(), 4))
    add(st, "evtime_median", sas_pctl(ev, 50))
    add(st, "evtime_sd", round(ev.std(ddof=1), 4))
    add(st, "corr_age_logreimb", round(np.corrcoef(prim.age, prim.logreimb)[0, 1], 4))
    add(st, "corr_chf_ckd", round(np.corrcoef(prim.chf, prim.ckd)[0, 1], 4))
    add(st, "corr_chf_ihd", round(np.corrcoef(prim.chf, prim.ihd)[0, 1], 4))

    # ---------------- Stage 3: 05_nonparametric (seminar section 3) ------------
    st = "05_nonparametric"
    for lab, sub in [("chf1", prim[prim.chf == 1]), ("chf0", prim[prim.chf == 0]), ("all", prim)]:
        k = km(sub.time, sub.event, at=[365, 730])
        h = nelson_aalen(sub.time, sub.event, at=[365, 730])
        for t in (365, 730):
            add(st, f"km_s{t}_{lab}", round(k[t][0], 6))
            add(st, f"km_se{t}_{lab}", round(k[t][1], 6))
            add(st, f"na_h{t}_{lab}", round(h[t], 6))
        add(st, f"rmst730_{lab}", round(rmst(sub.time, sub.event, 730), 4))
    add(st, "rmst730_diff_chf1_minus_chf0",
        round(rmst(prim.time[prim.chf == 1], prim.event[prim.chf == 1], 730)
              - rmst(prim.time[prim.chf == 0], prim.event[prim.chf == 0], 730), 4))
    for w in ("logrank", "wilcoxon", "tarone", "peto"):
        chi, pv = weighted_logrank(prim.time, prim.event, prim.chf, w)
        add(st, f"test_{w}_chisq", round(chi, 4))
        add(st, f"test_{w}_p", round(pv, 4))
    add(st, "distinct_event_times_primary", int(prim.time[prim.event == 1].nunique()))

    # ---------------- Stage 4: 06_cox (seminar section 5) -----------------------
    st = "06_cox"
    from scipy.stats import norm
    z = norm.ppf(0.975)
    T, E = prim.time.values, prim.event.values
    d = prim.assign(race_black=(prim.race == 2).astype(int), race_other=(prim.race == 3).astype(int),
                    race_hisp=(prim.race == 5).astype(int), agec=prim.age - 75)
    adj_cols = ["chf", "age", "female", "race_black", "race_other", "race_hisp", "esrd", "alzh", "ckd",
                "cancer", "copd", "depress", "diab", "ihd", "osteo", "raoa", "stroke", "logreimb"]

    def add_model(tag, cols, b, se, ll, keep):
        add(st, f"{tag}_m2ll", round(-2 * ll, 4))
        for c in keep:
            k = cols.index(c)
            add(st, f"{tag}_{c}_beta", round(b[k], 6))
            add(st, f"{tag}_{c}_se", round(se[k], 6))
            add(st, f"{tag}_{c}_hr", round(np.exp(b[k]), 6))
            add(st, f"{tag}_{c}_hrlcl", round(np.exp(b[k] - z * se[k]), 6))
            add(st, f"{tag}_{c}_hrucl", round(np.exp(b[k] + z * se[k]), 6))

    b, se, ll = cox(d[["chf"]].values, T, E, "efron")
    add_model("cox_unadj", ["chf"], b, se, ll, ["chf"])
    b, se, ll = cox(d[["chf"]].values, T, E, "breslow")
    add_model("cox_unadj_breslow", ["chf"], b, se, ll, ["chf"])
    b, se, ll = cox(d[adj_cols].values, T, E, "efron")
    add_model("cox_adj", adj_cols, b, se, ll, adj_cols)

    int_cols = ["chf", "agec"] + adj_cols[2:] + ["chf_agec", "chf_female"]
    di = d.assign(chf_agec=d.chf * d.agec, chf_female=d.chf * d.female)
    b, se, ll = cox(di[int_cols].values, T, E, "efron")
    add_model("cox_int", int_cols, b, se, ll, ["chf", "chf_agec", "chf_female"])
    # Hazard ratio for CHF at ages 70, 80, 90 for men and women, with Wald limits.
    bi = dict(zip(int_cols, b))
    for age_ in (70, 80, 90):
        for fem in (0, 1):
            lin = bi["chf"] + bi["chf_agec"] * (age_ - 75) + bi["chf_female"] * fem
            add(st, f"cox_int_hr_chf_age{age_}_female{fem}", round(np.exp(lin), 6))

    dt = d.assign(chf_lt=d.chf)
    b, se, ll = cox_tv(dt[["chf", "chf_lt"]].values, T, E, [(1, np.log)], "efron")
    add_model("cox_tv", ["chf", "chf_lt"], b, se, ll, ["chf", "chf_lt"])

    # ---------------- Stage 5: 06b_plasmode ----------------------------------
    # Every random draw is hashed from the beneficiary ID (see hash_u), so the
    # simulated data are identical to the SAS data set and every fitted value
    # can be checked.
    st = "06b_plasmode"
    from scipy.optimize import brentq
    pc = plasmode_config()
    seed = read_seed()
    agec = prim.age.values - 75
    eta0 = (pc["plas_chf_a"] * prim.chf.values + pc["plas_age1"] * agec + pc["plas_age2"] * agec**2
            + np.log(pc["plas_hr_diab"]) * prim.diab.values + np.log(pc["plas_hr_female"]) * prim.female.values)
    shape = pc["plas_shape"] + pc["plas_chf_b"] * prim.chf.values
    h730 = pc["plas_shape"] * np.exp(eta0) * ADMIN_DAYS**shape / shape
    loglam = brentq(lambda L: np.mean(1 - np.exp(-np.exp(L) * h730)) - pc["plas_target"], -20, 0, xtol=1e-15)
    pm = 1 - np.exp(-np.exp(loglam) * h730)
    add(st, "plas_n", len(prim))
    add(st, "plas_loglambda", loglam)
    add(st, "plas_expected_mort_all", pm.mean())
    add(st, "plas_expected_mort_chf1", pm[prim.chf.values == 1].mean())
    add(st, "plas_expected_mort_chf0", pm[prim.chf.values == 0].mean())

    # Self-test values: the five smallest IDs, draws 0, 1, 2.
    ids = np.asarray(prim.index)
    for i, pid in enumerate(sorted(ids)[:5], 1):
        for k in range(3):
            h = hash_h([pid], seed, k)[0]
            add(st, f"plas_hash_h_i{i}_k{k}", str(int(h)))
            add(st, f"plas_hash_u_i{i}_k{k}", repr(float((int(h) + 0.5) / 2**52)))

    pl = prim.assign(agec=agec, agec2=agec**2, race_black=(prim.race == 2).astype(int),
                     race_other=(prim.race == 3).astype(int), race_hisp=(prim.race == 5).astype(int))
    gen = lambda idx, ku, kc: plasmode_generate(ids[idx], eta0[idx], shape[idx], loglam, pc, seed, ku, kc)
    allidx = np.arange(len(prim))
    tday, tmon, ev, ctype = gen(allidx, 1, 2)
    add(st, "plas_deaths", int(ev.sum()))
    add(st, "plas_random_censored", int((ctype == 1).sum()))
    add(st, "plas_admin_censored", int((ctype == 2).sum()))
    add(st, "plas_event_times_monthly", len(np.unique(tmon[ev == 1])))
    add(st, "plas_event_times_daily", len(np.unique(tday[ev == 1])))
    add(st, "plas_sum_time", int(tmon.sum()))
    add(st, "plas_sum_time_day", int(tday.sum()))

    pl_cols = ["chf", "chf_lt", "agec", "agec2", "female", "race_black", "race_other", "race_hisp",
               "esrd", "alzh", "ckd", "cancer", "copd", "depress", "diab", "ihd", "osteo", "raoa",
               "stroke", "logreimb"]
    Xp = pl.assign(chf_lt=pl.chf)[pl_cols].values.astype(float)
    specs = {"DE": (tday, "efron", 0), "ME": (tmon, "efron", 0),
             "MB": (tmon, "breslow", 0), "MM": (tmon, "efron", 15)}
    for spec, (tt, ties, off) in specs.items():
        b, se, ll = cox_tv(Xp, tt, ev, [(1, lambda t, o=off: np.log(t + o))], ties)
        add(st, f"plas_{spec}_m2ll", round(-2 * ll, 4))
        for j, c in enumerate(pl_cols):
            add(st, f"plas_{spec}_{c}_beta", round(b[j], 6))
            add(st, f"plas_{spec}_{c}_se", round(se[j], 6))
        print(f"plasmode fit {spec} done", file=sys.stderr, flush=True)

    # Replicate study on the fixed subsample (smallest draw-0 uniforms).
    nrep, nsub = int(pc["plas_reps"]), int(pc["plas_subn"])
    if nrep > 0:
        u0 = (hash_h(ids, seed, 0) + 0.5) / 2**52
        sub = np.argsort(u0, kind="stable")[:nsub]
        truth = {"chf": pc["plas_chf_a"], "chf_lt": pc["plas_chf_b"], "agec": pc["plas_age1"],
                 "agec2": pc["plas_age2"], "female": np.log(pc["plas_hr_female"]),
                 "diab": np.log(pc["plas_hr_diab"])}
        z = 1.959963984540054
        est = {sp: [] for sp in ("ME", "MB", "MM")}
        for r in range(1, nrep + 1):
            _, tm_r, ev_r, _ = gen(sub, 2 * r + 1, 2 * r + 2)
            for sp, ties, off in (("ME", "efron", 0), ("MB", "breslow", 0), ("MM", "efron", 15)):
                b, se, _ = cox_tv(Xp[sub], tm_r, ev_r, [(1, lambda t, o=off: np.log(t + o))], ties)
                est[sp].append((b, se))
        for sp, lst in est.items():
            B = np.array([x[0] for x in lst]); S = np.array([x[1] for x in lst])
            for c, tv in truth.items():
                j = pl_cols.index(c)
                cov = np.mean((B[:, j] - z * S[:, j] <= tv) & (tv <= B[:, j] + z * S[:, j]))
                add(st, f"plas_rep_{sp}_{c}_mean", round(B[:, j].mean(), 6))
                add(st, f"plas_rep_{sp}_{c}_sd", round(B[:, j].std(ddof=1), 6))
                add(st, f"plas_rep_{sp}_{c}_se", round(S[:, j].mean(), 6))
                add(st, f"plas_rep_{sp}_{c}_cov", cov)
        print("replicates done", file=sys.stderr, flush=True)

    # ---------------- Stage 5: 07_functional_form and 08_ph_assumption ---------
    # Run on the real cohort and on the plasmode data (metric prefix plas_).
    knots = [sas_pctl(prim.age, q) for q in (5, 27.5, 50, 72.5, 95)]
    for i, k in enumerate(knots, 1):
        add("07_functional_form", f"ff_knot{i}", k)
    sb = rcs_basis(d.age.values, knots)
    ds = d.assign(s1=sb[:, 0], s2=sb[:, 1], s3=sb[:, 2])
    spl_cols = adj_cols + ["s1", "s2", "s3"]
    nochf = [c for c in adj_cols if c != "chf"]
    dp = d.assign(p1=d.chf, p2=d.chf, p3=d.chf)
    for pre, TT, EE in (("", T, E), ("plas_", tmon, ev)):
        st = "07_functional_form"
        _, _, ll_lin = cox(d[adj_cols].values, TT, EE, "efron")
        b, se, ll_spl = cox(ds[spl_cols].values, TT, EE, "efron")
        add(st, f"{pre}ff_linear_m2ll", round(-2 * ll_lin, 4))
        add(st, f"{pre}ff_spline_m2ll", round(-2 * ll_spl, 4))
        add(st, f"{pre}ff_lr_nonlin_chisq", round(2 * (ll_spl - ll_lin), 4))
        for c in ("age", "s1", "s2", "s3"):
            k = spl_cols.index(c)
            add(st, f"{pre}ff_spline_{c}_beta", round(b[k], 6))
            add(st, f"{pre}ff_spline_{c}_se", round(se[k], 6))

        st = "08_ph_assumption"
        for v in ("chf", "age", "female", "diab"):
            cols = adj_cols + [v]
            b, se, _ = cox_tv(d[cols].values, TT, EE, [(len(cols) - 1, np.log)], "efron")
            add(st, f"{pre}ph_tv_{v}_beta", round(b[-1], 6))
            add(st, f"{pre}ph_tv_{v}_se", round(se[-1], 6))
        cols = ["p1", "p2", "p3"] + nochf
        b, se, _ = cox_tv(dp[cols].values, TT, EE,
                          [(0, lambda t: float(t <= 181)), (1, lambda t: float(181 < t <= 365)),
                           (2, lambda t: float(t > 365))], "efron")
        for i in range(3):
            add(st, f"{pre}ph_piece_p{i + 1}_beta", round(b[i], 6))
            add(st, f"{pre}ph_piece_p{i + 1}_se", round(se[i], 6))

    # Plasmode only: the correctly specified model (true age and CHF forms) with
    # each covariate by log(time + 15) added in turn.
    st = "08_ph_assumption"
    f15 = lambda t: np.log(t + 15)
    for v in ("agec", "female", "diab"):
        cols = pl_cols + [v]
        Xc = pl.assign(chf_lt=pl.chf)[cols].values.astype(float)
        b, se, _ = cox_tv(Xc, tmon, ev, [(1, f15), (len(cols) - 1, f15)], "efron")
        add(st, f"plas_ph_tvc_{v}_beta", round(b[-1], 6))
        add(st, f"plas_ph_tvc_{v}_se", round(se[-1], 6))

    # ---------------- Stage 5: 09_influence (seminar section 8) ----------------
    st = "09_influence"
    Xa = d[adj_cols].values.astype(float)
    for pre, TT, EE in (("", T, E), ("plas_", tmon, ev)):
        b, se, dfb, ld = cox_influence(Xa, TT, EE)
        add(st, f"inf_{pre}max_ld", round(ld.max(), 8))
        add(st, f"inf_{pre}beta_chf_all", round(b[0], 6))
        top = np.argsort(-ld, kind="stable")[:20]
        keep = np.ones(len(TT), bool); keep[top] = False
        b2, _, _ = cox(Xa[keep], TT[keep], EE[keep], "efron")
        add(st, f"inf_{pre}beta_chf_without_top20", round(b2[0], 6))
        for c in ("chf", "age", "female", "diab"):
            k = adj_cols.index(c)
            add(st, f"inf_{pre}n_over_0_1_{c}", int((np.abs(dfb[:, k] / se[k]) > 0.1).sum()))

    # ---------------- Stage 6: 10_sensitivity ----------------------------------
    # The adjusted model of 06 (race as indicators, Efron ties) in each design.
    st = "10_sensitivity"
    x8 = x.assign(age_int08=age_int(x.birth_dt, pd.Timestamp("2008-01-01")),
                  age08=(pd.Timestamp("2008-01-01") - x.birth_dt).dt.days / 365.25)
    nv = x8[(x8.age_int08 >= 65) & (x8.ffs08 == 1)].copy()
    nv_ev = nv.death_dt.between(pd.Timestamp("2008-01-01"), pd.Timestamp("2010-12-31"))
    nv["event"] = nv_ev.astype(int)
    nv["time"] = np.where(nv_ev, (nv.death_dt - pd.Timestamp("2007-12-31")).dt.days, NAIVE_ADMIN).astype(int)
    nv["age"] = nv.age08
    add(st, "sens_naive_ltfu", int((~nv.in2009 & nv.death_dt.isna()).sum()
                                  + (nv.in2009 & ~nv.in2010 & nv.death_dt.isna()).sum()))
    covc = prim[~((prim.time_cov == 0) & (prim.event_cov == 0))].assign(time=lambda q: q.time_cov,
                                                                       event=lambda q: q.event_cov)
    noreimb = [c for c in adj_cols if c != "logreimb"]
    designs = [("primary", prim, adj_cols, []), ("naive_unadj", nv, ["chf"], []),
               ("naive_noreimb", nv, noreimb, []), ("naive", nv, adj_cols, []),
               ("cov", covc, adj_cols, []), ("full65", full65, adj_cols, []),
               ("allages", alive09.assign(under65=(alive09.age_int < 65).astype(int)), adj_cols, ["under65"])]
    for tag, df, base, extra in designs:
        dd = df.assign(race_black=(df.race == 2).astype(int), race_other=(df.race == 3).astype(int),
                       race_hisp=(df.race == 5).astype(int))
        cols = base + extra
        b, se, ll = cox(dd[cols].values.astype(float), dd.time.values, dd.event.values, "efron")
        add(st, f"sens_{tag}_n", len(dd))
        add(st, f"sens_{tag}_events", int(dd.event.sum()))
        add(st, f"sens_{tag}_m2ll", round(-2 * ll, 4))
        add(st, f"sens_{tag}_chf_beta", round(b[0], 6))
        add(st, f"sens_{tag}_chf_se", round(se[0], 6))
        add(st, f"sens_{tag}_chf_hr", round(np.exp(b[0]), 6))
        add(st, f"sens_{tag}_chf_hrlcl", round(np.exp(b[0] - z * se[0]), 6))
        add(st, f"sens_{tag}_chf_hrucl", round(np.exp(b[0] + z * se[0]), 6))
        if extra:
            add(st, f"sens_{tag}_under65_beta", round(b[-1], 6))
            add(st, f"sens_{tag}_under65_se", round(se[-1], 6))
        print(tag, len(dd), int(dd.event.sum()), round(np.exp(b[0]), 4),
              round(np.exp(b[0] - z * se[0]), 4), round(np.exp(b[0] + z * se[0]), 4), file=sys.stderr, flush=True)

    # Why the naive adjusted model moves: 2008 reimbursement by death in 2008.
    died08 = ((nv.event == 1) & (nv.time <= 366)).astype(int)
    for g, m in (("all", np.ones(len(nv), bool)), ("chf1", nv.chf.values == 1), ("chf0", nv.chf.values == 0)):
        for dv in (0, 1):
            sel = m & (died08.values == dv)
            add(st, f"sens_naive_mlogreimb_{g}_died{dv}", round(nv.logreimb.values[sel].mean(), 6))
    for cv in (0, 1):
        add(st, f"sens_naive_died08_rate_chf{cv}", round(died08[nv.chf == cv].mean(), 6))
    add(st, "sens_naive_corr_chf_logreimb", round(np.corrcoef(nv.chf, nv.logreimb)[0, 1], 6))

    out = pd.DataFrame(vals)
    # Strings are already formatted at full precision (hash self-test values).
    out["value"] = out.value.map(lambda v: v if isinstance(v, str) else f"{float(v):.10g}")
    OUT.parent.mkdir(exist_ok=True)
    out.to_csv(OUT, index=False, lineterminator="\n")
    att.to_csv(REPO / "validation" / "attrition_python.csv", index=False, lineterminator="\n")
    print(att.to_string(index=False))
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
