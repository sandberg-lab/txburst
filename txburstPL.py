#!/usr/bin/python3
from scipy import special
from scipy.stats import poisson,norm, chi2
from scipy.special import j_roots
from scipy.special import beta as beta_fun
from scipy.optimize import minimize
import argparse
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
from joblib import delayed,Parallel
import os
from scipy.interpolate import interp1d
from scipy.optimize import minimize, fsolve
from scipy import special
from scipy.stats import poisson,norm, chi2
from scipy.special import j_roots
from scipy.special import beta as beta_fun
import numpy as np
def dBP(at, alpha, bet, lam):
    at.shape = (len(at), 1)
    np.repeat(at, 50, axis = 1)
    def fun(at, m):
        if(max(m) < 1e6):
            return(poisson.pmf(at,m))
        else:
            return(norm.pdf(at,loc=m,scale=np.sqrt(m)))
    x,w = j_roots(50,alpha = bet - 1, beta = alpha - 1)
    gs = np.sum(w*fun(at, m = lam*(1+x)/2), axis=1)
    prob = 1/beta_fun(alpha, bet)*2**(-alpha-bet+1)*gs
    return(prob)
def ll_kon(x, l, obs):
    mu = x[0]
    v = x[1]
    return(-np.sum(np.log(dBP(obs,l,mu,v) + 1e-10)))
def ll_burst_size(x, burst_size, obs):
    l = x[0]
    v = x[1]
    mu = v/burst_size
    return(-np.sum(np.log(dBP(obs,l,mu,v) + 1e-10)))
def minll_kon(x, kon, obs):
    bnds = ((1e-3, 1e3), (1,1e10))
    res = minimize(ll_kon, x, args = (kon,obs), method='L-BFGS-B', bounds = bnds)
    ll = res.fun
    x0 = res.x
    return ll, x0
def minll_burst_size(x, burst_size, obs):
    bnds = ((1e-3, 1e3), (1,1e6))
    res = minimize(ll_burst_size, x, args = (burst_size,obs), method='L-BFGS-B', bounds = bnds)
    ll = res.fun
    x0 = res.x
    return ll, x0
def get_h(param):
    x = np.log10(param)
    h = 10**(x-2)
    return(h)

def kon_ll(vals,start_estim):
    cutoff = chi2.ppf(1-alph, 1)/2
    x0 = [start_estim[1], start_estim[2]]
    bnds = ((1e-3, 1e3), (1,1e10))
    N = 100
    start = start_estim[0]
    res = minimize(ll_kon, x0, args = (start,vals), method='L-BFGS-B', bounds = bnds)
    h = get_h(start)
    h = 5*h
    ll_p = res.fun
    ll_l = np.zeros(N)
    kon_l = np.array([])
    try:
        for i in range(N):
            kon = start - h*i
            if kon <= 0:
                break
            kon_l = np.append(kon_l, kon)
            try:
                res = minimize(ll_kon, res.x, args = (kon,vals), method='L-BFGS-B', bounds = bnds)
            except ValueError:
                break
            ll_l[i] = res.fun
            if (2*(ll_l[i] - min(ll_l[ll_l > 0])) > cutoff + 0.5) and (ll_l[i] > ll_l[i-1]):
                break
    except ValueError:
        return np.array([minimum,np.nan,np.nan]), np.nan

    ll_l = ll_l[:i+1]
    ll_l = ll_l[::-1]
    kon_l = kon_l[::-1]
    res = minimize(ll_kon, x0, args = (start,vals), method='L-BFGS-B', bounds = bnds)
    ll_u = np.zeros(N)
    kon_u = np.array([])
    try:
        for j in range(N):
            kon = start + h*j
            kon_u = np.append(kon_u, kon)
            try:
                res = minimize(ll_kon, res.x, args = (kon,vals), method='L-BFGS-B', bounds = bnds)
            except ValueError:
                break
            ll_u[j] = res.fun
            if (2*(ll_u[j] - min(ll_u[ll_u > 0])) > cutoff + 0.5) and (ll_u[j] > ll_u[j-1]):
                break
    except ValueError:
        return np.array([minimum,np.nan,np.nan]), np.nan

    ll_u = ll_u[:j+1]

    ll = np.concatenate((ll_l[:-1],np.array([ll_p]),ll_u[1:]))
    kon_indexed = np.concatenate((kon_l[:-1],np.array([start]),kon_u[1:])).squeeze()

    ll_ratio = 2*(ll - min(ll)).squeeze()

    minimum_idx = np.argmin(ll_ratio)

    ll_right_side = ll_ratio[minimum_idx:]
    ll_left_side = ll_ratio[:minimum_idx]

    minimum = kon_indexed[minimum_idx]
    kon_right_side = kon_indexed[minimum_idx:]
    kon_left_side = kon_indexed[:minimum_idx]
    try:
        f_1 = interp1d(ll_left_side,kon_left_side, kind='cubic')
        f_2 = interp1d(ll_right_side,kon_right_side , kind='cubic')
        res = np.array([minimum, f_1(cutoff), f_2(cutoff)])
        return res, kon_indexed, ll_ratio
    except (ValueError,np.linalg.linalg.LinAlgError, TypeError):
        return np.array([minimum,np.nan,np.nan]), kon_indexed, ll_ratio

def burst_size_ll(vals, start_estim):
    cutoff = chi2.ppf(1-alph, 1)/2
    x0 = [start_estim[0], start_estim[2]]
    bnds = ((1e-3, 1e3), (1,1e10))
    N = 100
    start = start_estim[2]/start_estim[1]
    res = minimize(ll_burst_size, x0, args = (start,vals), method='L-BFGS-B', bounds = bnds)
    h = get_h(start)
    h = 3*h
    ll_p = res.fun
    ll_l = np.zeros(N)
    burst_size_l = np.array([])
    try:
        for i in range(N):
            burst_size = start - h*i
            if burst_size <= 0:
                break
            burst_size_l = np.append(burst_size_l, burst_size)
            try:
                res = minimize(ll_burst_size, res.x, args = (burst_size,vals), method='L-BFGS-B', bounds = bnds)
            except ValueError:
                break
            ll_l[i] = res.fun
            if (2*(ll_l[i] - min(ll_l[ll_l > 0])) > cutoff + 0.5) and (ll_l[i] > ll_l[i-1]):
                break
    except ValueError:
        return np.array([minimum,np.nan,np.nan]), np.nan
    ll_l = ll_l[:i+1]
    ll_l = ll_l[::-1]
    burst_size_l = burst_size_l[::-1]
    res = minimize(ll_burst_size, x0, args = (start,vals), method='L-BFGS-B', bounds = bnds)
    ll_u = np.zeros(N)
    burst_size_u = np.array([])
    try:
        for j in range(N):
            burst_size = start + h*j
            burst_size_u = np.append(burst_size_u, burst_size)
            try:
                res = minimize(ll_burst_size, res.x, args = (burst_size,vals), method='L-BFGS-B', bounds = bnds)
            except ValueError:
                break
            ll_u[j] = res.fun
            if (2*(ll_u[j] - min(ll_u[ll_u > 0])) > cutoff + 0.5) and (ll_u[j] > ll_u[j-1]):
                break
    except ValueError:
        return np.array([minimum,np.nan,np.nan]), np.nan

    ll_u = ll_u[:j+1]

    ll = np.concatenate((ll_l[:-1],np.array([ll_p]),ll_u[1:]))
    burst_size_indexed = np.concatenate((burst_size_l[:-1],np.array([start]),burst_size_u[1:])).squeeze()

    ll_ratio = 2*(ll - min(ll)).squeeze()

    minimum_idx = np.argmin(ll_ratio)

    ll_right_side = ll_ratio[minimum_idx:]
    ll_left_side = ll_ratio[:minimum_idx]

    minimum = burst_size_indexed[minimum_idx]
    burst_size_right_side = burst_size_indexed[minimum_idx:]
    burst_size_left_side = burst_size_indexed[:minimum_idx]
    try:
        f_1 = interp1d(ll_left_side,burst_size_left_side, kind='linear')
        f_2 = interp1d(ll_right_side,burst_size_right_side , kind='linear')
        res = np.array([minimum, f_1(cutoff), f_2(cutoff)])
        return res, burst_size_indexed, ll_ratio
    except (ValueError,np.linalg.linalg.LinAlgError, TypeError):
        return np.array([minimum,np.nan,np.nan]), burst_size_indexed, ll_ratio
def PL2(vals, start_estim):
    vals_ = np.copy(vals)

    res_kon, kon_idx, kon_ll_ratio = kon_ll(vals_, start_estim)
    res_burst_size, bs_idx, bs_ll_ratio = burst_size_ll(vals_,start_estim)

    return start_estim, res_kon, res_burst_size

parser = argparse.ArgumentParser(description='Confidence intervals on parameters of bursting kinetics from scRNA-seq data')
parser.add_argument('--file', metavar='file', type=str, nargs=1,help='.csv file file with allelic-resolution transcript counts' )
parser.add_argument('--njobs', default=[50], nargs=1, type=int, help='Number of jobs for the parallelization, default 50')
parser.add_argument('--alpha', default=[0.05], nargs=1, type=float, help='Alpha significance level (default 0.05)')
parser.add_argument("--MLFile", default=None, nargs=1, type=str, help='Maximum Likelihood file (required)')
args = parser.parse_args()
filename = args.file[0]
njobs = args.njobs[0]
ML = args.MLFile[0]
alph = args.alpha[0]
print('Reading file ' + filename)
counts = pd.read_csv(filename, index_col=0)

print('Reading file ' + ML)
MLPickle = pd.read_pickle(ML)

MLPickle = MLPickle[MLPickle[1]][0]

MLPickle = MLPickle[~MLPickle.index.duplicated(keep='first')]

genes = np.intersect1d(MLPickle.index.values, counts.index.values)

PL_CIs = Parallel(n_jobs=njobs, verbose = 3)(delayed(PL2)(np.around(counts.loc[gene][pd.notnull(counts.loc[gene])]), MLPickle[gene]) for gene in genes)

base1 = os.path.splitext(os.path.basename(filename))[0]

base = base1 + '_PL.pkl'

print('Saving result to: {}'.format(base))

pd.to_pickle(pd.DataFrame(PL_CIs, index=counts.loc[genes].index), base)
