#!/usr/bin/python3
import argparse
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from joblib import delayed,Parallel
import os

def whichKeep(est_params):
    kon = np.array(est_params)[:,0]
    koff = np.array(est_params)[:,1]
    ksyn = np.array(est_params)[:,2]
    which_kon = ~(kon < 2*1e-3)*~(kon > 1e3 - 1)
    which_koff = ~(koff < 2*1e-3)*~(koff > 1e3 - 1)
    which_burst = ksyn/koff > 1
    which_ksyn = ksyn > 1
    which = which_burst*which_koff*which_kon*which_ksyn
    return which


def MaximumLikelihood(vals, export_asymp_ci = False, fix = 0, metod = 'L-BFGS-B'):
    from scipy.interpolate import interp1d
    from scipy.optimize import minimize
    from scipy import special
    from scipy.stats import poisson,norm
    from scipy.special import j_roots
    from scipy.special import beta as beta_fun
    import numpy as np
    if len(vals) == 0:
        return np.array([np.nan, np.nan, np.nan])
    def dBP(at, alpha, bet, lam):
        at.shape = (len(at), 1)
        np.repeat(at, 50, axis = 1)
        def fun(at, m):
            if(max(m) < 1e6):
                return(poisson.pmf(at,m))
            else:
                return(norm.pdf(at,loc=m,scale=sqrt(m)))

        x,w = j_roots(50,alpha = bet - 1, beta = alpha - 1)
        gs = np.sum(w*fun(at, m = lam*(1+x)/2), axis=1)
        prob = 1/beta_fun(alpha, bet)*2**(-alpha-bet+1)*gs
        return(prob)
    def LogLikelihood(x, vals):
        kon = x[0]
        koff = x[1]
        ksyn = x[2]
        return(-np.sum(np.log( dBP(vals,kon,koff,ksyn) + 1e-10) ) )
    x0 = MomentInference(vals)
    if np.isnan(x0).any() or any(x0 < 0):
        x0 = np.array([10,10,10])
    bnds = ((1e-3,1e3),(1e-3,1e3), (1, 1e4))
    vals_ = np.copy(vals) # Otherwise the structure is violated.
    try:
        ll = minimize(LogLikelihood, x0, args = (vals_), method=metod, bounds=bnds)
    except:
        return np.array([np.nan,np.nan,np.nan])
    #se = ll.hess_inv.todense().diagonal()
    estim = ll.x
    return estim

# moment-based inference
def MomentInference(vals, export_moments=False):
    # code from Anton Larsson's R implementation
    from scipy import stats # needs imports inside function when run in ipyparallel
    import numpy as np
    m1 = float(np.mean(vals))
    m2 = float(sum(vals*(vals - 1))/len(vals))
    m3 = float(sum(vals*(vals - 1)*(vals - 2))/len(vals))

    # sanity check on input (e.g. need at least on expression level)
    if sum(vals) == 0: return np.nan
    if m1 == 0: return np.nan
    if m2 == 0: return np.nan

    r1=m1
    r2=m2/m1
    r3=m3/m2

    if (r1*r2-2*r1*r3 + r2*r3) == 0: return np.nan
    if ((r1*r2 - 2*r1*r3 + r2*r3)*(r1-2*r2+r3)) == 0: return np.nan
    if (r1 - 2*r2 + r3) == 0: return np.nan

    lambda_est = (2*r1*(r3-r2))/(r1*r2-2*r1*r3 + r2*r3)
    mu_est = (2*(r3-r2)*(r1-r3)*(r2-r1))/((r1*r2 - 2*r1*r3 + r2*r3)*(r1-2*r2+r3))
    v_est = (2*r1*r3 - r1*r2 - r2*r3)/(r1 - 2*r2 + r3)

    if export_moments:
        return np.array([lambda_est, mu_est, v_est, r1, r2, r3])

    return np.array([lambda_est, mu_est, v_est])


parser = argparse.ArgumentParser(description='Maximum likelihood inference of bursting kinetics from scRNA-seq data')
parser.add_argument('file', metavar='file', type=str, nargs=1,help='.csv file with allelic-resolution transcript counts' )
parser.add_argument('--njobs', default=[50], nargs=1, type=int, help='Number of jobs for the parallelization, default 50')
args = parser.parse_args()
filename = args.file[0]
njobs = args.njobs[0]
print('Reading file ' + filename)
rpkm = pd.read_csv(filename, index_col=0)

print('Inferring kinetics:')
params = Parallel(n_jobs=njobs, verbose = 3)(delayed(MaximumLikelihood)(np.around(rpkm[pd.notnull(rpkm)])) for i,rpkm in rpkm.iterrows())
keep = whichKeep(params)

print('Inferred kinetics of {} genes out of {} total'.format(np.sum(keep), len(keep)))

base = os.path.splitext(os.path.basename(filename))[0]
base = base + '_ML.pkl'
print('Saving result to ' + base)

pd.to_pickle(pd.DataFrame([ params, list(keep)], columns=rpkm.index).T, base)
