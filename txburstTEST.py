#!/usr/bin/python3
import argparse
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from joblib import delayed,Parallel
import os

def fold_change_bf_and_bs(vals_1,vals_2, start_estim_1, start_estim_2, name):
    from scipy.optimize import minimize
    from scipy import special
    from scipy.stats import poisson,norm, chi2
    from scipy.special import j_roots
    from scipy.special import beta as beta_fun
    import numpy as np
    try:
        def dBP(at, alpha, bet, lam):
            at.shape = (len(at), 1)
            np.repeat(at, 40, axis = 1)
            def fun(at, m):
                if(max(m) < 1e6):
                    return(poisson.pmf(at,m))
                else:
                    return(norm.pdf(at,loc=m,scale=np.sqrt(m)))

            x,w = j_roots(40,alpha = bet - 1, beta = alpha - 1)
            gs = np.sum(w*fun(at, m = lam*(1+x)/2), axis=1)
            prob = 1/beta_fun(alpha, bet)*2**(-alpha-bet+1)*gs
            return(prob)
        def minll_bs(x, theta, obs_1, obs_2):
            theta = 2**theta
            l_1 = x[0]
            mu_1 = x[1]
            v_1 = x[2]
            l_2 = x[3]
            mu_2 = x[4]
            v_2 = (v_1*mu_2)/(mu_1*theta)
            return(- (np.sum(np.log(dBP(obs_1,l_1,mu_1,v_1) + 1e-10)) + np.sum(np.log(dBP(obs_2,l_2,mu_2,v_2) + 1e-10))) )
        def minll_bf(x, theta, obs_1, obs_2):
            theta = 2**theta
            l_1 = x[0]
            mu_1 = x[1]
            v_1 = x[2]
            l_2 = theta*l_1
            mu_2 = x[3]
            v_2 = x[4]
            return(- (np.sum(np.log(dBP(obs_1,l_1,mu_1,v_1) + 1e-10)) + np.sum(np.log(dBP(obs_2,l_2,mu_2,v_2) + 1e-10))) )


        bnds_bs = ((1e-4,1e4), (1e-4,1e4), (1e-4,1e4), (1e-4,1e4), (1e-4,1e4))
        ll_null_bs = minimize(minll_bs, x0 = [start_estim_1[0], start_estim_1[1],start_estim_1[2], start_estim_2[0], start_estim_2[1]], args = (0,vals_1,vals_2), method='TNC', bounds = bnds_bs).fun
        ll_alt_bs = minimize(minll_bs, x0 = [start_estim_1[0], start_estim_1[1],start_estim_1[2], start_estim_2[0], start_estim_2[1]], args = (np.log2((start_estim_1[2]*start_estim_2[1])/(start_estim_2[2]*start_estim_1[1])),vals_1,vals_2), method='TNC', bounds = bnds_bs).fun
        p_bs = 1-chi2.cdf(2*(ll_null_bs - ll_alt_bs),1)

        bnds_bf = ((1e-4,1e4), (1e-4,1e4), (1e-4,1e4), (1e-4,1e4), (1e-4,1e4))
        ll_null_bf = minimize(minll_bf, x0 = [start_estim_1[0], start_estim_1[1],start_estim_1[2], start_estim_2[1], start_estim_2[2]], args = (0,vals_1,vals_2), method='TNC', bounds = bnds_bf).fun
        ll_alt_bf = minimize(minll_bf, x0 = [start_estim_1[0], start_estim_1[1],start_estim_1[2], start_estim_2[1], start_estim_2[2]], args = (np.log2(start_estim_2[0]/start_estim_1[0]),vals_1,vals_2), method='TNC', bounds = bnds_bf).fun
        p_bf = 1-chi2.cdf(2*(ll_null_bf - ll_alt_bf),1)
    except ValueError:
        return start_estim_1, start_estim_2, np.nan, np.nan
    return start_estim_1, start_estim_2, p_bf, p_bs, name

parser = argparse.ArgumentParser(description='Hypothesis testing of differences in bursting kinetics from scRNA-seq data')
parser.add_argument("--file1", metavar='file1',type=str, nargs=1,help='.csv file 1 with allelic-resolution transcript counts' )
parser.add_argument("--file2", metavar='file2', type=str, nargs=1,help='.csv file 2 with allelic-resolution transcript counts' )
parser.add_argument("--njobs", default=[50], nargs=1, type=int, help='Number of jobs for the parallelization, default 50')
parser.add_argument("--ML1", default=None, nargs=1, type=str, help='Maximum Likelihood file 1 (required)')
parser.add_argument("--ML2", default=None, nargs=1, type=str, help='Maximum Likelihood file 2 (required)')

args = parser.parse_args()
filename1 = args.file1[0]
filename2 = args.file2[0]
njobs = args.njobs[0]
ML1 = args.ML1[0]
ML2 = args.ML2[0]

MLPickle1 = pd.read_pickle(ML1)
MLPickle2 = pd.read_pickle(ML2)
print('Reading in {}'.format(filename1))
rpkm1 = pd.read_csv(filename1, index_col=0)
print('Reading in {}'.format(filename2))
rpkm2 = pd.read_csv(filename2, index_col=0)

#MLPickle1 = MLPickle1[MLPickle1[1]][0]
#MLPickle2 = MLPickle2[MLPickle2[1]][0]

MLPickle1 = MLPickle1[0]
MLPickle2 = MLPickle2[0]

MLPickle1 = MLPickle1[~MLPickle1.index.duplicated(keep='first')]
MLPickle2 = MLPickle2[~MLPickle2.index.duplicated(keep='first')]

df = pd.DataFrame([MLPickle1, MLPickle2], index=['ML1', 'ML2']).T.dropna(how='any')

print('Calculating fold changes and hypothesis testing of {} genes: '.format(len(df.index.values)))

fold_changes = Parallel(n_jobs=njobs, verbose = 3)(delayed(fold_change_bf_and_bs)(np.array(np.around(rpkm1.loc[gene][pd.notnull(rpkm1.loc[gene])])), np.array(np.around(rpkm2.loc[gene][pd.notnull(rpkm2.loc[gene])])), df['ML1'][gene], df['ML2'][gene], gene) for gene in df.index.values)

base1 = os.path.splitext(os.path.basename(filename1))[0]
base2 = os.path.splitext(os.path.basename(filename2))[0]
base = base1 + '_vs_' + base2
base = base + '_TEST.pkl'

print('Saving result to: {}'.format(base))

pd.to_pickle(pd.DataFrame(fold_changes, index=df.index), base)
