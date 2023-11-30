import os
import sys
import numpy as np

path_lib = os.path.join( '/',*os.getcwd().split('/')[:-1], 'simulator', 'lib' )
if path_lib not in sys.path:
    sys.path.append( path_lib )

from simulate_experiment import simulate_experiment_multispecies_newintro


def log_likelihood_multispecies_SEEIR( data, Nsim, seed, beta, sigmas, sigmaExp, mu, pClusterInterv, pIntroInterv, pIntroControl, pIntroBulk, Nms, pSurv,
 checkPts=[12,36,84], dtBurnin = 15 * 24, dt01 = 60, Nc=5, Ni=5, cost = 10. ):

    '''
    Computes the log-likelihood function for a single type of chickens.
    
    In the main analysis, this function should be used once for broiler and once for backyard chickens
    
    Uses an external simulator to compute the probability of a chicken testing positive
    at different stages of the experiment
    '''

    # get data about +ve and -ve chickens

    n_pos  = data['npos']
    n_susc = data['nsusc']

    # make sure model is correct

    nEstages = 2
    nIstages = 1

    # simulate process
    dT = len( pSurv )
    n_pos_sim = simulate_experiment_multispecies_newintro( Nsim, seed, checkPts, dtBurnin, dt01, Nc, Ni, pSurv, Nms, beta, sigmas, sigmaExp, mu, nEstages, nIstages, 
    pClusterInterv, pIntroInterv, pIntroControl, pIntroBulk )

    n_pos_sim = np.array( n_pos_sim )

    # compute probabilities control/intervention group

    # for positive cases

    pc = np.mean( n_pos_sim[0], axis = 0 ) / Nc
    pi = np.mean( n_pos_sim[1], axis = 0 ) / Ni

    cumpc = np.cumsum( pc )
    cumpi = np.cumsum( pi )

    # for negative/dead cases

    tmp = np.concatenate( [ [0.], cumpc ] )
    psc = 1. - tmp 

    tmp = np.concatenate( [ [0.], cumpi ] )
    psi = 1. - tmp

    psc[psc < 0.] = 0.
    psi[psi < 0.] = 0. 

    # compute log-likelihood

    logl = 0.

    # (main) Individual chicken state likelihood
    
    for group in ['0','1']:
        if group == '0':
            p = pc
            ps = psc
        else:
            p = pi
            ps = psi
        
        # infected chicken likelihood

        for i, n in enumerate( n_pos[group] ):
            if p[i] == 0.:
                logl -= n * cost
            else:
                logl += n * np.log( p[i] )

        # susceptible/dead chicken likelihood

        for i, n in enumerate( n_susc[group] ):
            if ps[i] == 0.:
                logl -= n * cost
            else:
                logl += n * np.log( ps[i] )

    if ( np.isnan( logl ) ):
        return -np.infty
    
    return logl