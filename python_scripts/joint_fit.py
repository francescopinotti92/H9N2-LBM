import os
import sys
import json
import time
import numpy as np
import emcee
from inference import log_likelihood_multispecies_SEEIR
from multiprocessing import Pool, Value, Manager
from scipy.special import gammaincc
import utils

def init_globals(counter):
    # This allows to set a global counter to be shared between processes                                                                                     
    global rng_gen
    rng_gen = counter

def compute_prob_SEEIRR_gamma( lam, k, sigma, mu, u, p ):
    '''
    Computes introduction probability for a SEEIRR model.
    Uses a gamma prior on time since exposure.
    
    lam : inverse scale gamma prior
    k   : shape gamma prior
    sigma: incubation rate
    mu : recovery rate
    u : rate from R1 -> R2 (R2 chickens are considered -ve)
    p : probability of a chickens being exposed
    '''    

    epsi = sigma * 2 # rate from E1 -> E2 & from E2 -> I 
    
    res = np.array( [ 0., 0., 0., 0., 0., 0. ] )
    
    res[0] = 1. - p
    res[1] = p * np.power( lam / ( lam + epsi ), k )
    res[2] = p * k * epsi * np.power( lam, k ) / np.power( lam + epsi, k + 1 )
    res[3] = p * np.square( epsi ) * np.power( lam, k ) * ( np.power( lam + mu, -k ) - np.power( lam + epsi, -k ) * ( 1 - k * ( mu - epsi ) / ( lam + epsi ) ) ) / np.square( mu - epsi )

    tmp = 0.
   
    tmp += np.power( lam + u, -k ) * np.square( ( epsi - mu ) / ( epsi - u ) ) / ( mu - u )
    tmp += -np.power( lam + mu, -k ) / ( mu - u )
    tmp += np.power( lam + epsi, -k ) * ( 2 * epsi - u - mu ) / np.square( epsi - u )
    tmp += k * np.power( lam + epsi, -k-1 ) * ( epsi - mu ) / ( epsi - u )
    
    res[4] = p * mu * np.square( epsi ) * np.power( lam, k ) * tmp / np.square( mu - epsi )
    res[5] = 1. - res[:5].sum()
    
    # check all elems pos and sum to 1
    
    res[res < 0.] = 0.
    res = res / res.sum()
    
    return res


def mean_exp( a, lam, k, dt = 0. ):
    '''
    Computes the average of exp(-ax) over f(x+dt|lam,k)/P(x>dt|lam,k),
    where f is the pdf of a gamma distribution with parameters lam and k
    '''
    return np.exp( a * dt ) * np.power( lam / ( lam + a ), k ) * gammaincc( k, ( lam + a ) * dt ) / gammaincc( k, lam * dt ) 

def mean_xexp( a, lam, k, dt = 0. ):
    '''
    Computes the average of x*exp(-ax) over f(x+dt|lam,k)/P(x>dt|lam,k),
    where f is the pdf of a gamma distribution with parameters lam and k
    '''
    
    A = np.exp( a * dt ) * np.power( lam / ( lam + a ), k ) / gammaincc( k, lam * dt )
    B = ( k / ( lam + a ) ) * gammaincc( k + 1, ( lam + a ) * dt ) - dt * gammaincc( k, ( lam + a ) * dt )
    
    return A * B

def compute_prob_SEEIRR_gamma_dt( lam, k, sigma, mu, u, p, dt = 0.):
    
    '''
    Computes the probability P(X) that a chicken is in compartment X=S,E1,E2,I,R1 or R2
    when it enters the market.
    
    p is the probability of prior exposure, hence P(S)=1-p
    in contrast P(E1)+P(E2)+P(I)+P(R1)+P(R2)=1-P(S)=p
    
    sigma is the inverse of the incubation period
    
    mu is the inverse of the infectious period
    
    u is the rate at which recently recovered chickens (R1) move to the R2 compartment
    
    lam, k are parameters of a gamma distribution modelling the timing of prior exposures
    
    dt is the duration of time (in hours) spent in a safe space before T1. It is 0 for control chickens
    and 60 for intervention chickens.
    
    '''
    
    epsi = sigma * 2 # rate from E1 -> E2 & from E2 -> I 
    
    res = np.array( [ 0., 0., 0., 0., 0., 0. ] )
    
    res[0] = 1. - p
    res[1] = p * mean_exp( epsi, lam, k, dt )
    res[2] = p * epsi * mean_xexp( epsi, lam, k, dt )
    res[3] = p * np.square( epsi / ( epsi - mu ) ) * ( mean_exp( mu, lam, k, dt ) - mean_exp( epsi, lam, k, dt ) - ( epsi - mu ) * mean_xexp( epsi, lam, k, dt ) )
    
    #tmp =  mean_exp( u, lam, k, dt ) * ( 1./( mu - u ) - 1./( epsi - u ) - (epsi - mu)/np.square( epsi - u ) )
    
    tmp =  mean_exp( u, lam, k, dt ) * np.square( epsi - mu ) / ( ( mu - u ) * np.square( epsi - u ) )
    tmp += - mean_exp( mu, lam, k, dt ) / ( mu - u )
    tmp += mean_exp( epsi, lam, k, dt ) * (  2 * epsi - mu - u ) / np.square( epsi - u )
    tmp += mean_xexp( epsi, lam, k, dt ) * ( epsi - mu ) / ( epsi - u )

        
    res[4] = p * mu * np.square( epsi / ( epsi - mu ) ) * tmp
    res[5] = 1. - res[:5].sum()

     # check all elems pos and sum to 1
        
    return res
    
    res[res < 0.] = 0.
    res = res / res.sum()
    
    return res

def log_posterior_SEEIRR( theta, data, Nsim, pSurv, priorSheddingTime, priorR0Scale, checkPts=[12,36,84], dtBurnin = 25 * 24, dt01 = 60, Nc=100000, Ni=100000, Nms = [2852, 219], cost = 10.):

    """
    Computes posterior distribution using only broiler chicken data.

    theta is an array with all transformed parameters
    
    data is a json file with chicken counts stratified by chicken type and experimental arm
    
    pSurv is an array with the survival probability of chickens in the market
    
    priorSheddingTime is the a priori mean value of total shedding time (TS=TE+TI).
    It is used to set a strict prior on TE+TI.
    
    priorR0Scale is used to set a prior on the bare reproductive number beta/mu


    """

    # some intrusive output stuff
    global shared_list
    
    lpost = 0.

    # transform parameters bac
    
    lbeta, lsigma_br, lsigma_de, lmu, luu, llam_br, llam_de, lshape_br, lshape_de, lPc_br, lPi_br, lPc_de, lPi_de = theta

    beta     = utils.inverse_transform_positive( lbeta, 0. )
    sigma_br = utils.inverse_transform_positive( lsigma_br, 0. )
    sigma_de = utils.inverse_transform_positive( lsigma_de, 0. )
    mu       = utils.inverse_transform_positive( lmu, 0. )
    uu       = utils.inverse_transform_positive( luu, 0. )
    lam_br   = utils.inverse_transform_positive( llam_br, 0. )
    lam_de   = utils.inverse_transform_positive( llam_de, 0. )
    shape_br = utils.inverse_transform_positive( lshape_br, 0. )
    shape_de = utils.inverse_transform_positive( lshape_de, 0. ) 
    Pc_br    = utils.inverse_transform_bounded( lPc_br, 0., 1. )
    Pi_br    = utils.inverse_transform_bounded( lPi_br, 0., 1. )
    Pc_de    = utils.inverse_transform_bounded( lPc_de, 0., 1. )
    Pi_de    = utils.inverse_transform_bounded( lPi_de, 0., 1. )

    sigmas = [ sigma_br, sigma_de ]

    # compute introduction probabilities

    #shape = 1.5 # set shape gamma hyper-prior for exposure (shape = 1 corresponds to the exponential case)
    pIntroInterv_br  = compute_prob_SEEIRR_gamma_dt( lam_br, shape_br, sigma_br, mu, uu, Pi_br, dt = dt01 ).tolist()
    pIntroControl_br = compute_prob_SEEIRR_gamma_dt( lam_br, shape_br, sigma_br, mu, uu, Pc_br, dt = 0. ).tolist()    
    pIntroInterv_de  = compute_prob_SEEIRR_gamma_dt( lam_de, shape_de, sigma_de, mu, uu, Pi_de, dt = dt01 ).tolist()
    pIntroControl_de = compute_prob_SEEIRR_gamma_dt( lam_de, shape_de, sigma_de, mu, uu, Pc_de, dt = 0. ).tolist()

    pIntroBulk = [ pIntroControl_br, pIntroControl_de ]  

    # compute log-prior

    lprior = 0.
    lprior += utils.hard_penalty( sigma_br, 1./( 24 * 20 ), 1. )    
    lprior += utils.hard_penalty( sigma_de, 1./( 24 * 20 ), 1. )    
    lprior += utils.hard_penalty( mu, 1./(24 * 15), 1./( 12 ) )
    lprior += utils.hard_penalty( uu, 1./(24 * 10), 1./( 24 ) )
    lprior += utils.hard_penalty( lam_br, 0., np.inf )
    lprior += utils.hard_penalty( lam_de, 0., np.inf )
    lprior += utils.hard_penalty( shape_br, 0., np.inf ) # in main analysis, lower bound was 1.
    lprior += utils.hard_penalty( shape_de, 0., np.inf ) # in main analysis, lower bound was 1. 
    lprior += utils.hard_penalty( Pc_br, 0., 0.5 )
    lprior += utils.hard_penalty( Pi_br, 0., 0.5 )
    lprior += utils.hard_penalty( Pc_de, 0., 0.5 )
    lprior += utils.hard_penalty( Pi_de, 0., 0.5 )

    lprior += -lam_br * 0.5 + np.log( lam_br )
    lprior += -lam_de * 0.5 + np.log( lam_de )
    lprior += -shape_br * 0.5 + np.log( shape_br )
    lprior += -shape_de * 0.5 + np.log( shape_de )

    # add jacobians

    lprior += np.log( beta )
    lprior += np.log( sigma_br )
    lprior += np.log( sigma_de )
    lprior += np.log( mu )
    lprior += np.log( uu )
    lprior += np.log( lam_br )
    lprior += np.log( lam_de )
    lprior += np.log( shape_br )
    lprior += np.log( shape_de )
    lprior += np.log( Pc_br * ( 1. - Pc_br ) )
    lprior += np.log( Pi_br * ( 1. - Pi_br ) )
    lprior += np.log( Pc_de * ( 1. - Pc_de ) )
    lprior += np.log( Pi_de * ( 1. - Pi_de ) )
    
    # soft exp prior on beta, based on R0 not being too large
    lprior -= beta / ( priorR0Scale * sum( Nms ) * mu ) + np.log( mu ) # mid prior

    # soft prior on T_I + T_E
    Ti = 1. / mu # should be in h
    Te_br = 1. / sigma_br
    Te_de = 1. / sigma_de
    Te_mean = 0.5 * ( Te_br + Te_de )

    lprior -= 20. * ( Ti + Te_mean - 24 * priorSheddingTime )**2 

    if not np.isfinite( lprior ):
        return -np.infty
    else:
        lpost += lprior


    # update seed for c++ routine
    with rng_gen.get_lock():
        seed = rng_gen.value
        rng_gen.value += 1

    # compute broiler log-likelihood
   
    llikelihood_broiler = log_likelihood_multispecies_SEEIR( data["broiler"], Nsim, seed, beta, sigmas, sigma_br, mu, 1., pIntroInterv_br, pIntroControl_br, 
    pIntroBulk, Nms = Nms, pSurv=pSurv, checkPts=checkPts, dtBurnin = dtBurnin, dt01 = dt01, Nc=Nc, Ni=Ni, cost = cost)

    if not np.isfinite( llikelihood_broiler ):
        return -np.infty
    else:
        lpost += llikelihood_broiler

     # update seed for c++ routine
    with rng_gen.get_lock():
        seed = rng_gen.value
        rng_gen.value += 1

    # compute deshi log-likelihood
    llikelihood_deshi = log_likelihood_multispecies_SEEIR( data["deshi"], Nsim, seed, beta, sigmas, sigma_de, mu, 1., pIntroInterv_de, pIntroControl_de, 
    pIntroBulk, Nms = Nms, pSurv=pSurv, checkPts=checkPts, dtBurnin = dtBurnin, dt01 = dt01, Nc=Nc, Ni=Ni, cost = cost) 
    
    if not np.isfinite( llikelihood_deshi ):
        return -np.infty
    else:
        lpost += llikelihood_deshi
        
    return lpost


if __name__ == "__main__":

    #=== read script args

    nargs = len( sys.argv )
    if nargs == 8:
        nwalkers       = int( sys.argv[1] ) # must be twice ndim (see below)
        nprocesses     = int( sys.argv[2] ) # number of cores to be used
        n_iterations   = int( sys.argv[3] ) # number of mcmc steps
        seed0          = int( sys.argv[4] ) # initial seed
        priorSheddingT = float( sys.argv[5] ) # for prior on shedding time TS=TE+TI
        priorR0Scale   = float( sys.argv[6] ) # for prior on bare reproductive number R
        ctFit          = sys.argv[7] # which data to fit (ct40 or ct33)
    else:
        raise ValueError("Wrong number of arguments")

    Nsim               = 10 # number of simulations
    p_surv             = [1.0, 0.988392768, 0.971711118, 0.953331336, 0.933313355, 0.909279792, 0.880910998, 0.839356708, 0.797942264, 0.754749775, 0.708121067, 0.657616622, 0.6080911, 0.557446809, 0.506143242, 0.457356907, 0.410148836, 0.367375887, 0.322185596, 0.275556887, 0.231944861, 0.193007692, 0.163859754, 0.140885026, 0.120107881, 0.115413046, 0.112396364, 0.109859155, 0.107202078, 0.104604935, 0.101288583, 0.096833483, 0.09231845, 0.086844471, 0.08139047, 0.07613625, 0.070242733, 0.064908601, 0.058815303, 0.053381281, 0.048466687, 0.044151433, 0.03943662, 0.034142443, 0.029407652, 0.025332135, 0.022275497, 0.019538508, 0.017600639, 0.017121167, 0.016661672, 0.016362002, 0.015842573, 0.015423035, 0.014923584, 0.014364199, 0.013864749, 0.013365298, 0.012406353, 0.01162721, 0.010908001, 0.010048946, 0.009329737, 0.008510638, 0.007911298, 0.00697233, 0.006392968, 0.005653781, 0.005114374, 0.004495055, 0.004195385, 0.003975627, 0.003655978, 0.00353611, 0.003516132, 0.003396264, 0.003356308, 0.00333633, 0.00323644, 0.003156528, 0.003016682, 0.002956748, 0.002796923, 0.002537209, 0.002377385, 0.002277495, 0.002137649, 0.002057736, 0.001897912, 0.00171811, 0.00151833, 0.00141844, 0.001258616, 0.001038857, 0.000978923, 0.000938967, 0.000899011, 0.000899011, 0.000859055, 0.000839077, 0.000819099, 0.000779143, 0.000779143, 0.000739187, 0.000719209, 0.000699231, 0.000679253, 0.000659275, 0.000599341, 0.000559385, 0.000539407, 0.000539407, 0.000499451, 0.000459495, 0.000459495, 0.000439517, 0.00039956, 0.000379582, 0.000359604, 0.000319648, 0.000319648] # survival probability of chickens in market (length of stay)

    #=== read and process data
    if (ctFit == 'ct40'):
        path_data = os.path.join('/', *os.getcwd().split('/')[:-1], 'data/data_fit_counts_ct40.json' )
    elif (ctFit == 'ct33'):
        path_data = os.path.join('/', *os.getcwd().split('/')[:-1], 'data/data_fit_counts_ct33.json' )
    else:
        raise ValueError("Please select ct40 or ct33")
    
    
    with open( path_data, 'r' ) as readfile:
        data = json.load( readfile )

    npos_br  = data["npos"]["A"]["H9"]
    nsusc_br = data["nsusc"]["A"]["H9"]
    npos_de  = data["npos"]["B"]["H9"]
    nsusc_de = data["nsusc"]["B"]["H9"]

    data_processed = { "broiler": { "npos": npos_br, "nsusc": nsusc_br }, "deshi": { "npos": npos_de, "nsusc": nsusc_de } }
   
    #=== initialise explored parameters 

    # a) choose initial values
    # b) map values to (-infty, +infty)
    # c) make copies (one for each chain)
    # d) add some noise (initial conditions MUST BE DIFFERENT)
    
    lbeta0     = utils.transform_positive( 5., 0. )
    lsigma_br0 = utils.transform_positive( 1./( 10. ), 0. )
    lsigma_de0 = utils.transform_positive( 1./( 1.25 * 24 ), 0. )
    lmu0       = utils.transform_positive( 1./( priorSheddingT * 24 - 1.7 ), 0. )
    luu0       = utils.transform_positive( 1./( 1.2 * 24 ), 0. )
    llam_br0   = utils.transform_positive( 1./( 12 ), 0. )
    llam_de0   = utils.transform_positive( 1./( 10. * 24 ), 0. ) 
    lshape_br0 = utils.transform_positive( 1.5, 0. )
    lshape_de0 = utils.transform_positive( 1.5, 0. )
    lPc_br0    = utils.transform_bounded( 0.2 )
    lPi_br0    = utils.transform_bounded( 0.06 )                             
    lPc_de0    = utils.transform_bounded( 0.25 )                                                           
    lPi_de0    = utils.transform_bounded( 0.25 )
                                                                                                
    theta0_block = np.array( [ lbeta0, lsigma_br0, lsigma_de0, lmu0, luu0, llam_br0, llam_de0, lshape_br0, lshape_de0, lPc_br0, lPi_br0, lPc_de0, lPi_de0 ] )
    ndim = len( theta0_block ) # number of fit parameters

    theta0 = np.tile( theta0_block, nwalkers ).reshape( nwalkers, ndim )
    jitter = np.random.uniform( low = -0.05, high = 0.05, size = ( nwalkers, ndim ) ) * theta0_block[None,:]
    theta0 += jitter

    if nprocesses == 1:
        doParallel = False
    else:
        doParallel = True

    #=== initialise seed manager

    rng_gen = Value('i', seed0 )
    init_globals( rng_gen )

    manager = Manager()
    shared_list = manager.list()
    
    #=== set output stuff
    
    if (ctFit == 'ct40'):
        path_output = os.path.join( '/', *os.getcwd().split('/')[:-1], 'fit_results/data_fit_counts_ct40.json' )
    elif (ctFit == 'ct33'):
        path_output = os.path.join( '/', *os.getcwd().split('/')[:-1], 'fit_results/data_fit_counts_ct33.json' )
    else:
        raise ValueError("Please select ct40 or ct33")

    if not os.path.exists( path_output ):
        os.makedirs( path_output )
        
    output_name = 'samples_joint_Tshed{}_Rprior{}.h5'.format(
        str(priorSheddingT).replace( '.', 'p' ),
        str(priorR0Scale).replace( '.', 'p' )
    )

    filename = os.path.join( path_output, output_name )
    backend = emcee.backends.HDFBackend( filename )
    backend.reset( nwalkers, ndim )

    #=== run mcmc

    start = time.time()
    
    if doParallel: # parallel calculation
        with Pool( initializer=init_globals, initargs=( rng_gen, ) ) as pool:
            sampler = emcee.EnsembleSampler( nwalkers, ndim, log_posterior_SEEIRR,
            args=[ data_processed, Nsim, p_surv, priorSheddingT, priorR0Scale ], pool = pool, backend = backend  )

            for sample in sampler.sample( theta0, iterations = n_iterations, progress = False ):
                # pause every now and then and print statistics
                if sampler.iteration % 100 == 0:
                    tcheck = time.time()
                    print( "iter =", sampler.iteration, "time elapsed:", tcheck - start, flush = True )
                    print( sampler.acceptance_fraction, flush = True )
    else: # single core calculation
        sampler = emcee.EnsembleSampler( nwalkers, ndim, log_posterior_SEEIRR,
        args=[ data_processed, Nsim, p_surv, priorSheddingT, priorR0Scale ], pool = None, backend = backend  )
        x = sampler.run_mcmc( theta0, n_iterations, progress = True, store = True )
    # save results
    with open( os.path.join( path_output, 'acceptance_frac_joint_H9.npy' ), 'wb' ) as ff:
        np.save( ff, np.array( sampler.acceptance_fraction ) ) 
        
    stop = time.time()

    print( "time elapsed (total):", stop - start )
