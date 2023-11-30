//
//  simulate_experiment.cpp
//  pySEIRmarket
//
//  Created by user on 27/08/2021.
//

#include "simulate_experiment.hpp"

// as above, but introduction parameters are now fully custom
std::vector<std::vector<std::vector<int>>> simulate_experiment_multispecies_newintro( const int& Nsim, const int& seed, const std::vector<int>& checkPts, const int& dtBurnin, const int& dtT0T1, const int& Nc, const int& Ni, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& sigmaExperiment, const double& mu, const int& nEstages, const int& nIstages, const double& pFarmInterv, const std::vector<double>& pIntroInterv, const std::vector<double>& pIntroControl, const std::vector<std::vector<double>>& pIntroBulk ) {
    
    // res stores results
    std::vector<std::vector<std::vector<int>>> res = std::vector<std::vector<std::vector<int>>>( 2, std::vector<std::vector<int>>( Nsim, std::vector<int>( checkPts.size() + 2, 0 ) ) );
    
    Simulator simulator = Simulator( pSurv, Nms, beta, sigmas, mu, nEstages, nIstages, pIntroBulk );
    
    m_mt.seed( seed );
    
    int Tmax = checkPts.back() + 1;
    
    bool reuse_snapshot = false;
    
    int t = 0;
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        t = 0;
        
        if ( sim % 1 == 0 )
            reuse_snapshot = false;
        else
            reuse_snapshot = true;
        
        simulator.initializeExperiment( Nc, Ni, sigmaExperiment, 1., pFarmInterv, pIntroInterv, pIntroControl, dtBurnin, dtT0T1, reuse_snapshot );
        
        // get control chickens infectious at T1
        res[0][sim][1] = simulator.get_pos_c();
        
        // get intervention chickens infectious at T0 and T1
        res[1][sim][0] = simulator.get_pos_i_T0();
        res[1][sim][1] = simulator.get_pos_i_T1();
        
        int icheck = 0;
        int tcheck = checkPts[0];
        
        while ( t < Tmax ) {
            
            // update dynamics
            
            simulator.runStep();
            
            // check results
            
            if ( t == tcheck ) {
                auto new_pos = simulator.check_new_positive();
                res[0][sim][icheck+2] = new_pos.first; // control birds
                res[1][sim][icheck+2] = new_pos.second; // intervention birds

                ++icheck;
                if ( tcheck < Tmax )
                    tcheck = checkPts[icheck];
                else
                    tcheck = Tmax + 1;

            }
            
            ++t;
        }
        
    }
    return res;
    
    
    
}

// same as above but allows different susceptibilities across chicken types
std::vector<std::vector<std::vector<int>>> simulate_experiment_multispecies_newintro_susceptibility( const int& Nsim, const int& seed, const std::vector<int>& checkPts, const int& dtBurnin, const int& dtT0T1, const int& Nc, const int& Ni, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& sigmaExperiment, const double& mu, const std::vector<double> suscs, const double& suscExperiment, const int& nEstages, const int& nIstages, const double& pFarmInterv, const std::vector<double>& pIntroInterv, const std::vector<double>& pIntroControl, const std::vector<std::vector<double>> & pIntroBulk ) {
    
    
    // res stores results
    std::vector<std::vector<std::vector<int>>> res = std::vector<std::vector<std::vector<int>>>( 2, std::vector<std::vector<int>>( Nsim, std::vector<int>( checkPts.size() + 2, 0 ) ) );
    
    Simulator simulator = Simulator( pSurv, Nms, beta, sigmas, mu, nEstages, nIstages, pIntroBulk );
    
    simulator.setSusceptibility( suscs ); // set susceptibility
    
    m_mt.seed( seed );
    
    int Tmax = checkPts.back() + 1;
    
    bool reuse_snapshot = false;
    
    int t = 0;
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        t = 0;
        
        if ( sim % 1 == 0 )
            reuse_snapshot = false;
        else
            reuse_snapshot = true;
        
        simulator.initializeExperiment( Nc, Ni, sigmaExperiment, suscExperiment, pFarmInterv, pIntroInterv, pIntroControl, dtBurnin, dtT0T1, reuse_snapshot );
        
        // get control chickens infectious at T1
        res[0][sim][1] = simulator.get_pos_c();
        
        // get intervention chickens infectious at T0 and T1
        res[1][sim][0] = simulator.get_pos_i_T0();
        res[1][sim][1] = simulator.get_pos_i_T1();
        
        int icheck = 0;
        int tcheck = checkPts[0];
        
        while ( t < Tmax ) {
            
            // update dynamics
            
            simulator.runStep();
            
            // check results
            
            if ( t == tcheck ) {
                auto new_pos = simulator.check_new_positive();
                res[0][sim][icheck+2] = new_pos.first; // control birds
                res[1][sim][icheck+2] = new_pos.second; // intervention birds

                ++icheck;
                if ( tcheck < Tmax )
                    tcheck = checkPts[icheck];
                else
                    tcheck = Tmax + 1;

            }
            
            ++t;
        }
        
    }
    return res;
    
    
}


std::vector<std::vector<std::vector<std::vector<int>>>> simulate_experiment_multispecies_newintro_compartments( const int& Nsim, const int& seed, const int& Tmax, const int& dtBurnin, const int& dtT0T1, const int& Nc, const int& Ni, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& sigmaExperiment, const double& mu, const int& nEstages, const int& nIstages, const double& pFarmInterv, const std::vector<double>& pIntroInterv, const std::vector<double>& pIntroControl, const std::vector<std::vector<double>>& pIntroBulk ) {
    
    // res stores results
    std::vector<std::vector<std::vector<std::vector<int>>>> res = std::vector<std::vector<std::vector<std::vector<int>>>>( 2, std::vector<std::vector<std::vector<int>>>( Nsim, std::vector<std::vector<int>>( Tmax, std::vector<int>( 6, 0 ) ) ) );
    
    Simulator simulator = Simulator( pSurv, Nms, beta, sigmas, mu, nEstages, nIstages, pIntroBulk );

    
    m_mt.seed( seed );
        
    bool reuse_snapshot = false;
    
    int t = 0;
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        t = 0;
        
        if ( sim % 1 == 0 )
            reuse_snapshot = false;
        else
            reuse_snapshot = true;
        
        simulator.initializeExperiment( Nc, Ni, sigmaExperiment, 1., pFarmInterv, pIntroInterv, pIntroControl, dtBurnin, dtT0T1, reuse_snapshot );
        
        
        while ( t < Tmax ) {
            
            // update dynamics
            
            simulator.runStep();
            
            // check results
            
            simulator.printExperimentCompartments2Vec( res[0][sim][t], "control" );
            simulator.printExperimentCompartments2Vec( res[1][sim][t], "intervention" );

            ++t;
            
        }
        
    }
    return res;
    
    
    
}


// measure the proportion of chickens being susceptible
std::vector<std::vector<int>> measure_exposure_multispecies( const int& Nsim, const int& seed, const int& Tmax, const int& dtBurnin, const int& Nchickens, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& sigmaExperiment, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk ) {
    
    // res stores results
    std::vector<std::vector<int>> res =  std::vector<std::vector<int>>( Nsim, std::vector<int>( Tmax, 0 ) ) ;
    
    Simulator simulator = Simulator( pSurv, Nms, beta, sigmas, mu, nEstages, nIstages, pIntroBulk );
    
    m_mt.seed( seed );
        
    bool reuse_snapshot = false;
    
    int t = 0;
    
    std::vector<double> pIntervS = std::vector<double>( 3 + nEstages + nIstages, 0. );
    pIntervS[0] = 1.;
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        t = 0;
        
        simulator.initializeExperiment( 0, Nchickens, sigmaExperiment, 1., 0., pIntervS, pIntervS, dtBurnin, 0, reuse_snapshot );
        
        // get control chickens infectious at T1
        
        while ( t < Tmax ) {
            
            // update dynamics
            
            simulator.runStep();
            
            // get susceptible
            
            res[sim][t] = simulator.get_comp_i( 0 );
        
            
            ++t;
        }
        
    }
    return res;
    
    
    
}

// simulate market and return all compartments
std::vector<std::vector<std::vector<int>>> simulate_market( const int& Nsim, const int& seed, const int& Tmax, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk ) {
    
    
    std::vector<std::vector<std::vector<int>>> res =  std::vector<std::vector<std::vector<int>>>( Nsim, std::vector<std::vector<int>>( Tmax,  std::vector<int>( 3 + nEstages + nIstages, 0 ) ) ) ;

    
    Simulator simulator = Simulator( pSurv, Nms, beta, sigmas, mu, nEstages, nIstages, pIntroBulk );
    
    m_mt.seed( seed );
        
    bool reuse_snapshot = false;
    
    int t = 0;
    
    std::vector<double> pIntervS = std::vector<double>( 3 + nEstages + nIstages, 0. );
    pIntervS[0] = 1.;
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        t = 0;
        
        simulator.initializeExperiment( 0, 0, 0., 1., 0., pIntervS, pIntervS, 0, 0, reuse_snapshot );
        
        // get control chickens infectious at T1
        
        while ( t < Tmax ) {
            
            // update dynamics
            
            simulator.runStep();
            
            // store compartment content
            
            simulator.printCompartments2Vec( res[sim][t] );
        
            
            ++t;
        }
        
    }
    
    return res;
    
    
}

// Returns number of cumulative infections per batch
std::vector<std::vector<int>> simulate_market_cumulative_infections( const int& Nsim, const int& seed, const int& Tmax, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk ) {
    
    std::vector<std::vector<int>> res =  std::vector<std::vector<int>>( Nsim, std::vector<int>() );

    
    Simulator simulator = Simulator( pSurv, Nms, beta, sigmas, mu, nEstages, nIstages, pIntroBulk );
    
    m_mt.seed( seed );
        
    bool reuse_snapshot = false;
    
    int t = 0;
    
    std::vector<double> pIntervS = std::vector<double>( 3 + nEstages + nIstages, 0. );
    pIntervS[0] = 1.;
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        t = 0;
        
        simulator.initializeExperiment( 0, 0, 0., 1., 0., pIntervS, pIntervS, 0, 0, reuse_snapshot );
        
        // get control chickens infectious at T1
        
        while ( t < Tmax ) {
            
            // update dynamics
            
            simulator.runStep();
            
            
            ++t;
        }
        
        // store cumulative infections
        simulator.printBatchCumulativeInfections2Vec( res[sim] );

        
    }
    
    return res;
    
}

// same as simulate_market but prevents introduction of exposed chickens starting from a certain time
std::vector<std::vector<std::vector<int>>> simulate_market_nointro( const int& Nsim, const int& seed, const int& Tmin, const int& Tmax, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk ) {
    
    
    std::vector<std::vector<std::vector<int>>> res =  std::vector<std::vector<std::vector<int>>>( Nsim, std::vector<std::vector<int>>( Tmax,  std::vector<int>( 3 + nEstages + nIstages, 0 ) ) ) ;

    
    Simulator simulator = Simulator( pSurv, Nms, beta, sigmas, mu, nEstages, nIstages, pIntroBulk );
    
    m_mt.seed( seed );
        
    bool reuse_snapshot = false;
    
    int t = 0;
    
    std::vector<double> pIntervS = std::vector<double>( 3 + nEstages + nIstages, 0. );
    pIntervS[0] = 1.;
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        t = 0;
        
        simulator.initializeExperiment( 0, 0, 0.,1., 0., pIntervS, pIntervS, 0, 0, reuse_snapshot );
        
        // get control chickens infectious at T1
        
        while ( t < Tmax ) {
            
            // update dynamics
            
            if ( t < Tmin )
                simulator.runStep();
            else
                simulator.runStepSkipIntro(); // insert only susceptibles from now on
            
            // store compartment content
            
            simulator.printCompartments2Vec( res[sim][t] );
        
            
            ++t;
        }
        
    }
    
    return res;
    
    
}


std::vector<std::vector<std::vector<int>>> simulate_market_withclosure( const int& Nsim, const int& seed, const int& Tmax, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk, const int& t_close, const int& t_reopen )  {
    
    assert( t_reopen > t_close );
    
    std::vector<std::vector<std::vector<int>>> res =  std::vector<std::vector<std::vector<int>>>( Nsim, std::vector<std::vector<int>>( Tmax,  std::vector<int>( 3 + nEstages + nIstages, 0 ) ) ) ;

    
    Simulator simulator = Simulator( pSurv, Nms, beta, sigmas, mu, nEstages, nIstages, pIntroBulk );
    
    m_mt.seed( seed );
        
    bool reuse_snapshot = false;
    
    int t = 0;
    
    std::vector<double> pIntervS = std::vector<double>( 3 + nEstages + nIstages, 0. );
    pIntervS[0] = 1.;
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        t = 0;
        
        simulator.initializeExperiment( 0, 0, 0., 1., 0., pIntervS, pIntervS, 0, 0, reuse_snapshot );
        
        // get control chickens infectious at T1
        
        while ( t < Tmax ) {
            
            // update dynamics
            
            if ( ( t < t_close ) or ( t >= t_reopen ) ) {
            
                simulator.runStep();
                
            }
            else {
                
                if ( t == t_close ) {
                    
                    simulator.forceMarketClosure( 0. );

                    
                }
                
            }
            
            // store compartment content
            
            simulator.printCompartments2Vec( res[sim][t] );
        
            
            ++t;
        }
        
    }
    
    return res;
    
    
}


std::vector<std::vector<std::vector<int>>> simulate_market_environment( const int& Nsim, const int& seed, const int& Tmax, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta_environment, const double& decay_rate, const std::vector<double>& sigmas, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk, const double& disinfect ) {
    
    
    // last entry is the environment
    std::vector<std::vector<std::vector<int>>> res =  std::vector<std::vector<std::vector<int>>>( Nsim, std::vector<std::vector<int>>( Tmax,  std::vector<int>( 3 + nEstages + nIstages + 1, 0 ) ) ) ;

    
    Simulator simulator = Simulator( pSurv, Nms, beta_environment, decay_rate, sigmas, mu, nEstages, nIstages, pIntroBulk );
    
    m_mt.seed( seed );
        
    bool reuse_snapshot = false;
    
    int t = 0;
    
    std::vector<double> pIntervS = std::vector<double>( 3 + nEstages + nIstages + 1, 0. );
    pIntervS[0] = 1.;
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        t = 0;
        
        simulator.initializeExperiment( 0, 0, 0., 1., 0., pIntervS, pIntervS, 0, 0, reuse_snapshot );
        
        // get control chickens infectious at T1
        
        while ( t < Tmax ) {
            
            // update dynamics
            
            if ( t % 24 == 0 ) {
                // disinfect daily
                simulator.disinfectEnvironment( disinfect );
            }
            simulator.runStepEnvironment();
            
            // store compartment content
            
            simulator.printCompartments2Vec( res[sim][t] );
            res[sim][t][3 + nEstages + nIstages] = simulator.get_environment() ; // write environment
            
            ++t;
            
        }
        
    }
    
    return res;
    
}


std::vector<std::vector<std::vector<int>>> simulate_market_environment_withclosure( const int& Nsim, const int& seed, const int& Tmax, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta_environment, const double& decay_rate, const std::vector<double>& sigmas, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk, const int& t_close, const int& t_reopen, const double& pEnvReduction )  {
    
    
    // last entry is the environment
    std::vector<std::vector<std::vector<int>>> res =  std::vector<std::vector<std::vector<int>>>( Nsim, std::vector<std::vector<int>>( Tmax,  std::vector<int>( 3 + nEstages + nIstages + 1, 0 ) ) ) ;

    
    Simulator simulator = Simulator( pSurv, Nms, beta_environment, decay_rate, sigmas, mu, nEstages, nIstages, pIntroBulk );
    
    m_mt.seed( seed );
        
    bool reuse_snapshot = false;
    
    int t = 0;
    
    std::vector<double> pIntervS = std::vector<double>( 3 + nEstages + nIstages + 1, 0. );
    pIntervS[0] = 1.;
    
    for ( int sim = 0; sim < Nsim; ++sim ) {
        
        t = 0;
        
        simulator.initializeExperiment( 0, 0, 0., 1., 0., pIntervS, pIntervS, 0, 0, reuse_snapshot );
        
        // get control chickens infectious at T1
        
        while ( t < Tmax ) {
            
            // update dynamics
            
            if ( ( t < t_close ) or ( t >= t_reopen ) ) {
            
                simulator.runStepEnvironment();
                
            }
            else {
                
                if ( t == t_close ) {
                    
                    simulator.forceMarketClosure( pEnvReduction );

                }
                
            }
            
            // store compartment content
            
            simulator.printCompartments2Vec( res[sim][t] );
            res[sim][t][3 + nEstages + nIstages] = simulator.get_environment() ; // write environment
            
            ++t;
            
        }
        
    }
    
    return res;
    
}
