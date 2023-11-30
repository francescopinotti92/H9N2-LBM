//
//  simulate_experiment.hpp
//  pySEIRmarket
//
//  Created by user on 27/08/2021.
//

#ifndef simulate_experiment_hpp
#define simulate_experiment_hpp

#include "simulator.hpp"

// as above, but introduction parameters are now fully custom
std::vector<std::vector<std::vector<int>>> simulate_experiment_multispecies_newintro( const int& Nsim, const int& seed, const std::vector<int>& checkPts, const int& dtBurnin, const int& dtT0T1, const int& Nc, const int& Ni, const std::vector<double>& pSurv, const std::vector<int>& Nm, const double& beta, const std::vector<double>& sigmas, const double& sigmaExperiment, const double& mu, const int& nEstages, const int& nIstages, const double& pFarmInterv, const std::vector<double>& pIntroInterv, const std::vector<double>& pIntroControl, const std::vector<std::vector<double>> & pIntroBulk );

std::vector<std::vector<std::vector<int>>> simulate_experiment_multispecies_newintro_susceptibility( const int& Nsim, const int& seed, const std::vector<int>& checkPts, const int& dtBurnin, const int& dtT0T1, const int& Nc, const int& Ni, const std::vector<double>& pSurv, const std::vector<int>& Nm, const double& beta, const std::vector<double>& sigmas, const double& sigmaExperiment, const double& mu, const std::vector<double> suscs, const double& suscExperiment, const int& nEstages, const int& nIstages, const double& pFarmInterv, const std::vector<double>& pIntroInterv, const std::vector<double>& pIntroControl, const std::vector<std::vector<double>> & pIntroBulk );

std::vector<std::vector<std::vector<std::vector<int>>>> simulate_experiment_multispecies_newintro_compartments( const int& Nsim, const int& seed, const int& Tmax, const int& dtBurnin, const int& dtT0T1, const int& Nc, const int& Ni, const std::vector<double>& pSurv, const std::vector<int>& Nm, const double& beta, const std::vector<double>& sigmas, const double& sigmaExperiment, const double& mu, const int& nEstages, const int& nIstages, const double& pFarmInterv, const std::vector<double>& pIntroInterv, const std::vector<double>& pIntroControl, const std::vector<std::vector<double>>& pIntroBulk );

std::vector<std::vector<int>> measure_exposure_multispecies( const int& Nsim, const int& seed, const int& Tmax, const int& dtBurnin, const int& Nchickens, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& sigmaExperiment, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk );


std::vector<std::vector<std::vector<int>>> simulate_market( const int& Nsim, const int& seed, const int& Tmax, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk );

std::vector<std::vector<std::vector<int>>> simulate_market_nointro( const int& Nsim, const int& seed, const int& Tmin, const int& Tmax, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk );

std::vector<std::vector<std::vector<int>>> simulate_market_withclosure( const int& Nsim, const int& seed, const int& Tmax, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk, const int& t_close, const int& t_reopen );

std::vector<std::vector<int>> simulate_market_cumulative_infections( const int& Nsim, const int& seed, const int& Tmax, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta, const std::vector<double>& sigmas, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk );


//=== Functions using environmental transmission

std::vector<std::vector<std::vector<int>>> simulate_market_environment( const int& Nsim, const int& seed, const int& Tmax, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta_environment, const double& decay_rate, const std::vector<double>& sigmas, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk, const double& disinfect = 0. );

std::vector<std::vector<std::vector<int>>> simulate_market_environment_withclosure( const int& Nsim, const int& seed, const int& Tmax, const std::vector<double>& pSurv, const std::vector<int>& Nms, const double& beta_environment, const double& decay_rate, const std::vector<double>& sigmas, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>> & pIntroBulk, const int& t_close, const int& t_reopen, const double& pEnvReduction );


#endif /* simulate_experiment_hpp */
