//
//  simulator.hpp
//  SimulateSingleMarket
//
//  Created by user on 27/08/2021.
//

#ifndef simulator_hpp
#define simulator_hpp

#include "random.hpp"
#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>

void printVec( const std::vector<int>& v );

//===== Simulator class =====//

class Simulator {
public:
    
    /* Constructors */

    // multiple species, fully custom compartment initialisation
    Simulator( const std::vector<double>& pSurv, const std::vector<int>& Nm, const double& beta, const std::vector<double>& sigma, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>>& pIntro );
    
    Simulator( const std::vector<double>& pSurv, const std::vector<int>& Nm, const double& beta_environment, const double& decay_rate, const std::vector<double>& sigma, const double& mu, const int& nEstages, const int& nIstages, const std::vector<std::vector<double>>& pIntro );
    
    /* Functions to setup experiment chickens */
    
    // fully custom initial conditions (pE, pI is not split uniformly across ranges)
    void initializeExperiment( const int& Nc, const int& Ni, const double& sigmaExperiment, const double& birdSuscExperiment, const double& pClusterInterv, const std::vector<double>& pIntroInterv, const std::vector<double>& pIntroControl, const int& dt_burn, const int& dt_t0_t1, const bool& reuse_snaphot = false );
    
    /* Run the simulation */
    
    void runStep( const bool& doEpi = true );
    void runStepSkipIntro( const bool& doEpi = true );
    void runStepEnvironment( const bool& doEpi = true );

    /* Interventions */
    void forceMarketClosure( const double& pEnvReduction );
    void disinfectEnvironment( const double& disinfect );

    
    int get_Ie_c();
    int get_pos_c();
    int get_pos_i_T0() const { return pos_i_T0; }
    int get_pos_i_T1() const { return pos_i_T1; }
    int get_comp_i( const int& icomp ) { return interventionChickens[icomp]; }
    int get_environment() const { return envCont; }
    std::pair<int, int> check_new_positive( const bool& reset = true );
    void printCompartments2Vec( std::vector<int>& v );
    void printExperimentCompartments2Vec( std::vector<int>& v, const std::string& group );
    void printBatchCumulativeInfections2Vec( std::vector<int>& v );
    
    void setSusceptibility( std::vector<double> vals ) { susc = vals ; }
    
    // mostly used when computing R0 and experimenting with individual chickens
    
    /*
    int get_susc_size();
    int get_random_susc_chicken();
    int getTimeToSale( const int& batchID );
    int getTimeToI();
    int getTimeToR();
     */
    
private:
    
    // private methods
    
    std::pair<int,int> updateEpiBirdVec( std::vector<int>& v, const double& expsig, const double& birdSusc );
    void clearBirdVec( std::vector<int>& v );
    
    void fillBirdVecFixedProb( std::vector<int>& v, const int& n, const std::vector<double>& p );

    void fillBirdVecClustered( std::vector<int>& v, const int& n, const double& pFarm, const std::vector<double>& p );
    
    
    void runEpiStep( const bool& doBulkUpdate = true, const bool& preventInfection = false );
    void runEpiStepEnvironment( const bool& doBulkUpdate = true, const bool& preventInfection = false );

    void allocateDailyShipment( const bool& onlySusceptible = false );
    void sellChickens();
    void compute_foi();
    void compute_foi_environment();
    void updateEnvironment();
    
    // time-related variables
    
    int t;
    int hour;
    int day;
    int daymod;
    
    // miscellanea
    
    int nMaxBatches;
    int nCompartments;
    int nSpecies;
    std::vector<int> batchAge; // number of days a batch has been in the market
    
    // various indexes
    
    int ixS;
    int ixE0;
    int ixEn;
    int ixI0;
    int ixIn;
    int ixR1;
    int ixR2;
    
    int envCont;

    
    std::vector<double> pSellCorrected;
    
    // compartments
    
    std::vector<std::vector<std::vector<int>>> bulkChickens;
    std::vector<std::vector<std::vector<int>>> bulkChickensSnapshot;
    std::vector<int> controlChickens;
    std::vector<int> interventionChickens;
    
    std::vector<int> cumulNewInfections;
    std::vector<int> cumulNewInfectionsStore;
    
    int _Nc;
    int _Ni;
    
    // epi params
    
    double beta;
    std::vector<double> sigma;
    std::vector<double> expsigma;
    double expsigmaExperiment;
    double mu;
    double expmu;
    int nEstages;
    int nIstages;
    std::vector<double> susc ; // intrinsic susceptibility
    double birdSuscExperiment ;

    
    // environment
    double beta_environment;
    double decay_rate;
    double exp_decay_rate;
    
    std::vector<int> NmVec;
    std::vector<std::vector<double>> pIntroVec;

    double foi;
    double Nfoi; // scaling factor FOI
    
    // output
    
    int pos_i_T0;
    int pos_i_T1;
    
    int newPos_c;
    int newPos_i;
    
};



#endif /* simulator_hpp */
