//
//  main.cpp
//  SimulateSingleMarket
//
//  Created by user on 03/08/2021.
//

#include "simulate_experiment.hpp"
#include "random.hpp"

int main(int argc, const char * argv[]) {
    // insert code here...
    
    int Nsim = 1;
    int seed = 0;
    std::vector<int> checkPts = {12,36,84};
    int Nc = 5;
    int Ni = 5;
    int dtBurn = 24 * 50;
    int dtT0T1 = 60;
    
    /*
    std::vector<double> pSurv = {1.0,
        0.9, 0.81,0.7290000000000001,
        0.6561,0.5904900000000001,
        0.531441,0.4782969000000001,
        0.4304672100000001,0.3874204890000001,
        0.3486784401000001,0.31381059609000006,
        0.2824295364810001,0.2541865828329001,
        0.2287679245496101,0.20589113209464907,
        0.18530201888518416,0.16677181699666577};
     */
    
    std::vector<double> pSurv = {1.0, 0.988392768, 0.971711118, 0.953331336, 0.933313355, 0.909279792, 0.880910998, 0.839356708, 0.797942264, 0.754749775, 0.708121067, 0.657616622, 0.6080911, 0.557446809, 0.506143242, 0.457356907, 0.410148836, 0.367375887, 0.322185596, 0.275556887, 0.231944861, 0.193007692, 0.163859754, 0.140885026, 0.120107881, 0.115413046, 0.112396364, 0.109859155, 0.107202078, 0.104604935, 0.101288583, 0.096833483, 0.09231845, 0.086844471, 0.08139047, 0.07613625, 0.070242733, 0.064908601, 0.058815303, 0.053381281, 0.048466687, 0.044151433, 0.03943662, 0.034142443, 0.029407652, 0.025332135, 0.022275497, 0.019538508, 0.017600639, 0.017121167, 0.016661672, 0.016362002, 0.015842573, 0.015423035, 0.014923584, 0.014364199, 0.013864749, 0.013365298, 0.012406353, 0.01162721, 0.010908001, 0.010048946, 0.009329737, 0.008510638, 0.007911298, 0.00697233, 0.006392968, 0.005653781, 0.005114374, 0.004495055, 0.004195385, 0.003975627, 0.003655978, 0.00353611, 0.003516132, 0.003396264, 0.003356308, 0.00333633, 0.00323644, 0.003156528, 0.003016682, 0.002956748, 0.002796923, 0.002537209, 0.002377385, 0.002277495, 0.002137649, 0.002057736, 0.001897912, 0.00171811, 0.00151833, 0.00141844, 0.001258616, 0.001038857, 0.000978923, 0.000938967, 0.000899011, 0.000899011, 0.000859055, 0.000839077, 0.000819099, 0.000779143, 0.000779143, 0.000739187, 0.000719209, 0.000699231, 0.000679253, 0.000659275, 0.000599341, 0.000559385, 0.000539407, 0.000539407, 0.000499451, 0.000459495, 0.000459495, 0.000439517, 0.00039956, 0.000379582, 0.000359604, 0.000319648, 0.0};
    
    std::vector<double> sigmas = { 0.1, 0.1 };
    double sigmaExp = 0.1;
    double mu = 1./24;
    double beta = mu * 10.;
    int nEstages = 2;
    int nIstages = 1;
   

        
    //auto res = simulate_experiment( Nsim, seed, checkPts, dtBurn, dtT0T1, Nc, Ni, pSurv, Nm, beta, sigma, mu, nEstages, nIstages, pc, piCL, piInd );
    
    std::vector<double> pIntroInterv = {1.,0.,0.,0.,0.,0.};
    std::vector<double> pIntroControlBR = {0.9,0.,0.,0.1,0.,0.};
    std::vector<double> pIntroControlDE = {0.9,0.,0.,0.1,0.,0.};
    
    std::vector<int> Nms = {3000, 200};
    //simulate_experiment_multispecies_newintro( Nsim, seed, checkPts, dtBurn, dtT0T1, Nc, Ni, pSurv, Nms, beta, sigmas, sigmaExp, mu, nEstages, nIstages, 1., pIntroInterv, pIntroControlBR, { pIntroControlBR, pIntroControlDE } );
    int Tmax = 24 * 20;
    
    simulate_market_cumulative_infections(Nsim, seed, Tmax, pSurv, Nms, beta, sigmas, mu, nEstages, nIstages, { pIntroControlBR, pIntroControlDE } );
    
    
    
    /*
    int Tmax = 24 * 30;
    int dt = 1;
    auto res = simulate_market(Nsim, Tmax, dt, seed, pSurv, Nm, beta, sigma, mu, nEstages, nIstages, pc);
    */
    
    
    return 0;
}
