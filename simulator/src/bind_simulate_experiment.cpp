//
//  bind_simulate_experiment.cpp
//  SimulateSingleMarket
//
//  Created by user on 27/08/2021.
//

#include "simulate_experiment.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// module: is the name of the python module we are creating, should coincide with cpp file name?
// python complained when I used 'mymodule' instead of 'module'
// m: is an instance of py::module_, necessary to create bindings

namespace py = pybind11;

PYBIND11_MODULE(simulate_experiment, m) {
    m.doc() = "python binding for c++ code simulating SEIR dynamics in a market"; // optional module docstring

    //py::bind_vector<std::vector<int>>(m, "IntVector3D");
    py::bind_vector<std::vector<std::vector<std::vector<int>>>>(m, "IntVector3D");
    
    m.def("simulate_experiment_multispecies_newintro", &simulate_experiment_multispecies_newintro, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("checkPts"),
          py::arg("dtBurnin"),
          py::arg("dt01"),
          py::arg("Nc"),
          py::arg("Ni"),
          py::arg("pSurv"),
          py::arg("Nms"),
          py::arg("beta"),
          py::arg("sigmas"),
          py::arg("sigmaExp"),
          py::arg("mu"),
          py::arg("nEstages"),
          py::arg("nIstages"),
          py::arg("pFarmInterv"),
          py::arg("pIntroInterv"),
          py::arg("pIntroControl"),
          py::arg("pIntroBulk") );
    
    m.def("simulate_experiment_multispecies_newintro_susceptibility", &simulate_experiment_multispecies_newintro_susceptibility, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("checkPts"),
          py::arg("dtBurnin"),
          py::arg("dt01"),
          py::arg("Nc"),
          py::arg("Ni"),
          py::arg("pSurv"),
          py::arg("Nms"),
          py::arg("beta"),
          py::arg("sigmas"),
          py::arg("sigmaExp"),
          py::arg("mu"),
          py::arg("suscs"),
          py::arg("suscExp"),
          py::arg("nEstages"),
          py::arg("nIstages"),
          py::arg("pFarmInterv"),
          py::arg("pIntroInterv"),
          py::arg("pIntroControl"),
          py::arg("pIntroBulk") );
    
    m.def("simulate_experiment_multispecies_newintro_compartments", &simulate_experiment_multispecies_newintro_compartments, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("dtBurnin"),
          py::arg("dt01"),
          py::arg("Nc"),
          py::arg("Ni"),
          py::arg("pSurv"),
          py::arg("Nms"),
          py::arg("beta"),
          py::arg("sigmas"),
          py::arg("sigmaExp"),
          py::arg("mu"),
          py::arg("nEstages"),
          py::arg("nIstages"),
          py::arg("pFarmInterv"),
          py::arg("pIntroInterv"),
          py::arg("pIntroControl"),
          py::arg("pIntroBulk") );
    
    m.def("measure_exposure_multispecies", &measure_exposure_multispecies, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("dtBurnin"),
          py::arg("Nchickens"),
          py::arg("pSurv"),
          py::arg("Nms"),
          py::arg("beta"),
          py::arg("sigmas"),
          py::arg("sigmaExp"),
          py::arg("mu"),
          py::arg("nEstages"),
          py::arg("nIstages"),
          py::arg("pIntroBulk") );
    
    
    m.def("simulate_market", &simulate_market, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("pSurv"),
          py::arg("Nms"),
          py::arg("beta"),
          py::arg("sigmas"),
          py::arg("mu"),
          py::arg("nEstages"),
          py::arg("nIstages"),
          py::arg("pIntroBulk") );
    
    m.def("simulate_market_nointro", &simulate_market_nointro, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmin"),
          py::arg("Tmax"),
          py::arg("pSurv"),
          py::arg("Nms"),
          py::arg("beta"),
          py::arg("sigmas"),
          py::arg("mu"),
          py::arg("nEstages"),
          py::arg("nIstages"),
          py::arg("pIntroBulk") );
    
    m.def("simulate_market_withclosure", &simulate_market_withclosure, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("pSurv"),
          py::arg("Nms"),
          py::arg("beta"),
          py::arg("sigmas"),
          py::arg("mu"),
          py::arg("nEstages"),
          py::arg("nIstages"),
          py::arg("pIntroBulk"),
          py::arg("t_close"),
          py::arg("t_reopen") );
    
    m.def("simulate_market_cumulative_infections", &simulate_market_cumulative_infections, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("pSurv"),
          py::arg("Nms"),
          py::arg("beta"),
          py::arg("sigmas"),
          py::arg("mu"),
          py::arg("nEstages"),
          py::arg("nIstages"),
          py::arg("pIntroBulk") );
    
    
    m.def("simulate_market_environment", &simulate_market_environment, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("pSurv"),
          py::arg("Nms"),
          py::arg("beta_environment"),
          py::arg("decay_rate"),
          py::arg("sigmas"),
          py::arg("mu"),
          py::arg("nEstages"),
          py::arg("nIstages"),
          py::arg("pIntroBulk"),
          py::arg("disinfect") = 0.);
    
    m.def("simulate_market_environment_withclosure", &simulate_market_environment_withclosure, "Simulates spread",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("Tmax"),
          py::arg("pSurv"),
          py::arg("Nms"),
          py::arg("beta_environment"),
          py::arg("decay_rate"),
          py::arg("sigmas"),
          py::arg("mu"),
          py::arg("nEstages"),
          py::arg("nIstages"),
          py::arg("pIntroBulk"),
          py::arg("t_close"),
          py::arg("t_reopen"),
          py::arg("pEnvReduction") );
          
    
    
    /*
    m.def("simulate_market", &simulate_market, "Simulates spread",
          py::arg("Nsim"),
          py::arg("Tmax"),
          py::arg("dt"),
          py::arg("seed"),
          py::arg("pSurv"),
          py::arg("Nm"),
          py::arg("beta"),
          py::arg("sigma"),
          py::arg("mu"),
          py::arg("nEstages"),
          py::arg("nIstages"),
          py::arg("pEIc") );
    
    m.def("measure_R0", &measure_R0, "measure R0",
          py::arg("Nsim"),
          py::arg("seed"),
          py::arg("dtBurn"),
          py::arg("pSurv"),
          py::arg("Nm"),
          py::arg("beta"),
          py::arg("sigma"),
          py::arg("mu"),
          py::arg("nEstages"),
          py::arg("nIstages") );
     */
    
    // use py::arg("name_arg") to enable named arguments in python
    // use py::arg("name_arg") = x to enable named arguments with default parameters
    /*
    m.def("add", &add, "A function which adds two numbers",
    py::arg("i") = 1, py::arg("j") = 2);
     */
}
