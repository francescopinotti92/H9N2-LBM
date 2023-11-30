The current folder contains data and code to perform the analyses described in the accompanying manuscript and supplementary information. 

This text file describes each of the subfolders. 

- data: contains files with counts of chickens testing positive to H9N2 avian influenza at different time points (T0-T4) when the cycle threshold criterion is set to Ct = 33 or Ct = 40 (each list contains counts for T0-T4). This information is stored under the entry "npos" in each .json file. "nsusc" contains counts of chickens that were still susceptible when they left the experiment; this could be because the experiment terminated (at T4) or because the chicken died during the experiment. Broiler and backyard chickens are coded as "A" and "B" respectively. Control and intervention chickens are coded as "0" and "1" respectively.

- simulator: contains C++ code to run various experiments. The subfolder simulator/src also contains code to create a python module (simulate_experiment) that wraps C++ code. simulate_experiment must be compiled before it can be imported and requires pybind11. pybind11 is hosted at https://github.com/pybind/pybind11 and additional details on installation can be found at https://pybind11.readthedocs.io/en/stable/installing.html. To compile simulate_experiment, open the terminal and navigate to the simulator folder (where a file named 'makefile' is located). Once there, type 'make' in the command line to compile simulate_experiment. Compilation requires a working C++11 compiler, like g++ for Ubuntu or clang for macOS.
We were able to compile simulate_experiment on a macOS and an Ubuntu machine. No test has been carried out on Windows machines.

- fit_results: contains some files with MCMC samples. These can be used directly to analyse MCMC output without the need to perform MCMC first.

- python_scripts: contains the following files:
    - utils.py: contains some auxiliary functions.
    - inference.py: contains code to calculate the log-likelihood function.
    - joint_fit.py: contains code to run MCMC using the emcee module (https://emcee.readthedocs.io/en/stable/).
    - launch_fit.sh: used to launch joint_fit.py with default parameter values. The first line fits the model to Ct=40 data, while the second line uses Ct=33 data. Make sure to uncomment the desired line and comment the other. MCMC samples are then saved in the fit_results folder.
    - Analyse output.ipynb: a Jupiter notebook with minimal code to visualise MCMC output and to run further simulations (e.g. for posterior predictive checks). Shown examples also demonstrate how to simulate interventions with direct or indirect transmission. Additional dependencies include numpy, scipy, matplotlib, emcee and simulate_experiment (which has to be compiled before it can be imported, as described above).

