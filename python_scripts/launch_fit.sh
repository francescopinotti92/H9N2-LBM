#source activate /data/zool-cbm/zool2492/env_emcee

#== arguments
# 1) Number of walkers: must be even and at least twice the number of fitted parameters
# 2) Number of MCMC iterations
# 3) Random number generator seed
# 4) Prior hyperparam on duration of shedding time (i.e. incubation + infectious periods)
# 5) Prior hyperparam on basic reproductive number
# 6) Which data to fit (ct40 or ct33)

python3 joint_fit.py 26 1 100 0 5 0.005 "ct40"  # fitting to Ct=40 data
#python3 joint_fit.py 26 1 100 1 5 0.005 "ct33" # fitting to Ct=33 data
