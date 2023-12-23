#########################################################################
# model_ensemble.jl                                                     #
#                                                                       #
# Runs full model chain ensemble (looking at the full parameter         #
#   uncertainty) using MimiBRICK and SNEASY.                            # #                                                                       #
#                                                                       #
# This script requires the calibrated SNEASY-BRICK parameter set        #
#   (e.g. from https://zenodo.org/records/6626335). We used             #
#   parameters_subsample_sneasybrick.csv; using the full MCMC chain     #
#   would require some further subsampling.                             #
#                                                                       #
# Note: It takes about 25 minutes to run 10,000 samples with this       #
#    script                                                             #
#########################################################################

# activate the environment
using Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

# import necessary packages
using Random # ability to set seeds
using QuasiMonteCarlo # QMC sampling (e.g. Sobol sampling)
using DataFrames # tabular data structure
using Distributions # API for working with statistical distribution
using CSVFiles # need to keep for "load" function with CSVs

include(joinpath(@__DIR__, "functions.jl")) # include functions from other scripts
include(joinpath(@__DIR__, "gsa_functions.jl"))
include(joinpath(@__DIR__, "run_sneasy_brick.jl")) # include functions from other scripts

# initial set up
Random.seed!(1)
output_dir = "default"
start_year = 1850
end_year = 2300
num_samples = 100_000

# generate parameter samples
# set bounds for the 3 emissions parameters (growth, t_peak, decline)
lb = zeros(3) # lower bound (0)
ub = ones(3) # upper bound (1)

# create design matrices for emissions parameters (A1 and B1)
A1 = QuasiMonteCarlo.sample(num_samples, lb, ub, SobolSample()) # all values in range 0 to 1

# define distributions for emissions parameters
γ_g_dist     = truncated(Normal(0.004,0.0075), 0.001, Inf)  # growth parameter
t_peak_dist  = truncated(Normal(2070,25), 2030, 2200)       # peaking time
γ_d_dist     = truncated(Normal(0.07,0.05), 0.001, 0.2)     # decline parameter

# combine emissions paramaters' distributions into a matrix
emissions_dist = hcat(γ_g_dist, t_peak_dist, γ_d_dist)

# convert the 0-1 quantiles into the actual values for growth, peaking, and decline
for i in 1:3 # loop through the three emissions parameters
    A1[i,:] = quantile(emissions_dist[i], A1[i,:])
end

# truncate the peaking years from Float to Integer
A1[2,:] = trunc.(Int64, A1[2,:])

# read in subsample of MCMC parameters
mcmc_params = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "data", "calibrated_parameters", "parameters_subsample_sneasybrick.csv"))) # read in subsample
select!(mcmc_params, Not([:sigma_whitenoise_co2, :alpha0_CO2])) # drop these parameters related to CO2 flux postprocessing; we don't use them
num_params = size(mcmc_params, 2)

# initialize storage for A2 and B2 matrices
A2 = zeros(num_samples,num_params)

# sample values for each parameter to create A2 and B2 matrices
mcmc_idx = sample(1:size(mcmc_params, 1), num_samples)
for i in 1:num_params
    for j in 1:num_samples
        A2[j, i] = mcmc_params[mcmc_idx[j], i]
    end
end

# combine design matrices and add parameter names
A = hcat(A1', A2)
param_names = ["gamma_g", "t_peak", "gamma_d", names(mcmc_params)...]
A = DataFrame(A, param_names)
A[!, :lw_random_sample] = rand(Normal(0.0003, 0.00018), num_samples)
calibrated_params = Matrix(A)

run_model(A, start_year, end_year, output_dir) # run the model using this design matrix