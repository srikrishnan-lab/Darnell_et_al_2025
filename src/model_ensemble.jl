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
using DataFrames # tabular data structure
using Distributions # API for working with statistical distribution
using CSVFiles # need to keep for "load" function with CSVs

include(joinpath(@__DIR__, "functions.jl")) # include functions from other scripts
include(joinpath(@__DIR__, "run_sneasy_brick.jl")) # include functions from other scripts

# initial set up
Random.seed!(1)
scenarios = ["default", "optimistic", "pessimistic"]
scenario = scenarios[parse(Int64, ARGS[1])]
println("$scenario")
output_dir = scenario
start_year = 1850
end_year = 2300
num_samples = 100_000

# generate parameter samples
emis_params = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "data", "emissions", scenario, "parameters.csv")))

# sample parameters
emis_idx = sample(1:size(emis_params, 1), num_samples)
emissions_params = emis_params[emis_idx, :]

# read in subsample of MCMC parameters
mcmc_params = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "data", "calibrated_parameters", "parameters_subsample_sneasybrick.csv"))) # read in subsample
select!(mcmc_params, Not([:sigma_whitenoise_co2, :alpha0_CO2])) # drop these parameters related to CO2 flux postprocessing; we don't use them
num_params = size(mcmc_params, 2)

# sample values for each parameter to create A2 and B2 matrices
mcmc_idx = sample(1:size(mcmc_params, 1), num_samples)
brick_params = mcmc_params[mcmc_idx, :]

# combine design matrices and add parameter names
A = hcat(emissions_params, brick_params)
A[!, :lw_random_sample] = rand(Normal(0.0003, 0.00018), num_samples)

run_model(A, start_year, end_year, output_dir) # run the model using this design matrix