#########################################################################
# peaking_ensemble.jl                                                   #
#                                                                       #
# Runs model chain ensemble (focusing on peaking year uncertainty)      #
#    using MimiBRICK and SNEASY.                                        #
#                                                                       #
# This script requires the calibrated SNEASY-BRICK parameter set        #
#   (e.g. from https://zenodo.org/records/6626335). We used             #
#   parameters_subsample_sneasybrick.csv; using the full MCMC chain     #
#   would require some further subsampling.                             #
#                                                                       #
#########################################################################

# load environment and packages
using Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

using Random # random seeds
using DataFrames # tabular data structures
using Distributions # API for statistical distributions
using CSVFiles # API for working with CSV files

γ_g_dist     = truncated(Normal(0.004,0.0075), 0.001, Inf)  # growth parameter
t_peak_dist  = truncated(Normal(2070,25), 2030, 2200)       # peaking time
γ_d_dist     = truncated(Normal(0.07,0.05), 0.001, 0.2)     # decline parameter

include(joinpath(@__DIR__, "functions.jl")) # include functions from other scripts
include(joinpath(@__DIR__, "run_sneasy_brick.jl")) # include functions from other scripts

# experiment settings
output_dir = "peaking"
n_samp = 100_000 # number of samples per peaking time
start_year = 1850
end_year = 2300

# load MCMC output
mcmc_params = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "data", "calibrated_parameters", "parameters_subsample_sneasybrick.csv"))) # read in subsample
select!(mcmc_params, Not([:sigma_whitenoise_co2, :alpha0_CO2]))
num_params = size(mcmc_params, 2)

# generate emissions samples and combine into design matrix
samples = zeros(n_samp, 2 + num_params)

samples[:, 1] = rand(γ_g_dist, n_samp)
samples[:, 2] = rand(γ_d_dist, n_samp)
mcmc_idx = rand(1:size(mcmc_params, 1), n_samp)
samples[:, 3:end] = Matrix(mcmc_params[mcmc_idx, :])
yrs = collect(2050:10:2100)
peak_yrs = trunc.(Int64, repeat(yrs, inner=n_samp))
parameters = repeat(samples, outer = [length(yrs), 1])
A = hcat(parameters[:, 1], peak_yrs, parameters[:, 2:end])

param_names = ["gamma_g", "t_peak", "gamma_d", names(mcmc_params)...]
A = DataFrame(A, param_names)
A[:, :lw_random_sample] = repeat(rand(Normal(0.0003, 0.00018), n_samp), outer=length(yrs))

# run experiment
run_model(A, start_year, end_year, output_dir)