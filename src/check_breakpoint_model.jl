#########################################################################
# check_breakpoint_model.jl                                             #
#                                                                       #
# Checks whether the breakpoint transient sensitivity models are        #
#   are statistically significant.                                      #
#                                                                       #
# This script requires the ensemble output to be present in             #
#   `results/default` .                                                 #
#                                                                       #
#########################################################################

# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSVs
using DataFrames # data structure for indices

using Distributions
using StatsBase # get mean function and density
using GLM

# load ensemble
output_dir = "results/default"
temperature = DataFrame(CSVFiles.load(joinpath(output_dir, "temperature.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
gmslr = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))

# define function to normalize data relative to some normalization period mean
function normalize_data!(dat, norm_yrs=nothing)
    # normalize to relevant period  defined by norm_yrs
    idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
    norm_mean = map(mean, eachrow(dat[:, idx_norm]))
    for row in axes(dat, 1)
        foreach(col -> dat[row, col] -= norm_mean[row], axes(dat, 2))
    end
    return dat
end

function cumulative_sum(dat)
    cum_df = zeros(size(temp))
    for i in 1:ncol(dat)
        if i == 1
            cum_df[:, i] = dat[:, i]
        else
            cum_df[:, i] = dat[:, i] + cum_df[:, i-1]
        end
    end
    return cum_df
end

# find cumulative emissions from 2022--2100
idx2100 = findfirst(names(gmslr) .== "2100")
idx2000 = findfirst(names(gmslr) .== "2000")
idx1850 = findfirst(names(gmslr) .== "1850")
idx1900 = findfirst(names(gmslr) .== "1900")

normalize_data!(temperature, 1850:1900)
normalize_data!(gmslr, 1995:2014)

temp = temperature[:, idx2000:idx2100]
temp_integral = cumulative_sum(temp)



function fit_piecewise(dat, minbp, maxbp, step)

    function add_breakpoint(dat, bp)
        pred_var = names(dat)[1]
        dat[!, "after_bp"] = max.(0, dat[!, pred_var] .- bp)
        return (pred_var, dat)
    end

    min_deviance = Inf
    best_model = nothing
    best_bp = 0
    current_model = nothing
    
    for bp in minbp:step:maxbp
      pred_var, dat_bp = add_breakpoint(dat, bp)
      current_model = lm(@formula(slr ~ temp + after_bp), dat_bp)
      if deviance(current_model) < min_deviance
        min_deviance = deviance(current_model)
        best_model = current_model
        best_bp = bp
      end
    end
    
    return best_model, best_bp
  end

avg_temp_2100 = temperature[:, idx2100] - temperature[:, idx2000]  # find emissions average
temp_diff = temperature[:, idx2000+1:idx2100] .- temperature[:, idx2000]
temp_int = sum(Matrix(temp_diff), dims=2) ./ 100
#cum_emissions = reduce(hcat, map(cumsum, eachrow(emissions[:, idx2000+1:idx2100])))'
#emis_int = sum(Matrix(cum_emissions), dims=2) ./ 100

avg_gmslr_2100 = (gmslr[:, idx2100] - gmslr[:, idx2000]) * 1000 / (2100 - 2000 + 1)

gmslr_dat = DataFrame(temp=vec(temp_int), slr=avg_gmslr_2100)
#emis_dat = DataFrame(temp=vec(emis_int), slr=avg_gmslr_2100)


# check p-value for breakpoint models
function sim_pval_dist(data, null_model, min_bp, max_bp, step, n_runs)

    function simulate_data(null_model)
        resids = Normal(0, dispersion(null_model.model))
        pred_data = fitted(null_model)
        sim_data = pred_data + rand(resids, length(pred_data))
        return sim_data
    end

    function sim_pvalue(data, null_model, min_bp, max_bp, step)
        sim_dat = simulate_data(null_model)
        data[:, :slr] = sim_dat
        best_model = fit_piecewise(data, min_bp, max_bp, step)
        pval = coeftable(best_model[1]).cols[4][3]
        return pval
    end

    pvals = [sim_pvalue(data, null_model, min_bp, max_bp, step) for _ in 1:n_runs]
    return pvals
end

print("2100 Model...\n")

temp_lm_all = fit_piecewise(gmslr_dat, 0, 4, 0.05)
#emis_lm_all = fit_piecewise(emis_dat, 1000, 4000, 100)

# check p-values of simulated null (no-breakpoint) data to check for significance of breakpoint model
n_runs = 1_000
# use null (no-breakpoint model) to examine if significant breakpoint is found a null pattern
pvals_temp = sim_pval_dist(DataFrame(temp=gmslr_dat[!, :temp]), lm(@formula(slr ~ temp), gmslr_dat), 0, 4, 0.1, n_runs)
#pvals_emis = sim_pval_dist(DataFrame(temp=emis_dat[!, :temp]), lm(@formula(slr ~ temp), emis_dat), 1000, 4000, 100, n_runs)
# check proprotion of p-values lower than original fitted p-value
print(sum(pvals_temp .<= coeftable(temp_lm_all[1]).cols[4][3]) / n_runs)
#sum(pvals_emis .<= coeftable(emis_lm_all[1]).cols[4][3]) / n_runs
print("\n")
print("Integrated Model...\n")

gmslr_dat = DataFrame(temp=vec(temp_integral), slr=vec(Matrix(gmslr[:, idx2000:idx2100])))
temp_pred_range = round(minimum(gmslr_dat[:, 1]); digits=0):1:round(maximum(gmslr_dat[:, 1]); digits=0)

temp_lm_all = fit_piecewise(gmslr_dat, 0, 250, 10)
pvals_temp_int = sim_pval_dist(DataFrame(temp=gmslr_dat[!, :temp]), lm(@formula(slr ~ temp), gmslr_dat), 0, 4, 0.1, n_runs)
print(sum(pvals_temp_int .<= coeftable(temp_lm_all[1]).cols[4][3]) / n_runs)