#########################################################################
# regression_drivers.jl                                                 #
#                                                                       #
# Makes plot summarizing the output ensemble.                           # #                                                                       #
#                                                                       #
# This script requires the ensemble output to be present in             #
#   `results/default` .                                                 #
#                                                                       #
#########################################################################

# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Random # set seed (since this estimator is Monte Carlo-based) and sample
using CSVFiles # read ensemble output
using DataFrames # tabular data structure
using Statistics # compute mean
using MLJ
using EvoTrees
using ShapML
using StatsBase

Random.seed!(1)

EvoTreeRegressor = @load EvoTreeRegressor pkg=EvoTrees

# Assume that we call this script from the project folder
output_dir = joinpath(@__DIR__, "..", "results", "default")
slr_out = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))
ar6_out = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "output", "ar6", "gmslr.csv")))

# normalize relative to 2000
# function to find normalize relative to a year or mean over some normalization period)
function normalize_data!(dat, norm_yrs=nothing)
    # normalize to relevant period  defined by norm_yrs
    idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
    norm_mean = map(mean, eachrow(dat[:, idx_norm]))
    for row in axes(dat, 1)
        foreach(col -> dat[row, col] -= norm_mean[row], axes(dat, 2))
    end
    return dat
end

# normalize BRICK SLR Output
normalize_data!(slr_out, 1995:2014)

# sample random features from the overall set
idx = rand(1:size(parameters, 1), 10_000)
features = parameters[idx, :]
targets_brick = slr_out[idx, :]
yrs = 2050:10:2150
targets_ar6 = zeros(10_000, length(yrs))
targets_ar6 = DataFrame(targets_ar6, string.(yrs))
for yr in yrs
    brick_cdf = ecdf(slr_out[:, Symbol(yr)])
    targets_q = brick_cdf.(targets_brick[:, Symbol(yr)])
    targets_ar6[:, Symbol(yr)] = quantile(ar6_out[:, Symbol(yr)], targets_q)
end

# define function for parallelized Shapley calculation
function predict_slr(model, data)
    pred = DataFrame(slr_pred = MLJ.predict(model,data))
    return pred
end

# make regression trees and compute Shapley values
function shapley_reg(yrs, features, targets)
    shap_df = DataFrame(feature_name = names(features))
    for yr in yrs
        println("$yr")
        slr_reg_tree = EvoTreeRegressor(nrounds=200, max_depth=5)
        slr_reg_mach = machine(slr_reg_tree, features, targets[:, Symbol(yr)])
        MLJ.fit!(slr_reg_mach, force=true)

        explain = copy(features)
        reference = copy(features)
        shap_out = ShapML.shap(explain = explain, 
                                reference = reference,
                                model = slr_reg_mach,
                                predict_function = predict_slr,
                                sample_size = 100,
                                seed = 1)
        shap_grouped = groupby(shap_out, :feature_name)
        shap_summary = combine(shap_grouped, 
            :shap_effect => (x -> mean(abs.(x))))
        rename!(shap_summary, Dict(:shap_effect_function => Symbol("mean_$(yr)")))
        shap_df = innerjoin(shap_df, shap_summary, on=:feature_name)
    end
    return shap_df
end

shap_df = shapley_reg(yrs, features, targets_ar6)
save(joinpath(@__DIR__, "..", "output", "shapley", "shapley_indices_ar6.csv"), shap_df)