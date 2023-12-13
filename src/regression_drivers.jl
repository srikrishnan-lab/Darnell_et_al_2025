import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Random
using CSVFiles
using DataFrames
using Statistics
using MLJ
using EvoTrees
using ShapML

Random.seed!(1)

EvoTreeRegressor = @load EvoTreeRegressor pkg=EvoTrees

# Assume that we call this script from the project folder
output_dir = "results/default"
slr_out = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))

# normalize relative to 2000
idx_2000 = findfirst(names(slr_out) .== "2000")
for row in axes(slr_out,1 )
    foreach(col -> slr_out[row, col] -= slr_out[row, idx_2000], axes(slr_out, 2))
end

# sample random features from the overall set
idx = rand(1:size(parameters, 1), 10000)
features = parameters[idx, :]
targets = slr_out[idx, :]

# define function for parallelized Shapley calculation
function predict_slr(model, data)
    pred = DataFrame(slr_pred = MLJ.predict(model,data))
    return pred
end

# make regression trees and compute Shapley values
function shapley_reg(yrs, features, targets)
    shap_df = DataFrame(feature_name = names(features))
    for yr in yrs
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

yrs = 2050:5:2200
shap_df = shapley_reg(yrs, features, targets)
save("output/shapley/shapley_indices.csv", shap_df)