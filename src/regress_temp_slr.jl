# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSVs
using DataFrames # data structure for indices
using Makie # plotting library
using CairoMakie
using Measures # adjust margins with explicit measures
using StatsBase # get mean function and density
using Turing
using Distributions

@model function breakpoint_lr(X, Y)
    a₁ ~ Normal(6, 0.5)

    b ~ Normal(0, 1)
    σ² ~ truncated(Normal(0, 1); lower=0)
    μ = b .+  a₁ * X
    Y .~ arraydist([Normal(μ[i], σ²) for i in eachindex(X)])
end

function predict_lr(X, chain)
    a₁ = chain[:a₁].data
    b = chain[:b].data
    σ² = chain[:σ²].data 
    pred = zeros(length(X), size(chain)[1])
    for i in 1:size(chain)[1]
        err = rand(Normal(0, σ²[i]), length(X))
        pred[:, i] = a₁[i] * X .+ b[i] .+ err
    end
    return pred
end

mod = breakpoint_lr(avg_temp[avg_ais .< 0.25], avg_ais[avg_ais .< 0.25])
chain = sample(mod, NUTS(0.65), 1_000)

temp_test = 0:0.1:4
pred = predict_lr(temp_test, chain)
pred_q = mapslices(col -> quantile(col, [0.05, 0.5, 0.95]), pred, dims=2)
Plots.scatter(avg_temp, avg_ais, alpha=0.1, color=:lightblue)
Plots.plot!(temp_test, pred_q[:, 2], color=:red, ribbon=(pred_q[:, 2] .- pred_q[:, 1], pred_q[:, 3] .- pred_q[:, 1]))