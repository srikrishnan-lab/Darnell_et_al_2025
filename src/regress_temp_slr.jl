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

@model function breakpoint_lr(Y, X)
    a₁ ~ Normal(0, 0.1)
    a₂ ~ truncated(Normal(0, 0.05); lower=0)

    b ~ Normal(0, 0.1)
    ν ~ Uniform(1, 2)
    μ = b .+ ifelse.(X .> ν, a₁ .+ a₂ * (X .- ν), a₁) .* X
    for i in eachindex(Y)
      Y[i] ~ Cauchy(μ[i], 1)
    end
end

mod = breakpoint_lr(avg_gmslr, avg_temp)
chain = sample(mod, NUTS(), 100)
