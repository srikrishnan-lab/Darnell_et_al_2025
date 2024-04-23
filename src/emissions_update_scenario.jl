#########################################################################
# emissions_update_scenario.jl                                          #
#                                                                       #
# Samples from the emissions distribution based on the probability of   # 
#     CO2 emissions in 2100                                             #
#                                                                       #
#                                                                       #
# Note: It takes about 25 minutes to run 10,000 samples with this       #
#    script                                                             #
#########################################################################

# activate the environment
using Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

# load packages
using DataFrames
using CSVFiles
using XLSX
using Interpolations
using Distributions
using StatsBase
using AdaptiveMCMC
using MCMCChains
using Makie
using CairoMakie

include(joinpath(@__DIR__, "functions.jl")) # include functions from other scripts

n_samples = 500_000

# load CMIP scenarios to get 
## load CO2 emissions for  CMIP6 scenarios
cmip_df =  XLSX.readtable(joinpath(@__DIR__, "..", "data", "cmip6_co2.xlsx"), "data") |> DataFrame
# select only rows and columns with scenario names and data   
select!(cmip_df, Not([:Model, :Region, :Variable, :Unit, :Notes]))
cmip_df = cmip_df[1:7, :]
# reformat scenario names to SSPn-x.x
cmip_df[!, :Scenario] = replace.(cmip_df[!, :Scenario], r"(\d)(\d)" =>   s"\1.\2")
cmip_df[!, :Scenario] = [split(cmip_df[i, :Scenario], " ")[1] for i in 1:nrow(cmip_df)]
# sort by scenarios
sort!(cmip_df, :Scenario)
# convert emissions  from MtCO2/yr to GtCO2/yr
cmip_df[!, Not(:Scenario)] = cmip_df[!, Not(:Scenario)] ./ 1000
# function to compute the cumulative emissions from 2022-2100 from a linear interpolation of the SSP emissions data
ssp_yrs = parse.(Float64, names(cmip_df)[2:end])
function compute_cum_emissions_ssp(x, ssp_yrs)
    # linearly interpolate and compute sums
    itp_yrs = collect(2022:2100)
    itp = interpolate((ssp_yrs, ), Vector(x), Gridded(Linear()))
    return sum(itp[itp_yrs])
end
cum_ssp = DataFrame(scenario=cmip_df[!, :Scenario], emissions=map(d -> compute_cum_emissions_ssp(d, ssp_yrs), eachrow(cmip_df[!, Not(:Scenario)])))

# read in historical emissions observations (1850-2021)
historical_data = historical_emissions(start_year=1850, end_year=2100)

# define prior distributions for emissions parameters
# these are intentional quite wide to allow the cumulative emissions update to meaningfully update them
γ_g_pri    = truncated(Normal(0.002,0.15), 0.001, 0.05)  # growth parameter
t_peak_pri  = truncated(Normal(2070,20), 2030, 2200)       # peaking time
γ_d_pri     = truncated(Normal(0.002,0.15), 0.001, 0.25)     # decline parameter

# define scenarios
# these can be added to or modified for new scenarios

# default: just need to sample from the below distributions
function default_model(p)
    γ_g, t_peak, γ_d = p
    cb_2deg = 1220
    mean_cum = 3263
    # set priors
    lp = 0
    lp += logpdf(γ_g_pri, γ_g)
    lp += logpdf(t_peak_pri, t_peak)
    lp += logpdf(γ_d_pri, γ_d)

    if !isinf(lp)
        # simulate emissions trajectory
        t, gtco2 = emissions_curve(historical_data, γ_g=γ_g, t_peak=trunc(Int64, t_peak), γ_d=γ_d,end_year=2100) # years and emissions
        cum_emis = sum(gtco2[indexin(2022:2100, t)])
        lp += logpdf(MixtureModel([Uniform(0, cb_2deg), truncated(Normal(mean_cum, 2000), lower=cb_2deg, upper=8000)], [0.5, 0.5]), cum_emis)
    end
    return lp
end
out = adaptive_rwm([0.05, 2070, 0.005], default_model, 550_000; acc_sw=0.30,fulladapt=true, b=50_000)
c = Chains(out.X', [:γ_g, :t_peak, :γ_d]) 
df_default = DataFrame(c)[!, 3:5]
df_default[!, :t_peak] = trunc.(Int64, df_default[!, :t_peak])
save(joinpath(@__DIR__, "..", "data", "emissions", "default", "parameters.csv"), df_default)

# optimistic: a 50% probability of sticking to a carbon budget consistent with the Paris Agreement based on Lamboll et al (2023), linear probability decrease until RCP 8.5
# this is a 50% probability of cumulative emissions 
function optimistic_model(p)
    γ_g, t_peak, γ_d = p
    cb_2deg = 1220
    # set priors
    lp = 0
    lp += logpdf(γ_g_pri, γ_g)
    lp += logpdf(t_peak_pri, t_peak)
    lp += logpdf(γ_d_pri, γ_d)

    if !isinf(lp)
        # simulate emissions trajectory
        t, gtco2 = emissions_curve(historical_data, γ_g=γ_g, t_peak=trunc(Int64, t_peak), γ_d=γ_d,end_year=2100) # years and emissions
        cum_emis = sum(gtco2[indexin(2022:2100, t)])
        lp += logpdf(MixtureModel([Uniform(0, cb_2deg), truncated(Normal(cb_2deg, 1500), lower=cb_2deg, upper=8000)], [0.75, 0.25]), cum_emis)
    end
    return lp
end
out = adaptive_rwm([0.05, 2070, 0.005], optimistic_model, 550_000; acc_sw=0.30,fulladapt=true, b=50_000)
c = Chains(out.X', [:γ_g, :t_peak, :γ_d]) 
df_optimistic = DataFrame(c)[!, 3:5]
df_optimistic[!, :t_peak] = trunc.(Int64, df_optimistic[!, :t_peak])
save(joinpath(@__DIR__, "..", "data", "emissions", "optimistic", "parameters.csv"), df_optimistic)


# pessimistic: a 50% probability of sticking to a carbon budget consistent with the Paris Agreement based on Lamboll et al (2023), linear probability decrease until RCP 8.5
# this is a 50% probability of cumulative emissions 
function pessimistic_model(p)    
    γ_g, t_peak, γ_d = p
    cb_2deg = 1220
    mean_cum = 5231
    # set priors
    lp = 0
    lp += logpdf(γ_g_pri, γ_g)
    lp += logpdf(t_peak_pri, t_peak)
    lp += logpdf(γ_d_pri, γ_d)

    if !isinf(lp)
        # simulate emissions trajectory
        t, gtco2 = emissions_curve(historical_data, γ_g=γ_g, t_peak=trunc(Int64, t_peak), γ_d=γ_d,end_year=2100) # years and emissions
        cum_emis = sum(gtco2[indexin(2022:2100, t)])
        lp += logpdf(truncated(Normal(mean_cum, 2500), lower=cb_2deg, upper=8000), cum_emis)
    end
    return lp
end
out = adaptive_rwm([0.05, 2070, 0.005], pessimistic_model, 550_000; acc_sw=0.30,fulladapt=true, b=50_000)
c = Chains(out.X', [:γ_g, :t_peak, :γ_d]) 
df_pessimistic = DataFrame(c)[!, 3:5]
df_pessimistic[!, :t_peak] = trunc.(Int64, df_pessimistic[!, :t_peak])
save(joinpath(@__DIR__, "..", "data", "emissions", "pessimistic", "parameters.csv"), df_pessimistic)

# plot priors and scenario distributions
fig = Figure(size=(800, 400), fontsize=20, figure_padding=10)
ax1 = Axis(fig[1, 1], xlabel="γg", ylabel="Probability Density")
ax2 = Axis(fig[1, 2], xlabel="tpeak", ylabel="Probability Density")
ax3 = Axis(fig[1, 3], xlabel="γd", ylabel="Probability Density")

# first plot prior then scenario distributions
Makie.lines!(ax1, Makie.KernelDensity.kde(rand(γ_g_pri, n_samples),bandwidth=0.002), color=:gray, label="Prior", linewidth=3)
Makie.lines!(ax1, Makie.KernelDensity.kde(Vector(df_default[!, :γ_g]), bandwidth=0.002), color=:black, label="Baseline Scenario", linewidth=3)
Makie.lines!(ax1, Makie.KernelDensity.kde(Vector(df_optimistic[!, :γ_g]), bandwidth=0.002), color=:darkorange, label="Optimistic Scenario", linewidth=3)
Makie.lines!(ax1, Makie.KernelDensity.kde(Vector(df_pessimistic[!, :γ_g]), bandwidth=0.002), color=:teal, label="Pessimistic Scenario", linewidth=3)

Makie.lines!(ax2, Makie.KernelDensity.kde(rand(t_peak_pri, n_samples),bandwidth=2), color=:gray, label="Prior", linewidth=3)
Makie.lines!(ax2, Makie.KernelDensity.kde(Vector(df_default[!, :t_peak]), bandwidth=2), color=:black, label="Baseline Scenario", linewidth=3)
Makie.lines!(ax2, Makie.KernelDensity.kde(Vector(df_optimistic[!, :t_peak]), bandwidth=2), color=:darkorange, label="Optimistic Scenario", linewidth=3)
Makie.lines!(ax2, Makie.KernelDensity.kde(Vector(df_pessimistic[!, :t_peak]), bandwidth=2), color=:teal, label="Pessimistic Scenario", linewidth=3)

Makie.lines!(ax3, Makie.KernelDensity.kde(rand(γ_d_pri, n_samples),bandwidth=0.01), color=:gray, label="Prior", linewidth=3)
Makie.lines!(ax3, Makie.KernelDensity.kde(Vector(df_default[!, :γ_d]), bandwidth=0.01), color=:black, label="Baseline Scenario", linewidth=3)
Makie.lines!(ax3, Makie.KernelDensity.kde(Vector(df_optimistic[!, :γ_d]), bandwidth=0.01), color=:darkorange, label="Optimistic Scenario", linewidth=3)
Makie.lines!(ax3, Makie.KernelDensity.kde(Vector(df_pessimistic[!, :γ_d]), bandwidth=0.01), color=:teal, label="Pessimistic Scenario", linewidth=3)

Legend(fig[2, 1:3], ax1, orientation=:horizontal)

CairoMakie.save(joinpath(@__DIR__, "..", "figures", "emissions_updates.png"), fig)