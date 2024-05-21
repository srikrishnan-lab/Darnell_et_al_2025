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
using Distributions
using StatsBase
using Optim
using Makie
using CairoMakie
using LaTeXStrings
using Measures

include(joinpath(@__DIR__, "functions.jl")) # include functions from other scripts

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

# read in historical emissions observations (1850-2021)
historical_data = historical_emissions(start_year=1850, end_year=2100)

function sim_emissions(df, end_year)
    γ_g = df[!, :γ_g]
    γ_d = df[!, :γ_d]
    t_peak = df[!, :t_peak]
    t = zeros(end_year - 1850 + 1)
    gtco2 = zeros((nrow(df), end_year - 1850 + 1))
    for i in eachindex(γ_g)
        t[:], gtco2[i, :] = emissions_curve(historical_data, γ_g=γ_g[i], t_peak=trunc(Int64, t_peak[i]), γ_d=γ_d[i],end_year=end_year) # years and emissions
    end
    return t, gtco2
end

# develop samples for illustration
t_peak = [2050, 2070, 2090]
γ_g = [0.005, 0.01, 0.03]
γ_d = [0.01, 0.04, 0.1]

emis_samp = DataFrame((t_peak=t, γ_g = x, γ_d = y) for t in t_peak for x in γ_g for y in γ_d)

# get emissions trajectories
t_sim, emis_sim =  sim_emissions(emis_samp, 2200)

## fit trajectories to SSP scenarios
function ssp_fit(p, ssp, ssp_yrs)
    params = DataFrame(γ_g = p[1], t_peak=trunc.(Int64, p[2]), γ_d = p[3])
    t, emis_sim =  sim_emissions(params, 2100)
    yr_idx = indexin(ssp_yrs, t)
    mse = mean((emis_sim[yr_idx] .- ssp).^2)
    return mse
end

bounds = [0.0 2020.0 0.0;
          0.5 2300.0 0.5]
sim_out = zeros(size(cmip_df, 1), 2100 - 1850 + 1)
pfit = zeros(size(cmip_df, 1), 3)
for idx in axes(cmip_df, 1)
      result = Metaheuristics.optimize(p -> ssp_fit(p, Vector(cmip_df[idx, 2:end]), ssp_yrs), bounds, DE())
      pfit[idx, :] = minimizer(result)
      params_ssp = DataFrame(γ_g = pfit[idx, 1], t_peak=trunc.(Int64, pfit[idx, 2]), γ_d = pfit[idx, 3])
      t, sim_out[idx, :] = sim_emissions(params_ssp, 2100)
end
t_ssp = 1850:2100

# plot resulting cumulative distributions
fig = Figure(size=(800, 400), fontsize=16, figure_padding=20)

axpdf = Axis(fig[1, 1], xlabel="Year", ylabel="Annual CO₂ Emissions (GtCO₂/yr)", limits=(2020, 2200, nothing, 130), xticks=2020:40:2200)
leg_sim = let x
    for idx in axes(emis_sim, 1)
        x =  Makie.lines!(axpdf, t_sim, emis_sim[idx, :], label="Simulated Trajectory", color=:black, alpha=0.2)
    end
    x
end
axislegend(axpdf, [leg_sim], ["Simulated Trajectory"], position=:rt)

axssp = Axis(fig[1, 2], xlabel="Year", ylabel="Annual CO₂ Emissions (GtCO₂/yr)", limits=(2020, 2100, nothing, nothing))
cmip_colors = cgrad(:Dark2_7, 7, categorical=true)
for i in axes(cmip_df, 1)
    Makie.lines!(axssp, t_ssp, sim_out[i, :], color=cmip_colors[i], linewidth=2)
    Makie.lines!(axssp, ssp_yrs, Vector(cmip_df[i, 2:end]), color=cmip_colors[i], linestyle=:dot, linewidth=2, label=cmip_df[i, :Scenario])
end

# add legend
leg = Legend(fig[1, 3], axssp, "SSP Scenarios", framevisible=false, tellwidth=true, titleposition=:top, tellheight=false)

Label(fig[1, 1, TopLeft()], "a", fontsize=20, font=:bold, padding = (0, 35, 10, 0), halign=:right)
Label(fig[1, 2, TopLeft()], "b", fontsize=20, font=:bold, padding = (0, 35, 10, 0), halign=:right)


CairoMakie.save(joinpath(@__DIR__, "..", "figures", "sample-emissions.png"), fig)
