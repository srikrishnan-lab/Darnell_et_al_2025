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
function compute_cum_emissions_ssp(x, ssp_yrs)
    # linearly interpolate and compute sums
    itp_yrs = collect(2022:2100)
    itp = interpolate((ssp_yrs, ), Vector(x), Gridded(Linear()))
    return sum(itp[itp_yrs])
end
cum_ssp = DataFrame(scenario=cmip_df[!, :Scenario], emissions=map(d -> compute_cum_emissions_ssp(d, ssp_yrs), eachrow(cmip_df[!, Not(:Scenario)])))

# read in historical emissions observations (1850-2021)
historical_data = historical_emissions(start_year=1850, end_year=2100)

# read in scenarios
df_default = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "data", "emissions", "default", "parameters.csv")))
df_optimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "data", "emissions", "optimistic", "parameters.csv")))
df_pessimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "data", "emissions", "pessimistic", "parameters.csv")))

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

function calc_cum_emissions(df, t, start_year, end_year)
    cum_emis = [sum(r[indexin(start_year:end_year, t)]) for r in eachrow(df)]
    return cum_emis
end

# get emissions trajectories
t, gtco2_default = sim_emissions(df_default, 2200)
t, gtco2_optimistic = sim_emissions(df_optimistic, 2200)
t, gtco2_pessimistic = sim_emissions(df_pessimistic, 2200)

q_default =  mapslices(col -> quantile(col, [0.05, 0.5, 0.95]), gtco2_default; dims=1)
q_optimistic =  mapslices(col -> quantile(col, [0.05, 0.5, 0.95]), gtco2_optimistic; dims=1)
q_pessimistic =  mapslices(col -> quantile(col, [0.05, 0.5, 0.95]), gtco2_pessimistic; dims=1)


# get cumulative emissions
cum_default = calc_cum_emissions(gtco2_default, t, 2022, 2100)
cum_optimistic =  calc_cum_emissions(gtco2_optimistic, t, 2022, 2100)
cum_pessimistic = calc_cum_emissions(gtco2_pessimistic, t, 2022, 2100)

# plot resulting cumulative distributions
fig = Figure(size=(800, 600), fontsize=16, figure_padding=20)
gemissions = fig[1:3, 1:2] = GridLayout()

axpdf = Axis(gemissions[1, 1], xlabel="Cumulative CO₂ Emissions, 2022--2018 (GtCO₂)", ylabel="Density")
Makie.lines!(axpdf, Makie.KernelDensity.kde(cum_default, bandwidth=200), color=:black, linewidth=3, label="Default Scenario")
Makie.lines!(axpdf, Makie.KernelDensity.kde(cum_optimistic, bandwidth=200), color=:darkorange, linewidth=3, label="Optimistic Scenario")
Makie.lines!(axpdf, Makie.KernelDensity.kde(cum_pessimistic, bandwidth=200), color=:teal, linewidth=3, label="Pessimistic Scenario")
Makie.vlines!(axpdf, [1220], linestyle=:dash, label="2°C Budget (50%)", color=:green, linewidth=2)
cmip_colors = cgrad(:Dark2_7, 7, categorical=true)
for i = 1:nrow(cum_ssp)
    Makie.vlines!(axpdf, [cum_ssp[i, :emissions]], color=cmip_colors[i], linestyle=:dot, label=cum_ssp[i, :scenario], linewidth=2)
end

default_cdf = ecdf(cum_default)
default_range = 0:maximum(cum_default)
optimistic_cdf = ecdf(cum_optimistic)
optimistic_range = 0:maximum(cum_optimistic)
pessimistic_cdf = ecdf(cum_pessimistic)
pessimistic_range = 0:maximum(cum_pessimistic)
axcdf = Axis(gemissions[1, 2], xlabel="Cumulative CO₂ Emissions, 2022--2018 (GtCO₂)", ylabel="Cumulative Probability")
Makie.lines!(axcdf, default_range, default_cdf.(default_range), label="Baseline Scenario", color=:black, linewidth=3)
Makie.lines!(axcdf, optimistic_range, optimistic_cdf.(optimistic_range), label="Optimistic Scenario", color=:darkorange, linewidth=3)
Makie.lines!(axcdf, pessimistic_range, pessimistic_cdf.(pessimistic_range), label="Pessimistic Scenario", color=:teal, linewidth=3)
Makie.vlines!(axcdf, [1220], linestyle=:dash, label="2°C Budget (50%)", color=:green, linewidth=2)
for i = 1:nrow(cum_ssp)
    Makie.vlines!(axcdf, [cum_ssp[i, :emissions]], color=cmip_colors[i], linestyle=:dot, label=cum_ssp[i, :scenario], linewidth=2)
end

# Time trajectories of emissions
axemissions = Axis(gemissions[2, 1], xlabel="Year", ylabel="CO₂ Emissions (Gt CO₂/yr)")
leg_med = Makie.lines!(axemissions, 2000:2200, q_default[2, indexin(2000:2200, t)], color=:black, linewidth=3)
leg_ci = Makie.band!(axemissions, 2000:2200, q_default[1, indexin(2000:2200, t)], q_default[3, indexin(2000:2200, t)], color=(:black, 0.1))
Makie.lines!(axemissions, 2000:2200, q_optimistic[2, indexin(2000:2200, t)], color=:darkorange2, linewidth=3)
Makie.band!(axemissions, 2000:2200, q_optimistic[1, indexin(2000:2200, t)], q_optimistic[3, indexin(2000:2200, t)], color=(:darkorange2, 0.1))
Makie.lines!(axemissions, 2000:2200, q_pessimistic[2, indexin(2000:2200, t)], color=:teal, linewidth=3)
Makie.band!(axemissions, 2000:2200, q_pessimistic[1, indexin(2000:2200, t)], q_pessimistic[3, indexin(2000:2200, t)], color=(:teal, 0.1))
# add in the SSP emissions scenarios for context
cmip_yrs = parse.(Float64, string.(names(cmip_df[!, Not(:Scenario)])))
cmip_colors = cgrad(:Dark2_7, 7, categorical=true)
for i = 1:nrow(cmip_df)
    Makie.lines!(axemissions, cmip_yrs, Vector(cmip_df[i, Not(:Scenario)]), color=cmip_colors[i], linestyle=:dash, linewidth=2)
end
axislegend(axemissions, [leg_med, leg_ci], ["Median", "90% Interval"], position=:rt)

default_cdf = ecdf(gtco2_default[:, findfirst(==(2100), t)])
default_range = 0:maximum(gtco2_default[:, findfirst(==(2100), t)])
optimistic_cdf = ecdf(gtco2_optimistic[:, findfirst(==(2100), t)])
optimistic_range = 0:maximum(gtco2_optimistic[:, findfirst(==(2100), t)])
pessimistic_cdf = ecdf(gtco2_pessimistic[:, findfirst(==(2100), t)])
pessimistic_range = 0:maximum(gtco2_pessimistic[:, findfirst(==(2100), t)])
ax2100 = Axis(gemissions[2, 2], xlabel="CO₂ Emissions in 2100 (GtCO₂/yr)", ylabel="Cumulative Probability")
Makie.lines!(ax2100, default_range, default_cdf.(default_range), label="Baseline Scenario", color=:black, linewidth=3)
Makie.lines!(ax2100, optimistic_range, optimistic_cdf.(optimistic_range), label="Optimistic Scenario", color=:darkorange, linewidth=3)
Makie.lines!(ax2100, pessimistic_range, pessimistic_cdf.(pessimistic_range), label="Pessimistic Scenario", color=:teal, linewidth=3)
for i = 1:nrow(cmip_df)
    Makie.vlines!(ax2100, [cmip_df[i, end]], color=cmip_colors[i], linestyle=:dot, label=cmip_df[i, end], linewidth=2)
end


rowgap!(gemissions, 1, Relative(0.15))

Label(fig[1, 1, TopLeft()], "a", fontsize=20, font=:bold, padding = (0, 35, 10, 0), halign=:right)
Label(fig[1, 2, TopLeft()], "b", fontsize=20, font=:bold, padding = (0, 35, 10, 0), halign=:right)
Label(fig[2, 1, TopLeft()], "c", fontsize=20, font=:bold, padding = (0, 35, -40, 0), halign=:right)
Label(fig[2, 2, TopLeft()], "d", fontsize=20, font=:bold, padding = (0, 35, -40, 0), halign=:right)

# add legend
leg = Legend(gemissions[3, 1:2], axcdf, orientation=:horizontal, framevisible=false, tellheight=true, titleposition=:top, tellwidth=false, nbanks=3)

CairoMakie.save(joinpath(@__DIR__, "..", "figures", "scenario-emissions.png"), fig)
