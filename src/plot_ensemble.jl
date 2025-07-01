#########################################################################
# plot_ensemble.jl                                                      #
#                                                                       #
# Makes plot summarizing series of ensemble projections.                # #                                                                       #
#                                                                       #
# This script requires the ensemble output to be present in             #
#   `results/default` .                                                 #
#                                                                       #
#########################################################################

# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSV of Shapley indices
using XLSX
using DataFrames # data structure for indices
using Makie # plotting library
using CairoMakie
using Measures # adjust margins with explicit measures
using StatsBase # get mean

# load the results
# read in scenarios
params_default = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "default", "parameters.csv")))
params_optimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "optimistic", "parameters.csv")))
params_pessimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "pessimistic", "parameters.csv")))
emis_default = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "default", "emissions.csv")))
emis_optimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "optimistic", "emissions.csv")))
emis_pessimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "pessimistic", "emissions.csv")))
temp_default = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "default", "temperature.csv")))
temp_optimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "optimistic", "temperature.csv")))
temp_pessimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "pessimistic", "temperature.csv")))
gmslr_default = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "default", "gmslr.csv")))
gmslr_optimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "optimistic", "gmslr.csv")))
gmslr_pessimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "pessimistic", "gmslr.csv")))
ais_default = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "default", "antarctic.csv")))
ais_optimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "optimistic", "antarctic.csv")))
ais_pessimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "pessimistic", "antarctic.csv")))
gis_default = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "default", "greenland.csv")))
gis_optimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "optimistic", "greenland.csv")))
gis_pessimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "pessimistic", "greenland.csv")))

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
cmip_df = cmip_df[Not(cmip_df[:, :Scenario] .== "SSP4-3.4"), :] # drop SSP4-3.4 since it doesn't show up in AR6
cmip_df = cmip_df[Not(cmip_df[:, :Scenario] .== "SSP4-6.0"), :] # drop SSP4-6.0 since it doesn't show up in AR6

ar6_lines_927 = Dict("SSP1-1.9" => colorant"rgb( 35, 197, 226)",
                    "SSP1-2.6" => colorant"rgb( 21, 73, 135)",
                    "SSP2-4.5" => colorant"rgb(255, 150,  60)",
                    "SSP3-7.0" => colorant"rgb(255,  66,  48)",
                    "SSP5-8.5" => colorant"rgb(178,  44,  26)",
                    "SSP5-8.5 LC" => colorant"rgb(200, 35, 200)"
)
ar6_temp = DataFrame(
    year = repeat([2050, 2090], inner=5),
    scenario = repeat(["SSP1-1.9", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"], outer=2),
    min = [1.2, 1.3, 1.6, 1.7, 1.9, 1.0, 1.3, 2.1, 2.8, 3.3],
    median = [1.6, 1.7, 2.0, 2.1, 2.4, 1.4, 1.8, 2.7, 3.6, 4.4],
    max = [2.0, 2.2, 2.5, 2.6, 3.0, 1.8, 2.4, 3.5, 4.6, 5.7]
)
# add colors to AR6 dataframe
cmip_df.color = [ar6_lines_927[s] for s in cmip_df.Scenario]
ar6_temp.color = [ar6_lines_927[s] for s in ar6_temp.scenario]



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

function compute_ais_fd(temps, params)
    ais_fd_trigger = 15.42 .+ 0.8365 * params[:, :antarctic_temp_threshold]
    temp_mat = Matrix(temps)
    ais_fd = zeros(size(temps))
    for i in 1:size(temp_mat, 1)
        ais_fd[i, :] = (temp_mat[i, :] .> ais_fd_trigger[i]) * params[i, :antarctic_lambda]
    end
    ais_fd_cum = mapslices(cumsum, ais_fd; dims=2)
    return ais_fd_cum
end

# normalize data
normalize_data!(temp_default, 1850:1900)
normalize_data!(temp_optimistic, 1850:1900)
normalize_data!(temp_pessimistic, 1850:1900)
normalize_data!(gmslr_default, 1995:2014)
normalize_data!(gmslr_optimistic, 1995:2014)
normalize_data!(gmslr_pessimistic, 1995:2014)
normalize_data!(gis_default, 1995:2014)
normalize_data!(gis_optimistic, 1995:2014)
normalize_data!(gis_pessimistic, 1995:2014)

# compute quantiles
emis_q_default =  mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), emis_default)
emis_q_optimistic =  mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), emis_optimistic)
emis_q_pessimistic =  mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), emis_pessimistic)
temp_q_default =  mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), temp_default)
temp_q_optimistic =  mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), temp_optimistic)
temp_q_pessimistic =  mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), temp_pessimistic)

gmslr_default_stack = DataFrame(stack(gmslr_default[:, [:"2050", :"2100", :"2150", :"2200"]]))
gmslr_default_stack.scenario .= "Baseline"
gmslr_optimistic_stack = DataFrame(stack(gmslr_optimistic[:, [:"2050", :"2100", :"2150", :"2200"]]))
gmslr_optimistic_stack.scenario .= "Optimistic"
gmslr_pessimistic_stack = DataFrame(stack(gmslr_pessimistic[:, [:"2050", :"2100", :"2150", :"2200"]]))
gmslr_pessimistic_stack.scenario .= "Pessimistic"
gmslr_all = vcat(gmslr_default_stack, gmslr_optimistic_stack, gmslr_pessimistic_stack)

gis_default_stack = DataFrame(stack(gis_default[:, [:"2050", :"2100", :"2150", :"2200"]]))
gis_default_stack.scenario .= "Baseline"
gis_optimistic_stack = DataFrame(stack(gis_optimistic[:, [:"2050", :"2100", :"2150", :"2200"]]))
gis_optimistic_stack.scenario .= "Optimistic"
gis_pessimistic_stack = DataFrame(stack(gis_pessimistic[:, [:"2050", :"2100", :"2150", :"2200"]]))
gis_pessimistic_stack.scenario .= "Pessimistic"
gis_all = vcat(gis_default_stack, gis_optimistic_stack, gis_pessimistic_stack)

ais_default_stack = DataFrame(stack(ais_default[:, [:"2050", :"2100", :"2150", :"2200"]]))
ais_default_stack.scenario .= "Baseline"
ais_optimistic_stack = DataFrame(stack(ais_optimistic[:, [:"2050", :"2100", :"2150", :"2200"]]))
ais_optimistic_stack.scenario .= "Optimistic"
ais_pessimistic_stack = DataFrame(stack(ais_pessimistic[:, [:"2050", :"2100", :"2150", :"2200"]]))
ais_pessimistic_stack.scenario .= "Pessimistic"
ais_all = vcat(ais_default_stack, ais_optimistic_stack, ais_pessimistic_stack)

# compute AIS fast dynamics contributions quantiles
fd_default = DataFrame(compute_ais_fd(temp_default, params_default), string.(1850:2300))
fd_optimistic = DataFrame(compute_ais_fd(temp_optimistic, params_optimistic), string.(1850:2300))
fd_pessimistic = DataFrame(compute_ais_fd(temp_pessimistic, params_pessimistic), string.(1850:2300))

fd_default_stack = DataFrame(stack(fd_default[:, [:"2050", :"2100", :"2150", :"2200"]]))
fd_default_stack.scenario .= "Baseline"
fd_optimistic_stack = DataFrame(stack(fd_optimistic[:, [:"2050", :"2100", :"2150", :"2200"]]))
fd_optimistic_stack.scenario .= "Optimistic"
fd_pessimistic_stack = DataFrame(stack(fd_pessimistic[:, [:"2050", :"2100", :"2150", :"2200"]]))
fd_pessimistic_stack.scenario .= "Pessimistic"
fd_all = vcat(fd_default_stack, fd_optimistic_stack, fd_pessimistic_stack)

# compute years of fast dynamics triggers
function find_fd_year(fd, years)
    fd_yr= Vector{Union{Int64, Nothing}}(undef, nrow(fd)) 
    fd_mat = Matrix(fd)
    for i in 1:nrow(fd_default)
        fd_yr_idx = findfirst((.>)(0), fd_mat[i, :])
        if !isnothing(fd_yr_idx)
            fd_yr[i] = years[fd_yr_idx]
        end
    end
    return fd_yr
end

fd_yr_default = identity.(filter(x -> !isnothing(x), find_fd_year(fd_default, 1850:2300)))
fd_yr_optimistic = identity.(filter(x -> !isnothing(x), find_fd_year(fd_optimistic, 1850:2300)))
fd_yr_pessimistic = identity.(filter(x -> !isnothing(x), find_fd_year(fd_pessimistic, 1850:2300)))

## make plot
inch = 96
mmx = inch / 25.4
fig = Figure(size=(180mmx, 120mmx), fontsize=9, figure_padding=10)

gemissions = fig[1:2, 1] = GridLayout()
gslr = fig[1:2, 2:3] = GridLayout()

# plot 1: emissions distribution
axemissions = Axis(gemissions[1, 1], xlabel="Year", ylabel="CO₂ Emissions (GtCO₂/yr)", xgridvisible=false, ygridvisible=false)
Makie.hlines!(axemissions, [0], color=:black, linewidth=1, alpha=0.8)
leg_med = Makie.lines!(axemissions, parse.(Int64, names(emis_q_default)), Vector(emis_q_default[2, :]), color=:darkgrey, linewidth=2)
leg_ci = Makie.band!(axemissions, parse.(Int64, names(emis_q_default)), Vector(emis_q_default[1, :]), Vector(emis_q_default[3, :]), color=(:darkgrey, 0.2))
Makie.lines!(axemissions, parse.(Int64, names(emis_q_optimistic)), Vector(emis_q_optimistic[2, :]), color=:darkorange, linewidth=2)
Makie.band!(axemissions, parse.(Int64, names(emis_q_optimistic)), Vector(emis_q_optimistic[1, :]), Vector(emis_q_optimistic[3, :]), color=(:darkorange, 0.2))
Makie.lines!(axemissions, parse.(Int64, names(emis_q_pessimistic)), Vector(emis_q_pessimistic[2, :]), color=:teal, linewidth=2)
Makie.band!(axemissions, parse.(Int64, names(emis_q_pessimistic)), Vector(emis_q_pessimistic[1, :]), Vector(emis_q_pessimistic[3, :]), color=(:teal, 0.2))
Makie.xlims!(axemissions, 2000, 2200)

# add in the SSP emissions scenarios for context
cmip_legend = Vector{LineElement}(undef, nrow(cmip_df)) # initialize storage for legend
cmip_yrs = parse.(Float64, string.(names(cmip_df[!, Not([:Scenario, :color])])))
for i = 1:nrow(cmip_df)
    Makie.lines!(axemissions, cmip_yrs, Vector(cmip_df[i, Not([:Scenario, :color])]), color=cmip_df.color[i], linestyle=:dash, linewidth=2)
    cmip_legend[i] = LineElement(color = cmip_df.color[i], linestyle = :dash, linewidth=4)
end

# plot 2: Temperature distribution
axtemp = Axis(gemissions[2, 1], ylabel="GMT Anomaly (°C)", xlabel="Year", xgridvisible=false, ygridvisible=false)
Makie.lines!(axtemp, parse.(Int64, names(temp_q_default)), Vector(temp_q_default[2, :]), color=:darkgrey, linewidth=2)
Makie.band!(axtemp, parse.(Int64, names(temp_q_default)), Vector(temp_q_default[1, :]), Vector(temp_q_default[3, :]), color=(:darkgrey, 0.2))
Makie.lines!(axtemp, parse.(Int64, names(temp_q_optimistic)), Vector(temp_q_optimistic[2, :]), color=:darkorange, linewidth=2)
Makie.band!(axtemp, parse.(Int64, names(temp_q_optimistic)), Vector(temp_q_optimistic[1, :]), Vector(temp_q_optimistic[3, :]), color=(:darkorange, 0.2))
Makie.lines!(axtemp, parse.(Int64, names(temp_q_pessimistic)), Vector(temp_q_pessimistic[2, :]), color=:teal, linewidth=2)
Makie.band!(axtemp, parse.(Int64, names(temp_q_pessimistic)), Vector(temp_q_pessimistic[1, :]), Vector(temp_q_pessimistic[3, :]), color=(:teal, 0.2))

Legend(gemissions[3, 1], [[LineElement(color=:black), PolyElement(color=(:black, 0.2))], cmip_legend], [["Median", "90% Interval"], cmip_df[!, :Scenario]], ["Ensemble Summary", "SSP Scenario"], vertical=false, framevisible=false, tellheight=true, nbanks=3, titleposition=:top, tellwidth=false)


# add dodge values for ar6 temperatures
plt_dodge = -10.0:5.0:10.0
ar6_temp.plt_x = ar6_temp.year + repeat(plt_dodge, outer=2)
Makie.rangebars!(axtemp, ar6_temp.plt_x,  ar6_temp.min, ar6_temp.max, color=ar6_temp.color, linewidth=1)
Makie.scatter!(axtemp, ar6_temp.plt_x, ar6_temp.median, color=ar6_temp.color, markersize=5)

Makie.xlims!(axtemp, 2000, 2200)
Makie.ylims!(axtemp, 0, 6)

Label(gemissions[1, 1, TopLeft()], "a", fontsize=10, font=:bold, padding = (10, 50, 20, 0), halign=:right)
Label(gemissions[2, 1, TopLeft()], "b", fontsize=10, font=:bold, padding = (10, 50, 20, 0), halign=:right)

dodge = identity.(indexin(ais_all.scenario, ["Optimistic", "Baseline", "Pessimistic"]))
color = [:darkorange, :darkgrey, :teal][dodge]

axgmsl = Axis(gslr[1, 1], ylabel="GMSL Anomaly (m)", xlabel="Year", xticks=(1:4, ["2050", "2100", "2150", "2200"]), yminorticks=IntervalsBetween(2), yminorticksvisible = true, xgridvisible=false, ygridvisible=false)
Makie.hlines!(axgmsl, [0], color=:black, alpha=0.8, linewidth=1)
Makie.boxplot!(axgmsl, identity.(indexin(gmslr_all.variable, ["2050", "2100", "2150", "2200"])), Vector(gmslr_all.value), dodge=dodge, color=color, markersize=5)

axgis = Axis(gslr[1, 2], ylabel="GMSLR from GIS (m SLR-eq)", xlabel="Year", xticks=(1:4, ["2050", "2100", "2150", "2200"]), yminorticks=IntervalsBetween(2), yminorticksvisible = true, xgridvisible=false, ygridvisible=false) 
Makie.hlines!(axgis, [0], color=:black, alpha=0.8, linewidth=1)
Makie.boxplot!(axgis, identity.(indexin(gis_all.variable, ["2050", "2100", "2150", "2200"])), Vector(gis_all.value), dodge=dodge, color=color, markersize=5)

axais = Axis(gslr[2, 1], ylabel="GMSLR from AIS (m SLR-eq)", xlabel="Year", xticks=(1:4, ["2050", "2100", "2150", "2200"]), yminorticks=IntervalsBetween(2), yminorticksvisible = true, xgridvisible=false, ygridvisible=false) 
Makie.hlines!(axais, [0], color=:black, alpha=0.8, linewidth=1)
Makie.boxplot!(axais, identity.(indexin(ais_all.variable, ["2050", "2100", "2150", "2200"])), Vector(ais_all.value), dodge=dodge, color=color, markersize=5)
Makie.ylims!(axais, -1.5, 4.75)

axaisnofd = Axis(gslr[2, 2], ylabel="Non-FD GMSLR from AIS (m SLR-eq)", xlabel="Year", xticks=(1:4, ["2050", "2100", "2150", "2200"]), yminorticks=IntervalsBetween(2), yminorticksvisible = true, xgridvisible=false, ygridvisible=false) 
Makie.hlines!(axaisnofd, [0], color=:black, alpha=0.8, linewidth=1)
Makie.boxplot!(axaisnofd, identity.(indexin(ais_all.variable, ["2050", "2100", "2150", "2200"])), Vector(ais_all.value - fd_all.value), dodge=dodge, color=color, markersize=5)
Makie.ylims!(axaisnofd, -1.5, 4.75)

Legend(gslr[3, 1:2], [LineElement(color=:darkorange, linewidth=4), LineElement(color=:grey, linewidth=4), LineElement(color=:teal, linewidth=4)], ["Optimistic", "Baseline", "Pessimistic"], ["Emissions Scenario"], orientation=:horizontal, framevisible=false, tellheight=true, titleposition=:top, tellwidth=false)

Label(gslr[1, 1, TopLeft()], "c", fontsize=10, font=:bold, padding = (0, 50, 20, 0), halign=:right)
Label(gslr[2, 1, TopLeft()], "d", fontsize=10, font=:bold, padding = (0, 50, 20, 0), halign=:right)
Label(gslr[1, 2, TopLeft()], "e", fontsize=10, font=:bold, padding = (0, 50, 20, 0), halign=:right)
Label(gslr[2, 2, TopLeft()], "f", fontsize=10, font=:bold, padding = (0, 50, 20, 0), halign=:right)

colgap!(fig.layout, 1, Relative(0.05))

save(joinpath(@__DIR__, "..", "figures", "ensemble_projections.pdf"), fig)
