#########################################################################
# plot_slr_temperatures.jl                                              #
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

using CSVFiles # read CSVs
using DataFrames # data structure for indices
using Makie # plotting library
using CairoMakie
using Measures # adjust margins with explicit measures
using StatsBase # get mean function and density

# load ensemble
output_dir = "results/default"
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
temperature = DataFrame(CSVFiles.load(joinpath(output_dir, "temperature.csv")))
gmslr = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))

# define function to normalize data relative to some normalization period mean
function normalize_data!(dat, norm_yrs=nothing)
    # normalize to relevant period  defined by norm_yrs
    if !isnothing(norm_yrs)
        idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
        for row in axes(dat, 1)
            foreach(col -> dat[row, col] -= mean(dat[row, idx_norm]), axes(dat, 2))
        end
    end
    return dat
end

normalize_data!(temperature, 1850:1900)
normalize_data!(gmslr, 1995:2014)

# define function to compute quantiles relative to some normalization period
function compute_norm_quantiles(dat, norm_yrs=nothing)
    # normalize to relevant period  defined by norm_yrs
    if !isnothing(norm_yrs)
        idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
        for row in axes(dat, 1)
            foreach(col -> dat[row, col] -= mean(dat[row, idx_norm]), axes(dat, 2))
        end
    end
    # compute median and 95% prediction interval
    quantiles = mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), dat)
    return quantiles
end

emissions_q = compute_norm_quantiles(emissions)

idx1850 = findfirst(names(temperature) .== "1850")
idx1900 = findfirst(names(temperature) .== "1900")

# convert AIS threshold from local to global mean temperature
# uses the regression fit from the model calibration
ais_threshold = [15.42 .+ 0.8365 * parameters[i, :antarctic_temp_threshold] - mean(temperature[i, idx1850:idx1900]) for i in axes(temperature, 1)]

# find years in which AIS threshold is exceeded
ais_exceed_yr = Vector{Union{Float64, Nothing}}(undef, length(ais_threshold))
for j in eachindex(ais_threshold)
    exceed_ais = [temperature[j, k] > ais_threshold[j] for k in axes(temperature, 2)]
    exceed_ais_idx = findfirst(exceed_ais)
    if isnothing(exceed_ais_idx)
        ais_exceed_yr[j] = nothing
    else
        ais_exceed_yr[j] = parse(Int64, names(temperature)[exceed_ais_idx])
    end
end

ais_exceed = hcat(parameters[:, :t_peak], ais_exceed_yr)
ais_exceed[:, 2] = replace(ais_exceed[:, 2], nothing => 2305)

# find cumulative emissions from 2022--2100
idx2100 = findfirst(names(gmslr) .== "2100")
idx2022 = findfirst(names(gmslr) .== "2022")
cum_emissions = [sum(emissions[i, idx2022:idx2100]) for i in 1:nrow(gmslr)]


idx_paris = findall(temperature[!, idx2100] .< 2.0) # find Paris consistent SOWs

fig = Figure(size=(400, 600), fontsize=14, figure_padding=10)

# set up layout
ga = fig[1, 1] = GridLayout()
gb = fig[2, 1] = GridLayout()

ax_slr_emissions = Axis(ga[1, 1], xlabel="Cumulative CO₂ Emissions (Gt CO₂)", ylabel="Global Sea Level Anomaly (m)", alignmode=Inside())

colors_paris = cgrad(:vik100, [0.1, 0.5, 0.9], rev=true)
plt_slr_emissions =  Makie.scatter!(ax_slr_emissions, cum_emissions, gmslr[:, idx2100], color=ais_threshold, colormap=colors_paris, markersize=6)
cbscatter = Colorbar(ga[1, 2], plt_slr_emissions, label="AIS Temperature Threshold (°C)", vertical=true, flipaxis=false, tellwidth=true, tellheight=false)

Label(ga[1, 1, TopLeft()], "a", fontsize=18, font=:bold, padding = (0, 50, 0, 0), halign=:right)

ax_slr_paris = Axis(gb[1, 1], xlabel="Cumulative CO₂ Emissions (Gt CO₂)", ylabel="Global Sea Level Anomaly (m)", alignmode=Inside())

colors_paris = cgrad(:vik100, [0.1, 0.5, 0.9], rev=true)
plt_slr_paris =  Makie.scatter!(ax_slr_paris, cum_emissions[idx_paris], gmslr[idx_paris, idx2100], color=ais_threshold[idx_paris], colormap=colors_paris, markersize=6)
cbscatter = Colorbar(gb[1, 2], plt_slr_paris, label="AIS Temperature Threshold (°C)", vertical=true, flipaxis=false, tellwidth=true, tellheight=false)


Label(gb[1, 1, TopLeft()], "b", fontsize=18, font=:bold, padding = (0, 50, 0, 0), halign=:right)

CairoMakie.save("figures/slr_temps_all.png", fig)








