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
using Loess

# load ensemble
output_dir = "results/default"
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
temperature = DataFrame(CSVFiles.load(joinpath(output_dir, "temperature.csv")))
gmslr = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
ais = DataFrame(CSVFiles.load(joinpath(output_dir, "antarctic.csv")))

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

normalize_data!(temperature, 1850:1900)
normalize_data!(gmslr, 1995:2014)
normalize_data!(ais, 1995:2014)

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
idx2000 = findfirst(names(gmslr) .== "2000")

cum_emissions = [sum(emissions[i, idx2022:idx2100]) for i in 1:nrow(gmslr)]
avg_temp = temperature[:, idx2100] - temperature[:, idx2000]  # find emissions average
avg_emissions = cum_emissions / (2100 - 2000 + 1) # find emissions average
avg_gmslr = gmslr[:, idx2100] - gmslr[:, idx2000]
avg_ais = ais[:, idx2100] - ais[:, idx2000]

# compute densities
dens_gmslr = kde(hcat(avg_temp, avg_gmslr))
dens_ais = kde(hcat(avg_temp, avg_ais))
# regress temperature against the outputs
# bootstrap confidence intervals for the loess
function loess_bootstrap(xs, ys; nboot=1_000, xstep=0.01, ci_width=0.95)
    
    mod = loess(xs, ys)
    us = range(extrema(xs)...; step=xstep)
    vs_out = zeros(3, length(us))
    vs_out[2, :] .= predict(mod, us)

    vs_boot = zeros(nboot, length(us)) # will hold all bootstrap predictions here
    nx = length(xs)
    # bootstrap
    for iboot=1:nboot
        # sample
        idx_boot = sample( (1:nx) |> Array |> vec, nx-2, replace=true)
        push!(idx_boot, argmin(xs))
        push!(idx_boot, argmax(xs))
        # run model on bootstrap sample
        model = loess(xs[idx_boot], ys[idx_boot], span=0.5)

        vs_boot[iboot, :] .= Loess.predict(model, us)
    end
    # get confidence intervals
    for i = 1:length(us)
        vs_out[[1, 3], i] .= quantile(vs_boot[:, i], ((1 - ci_width) / 2, (1 + ci_width) / 2))
    end
    return (us, vs_out)
end

temp_grid, slr_loess = loess_bootstrap(avg_temp, avg_gmslr;nboot=200)
temp_grid, ais_loess = loess_bootstrap(avg_temp, avg_ais)

fig = Figure(size=(400, 600), fontsize=14, figure_padding=10)
# set up layout
ga = fig[1, 1] = GridLayout()
gb = fig[2, 1] = GridLayout()

ax_gmt = Makie.Axis(ga[1,1])
ax_main = Axis(ga[2, 1], xlabel="GMST Anomaly from 2000--2100 (°C)", ylabel="GSLR Rate (m)", alignmode=Inside())
ax_gmslr = Makie.Axis(ga[2,2])

# link axes and hide decorations for marginal plots
linkyaxes!(ax_main, ax_gmslr)
linkxaxes!(ax_main, ax_gmt)
hidedecorations!(ax_gmslr)
hidedecorations!(ax_gmt)
hidespines!(ax_gmslr)
hidespines!(ax_gmt)
# plot marginal densities
Makie.density!(ax_gmt, avg_temp)
Makie.density!(ax_gmslr, avg_gmslr, direction=:y)

colsize!(ga, 2, Auto(0.25))
rowsize!(ga, 1, Auto(0.25))
colgap!(ga, 1, Relative(-0.01))
rowgap!(ga, 1, Relative(-0.01))

Makie.contour!(ax_main, dens_gmslr.x, dens_gmslr.y, dens_gmslr.density, levels=0.05:0.2:2.15)
Makie.band!(ax_main, temp_grid, slr_loess[1, :], slr_loess[3, :], alpha=0.5)
Makie.lines!(ax_main, temp_grid, slr_loess[2, :], linewidth=1)

Makie.xlims!(ax_main, -0.04, 4)
Makie.ylims!(ax_gmslr, 0, 2)

Label(ga[1, 1, TopLeft()], "a", fontsize=18, font=:bold, padding = (0, 50, 0, 0), halign=:right)

ax_gmt = Makie.Axis(gb[1,1])
ax_main = Axis(gb[2, 1], xlabel="GMST Anomaly from 2000--2100 (°C)", ylabel="AIS Contribution (m SLE-eq)", alignmode=Inside())
ax_ais = Makie.Axis(gb[2,2])

# link axes and hide decorations for marginal plots
linkyaxes!(ax_main, ax_ais)
linkxaxes!(ax_main, ax_gmt)
hidedecorations!(ax_ais)
hidedecorations!(ax_gmt)
hidespines!(ax_ais)
hidespines!(ax_gmt)
# plot marginal densities
Makie.density!(ax_gmt, avg_temp)
Makie.density!(ax_ais, avg_ais, direction=:y)

colsize!(gb, 2, Auto(0.25))
rowsize!(gb, 1, Auto(0.25))
colgap!(gb, 1, Relative(-0.01))
rowgap!(gb, 1, Relative(-0.01))

Makie.contour!(ax_main, dens_ais.x, dens_ais.y, dens_ais.density, levels=0.05:0.5:4.35)
Makie.band!(ax_main, temp_grid, ais_loess[1, :], ais_loess[3, :], alpha=0.5)
Makie.lines!(ax_main, temp_grid, ais_loess[2, :], linewidth=1)


Makie.xlims!(ax_main, -0.04, 4)
Makie.ylims!(ax_ais, -0.1, 1)


Label(gb[1, 1, TopLeft()], "b", fontsize=18, font=:bold, padding = (0, 50, 0, 0), halign=:right)

CairoMakie.save("figures/slr_temps_all.png", fig)
