using NetCDF
using CSVFiles
using DataFrames
using Statistics
using KernelDensity
using Interpolations
using Makie
using CairoMakie
using StatsBase
using StatsPlots

# activate the environment
using Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()


temps_ipcc = [1.5, 2.0, 3.0, 4.0, 5.0]
quantiles = [0, 5, 50, 95, 100]
function ipcc_quantiles(temp, q)
    ar6_path = "$(dirname(@__DIR__))/data/ar6"
    ais_nc = "$(ar6_path)/global/confidence_output_files/medium_confidence/tlim$(temp)win0.25/AIS_tlim$(temp)win0.25_medium_confidence_values.nc"
    ais_slr = ncread(ais_nc, "sea_level_change")
    ais_q = 100 * ncread(ais_nc, "quantiles")
    ais_yr = ncread(ais_nc, "years")
    idx_yr_2100 = findfirst(ais_yr .== 2100)
    idx_q = findfirst(ais_q .== q)
    return ais_slr[:, idx_yr_2100, idx_q][1]
end

ipcc_q = [ipcc_quantiles(temp, q) for temp in temps_ipcc, q in quantiles] . / 1000

# get runs from our ensemble
temp_default = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "default", "temperature.csv")))
ais_default = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "default", "antarctic.csv")))

# normalize SLR relative to 1995-2014 and compute quantiles
years = parse.(Int64, names(ais_default))
slr_norm = 1995:2014
temp_norm = 1850:1900
# function to find the 66% prediction interval (optional: after normalizing relative to a year or mean over some normalization period)
function norm_df(dat, ret_yr, norm_yrs)
    # normalize to relevant period  defined by norm_yrs
    idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
    mean_val = map(row -> mean(row[idx_norm]), eachrow(dat))
    for row in axes(dat, 1)
        foreach(col -> dat[row, col] -= mean_val[row], axes(dat, 2))
    end
    idx_ret = findfirst((names(dat) .== string(ret_yr)))
    return dat[:, idx_ret]
end

ais_out = norm_df(ais_default, 2100, slr_norm)
temp_out = norm_df(temp_default, 2100, temp_norm)

fig = Makie.Figure(size=(600, 600))
ax_gmt = Makie.Axis(fig[1,1])
ax_main = Makie.Axis(fig[2,1], xlabel="Global Mean Temperature in 2100 (°C)", ylabel="AIS Contribution to Global Sea Level in 2100 (m)", xgridvisible=false, ygridvisible=false)
ax_ais = Makie.Axis(fig[2,2])
# link axes and hide decorations for marginal plots
linkyaxes!(ax_main, ax_ais)
linkxaxes!(ax_main, ax_gmt)
hidedecorations!(ax_ais)
hidedecorations!(ax_gmt)
hidespines!(ax_ais)
hidespines!(ax_gmt)

# plot marginal densities
Makie.density!(ax_gmt, temp_out)
Makie.density!(ax_ais, ais_out, direction=:y)

colsize!(fig.layout, 2, Auto(0.25))
rowsize!(fig.layout, 1, Auto(0.25))
colgap!(fig.layout, 1, Relative(-0.01))
rowgap!(fig.layout, 1, Relative(-0.01))

brick_15deg = ais_out[1.45 .<= temp_out .<= 1.55]
Makie.violin!(ax_main, 1.5 .+ zeros(length(brick_15deg)), brick_15deg, show_median=true, side=:left, color=:orange)

brick_2deg = ais_out[1.95 .<= temp_out .<= 2.05]
Makie.violin!(ax_main, 2 .+ zeros(length(brick_2deg)), brick_2deg, show_median=true, side=:left, color=:orange)

brick_3deg = ais_out[2.95 .<= temp_out .<= 3.05]
Makie.violin!(ax_main, 3 .+ zeros(length(brick_3deg)), brick_3deg, show_median=true, side=:left, color=:orange)

brick_4deg = ais_out[3.95 .<= temp_out .<= 4.05]
Makie.violin!(ax_main, 4 .+ zeros(length(brick_4deg)), brick_4deg, show_median=true, side=:left, color=:orange)

brick_5deg = ais_out[4.95 .<= temp_out .<= 5.05]
Makie.violin!(ax_main, 5 .+ zeros(length(brick_5deg)), brick_5deg, show_median=true, side=:left, color=:orange)


Makie.errorbars!(ax_main, temps_ipcc .+ 0.05, ipcc_q[:, 3], ipcc_q[:, 3] - ipcc_q[:, 1], ipcc_q[:, 5] - ipcc_q[:, 3], linewidth=1, whiskerwidth=5, color=:blue)
Makie.errorbars!(ax_main, temps_ipcc .+ 0.05, ipcc_q[:, 3], ipcc_q[:, 3] - ipcc_q[:, 2], ipcc_q[:, 4] - ipcc_q[:, 3], linewidth=3, whiskerwidth=5, color=:blue)
Makie.scatter!(ax_main, temps_ipcc .+ 0.05, ipcc_q[:, 3], color=:blue)

Makie.xlims!(ax_main, 1, 5.25)

# make legend
elem_ipcc_90 = LineElement(color = :blue, linestyle = nothing, linewidth=3)
elem_ipcc_full = LineElement(color = :blue, linestyle = nothing, linewidth=1)
elem_ipcc_med = MarkerElement(color=:blue, marker=:circle)
elem_brick = PolyElement(color=:orange)

leg = Legend(fig[3, 1:2], [elem_ipcc_med, elem_ipcc_90, elem_ipcc_full, elem_brick], ["AR6 Median", "AR6 90% CI", "AR6 Full Range", "BRICK Projection"], orientation=:horizontal, tellwidth=true)

CairoMakie.save(joinpath(@__DIR__, "..", "figures", "brick_ipcc_ais_compare.png"), fig)
