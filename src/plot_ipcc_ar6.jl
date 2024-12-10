# activate the environment
using Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

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


temps_ipcc = [1.5, 2.0, 3.0, 4.0, 5.0]
quantiles = [0, 5, 50, 95, 100]
function ipcc_quantiles(temp, q, variable)
    ar6_path = joinpath(@__DIR__, "data", "ar6")
    ipcc_nc = "$(ar6_path)/global/confidence_output_files/medium_confidence/tlim$(temp)win0.25/$(variable)_tlim$(temp)win0.25_medium_confidence_values.nc"
    ipcc_out = ncread(ipcc_nc, "sea_level_change")
    ipcc_q = 100 * ncread(ipcc_nc, "quantiles")
    ipcc_yr = ncread(ipcc_nc, "years")
    idx_yr_2100 = findfirst(ipcc_yr .== 2100)
    idx_q = findfirst(ipcc_q .== q)
    return ipcc_out[:, idx_yr_2100, idx_q][1]
end

ais_ipcc_q = [ipcc_quantiles(temp, q, "AIS") for temp in temps_ipcc, q in quantiles] ./ 1000
gis_ipcc_q = [ipcc_quantiles(temp, q, "GIS") for temp in temps_ipcc, q in quantiles] ./ 1000
gsic_ipcc_q = [ipcc_quantiles(temp, q, "glaciers") for temp in temps_ipcc, q in quantiles] ./ 1000
slr_ipcc_q = [ipcc_quantiles(temp, q, "total") for temp in temps_ipcc, q in quantiles] ./ 1000


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

function read_brick_output(variable, ret_yr=2100, norm_yrs=nothing)
    brick_raw = DataFrame(CSVFiles.load(joinpath(@__DIR__, "results", "default", "$(variable).csv")))
    brick_out = norm_df(brick_raw, ret_yr, norm_yrs)
    return brick_out
end


# normalize SLR relative to 1995-2014 and compute quantiles
slr_norm = 1995:2014
temp_norm = 1850:1900

ais_out = read_brick_output("antarctic", 2100, slr_norm)
gis_out = read_brick_output("greenland", 2100, slr_norm)
gsic_out = read_brick_output("gsic", 2100, slr_norm)
slr_out = read_brick_output("gmslr", 2100, slr_norm)
temp_out = read_brick_output("temperature", 2100, temp_norm)

## make plot
fig = Makie.Figure(size=(1200, 1200), fontsize=24)
function plot_comparison!(brick_out, slr_ipcc_q, temps, temps_ipcc, labels, subplot_label, title; f=Figure())
    ax_gmt = Makie.Axis(f[1,1])
    ax_main = Makie.Axis(f[2,1], xlabel=labels[1], ylabel=labels[2], xgridvisible=false, ygridvisible=false)
    ax_brick = Makie.Axis(f[2,2])
    # link axes and hide decorations for marginal plots
    linkyaxes!(ax_main, ax_brick)
    linkxaxes!(ax_main, ax_gmt)
    hidedecorations!(ax_brick)
    hidedecorations!(ax_gmt)
    hidespines!(ax_brick)
    hidespines!(ax_gmt)

    # plot marginal densities
    Makie.density!(ax_gmt, temps)
    Makie.density!(ax_brick, brick_out, direction=:y)

    colsize!(f, 2, Auto(0.25))
    rowsize!(f, 1, Auto(0.25))
    colgap!(f, 1, Relative(-0.01))
    rowgap!(f, 1, Relative(-0.01))

    Label(f[1, 1, Top()], title, font=:bold,fontsize=26,
    halign = :center)    
    Label(f[1, 1, TopLeft()], subplot_label, font=:bold,fontsize=26,
    padding = (0, 5, 5, 5),
    halign = :left)    

    # plot BRICK distributions
    slr_brick_15deg = brick_out[1.45 .<= temps .<= 1.55]
    Makie.violin!(ax_main, 1.5 .+ zeros(length(slr_brick_15deg)), slr_brick_15deg, show_median=true, side=:left, color=:orange)

    slr_brick_2deg = brick_out[1.95 .<= temps .<= 2.05]
    Makie.violin!(ax_main, 2 .+ zeros(length(slr_brick_2deg)), slr_brick_2deg, show_median=true, side=:left, color=:orange)

    slr_brick_3deg = brick_out[2.95 .<= temps .<= 3.05]
    Makie.violin!(ax_main, 3 .+ zeros(length(slr_brick_3deg)), slr_brick_3deg, show_median=true, side=:left, color=:orange)

    slr_brick_4deg = brick_out[3.95 .<= temps .<= 4.05]
    Makie.violin!(ax_main, 4 .+ zeros(length(slr_brick_4deg)), slr_brick_4deg, show_median=true, side=:left, color=:orange)

    slr_brick_5deg = brick_out[4.95 .<= temps .<= 5.05]
    Makie.violin!(ax_main, 5 .+ zeros(length(slr_brick_5deg)), slr_brick_5deg, show_median=true, side=:left, color=:orange)

    Makie.errorbars!(ax_main, temps_ipcc .+ 0.05, slr_ipcc_q[:, 3], slr_ipcc_q[:, 3] - slr_ipcc_q[:, 1], slr_ipcc_q[:, 5] - slr_ipcc_q[:, 3], linewidth=1, whiskerwidth=5, color=:blue)
    Makie.errorbars!(ax_main, temps_ipcc .+ 0.05, slr_ipcc_q[:, 3], slr_ipcc_q[:, 3] - slr_ipcc_q[:, 2], slr_ipcc_q[:, 4] - slr_ipcc_q[:, 3], linewidth=3, whiskerwidth=5, color=:blue)
    Makie.scatter!(ax_main, temps_ipcc .+ 0.05, slr_ipcc_q[:, 3], color=:blue)

    Makie.xlims!(ax_main, 1, 5.25)

    return f
end

f = Figure(size=(1200, 1200), fontsize=24)
ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()
gc = f[2, 1] = GridLayout()
gd = f[2, 2] = GridLayout()

plot_comparison!(slr_out, slr_ipcc_q, temp_out, temps_ipcc, ["Global Mean Temperature in 2100 (°C)", "Global Sea Level Anomaly in 2100 (m)"], "a)", "Global Mean Sea Level"; f=ga)
plot_comparison!(ais_out, ais_ipcc_q, temp_out, temps_ipcc, ["Global Mean Temperature in 2100 (°C)", "SLR Contribution in 2100 (m SLE-eq)"], "b)", "Antarctic Ice Sheet"; f=gb)
plot_comparison!(gis_out, gis_ipcc_q, temp_out, temps_ipcc, ["Global Mean Temperature in 2100 (°C)", "SLR Contribution in 2100 (m SLE-eq)"], "c)", "Greenland Ice Sheet"; f=gc)
plot_comparison!(gsic_out, gsic_ipcc_q, temp_out, temps_ipcc, ["Global Mean Temperature in 2100 (°C)", "SLR Contribution in 2100 (m SLE-eq)"], "d)", "Glaciers"; f=gd)


# make legend
elem_ipcc_90 = LineElement(color = :blue, linestyle = nothing, linewidth=3)
elem_ipcc_full = LineElement(color = :blue, linestyle = nothing, linewidth=1)
elem_ipcc_med = MarkerElement(color=:blue, marker=:circle)
elem_brick = PolyElement(color=:orange)

leg = Legend(f[3, 1:2], [elem_ipcc_med, elem_ipcc_90, elem_ipcc_full, elem_brick], ["AR6 Median", "AR6 90% CI", "AR6 Full Range", "BRICK Projection"], orientation=:horizontal, tellwidth=true)

CairoMakie.save(joinpath(@__DIR__, "..", "figures", "brick_ipcc_compare.png"), f)
