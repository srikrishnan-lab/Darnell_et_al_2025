import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSV of Shapley indices
using DataFrames # data structure for indices
using Makie # plotting library
using CairoMakie
using ColorSchemes
using Measures # adjust margins with explicit measures
using Statistics # get mean function

output_dir = "results/peaking"
slr_out = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
ais_out = DataFrame(CSVFiles.load(joinpath(output_dir, "antarctic.csv")))
gis_out = DataFrame(CSVFiles.load(joinpath(output_dir, "greenland.csv")))
temp_out = DataFrame(CSVFiles.load(joinpath(output_dir, "temperature.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))

# normalize relative to 2000
idx_2000 = findfirst(names(slr_out) .== "2000")
idx_2100 = findfirst(names(slr_out) .== "2100")
idx_1850 = findfirst(names(slr_out) .== "1850")
idx_1900 = findfirst(names(slr_out) .== "1900")
for row in axes(slr_out,1 )
    foreach(col -> slr_out[row, col] -= slr_out[row, idx_2000], axes(slr_out, 2))
    foreach(col -> ais_out[row, col] -= ais_out[row, idx_2000], axes(slr_out, 2))
    foreach(col -> gis_out[row, col] -= gis_out[row, idx_2000], axes(slr_out, 2))
    foreach(col -> temp_out[row, col] -= mean(temp_out[row, idx_1850:idx_1900]), axes(slr_out, 2))
end

function plot_peaking_ensemble(parameters, emissions, slr_out, temp_out, ais_out)

    function plot_quantiles(df, params, peak_yrs, idx, label, colors; f=Figure(), limits=(nothing, nothing))
        val_low_all = df[idx, :]
        param_low_all = params[idx, :]
        ax = Axis(f, xlabel="Year", ylabel=label)
        for (i, yr) in pairs(peak_yrs)
            q = mapcols(x -> quantile(x, [0.025, 0.5, 0.975]), val_low_all[param_low_all[:, :t_peak] .== yr, :])
            Makie.lines!(ax, parse.(Int64, names(df)), Vector(q[2, :]), color=colors[i], linewidth=2)
            Makie.lines!(ax, parse.(Int64, names(df)), Vector(q[1, :]), color=colors[i], linestyle=:dash)
            Makie.lines!(ax, parse.(Int64, names(df)), Vector(q[3, :]), color=colors[i],  linestyle=:dash)
        end
        Makie.ylims!(ax, limits[1], limits[2])
        Makie.xlims!(ax, 2020, 2100)
        return f
    end

    ais_threshold = 15.42 .+ 0.8365 * parameters[:, :antarctic_temp_threshold] .- map(mean, eachrow(temp_out[:, idx_1850:idx_1900]))
    ais_exceed_yr = Vector{Union{Float64, Nothing}}(undef, length(ais_threshold))
    temp_exceed_max = Vector{Union{Float64, Nothing}}(undef, length(ais_threshold))

    for j in eachindex(ais_threshold)
        exceed_ais = [temp_out[j, k] > ais_threshold[j] for k in axes(temp_out, 2)]
        exceed_ais_idx = findfirst(exceed_ais)
        if !isnothing(exceed_ais_idx)
            ais_exceed_yr[j] = parse(Int64, names(temp_out)[exceed_ais_idx])
            temp_exceed_max[j] = temp_out[j, exceed_ais_idx]
        end
    end

    ais_exceed = hcat(parameters[:, :t_peak], ais_exceed_yr)
    ais_exceed_idx = .!(isnothing).(ais_exceed[:, 2])
    ais_exceed = Float64.(ais_exceed[ais_exceed_idx, :])
    temp_all = hcat(parameters[:, :t_peak], temp_exceed_max)
    temp_exceed = Float64.(temp_all[ais_exceed_idx, :])

    pk_colors = ColorSchemes.seaborn_colorblind[1:6]

    fig = Figure(size=(800, 600), fontsize=14, figure_padding=(10, 20, 10, 10))

    gts = fig[1, 1:3] = GridLayout()
    ghist = fig[2, 1:3] = GridLayout()

    plot_quantiles(emissions, parameters, unique(parameters[!, :t_peak]), 1:length(ais_threshold), "CO₂ Emissions\n(GtCO2/yr)", pk_colors, f=gts[1, 1])
    plot_quantiles(slr_out, parameters, unique(parameters[!, :t_peak]), 1:length(ais_threshold), "Global Mean\nSea Level (m)", pk_colors, f=gts[1, 2], limits=(0, 1.5))
    #plot_quantiles(gis_out, parameters, unique(parameters[!, :t_peak]), low_all_idx, "GIS Contribution (m/yr)", pk_colors, f=fig[2, 1], limits=(0, 0.5))
    plot_quantiles(ais_out, parameters, unique(parameters[!, :t_peak]), 1:length(ais_threshold), "AIS Contribution (m)", pk_colors, f=gts[1, 3], limits=(-0.05, 0.75))

    # make legend
    leg_elem = Vector{LineElement}(undef, length(unique(parameters[!, :t_peak])))
    leg_label = Vector{String}(undef, length(unique(parameters[!, :t_peak])))
    for i in 1:length(unique(parameters[!, :t_peak]))
        leg_elem[i]= LineElement(color = pk_colors[i], linestyle = nothing, linewidth=4)
        leg_label[i] = string(Int64(unique(parameters[!, :t_peak])[i]))
    end

    ax_ais_hist = Axis(ghist[1, 1], xlabel="Year AIS Threshold Exceeded", ylabel="Emissions Peak Year", yticks=(1:6, string.(2050:10:2100)))
    rainclouds!(ax_ais_hist, string.(Int64.(ais_exceed[:, 1])), ais_exceed[:, 2], bins=25, color=pk_colors[indexin(ais_exceed[:, 1], unique(ais_exceed[:, 1]))], orientation=:horizontal, violin_limits=extrema, gap=0.1)
    Makie.xlims!(ax_ais_hist, 2030, 2300)
    Makie.yticks!

    ax_temp_hist = Axis(ghist[1, 2], xlabel="AIS Trigger Temperature (°C)", ylabel="Emissions Peak Year", yticks=(1:6, string.(2050:10:2100)))
    rainclouds!(ax_temp_hist, string.(Int64.(temp_exceed[:, 1])), temp_exceed[:, 2], bins=25, color=pk_colors[indexin(temp_exceed[:, 1], unique(temp_exceed[:, 1]))], orientation=:horizontal, violin_limits=extrema, gap=0.1)


    # ax_temp2_hist = Axis(fig[3, 1], xlabel="Maximum Temperature Reached (°C)", ylabel="Emissions Peak Year", yticks=(1:6, string.(2050:10:2100)))
    # rainclouds!(ax_temp2_hist, string.(Int64.(temp_no_exceed[:, 1])), temp_no_exceed[:, 2], bins=25, color=pk_colors[indexin(temp_no_exceed[:, 1], unique(temp_no_exceed[:, 1]))], orientation=:horizontal, violin_limits=extrema, gap=0.1)
    # Makie.xlims!(ax_temp2_hist, 1.3, 3.5)

    # add numbers

    text!(ax_ais_hist, [(2200, i + 0.1) for i in 1:6], text = "n = " .* string.([sum(temp_exceed[:, 1] .== yr) for yr in 2050:10:2100]) .* " (" .* string.([round(Int64, 100 * sum(temp_exceed[:, 1] .== yr) / sum(temp_all[:, 1] .== yr)) for yr in 2050:10:2100]) .* "%)")

    Label(gts[1, 1, TopLeft()], "a", fontsize=14, font=:bold, padding = (0, 50, 0, 0), halign=:right)
    Label(gts[1, 2, TopLeft()], "b", fontsize=14, font=:bold, padding = (0, 50, 0, 0), halign=:right)
    Label(gts[1, 3, TopLeft()], "c", fontsize=14, font=:bold, padding = (0, 50, 0, 0), halign=:right)
    Label(ghist[1, 1, TopLeft()], "d", fontsize=14, font=:bold, padding = (0, 50, 0, 0), halign=:right)
    Label(ghist[1, 2, TopLeft()], "e", fontsize=14, font=:bold, padding = (0, 50, 0, 0), halign=:right)

    Legend(fig[3, 1:3], leg_elem, leg_label, "Year Emissions Peak", orientation=:horizontal, nbanks=1)

    rowsize!(fig.layout, 2, Auto(1.5))

    return fig
end

# make main figure, top 10%
f_base = plot_peaking_ensemble(parameters, emissions, slr_out, temp_out, ais_out)
save("figures/peaking_threshold_all_supp.png", f_base)
