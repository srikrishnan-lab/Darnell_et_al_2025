import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSVs
using DataFrames # data structure for indices
using Makie # plotting library
using CairoMakie
using Measures # adjust margins with explicit measures
using StatsBase # get mean function and density

output_dir = "results/default"
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))
temperature = DataFrame(CSVFiles.load(joinpath(output_dir, "temperature.csv")))
conc_co2 = DataFrame(CSVFiles.load(joinpath("results", "ssp", "ssp119", "concentrations.csv")))
ohc = DataFrame(CSVFiles.load(joinpath(output_dir, "ocean_heat.csv")))
slr = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
ais = DataFrame(CSVFiles.load(joinpath(output_dir, "antarctic.csv")))
gis = DataFrame(CSVFiles.load(joinpath(output_dir, "greenland.csv")))
gsic = DataFrame(CSVFiles.load(joinpath(output_dir, "gsic.csv")))

q_temp = mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), temperature)
q_slr = mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), slr)
q_co2 = mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), conc_co2)
q_ohc = mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), ohc)
q_ais = mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), ais)
q_gis = mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), gis)
q_gsic = mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), gsic)

calibration_data = DataFrame(CSVFiles.load(joinpath(@__DIR__, "data", "calibration", "all_calibration_data_combined.csv")))

hind_yrs = 1850:2020
hind_idx = findall((!isnothing).(indexin(names(temperature), string.(hind_yrs))))

norm_ohc = 1960:1990
norm_temp = 1850:1900
norm_gmslr = 1995:2013
norm_idx_ohc = findall((in)(norm_ohc), calibration_data.year) 
norm_idx_temp = findall((in)(norm_temp), calibration_data.year) 
norm_idx_gmslr = findall((in)(norm_gmslr), calibration_data.year) 
temp_norm = calibration_data.hadcrut_temperature_obs .- mean(calibration_data.hadcrut_temperature_obs[norm_idx_temp])
ohc_norm = calibration_data.ocean_heat_obs .- mean(calibration_data.ocean_heat_obs[norm_idx_ohc])
gmslr_norm = calibration_data.gmsl_obs .- mean(calibration_data.gmsl_obs[norm_idx_gmslr])


# Calculate indices for each year that has an observation in calibration data sets.
indices_maunaloa_co2_data  = findall(x-> !ismissing(x), calibration_data.maunaloa_co2_obs)
indices_lawdome_co2_data   = findall(x-> !ismissing(x), calibration_data.lawdome_co2_obs)
indices_oceanco2_flux_data = findall(x-> !ismissing(x), calibration_data.oceanco2_flux_obs)
indices_temperature_data   = findall(x-> !ismissing(x), calibration_data.hadcrut_temperature_obs)
indices_oceanheat_data     = findall(x-> !ismissing(x), calibration_data.ocean_heat_obs)
indices_glaciers_data      = findall(x-> !ismissing(x), calibration_data.glaciers_obs)
indices_greenland_data     = findall(x-> !ismissing(x), calibration_data.merged_greenland_obs) # Use merged Greenland data.
indices_antarctic_data     = findall(x-> !ismissing(x), calibration_data.antarctic_imbie_obs)
indices_gmsl_data          = findall(x-> !ismissing(x), calibration_data.gmsl_obs)

# Combine CO₂ indices from Law Dome and Mauna Loa observations.
indices_co2_data = sort(vcat(indices_lawdome_co2_data, indices_maunaloa_co2_data))

# Calculate number of ice core observations for CO₂ (used for indexing).
n_lawdome_co2 = length(indices_lawdome_co2_data)

# plot constraints with posterior projections

fig = Figure(size=(700, 500), fontsize=16, figure_padding=10)

g_temp = fig[1, 1]
ax_temp = Axis(g_temp[1, 1], xlabel="Year", ylabel="GMT Anomaly (°C)", alignmode=Inside())
Makie.band!(ax_temp, hind_yrs, Matrix(q_temp)[1, hind_idx], Matrix(q_temp)[3, hind_idx], color=:blue, alpha=0.3, label=false)
proj = Makie.lines!(ax_temp, hind_yrs, Matrix(q_temp)[2, hind_idx], color=:blue, linewidth=2, label="Projection")
obs = Makie.scatter!(ax_temp, (1850:2020)[indices_temperature_data], temp_norm[indices_temperature_data], color=:black, markersize=8)

g_co2 = fig[1, 2]
ax_co2 = Axis(g_co2[1, 1], xlabel="Year", ylabel="CO₂ Concentration (ppm)", alignmode=Inside())
Makie.band!(ax_co2, hind_yrs, Matrix(q_co2)[1, hind_idx], Matrix(q_co2)[3, hind_idx], color=:blue, alpha=0.3, label=false)
Makie.lines!(ax_co2, hind_yrs, Matrix(q_co2)[2, hind_idx], color=:blue, linewidth=2, label="Projection")
Makie.scatter!(ax_co2, hind_yrs[indices_lawdome_co2_data], convert(Vector{Float64}, calibration_data.lawdome_co2_obs[indices_lawdome_co2_data]), color=:black, markersize=8, marker=:utriangle)
Makie.scatter!(ax_co2, hind_yrs[indices_maunaloa_co2_data], convert(Vector{Float64}, calibration_data.maunaloa_co2_obs[indices_maunaloa_co2_data]), color=:black, markersize=8)

g_ohc = fig[2, 1]
ohc_yrs = 1950:2000
ohc_idx = findall((!isnothing).(indexin(names(ohc), string.(ohc_yrs))))
ax_ohc = Axis(g_ohc[1, 1], xlabel="Year", ylabel="OHC (10²² J)", alignmode=Inside())
Makie.band!(ax_ohc, ohc_yrs, Matrix(q_ohc)[1, ohc_idx], Matrix(q_ohc)[3, ohc_idx], color=:blue, alpha=0.3, label=false)
Makie.lines!(ax_ohc, ohc_yrs, Matrix(q_ohc)[2, ohc_idx], color=:blue, linewidth=2, label="Projection")
Makie.scatter!(ax_ohc, hind_yrs[indices_oceanheat_data], convert(Vector{Float64},ohc_norm[indices_oceanheat_data]), color=:black, markersize=8)

g_gmsl = fig[2, 2]
ax_gmsl = Axis(g_gmsl[1, 1], xlabel="Year", ylabel="GMSL Anomaly (m)", alignmode=Inside())
Makie.band!(ax_gmsl, hind_yrs, Matrix(q_slr)[1, hind_idx], Matrix(q_slr)[3, hind_idx], color=:blue, alpha=0.3, label=false)
Makie.lines!(ax_gmsl, hind_yrs, Matrix(q_slr)[2, hind_idx], color=:blue, linewidth=2, label="Projection")
Makie.scatter!(ax_gmsl, (1850:2020)[indices_gmsl_data], convert(Vector{Float64}, gmslr_norm[indices_gmsl_data]), color=:black, markersize=8)

Makie.save("figures/observational_constraints.png", fig)