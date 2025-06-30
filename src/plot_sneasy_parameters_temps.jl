import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSVs
using DataFrames # data structure for indices
using Makie # plotting library
using CairoMakie
using Measures # adjust margins with explicit measures
using StatsBase # get mean function and density
using Interpolations

data_path = joinpath(@__DIR__, "..",  "data")

## load emissions scenarios from SSP database
emis_dat = DataFrame(CSVFiles.load(joinpath(data_path, "rcmip-emissions-annual-means-v5-1-0.csv")))
conc_dat = DataFrame(CSVFiles.load(joinpath(data_path, "rcmip-concentrations-annual-means-v5-1-0.csv")))
forc_dat = DataFrame(CSVFiles.load(joinpath(data_path, "rcmip-radiative-forcing-annual-means-v5-1-0.csv")))
filter!(:Mip_Era => isequal("CMIP6"), emis_dat)
filter!(:Mip_Era => isequal("CMIP6"), conc_dat)
filter!(:Mip_Era => isequal("CMIP6"), forc_dat)
scenario_regexp = r"^ssp[0-9]*$"
filter!(:Scenario => contains(scenario_regexp), emis_dat)
filter!(:Scenario => contains(scenario_regexp), conc_dat)
filter!(:Scenario => contains(scenario_regexp), forc_dat)
filter!(:Region => isequal("World"), emis_dat)
filter!(:Region => isequal("World"), conc_dat)
filter!(:Region => isequal("World"), forc_dat)
filter!(:Variable => isequal("Emissions|CO2"), emis_dat)
aerosol_forc = filter(:Variable => isequal("Effective Radiative Forcing|Anthropogenic|Aerosols"), forc_dat)

output_dir = joinpath(@__DIR__, ".."  "results")
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "default", "parameters.csv")))
temperature = DataFrame(CSVFiles.load(joinpath(output_dir, "default","temperature.csv")))
conc_co2 = DataFrame(CSVFiles.load(joinpath(output_dir, "default", "concentrations.csv")))
ohc = DataFrame(CSVFiles.load(joinpath(output_dir, "default", "ocean_heat.csv")))


idx2020 = findfirst(names(co2_forc) .== "2020")
idx2100 = findfirst(names(temperature) .== "2100")
ssp_idx2014 = findfirst(names(co2_forc) .== "2014")
ssp_idx1750 = findfirst(names(co2_forc) .== "1750")
ssp_idx1850 = findfirst(names(co2_forc) .== "1850")
idx2014 = findfirst(names(temperature) .== "2014")
idx1850 = findfirst(names(temperature) .== "1850")
idx1900 = findfirst(names(temperature) .== "1900")

# interpolate co2 emissions on an annual grid using a linear interpolation
xs_idx = findall((!ismissing).(collect(emis_dat[1, ssp_idx1750:ssp_idx2014])))
yrs = parse.(Int, names(emis_dat)[ssp_idx1750:ssp_idx2014])
xs = yrs[xs_idx]
interp = LinearInterpolation(xs, collect(emis_dat[1, xs_idx .+ ssp_idx1750 .- 1]))
emis_hist = interp.(yrs)


fig = Figure(size=(700, 600), fontsize=16, figure_padding=10)
# panel 1: SNEASY ECS distribution
ga = fig[1, 1]
ecs = parameters[:, :climate_sensitivity]
ecs_assessed = [2.3, 2.6, 3.1, 3.9, 4.7]
ecs_q = quantile(ecs, [0.05, 0.17, 0.5, 0.83, 0.95])
ax1 = Axis(ga[1, 1], xlabel="Equilibrium Climate Sensitivity (°C)", ylabel="Density", alignmode=Inside(), xticks=0:1:8, yticksvisible=false, yticklabelsvisible=false)
Makie.hist!(ax1, ecs, normalization=:pdf, strokecolor=:black, strokewidth=1, color=:blue)
Makie.rangebars!(ax1, [-0.1], [ecs_q[2]], [ecs_q[4]], direction=:x, linewidth=25, color=:blue)
Makie.rangebars!(ax1, [-0.1], [ecs_q[1]], [ecs_q[5]], direction=:x, linewidth=5, whiskerwidth=10, color=:blue)
Makie.rangebars!(ax1, [ecs_q[3]], [-0.01], [-0.15], color=:white)
Makie.rangebars!(ax1, [-0.25], [ecs_assessed[2]], [ecs_assessed[4]], direction=:x, linewidth=25, color=:black)
Makie.rangebars!(ax1, [-0.25], [ecs_assessed[1]], [ecs_assessed[5]], direction=:x, linewidth=5, whiskerwidth=10, color=:black)
Makie.rangebars!(ax1, [ecs_assessed[3]], [-0.2], [-0.3], color=:white)


# panel 2: SNEASY TCRE
gb = fig[1, 2]
tcre_assessed = [1.03, 1.4, 1.77, 2.14, 2.51]
temp_norm = temperature[:, idx2014] - mean.(eachrow(temperature[:, idx1850:idx1900]))
tcre = tcre = temp_norm ./ (sum(emis_hist) / 1e6 / 3.67) # need to convert to TtC instead of MtCO2
tcre_q = quantile(tcre, [0.05, 0.17, 0.5, 0.83, 0.95])
ax2 = Axis(gb[1, 1], xlabel="TCRE (°C/TtC)", ylabel="Density", alignmode=Inside(), xticks=0:0.25:2.75, yticksvisible=false, yticklabelsvisible=false)
Makie.hist!(ax2, tcre, normalization=:pdf, strokecolor=:black, strokewidth=1, color=:blue)
Makie.rangebars!(ax2, [-0.27], [tcre_q[2]], [tcre_q[4]], direction=:x, linewidth=25, color=:blue)
Makie.rangebars!(ax2, [-0.27], [tcre_q[1]], [tcre_q[5]], direction=:x, linewidth=5, whiskerwidth=10, color=:blue)
Makie.rangebars!(ax2, [tcre_q[3]], [-0.01], [-0.5], color=:white)
Makie.rangebars!(ax2, [-0.75], [tcre_assessed[2]], [tcre_assessed[4]], direction=:x, linewidth=25, color=:black)
Makie.rangebars!(ax2, [-0.75], [tcre_assessed[1]], [tcre_assessed[5]], direction=:x, linewidth=5, whiskerwidth=10, color=:black)
Makie.rangebars!(ax2, [tcre_assessed[3]], [-0.35], [-0.95], color=:white)


# panel 3: 2014 CO2 effective radiative forcing
gc = fig[2, 1]
co2_assessed = [1.69, 1.80, 1.91]
co2_diff = 5.35 * log.(conc_co2[:, idx2014] ./ conc_co2[:, idx1850])
co2_q = quantile(co2_diff, [0.17, 0.5, 0.83])
ax3 = Axis(gc[1, 1], xlabel="Effective CO₂ RF, 2014 (W/m²)", ylabel="Density", alignmode=Inside(), xticks=1.25:0.25:2.5, yticksvisible=false, yticklabelsvisible=false)
Makie.hist!(ax3, co2_diff, normalization=:pdf, strokecolor=:black, strokewidth=1, color=:blue)
Makie.rangebars!(ax3, [-0.6], [co2_q[1]], [co2_q[3]], direction=:x, linewidth=25, color=:blue)
Makie.rangebars!(ax3, [co2_q[2]], [-0.10], [-1.2], color=:white)
Makie.rangebars!(ax3, [-1.8], [co2_assessed[1]], [co2_assessed[3]], direction=:x, linewidth=25, color=:black)
Makie.rangebars!(ax3, [co2_assessed[2]], [-1.35], [-2.35], color=:white)

# panel 4: 2018 Ocean Heat Content Relative to 1971
gd = fig[2, 2]
α = parameters[:, :rf_scale_aerosol]
aer_rf_assessed = [-1.37, -1.01, -0.63]
aer_rf_q = quantile(aerosol_forc[1, ssp_idx2014] .* α, [0.17, 0.5, 0.83])
ax4 = Axis(gd[1, 1], xlabel="Effective Aerosol RF, 2014 (W/m²)", ylabel="Density", alignmode=Inside(), xticks=-2:0.5:-0.5, yticksvisible=false, yticklabelsvisible=false)
Makie.hist!(ax4, aerosol_forc[1, :"2014"] * α, normalization=:pdf, strokecolor=:black, strokewidth=1, color=:blue)
Makie.rangebars!(ax4, [-0.22], [aer_rf_q[1]], [aer_rf_q[3]], direction=:x, linewidth=25, color=:blue)
Makie.rangebars!(ax4, [aer_rf_q[2]], [-0.35], [-0.020], color=:white)
Makie.rangebars!(ax4, [-0.65], [aer_rf_assessed[1]], [aer_rf_assessed[3]], direction=:x, linewidth=25, color=:black)
Makie.rangebars!(ax4, [aer_rf_assessed[2]], [-0.5], [-0.85], color=:white)

Label(fig[1, 1, TopLeft()], "a", fontsize=20, font=:bold, padding = (0, 50, 10, 0), halign=:right)
Label(fig[1, 2, TopLeft()], "b", fontsize=20, font=:bold, padding = (0, 50, 10, 0), halign=:right)
Label(fig[2, 1, TopLeft()], "c", fontsize=20, font=:bold, padding = (0, 50, 10, 0), halign=:right)
Label(fig[2, 2, TopLeft()], "d", fontsize=20, font=:bold, padding = (0, 50, 10, 0), halign=:right)

CairoMakie.save(joinpath(@__DIR__, "..", "figures", "rcmip_proxy.png"), fig)