import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles
using DataFrames
using Statistics
using StatsBase
using Plots
using StatsPlots

conc_dat = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "data", "rcmip-concentrations-annual-means-v5-1-0.csv")))
filter!(:Mip_Era => isequal("CMIP6"), conc_dat)
scenario_regexp = r"^ssp[0-9]*$"
filter!(:Scenario => contains(scenario_regexp), conc_dat)
filter!(:Region => isequal("World"), conc_dat)
co2_conc = filter(:Variable => isequal("Atmospheric Concentrations|CO2"), conc_dat)

ssps = unique(co2_conc[!, :Scenario])
nsamples = 10_000 # number of samples per simulation
ssp_sim_conc = zeros(nsamples, length(ssps))
for (i, ssp) in pairs(ssps)
    ssp_out = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "ssp", ssp, "concentrations_nonoise.csv")))
    ssp_sim_conc[:, i] = ssp_out[!, :"2100"]
    Makie.boxplot!(ax, i, ssp_sim_conc[:, i])

end
ssp_sim_df = DataFrame(ssp_sim_conc, ssps)
ssp_sim_stack = stack(ssp_sim_df)

# make plot
p = StatsPlots.violin(ssp_sim_stack.variable, ssp_sim_stack.value, show_median=true, label="SNEASY Simulation", ylabel="CO₂ Concentration in 2100 (ppm)", xlabel="Shared Socioeconomic Pathway")
StatsPlots.scatter!(p, co2_conc.Scenario, co2_conc[!, :"2100"], label="RCMIP Value")
savefig(p, joinpath(@__DIR__, "..", "figures", "brick_ssp_conc_compare.png"))
