# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSV of Shapley indices
using DataFrames # data structure for indices
using Plots # plotting library
using ColorSchemes
using Measures # adjust margins with explicit measures
using Statistics # get mean function

# assume this is called from the project root directory
output_path = "output/shapley"
shap_ind = DataFrame(CSVFiles.load(joinpath(output_path, "shapley_indices.csv")))

# assign parameters to groups by subsystem
groups = Dict("Emissions"         => ["gamma_g", "t_peak", "gamma_d"],
            "Carbon Cycle"        => ["CO2_0", "N2O_0", "Q10",  "CO2_fertilization", "CO2_diffusivity"],
            "Climate System"      => ["temperature_0", "ocean_heat_0",  "heat_diffusivity", "rf_scale_aerosol", "climate_sensitivity"],
            "Thermal Expansion"   => ["thermal_alpha", "thermal_s0"],
            "Greenland Ice Sheet"           => ["greenland_a", "greenland_b", "greenland_alpha", "greenland_beta", "greenland_v0"],
            "Antarctic Ice Sheet"           => ["antarctic_s0", "antarctic_gamma", "antarctic_alpha", "antarctic_mu", "antarctic_nu",
                                        "antarctic_precip0", "antarctic_kappa", "antarctic_flow0", "antarctic_runoff_height0",
                                        "antarctic_c", "antarctic_bed_height0", "antarctic_slope", "antarctic_lambda",
                                        "antarctic_temp_threshold", "anto_alpha", "anto_beta"],             
            "Glaciers and SIC"    => ["glaciers_beta0", "glaciers_n", "glaciers_v0", "glaciers_s0"],
            "Land Water Storage"  => ["lw_random_sample"],            
            "Other"   => ["sd_temp", "sd_ocean_heat", "sd_glaciers", "sd_greenland", "sd_antarctic", "sd_gmsl",
            "rho_temperature", "rho_ocean_heat", "rho_glaciers", "rho_greenland", "rho_antarctic", "rho_gmsl"],
)
group_col = Array{String}(undef, size(shap_ind)[1]) # preallocate space
for (key,value) in groups # loop through dictionary and create vector with group labels
    indices = findall((in)(value), shap_ind.feature_name) # find indices for current group
    group_col[indices] .= key # fill the indices with the name of current group
end
shap_ind.group = group_col # add group column to dataframe   

# average shapley effects by group for each year
shap_group = groupby(shap_ind[!, Not(:feature_name)], :group)
shap_group = combine(shap_group, Not(:group) .=> sum)
# normalize so the grouped shapley sums equal 1
shap_norm = mapcols(x -> x / sum(x), shap_group[!, Not([:group])])
insertcols!(shap_norm, 1, :group => shap_group.group)
group_order = ["Emissions", "Carbon Cycle", "Climate System", "Greenland Ice Sheet", "Antarctic Ice Sheet", "Glaciers and SIC", "Thermal Expansion", "Land Water Storage", "Other"]
shap_norm = shap_norm[indexin(group_order, shap_group.group), :]
shap_permute = permutedims(shap_norm, 1)

group_colors = ColorSchemes.Set1_9[1:9]

# plot indices over time
yrs = 2050:5:2200
plt = areaplot(yrs, Matrix(shap_permute[!, Not(:group)]), label=permutedims(shap_norm.group), xlabel="Year", ylabel="Normalized Grouped Shapley Index", color_palette=group_colors, left_margin=5mm, right_margin=5mm, bottom_margin=10mm, legendfontsize=10, guidefontsize=10, tickfontsize=9, legend=:outerright)
Plots.xticks!(plt, 2050:25:2200)
Plots.xlims!((2050, 2200))
Plots.ylims!((0, 1))
plot!(size=(800, 400))

#shap_ind_norm = DataFrame(mapcols(x -> x ./ sum(x), shap_ind[!, Not([:feature_name, :group])]), names(shap_group)[2:end])
#shap_ind_norm.feature_name = shap_ind.feature_name
#shap_ind_norm.group = shap_ind.group
#shap_ind_2100 = sort(shap_ind_norm[!, [:feature_name, :mean_2100_sum, :group]], :mean_2100_sum, rev=true)[1:5, :]
#sort!(shap_ind_2100, :group, rev=true)
#shap_ind_2150 = sort(shap_ind_norm[!, [:feature_name, :mean_2150_sum, :group]], :mean_2150_sum, rev=true)[1:5, :]
#sort!(shap_ind_2150, :group, rev=true)

# p_ind2100 = bar(shap_ind_2100[!, :feature_name], shap_ind_2100[!, :mean_2100_sum],group=shap_ind_2100[!, :group], color_palette=group_colors[[1, 3, 5]], xlabel="Parameter", ylabel="Normalized Shapley Index", left_margin=5mm, right_margin=5mm, bottom_margin=30mm, legend=:false, guidefontsize=10, tickfontsize=10, xrotation = 45)

# p_ind2150 = bar(shap_ind_2150[!, :feature_name], shap_ind_2150[!, :mean_2150_sum],group=shap_ind_2150[!, :group], xlabel="Parameter", ylabel="Normalized Shapley Index", color_palette=group_colors[[1, 3, 5]], left_margin=5mm, right_margin=5mm, legend=:false, guidefontsize=10, tickfontsize=8, xrotation = 45)

savefig(plt, "figures/stacked-shapley-index.png")
