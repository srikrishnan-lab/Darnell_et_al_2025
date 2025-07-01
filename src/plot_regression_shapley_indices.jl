#########################################################################
# plot_regression_shapley_indices.jl                                    #
#                                                                       #
# Plots grouped Shapley indices over time.                              # #                                                                       #
#                                                                       #
# This script requires the Shapley output to be in `output/shapley`.    #
# If this is not there, or if a different ensemble was generated,       #
#    run `src/regression_drivers.jl` to recompute.                      #
#                                                                       #
#########################################################################

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
using StatsPlots
using CategoricalArrays

# assume this is called from the project root directory
output_path = joinpath(@__DIR__, "..", "output", "shapley")

yrs = 2050:5:2200
group_colors = ColorSchemes.glasbey_hv_n256[1:9]

shap_ind_default = DataFrame(CSVFiles.load(joinpath(output_path, "shapley_indices_default.csv")))

function normalize_shap_groups(shap_ind)
    # assign parameters to groups by subsystem
    groups = Dict("Emissions"         => ["γ_g", "t_peak", "γ_d"],
    "Carbon Cycle"        => ["CO2_0", "N2O_0", "Q10",  "CO2_fertilization", "CO2_diffusivity"],
    "Climate System"      => ["temperature_0", "ocean_heat_0",  "heat_diffusivity", "rf_scale_aerosol", "climate_sensitivity"],
    "Thermal Expansion"   => ["thermal_alpha", "thermal_s0"],
    "Greenland Ice Sheet"           => ["greenland_a", "greenland_b", "greenland_alpha", "greenland_beta", "greenland_v0"],
    "Antarctic Ice Sheet"           => ["antarctic_s0", "antarctic_gamma", "antarctic_alpha", "antarctic_mu", "antarctic_nu",
                                "antarctic_precip0", "antarctic_kappa", "antarctic_flow0", "antarctic_runoff_height0",
                                "antarctic_c", "antarctic_bed_height0", "antarctic_slope", "antarctic_lambda",
                                "antarctic_temp_threshold", "anto_alpha", "anto_beta"],             
    "Glaciers and Small Ice Caps"    => ["glaciers_beta0", "glaciers_n", "glaciers_v0", "glaciers_s0"],
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
    group_order = ["Emissions", "Carbon Cycle", "Climate System", "Greenland Ice Sheet", "Antarctic Ice Sheet", "Glaciers and Small Ice Caps", "Thermal Expansion", "Land Water Storage", "Other"]
    shap_norm = shap_norm[indexin(group_order, shap_group.group), :]
    shap_permute = permutedims(shap_norm, 1)
    return shap_permute
end
shap_default = normalize_shap_groups(shap_ind_default)

# plot indices over time
inch = 96
mmx = inch / 25.4
p_shap = areaplot(yrs, Matrix(shap_default[!, Not(:group)]), xlabel="Year", ylabel="Relative Group Importance", color_palette=group_colors, guidefontsize=7, tickfontsize=6, legend=:outerbottom, label=permutedims(names(shap_default)[2:end]), legendfontsize=7, fg_color_legend=false, rightmargin=5mm)
annotate!(p_shap, 2030, 1.05, text("a", :left, 8, "Helvetica Bold"))
Plots.xticks!(p_shap, 2050:25:2200)
Plots.xlims!(p_shap, (2050, 2200))
Plots.ylims!(p_shap, (0, 1))

# plot group importance differences betweend other scenarios and default
shap_ind_optimistic = DataFrame(CSVFiles.load(joinpath(output_path, "shapley_indices_optimistic.csv")))
shap_ind_pessimistic = DataFrame(CSVFiles.load(joinpath(output_path, "shapley_indices_pessimistic.csv")))
shap_optimistic = normalize_shap_groups(shap_ind_optimistic)
shap_pessimistic = normalize_shap_groups(shap_ind_pessimistic)
shap_opt_diff = shap_optimistic[:, 2:end] .- shap_default[:, 2:end]
shap_pess_diff = shap_pessimistic[:, 2:end] .- shap_default[:, 2:end]

function plot_differences(shap_opt_diff, shap_pess_diff, groupname, label)
    diff_mat = Matrix(hcat(shap_opt_diff[!, groupname], shap_pess_diff[!, groupname]))
    plt = Plots.plot(yrs, diff_mat, linewidth=2, xlabel="Year", ylabel="Difference From\nBaseline", title=groupname, legend=:false, colors=[:orange, :teal], grid=false,  guidefontsize=7, tickfontsize=6, legendfontsize=7, titlefontsize=8, rightmargin=5mm, bottommargin=5mm)
    annotate!(plt, 1950, maximum(diff_mat) + (maximum(diff_mat) - minimum(diff_mat)) / 8, text(label, :left, 8, "Helvetica Bold"))
    xticks!(plt, 2050:50:2200)
    hline!(plt, [0.0], color=:black, linestyle=:dash)
    return plt
end
p_emis = plot_differences(shap_opt_diff, shap_pess_diff, "Emissions", "b")
p_clim = plot_differences(shap_opt_diff, shap_pess_diff, "Climate System", "c")
p_ais = plot_differences(shap_opt_diff, shap_pess_diff, "Antarctic Ice Sheet", "d")
p_gis = plot_differences(shap_opt_diff, shap_pess_diff, "Greenland Ice Sheet", "e")
p_leg = plot((-2:-1)', (-2:-1)', lims=(0, 0.05), legendfontsize=8, legend=:bottom, fg_color_legend=false, labels=["Optimistic" "Pessimistic"], fc=[:orange :teal], frame=:none)
l = @layout [
            a{0.5w} [grid(2, 2)
                     b{0.005h}    ]
]
plt = Plots.plot(p_shap, p_emis, p_clim, p_ais, p_gis, p_leg, layout=l, dpi=300, size=(180mmx, 90mmx))

savefig(plt, joinpath(@__DIR__, "..", "figures", "stacked-shapley-index-scenarios.pdf"))

# make supplemental figure for barplots in selected years
yrlabel = repeat(["2050", "2060", "2075", "2100"], inner=ncol(shap_default)-1)
cat = CategoricalArray(repeat(names(shap_default)[2:end], outer=4))
levels!(cat, names(shap_default)[2:end])
yrplt = groupedbar(yrlabel, Matrix(shap_default[[1, 3, 6, 11], 2:end])', group=cat, barposition=:dodge, xlabel="Uncertainty Group", ylabel="Normalized Shapley Value", ylims=(0, 0.6), palette=group_colors, legend=:topleft)
savefig(yrplt, joinpath(@__DIR__, "..", "figures", "shapley-index-years.png"))