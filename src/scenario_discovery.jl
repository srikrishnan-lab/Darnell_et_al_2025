#########################################################################
# scenario_discovery.jl                                                 #
#                                                                       #
# Plots factor maps for probabilities of exceeding thresholds.          # #                                                                       #
#                                                                       #
# This script requires the ensemble output to be present in             #
#   `results/peaking`.                                                  #
#                                                                       #
#########################################################################

# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles
using Statistics
using Plots
using DataFrames
using MLJ
using EvoTrees
using CairoMakie

# define random forest classifier object
classifier = @load EvoTreeClassifier pkg=EvoTrees

# load ensemble output
output_dir = "results/peaking"
slr_out = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
ais_out = DataFrame(CSVFiles.load(joinpath(output_dir, "antarctic.csv")))
gis_out = DataFrame(CSVFiles.load(joinpath(output_dir, "greenland.csv")))
temp_out = DataFrame(CSVFiles.load(joinpath(output_dir, "temperature.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))
# normalize relative to 2000
idx_2000 = findfirst(names(slr_out) .== "2000")
for row in axes(slr_out,1 )
    foreach(col -> slr_out[row, col] -= slr_out[row, idx_2000], axes(slr_out, 2))
end

## Start with SLR outcomes in 2100
features = parameters
idx_1850 = findfirst(names(slr_out) .== "1850")
idx_1900 = findfirst(names(slr_out) .== "1900")
features.antarctic_temp_threshold = [15.42 .+ 0.8365 * parameters[i, :antarctic_temp_threshold] - mean(temp_out[i, idx_1850:idx_1900]) for i in axes(temp_out, 1)]


function plot_feature_importance(year, threshold, features, subplot_label; f=Figure())
    slr_labels = ifelse.(threshold .< slr_out[:, Symbol(year)], "high", "normal")

    slr_key_tree = EvoTreeClassifier(nrounds=200, max_depth=3)
    slr_key_mach = machine(slr_key_tree, features, slr_labels)
    MLJ.fit!(slr_key_mach)

    param_import = stack(DataFrame(feature_importances(slr_key_mach)))
    ax = Axis(f[1,1], xticks = (1:10, param_import.variable[1:10]),                
        xticklabelrotation = pi/4, ylabel = "Feature Importance")
    Makie.barplot!(ax, param_import.value[1:10])

    Label(f[1, 1, TopLeft()], subplot_label, font=:bold,fontsize=16,
    padding = (0, 5, 5, 0),
    halign = :left)    

    return f
end

function plot_sos_contours(year, threshold, features, key_params, stepsize, limits, labels, subplot_label, colors; f=Figure())
    #refit tree using just t_peak, antarctic temp threshold, and climate sensitivity
    slr_labels = ifelse.(threshold .< slr_out[:, Symbol(year)], "high", "normal")

    slr_key_tree = EvoTreeClassifier(nrounds=200, max_depth=3)
    slr_key_mach = machine(slr_key_tree, features[:, key_params], slr_labels)
    MLJ.fit!(slr_key_mach)
    # set up grid for predictions
    bds1 = minimum(features[:,  key_params[2]]):stepsize[1]:maximum(features[:,  key_params[2]])
    bds2 = minimum(features[:,  key_params[3]]):stepsize[2]:maximum(features[:,  key_params[3]])
    coords = Iterators.product(bds1, bds2)

    # plot dividing contour
    key_features = [[2050.0, grid[1], grid[2]] for grid in coords]
    key_feature_df = DataFrame(mapreduce(permutedims, vcat, key_features), key_params)
    predict_key = MLJ.predict(slr_key_mach, key_feature_df)
    predict_class_2050 = [pred.prob_given_ref[1] for pred in predict_key]

    key_features = [[2070.0, grid[1], grid[2]] for grid in coords]
    key_feature_df = DataFrame(mapreduce(permutedims, vcat, key_features), key_params)
    predict_key = MLJ.predict(slr_key_mach, key_feature_df)
    predict_class_2070 = Float64.([pred.prob_given_ref[1] for pred in predict_key])


    key_features = [[2090.0, grid[1], grid[2]] for grid in coords]
    key_feature_df = DataFrame(mapreduce(permutedims, vcat, key_features), key_params)
    predict_key = MLJ.predict(slr_key_mach, key_feature_df)
    predict_class_2090 = Float64.([pred.prob_given_ref[1] for pred in predict_key])

    axmain = Axis(f[2,1], xlabel=labels[1], ylabel=labels[2])
    axtop = Axis(f[1,1], limits=(limits[1], (0, nothing)))
    axright = Axis(f[2, 2], limits=((0, nothing), limits[2]))
    linkyaxes!(axmain, axright)
    linkxaxes!(axmain, axtop)

    fmap = Makie.contour!(axmain, key_feature_df[:, 2], key_feature_df[:, 3], predict_class_2050, color=colors[0.2], levels=0.0:0.5:1.0, linewidth=4, xlabel=labels[1], ylabel=labels[2])
    fmap3 =Makie.contour!(axmain, key_feature_df[:, 2], key_feature_df[:, 3], predict_class_2070, color=colors[0.5], levels=0.0:0.5:1.0, linewidth=4)
    fmap5 =Makie.contour!(axmain, key_feature_df[:, 2], key_feature_df[:, 3], predict_class_2090, color=colors[0.5], levels=0.0:0.5:1.0, linewidth=4)


    # plot marginal densities for the geophysical uncertainties
    Makie.density!(axtop, features[!, key_params[2]])
    Makie.density!(axright, features[!, key_params[3]],direction=:y)
    hidedecorations!(axtop)
    hidedecorations!(axright)
    hidespines!(axtop) 
    hidespines!(axright)

    Label(f[1, 1, TopLeft()], subplot_label, font=:bold,fontsize=16,
    padding = (0, 5, 5, 5),
    halign = :left)    

    colsize!(f, 2, Auto(0.25))
    rowsize!(f, 1, Auto(0.25))
    colgap!(f, 1, Relative(0.0))
    rowgap!(f, 1, Relative(0.0))

    return f
end

f_imp = Figure(resolution=(700, 400), fontsize=12, figure_padding=(50, 10, 0, 30))

ga = f_imp[1, 1] = GridLayout()
gb = f_imp[1, 2] = GridLayout()

fig_imp_2100 = plot_feature_importance(2100, 1.0, features, "a"; f=ga)
fig_imp_2150 = plot_feature_importance(2150, 1.5, features, "b"; f=gb)

CairoMakie.save("figures/feature_importance.png", f_imp)

contour_colors = cgrad(:Reds_5)

f_sd = Figure(resolution=(700, 400), fontsize=16)

ga = f_sd[1, 1] = GridLayout()
gb = f_sd[1, 2] = GridLayout()

fig_2100_1m = plot_sos_contours(2100, 1.0, features,  [:t_peak, :climate_sensitivity, :antarctic_temp_threshold], [0.001, 0.001], [(1.5, 6), (1.2, 3.8)], ["Equilibrium Climate Sensitivity (°C)", "AIS Temperature Threshold (°C)"], "a", contour_colors; f=ga)

fig_2150_2m = plot_sos_contours(2150, 2.0, features,  [:t_peak, :climate_sensitivity, :antarctic_temp_threshold], [0.001, 0.001], [(1.5, 6), (1.2, 3.8)], ["Equilibrium Climate Sensitivity (°C)", "AIS Temperature Threshold (°C)"], "b", contour_colors; f=gb)

# create legend
elem_2050 = LineElement(color = contour_colors[0.0], linestyle = nothing, linewidth=4)
elem_2070 = LineElement(color = contour_colors[0.4], linestyle = nothing, linewidth=4)
elem_2090 = LineElement(color = contour_colors[0.8], linestyle = nothing, linewidth=4)

leg = Legend(f_sd[2,1:2], [elem_2050, elem_2070, elem_2090, ], ["2050", "2070", "2090"], "Year Emissions Peak", orientation=:horizontal, tellwidth=false, tellheight=true, framevisible=false)

CairoMakie.save("figures/factor_map_all.png", f_sd)
