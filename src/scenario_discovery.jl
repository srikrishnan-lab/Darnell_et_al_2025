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
output_dir = joinpath(@__DIR__, "..", "results")
slr_default = DataFrame(CSVFiles.load(joinpath(output_dir, "default", "gmslr.csv")))
slr_optimistic = DataFrame(CSVFiles.load(joinpath(output_dir, "optimistic", "gmslr.csv")))
slr_pessimistic = DataFrame(CSVFiles.load(joinpath(output_dir, "pessimistic", "gmslr.csv")))
temp_out = DataFrame(CSVFiles.load(joinpath(output_dir, "default", "temperature.csv")))
parameters_default = DataFrame(CSVFiles.load(joinpath(output_dir, "default", "parameters.csv")))
parameters_optimistic = DataFrame(CSVFiles.load(joinpath(output_dir, "optimistic", "parameters.csv")))
parameters_pessimistic = DataFrame(CSVFiles.load(joinpath(output_dir, "pessimistic", "parameters.csv")))

# normalize relative to 2000
# function to find normalize relative to a year or mean over some normalization period)
function normalize_data!(dat, norm_yrs=nothing)
    # normalize to relevant period  defined by norm_yrs
    idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
    norm_mean = map(mean, eachrow(dat[:, idx_norm]))
    for row in axes(dat, 1)
        foreach(col -> dat[row, col] -= norm_mean[row], axes(dat, 2))
    end
    return dat
end

normalize_data!(slr_default, 1995:2014)
normalize_data!(slr_optimistic, 1995:2014)
normalize_data!(slr_pessimistic, 1995:2014)

## Start with SLR outcomes in 2100
idx_1850 = findfirst(names(slr_default) .== "1850")
idx_1900 = findfirst(names(slr_default) .== "1900")
parameters_default.antarctic_temp_threshold = [15.42 .+ 0.8365 * parameters_default[i, :antarctic_temp_threshold] - mean(temp_out[i, idx_1850:idx_1900]) for i in axes(temp_out, 1)]
parameters_optimistic.antarctic_temp_threshold = [15.42 .+ 0.8365 * parameters_optimistic[i, :antarctic_temp_threshold] - mean(temp_out[i, idx_1850:idx_1900]) for i in axes(temp_out, 1)]
parameters_pessimistic.antarctic_temp_threshold = [15.42 .+ 0.8365 * parameters_pessimistic[i, :antarctic_temp_threshold] - mean(temp_out[i, idx_1850:idx_1900]) for i in axes(temp_out, 1)]


function plot_feature_importance(slr_out, year, threshold, features, subplot_label, title; f=Figure())
    slr_labels = ifelse.(threshold .< slr_out[:, Symbol(year)], "high", "normal")

    slr_key_tree = EvoTreeClassifier(nrounds=200, max_depth=3)
    slr_key_mach = machine(slr_key_tree, features, slr_labels)
    MLJ.fit!(slr_key_mach)

    param_import = stack(DataFrame(feature_importances(slr_key_mach)))
    ax = Axis(f[1,1], xticks = (1:10, param_import.variable[1:10]),                
        xticklabelrotation = pi/4, ylabel = "Feature Importance", title=title)
    Makie.barplot!(ax, param_import.value[1:10])

    Label(f[1, 1, TopLeft()], subplot_label, font=:bold,fontsize=16,
    padding = (0, 5, 5, 0),
    halign = :left)    

    return f
end



function plot_sos_contours(slr_out, year, threshold, features, key_params, stepsize, limits, labels, subplot_label, colors, title; f=Figure())
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
    axmain.xticks = 2.0:1.0:6.0

    fmap = Makie.contour!(axmain, key_feature_df[:, 2], key_feature_df[:, 3], predict_class_2050, color=colors[0.2], levels=0.0:0.5:1.0, linewidth=4, xlabel=labels[1], ylabel=labels[2])
    fmap3 =Makie.contour!(axmain, key_feature_df[:, 2], key_feature_df[:, 3], predict_class_2070, color=colors[0.5], levels=0.0:0.5:1.0, linewidth=4)
    fmap5 =Makie.contour!(axmain, key_feature_df[:, 2], key_feature_df[:, 3], predict_class_2090, color=colors[0.8], levels=0.0:0.5:1.0, linewidth=4)

    # plot marginal densities for the geophysical uncertainties
    Makie.density!(axtop, features[!, key_params[2]])
    Makie.density!(axright, features[!, key_params[3]],direction=:y)
    hidedecorations!(axtop)
    hidedecorations!(axright)
    hidespines!(axtop) 
    hidespines!(axright)

    Label(f[1, 1, Top()], title, font=:bold,fontsize=16,
    halign = :center)    
    Label(f[1, 1, TopLeft()], subplot_label, font=:bold,fontsize=16,
    padding = (0, 5, 5, 5),
    halign = :left)    

    colsize!(f, 2, Auto(0.25))
    rowsize!(f, 1, Auto(0.25))
    colgap!(f, 1, Relative(0.0))
    rowgap!(f, 1, Relative(0.0))

    return f
end

f_imp = Figure(resolution=(1000, 400), fontsize=12, figure_padding=(50, 10, 0, 30))

ga = f_imp[1, 1] 
gb = f_imp[1, 2] 
gc = f_imp[1, 3]

fig_imp_default = plot_feature_importance(slr_default, 2100, 1.0, parameters_default, "a", "Baseline"; f=ga)
fig_imp_optimistic = plot_feature_importance(slr_default, 2100, 0.25, parameters_optimistic, "b", "Optimistic"; f=gb)
fig_imp_pessimistic = plot_feature_importance(slr_pessimistic, 2100, 0.5, parameters_pessimistic, "c", "Pessimistic"; f=gc)

CairoMakie.save("figures/feature_importance_scenarios.png", f_imp)

contour_colors = cgrad(:Reds_5)

f_sd = Figure(size=(1000, 400), fontsize=16)

ga = f_sd[1, 1] = GridLayout()
gb = f_sd[1, 2] = GridLayout()
gc = f_sd[1, 3] = GridLayout()

fig_2100_default = plot_sos_contours(slr_default, 2100, 0.25, parameters_default,  [:t_peak, :climate_sensitivity, :antarctic_temp_threshold], [0.001, 0.001], [(1.5, 6), (1.2, 3.8)], ["Equilibrium Climate Sensitivity (°C)", "AIS Temperature Threshold (°C)"], "a", contour_colors, "Baseline"; f=ga)

fig_2100_optimistic = plot_sos_contours(slr_optimistic, 2100, 0.5, parameters_optimistic,  [:t_peak, :climate_sensitivity, :antarctic_temp_threshold], [0.001, 0.001], [(1.5, 6), (1.2, 3.8)], ["Equilibrium Climate Sensitivity (°C)", "AIS Temperature Threshold (°C)"], "b", contour_colors, "Optimistic"; f=gb)

fig_2100_pessimistic = plot_sos_contours(slr_pessimistic, 2100, 1.0, parameters_pessimistic,  [:t_peak, :climate_sensitivity, :antarctic_temp_threshold], [0.001, 0.001], [(1.5, 6), (1.2, 3.8)], ["Equilibrium Climate Sensitivity (°C)", "AIS Temperature Threshold (°C)"], "c", contour_colors, "Pessimistic"; f=gc)

# create legend
elem_2050 = LineElement(color = contour_colors[0.0], linestyle = nothing, linewidth=4)
elem_2070 = LineElement(color = contour_colors[0.4], linestyle = nothing, linewidth=4)
elem_2090 = LineElement(color = contour_colors[0.8], linestyle = nothing, linewidth=4)

leg = Legend(f_sd[2,1:3], [elem_2050, elem_2070, elem_2090, ], ["2050", "2070", "2090"], "Year Emissions Peak", orientation=:horizontal, tellwidth=false, tellheight=true, framevisible=false)

CairoMakie.save(joinpath(@__DIR__, "..", "figures", "factor_map_all_scenarios.png"), f_sd)
