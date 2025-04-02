#########################################################################
# plot_slr_temperatures.jl                                              #
#                                                                       #
# Makes plot summarizing the output ensemble.                           # #                                                                       #
#                                                                       #
# This script requires the ensemble output to be present in             #
#   `results/default` .                                                 #
#                                                                       #
#########################################################################

# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSVs
using DataFrames # data structure for indices
using Makie # plotting library
using CairoMakie
using Measures # adjust margins with explicit measures
using StatsBase # get mean function and density
using GLM
using KernelDensity

# load ensemble
output_dir = "results/default"
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))
temperature = DataFrame(CSVFiles.load(joinpath(output_dir, "temperature.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
gmslr = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
ais = DataFrame(CSVFiles.load(joinpath(output_dir, "antarctic.csv")))

# define function to normalize data relative to some normalization period mean
function normalize_data!(dat, norm_yrs=nothing)
    # normalize to relevant period  defined by norm_yrs
    idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
    norm_mean = map(mean, eachrow(dat[:, idx_norm]))
    for row in axes(dat, 1)
        foreach(col -> dat[row, col] -= norm_mean[row], axes(dat, 2))
    end
    return dat
end

normalize_data!(temperature, 1850:1900)
normalize_data!(gmslr, 1995:2014)

# find cumulative emissions from 2022--2100
idx2100 = findfirst(names(gmslr) .== "2100")
idx2000 = findfirst(names(gmslr) .== "2000")
idx1850 = findfirst(names(gmslr) .== "1850")
idx1900 = findfirst(names(gmslr) .== "1900")

avg_temp_2100 = temperature[:, idx2100] - temperature[:, idx2000]  # find emissions average
temp_diff = temperature[:, idx2000+1:idx2100] .- temperature[:, idx2000]
temp_int = sum(Matrix(temp_diff), dims=2) ./ 100
cum_emissions = reduce(hcat, map(cumsum, eachrow(emissions[:, idx2000+1:idx2100])))'
emis_int = sum(Matrix(cum_emissions), dims=2) ./ 100

avg_gmslr_2100 = (gmslr[:, idx2100] - gmslr[:, idx2000]) * 1000 / (2100 - 2000 + 1)

# compute densities
dens_gmslr_2100 = kde(hcat(vec(temp_int), avg_gmslr_2100))
dens_emis_2100 = kde(hcat(vec(emis_int), avg_gmslr_2100))

ais_threshold = [15.42 .+ 0.8365 * parameters[i, :antarctic_temp_threshold] - mean(temperature[i, idx1850:idx1900]) for i in axes(temperature, 1)]

ais_exceed_yr = Vector{Union{Float64, Nothing}}(undef, length(ais_threshold))
for j in eachindex(ais_threshold)
    exceed_ais = [temperature[j, k] > ais_threshold[j] for k in axes(temperature, 2)]
    if sum(exceed_ais[idx1850:idx2100]) == 0
        ais_exceed_yr[j] = 0
    else
        ais_exceed_yr[j] = 1
    end
end

function add_breakpoint(dat, bp)
    pred_var = names(dat)[1]
    dat[!, "after_bp"] = max.(0, dat[!, pred_var] .- bp)
    return (pred_var, dat)
end

function fit_piecewise(dat, minbp, maxbp, step)

    min_deviance = Inf
    best_model = nothing
    best_bp = 0
    current_model = nothing
    
    for bp in minbp:step:maxbp
      pred_var, dat_bp = add_breakpoint(dat, bp)
      current_model = lm(@formula(slr ~ temp + after_bp), dat_bp)
      if deviance(current_model) < min_deviance
        min_deviance = deviance(current_model)
        best_model = current_model
        best_bp = bp
      end
    end
    
    return best_model, best_bp
  end

gmslr_dat = DataFrame(temp=vec(temp_int), slr=avg_gmslr_2100)
emis_dat = DataFrame(temp=vec(emis_int), slr=avg_gmslr_2100)

temp_pred_range = 0:0.1:4
emis_pred_range = 1000:100:4000

temp_lm_all = fit_piecewise(gmslr_dat, 0, 4, 0.05)
temp_predict_all = disallowmissing!(GLM.predict(temp_lm_all[1], add_breakpoint(DataFrame(temp=temp_pred_range), temp_lm_all[2])[2], interval=:confidence, level=0.95))
emis_lm_all = fit_piecewise(emis_dat, 1000, 4000, 100)
emis_predict_all = disallowmissing!(GLM.predict(emis_lm_all[1], add_breakpoint(DataFrame(temp=emis_pred_range), emis_lm_all[2])[2], interval=:confidence, level=0.95))

temp_bp_idx = findfirst(temp_pred_range .== temp_lm_all[2])
emis_bp_idx = findfirst(emis_pred_range .== emis_lm_all[2])

temp_lm_nonexceed = fit_piecewise(gmslr_dat[ais_exceed_yr .== 0, :], -15, 60, 0.05)
temp_predict_nonexceed = disallowmissing!(GLM.predict(temp_lm_nonexceed[1], add_breakpoint(DataFrame(temp=temp_pred_range), temp_lm_nonexceed[2])[2], interval=:confidence, level=0.95))
emis_lm_nonexceed = fit_piecewise(emis_dat[ais_exceed_yr .== 0, :], 1000, 4000, 100)
emis_predict_nonexceed = disallowmissing!(GLM.predict(emis_lm_nonexceed[1], add_breakpoint(DataFrame(temp=emis_pred_range), emis_lm_nonexceed[2])[2], interval=:confidence, level=0.95))

fig = Figure(size=(700, 500), fontsize=16, figure_padding=10)
# set up layout
ga = fig[1, 1] = GridLayout()
gb = fig[1, 2] = GridLayout()

ax_gmt = Makie.Axis(ga[1,1])
ax_main = Axis(ga[2, 1], xlabel="Avg GMST Time-Integral (°C)", ylabel="GSLR Rate (mm/yr)", alignmode=Inside())
ax_gmslr = Makie.Axis(ga[2,2])

# link axes and hide decorations for marginal plots
linkyaxes!(ax_main, ax_gmslr)
linkxaxes!(ax_main, ax_gmt)
hidedecorations!(ax_gmslr)
hidedecorations!(ax_gmt)
hidespines!(ax_gmslr)
hidespines!(ax_gmt)
# plot marginal densities
Makie.density!(ax_gmt, vec(temp_int))
Makie.density!(ax_gmslr, avg_gmslr_2100, direction=:y)

colsize!(ga, 2, Auto(0.25))
rowsize!(ga, 1, Auto(0.25))
colgap!(ga, 1, Relative(-0.01))
rowgap!(ga, 1, Relative(-0.01))

Makie.contour!(ax_main, dens_gmslr_2100.x, dens_gmslr_2100.y, dens_gmslr_2100.density, levels=10)
Makie.band!(ax_main, temp_pred_range, temp_predict_all[!, :lower], temp_predict_all[!, :upper], color=:orange, alpha=0.3, label=false)
lin_all = Makie.lines!(ax_main, 0:0.1:temp_lm_all[2], temp_predict_all[1:length(0:0.1:temp_lm_all[2]), :prediction], color=:orange, linewidth=2, label="All Simulations")
Makie.lines!(ax_main, temp_lm_all[2]:0.1:4, temp_predict_all[length(0:0.1:temp_lm_all[2]):end, :prediction], color=:orange, linewidth=2, label=false, linestyle=:dash)
Makie.band!(ax_main, temp_pred_range, temp_predict_nonexceed[!, :lower], temp_predict_nonexceed[!, :upper], color=:blue, alpha=0.3, label=false)
lin_nonexceed = Makie.lines!(ax_main, temp_pred_range, temp_predict_nonexceed[!, :prediction], color=:blue, linewidth=2, label="Neglecting Fast Dynamics")
Makie.scatter!(ax_main, [temp_lm_all[2]], [temp_predict_all[temp_bp_idx, :prediction]], color=:orange, label=false)


Makie.xlims!(ax_main, -0.04, 2)
Makie.ylims!(ax_gmslr, 0, 20)


ax_gmt = Makie.Axis(gb[1,1])
ax_main2 = Axis(gb[2, 1], xlabel="Avg CO₂ Emissions Time-Integral (GtCO₂)", ylabel="GMSLR Rate (mm/yr)", alignmode=Inside())
ax_ais = Makie.Axis(gb[2,2])

# link axes and hide decorations for marginal plots
linkyaxes!(ax_main2, ax_ais)
linkxaxes!(ax_main2, ax_gmt)
hidedecorations!(ax_ais)
hidedecorations!(ax_gmt)
hidespines!(ax_ais)
hidespines!(ax_gmt)
# plot marginal densities
Makie.density!(ax_gmt, vec(emis_int))
Makie.density!(ax_ais, avg_gmslr_2100, direction=:y)

colsize!(gb, 2, Auto(0.25))
rowsize!(gb, 1, Auto(0.25))
colgap!(gb, 1, Relative(-0.01))
rowgap!(gb, 1, Relative(-0.01))

Makie.contour!(ax_main2, dens_emis_2100.x, dens_emis_2100.y, dens_emis_2100.density, levels=10)
Makie.band!(ax_main2, emis_pred_range, emis_predict_all[!, :lower], emis_predict_all[!, :upper], color=:orange, alpha=0.3, label=false)
Makie.lines!(ax_main2, 1000:100:emis_lm_all[2], emis_predict_all[1:emis_bp_idx, :prediction], color=:orange, linewidth=2, label="All Simulations")
Makie.lines!(ax_main2, emis_lm_all[2]:100:4000, emis_predict_all[emis_bp_idx:end, :prediction], color=:orange, linewidth=2, label=false, linestyle=:dash)
Makie.band!(ax_main2, emis_pred_range, emis_predict_nonexceed[!, :lower], emis_predict_nonexceed[!, :upper], color=:blue, alpha=0.3, label=false)
Makie.lines!(ax_main2, emis_pred_range, emis_predict_nonexceed[!, :prediction], color=:blue, linewidth=2, label="Neglecting Fast Dynamics")
Makie.scatter!(ax_main2, [emis_lm_all[2]], [emis_predict_all[emis_bp_idx, :prediction]], color=:orange)

Makie.xlims!(ax_main2, 1000, 4000)
Makie.ylims!(ax_ais, 0, 20)

Label(ga[1, 1, TopLeft()], "a", fontsize=18, font=:bold, padding = (0, 50, 0, 0), halign=:right)
Label(gb[1, 1, TopLeft()], "b", fontsize=18, font=:bold, padding = (0, 50, 0, 0), halign=:right)

Legend(fig[2, 1:2], [lin_all, lin_nonexceed], ["All", "No Fast Dynamics"], "Simulations", framevisible=true, orientation=:horizontal)

rowgap!(fig.layout, 1, Relative(0.05))
resize_to_layout!(fig)

CairoMakie.save("figures/slr_temps.png", fig)
