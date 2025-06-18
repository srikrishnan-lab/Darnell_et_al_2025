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
using ColorSchemes
using Measures # adjust margins with explicit measures
using StatsBase # get mean function and density
using GLM
using KernelDensity

# load ensemble
output_dir = "results/default"
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))
temperature = DataFrame(CSVFiles.load(joinpath(output_dir, "temperature.csv")))
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

# define function to calculate cumulative sums of 
function cumulative_sum(dat)
    cum_df = zeros(size(temp))
    for i in 1:ncol(dat)
        if i == 1
            cum_df[:, i] = dat[:, i]
        else
            cum_df[:, i] = dat[:, i] + cum_df[:, i-1]
        end
    end
    return cum_df
end

normalize_data!(temperature, 1850:1900)
normalize_data!(gmslr, 1995:2014)

# find cumulative emissions from 2022--2100
idx2100 = findfirst(names(gmslr) .== "2100")
idx2000 = findfirst(names(gmslr) .== "2000")
idx1850 = findfirst(names(gmslr) .== "1850")
idx1900 = findfirst(names(gmslr) .== "1900")

# find temperature and GMSLR century increase
avg_temp_2100 = temperature[:, idx2100] - temperature[:, idx2000]  
temp_diff = temperature[:, idx2000+1:idx2100] .- temperature[:, idx2000]
temp_int = sum(Matrix(temp_diff), dims=2) ./ 100
avg_gmslr_2100 = (gmslr[:, idx2100] - gmslr[:, idx2000])

temp = temperature[:, idx2000:idx2100]
temp_integral = cumulative_sum(temp)

# compute densities
dens_gmslr_2100 = kde(hcat(vec(temp_int), avg_gmslr_2100))

# compute AIS threshold exceedance years
ais_threshold = [15.42 .+ 0.8365 * parameters[i, :antarctic_temp_threshold] - mean(temperature[i, idx1850:idx1900]) for i in axes(temperature, 1)]

exceed_ais = zeros(size(temp))
#ais_exceed_yr = Vector{Union{Int64, Nothing}}(undef, length(ais_threshold))
for j in eachindex(ais_threshold)
    exceed_ais_yr = [temp[j, k] > ais_threshold[j] for k in axes(temp, 2)]
    if sum(exceed_ais_yr) > 0
        yr_idx = findfirst(exceed_ais_yr .== 1)
        exceed_ais[j, yr_idx:end] .= 1
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

gmslr_dat = DataFrame(temp=vec(temp_int), slr=vec(avg_gmslr_2100), exceed=vec(sum(eachcol(exceed_ais)) .> 0))

function fit_and_predict(dat, pred_range, step)
    lm_fit = fit_piecewise(dat, pred_range[1], pred_range[end], step)
    pred = disallowmissing!(GLM.predict(lm_fit[1], add_breakpoint(DataFrame(temp=pred_range), lm_fit[2])[2], interval=:prediction, level=0.95))
    return (lm_fit, pred)
end

# fit model for 2100 temps/GMSLR
temp_pred_range = round(minimum(gmslr_dat[:, 1]); digits=0):0.05:round(maximum(gmslr_dat[:, 1]); digits=0)

temp_lm_all, temp_predict_all = fit_and_predict(gmslr_dat, temp_pred_range, 0.05)
temp_bp_idx = findfirst(temp_pred_range .== temp_lm_all[2])
temp_lm_nonexceed, temp_predict_nonexceed = fit_and_predict(gmslr_dat[gmslr_dat.exceed .== 0, :], 0:0.05:4, 0.05)

# fit model for time-integrated temperature
gmslr_dat = DataFrame(temp=vec(temp_integral), slr=vec(Matrix(gmslr[:, idx2000:idx2100])), exceed=vec(exceed_ais))

int_pred_range = round(minimum(gmslr_dat[:, 1]); digits=0):1:round(maximum(gmslr_dat[:, 1]); digits=0)
int_pred_nonexceed_range = round(minimum(gmslr_dat[gmslr_dat.exceed .== 0, 1]); digits=0):1:round(maximum(gmslr_dat[gmslr_dat.exceed .== 0, 1]); digits=0)
int_pred_exceed_range = round(minimum(gmslr_dat[gmslr_dat.exceed .== 1, 1]); digits=0):1:round(maximum(gmslr_dat[gmslr_dat.exceed .== 1, 1]); digits=0)

int_lm_all, int_predict_all = fit_and_predict(gmslr_dat, int_pred_range, 5)
int_bp_idx = findfirst(int_pred_range .== int_lm_all[2])
int_lm_nonexceed, int_predict_nonexceed = fit_and_predict(gmslr_dat[gmslr_dat.exceed .== 0, :],  int_pred_nonexceed_range, 5)
int_lm_exceed, int_predict_exceed = fit_and_predict(gmslr_dat[gmslr_dat.exceed .== 1, :],  int_pred_exceed_range, 5)

# make plot
colors = ColorSchemes.tol_bright[[6, 3, 2]]
fig = Figure(size=(700, 500), fontsize=16, figure_padding=10)
# set up layout
ga = fig[1, 1] = GridLayout()
gb = fig[1, 2] = GridLayout()

ax_main = Axis(ga[1, 1], xlabel="GMST Anomaly (°C/century)", ylabel="GMSLR (m/century)", alignmode=Inside())

Makie.contour!(ax_main, dens_gmslr_2100.x, dens_gmslr_2100.y, dens_gmslr_2100.density, levels=10)
Makie.band!(ax_main, temp_pred_range, temp_predict_all[!, :lower], temp_predict_all[!, :upper], color=colors[1], alpha=0.3, label=false)
Makie.lines!(ax_main, temp_pred_range, temp_predict_all[:, :prediction], color=colors[1], linewidth=2, label="All Simulations")
Makie.band!(ax_main, temp_pred_nonexceed_range, temp_predict_nonexceed[!, :lower], temp_predict_nonexceed[!, :upper], color=colors[2], alpha=0.3, label=false)
Makie.lines!(ax_main, temp_pred_nonexceed_range, temp_predict_nonexceed[!, :prediction], color=colors[2], linewidth=2, label="Neglecting Fast Dynamics")
Makie.scatter!(ax_main, [temp_lm_all[2]], [temp_predict_all[temp_bp_idx, :prediction]], color=colors[1], label=false)


#Makie.xlims!(ax_main, -0.04, 2)

ax_main2 = Axis(gb[1, 1], xlabel="GMST Time-Integral (°C-yr)", ylabel="GMSLR (m)", alignmode=Inside())

Makie.band!(ax_main2, int_pred_range, int_predict_all[!, :lower], int_predict_all[!, :upper], color=colors[1], alpha=0.2, label=false)
lin_all = Makie.lines!(ax_main2, int_pred_range, int_predict_all[:, :prediction], color=colors[1], linewidth=2, label="All Simulations")
Makie.band!(ax_main2, int_pred_nonexceed_range, int_predict_nonexceed[!, :lower], int_predict_nonexceed[!, :upper], color=colors[2], alpha=0.2, label=false)
lin_nonexceed = Makie.lines!(ax_main2, int_pred_nonexceed_range, int_predict_nonexceed[!, :prediction], color=colors[2], linewidth=2, label="No Fast Dynamics")
Makie.band!(ax_main2, int_pred_exceed_range, int_predict_exceed[!, :lower], int_predict_exceed[!, :upper], color=colors[3], alpha=0.2, label=false)
lin_exceed = Makie.lines!(ax_main2, int_pred_exceed_range, int_predict_exceed[!, :prediction], color=colors[3], linewidth=2, label="Triggering Fast Dynamics")
Makie.scatter!(ax_main2, [int_lm_all[2]], [int_predict_all[int_bp_idx, :prediction]], color=colors[1], label=false)


Label(ga[1, 1, TopLeft()], "a", fontsize=18, font=:bold, padding = (0, 50, 0, 0), halign=:right)
Label(gb[1, 1, TopLeft()], "b", fontsize=18, font=:bold, padding = (0, 50, 0, 0), halign=:right)

Legend(fig[2, 1:2], [lin_all, lin_nonexceed, lin_exceed], ["All", "No Fast Dynamics", "Triggered Fast Dynamics"], "Simulations", framevisible=true, orientation=:horizontal)

resize_to_layout!(fig)

CairoMakie.save("figures/slr_temps.png", fig)
