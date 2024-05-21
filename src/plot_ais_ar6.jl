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

using StatsBase
using CSVFiles
using DataFrames
using Colors
using Makie
using CairoMakie

# AR6 scenario SLR ranges (Ch.9, Table 9.9, 66% ranges, m relative to 1995-2014)

ar6_slr = DataFrame(
    year = repeat([2030, 2050, 2090, 2100, 2150], inner=6),
    scenario = repeat(["SSP1-1.9", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5", "SSP5-8.5 LC"], outer=5),
    min = [0.08, 0.08, 0.08, 0.08, 0.09, 0.09, 0.15, 0.16, 0.17, 0.18, 0.20, 0.20, 0.26, 0.30, 0.38, 0.46, 0.52, 0.52, 0.28, 0.32, 0.44, 0.55, 0.63, 0.63, 0.37, 0.46, 0.66, 0.89, 0.98, 0.98],
    max = [0.12, 0.12, 0.12, 0.12, 0.12, 0.15, 0.23, 0.25, 0.26, 0.27, 0.29, 0.40, 0.49, 0.54, 0.65, 0.74, 0.83, 1.30, 0.55, 0.62, 0.76, 0.90, 1.01, 1.60, 0.86, 0.99, 1.33, 1.65, 1.88, 4.82]
)

ar6_temp = DataFrame(
    year = repeat([2030, 2050, 2090], inner=5),
    scenario = repeat(["SSP1-1.9", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"], outer=3),
    min = [1.2, 1.2, 1.2, 1.2, 1.3, 1.2, 1.3, 1.6, 1.7, 1.9, 1.0, 1.3, 2.1, 2.8, 3.3],
    median = [1.5, 1.5, 1.5, 1.5, 1.6, 1.6, 1.7, 2.0, 2.1, 2.4, 1.4, 1.8, 2.7, 3.6, 4.4],
    max = [1.7, 1.8, 1.8, 1.8, 1.9, 2.0, 2.2, 2.5, 2.6, 3.0, 1.8, 2.4, 3.5, 4.6, 5.7]
)

# colors selected for consistency with AR6 Fig. 9.27
ar6_lines_927 = Dict("SSP1-1.9" => colorant"rgb( 35, 197, 226)",
                     "SSP1-2.6" => colorant"rgb( 21, 73, 135)",
                     "SSP2-4.5" => colorant"rgb(255, 150,  60)",
                     "SSP3-7.0" => colorant"rgb(255,  66,  48)",
                     "SSP5-8.5" => colorant"rgb(178,  44,  26)",
                     "SSP5-8.5 LC" => colorant"rgb(200, 35, 200)"
)
# add colors to AR6 dataframe
ar6_slr.color = [ar6_lines_927[s] for s in ar6_slr.scenario]
ar6_temp.color = [ar6_lines_927[s] for s in ar6_temp.scenario]

# add dodge values for plotting
plt_dodge = -2.5:1.0:2.5
ar6_slr.plt_x = ar6_slr.year + repeat(plt_dodge, outer=5)
plt_dodge = -0.5:1.0:3.5
ar6_temp.plt_x = ar6_temp.year + repeat(plt_dodge, outer=3)

# read in SLR and temperature output
slr_default = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "default", "gmslr.csv")))
slr_optimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results",  "optimistic", "gmslr.csv")))
slr_pessimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results",  "pessimistic", "gmslr.csv")))
temp_default = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results", "default", "temperature.csv")))
temp_optimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results",  "optimistic", "temperature.csv")))
temp_pessimistic = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "results",  "pessimistic", "temperature.csv")))

# normalize SLR relative to 1995-2014 and compute quantiles
years = parse.(Int64, names(slr_default))
slr_norm = 1995:2014
temp_norm = 1850:1900

# function to find the 66% prediction interval (optional: after normalizing relative to a year or mean over some normalization period)
function compute_norm_quantiles(dat; norm_yrs=nothing, mean_yrs=nothing)
    # normalize to relevant period  defined by norm_yrs
    if !isnothing(norm_yrs)
        idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
        mean_val = map(row -> mean(row[idx_norm]), eachrow(dat))
        for row in axes(dat, 1)
            foreach(col -> dat[row, col] -= mean_val[row], axes(dat, 2))
        end
    end

    if !isnothing(mean_yrs)
        quantiles = DataFrame(near=zeros(3), medium=zeros(3), far=zeros(3))
        for i in eachindex(mean_yrs)
            idx_mean = findall((!isnothing).(indexin(names(dat), string.(mean_yrs[i]))))
            quantiles[:, i] = quantile(map(mean, eachrow(dat[:, idx_mean])), [0.05, 0.5, 0.95])
        end
    else
            # compute median and 90% prediction interval
        quantiles = mapcols(col -> quantile(col, [0.17, 0.5, 0.83]), dat)
    end
    return quantiles
end

default_slr_q = compute_norm_quantiles(slr_default, norm_yrs=slr_norm)
optimistic_slr_q = compute_norm_quantiles(slr_optimistic, norm_yrs=slr_norm)
pessimistic_slr_q = compute_norm_quantiles(slr_pessimistic, norm_yrs=slr_norm)

default_temp_q = compute_norm_quantiles(temp_default, norm_yrs=temp_norm, mean_yrs=(2021:2040, 2041:2060, 2081:2100))
optimistic_temp_q = compute_norm_quantiles(temp_optimistic, norm_yrs=temp_norm, mean_yrs=(2021:2040, 2041:2060, 2081:2100))
pessimistic_temp_q = compute_norm_quantiles(temp_pessimistic, norm_yrs=temp_norm, mean_yrs=(2021:2040, 2041:2060, 2081:2100))


# make plot
plt_yrs = 2000:2150
plt_idx = indexin(plt_yrs, years)
fig = Figure(size=(1000, 800), fontsize=20, figure_padding=10)
ax_slr = Axis(fig[1,1:2], xlabel="Year", ylabel="Global Mean Sea Level Relative\n to 1995-2014 Mean (m)")
# plot our scenarios
Makie.band!(ax_slr, years[plt_idx], Vector(default_slr_q[1, plt_idx]), Vector(default_slr_q[3, plt_idx]), color=(:black, 0.2), label="Baseline Scenario")
Makie.lines!(ax_slr, years[plt_idx], Vector(default_slr_q[2, plt_idx]), color=:black, linewidth=2, label="Baseline Scenario")
Makie.band!(ax_slr, years[plt_idx], Vector(optimistic_slr_q[1, plt_idx]), Vector(optimistic_slr_q[3, plt_idx]), color=(:darkorange, 0.2), label="Optimistic Scenario")
Makie.lines!(ax_slr, years[plt_idx], Vector(optimistic_slr_q[2, plt_idx]), color=:darkorange, linewidth=2, label="Optimistic Scenario")
Makie.band!(ax_slr, years[plt_idx], Vector(pessimistic_slr_q[1, plt_idx]), Vector(pessimistic_slr_q[3, plt_idx]), color=(:teal, 0.1), label="Pessimistic Scenario")
Makie.lines!(ax_slr, years[plt_idx], Vector(pessimistic_slr_q[2, plt_idx]), color=:teal, linewidth=2, label="Pessimistic Scenario")
# plot AR6 ranges
Makie.rangebars!(ax_slr, ar6_slr.plt_x, ar6_slr.min, ar6_slr.max, color=ar6_slr.color, label=Vector(ar6_slr.scenario), linewidth=2)
ylims!(ax_slr, (0, 2.75))

ax_temp = Axis(fig[2,1], xlabel="Year", ylabel="Temperature Anomaly Relative\n to 1850-1990 Mean (°C)")
# plot our scenarios

# plot AR6 ranges
Makie.rangebars!(ax_temp, [2026.5, 2046.5, 2086.5], Vector(default_temp_q[1, :]), Vector(default_temp_q[3, :]), color=:black, linewidth=2)
Makie.scatter!(ax_temp, [2026.5, 2046.5, 2086.5], Vector(default_temp_q[2, :]), color=:black, markersize=8)
Makie.rangebars!(ax_temp, [2027.5, 2047.5, 2087.5], Vector(optimistic_temp_q[1, :]), Vector(optimistic_temp_q[3, :]), color=:darkorange, linewidth=2)
Makie.scatter!(ax_temp, [2027.5, 2047.5, 2087.5], Vector(optimistic_temp_q[2, :]), color=:darkorange, markersize=8)
Makie.rangebars!(ax_temp, [2028.5, 2048.5, 2088.5], Vector(pessimistic_temp_q[1, :]), Vector(pessimistic_temp_q[3, :]), color=:teal, linewidth=2)
Makie.scatter!(ax_temp, [2028.5, 2048.5, 2088.5], Vector(pessimistic_temp_q[2, :]), color=:teal, markersize=8)
Makie.rangebars!(ax_temp, ar6_temp.plt_x,  ar6_temp.min, ar6_temp.max, color=ar6_temp.color, linewidth=2)
ylims!(ax_temp, 0, 5.75)

Makie.poly!(ax_temp, Rect(2023, 0.7, 15, 0.2), color=:lightgray)
Makie.poly!(ax_temp, Point2f[(2021, 0.8), (2023, 0.95), (2023, 0.65)], color=:lightgray)
Makie.poly!(ax_temp, Point2f[(2040, 0.8), (2038, 0.95), (2038, 0.65)], color=:lightgray)

Makie.poly!(ax_temp, Rect(2043, 0.7, 15, 0.2), color=:lightgray)
Makie.poly!(ax_temp, Point2f[(2041, 0.8), (2043, 0.95), (2043, 0.65)], color=:lightgray)
Makie.poly!(ax_temp, Point2f[(2060, 0.8), (2058, 0.95), (2058, 0.65)], color=:lightgray)

Makie.poly!(ax_temp, Rect(2083, 0.7, 15, 0.2), color=:lightgray)
Makie.poly!(ax_temp, Point2f[(2081, 0.8), (2083, 0.95), (2083, 0.65)], color=:lightgray)
Makie.poly!(ax_temp, Point2f[(2100, 0.8), (2098, 0.95), (2098, 0.65)], color=:lightgray)
Makie.text!(ax_temp, Point.([2030, 2050, 2090], 0.325), text=["Near-Term\n (2021-2040)", "Mid-Term\n (2041-2060)", "Far-Term\n (2081-2100)"], color=:black, align=(:center, :center), fontsize=15)

# build legend
def_leg = [LineElement(color=:black, linewidth=2), PolyElement(color=(:black, 0.2)), MarkerElement(marker=:circle, color=:black, markersize=7)]
opt_leg = [LineElement(color=:darkorange, linewidth=2), PolyElement(color=(:darkorange, 0.2)), MarkerElement(marker=:circle, color=:darkorange, markersize=7)]
pes_leg = [LineElement(color=:teal, linewidth=2), PolyElement(color=(:teal, 0.2)), MarkerElement(marker=:circle, color=:teal, markersize=7)]
ar6_leg = Vector{LineElement}(undef, length(unique(ar6_slr.scenario)))
for i in eachindex(unique(ar6_slr.scenario))
    idx = findfirst(ar6_slr.scenario .== unique(ar6_slr.scenario)[i])
    ar6_leg[i] = LineElement(color=ar6_slr.color[idx])
end

Legend(fig[2, 2], reduce(vcat, ([def_leg, opt_leg, pes_leg], ar6_leg)), reduce(vcat, (["Baseline Scenario", "Optimistic Scenario","Pessimistic Scenario"], unique(ar6_slr.scenario))), orientation=:vertical, tellwidth=true, tellheight=false)


rowgap!(fig.layout, 1, Fixed(50))
colsize!(fig.layout, 2, Auto(0.35))

CairoMakie.save(joinpath(@__DIR__, "..", "figures", "slr_ar6.png"), fig)