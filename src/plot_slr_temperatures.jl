# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSVs
using XLSX # read XL files
using DataFrames # data structure for indices
using AlgebraOfGraphics
using AlgebraOfGraphics: density
using Makie # plotting library
using CairoMakie
using Measures # adjust margins with explicit measures
using StatsBase

output_dir = "results/default"
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
temperature = DataFrame(CSVFiles.load(joinpath(output_dir, "temperature.csv")))
gmslr = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
antarctic = DataFrame(CSVFiles.load(joinpath(output_dir, "antarctic.csv")))
gsic = DataFrame(CSVFiles.load(joinpath(output_dir, "gsic.csv")))
greenland = DataFrame(CSVFiles.load(joinpath(output_dir, "greenland.csv")))
lw_storage = DataFrame(CSVFiles.load(joinpath(output_dir, "lw_storage.csv")))
thermal_expansion = DataFrame(CSVFiles.load(joinpath(output_dir, "thermal_expansion.csv")))

function normalize_data!(dat, norm_yrs=nothing)
    # normalize to relevant period  defined by norm_yrs
    if !isnothing(norm_yrs)
        idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
        for row in axes(dat, 1)
            foreach(col -> dat[row, col] -= mean(dat[row, idx_norm]), axes(dat, 2))
        end
    end
    return dat
end

normalize_data!(temperature, 1850:1900)
normalize_data!(gmslr, [2000])
normalize_data!(antarctic, [2000])
normalize_data!(greenland, [2000])

function compute_norm_quantiles(dat, norm_yrs=nothing)
    # normalize to relevant period  defined by norm_yrs
    if !isnothing(norm_yrs)
        idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
        for row in axes(dat, 1)
            foreach(col -> dat[row, col] -= mean(dat[row, idx_norm]), axes(dat, 2))
        end
    end
    # compute median and 95% prediction interval
    quantiles = mapcols(col -> quantile(col, [0.05, 0.5, 0.95]), dat)
    return quantiles
end

emissions_q = compute_norm_quantiles(emissions)

## load CO2 emissions for  CMIP6 scenarios
cmip_df =  XLSX.readtable("data/cmip6_co2.xlsx", "data") |> DataFrame
# select only rows and columns with scenario names and data   
select!(cmip_df, Not([:Model, :Region, :Variable, :Unit, :Notes]))
cmip_df = cmip_df[1:7, :]
# reformat scenario names to SSPn-x.x
cmip_df[!, :Scenario] = replace.(cmip_df[!, :Scenario], r"(\d)(\d)" =>   s"\1.\2")
cmip_df[!, :Scenario] = [split(cmip_df[i, :Scenario], " ")[1] for i in 1:nrow(cmip_df)]
# sort by scenarios
sort!(cmip_df, :Scenario)
# convert emissions  from MtCO2/yr to GtCO2/yr
cmip_df[!, Not(:Scenario)] = cmip_df[!, Not(:Scenario)] ./ 1000

idx1850 = findfirst(names(temperature) .== "1850")
idx1900 = findfirst(names(temperature) .== "1900")

ais_threshold = [15.42 .+ 0.8365 * parameters[i, :antarctic_temp_threshold] - mean(temperature[i, idx1850:idx1900]) for i in axes(temperature, 1)]
ais_exceed_yr = Vector{Union{Float64, Nothing}}(undef, length(ais_threshold))

for j in eachindex(ais_threshold)
    exceed_ais = [temperature[j, k] > ais_threshold[j] for k in axes(temperature, 2)]
    exceed_ais_idx = findfirst(exceed_ais)
    if isnothing(exceed_ais_idx)
        ais_exceed_yr[j] = nothing
    else
        ais_exceed_yr[j] = parse(Int64, names(temperature)[exceed_ais_idx])
    end
end

ais_exceed = hcat(parameters[:, :t_peak], ais_exceed_yr)
ais_exceed[:, 2] = replace(ais_exceed[:, 2], nothing => 2305)

idx2100 = findfirst(names(gmslr) .== "2100")
idx2020 = findfirst(names(gmslr) .== "2020")
cum_emissions = [sum(emissions[i, idx2020:idx2100]) for i in 1:nrow(gmslr)]

slr_all = Vector(gmslr[:, idx2100])
dat = DataFrame(slr=slr_all, emissions=cum_emissions, ais_yr = ais_exceed[:, 2])
function assign_bin(yr)
    if yr < 2050
        return "Before 2050"
    elseif yr < 2075
        return "2050 -- 2075"
    elseif yr < 2100
        return "2075 -- 2100"
    else
        return "After 2100"
    end
end

dat.ais_bin = assign_bin.(dat.ais_yr)

fig = Figure(size=(600, 800), fontsize=14, figure_padding=10)

# set up layout
ga = fig[1, 1] = GridLayout()
gcd = fig[2, 1] = GridLayout()
gc = gcd[1, 1] = GridLayout()
gd = gcd[1, 2] = GridLayout()

axemissions = Axis(ga[1, 1], xlabel="Year", ylabel="CO₂ Emissions (Gt CO₂/yr)")
leg_med = Makie.lines!(axemissions, parse.(Int64, names(emissions_q)), Vector(emissions_q[2, :]), color=:black)
leg_ci = Makie.band!(axemissions, parse.(Int64, names(emissions_q)), Vector(emissions_q[1, :]), Vector(emissions_q[3, :]), color=(:gray, 0.2))
leg_ext = Makie.lines!(axemissions, parse.(Int64, names(emissions_q)), Vector((col -> maximum(col)).(eachcol(emissions))), color=:gray)
Makie.lines!(axemissions, parse.(Int64, names(emissions_q)), Vector((col -> minimum(col)).(eachcol(emissions))), color=:gray)
cmip_legend = Vector{LineElement}(undef, nrow(cmip_df)) # initialize storage for legend
cmip_yrs = parse.(Float64, string.(names(cmip_df[!, Not(:Scenario)])))
cmip_colors = cgrad(:Dark2_7, 7, categorical=true)
for i = 1:nrow(cmip_df)
    Makie.lines!(axemissions, cmip_yrs, Vector(cmip_df[i, Not(:Scenario)]), color=cmip_colors[i], linestyle=:dash)
    cmip_legend[i] = LineElement(color = cmip_colors[i], linestyle = :dash, linewidth=4)
end
Makie.xlims!(axemissions, 2000, 2300)

Legend(ga[2, 1], [[LineElement(color=:black), LineElement(color=:gray), PolyElement(color=:gray)], cmip_legend], [["Median", "Extrema", "95% Projection Interval"], cmip_df[!, :Scenario]], ["Ensemble Summary", "Emissions Scenario"], vertical=false, framevisible=false, tellheight=true, nbanks=4, titleposition=:top, tellwidth=false)

rowgap!(ga, 5)
rowsize!(ga, 1, Auto(0.85))
rowsize!(ga, 2, Auto(0.15))

# axcdf = Axis(fig[3, 1], ylabel="CDF", xlabel="CO₂ Emissions in 2100 (Gt CO₂/yr)")
# emis_range = 0:maximum(emissions[!, :"2100"])
# Makie.lines!(axcdf, emis_range, emis_cdf.(emis_range), color=:black)
# for i = 1:nrow(cmip_df)
#     Makie.vlines!(axcdf, cmip_df[i, end], color=cmip_colors[i], linestyle=:dash)
# end
# Label(fig[3, 1, TopLeft()], "b", fontsize=22, font=:bold, padding = (0, 50, 20, 0), halign=:right)

# axslr = Axis(gb[1,1], xlabel="Year", ylabel="Global Mean Sea Level (m)")
# color_slr = cgrad(:vik100, [0.1, 0.3, 0.6, 0.9])
idx_paris = findall(temperature[!, idx2100] .< 2.0) # find Paris consistent SOWs
# for i in idx_paris
#     Makie.lines!(axslr, parse.(Int64, names(gmslr)), Vector(gmslr[i, :]), color=cum_emissions[i], colorrange=(800, 6000), colormap=color_slr)
# end
# xlims!(2000, 2300)
# cbtrace = Colorbar(gb[1, 2], limits=(800, 6000), colormap=color_slr, label="Cumulative CO₂ Emissions (Gt CO₂)", vertical=true, flipaxis=false, tellheight=false)


Label(ga[1, 1, TopLeft()], "a", fontsize=18, font=:bold, padding = (0, 50, 10, 0), halign=:right)
# Label(gb[1, 1, TopLeft()], "b", fontsize=18, font=:bold, padding = (0, 50, 10, 0), halign=:right)

# rowsize!(gb, 1, Auto(0.9))
# rowsize!(gb, 2, Auto(0.1))

# rowsize!(gab, 1, Auto(1.3))
# rowgap!(gab, 5)

# axscatter = Axis(gc[1,1], xlabel="Global Temperature Anomaly (°C)", ylabel="Global Sea Level Anomaly (m)")
# colors = cgrad(:vik100, [0.1, 0.5, 0.9])
# plt =  Makie.scatter!(axscatter, temperature[:,idx2100  ], gmslr[:, idx2100], color=cum_emissions, colormap=colors, markersize=6)
# cbscatter = Colorbar(gc[1, 2], plt, label="Cumulative CO₂ Emissions (Gt CO₂)")


# axtrace = Axis(fig[3:4,2], xlabel="Year", ylabel="Global Sea Level Anomaly (m)")
# for i in idx_low_emissions
#     Makie.lines!(axtrace, 2020:2100, Vector(gmslr[i, idx2020:idx2100]), color=parameters[i, :antarctic_temp_threshold], colorrange=(2030, 2050), colormap=cgrad(:Reds))
# end
# xlims!(axtrace, (2020, 2100))
# axright = Axis(fig[3:4, 3], limits=((0, nothing), (0, 1.8)))
# linkyaxes!(axtrace, axright)
# Makie.density!(axright, gmslr[:, idx2100],direction=:y)
# hidedecorations!(axright)
# hidespines!(axright)

# cbtrace = Colorbar(fig[3:4, 4], limits=(2030, 2050), colormap=:Reds, label="Year In Which Emissions Peak")
# Label(fig[3, 2, TopLeft()], "d", fontsize=22, font=:bold, padding = (0, 50, 20, 0), halign=:right)

ais_labels = ["Before 2050", "2050 -- 2075", "2075 -- 2100", "After 2100"]
layers = density() * visual(Contour) + linear(interval=:confidence)
slr_gp = data(dat) *  mapping(:emissions => "Cumulative CO₂ Emissions (Gt CO₂)", :slr => "Global Sea Level Anomaly (m)") * layers  * mapping(color = :ais_bin => sorter(ais_labels) => "AIS Instability Trigger Year")
ais_trigger = draw!(gc[1, 1], slr_gp)
legend!(gc[2, 1], ais_trigger, tellwidth=false, tellheight=true, nbanks=2, framevisible=false)

Label(gc[1, 1, TopLeft()], "b", fontsize=18, font=:bold, padding = (0, 50, 0, 0), halign=:right)

axparis = Axis(gd[1, 1], xlabel="Cumulative CO₂ Emissions (Gt CO₂)", ylabel="Global Sea Level Anomaly (m)", alignmode=Inside())

colors_paris = cgrad(:vik100, [0.1, 0.5, 0.9], rev=true)
plt_paris =  Makie.scatter!(axparis, cum_emissions[idx_paris], gmslr[idx_paris, idx2100], color=ais_threshold[idx_paris], colormap=colors_paris, markersize=6)
cbscatter = Colorbar(gd[2, 1], plt_paris, label="AIS Temperature Threshold (°C)", vertical=false, flipaxis=true, tellwidth=false, tellheight=true)

Label(gd[1, 1, TopLeft()], "c", fontsize=18, font=:bold, padding = (0, 50, 0, 0), halign=:right)

#rowsize!(fig.layout, 1, Auto(0.5))
#rowsize!(fig.layout, 2, Auto(1.25))
colsize!(gcd, 1, Relative(0.5))

CairoMakie.save("figures/slr_temps.png", fig)





