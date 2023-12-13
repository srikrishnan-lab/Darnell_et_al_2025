# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles
using DataFrames
using Makie, CairoMakie

# Assume that we call this script from the project folder
output_dir = "output/default"
slr_out = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))
# normalize relative to 2000
idx_2000 = findfirst(names(slr_out) .== "2000")
for row in axes(slr_out,1 )
    foreach(col -> slr_out[row, col] -= slr_out[row, idx_2000], axes(slr_out, 2))
end

milestones = [0.5, 1.0, 1.5, 2.0]

yrs = DataFrame(milestone=String[], year=Int64[])
for milestone in milestones
    idx = findfirst.(>=(milestone), eachrow(slr_out))
    filter!(!isnothing, idx)
    append!(yrs, DataFrame(milestone="$(milestone)m", year = parse.(Int, String.(idx))))
end

fig = Figure(resolution=(450, 450), fontsize=18, figure_padding=1)
ax = Axis(fig[1,1], xlabel="Year", ylabel="Sea Level Milestone (m)")
colors = Makie.wong_colors()
plt = Makie.rainclouds!(ax, yrs.milestone, yrs.year, color=colors[indexin(yrs.milestone, unique(yrs.milestone))], orientation=:horizontal, clouds=hist)

save("figures/slr-milestone.png", fig)
