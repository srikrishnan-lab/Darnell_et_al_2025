# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSV of Shapley indices
using DataFrames # data structure for indices
using PrettyTables # table formatting
using Crayons # crayons for highlighting

ssp_path = joinpath(@__DIR__, "results", "ssp")
ssps = readdir(ssp_path)

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

ar6_temp = DataFrame(
    year = repeat([2030, 2050, 2090], inner=5),
    scenario = repeat(["SSP1-1.9", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"], outer=3),
    min = [1.2, 1.2, 1.2, 1.2, 1.3, 1.2, 1.3, 1.6, 1.7, 1.9, 1.0, 1.3, 2.1, 2.8, 3.3],
    median = [1.5, 1.5, 1.5, 1.5, 1.6, 1.6, 1.7, 2.0, 2.1, 2.4, 1.4, 1.8, 2.7, 3.6, 4.4],
    max = [1.7, 1.8, 1.8, 1.8, 1.9, 2.0, 2.2, 2.5, 2.6, 3.0, 1.8, 2.4, 3.5, 4.6, 5.7]
)

colnames = ["Scenario", "Years", "AR6", "SNEASY"]
temp_comp = DataFrame(Matrix(undef, 0, length(colnames)), colnames)
formatted_ssps = ["SSP1-1.9", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]
for (i, ssp) in pairs(["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"])
    brick_dat = DataFrame(CSVFiles.load(joinpath(ssp_path, ssp, "temperature.csv")))
    brick_q = round.(compute_norm_quantiles(brick_dat, norm_yrs=1850:1900, mean_yrs=(2021:2040, 2041:2060, 2081:2100)); digits=1)
    ar6_dat = filter(:scenario => isequal(formatted_ssps[i]), ar6_temp)
    push!(temp_comp, [formatted_ssps[i], "2021-2040", "$(ar6_dat[1, :median]) ($(ar6_dat[1, :min]), $(ar6_dat[1, :max]))", "$(brick_q[2, :near]) ($(brick_q[1, :near]), $(brick_q[3, :near]))"])
    push!(temp_comp, [formatted_ssps[i], "2041-2060", "$(ar6_dat[2, :median]) ($(ar6_dat[2, :min]), $(ar6_dat[2, :max]))", "$(brick_q[2, :medium]) ($(brick_q[1, :medium]), $(brick_q[3, :medium]))"])
    push!(temp_comp, [formatted_ssps[i], "2081-2100", "$(ar6_dat[3, :median]) ($(ar6_dat[3, :min]), $(ar6_dat[3, :max]))", "$(brick_q[2, :far]) ($(brick_q[1, :far]), $(brick_q[3, :far]))"])

end

PrettyTables.pretty_table(temp_comp, backend=Val(:latex))