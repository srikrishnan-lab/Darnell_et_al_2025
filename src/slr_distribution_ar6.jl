using Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

using NetCDF
using CSVFiles
using XLSX
using DataFrames
using Statistics
using Distributions
using StatsBase

# estimate quantiles for the RCPs
# function to find the 66% prediction interval (optional: after normalizing relative to a year or mean over some normalization period)
function norm_df(dat, ret_yr, norm_yrs)
    # normalize to relevant period  defined by norm_yrs
    idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
    mean_val = map(row -> mean(row[idx_norm]), eachrow(dat))
    for row in axes(dat, 1)
        foreach(col -> dat[row, col] -= mean_val[row], axes(dat, 2))
    end
    idx_ret = findfirst((names(dat) .== string(ret_yr)))
    return dat[:, idx_ret]
end

function read_brick_output(variable, scenario, ret_yr=2100, norm_yrs=nothing)
    brick_raw = DataFrame(CSVFiles.load(joinpath(@__DIR__, "results", scenario, "$(variable).csv")))
    brick_out = norm_df(brick_raw, ret_yr, norm_yrs)
    return brick_out
end

# normalize SLR relative to 1995-2014 and compute quantiles
slr_norm = 1995:2014
temp_norm = 1850:1900

# get slr simulation output
slr_out_2100 = read_brick_output("gmslr", "default", 2100, slr_norm)
slr_out_2150 = read_brick_output("gmslr", "default", 2150, slr_norm)
slr_out_2200 = read_brick_output("gmslr", "default", 2200, slr_norm)
temp_out_2100 = read_brick_output("temperature", "default", 2100, temp_norm)
temp_out_2150 = read_brick_output("temperature", "default", 2150, temp_norm)
temp_out_2200 = read_brick_output("temperature", "default", 2200, temp_norm)

temps_ipcc = [1.5, 2.0, 3.0, 4.0, 5.0]
temps_q = ecdf(temp_out).(temps_ipcc)
function ipcc_ensemble(nsample, yr, temps_q)
    # read in IPCC AR6 files
    ar6_path = joinpath(@__DIR__, "data", "ar6")
    ipcc_nc_15 = "$(ar6_path)/global/confidence_output_files/medium_confidence/tlim1.5win0.25/total_tlim1.5win0.25_medium_confidence_values.nc"
    ipcc_nc_20 = "$(ar6_path)/global/confidence_output_files/medium_confidence/tlim2.0win0.25/total_tlim2.0win0.25_medium_confidence_values.nc"
    ipcc_nc_30 = "$(ar6_path)/global/confidence_output_files/medium_confidence/tlim3.0win0.25/total_tlim3.0win0.25_medium_confidence_values.nc"
    ipcc_nc_40 = "$(ar6_path)/global/confidence_output_files/medium_confidence/tlim4.0win0.25/total_tlim4.0win0.25_medium_confidence_values.nc"
    ipcc_nc_50 = "$(ar6_path)/global/confidence_output_files/medium_confidence/tlim5.0win0.25/total_tlim5.0win0.25_medium_confidence_values.nc"
    ipcc_out_15 = ncread(ipcc_nc_15, "sea_level_change")
    ipcc_out_20 = ncread(ipcc_nc_20, "sea_level_change")
    ipcc_out_30 = ncread(ipcc_nc_30, "sea_level_change")
    ipcc_out_40 = ncread(ipcc_nc_40, "sea_level_change")
    ipcc_out_50 = ncread(ipcc_nc_50, "sea_level_change")
    
    # get quantile and year indices
    ipcc_q = round.(Float64.(ncread(ipcc_nc_15, "quantiles")); digits=2)a
    ipcc_yr = findfirst(ncread(ipcc_nc_15, "years") .== yr)

    # sample 
    ar6_out = zeros(nsample)
    u = rand(Uniform(0, 1), nsample) # CDF location
    v = rand(Uniform(0, 1), nsample) # AR6 projection quantile
    idx_q = [findfirst(ipcc_q .== round.(q; digits=2)) for q in v]
    for i = 1:nsample
        if u[i] < (temps_q[1] + temps_q[2]) / 2
            ar6_out[i] = ipcc_out_15[:, ipcc_yr, idx_q[i]][1]
        elseif u[i] < (temps_q[2] + temps_q[3]) / 2
            ar6_out[i] = ipcc_out_20[:, ipcc_yr, idx_q[i]][1]
        elseif u[i] < (temps_q[3] + temps_q[4]) / 2
            ar6_out[i] = ipcc_out_30[:, ipcc_yr, idx_q[i]][1]
        elseif u[i] < (temps_q[4] + temps_q[5]) / 2
            ar6_out[i] = ipcc_out_40[:, ipcc_yr, idx_q[i]][1]
        else
            ar6_out[i] = ipcc_out_50[:, ipcc_yr, idx_q[i]][1]
        end
    end

    return ar6_out ./ 1000
end

ipcc_yrs = [2050, 2060, 2070, 2080, 2090, 2100, 2110, 2120, 2130, 2140, 2150]

ar6_all = reduce(hcat, [ipcc_ensemble(100_000, yr, temps_q) for yr in ipcc_yrs])
CSVFiles.save(joinpath(@__DIR__, "output", "ar6", "gmslr.csv"), DataFrame(ar6_all, string.(ipcc_yrs)))