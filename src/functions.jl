using CSVFiles
using XLSX
using Distributions
using DataFrames
using Random

# -------------------------------------------------------------------------------------------------------------------------------------------- #

#= create a function that simulates a stationary AR(1) process to superimpose noise onto the climate model projections
Function Arguments:
    num_years = number of time periods (years) the model is being run for
    σ         = calibrated standard deviation
    ρ         = calibrated autocorrelation coefficient
=#
function simulate_ar1_noise(num_years, σ, ρ)
    y = zeros(num_years)
    y[1] = rand(Normal(0, σ/sqrt(1-ρ^2)))
    ϵ = rand(Normal(0,σ), num_years) # time-varying observation errors
    for t in 2:num_years
        y[t] = ρ * y[t-1] + ϵ[t]
    end
    return y
end

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# create a function that returns a dataframe of historical emissions (1850-2021)
function historical_emissions(;start_year=1850, end_year=2300)
    # create vector of model years
    model_years = collect(start_year:end_year)

    # read in older historical observations (1850-2005) (from RCP files)
    RCP60 = scenario_df("RCP60")[in(model_years).(scenario_df("RCP60").year), :] # just picked RCP 6.0 scenario, but all scenarios have identical historical data from beginning-2005
    older_data = RCP60[start_year .<= RCP60.year .<= 2005, [:year, :rcp_co2_emissions]] # get historical data from 1850-2005 (same for all RCP scenarios)

    # read in more recent historical observations (2006-2021) (from Global Carbon Budget)
    recent_data = DataFrame(XLSX.readtable(joinpath(@__DIR__, "..", "data", "Global_Carbon_Budget_2022v1.0.xlsx"), "Global Carbon Budget", first_row=21)) # read in data
    recent_data = recent_data[2006 .<= recent_data.Year .<= 2021, [1,2,3]] # get emissions from years 2006-2021
    recent_data = DataFrame(year=recent_data.Year, rcp_co2_emissions=(recent_data[:,2] .+ recent_data[:,3]) .* 3.67) # add fossil emissions + land-use change emissions, and convert to units of GtCO₂

    # combine older and recent historical data
    historical_data = vcat(older_data, recent_data) # older_data is 1850-2005; recent_data is 2006-2021

    return historical_data
end

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# create a function that stores extreme RCP scenario emissions in a matrix (for plotting)
function rcp_emissions()

    # anthropogenic total CO2 emissions (PgC/yr) (Source: IPCC, 2013: Annex II: Climate System Scenario Tables)
    rcp26 = [2000 8.03;
            2010 9.70;
            2020 9.97;
            2030 8.00;
            2040 5.30;
            2050 3.50;
            2060 2.10;
            2070 0.81;
            2080 0.16;
            2090 -0.23;
            2100 -0.42]

    # anthropogenic total CO2 emissions (PgC/yr) (Source: IPCC, 2013: Annex II: Climate System Scenario Tables)
    rcp85 = [2000 8.03;
            2010 9.98;
            2020 12.28;
            2030 14.53;
            2040 17.33;
            2050 20.61;
            2060 23.83;
            2070 26.17;
            2080 27.60;
            2090 28.44;
            2100 28.77]

    # multiply by 3.67 to convert PgC/yr to GtCO₂/yr (1 PgC = 1 GtC)
    rcp26[:,2] = rcp26[:,2] * 3.67
    rcp85[:,2] = rcp85[:,2] * 3.67

    return rcp26, rcp85 # units of GtCO₂
end

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# create a function that extracts and stores relevant information for each RCP scenario into a dataframe
function scenario_df(rcp_scenario)

    # load emissions and forcing data
    rcp_emissions      = DataFrame(load(joinpath(@__DIR__, "..", "data", rcp_scenario*"_emissions.csv"), skiplines_begin=36))
    rcp_concentrations = DataFrame(load(joinpath(@__DIR__, "..", "data", rcp_scenario*"_concentrations.csv"), skiplines_begin=37))
    rcp_forcing        = DataFrame(load(joinpath(@__DIR__, "..", "data", rcp_scenario*"_midyear_radforcings.csv"), skiplines_begin=58))

    # index for years we want to consider
    year_idx = [1:1:736;] # rows 1-736
    years = rcp_emissions.YEARS[year_idx] # years 1765-2500 (full span of years in csv files)

    # calculate CO₂ emissions
    rcp_co2_emissions = (rcp_emissions.FossilCO2 .+ rcp_emissions.OtherCO2)[year_idx] * 3.67 # multiply by 3.67 to convert GtC/yr to GtCO₂/yr

    # get RCP N₂O concentrations (used for CO₂ radiative forcing calculations)
    rcp_n2o_concentration = rcp_concentrations.N2O[year_idx] # units in ppb

    # calculate exogenous radiative forcings
    rcp_aerosol_forcing = (rcp_forcing.TOTAER_DIR_RF .+ rcp_forcing.CLOUD_TOT_RF)[year_idx] # units in W/m₂
    rcp_other_forcing   = (rcp_forcing.TOTAL_INCLVOLCANIC_RF .- rcp_forcing.CO2_RF .- rcp_forcing.TOTAER_DIR_RF .- rcp_forcing.CLOUD_TOT_RF)[year_idx] # units in W/m₂

    df = DataFrame(year = years, # years 1765-2500
                   rcp_co2_emissions = rcp_co2_emissions, # units in GtCO₂/yr
                   rcp_n2o_concentration = rcp_n2o_concentration, # units in ppb
                   rcp_aerosol_forcing = rcp_aerosol_forcing, # units in W/m₂
                   rcp_other_forcing = rcp_other_forcing) # units in W/m₂
    return df
end

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# create a function to match N₂O concentrations, aerosol forcing, and other forcing to closest RCP scenario
function match_inputs(RCP26, RCP45, RCP60, RCP85, emissions_df)

    years = emissions_df.Year

    # initialize variables
    n2o_concentration = zeros(length(years))
    aerosol_forcing = zeros(length(years))
    other_forcing = zeros(length(years))

    for year in 1:nrow(emissions_df)

        # isolate emission info for each scenario
        true_val  = emissions_df[year, :Emissions]
        RCP26_val = RCP26[year, :rcp_co2_emissions]
        RCP45_val = RCP45[year, :rcp_co2_emissions]
        RCP60_val = RCP60[year, :rcp_co2_emissions]
        RCP85_val = RCP85[year, :rcp_co2_emissions]

        # put four RCP emissions values in an array
        RCP_array = [RCP26_val, RCP45_val, RCP60_val, RCP85_val]
        loc = findmin(abs.(RCP_array .- true_val))[2] # find index of minimum distance

        # initialize variables
        n2o_concentration_val = 0
        aerosol_forcing_val = 0
        other_forcing_val = 0

        # assign values based on closest RCP Scenario
        if loc == 1
            # use RCP 2.6 values
            n2o_concentration_val = RCP26[year, :rcp_n2o_concentration]
            aerosol_forcing_val = RCP26[year, :rcp_aerosol_forcing]
            other_forcing_val = RCP26[year, :rcp_other_forcing]
        elseif loc == 2
            # use RCP 4.5 values
            n2o_concentration_val = RCP45[year, :rcp_n2o_concentration]
            aerosol_forcing_val = RCP45[year, :rcp_aerosol_forcing]
            other_forcing_val = RCP45[year, :rcp_other_forcing]
        elseif loc == 3
            # use RCP 6.0 values
            n2o_concentration_val = RCP60[year, :rcp_n2o_concentration]
            aerosol_forcing_val = RCP60[year, :rcp_aerosol_forcing]
            other_forcing_val = RCP60[year, :rcp_other_forcing]
        elseif loc == 4
            # use RCP 8.5 values
            n2o_concentration_val = RCP85[year, :rcp_n2o_concentration]
            aerosol_forcing_val = RCP85[year, :rcp_aerosol_forcing]
            other_forcing_val = RCP85[year, :rcp_other_forcing]
        end

        n2o_concentration[year] = n2o_concentration_val
        aerosol_forcing[year] = aerosol_forcing_val
        other_forcing[year] = other_forcing_val
    end

    return n2o_concentration, aerosol_forcing, other_forcing
end

#=
Below code is for plotting the matched inputs to nearest scenario:
# create plot with line for CO₂ emissions associated with each RCP scenario
plot(RCP26.year, RCP26.rcp_co2_emissions, xlabel="Year", ylabel="CO₂ emissions (GtCO₂/year)", label="RCP 2.6")
plot!(RCP45.year, RCP45.rcp_co2_emissions, label="RCP 4.5")
plot!(RCP60.year, RCP60.rcp_co2_emissions, label="RCP 6.0")
plot!(RCP85.year, RCP85.rcp_co2_emissions, label="RCP 8.5")
# add line for DICE scenario
years = final_co2_df[:,:Year]
emissions = final_co2_df[:,:Emissions]
plot!(years, emissions, label="DICE Scenario")
# add dots for closest RCP scenario to the DICE scenario in question
scatter!(years, match_inputs(), markersize=2, label="Closest RCP")
#savefig(joinpath(@__DIR__, "input_plot.pdf"))
=#

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# create a function that generates one emissions curve given 3 sampled parameters (γ_g, t_peak, γ_d)
function emissions_curve(historical_data::DataFrame; γ_g=0.004, t_peak=2070, γ_d=0.07, start_year=1850, end_year=2300)

    # allocate space for years and emissions
    model_years = collect(1850:2300) # full range of years, will slice at the end if needed
    num_years = length(model_years)
    t = zeros(Int64, num_years) # years
    gtco2 = zeros(num_years) # emissions
    gtco2_tpeak = 0 # initialize value for emissions at peaking time   

    # separate historical data df into years and emissions
    historical_years = historical_data[:, :year]
    historical_emissions = historical_data[:, :rcp_co2_emissions] # units of GtCO₂

    # add historical emissions to arrays
    idx = length(historical_years)
    t[1:idx] = historical_years
    gtco2[1:idx] = historical_emissions

    for i = idx+1:length(t) # loop through years

        t[i] = t[i-1] + 1 # fill in current year to the t array

        # before peaking time: quadratic increase
        if t[i] <= t_peak
            gtco2[i] = gtco2[i-1] + γ_g * (t_peak - t[i]) # calculate emissions for current year
            
            # at peaking time: save value
            if t[i] == t_peak
                gtco2_tpeak = gtco2[i-1] + γ_g * (t_peak - t[i]) # save value of emissions at t_peak for scaling in logistic function
            end

        # after peaking time: logistic decrease
        elseif t[i] > t_peak
            gtco2[i] = gtco2[i-1] - ((2 * gtco2_tpeak * γ_d * exp(γ_d * (t[i] - t_peak))) / (exp(γ_d * (t[i] - t_peak)) + 1)^2) # calculate emissions for current year
        end

    end

    # return specified start and end years
    indices = findall((in)(start_year:end_year), 1850:2300)
    return t[indices], gtco2[indices] # return years and emissions

end

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# create a NetCDF file
function create_netcdf(data) # note: data is output_df
    ds = Dataset(joinpath(@__DIR__, "../results/slr_output.nc"), "c") # create file

    # add dimensions
    defDim(ds, "time", length(data.time))
    defDim(ds, "run", length(data.time)) # FIX LENGTH FOR RUN

    # add global attributes
    ds.attrib["creator"] = "Chloe Darnell"
    ds.attrib["title"] = "Sea Level Rise Data"
    ds.attrib["date"] = string(today())

    # add variables
    co2_emissions = defVar(ds, "co2_emissions", Float64, ("time",))
    radiative_forcing = defVar(ds, "radiative_forcing", Float64, ("time",))
    temperature = defVar(ds, "temperature", Float64, ("time",))
    sea_level_rise = defVar(ds, "sea_level_rise", Float64, ("time",))
    slr_antarctic_icesheet = defVar(ds, "slr_antarctic_icesheet", Float64, ("time",))
    slr_glaciers_small_ice_caps = defVar(ds, "slr_glaciers_small_ice_caps", Float64, ("time",))
    slr_greenland_icesheet = defVar(ds, "slr_greenland_icesheet", Float64, ("time",))
    slr_landwater_storage = defVar(ds, "slr_landwater_storage", Float64, ("time",))
    slr_thermal_expansion = defVar(ds, "slr_thermal_expansion", Float64, ("time",))

    # add variable attributes
    co2_emissions.attrib["units"] = "GtCO₂/yr"
    radiative_forcing.attrib["units"] = "W/m₂"
    temperature.attrib["units"] = "°C"
    sea_level_rise.attrib["units"] = "m"
    slr_antarctic_icesheet.attrib["units"] = "m"
    slr_glaciers_small_ice_caps.attrib["units"] = "m"
    slr_greenland_icesheet.attrib["units"] = "m"
    slr_landwater_storage.attrib["units"] = "m"
    slr_thermal_expansion.attrib["units"] = "m"

    # add data
    co2_emissions = data[:, :co2_emissions]
    radiative_forcing = data[:, :radiative_forcing]
    temperature = data[:, :temperature]
    sea_level_rise = data[:, :sea_level_rise]
    slr_antarctic_icesheet = data[:, :slr_antartic_icesheet]
    slr_glaciers_small_ice_caps = data[:, :slr_glaciers_small_ice_caps]
    slr_greenland_icesheet = data[:, :slr_greeland_icesheet]
    slr_landwater_storage = data[:, :slr_landwater_storage]
    slr_thermal_expansion = data[:, :slr_thermal_expansion]

    close(ds) # close the file
end