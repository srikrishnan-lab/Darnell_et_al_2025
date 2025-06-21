# activate the environment
using Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

## import necessary packages
using Random # ability to set seeds
using DataFrames # tabular data structure
using Distributions # API for working with statistical distribution
using CSVFiles # need to keep for "load" function with CSVs
using Interpolations
using Mimi
using MimiBRICK

data_path = joinpath(@__DIR__, "..", "data")

## load emissions scenarios from SSP database
emis_dat = DataFrame(CSVFiles.load(joinpath(data_path, "rcmip-emissions-annual-means-v5-1-0.csv")))
conc_dat = DataFrame(CSVFiles.load(joinpath(data_path, "rcmip-concentrations-annual-means-v5-1-0.csv")))
forc_dat = DataFrame(CSVFiles.load(joinpath(data_path, "rcmip-radiative-forcing-annual-means-v5-1-0.csv")))
filter!(:Mip_Era => isequal("CMIP6"), emis_dat)
filter!(:Mip_Era => isequal("CMIP6"), conc_dat)
filter!(:Mip_Era => isequal("CMIP6"), forc_dat)
scenario_regexp = r"^ssp[0-9]*$"
filter!(:Scenario => contains(scenario_regexp), emis_dat)
filter!(:Scenario => contains(scenario_regexp), conc_dat)
filter!(:Scenario => contains(scenario_regexp), forc_dat)
filter!(:Region => isequal("World"), emis_dat)
filter!(:Region => isequal("World"), conc_dat)
filter!(:Region => isequal("World"), forc_dat)
filter!(:Variable => isequal("Emissions|CO2"), emis_dat)
co2_conc = filter(:Variable => isequal("Atmospheric Concentrations|CO2"), conc_dat)
n2o_conc = filter(:Variable => isequal("Atmospheric Concentrations|N2O"), conc_dat)
aerosol_forc = filter(:Variable => isequal("Effective Radiative Forcing|Anthropogenic|Aerosols"), forc_dat)
co2_forc = filter(:Variable => isequal("Effective Radiative Forcing|Anthropogenic|CO2"), forc_dat)
other_forc = filter(:Variable => isequal("Effective Radiative Forcing"), forc_dat)

# interpolate co2 emissions on an annual grid using a linear interpolation
for (i, row) in enumerate(eachrow(emis_dat))
    xs_idx = findall((!ismissing).(collect(emis_dat[i, 8:end])))
    yrs = parse.(Int, names(emis_dat)[8:end])
    xs = yrs[xs_idx]
    interp = LinearInterpolation(xs, collect(emis_dat[i, xs_idx .+ 7]))
    emis_dat[i, 8:end] = interp.(yrs)
end

ssps = unique(co2_conc[!, :Scenario])

# load posterior samples
params = DataFrame(CSVFiles.load(joinpath(@__DIR__, "..", "data", "calibrated_parameters/parameters_subsample_sneasybrick.csv")))

## set up MimiBRICK to run with SSP emissions
function run_model(A, co2_emis, n2o_concentration, aerosol_forcings, other_forcings, start_year, end_year, output_dir)
    # create MimiBRICK instance
    m = MimiBRICK.create_sneasy_brick(start_year=start_year, end_year=end_year) 
    num_samples = size(A, 1)
    model_years = collect(start_year:end_year)
    num_years = length(model_years)
    # add land-water storage rate samples
    A[!, :lw_random_sample] = rand(Normal(0.0003, 0.00018), num_samples)

    # pre-allocate arrays to store results
    parameters                  = zeros(Float64, num_samples, 36)
    co2_emissions               = zeros(Float64, num_samples, num_years)
    co2_concentrations          = zeros(Float64, num_samples, num_years)
    co2_concentrations_no_noise = zeros(Float64, num_samples, num_years)
    oceanco2                    = zeros(Union{Missing, Float64}, num_samples, num_years)
    radiative_forcing           = zeros(Float64, num_samples, num_years)
    temperature                 = zeros(Float64, num_samples, num_years)
    global_mean_sea_level_rise  = zeros(Float64, num_samples, num_years)
    slr_antarctic_icesheet      = zeros(Float64, num_samples, num_years)
    slr_glaciers_small_ice_caps = zeros(Float64, num_samples, num_years)
    slr_greenland_icesheet      = zeros(Float64, num_samples, num_years)
    slr_landwater_storage       = zeros(Float64, num_samples, num_years)
    slr_thermal_expansion       = zeros(Float64, num_samples, num_years)
    ocean_heat                  = zeros(Float64, num_samples, num_years)
    # loop over each sample input and evaluate the model
    for i = 1:num_samples
        gtco2 = co2_emis / 1000 / 3.67 # divide by (1000 * 3.67) to convert MtCO₂/yr to GtC/yr (SNEASY needs input of GtC/yr)

        # feed CO₂ emissions and forcings into SNEASY-BRICK
        update_param!(m, :ccm, :CO2_emissions, gtco2)
        update_param!(m, :rfco2, :N₂O, n2o_concentration)
        update_param!(m, :radiativeforcing, :rf_aerosol, aerosol_forcings)
        update_param!(m, :radiativeforcing, :rf_other, other_forcings)

        # ---------------------------------------------------- Create and Update SNEASY-BRICK Parameters --------------------------------------------------------- #

        # create parameters
        σ_temperature               = A[i,1]    # statistical noise parameter
        σ_ocean_heat                = A[i,2]    # statistical noise parameter
        σ_glaciers                  = A[i,3]    # statistical noise parameter
        σ_greenland                 = A[i,4]    # statistical noise parameter
        σ_antarctic                 = A[i,5]    # statistical noise parameter
        σ_gmsl                      = A[i,6]    # statistical noise parameter
        σ_co2                       = A[i,7]
        ρ_temperature               = A[i,8]   # statistical noise parameter
        ρ_ocean_heat                = A[i,9]   # statistical noise parameter
        ρ_glaciers                  = A[i,10]   # statistical noise parameter
        ρ_greenland                 = A[i,11]   # statistical noise parameter
        ρ_antarctic                 = A[i,12]   # statistical noise parameter
        ρ_gmsl                      = A[i,13]   # statistical noise parameter
        α₀_CO₂                      = A[i,14]
        CO2_0                       = A[i,15]
        N2O_0                       = A[i,16]
        temperature_0               = A[i,17]   # statistical noise parameter (will not use)
        ocean_heat_0                = A[i,18]   # statistical noise parameter
        thermal_s0                  = A[i,19]
        greenland_v0                = A[i,20]
        glaciers_v0                 = A[i,21]
        glaciers_s0                 = A[i,22]
        antarctic_s0                = A[i,23]
        Q10                         = A[i,24]
        CO2_fertilization           = A[i,25]
        CO2_diffusivity             = A[i,26]
        heat_diffusivity            = A[i,27]
        rf_scale_aerosol            = A[i,28]
        climate_sensitivity         = A[i,29]
        thermal_alpha               = A[i,30]
        greenland_a                 = A[i,31]
        greenland_b                 = A[i,32]
        greenland_alpha             = A[i,33]
        greenland_beta              = A[i,34]
        glaciers_beta0              = A[i,35]
        glaciers_n                  = A[i,36]
        anto_alpha                  = A[i,37]
        anto_beta                   = A[i,38]
        antarctic_gamma             = A[i,39]
        antarctic_alpha             = A[i,40]
        antarctic_mu                = A[i,41]
        antarctic_nu                = A[i,42]
        antarctic_precip0           = A[i,43]
        antarctic_kappa             = A[i,44]
        antarctic_flow0             = A[i,45]
        antarctic_runoff_height0    = A[i,46]
        antarctic_c                 = A[i,47]
        antarctic_bed_height0       = A[i,48]
        antarctic_slope             = A[i,49]
        antarctic_lambda            = A[i,50]
        antarctic_temp_threshold    = A[i,51]
        lw_random_sample            = A[i,52]

        # ----- Land Water Storage ----- #
        update_param!(m, :landwater_storage, :lws_random_sample, fill(lw_random_sample, num_years))

        # ----- Antarctic Ocean ----- #
        update_param!(m, :antarctic_ocean, :anto_α, anto_alpha)
        update_param!(m, :antarctic_ocean, :anto_β, anto_beta)

        # ----- Antarctic Ice Sheet ----- #
        update_param!(m, :antarctic_icesheet, :ais_sea_level₀, antarctic_s0)
        update_param!(m, :antarctic_icesheet, :ais_bedheight₀, antarctic_bed_height0)
        update_param!(m, :antarctic_icesheet, :ais_slope, antarctic_slope)
        update_param!(m, :antarctic_icesheet, :ais_μ, antarctic_mu)
        update_param!(m, :antarctic_icesheet, :ais_runoffline_snowheight₀, antarctic_runoff_height0)
        update_param!(m, :antarctic_icesheet, :ais_c, antarctic_c)
        update_param!(m, :antarctic_icesheet, :ais_precipitation₀, antarctic_precip0)
        update_param!(m, :antarctic_icesheet, :ais_κ, antarctic_kappa)
        update_param!(m, :antarctic_icesheet, :ais_ν, antarctic_nu)
        update_param!(m, :antarctic_icesheet, :ais_iceflow₀, antarctic_flow0)
        update_param!(m, :antarctic_icesheet, :ais_γ, antarctic_gamma)
        update_param!(m, :antarctic_icesheet, :ais_α, antarctic_alpha)
        update_param!(m, :antarctic_icesheet, :temperature_threshold, antarctic_temp_threshold)
        update_param!(m, :antarctic_icesheet, :λ, antarctic_lambda)

        # ----- Glaciers & Small Ice Caps ----- #
        update_param!(m, :glaciers_small_icecaps, :gsic_β₀, glaciers_beta0)
        update_param!(m, :glaciers_small_icecaps, :gsic_v₀, glaciers_v0)
        update_param!(m, :glaciers_small_icecaps, :gsic_s₀, glaciers_s0)
        update_param!(m, :glaciers_small_icecaps, :gsic_n, glaciers_n)

        # ----- Greenland Ice Sheet ----- #
        update_param!(m, :greenland_icesheet, :greenland_a, greenland_a)
        update_param!(m, :greenland_icesheet, :greenland_b, greenland_b)
        update_param!(m, :greenland_icesheet, :greenland_α, greenland_alpha)
        update_param!(m, :greenland_icesheet, :greenland_β, greenland_beta)
        update_param!(m, :greenland_icesheet, :greenland_v₀, greenland_v0)

        # ----- Thermal Expansion ----- #
        update_param!(m, :thermal_expansion, :te_α, thermal_alpha)
        update_param!(m, :thermal_expansion, :te_s₀, thermal_s0)

        # ----- SNEASY/DOECLIM Parameters ----- #
        update_param!(m, :doeclim, :t2co, climate_sensitivity)
        update_param!(m, :doeclim, :kappa, heat_diffusivity)
        update_param!(m, :radiativeforcing, :alpha, rf_scale_aerosol)
        update_param!(m, :model_CO₂_0, CO2_0)
        update_param!(m, :ccm, :Q10, Q10)
        update_param!(m, :ccm, :Beta, CO2_fertilization)
        update_param!(m, :ccm, :Eta, CO2_diffusivity)
        update_param!(m, :rfco2, :N₂O_0, N2O_0)

        # run SNEASY-BRICK
        run(m)

        # --------------------------------------------------- Retrieve and store output for current run -------------------------------------------------------- #

        # parameter values to be saved 
        param_vals = [CO2_0, N2O_0, thermal_s0, greenland_v0, glaciers_v0, glaciers_s0, antarctic_s0, Q10, CO2_fertilization, CO2_diffusivity, heat_diffusivity, rf_scale_aerosol, climate_sensitivity, 
                    thermal_alpha, greenland_a, greenland_b, greenland_alpha, greenland_beta, glaciers_beta0, glaciers_n, anto_alpha, anto_beta, antarctic_gamma, antarctic_alpha, antarctic_mu, antarctic_nu, 
                    antarctic_precip0, antarctic_kappa, antarctic_flow0, antarctic_runoff_height0, antarctic_c, antarctic_bed_height0, antarctic_slope, antarctic_lambda, antarctic_temp_threshold, lw_random_sample]

        # write current sample to respective array
        parameters[i,:]                     = param_vals                                                # values for each parameter
        co2_emissions[i,:]                  = m[:ccm, :CO2_emissions] .* 3.67                           # total CO₂ emissions (GtCO₂/yr)
        co2_concentrations[i, :]            = m[:ccm, :atmco2]
        co2_concentrations_no_noise[i, :]   = m[:ccm, :atmco2]
        oceanco2[i, :]                      = m[:ccm, :atm_oc_flux]
        radiative_forcing[i,:]              = m[:radiativeforcing, :rf]                                 # global radiative forcing (top of atmosphere) (W/m^2)
        temperature[i,:]                    = m[:ccm, :temp]                                            # global mean temperature anomaly (K), relative to preindustrial
        global_mean_sea_level_rise[i,:]     = m[:global_sea_level, :sea_level_rise]                     # total sea level rise from all components (m)
        slr_antarctic_icesheet[i,:]         = m[:global_sea_level, :slr_antartic_icesheet]              # sea level rise from the Antarctic ice sheet (m)
        slr_glaciers_small_ice_caps[i,:]    = m[:global_sea_level, :slr_glaciers_small_ice_caps]        # sea level rise from glaciers and small ice caps (m)
        slr_greenland_icesheet[i,:]         = m[:global_sea_level, :slr_greeland_icesheet]              # sea level rise from the Greenland ice sheet (m)
        slr_landwater_storage[i,:]          = m[:global_sea_level, :slr_landwater_storage]              # sea level rise from landwater storage (m)
        slr_thermal_expansion[i,:]          = m[:global_sea_level, :slr_thermal_expansion]              # sea level rise from thermal expansion (m)
        ocean_heat[i,:]                     = m[:doeclim, :heat_mixed] .+ m[:doeclim, :heat_interior]   # sum of ocean heat content anomaly in mixed layer and interior ocean (10²² J)

        # -------------------------------------------------- Incorporate statistical noise for current run ----------------------------------------------------- #
        
        # calculate the statistical noise
        noise_temperature   = MimiBRICK.simulate_ar1_noise(num_years, σ_temperature, ρ_temperature, zeros(end_year - start_year + 1))
        noise_ocean_heat    = MimiBRICK.simulate_ar1_noise(num_years, σ_ocean_heat, ρ_ocean_heat, zeros(end_year - start_year + 1))
        noise_glaciers      = MimiBRICK.simulate_ar1_noise(num_years, σ_glaciers, ρ_glaciers, zeros(end_year - start_year + 1))
        noise_greenland     = MimiBRICK.simulate_ar1_noise(num_years, σ_greenland, ρ_greenland, zeros(end_year - start_year + 1))
        noise_antarctic     = MimiBRICK.simulate_ar1_noise(num_years, σ_antarctic, ρ_antarctic, zeros(end_year - start_year + 1))
        noise_gmsl          = MimiBRICK.simulate_ar1_noise(num_years, σ_gmsl, ρ_gmsl, zeros(end_year - start_year + 1))
        normal_noise_oceanco2 = MimiBRICK.rand(Normal(0,0.4*sqrt(10)), num_years)
        car1_noise_co2 = MimiBRICK.simulate_car1_noise(num_years, α₀_CO₂, σ_co2, zeros(end_year - start_year + 1))

        # define baseline indices to normalize results
        temperature_norm_indices = findall((in)(1850:1900), start_year:end_year) # indices needed to normalize temperature anomalies relative to 1861-1880 mean
        sealevel_norm_indices_1995_2014 = findall((in)(1995:2014), start_year:end_year) # indices needed to normalize sea level rise sources relative to the 1992-2001 mean
        
        # subtract off the mean so we have results relative to a baseline
        @view(slr_glaciers_small_ice_caps[i,:]) .-= mean(slr_glaciers_small_ice_caps[i,:][sealevel_norm_indices_1995_2014]) .+ noise_glaciers
        @view(slr_greenland_icesheet[i,:])      .-= mean(slr_greenland_icesheet[i,:][sealevel_norm_indices_1995_2014]) .+ noise_greenland
        @view(slr_antarctic_icesheet[i,:])      .-= mean(slr_antarctic_icesheet[i,:][sealevel_norm_indices_1995_2014]) .+ noise_antarctic
        @view(global_mean_sea_level_rise[i,:])  .-= mean(global_mean_sea_level_rise[i,:][sealevel_norm_indices_1995_2014]) .+ noise_gmsl
        @view(temperature[i,:])                 .-= mean(temperature[i,:][temperature_norm_indices]) .+ noise_temperature .+ temperature_0
        @view(co2_concentrations[i, :]) .+= car1_noise_co2 
        @view(ocean_heat[i,:])                  .+= noise_ocean_heat .+ ocean_heat_0
        @view(oceanco2[i,:])                  .+= normal_noise_oceanco2
        
    end

    # transform matrices to dataframes to improve interpretability (add parameter names/years)
    co2_emissions_df                = DataFrame(co2_emissions, Symbol.([model_years...]))
    co2_concentrations_df           = DataFrame(co2_concentrations, Symbol.([model_years...]))
    co2_concentrations_no_noise_df  = DataFrame(co2_concentrations_no_noise, Symbol.([model_years...]))
    oceanco2_df                     = DataFrame(oceanco2, Symbol.([model_years...]))
    radiative_forcing_df            = DataFrame(radiative_forcing, Symbol.([model_years...]))
    temperature_df                  = DataFrame(temperature, Symbol.([model_years...]))
    global_mean_sea_level_rise_df   = DataFrame(global_mean_sea_level_rise, Symbol.([model_years...]))
    slr_antarctic_icesheet_df       = DataFrame(slr_antarctic_icesheet, Symbol.([model_years...]))
    slr_glaciers_small_ice_caps_df  = DataFrame(slr_glaciers_small_ice_caps, Symbol.([model_years...]))
    slr_greenland_icesheet_df       = DataFrame(slr_greenland_icesheet, Symbol.([model_years...]))
    slr_landwater_storage_df        = DataFrame(slr_landwater_storage, Symbol.([model_years...]))
    slr_thermal_expansion_df        = DataFrame(slr_thermal_expansion, Symbol.([model_years...]))
    ocean_heat_df                   = DataFrame(ocean_heat, Symbol.([model_years...]))

    # export output to .csv files
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "parameters.csv"), A)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "emissions.csv"), co2_emissions_df)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "concentrations.csv"), co2_concentrations_df)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "concentrations_nonoise.csv"), co2_concentrations_no_noise_df)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "oceanco2.csv"), oceanco2_df)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "radiative_forcing.csv"), radiative_forcing_df)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "temperature.csv"), temperature_df)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "gmslr.csv"), global_mean_sea_level_rise_df)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "antarctic.csv"), slr_antarctic_icesheet_df)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "gsic.csv"), slr_glaciers_small_ice_caps_df)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "greenland.csv"), slr_greenland_icesheet_df)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "lw_storage.csv"), slr_landwater_storage_df)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "thermal_expansion.csv"), slr_thermal_expansion_df)
    save(joinpath(@__DIR__, "..", "results", "ssp", output_dir, "ocean_heat.csv"), ocean_heat_df)

    return nothing
end    

for (i, ssp) in pairs(ssps)
    start_year = 1850
    end_year = 2100
    start_idx = findfirst(names(emis_dat) .== string(start_year))
    end_idx = findfirst(names(emis_dat) .== string(end_year))
    run_model(params, collect(emis_dat[i, start_idx:end_idx]), collect(n2o_conc[i, start_idx:end_idx]), collect(aerosol_forc[i, start_idx:end_idx]), collect(other_forc[i, start_idx:end_idx]) - (collect(aerosol_forc[i, start_idx:end_idx]) + collect(co2_forc[i, start_idx:end_idx])), start_year, end_year, ssp)
end
