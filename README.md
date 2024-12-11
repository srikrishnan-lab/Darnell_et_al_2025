# Darnell-etal_2024_inprep

**The interplay of future emissions and geophysical uncertainties for future sea level rise**

Darnell, Chloe<sup>1</sup>, Rennels, Lisa<sup>2,3</sup>, Errickson, Frank<sup>4</sup>, Wong, Tony<sup>5</sup>, Srikrishnan, Vivek<sup>1\*</sup>

<sup>1</sup> Department of Biological & Environmental Engineering, Cornell University, Ithaca, New York, USA  
<sup>2</sup> Doerr School of Sustainability, Stanford University, Palo Alto, California, USA
<sup>3</sup> Energy and Resources Group, University of California Berkeley, Berkeley, California, USA  
<sup>4</sup> School of Public and International Affairs, Princeton University, Princeton, New Jersey, USA  
<sup>5</sup> School of Mathematics and Statistics, Rochester Institute of Technology, Rochester, New York, USA

\* corresponding author:  viveks@cornell.edu

## Abstract

Uncertainty in future carbon dioxide (CO<sub>2</sub>) emissions, and the geophysical response to emissions, drives variability in future sea-level rise (SLR). However, the relative contribution of emissions and geophysical dynamics (e.g. Antarctic Ice Sheet (AIS) tipping points) to future sea-level projections is poorly understood.  Here, we disentangle their relative importance by propagating several ensembles of CO<sub>2</sub> emissions trajectories, representing relevant deep uncertainties, through a calibrated carbon cycle-climate-sea-level model chain. The CO<sub>2</sub> emissions trajectory, particularly the timing of when emissions are reduced, becomes the primary driver of sea-level variability only after 2075. The most extreme global mean SLR (exceeding 4m by 2200) is projected to occur regardless of optimism about limiting CO<sub>2</sub> emissions if accelerated AIS melting occurs. Further, delaying decarbonization reduces the “safe operating space” associated with the geophysical uncertainties. Our results highlight the potential that similar adaptation requirements may be needed regardless of optimism about future levels of CO<sub>2</sub> mitigation.

## Acknowledgements

CD was partially funded by the College of Agricultural \& Life Sciences, Cornell University. VS was partially funded by the U.S. Department of Energy, Office of Science, Biological and Environmental Research Program, Earth and Environmental Systems Modeling, MultiSector Dynamics as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.

## Journal Reference

* **Preprint:** <https://osf.io/preprints/osf/j47ts>
 
## Code Reference

## Data Reference

### Input Data

- Historical CO<sub>2</sub> emissions data was obtained from the [Global Carbon Project](https://www.globalcarbonproject.org/) 2022.
- RCP-SSP CO<sub>2</sub> emissions and non-CO<sub>2</sub> radiative forcing projections were obtained from the [CMIP6 SSP database](https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=10).
- SNEASY-BRICK calibration file, obtained from <https://zenodo.org/records/6626335>. The scripts here are designed to use the `parameters_subsample_sneasybrick.csv` file, which is a burned-in and thinned version; some modifications might be needed if working with the full chain (`parameters_full_chain_sneasybrick.csv`).
- **For BRICK benchmarking relative to AR6**: 
  - Data for RCMIP protocol v5.1, including [concentrations](https://gitlab.com/rcmip/rcmip/-/blob/master/data/protocol/rcmip-concentrations-annual-means-v5-1-0.csv?ref_type=heads), [emissions](https://gitlab.com/rcmip/rcmip/-/blob/master/data/protocol/rcmip-emissions-annual-means-v5-1-0.csv?ref_type=heads), and [radiative forcings](https://gitlab.com/rcmip/rcmip/-/blob/master/data/protocol/rcmip-radiative-forcing-annual-means-v5-1-0.csv?ref_type=heads).
  - [IPCC AR6 sea-level projections](https://zenodo.org/records/6382554), particularly `ar6.zip`.

### Output Data

- Output ensembles for the probabilistic ensembles can be found at <https://doi.org/10.5281/zenodo.10373089>. These include CSVs of all of the individual GMSLR components along with the total GMSLR across both the main/full ensemble (`default/`), and the optimistic (`optimistic/`) and pessimistic (`pessimistic/`) ensembles.
- Output for BRICK runs forced by SSP emissions can be found at <https://doi.org/10.5281/zenodo.14346559>

## Dependencies

This code is based on Julia 1.9.4. Relevant dependencies are in the `Project.toml` and `Manifest.toml` files (the `Manifest.jl` specifies the particular versions; this file should be kept as-is for perfect reproducibility but may need to be deleted and rebuilt with `Pkg.instantiate()` for different Julia versions).

## Reproduction

### Simulation 

1. After cloning the repository, install needed packages:
    ```julia
    import Pkg
    Pkg.activate(".") # from the cloned root directory
    Pkg.instantiate()
    ```
2. To re-calibrate the emissions scenario distributions, run `julia src/emissions_update_scenarios.jl`.
3. To re-simulate the main ensemble, run `julia src/model_ensemble.jl`. There should be an argument passed to the call with a 1 for the default/baseline scenario, a 2 for the optimistic scenario, and a 3 for the pessimistic scenario. This will write output into `results/<scenario>`. This should take about 10 hours for a given scenario on a typical computer. Reducing the size of the ensemble by modifying line 37 in `src/model_ensemble.jl` will speed this up.
4. To re-run the Shapley analysis, after the main ensemble is run, run `julia src/regression_drivers.jl`. There should be an argument passed to the call with a 1 for the default/baseline scenario, a 2 for the optimistic scenario, and a 3 for the pessimistic scenario. This will write output into `output/shapley/`. This can take 36-48 hours and may run into memory issues, but the ensemble size or number of years can be reduced for a quick check.
5. To force BRICK using the SSPs, run `julia src/ssp_emissions_concentrations.jl`, which will write output into `results/ssp`. This might take 10-24 hours to run depending on the computational environment.

### Figures

1. Run the simulations above or download the results and SSP-forced output from the Zenodo repositories. The Shapley output is provided with the GitHub repository or can be re-evaluated.  
2. Run the following scripts for each of the figures and tables (none should take longer than a half-hour on a typical computer):
    
    | Figure/Table | Script | How To Run | Output File |
    | --- | --- | --- | --- |
    | Fig. 1  | `src/plot_cumemissions_slr.jl` | `julia src/plot_cumemissions_slr.jl` | `figures/slr_temps_all.png` |
    | Fig. 2 | `src/plot_ensemble.jl` | `julia src/plot_ensemble.jl` | `figures/ensemble_projections.png` |
    | Fig. 3 | `src/plot_regression_shapley_indices.jl` | `julia src/regression_shapley_indices.jl` | `figures/stacked-shapley-index-scenarios.png` |
    | Fig. 4 | `src/scenario_discovery.jl` | `julia src/scenario_discovery.jl` | `figures/factor_map_all_scenarios.png` |
    | Supp. Fig. A1 | `src/plot_sample_emissions.jl`|  `julia src/plot_sample_emissions.jl`| `figures/sample-emissions.png` |
    | Supp. Fig. A2 | `src/plot_emissions.jl` | `julia src/plot_emissions.jl` | `figures/scenario_emissions.png` |
    | Supp. Fig. A3 | `src/plot_ipcc_ar6.jl` | `julia src/plot_ipcc_ar6.jl` | `figures/brick_ipcc_compare.png` |
    | Supp. Fig. A4 | `src/plot_ais_ar6.jl` | `julia src/plot_ais_ar6.jl` | `figures/slr_ar6.png` |
    | Supp. Fig. A5 | `src/emissions_update_scenario.jl` | `julia src/emissions_update_scenario.jl` | `figures/emissions_updates.png` |
    | Supp. Fig. A6 | `src/scenario_discovery.jl` | `julia src/scenario_discovery.jl` | `figures/feature_importance_scenarios.png` |

### Tables

To produce the information in Table A1, run `ssp_temperatures.jl`.

To produce Supp. Tables A2, A3, and A4, run `src/shapley_table.jl`.

