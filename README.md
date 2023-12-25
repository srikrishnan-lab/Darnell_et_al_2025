# Darnell-etal_2024_inprep

**Rapid decarbonization reduces but does not eliminate the risk of extreme sea level rise due to uncertain Antarctic Ice Sheet marine instability**

Darnell, Chloe<sup>1</sup>, Rennels, Lisa<sup>2</sup>, Errickson, Frank<sup>3</sup>, Wong, Tony<sup>4</sup>, Srikrishnan, Vivek<sup>1\*</sup>

<sup>1</sup> Department of Biological & Environmental Engineering, Cornell University, Ithaca, New York, USA
<sup>2</sup> Energy and Resources Group, University of California Berkeley, Berkeley, California, USA
<sup>3</sup> School of Public and International Affairs, Princeton University, Princeton, New Jersey, USA
<sup>4</sup> School of Mathematical Sciences, Rochester Institute of Technology, Rochester, New York, USA

\* corresponding author:  viveks@cornell.edu

## Abstract

Sea-level rise (SLR) poses risks to millions of people worldwide. Although these risks are partially driven by anthropogenic CO2 emissions, the relative contribution of uncertain societal choices (e.g. decarbonization) and Earth-system dynamics (e.g. Antarctic Ice Sheet (AIS) tipping points) on future sea-level projections is poorly understood.  This is in part due to a reliance on a handful of non-probabilistic and non-representative emissions scenarios, which are difficult to integrate with geophysical uncertainties. Here, we use a chained model framework to disentangle the relative impacts of these uncertainties on projected global mean sea levels. The emissions trajectory, particularly the timing of when emissions are reduced, becomes the primary driver of sea-level variability prior to 2100. Even under relatively optimistic assumptions about future emissions trajectories, there is a moderate probability of triggering AIS marine ice instabilities, resulting in high-end SLR. We quantify how delaying decarbonization by even a few decades can reduce the “safe operating space” associated with the geophysical uncertainties. Our results highlight the potential for rapid mitigation of CO<sub>2</sub>  emissions to reduce, but not eliminate, high-end SLR. However, the window to avoid future high-end SLR is closing and requires emissions to decrease rapidly and soon.

## Acknowledgements

CD was partially funded by the College of Agricultural \& Life Sciences, Cornell University. VS was partially funded by the U.S. Department of Energy, Office of Science, Biological and Environmental Research Program, Earth and Environmental Systems Modeling, MultiSector Dynamics as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.

## Journal Reference

* **Preprint:**
 
## Code Reference

## Data Reference

### Input Data

- Historical CO<sub>2</sub> emissions data was obtained from the [Global Carbon Project](https://www.globalcarbonproject.org/) 2022.
- RCP-SSP CO<sub>2</sub> emissions and non-CO<sub>2</sub> radiative forcing projections were obtained from the [CMIP6 SSP database](https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=10).
- SNEASY-BRICK calibration file, obtained from <https://zenodo.org/records/6626335>. The scripts here are designed to use the `parameters_subsample_sneasybrick.csv` file, which is a burned-in and thinned version; some modifications might be needed if working with the full chain (`parameters_full_chain_sneasybrick.csv`).

### Output Data

- Output ensembles can be found at <https://zenodo.org/records/10373090>. These include CSVs of all of the individual GMSLR components along with the total GMSLR across both the main/full experiment(`default/`) and the peaking uncertainty experiment (`peaking/`).

## Dependencies

This code is based on Julia 1.9.4. Relevant dependencies are in the `Project.toml` and `Manifest.toml` files.

## Reproduction

### Simulation 

1. After cloning the repository, install needed packages:
    ```julia
    import Pkg
    Pkg.activate(".") # from the cloned root directory
    Pkg.instantiate()
    ```
2. To re-simulate the main ensemble, run `julia src/model_ensemble.jl`. This will write output into `results/default/`.
3. To re-simulate the peaking ensemble, run `julia src/peaking_ensemble.jl`. This will write output into `results/peaking/`.
4. To re-run the Shapley analysis, after the main ensemble is run, run `julia src/regression_drivers.jl`.This will write output into `output/shapley/`.

### Figures

1. Run the simulations above or download the default and peaking results from the Zenodo repository. The Shapley output is provided with the GitHub repository or can be re-evaluated.  
2. Run the following scripts for each of the figures and tables:
    
    | Figure/Table | Script | How To Run |
    | --- | --- | --- | 
    | Fig. 1  | `plot_slr_temperatures.jl` | `julia plot_slr_temperatures.jl` |
    | Fig. 2 | `plot_regression_shapley_indices.jl` | `julia plot_regression_shapley_indices.jl` |
    | Fig. 3 | `plot_peaking_uncertainty.jl` | `julia plot_peaking_uncertainty.jl` |
    | Fig. 4 | `scenario_discovery.jl` | `julia scenario_discovery.jl` |
    | Supp. Fig. A1 | `plot_ensemble.jl`|  `julia plot_ensemble.jl`| 
    | Supp. Fig. A2 | `plot_peaking_uncertainty.jl` | `julia plot_peaking_uncertainty.jl` |
    | Supp. Fig. A3 | `plot_peaking_uncertainty_all.jl` | `julia plot_peaking_uncertainty_all.jl` |
    | Supp. Fig. A4 | `milestone_timing.jl` | `julia milestone_timing.jl` |
    | Supp. Fig. A5 | `scenario_discovery.jl` | `julia scenario_discovery.jl` |
