# Ecological recovery: the roles of space and food web complexity

This repository contains the data and code used in the manuscript [Gawecka, K.A, Barbour, M.A., Bullock, J.M. Bascompte, J. (2025) "The roles of space and food web complexity in mediating ecological recovery" Ecology Letters](https://www.biorxiv.org/content/10.1101/2025.05.13.653715v1).

All data and code have been archived on [Zenodo](https://zenodo.org/records/17242669), and are available under a CC-BY-4 license.

## `Code`
All code was created in R version 4.4.0.

`experiment_analysis.R` - analysis of the main experiment, including manuscript figures

`parameterisation.R` - analysis of parameterization experiments and determination of model parameters

`model_simulations.R` - metacommunity model functions, model simulations and model output analysis, including manuscript figures

## `Data`

`experiment_data.csv` - main experiment data

`parameterisation_data_X.csv` - parameterisation experiments data
- `_BRBR` - single aphid species (_Brevicoryne brassicae_) response
- `_LIER` - single aphid species (_Lipaphis erysimi_) response
- `_BRBR_LIER` - two aphid species (_Brevicoryne brassicae_ and _Lipaphis erysimi_) response
- `_DIRA_dispersal` - parasitoid (_Diaeretiella rapae_) dispersal
- `_DIRA_emergence` - parasitoid (_Diaeretiella rapae_) emergence from mummies
- `_DIRA_females` - parasitoid (_Diaeretiella rapae_) proportion of females
- `_DIRA_function` - parasitoid (_Diaeretiella rapae_) functional response
- `_DIRA_mortality` - parasitoid (_Diaeretiella rapae_) mortality
- `_DIRA_parasitism` - parasitoid (_Diaeretiella rapae_) maximum parasitism rate

`model_parameters.csv` - metacommunity model parameters determined in parameterisation procedure

## `Output`

`data_recovery_metapop_exp.csv` - metapopulation recovery credit computed on experimental data (output of `experiment_analysis.R`)

`data_recovery_pop_exp.csv` - population recovery credit computed on experimental data (output of `experiment_analysis.R`)

`M_land_50.csv` - simulated landscape adjacency matrix (output of `model_simulations.R`)

`out_X_Y.csv` - output of metacommunity model simulations as time series of species abundance (output of `model_simulations.R`); X - simulated community; Y - number of patches in simulated landscape
