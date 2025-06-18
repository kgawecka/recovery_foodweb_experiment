# Ecological recovery: the roles of space and food web complexity

This repository contains the data and code used in the manuscript [Gawecka et al. (2025) "The roles of space and food web complexity in mediating ecological recovery"](https://www.biorxiv.org/content/10.1101/2025.05.13.653715v1).

## `Code`
All code was created in R version 4.4.0.

`experiment_analysis.R` - analysis of the main experiment, including manuscript figures

`parameterization.R` - analysis of parameterization experiments and determination of model parameters

`model_simulations.R` - metacommunity model functions, model simulations and model output analysis, including manuscript figures

## `Data`

`experiment_data.xlsx` - main experiment data

`parameterization_data.xlsx` - parameterization experiments data

`model_parameters.csv` - metacommunity model parameters determined in parameterization procedure

## `Output`

`data_recovery_metapop_exp.csv` - metapopulation recovery credit computed on experimental data (output of `experiment_analysis.R`)

`data_recovery_pop_exp.csv` - population recovery credit computed on experimental data (output of `experiment_analysis.R`)

`M_land_50.csv` - simulated landscape adjacency matrix (output of `model_simulations.R`)
