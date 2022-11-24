# Theoretical_Analysis

Repo for the R code in "Quantifying the efficacy of taxation on antibiotics to control antimicrobial resistance in food animals using mathematical modelling". 

Code is seperated based on the different sections of the study. 1) Baseline Trajectory Plots (`Baseline`), 2) Uncertainty Analysis Model Code (`Euler_Run`), 3) ABC-SMC Model Fit (`Model_Fit`) and 4) Analysis of the Distribution of Uncertainty Analysis Runs (`Distribution_Analysis`). There is also a number of R files related to any extra analysis conducted (`Extra_Anal`). 

## Baseline Trajectory Plots 

Seperate trajectory plots were created for the Baseline Parameter Set, Model with Two Antibiotics and Model with Four Antibiotics. These can be found in the `Baseline` folder. 

## Uncertainty Analysis 

The .R file for the baseline uncertainty analysis can be found in the `Final_Runs` folder (`MDR_Sensitivity_interpolat_new.R`). This folder also includes model runs for a number of different scenario analyses. These include: 1) Four antibiotics (`MDR_sensitivity_four_v1.R`), 2) Two antibiotics (`MDR_Sensitivity_two_interpolat_v1.R`), 3) 25% effectiveness threshold (`MDR_Sensitivity_interpolat_25_v1.R`), 4) 75% effectiveness threshold (`MDR_Sensitivity_interpolat_75_v1.R`), 5) Biased PED matrix (`MDR_Sensitivity_interpolat_biasPED.R`) and 6) Realistic PED matrix (`MDR_Sensitivity_interpolat_realPED_v1.R`). An additional analysis was also included for the last section of the result - involving an uncertainty analysis comparing bans and optimal strategies for taxation (`MDR_Sensitivity_interpolat_v3_bans_new.R`). 

The resulting parameter combinations and absolute values from the uncertainty analyses, stored as .RDS files can be found in the `Model_Output` folder. 

## Analysis of Uncertainty Runs 

An .R file is also needed to present the information from the uncertainty analyses comparing the proportion of times each intervention "wins" and the distribution of absolute values across each intervention for each optimisation criteria. This can be found in the `Distribution_Analyses` folder. 

This figure was created for each of the different scenario analyses through importation of .RDS files into the main `CombHeatMap_Dist.R` file. Seperate .R files were used for the two and four antibiotic class examples and the comparison between bans and taxation. 

A seperate .R file was used to aggregate the information from the different uncertainty analysis into a single heatmap (`Scenario_Anal_Compare_v2.R`). 

## Model Fit 

The baseline model was fitted to ensure heterogeneity across the different antibiotic classes. The ABC-SMC files can be found in the `Model_Fit` folder and can be seperated into files for the analysis of the posterior dsitribution and running the algorithm. The `Model_Output` folder contains the .csv files for 10 generations obtained from the ABC-SMC model run. 
