# Theoretical_Analysis

Repo for the R code in "Quantifying the efficacy of taxation on antibiotics to control antimicrobial resistance in food animals using mathematical modelling". 

Code is seperated based on the different sections of the study. 1) Baseline Trajectory Plots (), 2) Uncertainty Analysis Model Code (), 3) ABC-SMC Model Fit () and 4) Analysis of the Distribution of Uncertainty Analysis Runs (). There is also a number of R files related to extra analysis conducted. 

Baseline Trajectory Plots 

Seperate trajectory plots were created for the Baseline Parameter Set, Model with Two Antibiotics and Model with Four Antibiotics. 

Uncertainty Analysis 

The .R file for the baseline uncertainty analysis can be found in the Final_Runs folder. This folder also includes model runs for a number of different scenario analyses. These include: 1) Four antibiotics (), 2) two antibiotics (), 3) 25% effectiveness threshold, 4) 75% effectiveness threshold, 5) Biased PED matrix and 6) Realistic PED matrix. An additional analysis was also included for the last section of the result - involving an uncertainty analysis comparing bans and optimal strategies for taxation (). 

The resulting parameter combinations and absolute values from the uncertainty analyses,s tored as .RDS files can be found in the Model_Output folder. 

Analysis of Uncertainty Runs 

An .R file is also needed to present the information from the uncertainty analyses comparin ghte proportion of times each intervention "wins" and the distribution of absolute values across each intervention for each optimisation criteria. This can be found in the Distribution_Analyses folder. 

This figure was created for each of the different scenario analyses through importation of .RDS files into the main CombHeatMap_Dist.R file. Seperate .R files wwere used for the two and four antibiotic class examples and the comparison between bans and taxation. 

A seperate .R file was used to aggregate the information from the different uncertainty analysis into a single heatmap (). 

Model Fit 

The baseline model was fitted to ensure heterogeneity across the different antibiotic classes. The ABC-SMC files can be found in the Model_Fit folder and can be seperated into files for the analysis of the poserior dsitribution and running the algorithm. The Model_Output folder contains the 10 generations obtained from the ABC-SMC model run. 

