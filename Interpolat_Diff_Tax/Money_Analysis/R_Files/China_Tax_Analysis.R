library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
rm(list=ls())

# Import in Resistance Dataset -------------------------------------------------------
res_data <- read.csv("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Money_Analysis/China_data.csv")

#Salmonella
salm_res_data <- res_data[res_data$Pathogens == "Salmonella" & res_data$PubDate > 2015 & res_data$SampleType == "Meat",]
table(salm_res_data$Class)

data_salm <- salm_res_data[salm_res_data$Class == "Polymyxins" & salm_res_data$Species == "Pig",]
mean(data_salm$Rescom, na.rm = T)

#Campylobacter
campy_res_data <- res_data[res_data$Pathogens == "Campylobacter" & res_data$PubDate > 2015 & res_data$SampleType == "Meat",]
table(campy_res_data$Class)

data_campy <- campy_res_data[campy_res_data$Class == "Lincosamides" & campy_res_data$Species == "Chicken",]
mean(data_campy$Rescom, na.rm = T)

# Import in Antibiotic Sales Data -----------------------------------------

sales_data <- read.csv("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Money_Analysis/Antimicrobial_Sales_data_2022_v6_tvb.csv")
china_sales_data <- sales_data[sales_data$CountryName == "China" & sales_data$Year == 2020,]
