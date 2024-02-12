library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
rm(list=ls())

# Import in Resistance Dataset -------------------------------------------------------
res_data <- read.csv("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Money_Analysis/Data/China_data.csv")

#Salmonella
salm_res_data <- res_data[res_data$Pathogens == "Salmonella" & res_data$PubDate > 2015 & res_data$SampleType == "Meat",]
table(salm_res_data$Class)

data_salm <- salm_res_data[salm_res_data$Class == "Polymyxins" & salm_res_data$Species == "Pig",]
mean(data_salm$Rescom, na.rm = T)

data_salm <- salm_res_data[salm_res_data$Class == grep("Ceph", salm_res_data$Class) & salm_res_data$Species == "Pig",]
mean(data_salm$Rescom, na.rm = T)

#Campylobacter
campy_res_data <- res_data[res_data$Pathogens == "Campylobacter" & res_data$PubDate > 2015 & res_data$SampleType == "Meat",]
table(campy_res_data$Class)

data_campy <- campy_res_data[campy_res_data$Class == "Lincosamides" & campy_res_data$Species == "Chicken",]
mean(data_campy$Rescom, na.rm = T)

# Cephalosporins ----------------------------------------------------------

#Salm
ceph_salm <- salm_res_data[grep("Ceph", salm_res_data$Class),]

data_ceph_salm <- ceph_salm[ceph_salm$Species == "Chicken",]
mean(data_ceph_salm$Rescom, na.rm = T)

data_ceph_salm <- ceph_salm[ceph_salm$Species == "Pig",]
mean(data_ceph_salm$Rescom, na.rm = T)

data_ceph_salm <- ceph_salm[ceph_salm$Species == "Cattle",]
mean(data_ceph_salm$Rescom, na.rm = T)

#Campy
ceph_campy <- campy_res_data[grep("Ceph", salm_res_data$Class),]

data_ceph_campy <- ceph_campy[ceph_campy$Species == "Chicken",]
mean(data_ceph_campy$Rescom, na.rm = T)

data_ceph_campy <- ceph_campy[ceph_campy$Species == "Pig",]
mean(data_ceph_campy$Rescom, na.rm = T)

data_ceph_campy <- ceph_campy[ceph_campy$Species == "Cattle",]
mean(data_ceph_campy$Rescom, na.rm = T)

# Import in Antibiotic Sales Data -----------------------------------------

sales_data <- read.csv("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Money_Analysis/Data/Antimicrobial_Sales_data_2022_v6_tvb.csv")
china_sales_data <- sales_data[sales_data$CountryName == "China" & sales_data$Year == 2020,]
