library("deSolve"); library("reshape2"); library("parallel")
rm(list=ls())

# Read in Sales Data and Remove CIs ---------------------------------------

sales_data <- read.csv("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Money_Analysis/UsebyCountrybyCLASS_2.csv")
sales_data$Total <- as.numeric(substr(sales_data$Total, 1, nchar(sales_data$Total)-8))
sales_data_usage <- sales_data[sales_data$Total != 0,]

sales_data_usage[] <- lapply(sales_data_usage, gsub, pattern=" (+NA%)", replacement="", fixed = TRUE)
sales_data_usage[,3:16] <- sapply(sales_data_usage[, c(3:16)], as.numeric)

# Remove NA Countries (Antibiotic Class) ----------------------------------

sales_data_usage_noNA <- sales_data_usage[rowSums(is.na(sales_data_usage)) != ncol(sales_data_usage[,4:16]),]

#What percentage of the total usage did we remove (no stratified class data)
percent_eff <- sum(sales_data_usage_noNA$Total) / sum(sales_data$Total)

sales_data_usage_noNA[,c(3:16)] <- sales_data_usage_noNA[,c(3:16)]*1000
sales_data_usage_noNA$Cephalosporins <- sales_data_usage_noNA$Cephalosporins_1_2_kg + sales_data_usage_noNA$Cephalosporins_3_4_kg
sales_data_usage_noNA <- sales_data_usage_noNA[,-c(7,8)]

# Assigning LMIC or HIC classification ------------------------------------

economy_data <- read.csv("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Money_Analysis/LMIC_class.csv")

merge_economy <- merge(sales_data_usage_noNA, economy_data, by.x = "Country", by.y = "Economy")

#Check if there are any countries missing from the merge when compared to the World Bank Rankings 

sales_data_usage_noNA[which(!sales_data_usage_noNA$Country %in%  merge_economy$Country ),]

#Using data from the World Bank for the 
merge_economy[merge_economy$Country == "Venezuela, RB", "Income.group"] = "Upper middle income"

#Collapse the Dataframe based on LMIC or HIC
merge_economy$Income.group

LIC_LMIC_UMIC <-  merge_economy[merge_economy$Income.group %in% c("Low income", "Lower middle income", "Upper middle income"),]
HIC <- merge_economy[merge_economy$Income.group %in% c("High income"),]

# Combined Dataframe ------------------------------------------------------

overall_antibiotics_sales <- data.frame(rbind(c("LIC_LMIC_UMIC", colSums(LIC_LMIC_UMIC[3:15])),
                                              c("HIC", colSums(HIC[3:15]))))
antibiotic_price <- data.frame(rbind(c(48.945, 34.2495, 68.031, 31.5825, 95.72475, 62.000025, 26.499, 55.899, 9.842893835, 33),
                                     c(1121.227778, 2502.533333, 534.1666667, 1221.28775, 19352.79206, 13592, 5318.8, 766.3333333, 292.477, 2396.78064)))
groupings <- data.frame(rbind(c(1, 1, 2, 1, 3, 2, 2, 2, 3, 3, 2),
                              c(1, 3, 2, 3, 1, 3, 3, 2, 3, 3, 2)))

colnames(antibiotic_price) <- c("Tetracyclines", "Amphenicols", "Penicillin", "Cephalosporins", "Sulphonamides", "Macrolides", "Aminoglycosides",
                                         "Quinolones", "Lincosamides", "Polymyxins", "Other")

colnames(groupings) <- c("Tetracyclines", "Amphenicols", "Penicillin", "Cephalosporins", "Sulphonamides", "Macrolides", "Aminoglycosides",
                                  "Quinolones", "Lincosamides", "Polymyxins", "Other")








HIC_LMIC_data <- data.frame()

merge_economy

all_sales_data <- sales_data[sales_data$Animal == "all", -c(26:31)]

unique(sales_data$CountryName)

for(i in 1:unique(all_sales_data$CountryName)) {
  
  
}
  