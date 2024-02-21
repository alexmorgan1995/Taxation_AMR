# Money Analysis


## Data 

This folder includes the data to perform the revenue assessment and the relevant .R files to conduct the Monte Carlo sampling and model runs. 

Data can be found in the `Data` folder. This includes information on the price of different antibiotics from Alibaba (obtained from Dr Cheng Zhao), `Alibaba.csv`, and for the US from websites like valleyvet.com (`Prices_import.csv`). These were then combined together into a .csv (`US_China_Group_Price_v1.csv`). 

Data on antibiotic sales was also available from Mulchandani et al, 2023, which provides estimates on antibiotic sales for every country globally in 2020. This data can be found in `UsebyCountrybyCLASS_2.csv`. 

Antibiotic resistance data is also available for the US and China. This was obtained from NARMS and from PPS surveys in China. This was combined together into an excel file `chinese+USresistance.xlsx`. 

Categorisation was also available to designate countries into LMIC or HIC brackets `LMIC_class.csv`. 

## Model Analysis Files 

Two sets of model analysis .R files are available for the revenue analysis. 

This first set involves Monte Carlo based sampling of model parameter space, and generates the change in antibiotic usage over time across the different interventions. This was done for 50% (baseline) (`Tax_Analysis_Sens_v2.R'`), 25% (`Tax_Analysis_Sens_v2_25.R`) and 75% (`Tax_Analysis_Sens_v2_75.R`) antibiotic taxation rates. 

These sets of files rely on `US_China_Group_Price_v1.csv` to provide information on the pricing of different antibiotic classes. 

The second set of files involves the analysis of the distribution of effect sizes or revenue from the revenue model output. 

