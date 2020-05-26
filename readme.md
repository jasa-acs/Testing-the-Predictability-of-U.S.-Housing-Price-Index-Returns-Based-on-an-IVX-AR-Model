# Testing the Predictability of U.S. Housing Price Index Returns Based on an IVX-AR Model

# Author Contributions Checklist Form

## Data

### Abstract

From the Federal Housing Finance Agency (FHFA), we collect the quarterly House Price Index (HPI, Index 1980:Q1=100) during 1975:Q1 -- 2018:Q2. This index is based on all-transactions of properties and calculated through a modified version of the weighted-repeat sales by Case and Shiller (1989). Ten macroeconomic variables are collected from the Federal Reserve Economic Data (FRED), and all data are quarterly between 1975:Q1 and 2018:Q2.




### Availability 

These data are publicly available on the Federal Reserve Bank of St. Louis (https://fred.stlouisfed.org/) and the Federal Housing Finance Agency (https://www.fhfa.gov/DataTools/Downloads/Pages/House-Price-Index.aspx).


### Description 

Civillian Unemployment Rate, Seasonally Adjusted (UNRATE)
Real Disposable Personal Income, Percent Change from Quarter One Year Ago, Seasonally Adjusted (A067RO1Q156NBEA)
Gross Domestic Product, Percent Change from Preceding Period, Seasonally Adjusted Annual Rate (A191RP1Q027SBEA)
Effective Federal Funds Rate, Percent, Not Seasonally Adjusted (FEDFUNDS)	
Shares of gross domestic product: Gross private domestic investment: Fixed investment: Residential, Percent, Not Seasonally Adjusted (A011RE1Q156NBEA)
Industrial Production Index, Index 2012=100, Seasonally Adjusted (INDPRO)
Consumer Price Index for All Urban Consumers: All items less shelter, Index 1982-1984=100, Seasonally Adjusted (CUSR0000SA0L2)
All-Transactions House Price Index for the United States, Index 1980:Q1=100, Not Seasonally Adjusted (USSTHPI)
30-Year Fixed Rate Mortgage Average in the United States, Percent, Not Seasonally Adjusted (MORTGAGE30US)
Gross Domestic Product: Implicit Price Deflator, Index 2012=100, Seasonally Adjusted (GDPDEF)
Total Reserve Balances Maintained with Federal Reserve Banks, Billions of Dollars, Not Seasonally Adjusted (RESBALNS)



## Code

### Abstract

The file, Simulation_README.txt, shows how to replicate our simulation studies, and the file, Empirical_README.txt, shows how to replicate the results in the empirical application.


### Description

1. The "Simulation" folder contains two subdirectories: "Univariate" and "Multivariate".
* The folder named "Univariate" contains R codes for the univariate model. 
- The R file named "uni_function.R" contains the R codes used in simulations for the univariate cases (IVX-AR by AIC and BIC, IVX-KMS, data generating processes);
- The R file named "univariate_size_tables.R" replicates Tables 1-2 in the main text;
- The R file named "univariate_power_plots.R" replicates Figures 1-3 in the main text.

* The folder named "Multivariate" contains R codes for the multivariate model. 
- The R file named "multi_function.R" contains the R codes used in simulations for the multivariate cases (IVX-AR by AIC and BIC, IVX-KMS, data generating processes);
- The R file named "multivariate_size_tables.R" replicates Table 3 in the main text;
- The R file named "multivariate_power_plots.R" replicates Figure 4 in the main text.

2. This "Empirical" folder contains two types of files:
* Eleven Excel spreadsheets (.csv) which contains the data used in the empirical application;
* "FUN_univariate.R" contains the functions used for the univariate predictive model;
* "FUN_multivariate.R" contains the functions used for the multivariate predictive model;
* "empirical_fullsample_replicate.R" contains codes to replicate Figures 5-8, Tables 4-6 and Panel 1 in Table 7 and 8;
* "empirical_subsample_replicate.R" contains codes to replicate Panel 2 in Table 7 and 8.


###Optional information 

R version: 3.5.0 (2018-04-23) – “Joy in Playing”
Packages: copula_0.999-18, ggplot2_2.2.1, lmtest_0.9-36, zoo_1.8-1, urca_1.3-0, forecast_8.4, tseries_0.10-45, doParallel_1.0.15, MASS_7.3-51.4


## Instructions for Use


### Reproducibility 

The simulation results are obtained by parallel calculation. 8 cores are registered, and the simulations need about 60 hours to complete.

The R codes, datasets and readme files can be downloaded from the public Google Drive folder:
https://drive.google.com/open?id=14XDjxhX00RInyqRp95gl3kzuy51NJl6x
