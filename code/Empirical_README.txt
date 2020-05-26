This "Empirical" folder contains two types of files:

* There are 11 Excel spreadsheets (.csv) which contains the data used in the empirical application. Specifically, they are:
 - Civillian Unemployment Rate, Seasonally Adjusted (UNRATE)
 - Real Disposable Personal Income, Percent Change from Quarter One Year Ago, Seasonally Adjusted (A067RO1Q156NBEA)
 - Gross Domestic Product, Percent Change from Preceding Period, Seasonally Adjusted Annual Rate (A191RP1Q027SBEA)
 - Effective Federal Funds Rate, Percent, Not Seasonally Adjusted (FEDFUNDS)	
 - Shares of gross domestic product: Gross private domestic investment: Fixed investment: Residential, Percent, Not Seasonally Adjusted (A011RE1Q156NBEA)
 - Industrial Production Index, Index 2012=100, Seasonally Adjusted (INDPRO)
 - Consumer Price Index for All Urban Consumers: All items less shelter, Index 1982-1984=100, Seasonally Adjusted (CUSR0000SA0L2)
 - All-Transactions House Price Index for the United States, Index 1980:Q1=100, Not Seasonally Adjusted (USSTHPI)
 - 30-Year Fixed Rate Mortgage Average in the United States, Percent, Not Seasonally Adjusted (MORTGAGE30US)
 - Gross Domestic Product: Implicit Price Deflator, Index 2012=100, Seasonally Adjusted (GDPDEF)
 - Total Reserve Balances Maintained with Federal Reserve Banks, Billions of Dollars, Not Seasonally Adjusted (RESBALNS)

* Of the 11 datasets, "USSTHPI" includes the quarterly US housing price index from 1975 to 2018. 
  This data can be downloaded from https://www.fhfa.gov/DataTools/Downloads/pages/house-price-index.aspx
  The remained 10 datasets are publicly available from https://fred.stlouisfed.org/


* "FUN_univariate.R" contains the functions used for the univariate predictive model.
  "FUN_multivariate.R" contains the functions used for the multivariate predictive model.

* "empirical_fullsample_replicate.R" contains codes to replicate Figures 5-8, Tables 4-6 and Panel 1 in Table 7 and 8. 
  The analysis is based on the full sample (1975:Q1-2018:Q2) observations.
  "empirical_subsample_replicate.R" contains codes to replicate Panel 2 in Table 7 and 8.
  The analysis is based on the subperiod (2000:Q1-2018:Q2) observations.