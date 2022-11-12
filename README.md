`main.R` reproduces the empirical results (fitted models, tables, plots) of 
Schnürch, Kleinow, Korn, Wagner (2022): The Impact of Mortality Shocks on Modeling and Insurance Valuation as Exemplified by COVID-19.
These results are not saved to disk but are available as R variables (data.frames and plot objects) after running the script.
Based on this, you can reproduce the results in the paper or run your own analyses.

Prerequisites:
- **Most importantly**, the following three data sources must be downloaded before being able to run the script (see also `settings.R`):
  - HMD: RData file containing a list with HMD data from mortality.org. The list should be called dataHMD and should have the entries mx, Dx, Ex containing death rates, deaths and exposures. This RData file can be generated with the helper function downloadData in /functions_data.R (HMD account required).
  - STMF: Folder containing STMF **input** files from mortality.org (https://mortality.org/File/GetDocument/Public/STMF/Inputs/STMFinput.zip). The folder must contain one csv file per country with at least the columns PopCode ("BEL", "FRATNP", ...), Year (integer, e.g. 2020), Week (1, 2, ..., 53), Sex ("m", "f", "b"), Age (0, 1, 5, 10, ...), Deaths (integer >= 0, e.g. 143).
  - European Standard Population: csv file with columns AgeGroup ("a-b years") and EuropeanStandardPopulation (integer, e.g. 5500), e.g. from https://www.opendata.nhs.scot/dataset/standard-populations/resource/edee9731-daf7-4e0d-b525-e4c1469b8f69.
- Data paths must be provided in `settings.R`, where you can also adjust further global variables of the analyses.
- Working directory must be source file location.
