######### PACKAGES #########
library(StMoMo)           # mortality models
library(MortalityLaws)    # download HMD
library(tidyverse)        # data preparation, evaluation, plotting
library(fst)              # fast serialization of data.frames
library(lubridate)        # date calculations
library(tsibble)          # for yearweek function
library(forecast)         # ARIMA
library(data.table)       # fread


######### DATA PATHS #########
# There are three data sources which this script requires to run. They must be downloaded in .RData or .csv
# format (as specified below) and put in one folder. The path of this folder relative to the working directory must be 
# specified below in DataPath. The paths of the single data sources relative to this folder must be specified in
# FileHMD, FolderSTMF, FileEuropeanStandardPopulation.
# 
#
# HMD: RData file containing a list with HMD data from mortality.org. 
#      The list should be called dataHMD and should have the entries mx, Dx, Ex containing death rates, deaths and exposures.
#      This RData file can be generated with the helper function downloadData in /functions_data.R (HMD account required).
#
# STMF: Folder containing STMF input files from mortality.org (https://www.mortality.org/Public/STMF/Inputs/STMFinput.zip). 
#       The folder must contain one csv file per country with at least the columns PopCode ("BEL", "FRATNP", ...), Year (integer, e.g. 2020), 
#       Week (1, 2, ..., 53), Sex ("m", "f", "b"), Age (0, 1, 5, 10, ...), Deaths (integer >= 0, e.g. 143).
#
# European Standard Population: csv file with columns AgeGroup ("a-b years") and EuropeanStandardPopulation (integer, e.g. 5500), e.g. from
#                               https://www.opendata.nhs.scot/dataset/standard-populations/resource/edee9731-daf7-4e0d-b525-e4c1469b8f69.
#

DataPath <- "data"  # path to data folder relative to working directory

# paths to .RData or .csv files in DataPath
FileHMD <- file.path(DataPath, "2021-04-14_10.43_dataHMD5x1")       
FolderSTMF <- file.path(DataPath, "2021-04-14_STMF_input")    
FileEuropeanStandardPopulation <- file.path(DataPath, "2021-02-02_european_standard_population")   

# path for data preparation results
FileMerged <- file.path(DataPath, "2021-04-14_Covid_Data_Merged") # .fst; contains results from a previous run of the data preprocessing
                                                                  # If this file exists, no new data preprocessing is run. Instead, the old data are loaded.


######### DATA SETTINGS #########
MinYear <- 1900
MaxYear <- 2021
MinAge <- 35
MaxAge <- Inf
ExtendHistoryGermany <- "Both"         # FALSE (no data extension for Germany); "West" (extend German data with West German data); "Both" (extend German data with West+East German data)
Sexes <- c("Male", "Female", "Total")
CappingAge <- 90                       # age above which data for all older age groups are summarized to the single age group "CappingAge+", e.g., 90+
LengthAgeRanges <- 5                   # integer or NULL; desired length of all ages ranges, e.g., 5 will create age ranges 0-5, 5-10, ...
NBaseYearsExposureExtrapolation <- 5   # number of previous years to use for exposure extrapolation


######### MODEL SETTINGS #########
CBDMinAge <- 60                        # minimal age to consider for CBD model
ModelMinYear <- 1991                   # earliest year to calibrate models on
ModelEvaluationYears <- 2021:2050      # forecast years for models
AnnuityStartingAge <- 65               # initial age of annuity receiver
LifeInsuranceStartingAge <- 35         # initial age of life insurance holder
AnnuityDuration <- 30                  # maturity annuity in years
LifeInsuranceDuration <- 30            # maturity life insurance in years
InterestRate <- 0.005                  # yearly interest rate 
PredictionIntervalProbability <- 95    # target prediction interval coverage probability (as a percentage)


######### PLOT SETTINGS #########
PlotAges <- c(40, 50, 60, 70, 80, 90)  # ages to plot

PlotTextSize <- 48

theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  text = element_text(size = PlotTextSize),
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                  legend.position = "bottom",
                  legend.key.width = unit(3, "cm")))
PlotLineSize <- 3
PlotPointSize <- 8


######### GLOBAL CONSTANTS ######### 
# country and sex labels in data.frame
CountryMapping <- c(
  "AUT"          =              "Austria",
  "BEL"          =              "Belgium",
  "CHE"          =              "Switzerland",
  "DEUTNP"       =              "Germany", 
  "DEUTE"        =              "East Germany",
  "DEUTW"        =              "West Germany",
  "ESP"          =              "Spain",
  "FRATNP"       =              "France", 
  "ITA"          =              "Italy",
  "POL"          =              "Poland",
  "SWE"          =              "Sweden"
)

CountryMapping <- data.frame(rawName = names(CountryMapping), 
                             fullName = CountryMapping,
                             stringsAsFactors = FALSE)
SexMappingSTMF <- data.frame(rawName = c("m", "f", "b"), 
                             fullName = c("Male", "Female", "Total"),
                             stringsAsFactors = FALSE)

# Covid pandemic year for which we perform special evaluations
COVID_YEAR <- 2020    