# Script to produce empirical results (fitted models, tables, plots) of 
# Schnürch, Kleinow, Korn, Wagner (2021): The Impact of Mortality Shocks on Modeling and Insurance Valuation as Exemplified by COVID-19.
# These results are not saved anywhere but are available as R variables (data.frames and plot objects) after running the script.
# Based on this, you can reproduce the results in the paper or run your own analyses.
#
# Working directory must be source file location.
# Data paths must be provided in ./settings.R.

# source settings and functions
source("settings.R")
source("functions_data.R")
source("functions_models.R")
source("functions_comparisons.R")


#############################
### LOAD AND PREPARE DATA ###
#############################
# European standard population 
dfStandardPopulation <- prepareDataESP(paste0(FileEuropeanStandardPopulation, ".csv"),
                                       cappingAge = CappingAge,
                                       lengthAgeRanges = LengthAgeRanges)

# HMD
if (!file.exists(paste0(FileHMD, "_prepared"))) {
  dfHMD <- prepareDataHMD(path = paste0(FileHMD, ".RData"), 
                          dfStandardPopulation = dfStandardPopulation,
                          countryMapping = CountryMapping,
                          minYear = MinYear, maxYear = MaxYear,
                          minAge = MinAge, maxAge = MaxAge,
                          sexes = Sexes,
                          cappingAge = CappingAge,
                          extendHistoryGermany = ExtendHistoryGermany)
  
  write_fst(dfHMD, path = paste0(FileHMD, "_prepared"))
} else {
  dfHMD <- read_fst(paste0(FileHMD, "_prepared"))
}

# STMF 
if (!file.exists(paste0(FolderSTMF, "_prepared"))) {
  dfSTMFInput <- prepareDataSTMFInput(folder = FolderSTMF,
                                      countryMapping = CountryMapping,    #  country codes are identical to HMD
                                      sexMapping = SexMappingSTMF,
                                      minYear = MinYear, maxYear = MaxYear, 
                                      minAge = MinAge, maxAge = MaxAge,  
                                      sexes = Sexes, 
                                      cappingAge = CappingAge,
                                      lengthAgeRanges = LengthAgeRanges)
  
  write_fst(dfSTMFInput, path = paste0(FolderSTMF, "_prepared"))
  
  # calculate excess death rates based on 4-year-average 
  dfSTMFInputExcess <- dfSTMFInput %>% 
    filter(Year == COVID_YEAR) %>%
    left_join(dfSTMFInput %>%
                filter(Year %in% (COVID_YEAR - 1:4)) %>%
                group_by(Population, Age, Week) %>%
                summarize(BaseAverage = mean(DeathsSTMFInput)),
              by = c("Population", "Age", "Week")) %>%
    group_by(Population, Age, Week) %>%
    mutate(ExcessOverAverage = DeathsSTMFInput / BaseAverage - 1)
  
  write_fst(dfSTMFInputExcess, path = paste0(FolderSTMF, "_excess"))
} else {
  dfSTMFInput <- read_fst(paste0(FolderSTMF, "_prepared"))
  dfSTMFInputExcess <- read_fst(paste0(FolderSTMF, "_excess"))
}

# dfTotal (merged)
if (!file.exists(FileMerged)) {
  dfTotal <- createDfTotal(dfHMD, dfSTMFInput, NBaseYearsExposureExtrapolation)
  write_fst(dfTotal, path = FileMerged)
} else {
  dfTotal <- read_fst(FileMerged)
}


####################################
### EXPLORATIVE PLOTS/STATISTICS ###
####################################
# plot: weekly ASDRs 2019 and 2020
dfTotalWeeklyASDR <- createDfWithAgeStandardizedRate(dfTotal %>% 
                                                       filter(Year %in% 2019:2020) %>%
                                                       select(Year, Week, YearWeek, Age, Population, Country, Sex, DeathRateSTMF) %>%
                                                       unique(),
                                                     dfStandardPopulation,  
                                                     colNamesGrouping = c("Sex", "Country", "Year", "Week"),  
                                                     colNameRate = "DeathRateSTMF",
                                                     colNameASR = "AgeStandardizedRateSTMF") 

plotTotalWeeklyASDRNoFacet <- ggplot(dfTotalWeeklyASDR %>% filter(Age == min(Age), Sex == "Total"), 
                                     aes(x = YearWeek, y = AgeStandardizedRateSTMF, color = Country)) +
  geom_line(size = PlotLineSize) +
  labs(x = "week", y = "age-standardized death rate", color = "country") +
  scale_x_yearweek(date_breaks = "12 weeks", date_labels = "%Y-%U") 

# plot: weekly excess death ratio
plotSTMFExcessDeathRate <- ggplot(dfSTMFInputExcess %>% filter(Sex == "Total", Week %in% 1:52, Age %in% PlotAges),
                                  aes_string(x = "Week", y = "ExcessOverAverage")) +
  geom_line(size = PlotLineSize) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = PlotLineSize * 0.75, color = "blue") +
  facet_grid(Country ~ Age) +
  labs(x = "week", y = "excess death ratio")

# plot: yearly ASDRs over last 15 years (2006 to 2020)
dfTotalYearlyASDR <- createDfWithAgeStandardizedRate(dfTotal %>% 
                                                       filter(Week == 1 | Year < 2019) %>%
                                                       group_by(Population, Country, Sex, Age, Year) %>%
                                                       summarize(DeathRate = unique(ifelse(is.na(DeathRateHMD), DeathRateSTMFYear, DeathRateHMD))),
                                                     dfStandardPopulation,
                                                     colNamesGrouping = c("Sex", "Country", "Year"),
                                                     colNameRate = "DeathRate",
                                                     colNameASR = "AgeStandardizedRate") %>%
  filter(Age == min(Age)) %>%
  mutate(Age = NULL)

plotTotalASDR_GER_POL_ESP <- ggplot(dfTotalYearlyASDR %>% filter(Year %in% 2006:2020, Country %in% c("Germany", "Poland", "Spain")), 
                                    aes(x = Year, y = AgeStandardizedRate, color = Sex, linetype = Sex)) +
  geom_line(size = PlotLineSize) +
  facet_wrap(~ Country, ncol = 3, scales = "fixed") +
  labs(x = "year", y = "age-standardized death rate", color = "sex", linetype = "sex") +
  scale_x_continuous(breaks = seq(0, 2100, by = 2)) +
  theme(text = element_text(size = 35))

# how many years of mortality development are lost due to Covid?
dfYearsLostDev <- dfTotalYearlyASDR %>% 
  filter(Year <= COVID_YEAR) %>% 
  arrange(Year) %>% 
  group_by(Population) %>% 
  summarize(jumpBackYear = max(Year[AgeStandardizedRate > tail(AgeStandardizedRate, 1)]), # latest year where ASDR was larger than in COVID year
            nYearsLostDev = COVID_YEAR - jumpBackYear)

minDevYears <- min(dfYearsLostDev$nYearsLostDev)
maxDevYears <- max(dfYearsLostDev$nYearsLostDev)
print(paste0("Least years (", minDevYears, ") have been lost by ", 
             paste(dfYearsLostDev %>% filter(nYearsLostDev == minDevYears) %>% pull(Population), collapse = ", ")))
print(paste0("Most years (", maxDevYears, ") have been lost by ", 
             paste(dfYearsLostDev %>% filter(nYearsLostDev == maxDevYears) %>% pull(Population), collapse = ", ")))

# table: top 10 worst years (improvement rates) by country
dfTotalYearlyASIR <- createDfWithAgeStandardizedRate(dfTotal %>% 
                                                       filter(Week == 1 | Year < COVID_YEAR - 1) %>%
                                                       select(Population, Country, Sex, Year, Age, ImprovementRate) %>%
                                                       unique(),
                                                     dfStandardPopulation,
                                                     colNamesGrouping = c("Population", "Year"),
                                                     colNameRate = "ImprovementRate",
                                                     colNameASR = "AgeStandardizedImprovementRate") 

dfWorstASIR <- dfTotalYearlyASIR %>% 
  filter(Age == min(Age), Sex == "Total") %>% 
  group_by(Country) %>%
  top_n(10, -AgeStandardizedImprovementRate) %>%
  arrange(Country, AgeStandardizedImprovementRate)


########################
### MORTALITY MODELS ###
########################
dfModelData <- dfTotal %>%
  select(Population, Age, Year, DeathRate, Deaths, Exposure = ExposureHMD) %>% 
  filter(!is.na(DeathRate), !is.na(Deaths), !is.na(Exposure)) %>%
  filter(Year <= COVID_YEAR) %>%
  unique() %>%
  mutate(Population = as.character(Population)) %>%
  arrange(Population, Year, Age)
mx <- dfModelData %>% select(Population, Age, Year, DeathRate)
Dx <- dfModelData %>% select(Population, Age, Year, Deaths)
Ex <- dfModelData %>% select(Population, Age, Year, Exposure)

for (modelName in c("LC", "CBD")) {
  if (modelName == "LC") {
    fitMethods <- c("SVD", "MLE")
  } else {
    fitMethods <- c("lm", "MLE")
  }
  
  for (fitMethod in fitMethods) {  
    if (modelName == "LC") {
      mxModel <- mx
      DxModel <- Dx
      ExModel <- Ex
    } else if (modelName == "CBD") {
      mxModel <- mx %>% filter(Age >= CBDMinAge)
      DxModel <- Dx %>% filter(Age >= CBDMinAge)
      ExModel <- Ex %>% filter(Age >= CBDMinAge)
    }
    
    # 1. model with real Covid data
    modelWithCovid <- fitIndividualModel(mx = mxModel, 
                                         Ex = ExModel,
                                         Dx = DxModel, 
                                         modelName = modelName,
                                         fitMethod = fitMethod,
                                         modelYearRange = ModelMinYear:max(mxModel$Year))
    modelWithCovid$name <- paste("real", COVID_YEAR)
    
    # 2. model before Covid data
    modelBeforeCovid <- fitIndividualModel(mx = mxModel %>% filter(Year < COVID_YEAR),
                                           Ex = ExModel %>% filter(Year < COVID_YEAR),
                                           Dx = DxModel %>% filter(Year < COVID_YEAR),
                                           modelName = modelName,
                                           fitMethod = fitMethod,
                                           modelYearRange = ModelMinYear:(COVID_YEAR - 1))
    modelBeforeCovid$name <- paste("before", COVID_YEAR)
    
    # 3. model where Covid data are replaced by LC forecast from 2019 for 2020
    mxForecast2020 <- mxModel %>%
      filter(Year < COVID_YEAR) %>%
      bind_rows(forecast(modelBeforeCovid, mxModel, COVID_YEAR)) %>%
      mutate(Name = NULL) %>%
      arrange(Population, Year, Age)
    DxForecast2020 <- DxModel %>% 
      filter(Year < COVID_YEAR) %>%
      bind_rows(mxForecast2020 %>%
                  filter(Year == COVID_YEAR) %>%
                  mutate(Rate = DeathRate * ExModel %>% filter(Year == COVID_YEAR) %>% pull(Exposure)) %>%
                  select(-DeathRate) %>%
                  rename(Deaths = Rate)) %>%
      arrange(Population, Year, Age)
    
    modelReplace2020 <- fitIndividualModel(mx = mxForecast2020, 
                                           Ex = ExModel,     
                                           Dx = DxForecast2020, 
                                           modelName = modelName,
                                           fitMethod = fitMethod,
                                           modelYearRange = ModelMinYear:max(mxForecast2020$Year))
    modelReplace2020$name <- paste(COVID_YEAR - 1, modelName, "estimate for", COVID_YEAR)
    
    # 4. model calibrated on real data up to 2020 and best estimate of LC-2019 for 2021 (assumption: exposure stays constant)
    mxWithCovid2021 <- mxModel %>%
      bind_rows(forecast(modelBeforeCovid, mxModel, COVID_YEAR + 1)) %>%
      mutate(Name = NULL) %>%
      arrange(Population, Year, Age)
    DxWithCovid2021 <- DxModel %>% 
      bind_rows(mxWithCovid2021 %>%
                  filter(Year == COVID_YEAR + 1) %>%
                  mutate(Rate = DeathRate * ExModel %>% filter(Year == COVID_YEAR) %>% pull(Exposure)) %>%  # an alternative would be to estimate an Ex for COVID_YEAR + 1
                  select(-DeathRate) %>%
                  rename(Deaths = Rate)) %>%
      arrange(Population, Year, Age)
    ExWith2021 <- ExModel %>%
      bind_rows(ExModel %>% 
                  filter(Year == COVID_YEAR) %>% 
                  mutate(Year = COVID_YEAR + 1))
    
    modelWithCovid2021 <- fitIndividualModel(mx = mxWithCovid2021, 
                                             Ex = ExWith2021,     
                                             Dx = DxWithCovid2021, 
                                             modelName = modelName,
                                             fitMethod = fitMethod,
                                             modelYearRange = (ModelMinYear + 1):max(mxWithCovid2021$Year))
    modelWithCovid2021$name <- paste0("real ", COVID_YEAR, ", ", COVID_YEAR - 1, " ", modelName, " est. for ", COVID_YEAR + 1)
    
    # 5. model calibrated on real data up to 2019 and best estimates of 2019 for 2020 and 2021 (assumption: exposure stays constant)
    mxForecast2021 <- mxForecast2020 %>%
      bind_rows(forecast(modelBeforeCovid, mxModel, COVID_YEAR + 1)) %>%
      mutate(Name = NULL) %>%
      arrange(Population, Year, Age)
    DxForecast2021 <- DxForecast2020 %>% 
      bind_rows(mxForecast2021 %>%
                  filter(Year == COVID_YEAR + 1) %>%
                  mutate(Rate = DeathRate * ExModel %>% filter(Year == COVID_YEAR) %>% pull(Exposure)) %>%  # an alternative would be to estimate an Ex for COVID_YEAR + 1
                  select(-DeathRate) %>%
                  rename(Deaths = Rate)) %>%
      arrange(Population, Year, Age)
    
    modelTo2021 <- fitIndividualModel(mx = mxForecast2021, 
                                      Ex = ExWith2021,     
                                      Dx = DxForecast2021, 
                                      modelName = modelName,
                                      fitMethod = fitMethod,
                                      modelYearRange = (ModelMinYear + 1):max(mxForecast2021$Year))
    modelTo2021$name <- paste(COVID_YEAR - 1, "LC est. up to", COVID_YEAR + 1)
    
    # Perform comparisons
    assign(paste0("comparisonRealVsReplace_", modelName, "_", fitMethod), 
           compareMortalityModels(modelWithCovid, modelReplace2020, mxModel, ModelEvaluationYears))
    assign(paste0("comparisonReal2021VsReplace_", modelName, "_", fitMethod),
           compareMortalityModels(modelWithCovid2021, modelTo2021, mxModel, ModelEvaluationYears + 1))
  }
}