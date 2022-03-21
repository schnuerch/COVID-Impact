downloadData <- function(username, 
                         password, 
                         targetFolder,
                         quantities = c("Dx", "Ex", "mx"), 
                         countries = NULL, 
                         interval = "5x1") {
  #' Downloads specified quantities for specified countries from human mortality database 
  #' and writes them to disk; working directory has to be source file location
  #'
  #' @param username HMD user name (registration on www.mortality.org required)
  #' @param password HMD password (registration on www.mortality.org required)
  #' @param targetFolder folder name relative to working directory where data should be saved
  #' @param quantities e.g. c("Dx", "Ex", "mx") for numbers of deaths, exposures and death rates 
  #' @param countries e.g. NULL for all available countries
  #' @param interval "1x1" by age and year; "5x1" by 5-year age group and year
  #' 
  #' @return void
  
  # initialize
  dataHMD <- list()
  timestamp <- Sys.time()
  
  # read all data from HMD 
  for(quantity in quantities) {
    dataHMD[[quantity]] <- (ReadHMD(what = quantity, 
                                    countries = countries,
                                    interval = interval,
                                    username = username, 
                                    password = password))$data
  }
  
  # write to disk
  save(dataHMD, file = paste0(file.path(targetFolder, paste0(format(timestamp, "%Y-%m-%d_%H.%M"), "_dataHMD")), 
                              ifelse(interval == "1x1", "", interval), ".RData"))
  
  return()
}


prepareDataESP <- function(filePath, cappingAge = Inf, lengthAgeRanges = NULL) {
  #' Prepares European standard population data 
  #'
  #' @param filePath path to data source
  #' @param cappingAge integer or Inf; age above which data for all older age groups are summarized 
  #'                   to the single age group "CappingAge+", e.g., 90+. If it is Inf, age groups are left as they come.
  #' @param lengthAgeRanges integer or NULL; desired length of all ages ranges, e.g., 5 will create age ranges 0-5, 5-10, ...
  #' 
  #' @return data.frame with prepared European standard population data
  
  # read and rename
  df <- fread(filePath) %>%
    rename(Age = "AgeGroup", Exposure = "EuropeanStandardPopulation") 
  
  # adjust age column
  df <- suppressWarnings(df %>% 
                           separate(Age, c("Age", NA))) %>%
    mutate(Age = as.numeric(gsub("([0-9]+).*$", "\\1", Age))) %>%   # e.g., 110+ --> 110
    arrange(Age)
  
  # age capping
  if (!is.infinite(cappingAge) && max(df$Age) >= cappingAge) {
    df <- df %>% 
      group_by(Age = pmin(Age, cappingAge)) %>%   # new Age column: capped at capping Age
      summarize(Exposure = sum(Exposure)) %>%
      ungroup()
  }
  
  # (try to) force age ranges of specified length
  if (!is.null(lengthAgeRanges)) {
    targetLengthsAgeRanges <- seq(min(df$Age), max(df$Age), by = lengthAgeRanges)
    
    if (length(setdiff(targetLengthsAgeRanges, unique(df$Age))) > 0) {
      warningMessage <- paste("prepareDataESP: Specified age range length of", lengthAgeRanges, "could not be implemented",
                              "because not all target ages from", min(df$Age), "to", 
                              max(df$Age), "are contained in the ESP data.")
      warning(warningMessage)
    } else {
      df <- df %>%  
        # sum up all deaths belonging to the age groups defined by lengthAgeRanges; create new age groups
        group_by(Age = floor(Age / lengthAgeRanges) * lengthAgeRanges) %>%
        summarize(Exposure = sum(Exposure))
    }
  }
  
  # share in total
  df$ExposureFraction <- df$Exposure / sum(df$Exposure)
  
  return(df)
}


prepareDataHMD <- function(path, 
                           dfStandardPopulation, 
                           countryMapping, 
                           minYear = 0, 
                           maxYear = Inf, 
                           minAge = 0 ,
                           maxAge = Inf, 
                           sexes = c("Male", "Female", "Total"), 
                           cappingAge = Inf, 
                           extendHistoryGermany = "Both") {
  #' Loads and preprocesses HMD data
  #'
  #' @param path path to .RData file containing a list "dataHMD" of HMD data
  #' @param dfStandardPopulation European standard population data.frame created by prepareDataESP
  #' @param countryMapping data.frame to join full names to df (with columns rawName and fullName)
  #' @param minYear Minimal year to include
  #' @param maxYear Maximal year to include
  #' @param minAge Minimal age to include
  #' @param maxAge Maximal age to include
  #' @param sexes Sexes to include (Male, Female, Total). 
  #' @param cappingAge integer or Inf; age above which data for all older age groups are summarized 
  #'                   to the single age group "CappingAge+", e.g., 90+. If it is Inf, age groups are left as they are.
  #' @param extendHistoryGermany FALSE: no data extension for Germany
  #'                             "West": extend German data with West German data
  #'                             "Both": extend German data with West+East German data
  #' 
  #' @return preprocessed HMD data.frame
  
  load(path)
  
  mappingTableHMD <- data.frame(rawName = c("mx", "Ex", "Dx"), 
                                displayName = c("DeathRateHMD", "ExposureHMD", "DeathsHMD"),
                                stringsAsFactors = FALSE)

  for (data in c("mx", "Ex", "Dx")) {
    dfNew <- filterDataHMD(df = dataHMD[[data]], 
                           countryMapping = countryMapping,
                           minYear = minYear, maxYear = maxYear, 
                           minAge = minAge, maxAge = maxAge, 
                           sexes = sexes, 
                           imposeNonZeroRates = (data == "mx")) %>%
      rename(!!mappingTableHMD$displayName[mappingTableHMD$rawName == data] := Rate) %>%
      mutate(Population = factor(if(length(Sexes) == 1) Country else paste(Country, substr(Sex, 1, 1), sep = "_")))
    
    if (data == "mx") {
      df <- dfNew
    } else {
      df <- df %>% 
        left_join(dfNew, by = c("Year", "Age", "Country", "Sex", "Population"))
    }
  }
  
  # extend data for Germany backwards
  if (is.character(extendHistoryGermany)) {
    minYearGermany <- min(df %>% filter(Country == "Germany") %>% pull(Year))
    
    if (extendHistoryGermany == "West") {
      df <- df %>% 
        bind_rows(df %>%
                    filter(Country == "West Germany", Year < minYearGermany) %>%
                    mutate(Country = "Germany",
                           Population = paste(Country, substr(Sex, 1, 1), sep = "_")))
    } else {
      df <- df %>%
        bind_rows(df %>%
                    filter(Country %in% c("West Germany", "East Germany"), Year < minYearGermany) %>%
                    mutate(Country = "Germany",
                           Population = paste(Country, substr(Sex, 1, 1), sep = "_")) %>%
                    group_by(Population, Country, Sex, Age, Year) %>%
                    summarize(DeathsHMD = sum(DeathsHMD), ExposureHMD = sum(ExposureHMD)))
    }
    
    df <- df %>%
      mutate(DeathRateHMD = DeathsHMD / ExposureHMD)
  }
  
  df <- df %>% 
    filter(Country != "West Germany", Country != "East Germany")
  
  # age capping
  if (!is.infinite(cappingAge) && max(df$Age) >= cappingAge) {
    df <- df %>% 
      group_by(Population, Country, Sex, Year, Age = pmin(Age, cappingAge)) %>%   # new Age column: capped at capping Age
      summarize(DeathsHMD = sum(DeathsHMD), ExposureHMD = sum(ExposureHMD), 
                DeathRateHMD = pmax(DeathsHMD / ExposureHMD, sqrt(.Machine$double.eps))) %>%
      ungroup()
  }
  
  df <- df %>%
    mutate(Sex = as.factor(Sex), Country = as.factor(Country), Population = as.factor(Population)) %>%
    group_by(Population, Age) %>%
    mutate(ImprovementRateHMD = 1 - DeathRateHMD / data.table::shift(DeathRateHMD, n = 1, fill = NA, type = "lag")) %>%
    ungroup()
  
  # add age standardized rates
  df <- createDfWithAgeStandardizedRate(df,
                                        dfStandardPopulation,
                                        colNamesGrouping = c("Population", "Year"),
                                        colNameRate = "DeathRateHMD",
                                        colNameASR = "AgeStandardizedRateHMD")
  
  df <- createDfWithAgeStandardizedRate(df,
                                        dfStandardPopulation,
                                        colNamesGrouping = c("Population", "Year"),
                                        colNameRate = "ImprovementRateHMD",
                                        colNameASR = "AgeStandardizedImprovementRateHMD")
  
  return(df)
}


prepareDataSTMFInput <- function(folder,
                                 countryMapping, 
                                 sexMapping, 
                                 minYear = 0, 
                                 maxYear = Inf, 
                                 minAge = 0, 
                                 maxAge = Inf, 
                                 sexes = c("Male", "Female", "Total"),
                                 cappingAge = Inf,
                                 lengthAgeRanges = NULL) {
  #' Filters and prepares short-term mortality fluctuation data ("input" data from STMF, containing no death rates)
  #'
  #' @param folder string, folder containing input csv data (relative to source file location)
  #' @param countryMapping data.frame to join full names to df (with columns rawName and fullName)
  #' @param sexMapping data.frame to join full names to df (with columns rawName and fullName)
  #' @param minYear Minimal year to include
  #' @param maxYear Maximal year to include
  #' @param minAge Minimal age to include
  #' @param maxAge Maximal age to include
  #' @param sexes Sexes to include (Male, Female, Total). 
  #' @param cappingAge integer or Inf; age above which data for all older age groups are summarized 
  #'                   to the single age group "CappingAge+", e.g., 90+. If it is Inf, age groups are left as they come.
  #' @param lengthAgeRanges integer or NULL; desired length of all ages ranges, e.g., 5 will create age ranges 0-5, 5-10, ...
  #' 
  #' @return data.frame with filtered and prepared short-term mortality data
  
  # read in files corresponding to relevant countries
  countryCodes <- countryMapping$rawName
  if (length(countryCodes) == 0) {   # read in all countries
    countryPattern <- NULL
  } else {
    countryPattern <- paste0("^(", paste(countryCodes, collapse = "|"), ")stmf.csv")
  }
 
  df <- list.files(path = folder, pattern = countryPattern, full.names = TRUE) %>%
    map_df(~ fread(., colClasses = c(Week = "character", Deaths = "character"))) %>%
    unique %>%
    mutate(Deaths = as.numeric(Deaths))
  
  # drop unnecessary columns
  df <- df %>%
    select(-one_of(c("Area", "AgeInterval", "Type", "Access"))) 
  
  # filter by year
  df <- df %>%
    filter(Year >= minYear, Year <= maxYear) 
  
  # country and sex mapping
  df <- df %>% 
    mutate(PopCode = as.character(PopCode)) %>%   # avoid join warning
    left_join(countryMapping %>% mutate(rawName = as.character(rawName)), by = c(PopCode = "rawName")) %>%
    select(-PopCode) %>%
    rename(Country = fullName)
  
  df <- df %>%
    mutate(Sex = as.character(Sex)) %>%
    left_join(sexMapping %>% mutate(rawName = as.character(rawName)), by = c(Sex = "rawName")) %>%
    select(-Sex) %>%
    rename(Sex = fullName)
  
  # filter sexes
  df <- df %>% 
    filter(Sex %in% sexes)
  
  # print a warning if there are multiple entries per Country, Sex, Year, Week, Age combination  
  dfMultiEntries <- df %>% 
    group_by(Country, Sex, Year, Week, Age) %>% 
    summarize(n = n()) %>% 
    filter(n > 1) %>%
    mutate(combination = paste(Country, Sex, Year, Week, Age))
  
  if (nrow(dfMultiEntries) > 0) {
    warningMessage <- paste0("prepareDataSTMFInput: ", nrow(dfMultiEntries), " combinations of Country, Sex, Year, Week, Age ",
                             "have multiple death counts associated to them. Here are the first such combinations: ",
                             paste(dfMultiEntries[1:min(5, nrow(dfMultiEntries)), ] %>% pull(combination), collapse = ", "),
                             ". Please double check the input data!")
    warning(warningMessage)
  }
  
  # rescale deaths per age group to total count, see https://timriffe.github.io/covid_age/GettingStarted.html (merge with STMF)
  # "redistribute deaths at unknown ages"
  df <- df %>% 
    group_by(Country, Sex, Year, Week) %>%
    mutate(TOT = min(Deaths[Age == "TOT"])) %>%    # column containing total deaths per age group (if we have multiple death counts, 
                                                   # this is peculiar, but we then take the smallest.)
    filter(Age != "TOT", Age != "UNK") %>%    # remove total and unknown counts from age column
    mutate(dist = Deaths / sum(Deaths),       # redistribute total death counts to age groups
           Deaths = dist * TOT) %>%
    select(-one_of(c("dist", "TOT"))) %>%
    ungroup()
  
  # also redistribute deaths at unknown weeks equally over the whole year
  df <- df %>%
    group_by(Country, Sex, Year, Age) %>%
    mutate(TOT = sum(Deaths)) %>%
    filter(Week != "UNK") %>%
    mutate(dist = Deaths / sum(Deaths),
           Deaths = dist * TOT) %>%
    select(-one_of(c("dist", "TOT"))) %>%
    mutate(Week = as.numeric(Week)) %>%
    ungroup()
  
  # prepare age column
  df <- df %>%
    mutate(Age = as.numeric(gsub("([0-9]+).*$", "\\1", Age))) %>%   # e.g., 85p --> 85
    arrange(Country, Year, Week, Age)
  
  # filter by age
  # if minAge is contained in some age interval (and not its lower bound), remove this age interval
  # if maxAge is contained in some age interval, remove this age interval
  df <- df %>%
    filter(Age >= minAge, c(unique(df$Age), Inf)[match(df$Age, unique(df$Age)) + 1] <= maxAge)  # filter "right-shifted" ages
  
  # age capping
  if (!is.infinite(cappingAge) && max(df$Age) >= cappingAge) {
    df <- df %>% 
      group_by(Country, Sex, Year, Week, Age = pmin(Age, cappingAge)) %>%   # new Age column: capped at capping Age
      summarize(Deaths = sum(Deaths)) %>%
      ungroup()
  }
  
  # (try to) force age ranges of specified length
  if (!is.null(lengthAgeRanges)) {
    targetLengthsAgeRanges <- seq(min(df$Age), max(df$Age), by = lengthAgeRanges)
    
    dfFailedTarget <- df %>%
      group_by(Country) %>%
      summarize(Failed = length(setdiff(targetLengthsAgeRanges, unique(Age))) > 0)
    
    if (any(dfFailedTarget$Failed)) {
      warningMessage <- paste("prepareDataSTMFInput: Specified age range length of", lengthAgeRanges, "could not be implemented",
                              "for", sum(dfFailedTarget$Failed), "countries because not all target ages from", min(df$Age), "to", 
                              max(df$Age), "are contained in the STMF data.")
      warning(warningMessage)
    } 
    
    # bind unchanged rows for countries where forced age ranges are impossible with changed rows for countries where they are possible
    df <- bind_rows(df %>% filter(Country %in% dfFailedTarget$Country[dfFailedTarget$Failed]),
                    df %>%  
                      filter(Country %in% dfFailedTarget$Country[!dfFailedTarget$Failed]) %>%
                      # sum up all deaths belonging to the age groups defined by lengthAgeRanges; create new age groups
                      group_by(Country, Sex, Year, Week, Age = floor(Age / lengthAgeRanges) * lengthAgeRanges) %>%
                      summarize(Deaths = sum(Deaths)))
  }
  
  # population column
  df$Population <- paste(df$Country, substr(df$Sex, 1, 1), sep = "_")
  
  # exclude populations where ages do not span the full range
  fullRangePopulations <- df %>% 
    group_by(Population) %>% 
    summarize(nAges = n_distinct(Age)) %>%
    filter(nAges == n_distinct(df$Age)) %>% 
    pull(Population)
  
  if (n_distinct(fullRangePopulations) < n_distinct(df$Population)) {
    warningMessage <- paste("prepareDataSTMFInput: Populations", 
                            paste(setdiff(unique(df$Population), unique(fullRangePopulations)), collapse = ", "),
                            "were removed from the data because they did not contain the full age range.")
    warning(warningMessage)
    
    df <- df %>% filter(Population %in% fullRangePopulations)
  }

  # add YearWeek column
  df <- df %>% 
    mutate(YearWeek = yearweek(paste(Year, "week", Week)))
  
  # rename Deaths column
  df <- df %>%
    rename(DeathsSTMFInput = "Deaths")
  
  # add yearly deaths column
  df <- df %>% 
    ungroup() %>%
    group_by(Country, Sex, Year, Age) %>%
    mutate(DeathsSTMFInputYear = sum(DeathsSTMFInput, na.rm = TRUE))
  
  # make Population, Country and Sex factors
  df <- df %>%
    ungroup() %>%
    mutate(Population = as.factor(Population), Country = as.factor(Country), Sex = as.factor(Sex))
  
  return(df)
}


createDfTotal <- function(dfHMD, 
                          dfSTMFInput, 
                          nBaseYearsExposureExtrapolation) {
  #' Merges and further prepares HMD and STMF data
  #'
  #' @param dfHMD HMD data.frame created by prepareDataHMD
  #' @param dfSTMFInput STMF data.frame created by prepareDataSTMFInput
  #' @param nBaseYearsExposureExtrapolation number of previous years to use for exposure extrapolation
  #' 
  #' @return preprocessed data.frame
  
  dfTotal <- dfSTMFInput %>%
    full_join(dfHMD, by = c("Population", "Country", "Sex", "Year", "Age")) %>%
    mutate(Country = factor(Country, levels = sort(levels(Country))))
  
  # extrapolate HMD exposures, calculate STMF exposures and death rates
  dfTotal <- dfTotal %>% 
    group_by(Population, Age) %>%
    arrange(Year) %>%
    mutate(ExposureHMD = extrapolateLinearReg(x = Year, y = ExposureHMD, nBaseY = nBaseYearsExposureExtrapolation)) %>%   # linear extrapolation (based on last nBaseY values)
    group_by(Year) %>%
    mutate(ExposureWeek = ExposureHMD / 52,   # 52 weeks per year
           DeathRateSTMF = DeathsSTMFInput / ExposureWeek,
           DeathRateSTMFYear = DeathsSTMFInputYear / ExposureHMD)  
  
  # add improvement rate and death rate according to both HMD and STMF
  dfTotal <- dfTotal %>% 
    left_join(dfTotal %>%
                filter(Week == 1 | Year < (COVID_YEAR - 1)) %>%
                select(Population, Age, Year, DeathRateHMD, DeathRateSTMFYear, DeathsHMD, DeathsSTMFInputYear) %>%
                unique() %>%
                group_by(Population, Age, Year) %>%       
                mutate(DeathRate = unique(ifelse(is.na(DeathRateHMD), DeathRateSTMFYear, DeathRateHMD)),
                       Deaths = unique(ifelse(is.na(DeathsHMD), DeathsSTMFInputYear, DeathsHMD))) %>%
                group_by(Population, Age) %>%
                mutate(ImprovementRate = 1 - DeathRate / data.table::shift(DeathRate, n = 1, fill = NA, 
                                                                           type = "lag")) %>%
                select(-one_of(c("DeathRateHMD", "DeathRateSTMFYear", "DeathsHMD", "DeathsSTMFInputYear"))), 
              by = c("Population", "Age", "Year"))
  
  return(dfTotal)
}


filterDataHMD <- function(df, 
                          countryMapping, 
                          minYear = 0, 
                          maxYear = Inf, 
                          minAge = 0, 
                          maxAge = Inf, 
                          sexes = c("Male", "Female", "Total"),
                          imposeNonZeroRates = TRUE) {
  #' Filters and prepares data from human mortality database
  #'
  #' @param df data.frame from HMD
  #' @param countryMapping data.frame to join full names to df (with columns rawName and fullName)
  #' @param minYear Minimal year to include
  #' @param maxYear Maximal year to include
  #' @param minAge Minimal age to include
  #' @param maxAge Maximal age to include
  #' @param sexes Sexes to include (Male, Female, Total). 
  #' @param imposeNonZeroRates boolean: replace zero values by small positive value
  #' 
  #' @return data.frame with filtered and prepared HMD mortality data
  
  # filter by year, countries
  df <- df %>%
    filter(Year >= minYear, Year <= maxYear) %>% 
    filter(country %in% countryMapping$rawName)
  
  # filter by age
  if (is.factor(df$Age)) {              # if age is factor (= age groups), map it to lower bound and adjust filtering
    df <- suppressWarnings(df %>% 
                             separate(Age, c("Age", NA))) %>%
      mutate(Age = as.numeric(gsub("([0-9]+).*$", "\\1", Age))) %>%   # e.g., 110+ --> 110
      arrange(country, Year, Age)
    
    # if minAge is contained in some age interval (and not its lower bound), remove this age interval
    # if maxAge is contained in some age interval, remove this age interval
    df <- df %>%
      filter(Age >= minAge, c(unique(df$Age), Inf)[match(df$Age, unique(df$Age)) + 1] <= maxAge)  # filter "right-shifted" ages
  } else {
    df <- df %>%
      filter(Age >= minAge, Age <= maxAge)
  }
  
  # countries get full names
  df <- df %>% 
    mutate(country = as.character(country)) %>%   # avoid join warning
    left_join(countryMapping %>% mutate(rawName = as.character(rawName)), by = c(country = "rawName")) %>%
    select(-country) %>%
    rename(Country = fullName)
  
  # filter sexes
  df <- df %>% 
    gather(sex, Rate, !!sexes) %>%
    rename(Sex = sex) %>% 
    filter(Sex %in% sexes)
  
  # cast Rate to numeric and exclude rows with NA rates
  df <- df %>% 
    mutate(Rate = as.numeric(Rate)) %>%
    filter(!is.na(Rate))
  
  # replace zero death rates with small positive value (many models are based on log death rates)
  if (imposeNonZeroRates) {
    df[df$Rate == 0, "Rate"] <- sqrt(.Machine$double.eps)
  }
  
  return(df)
}


createDfWithAgeStandardizedRate <- function(df, 
                                            dfStandardPopulation,
                                            colNamesGrouping,
                                            colNameRate = "DeathRate", 
                                            colNameASR = "ASR") {
  #' Adds column to df with age-standardized rates (not necessarily death rates). 
  #' Reference population is given in dfStandardPopulation and could for example be the European Standard Population.
  #'
  #' @param df data.frame containing columns given by colNameDeathRate, colNameExposure
  #' @param dfStandardPopulation data.frame containing columns Age, ExposureFraction
  #' @param colNamesGrouping character vector containing names of grouping columns in df except for Age
  #' @param colNameRate character containing name of rate column in df
  #' @param colNameASR character containing desired name of age-standardized rate column in result df
  #' 
  #' @return data.frame containing age standardized rates
  
  # Adjust fractions in standard population
  dfStandardPopulation <- dfStandardPopulation %>% 
    filter(Age %in% df$Age) %>%
    mutate(ExposureFraction = ExposureFraction / sum(ExposureFraction))
  
  # Calculate ASR
  df <- df %>% 
    left_join(dfStandardPopulation %>% select(Age, ExposureFraction), by = "Age") %>%
    mutate(!!colNameASR := ExposureFraction * !!rlang::sym(colNameRate)) %>%
    group_by_at(colNamesGrouping) %>%
    mutate(!!colNameASR := sum(!!rlang::sym(colNameASR))) %>%
    ungroup() %>%
    select(-ExposureFraction)
  
  return(df)
}


extrapolateLinearReg <- function(x, y, nBaseY) {
  #' Extrapolates y where it is NA from a regression of y on x (with last nBaseY available years)
  #'
  #' @param x vector, regression independent variable
  #' @param y vector of same length as x, regression dependent variable, contains NA at its end
  #' @param nBaseY number of past observations of y to use for regression
  #' 
  #' @return y with imputed values
  if (length(x) != length(y)) {
    errorMessage <- paste0("extrapolateLinearReg: lengths of x (", length(x), ") and y (", length(y), ") differ.")
    logerror(errorMessage)
    stop(errorMessage)
  }
  
  df <- unique(data.frame(x = x, y = y))
  linearModel <- lm(y ~ x, df %>% filter(!is.na(y)) %>% tail(nBaseY))
  y[is.na(y)] <- predict(linearModel, data.frame(x = x[is.na(y)]))
  
  return(y)
}
