compareMortalityModels <- function(model1, 
                                   model2, 
                                   mx, 
                                   evaluationYears) {
  #' Compares parameters, forecasts, annuity values and life insurance values of two mortality models
  #'
  #' @param model1 first model object
  #' @param model2 second model object
  #' @param mx data.frame of death rates, columns Population, Age, Year, DeathRate
  #' @param evaluationYears numeric vector of years to evaluate both models on
  #' 
  #' @return list of two lists:
  #'         "data": data.frames dfParameters, dfForecasts, dfAnnuity, dfLifeInsurance
  #'         "plots": annuity and life insurance plot objects
  
  if (class(model1) == "LC" && class(model2) == "LC") {   # available parameters depend on model
    pars <- c("ax", "bx", "kt")
  } else if (class(model1) == "CBD" && class(model2) == "CBD") {
    pars <- c("kt1", "kt2")
  } else {
    pars <- c()
  }
  
  # parameters
  dfParameters <- bind_rows(modelToDf(model1, pars), modelToDf(model2, pars)) %>%
    group_by(x, Population, parameter) 
  
  # forecasts
  dfForecasts <- bind_rows(forecast(model1, mx, evaluationYears, PredictionIntervalProbability),
                           forecast(model2, mx, evaluationYears, PredictionIntervalProbability))
  
  # annuity value and life insurance value
  # We only have death rates for 5-year age groups. For annuity/life insurance values we need death rates per age.
  # Therefore, we assume constant death rates over these 5-year age groups.
  dfAnnuity <- calcProductFromForecasts(dfForecasts, evaluationYears, "annuity", AnnuityDuration, AnnuityStartingAge,
                                        InterestRate)
  dfLifeInsurance <- calcProductFromForecasts(dfForecasts, evaluationYears, "life insurance",
                                              LifeInsuranceDuration, LifeInsuranceStartingAge, InterestRate)
  
  # plots
  plotAnnuity <- plotAnnuityValue(dfAnnuity, sex = "Male")
  if (!is.null(dfLifeInsurance)) {
    plotLifeInsurance <- plotLifeInsuranceValue(dfLifeInsurance, sex = "Male")
  } else {
    plotLifeInsurance <- NULL
  }
  
  plotsPar <- NULL
  
  for (par in pars) {
    plotsPar[[par]] <- .plotModelParameterBySex(df = dfParameters, 
                                                par = par, 
                                                sex = "Male", 
                                                countryNamesOnly = TRUE)
  }
  
  return(list(data = list(dfParameters = dfParameters, 
                          dfForecasts = dfForecasts, 
                          dfAnnuity = dfAnnuity, 
                          dfLifeInsurance = dfLifeInsurance),
              plots = list(plotAnnuity = plotAnnuity,
                           plotLifeInsurance = plotLifeInsurance,
                           plotsPar = plotsPar)))
}


modelToDf <- function(model, pars) {
  #' Writes parameters from a mortality model object into a data.frame
  #'
  #' @param model A mortality model object, i.e., a list containing parameters as specified by pars
  #' @param pars character vector containing parameter names
  #' 
  #' @return data.frame with columns Population, Age, Year, parameter, value
  
  result <- data.frame()
  
  for (par in pars) {
    result <- bind_rows(result,
                        model[[par]] %>%
                          as.data.frame() %>%
                          rownames_to_column(var = ifelse(substr(par, 1, 2) == "kt", "Population", "Age")) %>%
                          pivot_longer(-one_of(ifelse(substr(par, 1, 2) == "kt", "Population", "Age")), 
                                       names_to = ifelse(substr(par, 1, 2) == "kt", "Year", "Population"), 
                                       values_to = "value") %>%
                          mutate(parameter = par))
  }
  
  if (!"Age" %in% colnames(result)) {
    result$Age <- NA
  }
  
  result <- result %>%
    mutate(x = ifelse(is.na(Age), Year, Age)) %>%
    select(-one_of(c("Age", "Year"))) %>% 
    mutate(Name = model$name)
  
  return(result)
}


calcProductFromForecasts <- function(dfForecasts, 
                                     evaluationYears,
                                     productType,
                                     duration,
                                     startingAge,
                                     interestRate) {
  #' Calculates annuity/life insurance present value from forecasts
  #'
  #' @param dfForecasts data.frame containing Year, Age, Population, Name, DeathRate, LowerBound, UpperBound
  #' @param evaluationYears integer vector; years to evaluate  
  #' @param productType "annuity" or "life insurance"
  #' @param duration integer, duration in years
  #' @param startingAge integer
  #' @param interestRate double
  #' 
  #' @return data.frame with annuity/life insurance PVs
  
  if (min(dfForecasts$Age) > startingAge) {
    return(NULL)
  }
  
  df <- dfForecasts %>% 
    filter(Year >= min(evaluationYears), Year <= min(evaluationYears) + duration - 1, 
           Age >= startingAge, Age < startingAge + duration - 1) %>%
    group_by_at(intersect(c("Population", "Name", "Year"), colnames(dfForecasts))) %>%
    complete(Age = startingAge:(startingAge + duration - 1)) %>%
    fill(DeathRate, LowerBound, UpperBound) %>%
    group_by_at(intersect(c("Population", "Name"), colnames(dfForecasts))) %>%
    summarize(getProductValueDf(df = data.frame(Population, Year, Age, DeathRate, LowerBound, UpperBound), 
                                productType = productType,
                                duration = duration, 
                                interestRate = interestRate)) 
  
  if ("Name" %in% colnames(dfForecasts)) {
    df <- df %>%
      arrange(desc(Name)) %>%
      group_by(Population) %>%
      separate(col = Population, into = c("Country", "Sex"), sep = "_", remove = FALSE) %>%
      mutate(Sex = ifelse(Sex == "M", "Male", ifelse(Sex == "F", "Female", "Total")))
  }
  
  return(df)
}



getProductValueDf <- function(df, productType, duration, interestRate) {
  #' filters data.frame such that cohort death rates for calculating maximum-length annuities/term assurance remain 
  #' and then returns data.frame with one row per population containing corresponding cohort annuity/term assurance 
  #' value and possibly also annuity/term assurance value at lower and upper bound
  #'
  #' @param df data.frame containing columns Population, Year, Age, DeathRate and potentially LowerBound, UpperBound
  #' @param productType "annuity" or "life insurance"
  #' @param duration duration of the annuity/term assurance
  #' @param interestRate interest rate to use for annuity/term assurance value calculation
  #' 
  #' @return data.frame with best estimate cohort annuity/term assurance value and possibly also annuity/term assurance 
  #'         value at lower and upper bound
  
  # ensure that enough ages are contained in the data
  if (duration > n_distinct(df$Age)) {
    warning(paste0("getProductValueDf: Not enough ages (", n_distinct(df$Age), ") for duration ", duration))
    duration <- n_distinct(df$Age)
  }
  
  if (!("LowerBound" %in% colnames(df))) {
    df$LowerBound <- NA
    df$UpperBound <- NA
    dropBounds <- TRUE 
  } else {
    dropBounds <- FALSE
  }
  
  df <- df %>%
    filter(Year < min(Year) + duration & Age > max(Age) - duration) %>%
    filter((Year - Age) == (min(Year) - min(Age))) %>%
    group_by(Population) 
  
  # calculate product value
  if (productType == "annuity") {
    # annuity-due, PV = sum_{k=1}^n kpx(t) * v^k
    # cohort survival probabilities px(t), 2px(t) = px(t) * p_{x+1}(t+1), ..., npx(t) as function of cohort death rates mx(t), m_{x+1}(t+1), ..., m_{x+n}(t+n)
    p <- function(m) {    
      exp(-cumsum(m))
    }
    
    df <- df %>%
      dplyr::summarize(AnnuityValue = sum((1/(1 + interestRate))^(1:duration) * p(DeathRate)),
                       AnnuityValueLower = sum((1/(1 + interestRate))^(1:duration) * p(UpperBound)),
                       AnnuityValueUpper = sum((1/(1 + interestRate))^(1:duration) * p(LowerBound))) %>% 
      ungroup()
    
    if (dropBounds) {
      df$AnnuityValueLower <- dfAnnuity$AnnuityValueUpper <- NULL
    }
  } else {
    # cohort survival probabilities 1, px(t), 2px(t) = px(t) * p_{x+1}(t+1), ..., {n-1}px(t) as function of cohort death rates mx(t), m_{x+1}(t+1), ..., m_{x+n}(t+k)
    # and cohort death probabilities qx(t), q{x+1}(t+1), ..., q{x+n}(t+n) as function of these cohort death rates
    p <- function(m) {    
      c(1, (exp(-cumsum(m)))[1:(length(m) - 1)])
    }
    
    q <- function(m) {
      1 - exp(-m)
    }
    
    df <- df %>%
      dplyr::summarize(LifeInsuranceValue = sum((1/(1 + interestRate))^(1:duration) * p(DeathRate) * q(DeathRate)),
                       LifeInsuranceValueLower = sum((1/(1 + interestRate))^(1:duration) * p(LowerBound) * q(LowerBound)),
                       LifeInsuranceValueUpper = sum((1/(1 + interestRate))^(1:duration) * p(UpperBound) * q(UpperBound))) %>% 
      ungroup()
    
    if (dropBounds) {
      df$LifeInsuranceValueLower <- df$LifeInsuranceValueUpper <- NULL
    }
  }
  
  return(df)
}


plotAnnuityValue <- function(dfAnnuity, 
                             sex) {
  #' Plots annuity values of two different models 
  #'
  #' @param dfAnnuity data.frame containing columns Country, Sex, AnnuityValue, Name
  #' @param sex string containing sex to filter on (or NULL)
  #' 
  #' @return plot object
  
  if (!is.null(sex)) {
    dfAnnuity <- dfAnnuity %>% filter(Sex == sex)
  }
  
  plotAnnuity <- ggplot(dfAnnuity, aes(x = Country, y = AnnuityValue, color = Name, shape = Name)) +
    geom_point(size = PlotPointSize, position = position_dodge(width = ifelse(PredictionIntervalProbability > 0, 0.75, 0))) +
    labs(x = "country", y = "annuity value", color = "model", shape = "model") 
  
  if (is.null(sex)) {
    plotAnnuity <- plotAnnuity + facet_wrap(. ~ Sex, scales = "fixed") 
  }
  
  # plot annuity values with prediction intervals
  if (PredictionIntervalProbability > 0) {
    plotAnnuity <- plotAnnuity +
      geom_errorbar(aes(ymin = AnnuityValueLower, ymax = AnnuityValueUpper), 
                    size = PlotLineSize,
                    position = position_dodge(width = 0.75))   
  }
  
  return(plotAnnuity)
}


plotLifeInsuranceValue <- function(dfLifeInsurance,
                                   sex) {
  #' Plots life insurance values of two different models 
  #'
  #' @param dfLifeInsurance data.frame containing columns Country, Sex, LifeInsuranceValue, Name
  #' @param sex string containing sex to filter on (or NULL)
  #' 
  #' @return plot object
  
  if (!is.null(sex)) {
    dfLifeInsurance <- dfLifeInsurance %>% filter(Sex == sex)
  }
  
  plotLifeInsurance <- ggplot(dfLifeInsurance, aes(x = Country, y = LifeInsuranceValue, color = Name, shape = Name)) +
    geom_point(size = PlotPointSize + 2, position = position_dodge(width = ifelse(PredictionIntervalProbability > 0, 0.75, 0))) +
    labs(x = "country", y = "life insurance value", color = "model", shape = "model") 
  
  if (is.null(sex)) {
    plotLifeInsurance <- plotLifeInsurance + facet_wrap(. ~ Sex, scales = "fixed") 
  }
  
  # plot life insurance values with prediction intervals
  if (PredictionIntervalProbability > 0) {
    plotLifeInsurance <- plotLifeInsurance +
      geom_errorbar(aes(ymin = LifeInsuranceValueLower, ymax = LifeInsuranceValueUpper), 
                    size = PlotLineSize,
                    position = position_dodge(width = 0.75))  
  }
  
  return(plotLifeInsurance)
}


.plotModelParameterBySex <- function(dfParameters, 
                                     par, 
                                     sex, 
                                     countryNamesOnly) {
  #' Plots specified parameter of two different models for a particular sex
  #'
  #' @param dfParameters data.frame containing columns Population, parameter, x, value, Name
  #' @param par string specifying parameter name
  #' @param sex string specifying sex
  #' @param countryNamesOnly boolean, if TRUE, no population names but only country names are shown
  #' 
  #' @return plot object
  
  dfPlot <- dfParameters %>% 
    filter(str_sub(Population, -1, -1) == substr(sex, 1, 1),
           parameter == par) %>%
    mutate(Population = ifelse(countryNamesOnly, gsub(paste0("_", substr(sex, 1, 1)), "", Population), Population))
  
  plotParameters <- ggplot(dfPlot,
                           aes(x = as.numeric(x), y = value, color = Name, shape = Name)) +
    geom_point(size = PlotPointSize) +
    facet_wrap(. ~ Population, scales = "fixed") +
    labs(x = ifelse(substr(par, 1, 2) == "kt", "year", "age"), y = par, color = "model", shape = "model") 
  
  return(plotParameters)
}