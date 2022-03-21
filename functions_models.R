fitIndividualModel <- function(mx, 
                               Ex, 
                               Dx, 
                               modelName, 
                               fitMethod,
                               modelYearRange = NULL) {
  #' Fits individual (population-wise) stochastic mortality models to data.frame of death rates
  #'
  #' @param mx data.frame of death rates containing 4 columns: Population, Year, Age, Rate or DeathRate
  #' @param Ex data.frame of central exposures containing 4 columns: Population, Year, Age, Exposure
  #' @param Dx data.frame of deaths, contains Population, Year, Age, Deaths
  #' @param modelName string, "LC" or "CBD"
  #' @param fitMethod string, "SVD" or "MLE"
  #' @param modelYearRange integer vector, years on which the actual model should be calibrated. NULL for no filtering.
  #' 
  #' @return fitted model (list of e.g. ax [AxP], bx [Axk], kt [PxY], tsModel)
  
  if (!is.null(modelYearRange)) {
    mxModel <- mx %>% filter(Year %in% modelYearRange)
    DxModel <- Dx %>% filter(Year %in% modelYearRange)
    ExModel <- Ex %>% filter(Year %in% modelYearRange)
  }
  
  if (modelName == "LC") {
    model <- fitIndividualLC(mx = mxModel, Ex = ExModel, Dx = DxModel, fitMethod = fitMethod)
  } else if (modelName == "CBD") {
    model <- fitIndividualCBD(mx = mxModel, Ex = ExModel, Dx = DxModel, fitMethod = fitMethod)
  } 
  
  return(model)
}


fitIndividualLC <- function(mx,
                            Ex, 
                            Dx,
                            fitMethod) {
  #' Fits individual (population-wise) Lee-Carter models to data.frame of death rates and returns list of parameters in matrices
  #' Identifiability constraints: sum(bx) = 1, k_{t_1} = 0
  #'
  #' @param mx data.frame of death rates containing 4 columns: Population, Year, Age, Rate or DeathRate
  #' @param Ex data.frame of central exposures containing 4 columns: Population, Year, Age, Exposure
  #' @param Dx data.frame of deaths, contains Population, Year, Age, Deaths
  #' @param fitMethod string, "SVD" or "MLE"
  #' 
  #' @return fitted LC model (list of ax [AxP], bx [Axk], kt [PxY], random walk parameters)
  #' @export
  #'
  #' @examples
  
  populations <- sort((unique(mx$Population)))
  
  if ("DeathRate" %in% colnames(mx)) {
    colnames(mx)[colnames(mx) == "DeathRate"] <- "Rate"  # unify column name
  } 
  
  # initialize
  ax <- bx <- matrix(NA, ncol = length(populations), nrow = n_distinct(mx$Age), 
                     dimnames = list(unique(mx$Age), populations))
  kt <- matrix(NA, ncol = n_distinct(mx$Year), nrow = length(populations), 
               dimnames = list(populations, sort(unique(mx$Year))))
  
  if (fitMethod == "SVD") {
    for (pop in populations) {
      mxPop <- mx %>% filter(Population == pop)
      
      temp <- calcCentralisedLogDeathRates(mxPop)
      ax[, pop] <- temp$ax
      mx_prepared_svd <- as.matrix(temp$mx_log_centralised %>% dplyr::select(-c("Population", "Age")))
      
      # population-wise singular value decomposition
      decomposition <- svd(mx_prepared_svd %*% t(mx_prepared_svd))
      bx[, pop] <- decomposition$u[, 1] / sum(decomposition$u[, 1])
      
      popAvailableYears <- as.character(unique(mx %>% filter(Population == pop) %>% arrange(Year) %>% pull(Year)))
      kt[pop, popAvailableYears] <- solve(t(bx[, pop]) %*% bx[, pop]) %*% t(bx[, pop]) %*% mx_prepared_svd
      
      # implement identifiability constraint k_{t_1} = 0 and sum(b_x) = 1
      ax[, pop] <- ax[, pop] + bx[, pop] * kt[pop, min(which(!is.na(kt[pop, ])))]
      kt[pop, ] <- kt[pop, ] - kt[pop, min(which(!is.na(kt[pop, ])))]
    }
  } else if (fitMethod == "MLE"){
    loglik <- 0
    
    for (pop in populations) {
      # prepare input (age-year matrices)
      Dxt <- Dx %>% 
        filter(Population == pop) %>%
        mutate(Population = NULL) %>%
        pivot_wider(names_from = "Year", values_from = "Deaths") %>%
        column_to_rownames(var = "Age")
      
      Ext <- Ex %>% 
        filter(Population == pop) %>%
        mutate(Population = NULL) %>%
        pivot_wider(names_from = "Year", values_from = "Exposure") %>%
        column_to_rownames(var = "Age")
      
      model <- StMoMo::fit(lc(), Dxt = Dxt, Ext = Ext, ages = unique(mx$Age), years = unique(mx$Year), verbose = FALSE)
      
      popAvailableYears <- as.character(unique(mx %>% filter(Population == pop) %>% arrange(Year) %>% pull(Year)))
      ax[, pop] <- model$ax
      bx[, pop] <- model$bx
      kt[pop, popAvailableYears] <- model$kt
      loglik <- loglik + model$loglik
    }
    
    model <- list(ax = ax, bx = bx, kt = kt, loglik = loglik)
  } else {
    stop(paste0("fitIndividualLC: invalid fit method '", fitMethod, "' specified"))
  }
  
  # fit rwd 
  tsModel <- fitTimeSeries(kt, arimaOrder = c(0, 1, 0))
  
  model <- list(ax = ax, bx = bx, kt = kt, tsModel = tsModel)
  class(model) <- "LC"
  model$fitMethod <- fitMethod
  
  return(model)
}


fitIndividualCBD <- function(mx, 
                             Ex,
                             Dx,
                             fitMethod) {
  #' Fits individual CBD models to data.frame of death rates and returns list of parameters in matrices
  #'
  #' @param mx data.frame of death rates containing 4 columns: Population, Year, Age, Rate or DeathRate
  #' @param Ex data.frame of central exposures containing 4 columns: Population, Year, Age, Exposure
  #' @param Dx data.frame of deaths, contains Population, Year, Age, Deaths
  #' @param fitMethod string, "SVD" or "MLE"
  #' 
  #' @return list of fitted CBD models (StMoMo objects)
  #' @export
  #'
  #' @examples
  
  if (!fitMethod %in% c("MLE", "lm")) {
    errorMessage <- paste("fitIndividualCBD: fitMethod", fitMethod, "not implemented.")
    stop(errorMessage)
  }
  
  # initialize
  populations <- sort((unique(mx$Population)))
  
  if ("DeathRate" %in% colnames(mx)) {
    colnames(mx)[colnames(mx) == "DeathRate"] <- "Rate"  # unify column name
  } 
  
  models <- list()
  kt1 <- kt2 <- matrix(NA, nrow = length(populations), ncol = n_distinct(mx$Year), 
                       dimnames = list(populations, unique(mx$Year)))
  
  for (pop in populations) {
    # prepare input (age-year matrices)
    Dxt <- Dx %>% 
      filter(Population == pop) %>%
      mutate(Population = NULL) %>%
      pivot_wider(names_from = "Year", values_from = "Deaths") %>%
      column_to_rownames(var = "Age")
    
    Ext <- Ex %>% 
      filter(Population == pop) %>%
      mutate(Population = NULL) %>%
      pivot_wider(names_from = "Year", values_from = "Exposure") %>%
      column_to_rownames(var = "Age")
    
    # initial fit with StMoMo (to get model structure)
    # warning handler
    handler <- function(w) {
      if(w$message == "non-integer #successes in a binomial glm!" || grepl("non-integer x", w$message)) {
        invokeRestart("muffleWarning")
      }
    }
    models[[pop]] <- withCallingHandlers(fit(cbd(link = "log"), Dxt = Dxt, Ext = Ext, ages = unique(mx$Age), 
                                             years = unique(mx$Year), verbose = FALSE),
                                         warning = handler)
    if (fitMethod == "lm") {
      mxPop <- mx %>% 
        filter(Population == pop) %>%
        mutate(Age_centr = Age - mean(Age))
      
      coefficients <- sapply(unique(mxPop$Year), function(year) coef(lm(log(Rate) ~ Age_centr, mxPop %>% filter(Year == year))))
      models[[pop]]$kt[1, ] <- coefficients[1, ]
      models[[pop]]$kt[2, ] <- coefficients[2, ]
    }
    
    kt1[pop, ] <- models[[pop]]$kt[1, ]
    kt2[pop, ] <- models[[pop]]$kt[2, ]
  }
  
  models$kt1 <- kt1    # make period effect parameters available in the model list
  models$kt2 <- kt2
  class(models) <- "CBD"
  models$fitMethod <- fitMethod
  
  return(models)
}


calcCentralisedLogDeathRates <- function(mx) {
  #' calculates age-year table mx of centralised log death rates and mean matrix ax of log death rates
  #'
  #' @param mx death rates, data.frame with 4 columns: Population, Year, Age, Rate
  #' 
  #' @return list of average and centered log death rates
  
  mx_log <- spread(mx, Year, Rate) %>%
    mutate_at(vars(-Population, -Age), log)
  
  ax <- mx_log %>% 
    dplyr::select(-c("Population", "Age")) %>%
    rowMeans()
  
  mx_log_centralised <- mx_log %>%
    mutate_at(vars(-Population, -Age), function(x) x - ax)
  
  return(list(ax = matrix(ax, nrow = n_distinct(mx$Age), ncol = n_distinct(mx$Population), 
                          dimnames = list(unique(mx$Age), unique(mx$Population))), 
              mx_log_centralised = mx_log_centralised))
}


fitTimeSeries <- function(kt, arimaOrder, include.constant = TRUE) {
  #' Takes matrix of kappa parameters and fits row-wise ARIMA models with drift 
  #' via CSS-ML or CSS if CSS-ML is not possible
  #'
  #' @param kt kappa parameter, period effect as matrix (one row per population)
  #' @param arimaOrder (p,d,q) parameters to use for ARIMA model
  #' @param include.constant whether to include a drift term
  #' 
  #' @return list of fitted ARIMA models
  
  models <- list()
  
  for (i in 1:nrow(kt)) {
    if (sum(abs(kt[i, ])) == 0) {    # if kt identically equals zero, the ARIMA model should reproduce this exactly
      models[[i]] <- Arima(0, include.constant = FALSE)
    } else {
      tryCatch(models[[i]] <- Arima(as.numeric(kt[i, ]), order = arimaOrder, include.constant = include.constant,
                                    method = "CSS-ML"),
               error = function(e) {
                 tryCatch(models[[i]] <<- Arima(as.numeric(kt[i, ]), order = arimaOrder, 
                                                include.constant = include.constant, method = "CSS"),
                          error = function(e) {
                            message(paste("fitTimeSeries: ARIMA model of order", 
                                          paste(arimaOrder, collapse = ", "), 
                                          "could not be fit for population", rownames(kt)[i]))
                            models[[i]] <<- list()
                          }
                 )
               })
    }
    
    # keep arima order in output
    models[[i]]$arimaOrder <- arimaOrder
    models[[i]]$p <- arimaOrder[1]
    models[[i]]$d <- arimaOrder[2]
    models[[i]]$q <- arimaOrder[3]
  }
  
  names(models) <- rownames(kt)
  class(models) <- "ARIMA"
  
  return(models)
}


forecast.LC <- function(model,
                        mx, 
                        evaluationYears, 
                        intervalProbability = 0) {
  #' Forecasts from an ILC using an RWD and adds name column to result
  #'
  #' @param model An ILC model object (list containing ax, bx, kt)
  #' @param mx data.frame with columns Population, Age, Year
  #' @param evaluationYears Vector of years for evaluation
  #' @param intervalProbability Integer between 0 and 100; percentage for prediction interval. If == 0, no PIs are calculated
  #' 
  #' @return data.frame with columns Population, Age, Year, (predicted) Rate, Name and potentially LowerBound, UpperBound
  
  if (length(evaluationYears) > 1 && unique(diff(evaluationYears)) != 1) {
    stop("forecast.LC: evaluationYears have to be a sequence with step size 1.")
  }
  
  # the jump-off year of the model might not be directly one year before min(evaluationYears), so we might need to evaluate the years in between as well
  realEvaluationYears <- (max(as.numeric(colnames(model$kt))) + 1):max(evaluationYears)
  
  mxForecast <- data.frame(Population = rep(unique(mx$Population), length(realEvaluationYears) * n_distinct(mx$Age)), 
                           Year = rep(rep(realEvaluationYears, n_distinct(mx$Age)), each = n_distinct(mx$Population)), 
                           Age = rep(unique(mx$Age), each = n_distinct(mx$Population) * length(realEvaluationYears)),
                           Rate = NA,
                           stringsAsFactors = FALSE) 
  
  predictionDf <- createDfWithPredictedRate(mx = mxForecast, ax = model$ax, bx = model$bx, 
                                            tsModels = model$tsModel, 
                                            nSteps = max(realEvaluationYears) - min(realEvaluationYears) + 1, 
                                            intervalProbability = intervalProbability) %>%
    filter(Year %in% evaluationYears) %>%
    select(one_of(c("Population", "Year", "Age", "PredictedRate", 
                    if (intervalProbability > 0) c("LowerBound", "UpperBound") else c()))) %>% 
    rename(DeathRate = "PredictedRate") %>%
    mutate(Name = model$name)
  
  return(predictionDf)
}


forecast.CBD <- function(models, 
                         mx, 
                         evaluationYears, 
                         intervalProbability = 0) {
  #' Forecasts from an ICBD using a multivariate RWD and adds name column to result
  #'
  #' @param models A CBD model object (list containing StMoMo CBD model objects, one per population)
  #' @param mx only exists for consistency; no value necessary
  #' @param evaluationYears Vector of years for evaluation
  #' @param intervalProbability Integer between 0 and 100; percentage for prediction interval. If == 0, no PIs are calculated
  #' 
  #' @return data.frame with columns Population, Age, Year, (predicted) Rate, Name and potentially LowerBound, UpperBound
  
  if (length(evaluationYears) > 1 && unique(diff(evaluationYears)) != 1) {
    errorMessage <- "forecast.CBD: evaluationYears have to be a sequence with step size 1."
    stop(errorMessage)
  }
  
  populations <- setdiff(names(models), c("fitMethod", "name", "kt1", "kt2"))
  
  # the jump-off year of the model might not be directly one year before min(evaluationYears), so we might need to evaluate the years in between as well
  realEvaluationYears <- (max(as.numeric(colnames(models[[populations[1]]]$kt))) + 1):max(evaluationYears)
  
  predictionDf <- data.frame() 
  
  # call CBD forecast method from StMoMo
  for (pop in populations) {
    popForecast <- (forecast(models[[pop]], h = length(realEvaluationYears)))[["rates"]] %>%
      matrix(ncol = length(realEvaluationYears), dimnames = list(NULL, realEvaluationYears)) %>%
      as.data.frame() %>%
      mutate(Age = models[[pop]]$ages, Population = pop) %>%
      pivot_longer(cols = -one_of(c("Age", "Population")),
                   names_to = "Year", values_to = "DeathRate")
    
    # uncertainty (use simulation --> parameter uncertainty is ignored here!)
    if (intervalProbability > 0) {
      simulationResult <- (simulate(models[[pop]], h = max(2, length(realEvaluationYears)),  
                                    seed = 42))[["rates"]][, 1:length(realEvaluationYears), , drop = FALSE]          # this is a workaround bc h = 1 gives a bug of simulate.StMoMoFit
      
      for (bound in c("LowerBound", "UpperBound")) {
        if (bound == "LowerBound") {
          probs <- (1 - intervalProbability / 100) / 2
        } else {
          probs <- (1 + intervalProbability / 100) / 2
        }
        
        popSimulation <- apply(simulationResult, 1:2, quantile, probs = probs, type = 2) %>%
          as.data.frame() %>%
          mutate(Age = as.numeric(rownames(.)), Population = pop) %>%
          pivot_longer(cols = -one_of(c("Age", "Population")),
                       names_to = "Year", values_to = bound)
        popForecast <- popForecast %>% 
          left_join(popSimulation, by = c("Population", "Year", "Age"))
      }
    }
    
    predictionDf <- predictionDf %>%
      bind_rows(popForecast)
  }
  
  predictionDf <- predictionDf %>%
    mutate(Year = as.numeric(Year)) %>%
    filter(Year %in% evaluationYears) %>%
    select(one_of(c("Population", "Year", "Age", "DeathRate", 
                    if (intervalProbability > 0) c("LowerBound", "UpperBound") else c()))) %>% 
    mutate(Name = models$name)
  
  rownames(predictionDf) <- NULL
  
  return(predictionDf)
}



createDfWithPredictedRate <- function(mx, 
                                      ax, 
                                      bx, 
                                      tsModels, 
                                      nSteps = 1,
                                      intervalProbability = 95) {
  #' Takes LC model parameters and fitted arima models and creates data.frame with predicted mx values, prediction interval bounds
  #'
  #' @param mx data.frame with Population, Year, Age, Rate (= ground truth or NA if ground truth unknown)
  #' @param ax alpha parameter, mean log rates as data.frame (one column per population)
  #' @param bx beta parameter, age effect as data.frame (one column per population)
  #' @param tsModel named list of time series (ARIMA) models, one per population 
  #' @param nSteps number of values to forecast
  #' @param intervalProbability probability level for prediction intervals (in %)  
  #' 
  #' @return mx data.frame with added predicted rates and lower/upper bounds (Population, Year, Age, Rate, PredictedRate, LowerBound, UpperBound) 
  
  # no prediction is made
  if (nSteps <= 0) {
    return(data.frame(Rate = NA, PredictedRate = NA, LowerBound = NA, UpperBound = NA, Residual = NA, RateInInterval = NA))
  }
  
  # get aggregated values out of ARIMA model as a list of matrices (colnames: predicted years, rownames: populations)
  populationList <- lapply(1:length(tsModels),
                           function(i) {
                             forecastARIMA(tsModels[[i]], 
                                           nSteps = nSteps, 
                                           intervalProbability = intervalProbability)
                           })
  
  ktAggregated <- list()
  
  for (val in names(populationList[[1]])) {
    ktAggregated[[val]] <- do.call(rbind, lapply(populationList, function(elem) elem[[val]]))
    
    colnames(ktAggregated[[val]]) <- (sort(unique(mx$Year)))[1:ncol(ktAggregated[[val]])]
    rownames(ktAggregated[[val]]) <- names(tsModels)
  }
  
  # get mx with added columns
  mx <- addModelEvaluationToDf(mx, ax, bx, ktAggregated)
  
  return(mx)
}


forecastARIMA <- function(arimaModel, 
                          nSteps, 
                          intervalProbability = 95) {
  #' Wrapper. Takes ARIMA model and returns point forecasts along with prediction intervals as named list
  #'
  #' @param arimaModel model returned by Arima() function
  #' @param nSteps number of values to forecast
  #' @param intervalProbability probability level for prediction intervals (in %) 
  #' 
  #' @return list of time series forecasts including lower and upper bounds
  
  result <- list()
  
  forecastResult <- suppressWarnings(as.data.frame(forecast(arimaModel, h = nSteps, level = intervalProbability)))
  
  result[["Point Forecast"]] <- forecastResult[, "Point Forecast"]
  result[["Variance"]] <- arimaModel$sigma2
  result[["ArimaOrder"]] <- arimaModel$arimaOrder
  result[[paste("Lo", intervalProbability)]] <- forecastResult[, paste("Lo", intervalProbability)]
  result[[paste("Hi", intervalProbability)]] <- forecastResult[, paste("Hi", intervalProbability)]
  
  return(result)
}


addModelEvaluationToDf <- function(mx, 
                                   ax, 
                                   bx, 
                                   kt) {
  #' Adds fitted/predicted rates and possibly prediction interval bounds as new columns to mx
  #'
  #' @param mx data.frame with Population, Year, Age, Rate (= ground truth or NA if ground truth unknown)
  #' @param ax alpha parameter, mean log rates as data.frame (one column per population)
  #' @param bx beta parameter, age effect as data.frame (one column per population)
  #' @param kt kappa parameter, period effect as data.frame (one row per population) for fitted rates 
  #'           or as list of matrices (with time series aggregate values like mean, lower/upper bound) for predicted rates
  #' 
  #' @return mx data.frame with added predicted/fitted rates and possibly lower/upper bounds 
  
  ## preparations
  ages <- sort(unique(mx$Age))
  years <- sort(unique(mx$Year))
  populations <- sort((unique(as.character(mx$Population))))
  
  mx <- mx %>% arrange(Population, Year, Age)
  mx$PredictedRate <- as.numeric(NA)
  mx$UpperBound <- mx$LowerBound <- as.numeric(NA)
  
  ## evaluate model to obtain fitted/predicted rates
  for (i in 1:length(populations)) {
    # determine prediction interval probability
    intervalProbability <- as.numeric(readr::parse_number(names(kt)[substr(names(kt), 1, 2) == "Lo"]))
    
    # find prediction interval boundaries and point forecast by plugging in values from tsModel
    aggregateValues <- c("Point Forecast", paste("Lo", intervalProbability), paste("Hi", intervalProbability))
    
    populationDf <- list()
    
    for (val in aggregateValues) {
      populationDf <- c(populationDf,
                        gather(as.data.frame(evaluateLC(ax[, i, drop = FALSE], 
                                                        bx[, populations[i], drop = FALSE], 
                                                        matrix(kt[[val]][populations[i], ], nrow = length(populations[i])))), 
                               dummyColumn, !!val) %>%
                          dplyr::select(-dummyColumn))
    }
    
    populationDf <- as.data.frame(populationDf) 
    colnames(populationDf) <- c("PredictedRate", "LowerBound", "UpperBound")
    
    # lower and upper bound might not be ordered correctly because beta might be negative -> reorder
    populationDf[populationDf$LowerBound > populationDf$UpperBound, c("LowerBound", "UpperBound")] <- 
      populationDf[populationDf$LowerBound > populationDf$UpperBound, c("UpperBound", "LowerBound")]
    
    # write fitted values to data.frame 
    mx[((i-1) * nrow(populationDf) + 1):(i * nrow(populationDf)), 
       intersect(colnames(populationDf), c("PredictedRate", "LowerBound", "UpperBound"))] <- populationDf 
  }
  
  return(mx)
}


evaluateLC <- function(ax, 
                       bx, 
                       kt) {
  #' Takes LC-type model parameters and calculates mx evaluations (format: matrix with rows = ages, columns = years)
  #' Can be used for in-sample evaluation and for prediction
  #'
  #' @param ax alpha parameter, base mortality as column vector
  #' @param bx beta parameter, age effect as column vector (multiple age effects as columns of a matrix)
  #' @param kt kappa parameter, period effect as row vector (multiple period effects as rows of a matrix)
  #' 
  #' @return matrix of LC predicted rates
  #' @export
  #'
  #' @examples
  
  bx <- as.matrix(bx)
  kt <- as.matrix(kt)
  
  # every column in bx and corresponding row in kt are corresponding age-period effects. All period effects are added up.
  rates <- Reduce('+', lapply(1:ncol(bx), function(i) bx[, i, drop = FALSE] %*% kt[i, , drop = FALSE]))
  
  # add base level and exp
  rates <- exp(as.matrix(ax) %*% matrix(1, nrow = 1, ncol = ncol(kt)) + rates)
  
  return(rates)
}