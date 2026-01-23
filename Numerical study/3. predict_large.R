setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library("data.table")
library("Matrix")
library("dplyr")
library("xgboost")
library("survival")

options(scipen = 999)


dataSurv <- readRDS("Results/dataSurv_large.rds")

tstop <- 1
h <- 1/100

coxNaive <- readRDS("Results/coxNaiveModels_large.rds")
coxStack <- readRDS("Results/coxStackModels_large.rds")
xgbNaive <- readRDS("Results/xgbNaiveModels_large.rds")
xgbStack <- readRDS("Results/xgbStackModels_large.rds")

dataSurvPred <- dataSurv

for (modelName in names(coxNaive)) {
  
  fit <- coxNaive[[modelName]]$model
  
  base_surv <- survfit(fit)
  S0_fun <- stepfun(base_surv$time, c(1, base_surv$surv))
  
  risk_scores <- predict(fit, newdata = dataSurv, type = "risk")
  
  S0_at_tstop <- S0_fun(tstop)
  S0_at_entry <- S0_fun(dataSurv$S)
  
  surv_prob_tstop <- S0_at_tstop ^ risk_scores
  surv_prob_entry <- S0_at_entry ^ risk_scores
  
  dataSurvPred[[modelName]] <- surv_prob_tstop / surv_prob_entry
}

for (modelName in names(coxStack)) {
  
  fit <- coxStack[[modelName]]$model
  
  base_surv <- survfit(fit)
  S0_fun <- stepfun(base_surv$time, c(1, base_surv$surv))
  
  risk_scores <- predict(fit, newdata = dataSurv, type = "risk")
  
  S0_at_tstop <- S0_fun(tstop)
  S0_at_entry <- S0_fun(dataSurv$S)
  
  surv_prob_tstop <- S0_at_tstop ^ risk_scores
  surv_prob_entry <- S0_at_entry ^ risk_scores
  
  dataSurvPred[[modelName]] <- surv_prob_tstop / surv_prob_entry
}


for (i in 1:nrow(dataSurv)){
  print(i)
  datai <- dataSurv[i,]

  #xgb common precompute
  grid.pred <- seq(from=datai$S, to=tstop, by=h)
  datai <- bind_rows(replicate(length(grid.pred), datai, simplify = FALSE))
  datai$time <- grid.pred
  datai$exposure <- h
      
  #xgb naive
  for (modelName in names(xgbNaive)){
    fit <- xgbNaive[[modelName]]$model
    F0  <- xgbNaive[[modelName]]$F0
    X <- model.matrix(~ . - id - surv - exposure - S, data = datai)
    xgbDatai = xgb.DMatrix(X)
    setinfo(xgbDatai, "base_margin", log(datai$exposure)+log(F0))
    haz <- predict(object=fit, newdata=xgbDatai, type="response")
    surv <- exp(-sum(haz))
    
    dataSurvPred[i,modelName] <- surv
  }
  
  #xgb stacked landmarked
  for (modelName in names(xgbStack)){
    fit <- xgbStack[[modelName]]$model
    F0  <- xgbStack[[modelName]]$F0
    X <- model.matrix(~ . - id - surv - exposure, data = datai)
    xgbDatai = xgb.DMatrix(X)
    setinfo(xgbDatai, "base_margin", log(datai$exposure)+log(F0))
    haz <- predict(object=fit, newdata=xgbDatai, type="response")
    surv <- exp(-sum(haz))
    
    dataSurvPred[i,modelName] <- surv
  }
  
}

saveRDS(dataSurvPred, file = "Results/dataSurvPred_large.rds")