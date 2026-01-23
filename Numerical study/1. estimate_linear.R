setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./Simulation/data_sim_linear.R")
library("data.table")
library("Matrix")
library("dplyr")
library("xgboost")
library("survival")
options(na.action='na.pass')
options(scipen = 999)

##Initialize parameters
h <- 1/100
tstop <- 1
cens <- 0.2
grid.times <- seq(from=0,to=tstop,by=h)


###Simulate
seed <- 1
nsim <- 100000
nFit <- c(nsim,nsim*0.1,nsim*0.01,nsim*0.001)

set.seed(seed)
X1 <- rbinom(nsim,1,0.5)
X2 <- rbinom(nsim,1,0.5)
X3 <- rnorm(nsim,0.5,0.5)
cov0 <- data.frame(X1=X1,X2=X2,X3=X3)

dataRaw <- simulator_n(0, 0, cov0)

data <- dataRaw %>% 
  dplyr::filter(!is.na(id)) %>% 
  dplyr::filter(yn==0) %>% 
  dplyr::mutate(ynNext = ifelse(is.na(ynNext),0,ynNext),
                tnNext = ifelse(is.na(tnNext),tstop,tnNext)) %>%
  dplyr::mutate(ynNext = ifelse(ynNext==2,0,ynNext))

###Dynamic Landmarking with uniform sampling
dataStack <- tibble(
  id = character(),
  yn = numeric(),
  ynNext = numeric(),
  tn = numeric(),
  tnNext = numeric(),
  X1 = numeric(),
  X2 = numeric(),
  X3 = numeric(),
  S = numeric()
)

set.seed(seed)
nStack <- 10
for(i in 1:nStack){
  datai <- data %>%
    group_by(id) %>% 
    mutate(S = runif(1,0,tstop)) %>%
    dplyr::mutate(tn=pmax(tn,S), tnNext=pmax(tnNext,S)) %>% 
    dplyr::filter(!tn==tnNext) %>%
    arrange(id,tn) %>%
    mutate(X1=X1[1],X2=X2[1],X3=X3[1]) %>%
    ungroup() %>% 
    mutate(id=paste0(id,"_",i))
  
  dataStack <- union(dataStack, datai)
}

###Construct O/E data
oeStack <- survSplit(
  formula = Surv(tn, tnNext, ynNext) ~ ., 
  data = dataStack, 
  cut = grid.times,
  start = "start",
  end   = "stop"
) 

oe <- survSplit(
  formula = Surv(tn, tnNext, ynNext) ~ ., 
  data = data, 
  cut = grid.times,
  start = "start",
  end   = "stop"
)

#create folds for cross validation - ensure folds are chosen at the id-level (important when nStack > 1)
xgbFolds <- function(df){
  df$id <- sub("_.*", "", df$id)
  group_id <- df$id 
  groups <- unique(group_id)
  K <- 5
  group_folds <- sample(rep(1:K, length.out = length(groups)))
  fold_id <- group_folds[match(group_id, groups)]
  folds <- lapply(1:K, function(k) {
    which(fold_id == k)
  })
  return(folds)
}

xgbInput <- function(df){
  df_input <- df %>%
    dplyr::mutate(exposure=stop-start) %>%
    dplyr::select(-c("id","yn","stop")) %>%
    dplyr::rename("occ"="ynNext", "time"="start")
  
  return(df_input)
}

nQfilter <- function(n,Q,df){
  dfFilter <- df %>% dplyr::filter(
    as.numeric(sub("_.*", "", id)) <= n,
    as.numeric(sub(".*_", "", id)) <= Q
  )
  return(dfFilter)
}

nfilter <- function(n,df){
  dfFilter <- df %>% dplyr::filter(id <= n)
  return(dfFilter)
}



##Estimate: Cox Naive

coxNaiveModels <- list()
for (n in nFit){
  
  startTime <- Sys.time()
  coxTrain <- coxph(Surv(tn, tnNext, ynNext) ~ X1 + X2 + X3, 
                    data = nfilter(n,data),
                    model = TRUE)
  trainTime <- Sys.time() - startTime
  
  modelName <- paste0("coxNaive_", n)
  coxNaiveModels[[modelName]] <- list(model = coxTrain, trainTime = trainTime)
}

saveRDS(coxNaiveModels, file = "Results/coxNaiveModels_linear.rds")


##Estimate: Stacked landmark Cox

coxStackModels <- list()
for (q in c(1,2,5,10)){
  for (n in nFit){
    
    startTime <- Sys.time()
    coxTrain <- coxph(Surv(tn, tnNext, ynNext) ~ X1 + X2 + X3 + S, 
                      data = nQfilter(n,q,dataStack),
                      model = TRUE)
    trainTime <- Sys.time() - startTime
    
    modelName <- paste0("coxStack_", n, "_", q)
    coxStackModels[[modelName]] <- list(model = coxTrain, trainTime = trainTime)
    
  }
}  

saveRDS(coxStackModels, file = "Results/coxStackModels_linear.rds")
  
##Estimate: XGBoost Naive

xgbNaiveParam <- list(
  objective  = "count:poisson", 
  eval_metric = "poisson-nloglik",
  eta = 0.1, 
  max_depth = 1, 
  min_child_weight = 20, 
  subsample = 0.9,
  colsample_bytree = 0.7
)

xgbNaiveModels <- list()
for (n in nFit){
  data_input <- nfilter(n,oe)
  folds <- xgbFolds(data_input)
  data_input <- data_input %>% xgbInput()
  X <- model.matrix(occ ~. - exposure, data = data_input)
  xgtrain = xgb.DMatrix(X, label = data_input$occ)
  F0 <- sum(data_input$occ)/sum(data_input$exposure)
  setinfo(xgtrain, "base_margin", log(data_input$exposure)+log(F0))
  
  set.seed(seed)
  startTime <- Sys.time()
  xgbCv <- xgb.cv(params=xgbNaiveParam,
                  data=xgtrain,
                  folds=folds,
                  nrounds=2000, 
                  early_stopping_rounds = 5)
  
  xgbTrain <- xgb.train(
    params  = xgbNaiveParam,
    data    = xgtrain,
    nrounds = xgbCv$best_iteration
  )
  trainTime <- Sys.time() - startTime
  
  modelName <- paste0("xgbNaive_", n)
  xgbNaiveModels[[modelName]] <- list(model = xgbTrain, F0 = F0, trainTime = trainTime)
}

saveRDS(xgbNaiveModels, file = "Results/xgbNaiveModels_linear.rds")

##Estimate: Stacked landmark XGBoost

xgbStackParam <- list(
  objective  = "count:poisson", 
  eval_metric = "poisson-nloglik",
  eta = 0.1, 
  max_depth = 1, 
  min_child_weight = 20, 
  subsample = 0.9,
  colsample_bytree = 0.7
)

xgbStackModels <- list()
for (q in c(1,2,5,10)){
  for (n in nFit){
        data_input <- nQfilter(n,q,oeStack)
        folds <- xgbFolds(data_input)
        data_input <- data_input %>% xgbInput()
        X <- model.matrix(occ ~. - exposure, data = data_input)
        xgtrain = xgb.DMatrix(X, label = data_input$occ)
        F0 <- sum(data_input$occ)/sum(data_input$exposure)
        setinfo(xgtrain, "base_margin", log(data_input$exposure)+log(F0))
        
        set.seed(seed)
        startTime <- Sys.time()
        xgbCv <- xgb.cv(params=xgbStackParam,
                        data=xgtrain,
                        folds=folds,
                        nrounds=2000, 
                        early_stopping_rounds = 5)
        
        xgbTrain <- xgb.train(
          params  = xgbStackParam,
          data    = xgtrain,
          nrounds = xgbCv$best_iteration
        )
        trainTime <- Sys.time() - startTime
        
        modelName <- paste0("xgbStack_", n, "_", q)
        xgbStackModels[[modelName]] <- list(model = xgbTrain, F0 = F0, trainTime = trainTime)
        
        gc()
  }
}

saveRDS(xgbStackModels, file = "Results/xgbStackModels_linear.rds")






