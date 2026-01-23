setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./Simulation/data_sim_large.R")
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
nNoise <- 47
set.seed(seed)
A <- matrix(rnorm(nNoise * nNoise), nrow = nNoise, ncol = nNoise)
Sigma <- A %*% t(A) 
L_sigma <- t(chol(Sigma))
Z <- matrix(rnorm(nsim * nNoise), nrow = nsim, ncol = nNoise)

set.seed(seed)
X1 <- rbinom(nsim,1,0.5)
X2 <- rbinom(nsim,1,0.5)
X3 <- rnorm(nsim,0.5,0.5)
XNoise <- Z %*% t(L_sigma)
cov0 <- as.data.frame(cbind(X1,X2,X3,XNoise))
names(cov0) <- paste0("X", 1:(3+nNoise))

dataRaw <- simulator_n(0, 0, cov0)

data <- dataRaw %>% 
  dplyr::filter(!is.na(id)) %>% 
  dplyr::filter(yn==0) %>% 
  dplyr::mutate(ynNext = ifelse(is.na(ynNext),0,ynNext),
                tnNext = ifelse(is.na(tnNext),tstop,tnNext)) %>%
  dplyr::mutate(ynNext = ifelse(ynNext==2,0,ynNext))

###Dynamic Landmarking with uniform sampling
set.seed(seed)
nStack <- 10
dataStack_list <- lapply(1:nStack, function(i) {
  data %>%
    group_by(id) %>% 
    mutate(S = runif(1, 0, tstop)) %>%
    dplyr::mutate(tn = pmax(tn, S), tnNext = pmax(tnNext, S)) %>% 
    dplyr::filter(!tn==tnNext) %>%
    arrange(id, tn) %>%
    mutate(across(starts_with("X"), ~ .x[1])) %>%
    ungroup() %>% 
    mutate(id = paste0(id, "_", i))
})
dataStack <- bind_rows(dataStack_list)
rm(dataStack_list)

###Construct O/E data
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

xgbFoldsK3 <- function(df){
  df$id <- sub("_.*", "", df$id)
  group_id <- df$id 
  groups <- unique(group_id)
  K <- 3
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
  coxTrain <- coxph(Surv(tn, tnNext, ynNext) ~ ., 
                    data = nfilter(n,data) %>% dplyr::select(-c("id","yn")),
                    model = TRUE)
  trainTime <- Sys.time() - startTime
  
  modelName <- paste0("coxNaive_", n)
  coxNaiveModels[[modelName]] <- list(model = coxTrain, trainTime = trainTime)
}

saveRDS(coxNaiveModels, file = "Results/coxNaiveModels_large.rds")


##Estimate: Stacked landmark Cox

coxStackModels <- list()
for (q in c(1,2,5,10)){
  for (n in nFit){
    
    startTime <- Sys.time()
    
    coxTrain <- coxph(Surv(tn, tnNext, ynNext) ~ ., 
                      data = nQfilter(n,q,dataStack) %>% dplyr::select(-c("id","yn")),
                      model = TRUE)
    trainTime <- Sys.time() - startTime
    
    modelName <- paste0("coxStack_", n, "_", q)
    coxStackModels[[modelName]] <- list(model = coxTrain, trainTime = trainTime)
    
  }
}  

saveRDS(coxStackModels, file = "Results/coxStackModels_large.rds")
  
##Estimate: XGBoost Naive

xgbNaiveParam <- list(
  objective  = "count:poisson", 
  eval_metric = "poisson-nloglik",
  eta = 0.1, 
  max_depth = 1, 
  min_child_weight = 20,  
  subsample = 0.7,
  alpha = 100
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
  
  gc()
  
  xgbTrain <- xgb.train(
    params  = xgbNaiveParam,
    data    = xgtrain,
    nrounds = xgbCv$best_iteration
  )
  trainTime <- Sys.time() - startTime
  
  modelName <- paste0("xgbNaive_", n)
  xgbNaiveModels[[modelName]] <- list(model = xgbTrain, F0 = F0, trainTime = trainTime)
}

saveRDS(xgbNaiveModels, file = "Results/xgbNaiveModels_large.rds")

##Estimate: Stacked landmark XGBoost

rm(oe,cov0)
gc()

xgbStackParam <- list(
  objective  = "count:poisson", 
  eval_metric = "poisson-nloglik",
  eta = 0.1, 
  max_depth = 1, 
  min_child_weight = 100,  
  subsample = 0.7,
  alpha = 100
)

xgbStackModels <- list()
for (q in rev(c(1,2,5,10))){
  for (n in nFit){
    idx_subset <- which(as.numeric(sub("_.*", "", dataStack$id)) <= n & 
                          as.numeric(sub(".*_", "", dataStack$id)) <= q)

    data_input <- survSplit(
      formula = Surv(tn, tnNext, ynNext) ~ ., 
      data = dataStack[idx_subset,], 
      cut = grid.times,
      start = "start",
      end   = "stop"
    )
    if(q == 10 & n == nsim){
      folds <- xgbFoldsK3(data_input) #due to memory limitations
    } else{
      folds <- xgbFolds(data_input)
    }
    data_input <- data_input %>% xgbInput()
    X <- model.matrix(occ ~. - exposure, data = data_input)
    occVec <- data_input$occ
    expoVec <- data_input$exposure
    
    rm(data_input)
    gc()
    
    xgtrain = xgb.DMatrix(X, label = occVec)
    F0 <- sum(occVec)/sum(expoVec)
    setinfo(xgtrain, "base_margin", log(expoVec)+log(F0))
    
    rm(X)
    gc()
    
    set.seed(seed)
    startTime <- Sys.time()
    
    xgbCv <- xgb.cv(params=xgbStackParam,
                    data=xgtrain,
                    folds=folds,
                    nrounds=2000, 
                    early_stopping_rounds = 5)
    best_iteration <- xgbCv$best_iteration
    
    rm(xgbCv)
    gc()
    
    xgbTrain <- xgb.train(
      params  = xgbStackParam,
      data    = xgtrain,
      nrounds = best_iteration
    )
    trainTime <- Sys.time() - startTime
    
    modelName <- paste0("xgbStack_", n, "_", q)
    xgbStackModels[[modelName]] <- list(model = xgbTrain, F0 = F0, trainTime = trainTime)
    
    rm(xgbTrain,xgtrain)
    gc()
  }
}

saveRDS(xgbStackModels, file = "Results/xgbStackModels_large.rds")






