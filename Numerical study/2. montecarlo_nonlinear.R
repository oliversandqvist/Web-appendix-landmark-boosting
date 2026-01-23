setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./Simulation/data_sim_nonlinear.R")
library("data.table")
library("Matrix")
library("dplyr")
library("xgboost")
library("survival")
options(na.action='na.pass')
options(scipen = 999)

##Initialize parameters
tstop <- 1
cens <- 0

###Simulate test sample
seed <- 1337
nsim <- 1000

set.seed(seed)
nsimWithBuffer <- nsim*3 #simulate more than nsim to end up with nsim subjects, as some are lost in the landmarking
X1 <- rbinom(nsimWithBuffer,1,0.5)
X2 <- rbinom(nsimWithBuffer,1,0.5)
X3 <- rnorm(nsimWithBuffer,0.5,0.5)
X4 <- rep(0,nsimWithBuffer) #not used for fitting, only for simulation
cov0 <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4)

dataRaw <- simulator_n(0, 0, cov0)

data <- dataRaw %>% 
  dplyr::filter(!is.na(id)) %>% 
  dplyr::filter(yn==0) %>% 
  dplyr::mutate(ynNext = ifelse(is.na(ynNext),0,ynNext),
                tnNext = ifelse(is.na(tnNext),tstop,tnNext)) %>%
  dplyr::mutate(ynNext = ifelse(ynNext==2,0,ynNext))

set.seed(seed)
dataStack <- data %>%
  group_by(id) %>% 
  mutate(S = runif(1,0,tstop),
         X4 = lag(X3, default = 0)) %>%
  dplyr::filter(tn <= S & S < tnNext) %>%
  arrange(id,tn) %>%
  mutate(X1=X1[1],X2=X2[1],X3=X3[1],X4=X4[1]) %>%
  ungroup() %>%
  dplyr::select(-c("id","yn","ynNext","tn","tnNext")) %>%
  head(n=nsim)

###Calculate survival using Monte Carlo
namesDataAll <- c("id",paste0("X", 1:3),"S","surv")
nColAll <- length(namesDataAll)
dataSurv <- setNames(data.frame(matrix(nrow = nsim, ncol = nColAll)), namesDataAll)

set.seed(seed)
nMC <- 100000
for (i in 1:nsim){
  print(i)
  
  X1i <- as.numeric(dataStack[i,"X1"])
  X2i <- as.numeric(dataStack[i,"X2"])
  X3i <- as.numeric(dataStack[i,"X3"])
  X4i <- as.numeric(dataStack[i,"X4"])
  Si <- as.numeric(dataStack[i,"S"])
  covi <- data.frame(X1 = rep(X1i, nMC), 
                     X2 = rep(X2i, nMC),
                     X3 = rep(X3i, nMC),
                     X4 = rep(X4i, nMC))
  dataMC <- simulator_n(0, Si, covi)
  
  dataMC <- dataMC %>% 
    dplyr::filter(yn==0) %>% 
    dplyr::mutate(ynNext = ifelse(is.na(ynNext),0,ynNext))
  
  dataMC <- dataMC %>% 
    group_by(id) %>%
    summarize(ynNextMax=max(ynNext))
  
  survi <- 1-sum(dataMC$ynNextMax)/nMC 
  
  dataSurv[i,] <- c(i,X1i,X2i,X3i,Si,survi)
  
  gc()
}

saveRDS(dataSurv, file = "Results/dataSurv_nonlinear.rds")






