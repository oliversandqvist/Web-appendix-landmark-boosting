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
S <- runif(nsimWithBuffer,0,tstop)
nNoise <- 47
A <- matrix(rnorm(nNoise * nNoise), nrow = nNoise, ncol = nNoise)
Sigma <- A %*% t(A)
L_sigma <- t(chol(Sigma))
Z <- matrix(rnorm(nsimWithBuffer * nNoise), nrow = nsimWithBuffer, ncol = nNoise)
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

set.seed(seed)
dataStack <- data %>%
  group_by(id) %>% 
  mutate(S = runif(1, 0, tstop)) %>%
  dplyr::filter(tn <= S & S < tnNext) %>%
  arrange(id, tn) %>%
  mutate(across(starts_with("X"), ~ .x[1])) %>%
  ungroup() %>%
  dplyr::select(-c("id","yn","ynNext","tn","tnNext")) %>%
  head(n=nsim)


###Calculate survival using Monte Carlo
namesDataAll <- c("id",paste0("X", 1:50),"S","surv")
nColAll <- length(namesDataAll)
dataSurv <- setNames(data.frame(matrix(nrow = nsim, ncol = nColAll)), namesDataAll)

set.seed(seed)
nMC <- 100000
for (i in 1:nsim){
  print(i)
  
  covi <- dataStack[i,1:(3+nNoise)]
  Si <- as.numeric(dataStack[i,"S"])
  dataMC <- simulator_n(0, Si, as.data.frame(covi)[rep(1, nMC), ] )
  
  dataMC <- dataMC %>% 
    dplyr::filter(yn==0) %>% 
    dplyr::mutate(ynNext = ifelse(is.na(ynNext),0,ynNext))
  
  dataMC <- dataMC %>% 
    group_by(id) %>%
    summarize(ynNextMax=max(ynNext))
  
  survi <- 1-sum(dataMC$ynNextMax)/nMC 
  
  dataSurv[i,] <- c(i,as.numeric(covi),Si,survi)
  
  gc()
}

saveRDS(dataSurv, file = "Results/dataSurv_large.rds")






