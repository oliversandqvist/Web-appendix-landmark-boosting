#install.packages("HQM")
library("HQM")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("data.table")
library("data.table")
#install.packages("pdp")
library("pdp")
#install.packages("xgboost")
library("xgboost")
#install.packages("ggplot2")
library("ggplot2")
#install.packages("SHAPforxgboost")
library("SHAPforxgboost")

data(pbc2)

###Summary statistics
statusVec <- pbc2 %>% 
  group_by(id) %>% 
  dplyr::slice(c(1)) %>% 
  ungroup() %>%
  dplyr::select(status) %>% 
  unlist() 
ftable(statusVec)

png("Plots/visitingHistogram.png", width = 6, height = 4, units = 'in', res = 300)
pbc2 %>% 
  group_by(id) %>% 
  dplyr::mutate(normalizedYear = year/years) %>% 
  ungroup() %>% 
  dplyr::filter(normalizedYear > 0) %>%
  dplyr::select(normalizedYear) %>%
  unlist() %>%
  hist(, main="Histogram of normalized visiting times", xlab="Normalized visiting time")
dev.off()

pbc2 %>% dplyr::filter(year > 0) %>% nrow() / length(unique(pbc2$id)) #=5.233974

pbc2 %>% group_by(id) %>% summarise(years=mean(years)) %>% dplyr::select(years) %>% unlist() %>% as.numeric() %>% mean() #=6.411239 

pbc2 %>% dplyr::filter(year > 0) %>% nrow()
pbc2 %>% dplyr::filter(year > 0) %>% summarise(across(everything(), ~sum(is.na(.))))

###Create landmarked data
pbc2_landmark <- pbc2 %>% dplyr::mutate(sUnif=0,S=0,yearNext=0) %>% head(,n=0) #initialize empty dataset

nStack <- 10 
seed <- 1
T <- max(pbc2$years)
set.seed(seed)

for(i in 1:nStack){
datai <- pbc2 %>%
  dplyr::filter(year > 0) %>%
  group_by(id) %>% 
  arrange(id,year) %>%
  mutate(sUnif = runif(1,0,T), S = sample(year,1), yearNext=lead(year)) %>%
  mutate(yearNext=ifelse(is.na(yearNext),years,yearNext)) %>%
  dplyr::filter(year >= S) %>%
  mutate(ascites=ascites[1],
         hepatomegaly=hepatomegaly[1],
         spiders=spiders[1],
         edema=edema[1],
         serBilir=serBilir[1],
         serChol=serChol[1],
         albumin=albumin[1],
         alkaline=alkaline[1],
         SGOT=SGOT[1],
         platelets=platelets[1],
         prothrombin=prothrombin[1],
         histologic=histologic[1]) %>%
  mutate(age=age+S) %>%
  ungroup() %>% 
  dplyr::filter(sUnif <= years) %>% 
  mutate(id=paste0(id,"_",i)) #keep ids seperate to make exposure and occurrence transform simpler

pbc2_landmark <- union(pbc2_landmark,datai)
}

h <- 1/12
tmin <- 0
tmax <- 15
grid.times <- seq(from=tmin,to=tmax,by=h)
pbc2_landmark_OE <- data.table(pbc2_landmark)
pbc2_landmark_OE <- do.call("rbind", lapply(1:length(grid.times), function(t) {
  data_t <- copy(pbc2_landmark_OE)[, grid.left:=grid.times[t]] 
}))
pbc2_landmark_OE <- pbc2_landmark_OE[grid.left+h >= year]
pbc2_landmark_OE <- pbc2_landmark_OE[grid.left<yearNext]
pbc2_landmark_OE$grid.right <- pbc2_landmark_OE$grid.left+h
pbc2_landmark_OE <- pbc2_landmark_OE[yearNext<=grid.right, grid.right:=yearNext] 
pbc2_landmark_OE <- pbc2_landmark_OE[, expo.time:=grid.right-c(0, grid.right[-.N]), by="id"]
pbc2_landmark_OE <- pbc2_landmark_OE[grid.left<=S, expo.time:=expo.time-S]
pbc2_landmark_OE <- pbc2_landmark_OE[grid.left<=year, grid.left:=year]
pbc2_landmark_OE$occ <- 0
pbc2_landmark_OE <- pbc2_landmark_OE[years<=grid.right, occ:=ifelse(status=="alive",0,1)]

pbc2_landmark_OE_input <- pbc2_landmark_OE %>% 
  mutate(time=grid.left, 
         age=age+time) 

pbc2_landmark_OE_input$id <- sub("_.*", "", pbc2_landmark_OE_input$id)

#create folds for cross validation - ensure folds are chosen at the id-level (important when nStack > 1)
group_id <- pbc2_landmark_OE_input$id 
groups <- unique(group_id)
K <- 5
group_folds <- sample(rep(1:K, length.out = length(groups)))
fold_id <- group_folds[match(group_id, groups)]
folds <- lapply(1:K, function(k) {
  which(fold_id == k)
})

pbc2_landmark_OE_input <- pbc2_landmark_OE_input %>% 
  dplyr::select(-c("years","status","year","status2","yearNext","grid.left","grid.right","id","sUnif"))

#reparameterize for easier learning that time since landmarking is an important time axis
pbc2_landmark_OE_input$time <- pbc2_landmark_OE_input$time-pbc2_landmark_OE_input$S


##XGBoost fitting
F0 <- sum(pbc2_landmark_OE_input$occ)/sum(pbc2_landmark_OE_input$expo.time)

options(na.action='na.pass')
X <- model.matrix(occ ~. - expo.time, data = pbc2_landmark_OE_input)

xgtrain = xgb.DMatrix(X, label = pbc2_landmark_OE_input$occ)
setinfo(xgtrain, "base_margin", log(pbc2_landmark_OE_input$expo.time)+log(F0))

xgparam <- list(
  objective  = "count:poisson", 
  eval_metric = "logloss",
  eta = 0.01,
  max_depth = 4,    
  subsample = 0.8,
  lambda= 2,
  colsample_bytree= 0.5
)

set.seed(seed)
xgbCv <- xgb.cv(params=xgparam,
                data=xgtrain,
                folds=folds,
                nrounds=2000, 
                early_stopping_rounds = 25) 
#[638]	train-logloss:0.035914+0.001240	test-logloss:0.039814+0.004650

start.time <- Sys.time()
xgbTrain = xgb.train(
  params  = xgparam,
  data    = xgtrain,
  nrounds = 638
)
end.time <- Sys.time()
print(end.time - start.time)

#sanity check
pred <- predict(xgbTrain, newdata = xgtrain)
sum(pred) #495.7467
sum(pbc2_landmark_OE_input$occ) #523

#importance
importance_matrix <- xgb.importance(colnames(X), model = xgbTrain)
png("Plots/importance.png", width = 6, height = 4, units = 'in', res = 300)
xgb.plot.importance(importance_matrix, rel_to_first = TRUE, xlab = "Relative importance")
dev.off()

#shapley values
shap_val <- shap.prep(xgb_model = xgbTrain, X_train = X)
png("Plots/shapley.png", width = 6, height = 4, units = 'in', res = 300)
shap.plot.summary(shap_val)
dev.off()

#partial dependence plots
df_serBilirPred <- data.frame(pbc2_landmark_OE_input[,"serBilir"], 
                              pbc2_landmark_OE_input[,"expo.time"], 
                              pred)

df_serAlbuminPred <- data.frame(pbc2_landmark_OE_input[,"albumin"], 
                                pbc2_landmark_OE_input[,"expo.time"], 
                                pred)

pdp_serBilir <- pdp::partial(
  object = xgbTrain,
  pred.var = "serBilir",
  train = X,
  plot = FALSE
)
png("Plots/pdp serBilir.png", width = 6, height = 4, units = 'in', res = 300)
plot(pdp_serBilir$serBilir,F0*pdp_serBilir$yhat,type='l', xlab="Level of bilirubin", ylab="future conditional hazard")
rug(df_serBilirPred$serBilir, col = "#00000040")  
dev.off()

pdp_albumin <- pdp::partial(
  object = xgbTrain,
  pred.var = "albumin",
  train = X,
  plot = FALSE
)
png("Plots/pdp albumin.png", width = 6, height = 4, units = 'in', res = 300)
plot(pdp_albumin$albumin,F0*pdp_albumin$yhat,type='l', xlab="Level of albumin", ylab="future conditional hazard")
rug(df_serAlbuminPred$albumin, col = "#00000040")
dev.off()


#marginal plot
df_binned_serBilir <- df_serBilirPred %>%
  mutate(bin = cut(serBilir, breaks = 10)) %>%   
  group_by(bin) %>%
  summarize(mean_pred = sum(pred)/sum(expo.time), 
            mid_x = mean(serBilir))

png("Plots/mplot serBilir.png", width = 6, height = 4, units = 'in', res = 300)
plot(df_binned_serBilir$mid_x, df_binned_serBilir$mean_pred, type="l",
     xlab="Level of bilirubin (binned)", ylab="Average future conditional hazard")
rug(df_serBilirPred$serBilir, col = "#00000040")  
dev.off()

df_binned_serAlbumin <- df_serAlbuminPred %>%
  mutate(bin = cut(albumin, breaks = 10)) %>%  
  group_by(bin) %>%
  summarize(mean_pred = sum(pred)/sum(expo.time), 
            mid_x = mean(albumin))

png("Plots/mplot albumin.png", width = 6, height = 4, units = 'in', res = 300)
plot(df_binned_serAlbumin$mid_x, df_binned_serAlbumin$mean_pred, type="l",
     xlab="Level of albumin (binned)", ylab="Average future conditional hazard")
rug(df_serAlbuminPred$albumin, col = "#00000040")  
dev.off()

#survival function 10 years for selected id's
id1 <- "128_6"
id2 <- "25_2"

tmaxSurv <- 10
grid.times.surv <- seq(from=tmin,to=tmaxSurv,by=h)

dfX1 <- pbc2_landmark_OE %>% dplyr::filter(id==id1 & grid.left==S) %>% 
  mutate(time=0, age=age+S) %>%
  dplyr::select(-c("years","status","year","status2","yearNext","grid.left","grid.right","id","sUnif")) %>% as.data.frame()
dfX1 <- dfX1[rep(1, length(grid.times.surv)), ] 
dfX1$time <- grid.times.surv
dfX1$expo.time <- h
print(dfX1[1,])

X1 <- model.matrix(occ ~. - expo.time, data = dfX1)
xgtrain_X1 = xgb.DMatrix(X1, label = dfX1$occ)
setinfo(xgtrain_X1, "base_margin", log(dfX1$expo.time)+log(F0))
predX1 <- predict(xgbTrain, newdata = xgtrain_X1)
png("Plots/surv id128.png", width = 6, height = 4, units = 'in', res = 300)
plot(grid.times.surv,exp(-cumsum(predX1)), type="l", xlab="time", ylab="survival")
dev.off()

dfX2 <- pbc2_landmark_OE %>% dplyr::filter(id==id2 & grid.left==S) %>% 
  mutate(time=0, age=age+S) %>%
  dplyr::select(-c("years","status","year","status2","yearNext","grid.left","grid.right","id","sUnif")) %>% as.data.frame()
dfX2 <- dfX2[rep(1, length(grid.times.surv)), ] 
dfX2$time <- grid.times.surv
dfX2$expo.time <- h
print(dfX2[1,])


X2 <- model.matrix(occ ~. - expo.time, data = dfX2)
xgtrain_X2 = xgb.DMatrix(X2, label = dfX2$occ)
setinfo(xgtrain_X2, "base_margin", log(dfX2$expo.time)+log(F0))
predX2 <- predict(xgbTrain, newdata = xgtrain_X2)
png("Plots/surv id25.png", width = 6, height = 4, units = 'in', res = 300)
plot(grid.times.surv,exp(-cumsum(predX2)), type="l", xlab="time", ylab="survival")
dev.off()
