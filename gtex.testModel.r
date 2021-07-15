#!/usr/local/package/r/3.6.0/bin/Rscript --vanilla


if (!requireNamespace("doMC", quietly = TRUE)) {
  install.packages("doMC")
  library(doMC) } else {library(doMC)}

if (!requireNamespace("e1071", quietly = TRUE)) {
  install.packages("e1071")
  library(e1071) } else {library(e1071)}

if (!requireNamespace("LiblineaR", quietly = TRUE)) {
  install.packages("LiblineaR")
  library(LiblineaR) } else {library(LiblineaR)}

set.seed(999)
setwd("/archive/data/hgc1074/yingchen/Hiseq_2rounds")
tpm_active_clean <- readRDS("tpm_active.rds")
tpm_active_myiPS <- t(tpm_active_clean)
gtex_active <- readRDS("gtex_ips_active.rds")
#only gtex
gtex_active <- gtex_active[,0:17382]


gtex_anno_tissue <- readRDS("gtex_ips_anno.rds")
#only gtex
gtex_anno_tissue <- gtex_anno_tissue[0:17382,]

#tpm_myiPS <- readRDS("tpm_myiPS.rds")

train_data <- t(gtex_active)
label_data <- as.factor(gtex_anno_tissue$Tissue)
#test_data <-t(tpm_myiPS)

#center and scale data

scaled_train <- scale(train_data,center=TRUE,scale=TRUE)
scaled_test <- scale(test_data,attr(scaled_train,"scaled:center"),attr(scaled_train,"scaled:scale"))
scaled_test[is.na(scaled_test)] <- 0


#L2-regularized logistic regression
t=7 #L2-regularized logistic regression
# Tune the cost parameter of a logistic regression according to the Joachim’s heuristics
co=heuristicC(scaled_train)
m_LR=LiblineaR(data=scaled_train,labels=label_data,type=t,cost=co,bias=TRUE,verbose=FALSE)
save.image(file="gtex_prediciton_model_L2.RData")



# Make prediction

#p_LR=predict(m_LR,scaled_test)

#saveRDS(p_LR,"gtex_predicitin_liblinear.rds")


tryTypes=c(7:7)
tryCosts=c(1000,100,10,1,0.1,0.01,0.001)

bestCost=NA
bestAcc=0
bestType=NA
for(ty in tryTypes){
  for(co in tryCosts){
    acc=LiblineaR(data=scaled_train,labels=label_data,type=ty,cost=co,bias=TRUE,cross=10,verbose=FALSE)
    cat("Results for C=",co," : ",acc," accuracy.\n",sep="",file="Liblinear.result.txt",append=TRUE)
    if(acc>bestAcc){
      bestCost=co
      bestAcc=acc
      bestType=ty
    }
  }
}
sink("Liblinear.result.txt",append=TRUE)
cat("Best model type is:",bestType,"\n")
cat("Best cost is:",bestCost,"\n")
cat("Best accuracy is:",bestAcc,"\n")
sink()


m_best=LiblineaR(data=scaled_train,labels=label_data,type=bestType,cost=bestCost,bias=TRUE,verbose=FALSE)

pr=FALSE
if(bestType==0 | bestType==7) pr=TRUE
p_best=predict(m_best,scaled_test,proba=pr,decisionValues=TRUE)
save.image(file='Liblinear_result_L2.RData')


