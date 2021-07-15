#!/usr/local/package/r/3.6.0/bin/Rscript --vanilla
setwd("/archive/data/hgc1074/yingchen/Hiseq_2rounds")
source("./gtex_ips.ELNET.function.r")
if (!requireNamespace("doMC", quietly = TRUE)) {
  install.packages("doMC")
  library(doMC) } else {library(doMC)}

if (!requireNamespace("glmnet", quietly = TRUE)) { 
  install.packages("glmnet", dependencies = T)
  library(glmnet) 
} else {library(glmnet)}

if (!requireNamespace("c060", quietly = TRUE)) { 
  install.packages("c060")
  library(c060) } else {library(c060)}


run_nestedcv_GLMNET()


