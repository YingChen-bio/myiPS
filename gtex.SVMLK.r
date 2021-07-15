#!/usr/local/package/r/3.6.0/bin/Rscript --vanilla

setwd("/archive/data/hgc1074/yingchen/myiPS")
source("./gtex.SVMLK.function.r")

if (!requireNamespace("doMC", quietly = TRUE)) {
  install.packages("doMC")
  library(doMC) } else {library(doMC)}

if (!requireNamespace("e1071", quietly = TRUE)) { 
  install.packages("e1071")
  library(e1071) } else {library(e1071)}

if (!requireNamespace("LiblineaR", quietly = TRUE)) {
  install.packages("LiblineaR")
  library(LiblineaR) } else {library(LiblineaR)}

run_nestedcv_SVM_LiblineaR()



