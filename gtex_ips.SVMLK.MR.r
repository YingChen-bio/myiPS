#!/usr/local/package/r/3.6.0/bin/Rscript --vanilla
setwd("/archive/data/hgc1074/yingchen/Hiseq_2rounds")
source("./gtex_ips.SVMLK.MR.function.r")

if (!requireNamespace("caret", quietly = TRUE)) {
  install.packages("caret", dependencies = T)
  library(caret)
} else {library(caret)}

if (!requireNamespace("doMC", quietly = TRUE)) {
  install.packages("doMC")
  library(doMC) } else {library(doMC)}

if (!requireNamespace("e1071", quietly = TRUE)) {
  install.packages("e1071")
  library(e1071) } else {library(e1071)}

run_nestedcv_SVM_e1071()

