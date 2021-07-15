#!/usr/local/package/r/3.6.0/bin/Rscript --vanilla
setwd("/archive/data/hgc1074/yingchen/Hiseq_2rounds")
tpm_active_clean <- readRDS("tpm_active.rds")
tpm_active_myiPS <- t(tpm_active_clean)
gtex_active <- readRDS("gtex_ips_active.rds")
gtex_anno_tissue <- readRDS("gtex_ips_anno.rds")
gtex_ips_tpm <- merge(gtex_tpm,ips_tpm,by=0,all=TRUE)



