#!/usr/local/package/r/3.6.0/bin/Rscript --vanilla

#############function########


cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}


setwd("/archive/data/hgc1074/yingchen/myiPS")
load("gtex_prediciton_model_L2.RData")
gtex_weight_l2 <-as.data.frame(m_LR$W)
saveRDS(gtex_weight_l2,"gtex_weight_l2.rds")

gtex_weight_l2 <- gtex_weight_l2[,-ncol(gtex_weight_l2)]
myips <-readRDS("tpm_myiPS.rds")
#exp of the coeff of the weight 
gtex_weight_l2 <- exp(gtex_weight_l2)

gtex_weight_l2<-t(gtex_weight_l2)

#check if the rownames are exactly same
if(identical(rownames(gtex_weight_l2),rownames(myips))){
message("rownames are paired")
} else{
message("rownames are different!")
}
myips_tissue <-  data.frame()

###################### not exp the weight 
for(i in colnames(gtex_weight_l2)){
 index <- match(i, colnames(gtex_weight_l2))
 tissue_allips <- myips * rep(gtex_weight_l2[,index],ncol(myips))
 colnames(tissue_allips) <- paste(as.character(i),colnames(tissue_allips),sep="_")
 myips_tissue <- cbind.fill(myips_tissue,tissue_allips)
}

saveRDS(myips_tissue,"myips_tissue_1380.rds")

########################### exp the weight 
for(i in colnames(gtex_weight_l2)){
 index <- match(i, colnames(gtex_weight_l2))
 tissue_allips <- myips * rep(gtex_weight_l2[,index],ncol(myips))
 colnames(tissue_allips) <- paste(as.character(i),colnames(tissue_allips),sep="_")
 myips_tissue <- cbind.fill(myips_tissue,tissue_allips)
}

saveRDS(myips_tissue,"myips_exptissue_1380.rds")

