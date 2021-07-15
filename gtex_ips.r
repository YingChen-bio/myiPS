#!/usr/local/package/r/3.6.0/bin/Rscript --vanilla
library("zFPKM")
library("biomaRt")
library("ggplot2")
library("pheatmap")
library("ggfortify")
library("M3C")
library("Rtsne")
setwd("/archive/data/hgc1074/yingchen/Hiseq_2rounds")
source("./gtex_ips.function.r")

if(FALSE){

ips_tpm <- read.table("all.counts.ips_tpm.txt",sep="\t",header=TRUE)[,-c(2:6)]
rownames(ips_tpm) <-as.character(ips_tpm$Geneid)
ips_tpm <- ips_tpm[,-1]
#gtex_tpm <- readRDS("gtex_tpm.rds")
#gtex_tpm <- as.data.frame(gtex_tpm)

#get stable id of ENSG of gtex
#rownames(gtex_tpm) <- make.names(gsub("\\..*", "", rownames(gtex_tpm)),unique=TRUE)

#gtex_tmp_sample <- gtex_tpm[1:100,1:100]
#write.csv(gtex_tmp_sample,"gtex_tmp_sample.csv")
#write.csv(gtex_tpm,"gtex_tpm_noversion.csv")
gtex_tpm <- read.csv("gtex_tpm_noversion.csv")
rownames(gtex_tpm) <-as.character(gtex_tpm$X)
gtex_tpm <- gtex_tpm[,-1]
gtex_ips_tpm <- merge(gtex_tpm,ips_tpm,by=0,all=TRUE)

#write.csv(gtex_ips_tpm[1:100,1:100],"gtex_ips_tpm_sample.csv")

gtex_ips_tpm[is.na(gtex_ips_tpm)] <- 0

saveRDS(gtex_ips_tpm[1:100,17300:17450],"merged_gtex_ips.rds")


gtex_ips_tpm <- readRDS("merged_17838samples.rds")

gtex_ips_tpm$sum <- rowSums(gtex_ips_tpm)
gtex_cut <- subset(gtex_ips_tpm,sum>0)
gtex_clean <- gtex_cut[,-ncol(gtex_cut)]
zTPM <- zFPKM(gtex_clean,assayName ="tpm")
#here
activeGenes <- which(rowMeans(zTPM) > -3)
gtex_active <- gtex_clean[activeGenes,]
saveRDS(gtex_active,"gtex_ips_active.rds")


gtex_active <- readRDS("gtex_ips_active.rds")
gtex_active$var <- apply(gtex_active,1,var)
gtex_active_clean<- gtex_active[,-ncol(gtex_active)]
saveRDS(gtex_active,"gtex_ips_var.rds")
selectGenes <- rownames(gtex_active[order(gtex_active$var,decreasing = T),][1:1000,])


gtex_anno <- read.table("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",sep="\t",head=TRUE,quote = "",fill=NA)
gtex_anno$SAMPID <- gsub("-",".",gtex_anno$SAMPID)
gtex_anno <- subset(gtex_anno,gtex_anno$SAMPID %in% colnames(gtex_active))
gtex_anno_tissue <- cbind(gtex_anno$SAMPID, as.character(gtex_anno$SMTS))
colnames(gtex_anno_tissue) <- c("SampleID","Tissue")
gtex_anno_tissue <- as.data.frame(gtex_anno_tissue)
gtex_anno_tissue$Tissue <- as.factor(gtex_anno_tissue$Tissue)
ips_anno <- as.data.frame(colnames(gtex_ips_tpm)[17383:17838])
ips_anno$Tissue <- "ips"
colnames(ips_anno) <- c("SampleID","Tissue")
gtex_ips_anno <- rbind(gtex_anno_tissue,ips_anno)
gtex_ips_anno$Tissue <- as.factor(gtex_ips_anno$Tissue)



rownames(gtex_ips_anno) <- gtex_ips_anno[,1]
saveRDS(gtex_ips_anno,"gtex_ips_anno.rds")





pca_data <- t(gtex_active_clean[selectGenes,])
pca_data <- log2(pca_data+1)
pcaResults <- prcomp(pca_data)
saveRDS(pcaResults,"gtex_pca_1000.rds")

#pcaResults <- readRDS("gtex_pca_1000.rds")

pdf("gtex.tissue.pca.pdf",paper="a4r")
autoplot(pcaResults,data=gtex_anno_tissue,colour="Tissue",label.show.legend = FALSE)
dev.off()
}
#tsne
set.seed(999)
gtex_active <- readRDS("gtex_ips_active.rds")
gtex_anno_tissue <- readRDS("gtex_ips_anno.rds")
colnames(gtex_anno_tissue) <- c("ID","Tissue")

if(FALSE){
gtex_active$var <- apply(gtex_active,1,var)
gtex_active_clean<- gtex_active[,-ncol(gtex_active)]
selectGenes <- rownames(gtex_active[order(gtex_active$var,decreasing = T),][1:1000,])
tsne_data <- gtex_active_clean[selectGenes,]
tsne_data <- log2(tsne_data+1)saveRDS(tsne_data,"merged_tsne_data.rds")

#res <- M3C(tsne_data, des=gtex_anno_tissue, removeplots=TRUE, iters=25,fsize=8,maxK=20,analysistype="chi")
saveRDS(tsne_data,"merged_tsne_data.rds")
#saveRDS(res,"M3Cresults.rds")

#res <- M3C(tsne_data, des=gtex_anno_tissue,method=2,iters=25,maxK=20,analysistype="chi")
#saveRDS(res,"M3Cresults_fast.rds")
tsne <- Rtsne(t(tsne_data),dims=2,perplexity=30,verbose=TRUE,max_iter=500,num_threads=0)
saveRDS(tsne,"rtsne_result.rds")

embedding <- as.data.frame(tsne$Y)
embedding$Tissue <- as.factor(gtex_anno_tissue$Tissue)

#embedding$Tissue <- as.factor(ifelse(embedding$Tissue == "ips","ips","others"))

palette_64 <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
"#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
"#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
"#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
"#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
"#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
"#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
"#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
"#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
"#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
"#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
"#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

p <- ggplot(embedding, aes(x=V1, y=V2, color=Tissue)) +
     geom_point(size=1.25) +
     scale_colour_manual(values=palette_64)+
     guides(colour = guide_legend(override.aes = list(size=6))) +
     xlab("") + ylab("") +
     ggtitle("t-SNE Gtex (17382) and iPS(456)") +
     theme_light(base_size=20) +
     theme(strip.background = element_blank(),
           strip.text.x     = element_blank(),
           axis.text.x      = element_blank(),
           axis.text.y      = element_blank(),
           axis.ticks       = element_blank(),
           axis.line        = element_blank(),
           panel.border     = element_blank())

ggsave("rtsne_ipsgtex_color.pdf", p, width=24, height=18, units="in")
}

#colors <- rainbow(length(unique(gtex_anno_tissue$Tissue)))
#names(colors) <- unique(gtex_anno_tissue$Tissue)
#pdf("gtex_rtsne.pdf",paper="a4r")
#par(mgp=c(2.5,1,0))
#plot(tsne$Y, t='n', main="tSNE",xlab="tSNE dimension 1", ylab="tSNE dimension 2", "cex.main"=2,"cex.lab"=1.5)
#text(tsne$Y, labels=gtex_anno_tissue$Tissue,col=colors[gtex_anno_tissue$Tissue])
#dev.off()


#pdf("gtex.tissue.tsne.pdf",paper="a4r")
#tsne(tsne_data)
#dev.off()

## machine learning for establishing reference cohort data

#library("e1071")
#library("caret")
#library("e1071")
#library("glmnet")
#library("LiblineaR")
#library("HandTill2001")
#sampling and fold = y is label data 
makefolds <- function(y, cv.fold = 5){   
  n <- length(y)
  nlvl <- table(y)
  idx <- numeric(n)
  folds <- list()
  for (i in 1:length(nlvl)) {
    idx[which(y == levels(y)[i])] <- sample(rep(1:cv.fold,length = nlvl[i]))
  }
  for (i in 1:cv.fold){
    folds[[i]] <- list(train = which(idx!=i),
                       test =  which(idx==i)) 
  }  
  return(folds)
}

## Inner/Nested 5-fold CV loops ==> predict ==> score 1-5 (nested folds)  ==> train calibration model ==> PREDICT
makenestedfolds <- function(y, cv.fold = 5){
  nfolds <- list()
  folds <- makefolds(y,cv.fold)
  names(folds) <- paste0("outer",1:length(folds))
  for(k in 1:length(folds)){
    inner = makefolds(y[folds[[k]]$train],cv.fold)
    names(inner) <- paste0("inner",1:length(folds))
    for(i in 1:length(inner)){
      inner[[i]]$train <- folds[[k]]$train[inner[[i]]$train]
      inner[[i]]$test <- folds[[k]]$train[inner[[i]]$test]
    }
    nfolds[[k]] <- list(folds[k],inner) 
  }
  names(nfolds) <- paste0("outer",1:length(nfolds))
  return(nfolds)
}


nfolds <- makenestedfolds(as.factor(gtex_anno_tissue$Tissue))
saveRDS(nfolds,"./data/nfolds.rds")


#function `subfunc_load_betashdf5_subset_filter_match_save_betasKk()` -------------------------------------------------------------------------------------------------------

# 1. Load betas_v11.h5 (betashdf5) 
# 2. Subset K.k..train 
# 3. Unsupervised variance filtering (p = 10 000) 
# 4. Subset K.k.test 
# 5. Match CpG probes of filtered K.k.train to K.k.test 
# 6. Save also the full betas (2801 * 10000) but 10k CpG (are based on the respective K.k train set) 
# 7. Save betas.p.filtered.K.k into a separate folder

# Define utility / loader function ----------------------------------------------------------------------------------------------------------------


subfunc_load_betashdf5_subset_filter_match_save_betasKk <- function(K.start = 1, k.start = 0, n.cv.folds = 5, 
                                                                    nfolds.. = NULL,
                                                                    fpath.betasv11.hdf5 = NULL,
                                                                    p.var.filtering = 5000,
                                                                    out.path = "gtex_ips.varfilt.5k", out.fname = "gtex_ips.K.k"){
  
  # Check whether nfolds.. is provided
  nfolds.. <- nfolds
  # Check whether gtex_ips data existed 
  # rows: 9333 genes # cols: 17838 cell types
    message("Dimensions of loaded `gtex_ips_data` data file nrows: ", nrow(gtex_active), "; ncols: ", ncol(gtex_active))
 
  
  # Run CV scheme
  message("\nNested cross validation (CV) scheme starts ... ", Sys.time())
  for(K in K.start:n.cv.folds){  
    # Nested loop
    for(k in k.start:n.cv.folds){ 
      
      if(k > 0){ message("\n Subsetting & filtering inner/nested fold ", K,".", k,"  ... ",Sys.time())  
        fold <- nfolds..[[K]][[2]][[k]]  ### [[]][[2]][[]] means inner loop # Inner CV loops 1.1-1.5 (Fig. 1.)
      } else{                                                                          
        message("\n \nSubsetting & filtering outer fold ", K,".0  ... ",Sys.time()) 
        fold <- nfolds..[[K]][[1]][[1]]   ### [[]][[1]][[]] means outer loop # Outer CV loops 1.0-5.0 (Fig. 1.)
      }
      
 


     # Subset K.k$train
      message(" Step 2. Subsetting cases/columns: " , K, ".", k, " training set @ ", Sys.time())
      gtex.K.k.train <- gtex_active[ , fold$train] # rows genes # columns are cells types
      message(" Step 3. Unsupervised variance filtering of p = ", p.var.filtering,  
              " CpG probes on " , K, ".", k, " training set @ ", Sys.time(),
              "\n  It can take up to 1-2mins to finish.")
      # sd is calculated over all cols (i.e. cell_typpe) for each row (i.e. CpG probe) 
      gtex.p.filtered.K.k.train <- gtex.K.k.train[order(apply(gtex.K.k.train, 1, sd), decreasing = T)[1:5000], ] 
      message("  Dimension of `gtex.p.filtered.K.k.train` nrows: ", 
              nrow(gtex.p.filtered.K.k.train), " ncols: ", ncol(gtex.p.filtered.K.k.train),
              "\n  Variance filtering finished @ ", Sys.time()) # Duration @ single core ca. 1.25-1.5mins
      message(" \n  Check whether there is NA in train set : ", sum(is.na(gtex.p.filtered.K.k.train) == T))
      
      # Transposed afterwards!
      gtex.p.filtered.K.k.train <- t(gtex.p.filtered.K.k.train)
      # gtex.p.filtered.K.k.train # matrix 
      # fold$train (ca. 1700-2204) rows cases/patients (sentrixIDs) 
      # rows: patients; cols 10k most variable CpGs
      message("  Transposing `gtex.p.filtered.K.k.train` finished @ ", Sys.time())

      # Garbage collector (note: gc is not absolutely necessary)
      message("  Clean up memory (garbage collector) @ ", Sys.time())
      gc()
      
      # Subset genes of the corresping test set 
      # Select only 5000  CpG (`p.var.filtering`) genes (i.e. rows of betasv11.h5) that are filtered based on
      # the training (sub)fold (i.e. columns of betas.p.varfilt.train)
      message(" Step 4. Subsetting `betas_v11.h5` cases/columns: " , K, ".", k, " test/calibration set @ ", Sys.time())
      gtex.K.k.test <- gtex_active[ , fold$test]
      
      message(" Step 5. Matching variance filtered p = ", p.var.filtering,  
              " CpG probes corresponding to the " , K, ".", k, " training set @ ", Sys.time(),
              "\n  It can take up to 1-2mins to finish.")
      gtex.p.filtered.K.k.test <- gtex.K.k.test[match(colnames(gtex.p.filtered.K.k.train), rownames(gtex.K.k.test)), ]
      # Transpose $test
      # Note: wrappend in t() => rows: patients ; cols: CpG probenames
      gtex.p.filtered.K.k.test <- t(gtex.p.filtered.K.k.test)
      message("  Transposing `betas.p.filtered.K.k.test` finished @ ", Sys.time())
      message("  Dimension of `betas.p.filtered.K.k.test`  nrows: ", 
              nrow(gtex.p.filtered.K.k.test), " ncols: ", ncol(gtex.p.filtered.K.k.test),
              "\n CpG matching finished @ ", Sys.time())
      
      # Save also betas.K.k (2801 * 10k CpG selected on the training set)
      message(" Step 6. Matching variance filtered p = ", p.var.filtering, 
              " CpG probes corresponding to the " , K, ".", k, " training set @ ", Sys.time(),
              "\n  On the full `betas_v11.h5` data. It can take up to 1-2mins to finish.")
      gtex.K.k <- gtex_active[match(colnames(gtex.p.filtered.K.k.train), rownames(gtex_active)), ] 
      # rows = gens # columns = cell types # no subsetting
      message("  Transposing `betas.K.k` finished @ ", Sys.time())
      gtex.K.k <- t(gtex.K.k) # rows patients # cols CpGs
      

#########here
      # Security check
      message("\nAre column names (CpG probes) of $train and $test and full betas identical? ", 
              identical(colnames(gtex.p.filtered.K.k.train), colnames(gtex.p.filtered.K.k.test)))
      message("Are column names (CpG probes) of $train and full `betas.K.k`` identical? ", 
              identical(colnames(gtex.p.filtered.K.k.train), colnames(gtex.K.k)))
      message("Are column names (CpG probes) of $test and full `betas.K.k`` identical? ", 
              identical(colnames(gtex.p.filtered.K.k.train), colnames(gtex.K.k)))
      
      # Create output directory  
      folder.path <- file.path(getwd(), "data", out.path)
      dir.create(folder.path, showWarnings = F, recursive = T)
      #RData.path <- file.path(folder.path, paste(out.fname, K, k, "RData", sep = "."))
      
      # Save unsupervised variance filtered $train and $test sets
      save(gtex.p.filtered.K.k.train, 
           gtex.p.filtered.K.k.test,
           gtex.K.k,
           fold,
           file = file.path(folder.path, paste(out.fname, K, k, "RData", sep = "."))
      )  
      
    }
  }
}
subfunc_load_betashdf5_subset_filter_match_save_betasKk()


