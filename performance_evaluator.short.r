###     Performance evaluator   ###
###     2021.05.20              ###
###     YingChen                ###

# Includind all evuluated matrics
#       -Brier score(BS)
#       -Misclassification Error (ME)
#       -Multiclass AUC (AUC)
#       -Multiclass log loss (LL)

# Brier score (BS) -------------------------------------------------------------------------------------------------------------------------------------------------------------------

brier <- function(scores,y){
  ot <- matrix(0,nrow=nrow(scores),ncol=ncol(scores))
  # It can cause error: 
  # encountered errors in user code, all values of the jobs will be affectedError in matrix(0, nrow = nrow(scores), ncol = ncol(scores))
  arr.ind <- cbind(1:nrow(scores),match(y,colnames(scores)))
  ot[arr.ind] <- 1
  sum((scores - ot)^2)/nrow(scores)
}



#  Multiclass log loss (LL) ----------------------------------------------------------------------------------------------------------------------------------------------------------

mlogloss <- function(scores,y){
  N <- nrow(scores)
  y_true <- matrix(0,nrow=nrow(scores),ncol=ncol(scores))
  arr.ind <- cbind(1:nrow(scores),match(y,colnames(scores)))
  y_true[arr.ind] <- 1
  eps <- 1e-15 # we use Kaggle's definition of multiclass log loss with this constrain (eps) on extremly marginal scores (see reference below)  
  scores <- pmax(pmin(scores, 1 - eps), eps)
  (-1 / N) * sum(y_true * log(scores))
}

# Reference 
# <https://web.archive.org/web/20160316134526/https://www.kaggle.com/wiki/MultiClassLogLoss>



# Misclassification Error (ME) -------------------------------------------------------------------------------------------------------------------------------------------------------

# Subfunction for misclassification error (ME)
subfunc_misclassification_rate <- function(y.true.class, y.predicted){
  error_misclass <- sum(y.true.class != y.predicted)/length(y.true.class)
  return(error_misclass)
}



# Multiclass AUC after (Hand & Till 2001) --------------------------------------------------------------------------------------------------------------------------------------------- 

# Check & install package HandTill2001 if not in namespace & load
if (!requireNamespace("HandTill2001", quietly = TRUE)) {
  install.packages("HandTill2001")
  library(HandTill2001) } else {library(HandTill2001)}

# Subfunction for multiclass AUC & ROC by Hand & Till 2001

# Note: sum of row scores/probabilities must be scaled to 1 
subfunc_multiclass_AUC_HandTill2001 <- function(y.true.class, y.pred.matrix.rowsum.scaled1){
  auc_multiclass <- HandTill2001::auc(multcap(response = as.factor(y.true.class),
                                              predicted = y.pred.matrix.rowsum.scaled1))
  return(auc_multiclass)


}

 # 1. Checks
  # Check whether nfolds.. is available
     nfolds.. <- readRDS("./data/nfolds.rds")
     nfolds <- nfolds..


  # Check whether y true outcome is available
   
   gtex_anno_tissue <- readRDS("gtex_ips_anno.rds")
   gtex_anno_tissue <- gtex_anno_tissue[0:17382,]
  
   colnames(gtex_anno_tissue) <- c("ID","Tissue")
   gtex_anno_tissue$Tissue <- factor(gtex_anno_tissue$Tissue) 
   y.. <- as.factor(gtex_anno_tissue$Tissue)
   y <- y..
   anno.. <- gtex_anno_tissue
  



 if(length(anno..$ID) == length(y..)) {
    message("Checked: `anno..` and `y..` have the same corresponding dimensions. OK. \n")
  } else {
    stop("Error: anno.. and y.. have different dimensions.")}



  load.path.folder = "./SVM-LiblineaR"
  load.fname.stump = "CVfold"
  name.of.obj.to.load = NULL
  verbose=T
  scale.rowsum.to.1 = T
  misc.err = T
  multi.auc.HandTill2001 = T
 # brier = T
  mlogLoss = T

  probs.l <- list()
  for(i in 1:length(nfolds..)){
    # Create loading environment for safe loading the RData files of outerfolds
    env2load.outer <- environment()
    # Load
    if(verbose) message("\nLoad the test set from the ", i, ".0 outer fold @ ", Sys.time())
    path2load <- file.path(load.path.folder)
    fname2load <- file.path(path2load, paste(load.fname.stump, i, 0, "RData", sep = "."))
    load(fname2load, envir = env2load.outer)
    #probs.loaded <- get(as.character(name.of.obj.to.load), envir = env2load.outer)
    probs.loaded <- scores.pred.svm.liblinear.mod.type$probabilities
    # edit the row names of the sample
    test.sampleID <- fold$test

    rownames(probs.loaded) <- test.sampleID
    probs.loaded <- probs.loaded[ , levels(y..)]
    probs.l[[i]] <-  probs.loaded
  }
  if(verbose) message("\nAre predicted class labels (colnames) identically ordered as levels of the true outcome (y..): ", identical(colnames(probs), levels(y..)))
 if(verbose) message("\nCombine loaded outer fold test sets  1.0 ... ", length(nfolds..), ".0")
  probs <- do.call(rbind, probs.l)
  if(verbose) message("Dimensions of the combined matrix: ", nrow(probs), " x ", ncol(probs))
  probs <- probs[order(as.numeric(row.names(probs))),]

  # Rowsum check
  rowsum.p <- apply(probs, 1, sum)
  if(verbose) message("\nRow sum of first 10 cases :", paste(as.matrix(rowsum.p[1:10]), collapse = " ; "))
  # Show scores / probs ranges:
  p.range <- range(probs)
  if(verbose) message("The range of raw and/or calibrated object", " '", as.character(name.of.obj.to.load), "' ",
                      "is : ", "\n", p.range[[1]], " - ", p.range[[2]])

 # Predicted class column with highest prob for each row (sample))
  y.p.rowsum1 <- colnames(probs.rowsum1)[apply(probs.rowsum1,1, which.max)]
  y.l <- list(y.p.rowsum1)
  # Misclassification Error
  if(misc.err){
    err.misc.l <- lapply(y.l, subfunc_misclassification_rate, y.true.class = y..)
    err.misc <- unlist(err.misc.l)
    message("\nMisclassification Error: ", err.misc)
  }
  # AUC HandTIll2001
  if(multi.auc.HandTill2001){
    results.sc.p.rowsum1.l <- list(probs.rowsum1)
    message("Calculating multiclass AUC (Hand&Till 2001) ... ", Sys.time())
    auc.HT2001.l <- lapply(results.sc.p.rowsum1.l, subfunc_multiclass_AUC_HandTill2001, y.true.class = y..)
    auc.HT2001 <- unlist(auc.HT2001.l)
    message("Multiclass AUC (Hand&Till 2001): ", auc.HT2001)
  } else{
    message("Multiclass AUC (Hand&Till 2001) is set to FALSE =>",
            "Calculating RAW probs (row sum != 1). auc.HT2001 is set to NA.")
    auc.HT2001 <- NA
  }
  brierp.rowsum1 <- brier(scores = probs.rowsum1, y = y..)
  message("Brier score (BS): ", brierp.rowsum1)
  
 if(mlogLoss){
    loglp.rowsum1 <- mlogloss(scores = probs.rowsum1, y = y..)
    message("Multiclass log loss (LL): ", loglp.rowsum1)
  }

 # Results
  res <- list(misc.error = err.misc,
              auc.HandTill = auc.HT2001,
              brier = brierp.rowsum1,
              mlogloss = loglp.rowsum1)






