#--------------------------------------------------------------------------------
# ml4calibrated450k - Support Vector Machines (SVM) -  Utility/subfunctions  
# 
#                   - Linear Kernel SVM (e1071)
# 2021.05.12
#-------------------------------------------------------------------------------

# Check, install|load recquired packages -----------------------------------------------------------------------------------------------
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

# if (!requireNamespace("LiblineaR", quietly = TRUE)) {
#   install.packages("LiblineaR")
#   library(LiblineaR) } else {library(LiblineaR)}


# Utility functions ---------------------------------------------------------------------------------------------------------------------------

## Nested CV scheduler / stopper ---------------------------------------------------------------------------------------------------------------------------
subfunc_nestedcv_scheduler <- function(K, K.start, K.stop, k.start, k.stop, n.cv.folds, n.cv.inner.folds){
  
  # Correctly stop nested CV
  if(K > K.start && K < K.stop) {
    k.start <- 0
    n.cv.inner.folds <- n.cv.folds 
    
  } else {
    if(K == K.start) {
      if(K.start != K.stop) { # && k.start != 0){
        n.cv.inner.folds <- n.cv.folds # k inner goes to .5
      } else { # K == K.start == K.stop && k.start == 0  
        n.cv.inner.folds <- k.stop
      }
    } else { # K == K.stop
      if(k.start != 0) k.start <- 0
      n.cv.inner.folds <- k.stop
    } 
  }
  res <- list(k.start = k.start, 
              n.cv.inner.folds = n.cv.inner.folds
  )
  return(res)
}

## Cost tuner subfunction ---------------------------------------------------------------------------------------------------------------------------
subfunc_svm_e1071_linear_train_tuner_mc <- function(data.xTrain, target.yTrain, 
                                                    mod.type = "C-classification", 
                                                    kernel. = "linear", 
                                                    scale. = T,
                                                    C.base = 10, C.min = -3, C.max = 3, 
                                                    n.CV = 5, verbose = T, 
                                                    seed = 1234, 
                                                    parallel = T, 
                                                    mc.cores = 4L){ 
  
  # Cost C grid + give feedback and Sys.time                                     
  Cost.l <- as.list(C.base^(C.min:C.max))
  message("\nCost (C) = ", paste(simplify2array(Cost.l), sep = " ", collapse = " ; "), 
          " ; \nNr. of iterations: ", length(Cost.l),
          "\nStart at ", Sys.time())
  # Predefine empty list for results 
  cvfit.e1071.linear.C.tuner <- list() 
  
  # Parallel 
  if(parallel){
    cvfit.e1071.linear.C.tuner <- mclapply(seq_along(Cost.l), function(i){   # uses only lenght(Cost.l) workers
      set.seed(seed+1, kind ="L'Ecuyer-CMRG")
      svm(x = data.xTrain, y = target.yTrain, scale = scale., type = mod.type,
          kernel = kernel., 
          cost = Cost.l[[i]], 
          cross = n.CV, 
          probability = T, 
          fitted = T) # gamma: default 1/n.features; tolerance: default: 0.001) 
    }, mc.preschedule = T, mc.set.seed = T, mc.cores = mc.cores)
    print(Sys.time()) 
    return(cvfit.e1071.linear.C.tuner)
    
    # Sequential                                       
  } else {  
    cvfit.e1071.linear.C.tuner <- lapply(seq_along(Cost.l), function(i){
      cat("\nTuning e1071 with linear kernel function C = ", Cost.l[[i]], " @ ", Sys.time())   
      set.seed(seed+1, kind ="L'Ecuyer-CMRG")
      svm(x = data.xTrain, y = target.yTrain, scale = scale., type = mod.type, 
          kernel = kernel., cost = Cost.l[[i]], cross = n.CV, probability = T, fitted = T)
    })
    print(Sys.time())
    return(cvfit.e1071.linear.C.tuner)
  }
}    

## Cost selector subfunction ---------------------------------------------------------------------------------------------------------------------------
subfunc_svm_e1071_linear_C_selector <- function(results.cvfit.e1071.linear.C.tuner, 
                                                C.base = 10, C.min = -3, C.max = 3, 
                                                n.CV = 5, verbose = T){
  
  Costs.l <- as.list(C.base^(C.min:C.max))
  # Print simplified version of each crossvalidated fold accuracy for eventual manual selection
  res.cvfit.svm.accuracies.nCV <- sapply(seq_along(C.base^(C.min:C.max)), function(i){
    simplify2array(results.cvfit.e1071.linear.C.tuner[[i]]$accuracies)})
  colnames(res.cvfit.svm.accuracies.nCV) <- paste0("Cost_", Costs.l)
  rownames(res.cvfit.svm.accuracies.nCV) <- paste0("nCV", seq(1, n.CV, 1))
  # Print matrix of all CV accuracies
  if(verbose){
    message("\nMatrix of all CV accuracies:")
    print(res.cvfit.svm.accuracies.nCV)
  }
  
  # Average accuracy 
  res.cvfit.svm.accuracies.mean <- sapply(seq_along(C.base^(C.min:C.max)), function(i){
    simplify2array(results.cvfit.e1071.linear.C.tuner[[i]]$tot.accuracy)})
  names(res.cvfit.svm.accuracies.mean) <- paste0("Cost_", Costs.l)
  # Same as: res.cvfit.svm.accuracies.mean <- apply(res.cvfit.svm.accuracies.nCV, 2, mean)
  # Print list of average CV accuracies/ $tot.accuracy
  if(verbose){
    message("\nMean CV accuracies:")
    print(res.cvfit.svm.accuracies.mean)
  }
  
  # Selection
  # Chooses the smallest C with highest 5-fold cross-validated accuracy among possible choices 
  # => if C is large enough anyway (see Appendix-N.5.) doesnt make a difference 
  # => saves also computation time if C is smaller # => error-margin/Nr of supp.vecs.
  C.selected <- Costs.l[[which.max(res.cvfit.svm.accuracies.mean)]] 
  message("\nCost parameter with highest ", n.CV, "-fold CV accuracy : C = ", C.selected, " ; ", 
          "\n Note: If more than one maximal accuracy exists, C returns the smallest cost parameter with highest accuracy.", 
          "\n Once C is large than a certain value, the obtained models have similar performances", 
          " (for theoretical proof see Theorem 3 of Keerthi and Lin, 2003)")
  res <- list(C.selected = C.selected, 
              mtx.accuracies.nCV = res.cvfit.svm.accuracies.nCV, 
              mtx.accuracies.mean = res.cvfit.svm.accuracies.mean)
  return(res)
  
  # Literature:
  # Important # A practical guide to LIBLINEAR - Fan, Chang, Hsiesh, Wang and Lin 2008 
  # <https://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf>
  # Appendix-N.5. Parameter Selection:
  # 1. Solvers in LIBLINEAR are not very sensitive to C. Once C is large than a certain value, 
  # the obtained models have similar performances. 
  # Theoretical proof: Theorem 3 of Keerthi and Lin (2003)
}

## Re-fit training data subfunction ---------------------------------------------------------------------------------------------------------------------------

subfunc_svm_e1071_linear_modfit_train <- function(C.tuned, 
                                                  data.xTrain, target.yTrain,
                                                  results.cvfit.e1071.linear.C.tuner, 
                                                  C.selector.accuracy.mean, 
                                                  use.fitted = T){  #res.svm.C.tuner.l
  
  message("\n\nRe-fitting training data ... ", Sys.time())
  #Costs.l <- as.list(C.base^(C.min:C.max))
  i <- which.max(C.selector.accuracy.mean)
  
  # Ver.1 - Use predict function to refit training data
  # Note If the training set was scaled by svm (done by default), 
  # the new data is scaled accordingly using scale and center of the training data.
  modfit.train.svm.lin.pred <- predict(object = results.cvfit.e1071.linear.C.tuner[[i]], 
                                       newdata =  data.xTrain,
                                       decision.values = T,   
                                       # decision values of all binary classif. in multiclass setting are returned.
                                       probability = T)
  
  # Ver.2 - Use fitted() - see ??svm - Examples
  if(use.fitted){modfit.train.svm.lin.fitted <- fitted(results.cvfit.e1071.linear.C.tuner[[i]])}
  #message("\nBoth predict() & fitted() are ready @ ", Sys.time())
  message("\nPrediction is ready @ ", Sys.time())
  
  # Output
  res <- list(svm.e1071.model.object = results.cvfit.e1071.linear.C.tuner[[i]], 
              trainfit.svm.lin1 = modfit.train.svm.lin.pred, 
              trainfit.svm.lin2 = modfit.train.svm.lin.fitted) # => output file is rel. large / lot of large mtx.
  return(res)
}


#--------------------------------------------------------------------------------
# ml4calibrated450k - Support Vector Machines (SVM) -  Train function  
# 
#                   - Linear Kernel SVM (e1071)
#
# 2019-04-27 
#-------------------------------------------------------------------------------


### Training & tuning function - integrating the utility/subfuctions from `subfunctions_SVM_e1071.R`

train_SVM_e1071_LK <- function(y, betas.Train, 
                               seed, 
                               mc.cores, 
                               nfolds = 5, 
                               C.base = 10, C.min = -3, C.max = 3, 
                               scale.internally.by.e1071.svm = T,
                               mod.type = "C-classification"){  
  
  ## 1. Crossvalidate SVM/LiblineaR - Cost parameter for optimal
  set.seed(seed, kind = "L'Ecuyer-CMRG") 
  message("seed: ", seed)
  message("n: ", nrow(betas.Train))  # n_patients
  
  message("\nTuning SVM (e1071) linear kernel: hyperparameter C (cost) ... ", Sys.time())
  t1 <- system.time(
    cvfit.svm.e1071.linear.C.tuner <- subfunc_svm_e1071_linear_train_tuner_mc(data.xTrain = betas.Train, 
                                                                              target.yTrain = y,
                                                                              mod.type = mod.type, 
                                                                              kernel. = "linear", 
                                                                              scale. = scale.internally.by.e1071.svm,
                                                                              C.base = C.base, 
                                                                              C.min = C.min, C.max = C.max,
                                                                              n.CV = 5, 
                                                                              verbose = T, 
                                                                              seed = seed,
                                                                              parallel = T, 
                                                                              mc.cores = mc.cores)
  )
  
  # Extract optimal C or smallest C with highest accuracy 
  C.tuned.cv <-  subfunc_svm_e1071_linear_C_selector(results.cvfit.e1071.linear.C.tuner = cvfit.svm.e1071.linear.C.tuner, 
                                                     C.base = C.base, C.min = C.min, C.max = C.max, n.CV = nfolds, verbose = T)
  # C.tuned.cv = list of 3: $C.selected $mtx.accuracies.nCV $mtx.accuracies.mean
  # Provide message with value
  message(paste0("Optimal cost (C) parameter: ", C.tuned.cv$C.selected))
  
  # Refit models on s.xTrain L2R_LR (type0 - ca. 15 min/refit) and optionally only for classes Crammer & Singer (type4 - it takes just ca. +35s)
  message("\n(Re)Fitting optimal/tuned model on training data ... ", Sys.time())
  t2 <-  system.time(
    modfit.svm.linear.train <- subfunc_svm_e1071_linear_modfit_train(C.tuned = C.tuned.cv$C.selected, 
                                                                     data.xTrain = betas.Train, 
                                                                     target.yTrain = y, 
                                                                     results.cvfit.e1071.linear.C.tuner = cvfit.svm.e1071.linear.C.tuner, 
                                                                     C.selector.accuracy.mean = C.tuned.cv$mtx.accuracies.mean, 
                                                                     use.fitted = T)
    # uses predict supposed to scale data.xTrain / betas.Train automatically
    
  )
  # CAVE conames order is not the same as in levels(y) !!!
  pred.scores.trainfit.svm.lin1 <- attr(modfit.svm.linear.train$trainfit.svm.lin1, "probabilities") 
  
  # Results
  res <- list(modfit.svm.linear.train$svm.e1071.model.object,      # svm linear model object used in predict/fitted
              modfit.svm.linear.train$trainfit.svm.lin1,           # fitting train with predict.svm()/predict()
              pred.scores.trainfit.svm.lin1,                       # CAVE: colnames order is not the same as in levels(y) !!!
              modfit.svm.linear.train$trainfit.svm.lin2,           # fitting train with fitted() - as in the examples of svm{e1071} help
              cvfit.svm.e1071.linear.C.tuner,                      # list of 7: model objects of svm based on Cost tuner grid 10^(-3:3)
              C.tuned.cv,                                          # list of 3: $C.selected $mtx.accuracies.nCV $mtx.accuracies.mean 
              t1, t2) 
  return(res)
}

#--------------------------------------------------------------------------------
# ml4calibrated450k - Support Vector Machines (SVM) -  Run the nested CV scheme  
# 
#                   - Linear Kernel SVM (e1071)
#
# Matt Maros
# maros@uni-heidelberg.de
#
# 2019-04-27 
#-------------------------------------------------------------------------------

### Fit linear kernel SVM (SVM-LK) of the e1071 package in the integrated nested CV scheme 

run_nestedcv_SVM_e1071 <- function(y.. = NULL, 
                                   betas.. = NULL, 
                                   path.betas.var.filtered = "./data/gtex_ips.varfilt.5k/",
                                   fname.betas.p.varfilt = "gtex_ips.K.k",
                                   subset.CpGs.1k = F, 
                                   n.cv.folds = 5, 
                                   nfolds.. = NULL,   
                                   K.start = 1, k.start = 0,
                                   K.stop = NULL, k.stop = NULL, 
                                   n.CV. = 5,  # extra nested tuning in training loop
                                   C.base = 10, C.min = -3, C.max = 3, 
                                   n.mc.cores = 8L,   # standard 4 cores/8 threads machines
                                   seed = 1234, 
                                   out.path = "SVM-e1071-10k", 
                                   out.fname = "CVfold"){
  
  # Check:
  # Check whether y.. is provided
  
  # Check whether nfolds.. is provided
  gtex_anno_tissue <- readRDS("gtex_ips_anno.rds")
  colnames(gtex_anno_tissue) <- c("ID","Tissue")
  y.. <- as.factor(gtex_anno_tissue$Tissue)
  y <- y..

  # Check whether nfolds.. is provided
  nfolds.. <- readRDS("./data/nfolds.rds")
  nfolds <- nfolds.. 
  # <NOTE> if one of K|k.stop = NULL this gives an error. # Error in 1:n.cv.inner.folds : argument of length 0
  if(is.null(K.stop) && is.null(k.stop)) {
    message("\nK.stop & k.stop are at default NULL => the full `n.cv.folds` = ", n.cv.folds, " nested CV will be performed.")
    K.stop <- n.cv.outer.folds <- n.cv.folds
    k.stop <- n.cv.inner.folds <- n.cv.folds  
} else { # !is.null() # NOT NULL
    message("K.stop & k.stop are provided & nested CV is limited accordingly.")
    n.cv.outer.folds <-  K.stop
    n.cv.inner.folds <- k.stop  
  }
  
  # Start nested CV scheme:
  # Outer loop
  for(K in K.start:n.cv.outer.folds){
    
    # Schedule/Stop nested CV
    ncv.scheduler  <- subfunc_nestedcv_scheduler(K = K, 
                                                 K.start = K.start, K.stop = K.stop, 
                                                 k.start = k.start, k.stop = k.stop, 
                                                 n.cv.folds = n.cv.folds, 
                                                 n.cv.inner.folds = n.cv.inner.folds)
    k.start <- ncv.scheduler$k.start
    n.cv.inner.folds <- ncv.scheduler$n.cv.inner.folds
    
    # Inner/Nested loop
    for(k in k.start:n.cv.inner.folds){
      
      if(k > 0){ message("\n \nCalculating inner/nested fold ", K,".", k,"  ... ",Sys.time())  # Inner CV loops 1.1-1.5 (Fig. 1.)
        fold <- nfolds..[[K]][[2]][[k]]  ### [[]][[2]][[]] means inner loop
      } else{                                                                          
        message("\n \nCalculating outer fold ", K,".0  ... ",Sys.time()) # Outer CV loops 1.0-5.0 (Fig. 1.)
        fold <- nfolds..[[K]][[1]][[1]]   ### [[]][[1]][[]] means outer loop 
      }
      
      # 1. Load `betas.. 
      # Default is betas.. = NULL => Loads data from path
      if(is.null(betas..)) {
        # Load pre-filtered normalized but non-batchadjusted betas for fold K.k
        message("Setp 1. Loading pre-filtered normalized but non-batchadjusted betas for (sub)fold ", K, ".", k)
        # Safe loading into separate env
        env2load <- environment()
        # Define path (use defaults)
        path2load <- file.path(path.betas.var.filtered) # file.path("./data/betas.var.filtered/") # default 
        fname2load <- file.path(path2load, paste(fname.betas.p.varfilt, K, k, "RData", sep = "."))  
        # Load into env
        load(file = fname2load, envir = env2load)
        # Get 
        betas.train <- get(x = "gtex.p.filtered.K.k.train", envir = env2load)
        betas.test <- get(x = "gtex.p.filtered.K.k.test", envir = env2load)
        # Note that betas.train and betas.test columns/CpGs are ordered in deacreasing = T => simply subset [ , 1:1000] => 1k most variable
        if(subset.CpGs.1k) {
          betas.train <- betas.train[ , 1:1000] # matrix
          betas.test <- betas.test[ , 1:1000]   # matrix
        }
      } else { 
        # User provided `betas..` (e.g. `betas2801.1.0.RData`) including (`betas2801`) both $train & $test 
        # (only for a given fold => set K.start, k.start accordingly)
        message("\n<NOTE>: This is a legacy option. The `betas.. object should contain variance filtered cases of 2801 cases", 
                " according to the respective training set of (sub)fold ", K, ".", k, ".", 
                "\nThis option should be used only for a single fold corresponding to ", K.start, ".", k.start)
        betas.train <- betas..[fold$train, ] 
        betas.test <- betas..[fold$test, ]
        # If subset to 1k TRUE 
        if(subset.CpGs.1k) {
          betas.train <- betas.train[ , 1:1000]
          betas.test <- betas.test[ , 1:1000]
        }
      }
      
      message("\nStep 2. Pre-processing: scaling and centering training set ... ", Sys.time())
      s.betas.Train <- scale(betas.train, center = TRUE, scale = TRUE)
      
      # # Tune & train on training set # use internal scaling of e1071::svm()
      message("\nStart tuning on training set (scale & center happens internally within e1071::svm() ... ", Sys.time())
      message("\nScaling is performed internally within e1071::svm() before training,", 
              "and also using predict.svm() on the test/validation set according to the attributes of the training data")
      svm.linearcv <- train_SVM_e1071_LK(y = y..[fold$train], 
                                         betas.Train = betas.train, 
                                         seed = seed + 1, 
                                         nfolds = n.CV., 
                                         mc.cores = n.mc.cores,
                                         C.base = C.base, C.min = C.min, C.max = C.max)
      # Predict test set
      message("\nPredict SVM linear kernel {e1071} model with tuned cost (C) parameter ... ", 
              "\n Note predict.svm(): \'If the training set was scaled by svm (done by default),", 
              " the new data is scaled accordingly using scale and center of the training data.\' ...", Sys.time())
      scores.pred.svm.e1071.obj <- predict(object = svm.linearcv[[1]], 
                                           newdata = betas.test, 
                                           probability = T, 
                                           decisionValues = TRUE) 
      # probs.pred.SVM.e1071.obj => is a factor with attributes 
      # Get probabilities
      scores.pred.svm.e1071.mtx <- attr(scores.pred.svm.e1071.obj, "probabilities")
      # !!!CAVE: colnames() order might not be the same as in levels(y) originally!!!
      
      # Calculate Error                                                                                 
      # SVM Linear e1071: 
      err.svm.e1071.probs <- sum(colnames(scores.pred.svm.e1071.mtx)[apply(scores.pred.svm.e1071.mtx, 1, which.max)] != y..[fold$test]) / length(fold$test) 
      message("\nMisclassification error on test set estimated using [probabilities matrix] output: ",
              err.svm.e1071.probs, " ; ", Sys.time())
      
      # Control Steps
      message("\nControl step: whether rownames are identical (betas$K.k.fold$test):", 
              identical(rownames(scores.pred.svm.e1071.mtx), rownames(betas.test))) 
      message("Control step: whether colnames are identical (betas$K.k.fold$test):", 
              identical(colnames(scores.pred.svm.e1071.mtx), levels(y..))) 
      if(identical(colnames(scores.pred.svm.e1071.mtx), levels(y..[fold$test])) == FALSE){
        message("CAVE: Order of levels(y) and colnames(probs.pred.SVM.e1071.mtx)", 
                " => needs matching during performance evaluation!")
      }
      
      message("\nStep 5. Saving output objects & creating output folder (if necessary) @ ", Sys.time())
      # Create output directory  
      folder.path <- file.path(getwd(), out.path)
      dir.create(folder.path, recursive = T, showWarnings = F)
      #RData.path <- file.path(folder.path, paste(out.fname, K, k, "RData", sep = "."))
      
      # Save scores, SVM-LIBLINEAR-Modell, fold
      save(scores.pred.svm.e1071.mtx, 
           scores.pred.svm.e1071.obj, 
           svm.linearcv, 
           fold, 
           file = file.path(folder.path, paste(out.fname, K, k, "RData", sep = "."))
      )
    }
  }
  message("Full run finished @ ", Sys.time())
}



