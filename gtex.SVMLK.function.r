# LR  with linear kernels (LK) using the `LiblineaR`package 

# 2021.05.10
# yingchen

if (!requireNamespace("doMC", quietly = TRUE)) {
  install.packages("doMC")
  library(doMC) } else {library(doMC)}

if (!requireNamespace("e1071", quietly = TRUE)) { 
  install.packages("e1071")
  library(e1071) } else {library(e1071)}

if (!requireNamespace("LiblineaR", quietly = TRUE)) {
  install.packages("LiblineaR")
  library(LiblineaR) } else {library(LiblineaR)}




##########################################
###  Nested CV scheduler / stopper     ###
##########################################

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




#####################################################################
###             SUBFUNCTION - LIBLINEAR - PREDICTOR               ###
### FOR VERSION 1 - DOUBLE FOREACH - Both MODEL TYPE & COST TUNER ###
#####################################################################
# 
# subfunc_SVM_LiblineaR_predictor <- function(selected.bestType.l, selected.bestCost.l, 
#                                             xTrain.scaled.cent, yTrain, xTest, yTest, scale.center.xTest.as.xTrain = T){
#   # Can internally scale & center data 
#   if(scale.center.xTest.as.xTrain) {s2 = scale(xTest, attr(xTrain.scaled.cent, "scaled:center"), attr(xTrain.scaled.cent, "scaled:scale"))} 
#   m <- as.list()
#   mclapply(seq_along(selected.bestType.l), function(i) {
#     bestType.i <- selected.bestType.l[[i]]
#     
#     lapply(seq_along(selected.bestCost.l), fuction(j) {
#      m = LiblineaR(data = xTrain.scaled.cent, 
#                    target = yTrain, 
#                    type = bestType.i, 
#                    cost = bestCost[[j]]) #, bias = 1, verbose = FALSE) # are defaults
#       # Set probability output to FALSE
#       pr = FALSE
#       if(bestType.i == 0 || bestType.i == 7) pr = TRUE # only type 0 & type 2 models can generate probabilities
#       p = predict(m, s2, proba = pr, decisionValues = TRUE)
#       res <- list(mod.fit = m, probs.pred = p)
#       return(res)
#     })
#   })
# }

#############################################################################
### SUBFUNCTION TO RUN LIBLINEAR SVM - COST TUNER with MCLAPPLY           ###
###   for a given model type (type 0)                                     ###
#############################################################################

subfunc_svm_liblinear_train_tuner_mc <- function(data.scaled.xTrain, 
                                                 target.yTrain, 
                                                 mod.type = 0,
                                                 C.base = 10, C.min = -3, C.max = 3,
                                                 bias = 1, #default liblinear setting
                                                 n.CV = 5, 
                                                 verbose = T, 
                                                 seed = 1234, 
                                                 parallel = T, 
                                                 mc.cores){ 
  
  # Cost C grid + give feedback and Sys.time                                     
  Cost.l <- as.list(C.base^(C.min:C.max))
  message("Cost (C) = ", paste(simplify2array(Cost.l), sep = " ", collapse = " ; "),
          " ; \nNr. of CV folds: ", n.CV,
          "\nStart at ", Sys.time())
  # Predefine empty list for results 
  cvfit.liblineaR.C.tuner <- list() 
  
  # Parallelized for Cost tuning
  if(parallel){
    cvfit.liblineaR.C.tuner <- mclapply(seq_along(Cost.l), function(i){
      #message("Tuning LIBLINEAR type 0 - L2-reg. LR with C = ", Cost.l[[i]], " @ ", Sys.time())   
      set.seed(seed+1, kind ="L'Ecuyer-CMRG")
      LiblineaR(data = data.scaled.xTrain, target = target.yTrain, 
                cost = Cost.l[[i]], 
                type = mod.type, 
                bias = bias, cross = n.CV, verbose = verbose)
    }, mc.preschedule = T, mc.set.seed = T, mc.cores = mc.cores)            
    print(Sys.time())                                                       
    return(cvfit.liblineaR.C.tuner)                                         
  } else { # Sequential 
    cvfit.liblineaR.C.tuner <- lapply(seq_along(Cost.l), function(i){
      message("Tuning LIBLINEAR type 0 - L2-reg. LR with C = ", Cost.l[[i]], " @ ", Sys.time())   
      set.seed(seed+1, kind ="L'Ecuyer-CMRG")
      LiblineaR(data = data.scaled.xTrain, target = target.yTrain, type = mod.type, 
                cost = Cost.l[[i]], bias = 1, cross = n.CV, verbose = verbose) 
      # cross = LiblineaR help: if an integer value k>0 is specified, a k-fold cross validation on data is performed 
      #         to assess the quality of the model via a measure of the accuracy. 
      #         Note that this metric might not be appropriate if classes are largely unbalanced. Default is 0.
    })
    print(Sys.time())
    return(cvfit.liblineaR.C.tuner)
  }
}    

#############################################################
### SUBFUNCTION - SVM-LIBLINEAR - Cost parameter SELECTOR ###
#############################################################

# Literature: A practical guide to LIBLINEAR - Fan, Chang, Hsiesh, Wang and Lin 2008 
# <https://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf>
# Appendix - N.5. Parameter Selection:
# 1. Solvers in LIBLINEAR are not very sensitive to C. Once C is large than a certain value, the obtained models have similar performances. 
# Theoretical proof: Theorem 3 of Keerthi and Lin (2003)
  

subfunc_svm_liblinear_C_selector <- function(results.cvfit.liblineaR.C.tuner, 
                                             C.base = 10, C.min = -3, C.max = 3, 
                                             n.CV = 5){
    
  Costs.l <- as.list(C.base^(C.min:C.max))
  # Print simplified version for eventual manual selection
  names(results.cvfit.liblineaR.C.tuner) <- paste0("Cost_", Costs.l)
  print(simplify2array(results.cvfit.liblineaR.C.tuner))
  # Chooses the smallest C with highest 5-fold cross-validated accuracy among possible choices
  C.selected <- Costs.l[[which.max(results.cvfit.liblineaR.C.tuner)]] 
  message("\nCost parameter with highest ", n.CV, "-fold CV accuracy : C = ", C.selected, " ; ", 
          "\n If more than one maximal accuracy exists, C returns the smallest cost parameter with highest accuracy")
  return(C.selected)
}

####################################################################
### RE-FIT LINEAR SVM/LR MODEL with TUNED COST on s.xTRAIN       ###
###    THUS get model object to predict on test set              ###
####################################################################

subfunc_svm_liblinear_refit_train <- function(C.tuned, 
                                              data.scaled.xTrain, 
                                              target.yTrain, 
                                              mod.type = 0, 
                                              type4.CramSing = T, 
                                              verbose = T){
  if(type4.CramSing){ 
    # Crammer & Singer 
    # ca. 30 sec # tolerance f type is 1, 3, 4, 7, 12 or 13 # epsilon=0.1
    m.fit.ty4 <- LiblineaR(data = data.scaled.xTrain, target = target.yTrain, type = 4, cost = C.tuned, bias = 1, verbose = verbose)
    
    # L2-regularized Logsitic Regression - L2R_LR by default otherwise mod.type according to LiblineaR
    # 11mins # tolerance if type is 0, 2, 5 or 6 # epsilon=0.01
    m.fit <- LiblineaR(data = data.scaled.xTrain, target = target.yTrain, type = mod.type, cost = C.tuned, bias = 1, verbose = verbose) 
    res <- list(m.fit.mod.type = m.fit, m.fit.ty4.cs = m.fit.ty4)
  } else {
    m.fit <- LiblineaR(data = data.scaled.xTrain, target = target.yTrain, type = mod.type, cost = C.tuned, bias = 1, verbose = verbose)
    res <- list(m.fit.mod.type = m.fit, m.fit.ty4.cs = NULL)
  }
  return(res)
}




### Training & tuning function - integrating the utility/subfuctions from `subfunctions_SVM_LiblineaR.R`

train_SVM_LiblineaR <- function(y, 
                                s.betas.Train, 
                                seed, 
                                n.CV = 5, 
                                multicore = T, 
                                mc.cores, 
                                C.base = 10, C.min = -3, C.max = 3,
                                mod.type = 0, # defaults to type0 L2R_LR
                                type4.CramSing = T, # Crammer & Singer SVC is also fitted 
                                verbose = T){  
  
  ## 1. Crossvalidate SVM/LiblineaR - Cost parameter for optimal - L2R_LR (and Crammer & Singer models) 
  set.seed(seed, kind = "L'Ecuyer-CMRG") 
  message("seed: ", seed)
  message("n: ", nrow(s.betas.Train))
  if(multicore == T && mc.cores > 1) message("Parrallel computing with multicore: ", mc.cores)
  else message("Parralel computing (multicore = F) and/or mc.cores = 1 => Please revise (e.g. define backend)")
  
  model.names.4multiclass.LiblineaR <-c("0 – L2-regularized logistic regression (primal)", 
                                        "1 – L2-regularized L2-loss support vector classification (dual)", 
                                        "2 – L2-regularized L2-loss support vector classification (primal)", 
                                        "3 – L2-regularized L1-loss support vector classification (dual)", 
                                        "4 – support vector classification by Crammer and Singer", 
                                        "5 – L1-regularized L2-loss support vector classification", 
                                        "6 – L1-regularized logistic regression", 
                                        "7 – L2-regularized logistic regression (dual)")
  
  message("\nTuning linear kernel SVM hyperparameter C (cost) for `", 
          model.names.4multiclass.LiblineaR[[mod.type+1]], "`  ... ", Sys.time())
  if(type4.CramSing) message("\nAlso tuning for `", 
                             model.names.4multiclass.LiblineaR[[5]], "`  ... ", Sys.time())
  t1 <- system.time(
    cvfit.svm.liblinear.tuning <- subfunc_svm_liblinear_train_tuner_mc(data.scaled.xTrain = s.betas.Train, 
                                                                       target.yTrain = y, 
                                                                       mod.type = mod.type,  
                                                                       C.base = C.base, C.min = C.min, C.max = C.max, 
                                                                       bias = 1,
                                                                       n.CV = n.CV, 
                                                                       verbose = verbose,
                                                                       seed = seed,
                                                                       parallel = multicore, 
                                                                       mc.cores = mc.cores)
  )
  
  # Extract optimal C (i.e. smallest C with highest accuracy)
  C.tuned.cv <-  subfunc_svm_liblinear_C_selector(results.cvfit.liblineaR.C.tuner = cvfit.svm.liblinear.tuning, 
                                                  C.base = C.base, C.min = C.min, C.max = C.max, 
                                                  n.CV = n.CV)
  
  # Give message with values
  message(paste0("\nOptimal cost (C) parameter: ", C.tuned.cv))
  
  # Refit models on s.xTrain: 
  # - L2R_LR (type0 - ca. 15 min/refit) - gives prob output
  # - CS Crammer & Singer only for class output (ca. +35s)
  message("\nFitting optimal/tuned model on training data ... ", Sys.time())
  t2 <-  system.time(
    modfit.liblinear.train <- subfunc_svm_liblinear_refit_train(C.tuned = C.tuned.cv, 
                                                                data.scaled.xTrain = s.betas.Train, 
                                                                target.yTrain = y,
                                                                mod.type = mod.type, 
                                                                type4.CramSing = type4.CramSing)  
    # default is T for Cram&Sing # our prototyping showed that CS SVC is for all C value better/more accurate than type0
  )                                                                         
  
  # Results
  res <- list(modfit.liblinear.train$m.fit.mod.type, 
              modfit.liblinear.train$m.fit.ty4.cs, # if CS = F => NULL object
              cvfit.svm.liblinear.tuning, 
              C.tuned.cv, 
              t1, t2) 
  return(res)
}



### Fit linear kernel SVM (SVM-LK) of the `LiblineaR`` package in the integrated nested CV scheme 

## 3. Fit tuned `SVM LiblineaR` models in the integrated nested CV scheme

# 1. type 0: L2-regularized logistic regression (L2LR) and 
# 2. type 4: Crammer & Singer model 

### `run_nestedcv_SVM_LiblineaR` - hyperparameter (C, cost) tuning using 5-fold extra nested CV within the training loop 

run_nestedcv_SVM_LiblineaR <- function(y.. = NULL, 
                                       betas.. = NULL, 
                                       path.betas.var.filtered = "./data/gtex.varfilt.5k/",
                                       
                                       fname.betas.p.varfilt = "gtex.K.k",
                                       subset.CpGs.1k = F, 
                                       n.cv.folds = 5, 
                                       nfolds.. = NULL,   
                                       K.start = 1, k.start = 0,
                                       K.stop = NULL, k.stop = NULL, 
                                       n.CV. = 5,  # extra nested tuning in training loop
                                       C.base = 10, C.min = -3, C.max = 3, 
                                       mod.type = 0, # L2-LR with probability output
                                       type4.CramSing = T, # only class estimates no probability estimates
                                       verbose = T,
                                       parallel = T,
                                       n.mc.cores = 8L,   # standard 4 cores/8 threads machines
                                       seed = 1234, 
                                       out.path = "SVM-LiblineaR", 
                                       out.fname = "CVfold"){
  
  # Check:
  # Check whether y.. is provide y.. should be the label of the tissue
  gtex_anno_tissue <- readRDS("gtex_ips_anno.rds")
  #only gtex
  gtex_anno_tissue <- gtex_anno_tissue[0:17382,]
  colnames(gtex_anno_tissue) <- c("ID","Tissue")
  y.. <- as.factor(gtex_anno_tissue$Tissue)
  y <- y..
  
  # Check whether nfolds.. is provided
  nfolds.. <- readRDS("./data/nfolds.rds")
  nfolds <- nfolds..

  # Feedback messages
  # Check if K.stop & k.stop was provided
  if(is.null(K.stop) && is.null(k.stop)) {
    message("\nK.stop & k.stop are at default NULL => the full `n.cv.folds` = ", n.cv.folds, 
            " nested CV will be performed.")
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
      
      if(k > 0){ message("\n \nCalculating inner/nested fold ", K,".", k,"  ... ",Sys.time())  
        fold <- nfolds..[[K]][[2]][[k]]  ### Inner CV loops 1.1-1.5 (Fig. 1.)
      } else{                                                                          
        message("\n \nCalculating outer fold ", K,".0  ... ",Sys.time()) 
        fold <- nfolds..[[K]][[1]][[1]]   ### Outer CV loops 1.0-5.0 (Fig. 1.)
      }
      
      # 1. Load `betas.. 
      # Default is betas.. = NULL => Loads data from path
      if(is.null(betas..)) {
        # Load pre-filtered normalized but non-batchadjusted betas for fold K.k
        message("Step 1. Loading pre-filtered normalized but non-batchadjusted betas for (sub)fold ", K, ".", k)
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
        # Note that gtex.train and gtex.test columns/genes are ordered in deacreasing = T 
        # "Fast track" => simply subset [ , 1:1000] => 1k most variable 
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
      
      message("\nStep 3. Start tuning on training set ... ", Sys.time())
      # Run train- on scaled gtex TRAIN sets both Outer and Inner folds => "Train the model on train set (either nested or outer)"
      liblinearcv <- train_SVM_LiblineaR(y = y..[fold$train], 
                                         s.betas.Train = s.betas.Train, 
                                         n.CV = n.CV., 
                                         seed = seed+1, 
                                         multicore = parallel,  # if parallel=F (sequential) mc.cores is ignored
                                         mc.cores = n.mc.cores, 
                                         C.base = C.base, C.min = C.min, C.max = C.max, 
                                         mod.type = mod.type, 
                                         type4.CramSing = type4.CramSing, verbose = verbose)
      
      # Scale test set according to ranges of training set above
      message("\nStep 4. Scaling and centering test set using attributes of the trainin set ... ", Sys.time())
      s.betas.Test <- scale(betas.test, attr(s.betas.Train, "scaled:center"), attr(s.betas.Train, "scaled:scale"))  
      
      #message("\nCreating output folder (if necessary) @ ", Sys.time())
      # Create output directory  
      folder.path <- file.path(getwd(), out.path)
      dir.create(folder.path, recursive = T, showWarnings = F)
      
      # Use fitted model to predict the corresponding test set 
      message("\nStep 5. Predict on scaled `test set` using tuned model object ... ", Sys.time())
      # PROBABILITY OUTPUTS
      # only possible if the model was fitted with type=0, type=6 or type=7, i.e. a Logistic Regression. 
      # Default is FALSE
      if(mod.type == 0 || mod.type == 6 || mod.type == 7 && type4.CramSing == T){
        message("Probability output is only possible for model type = ", paste(c(0, 6, 7), sep = " ; "))
        scores.pred.svm.liblinear.mod.type <- predict(liblinearcv[[1]], s.betas.Test, proba = T, decisionValues=TRUE) 
        # [[1]] => modfit.liblinear.strain$m.fit.mod.type
        
        # Type 4 - CS model object 
        class.pred.svm.liblinear.ty4.CraSi <- predict(liblinearcv[[2]], s.betas.Test, proba = F, decisionValues=TRUE) 
        # [[2]] => modfit.liblinear.strain$m.fit.ty4.cs # if CS=T then exists otherwise NULL
        
        # Calculate Error                                                                                 
        # Specified model type: > default L2R-LogReg (type0)
        err.svm.mod.type <- sum(colnames(scores.pred.svm.liblinear.mod.type$probabilities)[apply(scores.pred.svm.liblinear.mod.type$probabilities,
                                                                                                 1, which.max)] != y..[fold$test])/length(fold$test) 
        message("\nMisclassification error of type (", mod.type, ") model on fold ", K, ".", k, " : ",
                err.svm.mod.type, " ; @ ", Sys.time())
        
        # Save scores, SVM-LIBLINEAR-Modell, fold
        message("\nStep 6. Saving output objects & creating output folder (if necessary) @ ", Sys.time())
        save(scores.pred.svm.liblinear.mod.type, 
             class.pred.svm.liblinear.ty4.CraSi, 
             liblinearcv, 
             fold, 
             file = file.path(folder.path, paste(out.fname, K, k, "RData", sep = ".")))
        
      } else {
        # Model types with only CLASS OUTPUTS                                                                   
        class.pred.svm.liblinear.ty4.CraSi <- predict(liblinearcv[[1]], s.betas.Test, proba = F, decisionValues=TRUE)
        # Calculate Error                                                                                 
        # Type 4 - Crammer & Singer
        err.svm.ty4 <- sum(class.pred.svm.liblinear.ty4.CraSi$predictions != y..[fold$test])/length(fold$test)
        message("\nMisclassification error fold ", K, ".", k, " : ", err.svm.ty4, " ; @ ", Sys.time())
        
        # Save scores, SVM-LIBLINEAR-Modell, fold
        message("\nStep 6. Saving output objects & creating output folder (if necessary) @ ", Sys.time())
        save(class.pred.svm.liblinear.ty4.CraSi, 
             liblinearcv, 
             fold, 
             file = file.path(folder.path, paste(out.fname, K, k, "RData", sep = ".")))
        
      }
    }
  }
  message("Full run finished ...", Sys.time())
}





























