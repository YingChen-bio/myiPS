#!/usr/local/package/r/3.6.0/bin/Rscript --vanilla

#function to make fold for CV
###	ying.chen	###
###	2021.05.27	###


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

## Inner/Nested 5-fold CV loops ==> predict ==> score 1-5 (nested folds)  ==><train calibration model>notfinishd  ==> PREDICT
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


## Command to excute functions

gtex_anno_tissue <- readRDS("gtex_anno_tissue.rds")

#only gtex
gtex_active <- readRDS("gtex_active.rds")
nfolds <- makenestedfolds(as.factor(gtex_anno_tissue$Tissue))
saveRDS(nfolds,"./data/nfolds.rds")

subfunc_load_betashdf5_subset_filter_match_save_betasKk <- function(K.start = 1, k.start = 0, n.cv.folds = 5,
                                                                    nfolds.. = NULL,
                                                                    fpath.betasv11.hdf5 = NULL,
                                                                    p.var.filtering = 5000,
                                                                    out.path = "gtex.varfilt.5k", out.fname = "gtex.K.k"){

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
              " genes  on " , K, ".", k, " training set @ ", Sys.time(),
              "\n  It can take up to 1-2mins to finish.")
      # sd is calculated over all cols (i.e. tissue_type) for each row (i.e. ENSG genes) 
      gtex.p.filtered.K.k.train <- gtex.K.k.train[order(apply(gtex.K.k.train, 1, sd), decreasing = T)[1:5000], ]
      message("  Dimension of `gtex.p.filtered.K.k.train` nrows: ",
              nrow(gtex.p.filtered.K.k.train), " ncols: ", ncol(gtex.p.filtered.K.k.train),
              "\n  Variance filtering finished @ ", Sys.time()) # Duration @ single core ca. 1.25-1.5mins
      message(" \n  Check whether there is NA in train set : ", sum(is.na(gtex.p.filtered.K.k.train) == T))

      # Transposed afterwards!
      gtex.p.filtered.K.k.train <- t(gtex.p.filtered.K.k.train)
      # gtex.p.filtered.K.k.train # matrix 
      # fold$train (ca. 1700-2204) rows cases/patients (sentrixIDs) 
      # rows: gtex_tissue_samples; cols 10k most variable genes
      message("  Transposing `gtex.p.filtered.K.k.train` finished @ ", Sys.time())

      # Garbage collector (note: gc is not absolutely necessary)
      message("  Clean up memory (garbage collector) @ ", Sys.time())
      gc()

      # Subset genes of the corresping test set 
      # Select only 5000  CpG (`p.var.filtering`) genes (i.e. rows of betasv11.h5) that are filtered based on
      # the training (sub)fold (i.e. columns of betas.p.varfilt.train)
      message(" Step 4. Subsetting `gtex_17382` cases/columns: " , K, ".", k, " test/calibration set @ ", Sys.time())
      gtex.K.k.test <- gtex_active[ , fold$test]

      message(" Step 5. Matching variance filtered p = ", p.var.filtering,
              " CpG probes corresponding to the " , K, ".", k, " training set @ ", Sys.time(),
              "\n  It can take up to 1-2mins to finish.")
      gtex.p.filtered.K.k.test <- gtex.K.k.test[match(colnames(gtex.p.filtered.K.k.train), rownames(gtex.K.k.test)), ]
      # Transpose $test
      # Note: wrappend in t() => rows: tissue_samples ; cols: genes
      gtex.p.filtered.K.k.test <- t(gtex.p.filtered.K.k.test)
      message("  Transposing `betas.p.filtered.K.k.test` finished @ ", Sys.time())
      message("  Dimension of `betas.p.filtered.K.k.test`  nrows: ",
              nrow(gtex.p.filtered.K.k.test), " ncols: ", ncol(gtex.p.filtered.K.k.test),
              "\n Gene  matching finished @ ", Sys.time())

      # Save also betas.K.k (2801 * 10k CpG selected on the training set)
      message(" Step 6. Matching variance filtered p = ", p.var.filtering,
              " CpG probes corresponding to the " , K, ".", k, " training set @ ", Sys.time(),
              "\n  On the full `gtex` data. It can take up to 1-2mins to finish.")
      gtex.K.k <- gtex_active[match(colnames(gtex.p.filtered.K.k.train), rownames(gtex_active)), ]
      # rows = gens # columns = tissue types # no subsetting
      message("  Transposing `betas.K.k` finished @ ", Sys.time())
      gtex.K.k <- t(gtex.K.k) # rows tissue_samples # cols genes


#########here
      # Security check
      message("\nAre column names (gene names) of $train and $test and full betas identical? ",
              identical(colnames(gtex.p.filtered.K.k.train), colnames(gtex.p.filtered.K.k.test)))
      message("Are column names (gene names) of $train and full `betas.K.k`` identical? ",
              identical(colnames(gtex.p.filtered.K.k.train), colnames(gtex.K.k)))
      message("Are column names (gene names) of $test and full `betas.K.k`` identical? ",
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

