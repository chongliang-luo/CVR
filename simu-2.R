
BiocManager::install("mixOmics")
help(package='mixOmics')

## DIABLO or sPLS-DA
# https://bioconductor.org/packages/release/bioc/vignettes/mixOmics/inst/doc/vignette.html#07
library(mixOmics)
# data(breast.TCGA)

# Extract training data and name each data frame
# Store as list
X <- list(X1 = X1, X2 = X2)

# Outcome: use event indicator as outcome as DIABLO can't handle survival outcome
# Y <- event
summary(event)

design <- matrix(0.1, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0
design 

# Number of components
diablo.tcga <- block.plsda(X, event, ncomp = 5, design = design)

set.seed(123) # For reproducibility, remove for your analyses
perf.diablo.tcga = perf(diablo.tcga, validation = 'Mfold', folds = 10, nrepeat = 10)

#perf.diablo.tcga$error.rate  # Lists the different types of error rates

# Plot of the error rates based on weighted vote
plot(perf.diablo.tcga)
ncomp <- perf.diablo.tcga$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]


# Number of variables to select 
# chunk takes about 2 min to run
set.seed(123) # for reproducibility
test.keepX <- list(X1 = c(5:19, seq(20, 50, 10)),
                   X2 = c(5:19, seq(20, 70, 10)) )

tune.diablo.tcga <- tune.block.splsda(X, event, ncomp = ncomp, 
                                      test.keepX = test.keepX, design = design,
                                      validation = 'Mfold', folds = 10, nrepeat = 1, 
                                      BPPARAM = BiocParallel::SnowParam(workers = 2),
                                      dist = "centroids.dist")
list.keepX <- tune.diablo.tcga$choice.keepX

# final model
diablo.tcga <- block.splsda(X, event, ncomp = ncomp, keepX = list.keepX, design = design)
diablo.tcga

# selected vars in each view
selectVar(diablo.tcga, block = 'X1', comp = 1)
selectVar(diablo.tcga, block = 'X2', comp = 1)

diablo.tcga$loadings$X1
# diablo.tcga$variates 
# diablo.tcga$names 
diablo.tcga$weights$comp1

# linear predictor
lptest = cbind(X1test %*% diablo.tcga$loadings$X1, X2test %*% diablo.tcga$loadings$X2) %*% diablo.tcga$weights$comp1
# predict survival outcome of test set
pred_diablo = UnoC(Surv(Y, event), Surv(Ytest, eventtest), lpnew=lptest)
