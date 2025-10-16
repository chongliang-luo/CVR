
## run CVR  with simulated data with survival outcome
require(data.table)
require(survival)
require(survAUC)
require(glmnet)
require(PMA) # sparse CCA

Rcpp::sourceCpp("src/cvrsolver.cpp")
source("R/myfunctions.R")
dd0 = readRDS('data/CVR_Cox_simu_data.RDS')

fit0 = CVR(Y=dd0$y[,1], Xlist=list(dd0$X1, dd0$X2), event=dd0$y[,2], rankseq=2, neta =10, nlam = 50,  
           Lamseq = NULL, family = "cox", Wini = NULL, penalty = "GL1" )