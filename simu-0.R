
## run CVR  with simulated data with survival outcome
require(data.table)
require(survival)
require(survAUC)
require(glmnet)
require(PMA) # sparse CCA

Rcpp::sourceCpp("src/cvrsolver.cpp")
source("R/myfunctions.R")
dd0 = readRDS('data/CVR_Cox_simu_data.RDS')

# fit CVR-Cox
fit0 = CVR(Y=dd0$y[,1], Xlist=list(dd0$X1, dd0$X2), event=dd0$y[,2], 
           rankseq=2, # use fixed rank=2
           neta =10,  # eta to balance correlation and prediction (CVR paper eq 2.5)
           nlam = 50, # lam to induce sparsity of X's coef
           Lamseq = NULL, family = "cox", Wini = NULL, penalty = "GL1" )

# predict test data, output Uno' concordance  
ff = fit0$cvout$cvr.fit
pred1 = PredCVR(dd0$ytest[,1], dd0$X1test, dd0$X2test, dd0$ytest[,2], 
                ff$W[[1]], ff$W[[2]], alpha = ff$alpha, beta =ff$beta, 
                Y=dd0$y[,1], X1=dd0$X1, X2=dd0$X2, event=dd0$y[,2], 
                family = "cox", refit = T, combine = T, type.measure = 'UnoC')

# CV_1(protein*W1, bb*W2)
# predict y using protein*W1 * b1 + bb*W2 * b2

