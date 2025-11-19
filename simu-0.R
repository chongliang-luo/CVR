
## run CVR  with simulated data with survival outcome
require(data.table)
require(survival)
require(survAUC)
require(glmnet)
require(PMA) # sparse CCA
setwd('CVR/')


Rcpp::sourceCpp("src/cvrsolver.cpp")
source("R/myfunctions.R") 


########## X1 X2 share common canonical variates which also predict Y
set.seed(42)   
n= 400
mydata <- SimulateCVR(family = "cox", n =n, rank = 4, p1 = 200, p2 = 300,  
                      pnz = 30, beta = c(2, 1, 0, 0))
X1 <- mydata$X1 
X2 <- mydata$X2
Xlist <- list(X1 = X1, X2 = X2); 
Y <- mydata$y[,1]
event = mydata$y[,2]
Ytest=mydata$ytest[,1]
X1test=mydata$X1test 
X2test=mydata$X2test 
eventtest=mydata$ytest[,2]
 
## CVR fit, fix rank = 4, tune eta and lambda 
out_cvr <- CVR(Y, Xlist, event=event, rankseq= 4, neta= 5, nlam = 25, family = "cox", nfold = 5)
# out_cvr$solution$W[[1]]
ff = out_cvr$cvout$cvr.fit # refit
round(cor(X1 %*% ff$W[[1]], X2 %*% ff$W[[2]]), 2)

# predict test data, output Uno's concordance 
pred1 = PredCVR(mydata$ytest[,1], mydata$X1test, mydata$X2test, mydata$ytest[,2], 
                ff$W[[1]], ff$W[[2]], alpha = ff$alpha, beta =ff$beta, 
                Y=mydata$y[,1], X1=mydata$X1, X2=mydata$X2, event=mydata$y[,2], 
                family = "cox", refit = T, combine = T, type.measure = 'UnoC')
pred1 # 0.846

# individual prediction:
# linear predictor (rank): ff$alpha + (X1test %*% ff$W[[1]], X2test %*% ff$W[[2]]) * ff$beta 
# also need to estimate baseline haz
 
########## X1 X2 have common+individual canonical variates, both predict Y 
set.seed(42)   
n = 400
# a different data generation...
mydata <- SimulateCVR(family = "cox", n = n, rank = 4, rank.common=2,
                      p1 = 200, p2 = 300, pnz = 30,  
                      beta = c(0,1,2,0, 0,0), cc = c(0.9, 0.3), 
                      param.cox=list(weibull_scale=200, weibull_shape=20, event.rate=0.3, tie=F), 
                      standardization = F)
X1 <- mydata$X1 # 400*200
X2 <- mydata$X2 # 400*300
Xlist <- list(X1 = X1, X2 = X2); 
Y <- mydata$y[,1]
event = mydata$y[,2]
# CVR fit, fix rank = 4, tune eta and lambda   
out_cvr <- CVR(Y, Xlist, event=event, rankseq= 4, neta= 5, nlam=25, family = "cox", nfold = 5)
# out_cvr$solution$W[[1]]
ff = out_cvr$cvout$cvr.fit # refit
pred1 = PredCVR(mydata$ytest[,1], mydata$X1test, mydata$X2test, mydata$ytest[,2], 
                ff$W[[1]], ff$W[[2]], alpha = ff$alpha, beta =ff$beta, 
                Y=mydata$y[,1], X1=mydata$X1, X2=mydata$X2, event=mydata$y[,2], 
                family = "cox", refit = T, combine = T, type.measure = 'UnoC')

# predict test data, output Uno' concordance  
ff = fit0$cvout$cvr.fit
pred1 = PredCVR(dd0$ytest[,1], dd0$X1test, dd0$X2test, dd0$ytest[,2], 
                ff$W[[1]], ff$W[[2]], alpha = ff$alpha, beta =ff$beta, 
                Y=dd0$y[,1], X1=dd0$X1, X2=dd0$X2, event=dd0$y[,2], 
                family = "cox", refit = T, combine = T, type.measure = 'UnoC')




######### CVR sequential estimation method (sCVR) 
source("R/scvr-engine.R")

offsetk=list(W1k=NULL, W2k=NULL, beta1k=NULL, beta2k=NULL) 

# sCVR fit: 1st layer, tune eta and lambda
seq1 = TuneCVR_seq(Y, X1, X2, event, offsetk = offsetk, etaseq=seq(0.1,0.5,0.1), family="cox", warm=F)  #  1 min 
# seq1$msparse # sparsity
# plot(seq1$cvm[1,])   # mean cross-validated concordance, neta*nlambda
# lines(seq1$cvm[2,])
# lines(seq1$cvm[3,])
# lines(seq1$cvm[4,])
# lines(seq1$cvm[5,])
fit1 = seq1$refit$fit    # refit result
# sum(fit1$w1!=0); sum(fit1$w2!=0) # selected variables
X1w1 = as.numeric(X1%*%fit1$w1 )
X2w2 = as.numeric(X2%*%fit1$w2 ) 
# predict test data using 1st layer of extracted CVs
preds1 = PredCVR(Ytest, X1test, X2test, eventtest, W1=as.matrix(fit1$w1), W2=as.matrix(fit1$w2), 
                alpha = 0, beta = fit1$beta, Y=Y, X1=X1, X2=X2, event=event, family = "cox") # 0.704

# sCVR fit: 2nd layer, tune eta and lambda 
# eta smaller -> more weight on prediction 
offset1 = fit1$offsetk # previous layer is included as offset
seq2 = TuneCVR_seq(Y, X1, X2, event, offsetk=offset1, etaseq=seq1$etaseq/2, family="cox", warm=F)  # 5 min 
fit2 = seq2$refit$fit
# sum(fit2$w1!=0); sum(fit2$w2!=0) #  
X1w12 = as.numeric(X1%*%fit2$w1)
X2w22 = as.numeric(X2%*%fit2$w2)   
W1=as.matrix(cbind(fit1$w1,fit2$w1))
W2=as.matrix(cbind(fit1$w2,fit2$w2))
preds2 = PredCVR(Ytest, X1test, X2test, eventtest, W1=W1, W2=W2, 
                 alpha = 0, beta = fit2$beta, Y=Y, X1=X1, X2=X2, event=event, family = "cox") # 0.747


