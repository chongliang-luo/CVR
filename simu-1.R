

# N=500, p1=60 blood biomarkers, p2=3000 proteins 
# training data: X1, X2, Y, event
# testing data: X1test, X2test, Ytest, eventtest


require(sup.r.jive)
# remotes::install_github("enorthrop/sup.r.jive")
# install.packages("r.jive")
require(r.jive)
require(PMA)
require(survival)
require(survAUC)
source("R/scvr-engine.R") # sequential CVR

###  X1 X2 have common+individual canonical variates, both predict Y 
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
Ytest=mydata$ytest[,1]
X1test=mydata$X1test 
X2test=mydata$X2test 
eventtest=mydata$ytest[,2]


#(1)# do separate prediction coxnet(Y~X1), coxnet(Y~X2)


#(2)# do combined prediction coxnet(Y~X1+X2)


#(3)# do 2-step: first sparse CCA(X1,X2), then Y~X1W1+X2W2 (correlated, but may not be predictive)
r=2
scca = SparseCCA(X1, X2, rank=r)
b.scca = coxph(Surv(Y, event)~cbind(X1%*%scca$W1, X2%*%scca$W2))$coef
pred_scca = UnoC(Surv(Y, event), Surv(Ytest, eventtest),
                 lpnew=cbind(X1test%*%scca$W1, X2test%*%scca$W2)%*%b.scca)  # lpnew is test data linear predictor


#(4)# JIVE (Joint and Individual Variation Explained), factorization based
# JIVE: https://pmc.ncbi.nlm.nih.gov/articles/PMC3671601/ 
j1 <- jive(list(X1=t(X1),X2=t(X2)), rankJ = 1, rankA = c(1,1),method = 'given') # 

pred.jive = jive.predict(list(X1=t(X1),X2=t(X2)), j1)
Uj = t(pred.jive$joint.scores)
Ui1= t(pred.jive$indiv.scores[[1]])
Ui2= t(pred.jive$indiv.scores[[2]])
cox_jive = coxph(Surv(Y,event) ~ Uj+Ui1+Ui2)
pred.jive.te = jive.predict(list(X1=t(X1test),X2=t(X2test)), j1)
Utej = t(pred.jive.te$joint.scores)
Utei1= t(pred.jive.te$indiv.scores[[1]])
Utei2= t(pred.jive.te$indiv.scores[[2]])
pred_jive = UnoC(Surv(Y, event), Surv(Ytest, eventtest),
                 lpnew=cbind(Utej, Utei1, Utei2)%*%cox_jive$coef) 

# supervised JIVE: can't handle time-to-event outcome
# sJIVE: https://pmc.ncbi.nlm.nih.gov/articles/PMC9481062/ 
# sj1 = sJIVE(list(t(X1), t(X2)), Y, rankJ = 1, rankA = c(1,1),method = 'given')
# pred.sjive = predict(sj1, list(t(X1),t(X2)))
# Uj = t(pred.sjive$Sj)
# Ui1= t(pred.sjive$Si[[1]])
# Ui2= t(pred.sjive$Si[[2]])
# pred.sjive.te = predict(sj1, list(t(X1te),t(X2te)))
# Utej = t(pred.sjive.te$Sj)
# Utei1= t(pred.sjive.te$Si[[1]])
# Utei2= t(pred.sjive.te$Si[[2]]) 


#(5)# CVR fit, fix rank = 4, tune eta and lambda 
out_cvr <- CVR(Y, Xlist, event=event, rankseq=4, neta= 5, nlam = 25, family = "cox", nfold = 5)
# out_cvr$solution$W[[1]]
ff = out_cvr$cvout$cvr.fit # refit
round(cor(X1 %*% ff$W[[1]], X2 %*% ff$W[[2]]), 2)

# predict test data, output Uno's concordance 
pred1 = PredCVR(mydata$ytest[,1], mydata$X1test, mydata$X2test, mydata$ytest[,2], 
                ff$W[[1]], ff$W[[2]], alpha = ff$alpha, beta =ff$beta, 
                Y=mydata$y[,1], X1=mydata$X1, X2=mydata$X2, event=mydata$y[,2], 
                family = "cox", refit = T, combine = T, type.measure = 'UnoC')
pred1 # 0.846


#(6)# sequential CVR (sCVR) 
offsetk=list(W1k=NULL, W2k=NULL, beta1k=NULL, beta2k=NULL)  
# sCVR fit: 1st layer, tune eta and lambda
Lamseq = seq1$Lamseq; 
Lamseq[,2] = Lamseq[,1] * 5;
seq1 = TuneCVR_seq(Y, X1, X2, event, offsetk = offsetk, etaseq=seq(0.1,0.5,0.1), Lamseq=Lamseq, family="cox", warm=F)  #  1 min 

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
                 alpha = 0, beta = fit2$beta, Y=Y, X1=X1, X2=X2, event=event, family = "cox") # 


# DIABLO
# https://mixomics.org/mixdiablo/