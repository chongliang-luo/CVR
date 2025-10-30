
## TODO: 
# adjust for covariate
# prediction: y~CV1+CV2
# rcpp_coxph_deriv_lp also in cvrsolver.cpp? 
# no need to use UnoC, just concordance(timewe='')
 
# logl
# vnorm
# solve_w
# cvrsolver_seq_1
# TuneCVR_seq


# solve_v <- function(eta, hv, hr, hessian_D, R2k, mu){ 
#   # M = (diag((1-eta)*hessian_D + eta + hv) + hr*R2k%*%t(R2k))/2
#   # mu = eta*v1 - (1-eta)*gradient + (1-eta)*diag(hessian_D)*v2+hv*v2+Gv - R2k*Gr
#   # find c such that ||solve(M+c*I_n, mu)|| = 1 
#   n = length(hessian_D)
#   k = ncol(R2k)
#   
#   Mcinv_mu = function(cc){
#     dinv = 1/((1-eta)*hessian_D + eta+hv+cc)  # n
#     # Mcinv = 2*(diag(dinv)-hr*(dinv*R2k)%*%solve(diag(k)+hr*t(R2k)%*%(R2k*dinv), t(dinv*R2k))) # n*n
#     # return(Mcinv%*%mu)  # n*1
#     tt = 2*(dinv*mu-hr*(dinv*R2k)%*%solve(diag(k)+hr*t(R2k)%*%(R2k*dinv), colSums(mu*dinv*R2k)))
#     return(tt)  # n*1
#   }
#   Minvmu = Mcinv_mu(0)
#   
#   ## check RcppNumerical to do optim in rcpp: optim_lbfgs()
#   # https://cran.r-project.org/web/packages/RcppNumerical/vignettes/introduction.html#numerical-optimization
#   chat = optimize(function(cc) (sum(Mcinv_mu(cc)^2) - 1)^2, c(-10,10))$minimum 
#   return(list(Minvmu=Minvmu, Minvmu_norm = sqrt(sum(Minvmu^2)), Mcinvmu = Mcinv_mu(chat), chat = chat))
# } 

## rcpp functions for cox reg
Rcpp::sourceCpp("src/rcpp_coxph.cpp")
Rcpp::sourceCpp("src/rcpp_coxph_deriv_lp.cpp")

## neg log (partial) likelihood fun of glm or coxph
logl = function(y, lp, event, family){
  if(family=='gaussian'){
    ll = sum((y-lp)^2)
  } else if(family=='binomial'){
    ll = -sum(y*lp - log(1+exp(lp)))
  } else if(family=='poisson'){
    ll = sum(exp(lp) - y*lp )
  } else if(family=='cox'){
    # ll = coxph_lpl(event, lp)
    ll = rcpp_coxph_logL(c(1), y, event, as.matrix(lp))
  }
  return(ll)
}

## vector norm
vnorm <- function(a) sqrt(sum(a^2))

## sequential fitting CVR
# for the (k+1)-th layer
# 1/2*eta*||X1w1-X2w2||^2 + (1-eta)*\sum_j=1^2 logL(Xjwj,beta|s^k) + \sum_j=1^2 P(wj,\lambda_j),
# s.t. wj^TXj^TXjwj = 1, Rj^k^TXjwj = 0_k, j=1,2
# ALM augmented lagrangian multiplier:
# for fixed w1, 
# w2 = argmin 1/2*eta*||v1hat-v2||^2 + (1-eta)*logL(v2,beta2|s^k) + P(w2,\lambda_2) 
#               + ..*||v2-X2w2-..||^2 + ..*||R2^k^Tv2 - ..||^2
#  s.t. v2=X2w2, v2^Tv2=1, R2^k^Tv2 = 0_k
solve_w <- function(w1hat, w2hat, X1, X2, Y, event, beta2, os2, eta, lam2, R2k, family='cox',
                    control=list(maxit=3, tol=1e-4,hh=1.01)){
  # for fixed w1, beta1, find w2 by ALM
  # w2 = argmin 1/2*eta*||v1hat-v2||^2 + (1-eta)*logL(v2,beta2|os) + P(w2, lam2) 
  #               + hv/2*||v2-X2w2-Gv/hv||^2 + hr/2*||R2k^Tv2-Gr/hr||^2
  #  s.t. v2=X2w2, v2^Tv2=1, R2k^Tv2 = 0_k
  ## for Cox reg, time Y needs to be ordered!!!
  # beta does not include intercept
  # switch 1 and 2 for finding w1 given w2 ...
  n = nrow(R2k)
  k = ncol(R2k)
  v1 = X1%*%w1hat
  v2 = X2%*%w2hat
  Gv = rep(0,n)
  Gr = rep(0,k)
  hv = hr = 1
  Obj = rep(NA, control$maxit+1)
  Obj[1] = 1/2*eta*sum((v1-v2)^2) + (1-eta)*logl(Y,as.vector(v2*beta2+os2),event,family)/n + lam2*sum(abs(w2hat)) 
  for(iter in 1:control$maxit){ 
    # update w2: hv/2*||X2w2-v2+Gv/hv||^2 + P(w2,lam2) 
    # cat(summary(v2-Gv/hv))
    # cat(v2)
    # cat(Gv)
    w2 = glmnet(X2, v2-Gv/hv, family='gaussian', lambda=lam2/hv/n)$beta # check: + t(R2k)%*%X2 ?
    
    # Do not use: update v2: min(v2^TM2v2 - mu^Tv2) 
    # loglik / lp: gradient vector, and hessian vector (diagonal)
    lp2 = as.vector(v2)*beta2 + os2
    elp2 = exp(lp2)
    if(family=='gaussian'){
      dz = -(Y - lp2) * beta2
      Hz = rep(beta2^2, n)
    } else if(family=='binomial'){
      dz = -(Y - elp2/(1+elp2)) * beta2 
      Hz = elp2/(1+elp2)^2 * beta2^2 
    } else if(family=='poisson'){
      dz = -(Y - elp2) * beta2 
      Hz = elp2 * beta2^2 
    } else if(family=='cox'){
      # d12 = rcpp_coxph_deriv_lp(event, as.vector(os2), beta2, as.vector(v2))
      # dz = d12$gradient
      # Hz = d12$hessianD
      tt = rcpp_coxph_deriv_lp(c(1), Y, event, lp2) # Y need to be sorted!!
      dz = tt[1:n]
      Hz = tt[-c(1:n)]
    }
    
    mu = eta*v1 - (1-eta)*dz/n + (1-eta)*Hz/n*v2 + hv*v2 + Gv - R2k%*%Gr
    # tt = solve_v(eta, hv, hr, Hz/n, R2k, mu)
    # v2 = tt$Mcinvmu 
    # a standard orthogonal Procrustes problem with rank=1
    v2 = mu / vnorm(mu)
    
    
    # update Gv, Gr
    Gv = Gv + hv*(X2%*%w2 - v2)  # n x 1
    Gr = Gr + hr*(t(R2k)%*%v2)   # k x 1
    
    # update hv, hr
    hv = hv*control$hh
    hr = hr*control$hh
    
    # obj 
    Obj[iter+1] = 1/2*eta*sum((v1-v2)^2) + (1-eta)*logl(Y,as.vector(v2*beta2+os2),event,family)/n + lam2*sum(abs(w2)) 
    if(is.na(Obj[iter+1])) break  # need debug, d12$hessianD might be singular?
    if(abs(Obj[iter]-Obj[iter+1])/abs(Obj[iter]) < control$tol) break
  } # end for  
  
  return(w2)
}


cvrsolver_seq_1 <- function(Y, X1, X2, event, offsetk,   # =list(W1k, W2k, beta1k, beta2k),  
                            eta, Lam, family, wini=NULL, control=list(glmnet.alpha=1, maxit=20, tol=1e-4, hh=1.01))
{
  # X1=scale(X1)
  # X2=scale(X2)
  # offsetk
  # eta=0.1
  # Lam=c(1,1)
  # family='cox'
  # X1=Xlist[[1]]
  # X2=Xlist[[2]]
  n = nrow(X1)
  p1 = ncol(X1)
  p2 = ncol(X2)
  k=length(offsetk$beta1k )
  if(is.null(offsetk$beta1k)){ # the 1st seq
    offsetk$W1k=matrix(0,p1,1)
    offsetk$W2k=matrix(0,p2,1)
    offsetk$beta1k=offsetk$beta2k=0
  }
  R1k = X1%*%offsetk$W1k
  R2k = X2%*%offsetk$W2k 
  lam1 = Lam[1]
  lam2 = Lam[2]
  os1 = as.vector(R1k%*%offsetk$beta1k)
  os2 = as.vector(R2k%*%offsetk$beta2k)
  
  # init 
  if(is.null(wini)){ 
    # ridge as init
    if(family=='cox'){
      tt = cv.glmnet(X1,Surv(Y,event),offset=os1, alpha=0,nlambda=20, nfolds = 5, family='cox')
      w1.ini = as.vector(coef(tt, s = 'lambda.min'))
      tt = cv.glmnet(X2,Surv(Y,event),offset=os2, alpha=0,nlambda=20, nfolds = 5, family='cox')
      w2.ini = as.vector(coef(tt, s = 'lambda.min'))
      # w1.ini = coxph(Surv(Y,event)~X1+offset(os1))$coef
      # w2.ini = coxph(Surv(Y,event)~X2+offset(os2))$coef
      alpha.ini=c(0,0)
    } else{
      tt = cv.glmnet(X1,Y,offset=os1, alpha=0,nlambda=20, nfolds = 5, family=family)
      w1.ini = as.vector(coef(tt, s = 'lambda.min'))
      tt = cv.glmnet(X2,Y,offset=os2, alpha=0,nlambda=20, nfolds = 5, family=family)
      w2.ini = as.vector(coef(tt, s = 'lambda.min'))
      # w1.ini = glm(Y~X1, offset=os1, family=family)$coef
      # w2.ini = glm(Y~X2, offset=os2, family=family)$coef 
      alpha.ini=c(w1.ini[1], w2.ini[1])
      w1.ini = w1.ini[-1]
      w2.ini = w2.ini[-1]
      # cat(w1.ini)
      # cat(w2.ini)
    }
    beta.ini = c(sqrt(sum((X1%*%w1.ini)^2)), sqrt(sum((X2%*%w2.ini)^2)))
    w1.ini = w1.ini/beta.ini[1]
    w2.ini = w2.ini/beta.ini[2]
    lp1.ini = X1%*%w1.ini
    lp2.ini = X2%*%w2.ini
  }else{
    w1.ini = wini[[1]]
    w2.ini = wini[[2]]
    lp1.ini = X1%*%w1.ini
    lp2.ini = X2%*%w2.ini
    if(family=='cox'){
      beta.ini = c(coxph(Surv(Y, event)~lp1.ini+offset(os1))$coef, 
                   coxph(Surv(Y, event)~lp2.ini+offset(os2))$coef)
      alpha.ini=c(0,0)
    }else{
      beta.ini = c(glm(Y~lp1.ini, offset=os1, family=family)$coef,
                   glm(Y~lp2.ini, offset=os2, family=family)$coef)
      alpha.ini=beta.ini[c(1,3)]
      beta.ini=beta.ini[c(2,4)]
    }
  }
  
  w1n = w1 = w1.ini
  w2n = w2 = w2.ini
  betan = beta = beta.ini
  alphan = alpha = alpha.ini
  Obj = matrix(NA, control$maxit+1, 6)
  ## use Y~X1W1 and Y~X2W2 is from CVR
  ## TO-DO: check if using Y~X1W1+X2W2 is better in the sCVR
  Obj[1,1:5] = c(1/2*sum((lp1.ini-lp2.ini)^2), 
                 logl(Y,alphan[1]+lp1.ini*betan[1]+os1,event,family)/n, 
                 logl(Y,alphan[2]+lp2.ini*betan[2]+os2,event,family)/n,
                 sum(abs(w1n)), 
                 sum(abs(w2n))) 
  Obj[1,6] = sum(Obj[1,1:5] * c(eta, 1-eta, 1-eta, lam1, lam2))
  
  for(iter in 1:control$maxit){ 
    # cat('iter=', iter)
    # cat(Obj[iter], ',')
    # for fixed w1, beta1, find w2
    w2n = solve_w(w1n, w2n, X1, X2, Y, event, betan[2], os2, eta, lam2, R2k, family)
    
    # for fixed w2, beta2, find w1
    w1n = solve_w(w2n, w1n, X2, X1, Y, event, betan[1], os1, eta, lam1, R1k, family)
    if(sum(abs(w2n))==0 | sum(abs(w1n))==0) break
    
    # for fixed w1, w2, find beta 
    lp1 = as.vector(X1%*%w1n)
    lp2 = as.vector(X2%*%w2n)
    if(family=='cox'){
      betan = c(coxph(Surv(Y, event)~lp1+offset(os1))$coef,
                coxph(Surv(Y, event)~lp2+offset(os2))$coef)
      alphan = c(0,0)
    }else{ 
      tt = c(glm(Y~lp1, offset=os1, family=family)$coef,
             glm(Y~lp2, offset=os2, family=family)$coef)
      betan = tt[c(2,4)] 
      alphan = tt[c(1,3)]  # need work
    }
    
    # obj 
    Obj[iter+1,1:5] = c(1/2*sum((lp1-lp2)^2), 
                        logl(Y,alphan[1]+lp1*betan[1]+os1,event,family)/n, 
                        logl(Y,alphan[2]+lp2*betan[2]+os2,event,family)/n,
                        sum(abs(w1n)),
                        sum(abs(w2n)))   
    # using lasso penalty for now, can be enet (glmnet.alpha=0.5)
    Obj[iter+1,6] = sum(Obj[iter+1,1:5] * c(eta, 1-eta, 1-eta, lam1, lam2))
    # cat(Obj[1:(iter+1)], '...', abs(Obj[iter+1]-Obj[iter])/abs(Obj[iter]), '\n')
    if(abs(Obj[iter+1,6]-Obj[iter,6])/abs(Obj[iter,6])<control$tol) break
  } # end for  
  
  # update offsetk for next seq estimation...
  if(k==0){
    offsetk$W1k = w1n 
    offsetk$W2k = w2n
    offsetk$beta1k = betan[1]
    offsetk$beta2k = betan[2]
  }else{
    offsetk$W1k = cbind(offsetk$W1k, w1n)
    offsetk$W2k = cbind(offsetk$W2k, w2n)
    offsetk$beta1k = c(offsetk$beta1k, betan[1])
    offsetk$beta2k = c(offsetk$beta2k, betan[2])
  }
  
  return(list(w1=w1n, w2=w2n, beta=betan, alpha=alphan, wini=list(w1.ini, w2.ini), 
              iter=iter, Obj=Obj, offsetk=offsetk))
}

# PredCVR_seq <- function()

## cross validate to select eta and lambda (one-layer)
TuneCVR_seq <- function(Y, X1, X2, event=NULL, offsetk, 
                        etaseq=NULL, Lamseq = NULL, neta=5, nlam=20,
                        family = c("gaussian", "binomial", "poisson", "cox"), 
                        wini = NULL, penalty = c("L1", "enet"),   # only L1 now # opts, 
                        nfold = 5, foldid = NULL, type.measure = NULL, warm=T, trace=F){    
  ##  cross-validate to select lam's and eta's
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  Y <- as.matrix(Y)
  if(family=='cox') event=as.matrix(event)
  n <- dim(Y)[1]
  K <- 2 # length(Xlist)
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  p1 <- ncol(X1) 
  p2 <- ncol(X2)
  
  k=length(offsetk$beta1k )
  if(is.null(offsetk$beta1k)){ # the 1st seq
    offsetk$W1k=matrix(0,p1,1)
    offsetk$W2k=matrix(0,p2,1)
    offsetk$beta1k=offsetk$beta2k=0
  }
  R1k = X1%*%offsetk$W1k
  R2k = X2%*%offsetk$W2k  
  os1 = as.vector(R1k%*%offsetk$beta1k)
  os2 = as.vector(R2k%*%offsetk$beta2k)
  
  # init 
  if(is.null(wini)){ 
    if(family=='cox'){
      # w1.ini = coxph(Surv(Y,event)~X1+offset(os1))$coef
      # w2.ini = coxph(Surv(Y,event)~X2+offset(os2))$coef
      # a quick glmnet
      tt = cv.glmnet(X1,Surv(Y,event),offset=os1, alpha=0,nlambda=20, nfolds = 5, family='cox')
      w1.ini = as.vector(coef(tt, s = 'lambda.min'))
      tt = cv.glmnet(X2,Surv(Y,event),offset=os2, alpha=0,nlambda=20, nfolds = 5, family='cox')
      w2.ini = as.vector(coef(tt, s = 'lambda.min'))
      alpha.ini=c(0,0)
    } else{
      # w1.ini = glm(Y~X1, offset=os1, family=family)$coef
      # w2.ini = glm(Y~X2, offset=os2, family=family)$coef 
      tt = cv.glmnet(X1,Y,offset=os1, alpha=0,nlambda=20, nfolds = 5, family=family)
      w1.ini = as.vector(coef(tt, s = 'lambda.min'))
      tt = cv.glmnet(X2,Y,offset=os2, alpha=0,nlambda=20, nfolds = 5, family=family)
      w2.ini = as.vector(coef(tt, s = 'lambda.min'))
      alpha.ini=c(w1.ini[1], w2.ini[1])
      w1.ini = w1.ini[-1]
      w2.ini = w2.ini[-1]
    }
    beta.ini = c(sqrt(sum((X1%*%w1.ini)^2)), sqrt(sum((X2%*%w2.ini)^2)))
    w1.ini = w1.ini/beta.ini[1]
    w2.ini = w2.ini/beta.ini[2]
    lp1.ini = X1%*%w1.ini
    lp2.ini = X2%*%w2.ini
    wini = list(w1.ini, w2.ini)
  }else{
    w1.ini = wini[[1]]
    w2.ini = wini[[2]]
    lp1.ini = X1%*%w1.ini
    lp2.ini = X2%*%w2.ini
    if(family=='cox'){
      beta.ini = c(coxph(Surv(Y, event)~lp1.ini+offset(os1))$coef, 
                   coxph(Surv(Y, event)~lp2.ini+offset(os2))$coef)
      alpha.ini=c(0,0)
    }else{
      beta.ini = c(glm(Y~lp1.ini, offset=os1, family=family)$coef,
                   glm(Y~lp2.ini, offset=os2, family=family)$coef)
      alpha.ini=beta.ini[c(1,3)]
      beta.ini=beta.ini[c(2,4)]
    }
  }
  
  if (is.null(type.measure)){
    if (family == "gaussian") {
      type.measure <- "mse";
    } else if (family == "binomial") {
      if (n / nfold > 10 & is.null(type.measure)) type.measure <- "auc";
      if (n / nfold < 11 & is.null(type.measure)) type.measure <- "deviance";   
    } else if (family == "poisson") {
      type.measure <- "deviance";   
    } else if (family == "cox") {
      type.measure <- "UnoC";   
    }
  }
  
  # warm <- TRUE;
  # opts$nrank <- rank
  # opts$family <- family
  # opts$penalty <- penalty
  
  if (is.null(etaseq)) {
    if(is.null(neta)) neta <- 5
    # etaseq <- seq(0,0.4,by=0.1)
    etaseq <- c(0.01, seq(0.1,0.4,by=0.1))
  } else {neta = length(etaseq)}
  if (is.null(Lamseq)) {
    if (is.null(nlam)) nlam <- 20
    lamseq <- 10^(seq(-2, 0.5, len = nlam)) # * 5
    Lamseq <- cbind(lamseq, lamseq * p2 / p1)
  } else {nlam = nrow(Lamseq)} 
  
  pred <- array(NA, c(neta, nlam, nfold))
  sparse <- array(NA, c(neta, nlam, nfold))    # monitor sparsity
  niter <- array(NA, c(neta, nlam, nfold)) 
  
  if (is.null(foldid)) {  
    tmp <- ceiling(c(1:n) / (n / nfold));
    foldid <- sample(tmp, n)   
  } 
  
  for (i in 1:nfold) {
    X1train <- X1[foldid != i, ];
    X2train <- X2[foldid != i, ];
    # Xtrain  <- list(X1 = X1train, X2 = X2train);
    Ytrain  <- as.matrix(Y[foldid != i,  ]);
    X1test  <- X1[foldid == i, ];
    X2test  <- X2[foldid == i, ];
    Ytest   <- as.matrix(Y[foldid == i, ]);
    eventtrain = eventtest = matrix(NA)
    if(family=='cox'){
      eventtrain = event[foldid != i,  ]
      eventtest = event[foldid == i,  ]
    } 
    # wini = SparseCCA(X1_[foldid == i, ], X2_[foldid == i, ], 1, 0.7, 0.7)
    
    for(ie in 1:neta){
      for (j in 1:nlam){  
        wini.j = lapply(wini, function(a) a*nfold/(nfold-1)) # X1w1 need to be norm-1
        if(warm==T & j>1) wini.j = fiti[1:2]
        fiti = cvrsolver_seq_1(Ytrain, X1train, X2train, eventtrain, offsetk, etaseq[ie], Lamseq[j,], family, wini.j)
        # refit to get coef       
        lptrain1 =  as.matrix(X1train%*%fiti$offsetk$W1k)   
        lptrain2 =  as.matrix(X2train%*%fiti$offsetk$W2k) 
        if (family == "cox") {   
          refit = coxph(Surv(Ytrain, eventtrain)~ lptrain1+lptrain2)$coef
        } else {
          refit <- coef(glm(Ytrain ~ lptrain1+lptrain2, family = family))[-1] 
        }
        
        if(any(!is.na(refit))==F) break ## all are degenerated
        refit[is.na(refit)] = 0  # when some lptrain is all zero
        lptest = as.vector(cbind(X1test%*%fiti$offsetk$W1k, X2test%*%fiti$offsetk$W2k) %*% refit)
        if(family=='cox')
          pred[ie, j, i] <- UnoC(Surv(Ytrain, eventtrain), Surv(Ytest, eventtest), lptest)  
        else if(family=='gaussian')
          pred[ie, j, i] <- 1-summary(lm(Ytest~lptest))$r.squ
        
        ## sparse = proportion of nonzero elements among W1 and W2
        sparse[ie, j, i] <- (sum(fiti$w1 != 0) + sum(fiti$w2 != 0)) / (p1 + p2)
        niter[ie, j, i] <- fiti$iter
        if(trace==T) cat('.')
      }
      if(trace==T) cat('\n')
    }
    if(trace==T) cat('fold=', i, as.character(Sys.time()))  # 20 min each
  }
  
  mpred <- apply(pred, c(1,2), mean, na.rm=T) # 
  msparse <- apply(sparse, c(1,2), mean, na.rm=T)
  # spthresh <- opts$spthresh
  if (type.measure == "mse" | type.measure == "deviance") {
    ## search the best lam only among sparse models
    ieta = which.min(apply(mpred,1,min,na.rm=T))
    ilam = which.min(apply(mpred,2,min,na.rm=T))
    etahat = etaseq[ieta]
    Lamhat = Lamseq[ilam, ]
  } else if (type.measure %in% c("auc", "UnoC", "C")) {
    ieta = which.max(apply(mpred,1,max,na.rm=T))
    ilam = which.max(apply(mpred,2,max,na.rm=T))
    etahat = etaseq[ieta]
    Lamhat = Lamseq[ilam, ]
    # mpred[ieta, ilam]
  }
  
  ## refit all the data with etahat and Lamhat   
  fit = cvrsolver_seq_1(Y, X1, X2, event, offsetk, etahat, Lamhat, family, wini = wini) 
  ## refit to get coef, for prediction
  XW <- as.matrix(cbind(X1 %*% fit$offsetk$W1k, X2 %*% fit$offsetk$W2k)  )
  if (family == "cox") {   
    ab <- coef(coxph(Surv(Y, event) ~ XW)); 
    ab[is.na(ab)] <- 0
    alpha.refit <- 0
    beta.refit <- ab
  } else {
    ab <- coef(glm(Y ~ XW, family = family))
    ab[is.na(ab)] <- 0
    alpha.refit <- ab[1]
    beta.refit <- ab[-1]
  }    
  
  refit <- list(alpha = alpha.refit, beta = beta.refit, fit = fit)
  cv.out <- list(Lamseq = Lamseq, etaseq = etaseq, # rank = rank, 
                 cverror = pred, cvm = mpred, msparse = msparse, 
                 ieta=ieta, ilam=ilam, Lamhat = Lamhat, etahat=etahat,  
                 refit = refit, type.measure = type.measure)
  class(cv.out) <- "TuneCVR_seq"
  return(cv.out) 
}



# CVR_seq <- function()
# TuneCVR_seq <- function(Y, X1, X2, event=NULL, etaseq, Lamseq = NULL, offsetk, 
#                         family = c("gaussian", "binomial", "poisson", "cox"), 
#                         Wini = NULL, penalty = c("L1", "enet"),   # opts, 
#                         nfold = 5, foldid = NULL, type.measure = NULL, warm=T)

