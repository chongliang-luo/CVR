## PMA package is needed to install CVR
## PMA depends on impute, which is removed from CRAN to bioconuctor
## Use the following to install impute 
## (try http:// if https:// URLs are not supported)
## source("https://bioconductor.org/biocLite.R")
## biocLite("impute")

## useDynLib(CVR, .registration = TRUE)

## row-wise norm of W
norm21 <-function(W) apply(W, 1, function(a) sqrt(sum(a^2)))

## projection matrix of U
Pj <- function(U) U%*%solve(t(U)%*%U, t(U))

## distance of Pj(U) and Pj(X1%*%W1)
Dspace <- function(U, X1W1) sum((Pj(U)-Pj(X1W1))^2)/(dim(as.matrix(U))[1]^2)

## precision, recall
PrecRec <- function(W, V, wise="row"){
  ## Precision=TP/(TP+FP) and Recall=TP/(TP+FN), 
  ## the larger the better for both
  ## V is the true sparse matrix
  if(dim(W)[1]!=dim(V)[1]) print(paste("Row dims don't match:",dim(W)[1],dim(V)[1]))
  if(sum(abs(W))==0) return(c(precision=0, recall=0))
  if(wise=="row"){
    w = apply(W, 1, function(a) sqrt(sum(a^2)))
    v = apply(V, 1, function(a) sqrt(sum(a^2)))
    tp = sum((w!=0) * (v!=0))
    fp = sum((w!=0) * (v==0))
    fn = sum((w==0) * (v!=0))
    
  }else if(wise=="entry"){
    if(  dim(W)[2]!=dim(V)[2]) print(paste("Column dims don't match:",dim(W)[2],dim(V)[2]))   
    tp = sum((W!=0) * (V!=0))
    fp = sum((W!=0) * (V==0))
    fn = sum((W==0) * (V!=0))
  }    
  return(c(precision=tp/(tp+fp), recall=tp/(tp+fn), mcr=(fp+fn)/dim(W)[1]))
}


# Rcpp::sourceCpp('coxph.cpp') 
# documentation to be added
# the  gradient / hessian diagonal of Cox PH log-partial likelihood function w.r.t. the linear predictor
# coxph(Surv(time, event) ~ os+lp*theta) # os=offset, lp=linear predictor
coxph_deriv_lp <- function(lp, time, event, os=rep(0,length(time)), theta=1, sd.jitter=0.1, sorted=F){
  # jitter time to get rid of ties
  if(anyDuplicated(time)>0) time0 = time+rnorm(length(time),0,sd.jitter)
  if(sorted==F){
    idx = order(time0)
    time0 = time0[idx]
    event0 = event[idx]
    os0 = os[idx]
    lp0 = lp[idx]
  }
  res = rcpp_coxph_deriv_lp(event=event0, os=os0, theta=theta, lp=lp0)
  return(res)
}


S0inv_weibull <- function(u, p=20, lam=1/200){ # inv surv of weibull 
  lam^(-1) * (-log(1-u))^(1/p) # exp(-(lam*t)^p)
}


## generate data of survival time, using Cox PH with Weibull basehaz 
generate_Cox_data <- function( 
  param = list(N = 1000,
               px = 2,
               beta = c(1:px),
               beta0 = -1,
               S0inv = S0inv_weibull,
               weibull_scale = 200,
               weibull_shape = 20,
               censor = F,
               event.rate = 1,
               tie = F), X = NULL, verbose=F){ 
  N = param$N
  # K = param$K
  px = param$px
  beta = param$beta
  beta0 = param$beta0
  S0inv = param$S0inv
  weibull_scale = param$weibull_scale
  weibull_shape = param$weibull_shape
  censor = param$censor
  event.rate = param$event.rate
  tie = param$tie 
  
  # Covariates X, 
  if(is.null(X)) {
    X <- matrix(runif(N*px, -3, 3), N, px)
    # X <- matrix(rbinom(N*px, 1, 0.5), N, px)
    colnames(X) <- paste0('X_', 1:px)
  }
  
  Xb <- c(X %*% beta + beta0)
  expXb <- exp(Xb)
  # plot(density(Xb))
  
  # generate true event Time with Weibull baseline distribution
  U <- runif(N, 0, 1)
  t_event <- S0inv(U, p=weibull_shape, lam=1/weibull_scale) / expXb^(1/weibull_shape)  
  
  ################################################################################################
  # NOTE: censoring has to be non-informative, 
  # and how to generating censoring time for given event.rate is tricky
  # better not to use: 
  # 1. one censoring time, e.g. quantile(t_event, event.rate)
  # 2. censoring time within a range, e.g. runif(N, 0.1, quantile(t_event, min(0.99,event.rate*2)))
  # suggest use:
  # Weibull r.v. with the same shape but scale weibull_scale_c = weibull_scale * ratio
  # notice ith sub has hazard w scale = weibull_scale / expXb[i]^(1/weibull_shape)
  # do some integral, know: 
  # censoring rate = mean(1/(1+(ratio)^weibull_shape*expXb))
  # thus do 1-d search for opt ratio
  ################################################################################################
  if(event.rate < 1){
    ratio = ((1/(1-event.rate)-1)/mean(expXb))^(1/weibull_shape)
    fn <- function(rr) abs(mean(1/(1+rr^weibull_shape*expXb)) - (1-event.rate) )
    ratio <- optimize(fn, interval = c(0.5, 1.5)*ratio)$minimum
    weibull_scale_c <- weibull_scale * ratio
    # weibull_scale_c <- weibull_scale * ((1/(1-event.rate)-1)/mean(expXb))^(1/weibull_shape)
    
    t_censor <- rweibull(N, weibull_shape, weibull_scale_c)
    # t_censor <- runif(N, 0.1, quantile(t_event, min(0.99,event.rate*2)))
    # t_censor <- quantile(t_event, event.rate)
    ind_event <- t_event < t_censor
    t_surv <- pmin(t_event, t_censor)
  } else {
    ind_event <- rep(TRUE, length(t_event))
    t_surv <- t_event
  }
  if(verbose) cat('censoring rate = ', 1-mean(ind_event), '\n')
  
  if(verbose) plot(density(t_surv))
  
  if(tie) t_surv <- ceiling(t_surv)  # pmax(round(t_surv), 1)  # 
  
  all_data <- list(t_surv=t_surv, 
                   t_event = t_event,
                   ind_event=ind_event,
                   X=X)
  
  return(list(all_data = all_data, param = param))
}




#' @useDynLib CVR
#' @import Rcpp 
#' @import graphics  
#' @import stats  
#' @title  Sparse canonical correlation analysis.
#'  
#' @description Get sparse CCA solutions of X1 and X2. Use \code{CCA} and \code{CCA.permute} from PMA package. 
#'                See PMA package for details.
#'                
#' @usage SparseCCA(X1, X2, rank = 2, penaltyx1 = NULL, penaltyx2 = NULL, 
#'              nperms = 25, ifplot = 0)
#' @param X1,X2 Numeric matrices representing the two sets of covariates. They should have the same number of rows and 
#'              cannot contain missing values.
#' @param rank The number of canonical variate pairs wanted. The default is 2.
#' @param penaltyx1,penaltyx2 Numeric vectors as the penalties to be applied to X1 and X2. The defaults are seq(0.1, 0.7, len = 20).
#'                            See PMA package for details.
#' @param nperms Number of times the data should be permuted. The default is 25. Don't use too small value. See PMA package for details.
#' @param ifplot 0 or 1.  The default is 0 which means don't plot the result. See PMA package for details.
#' @details This function is generally used for tuning the penalties in sparse CCA if \code{penaltyx1} and \code{penaltyx2} are supplied as vectors.
#'          The CCA solution is based on PMD algorithm and the tuning is based on permutation.
#'          The fitted  W1 and W2 are scaled so that the diagonals of W1'X1'X1W1 and W2'X2'X2W2 are all 1's.
#'          
#'          Specifically, if a single value of mild penalty is provided for both \code{penaltyx1} and \code{penaltyx2}, the result can be used as 
#'          initial values of W's in CVR. For instance, with \code{penaltyx1 = 0.7} and \code{penaltyx2 = 0.7},  
#'          the fitted W1 and W2 are only shrinked, but mostly not zero yet. 
#' @return \item{W1}{Loading matrix corresponding to X1. X1*W1 gives the canonical variates from X1.}
#' @return \item{W2}{Loading matrix corresponding to X2. X2*W2 gives the canonical variates from X2.}
#' @author Chongliang Luo, Kun Chen.
#' @references Daniela M. Witten,  Robert Tibshirani and Trevor Hastie (2009) A penalized matrix decomposition, with applications to 
#'                sparse principal components and canonical correlation analysis. Biostatistics 10(3), 515-534.
#'              
#'      Chongliang Luo, Jin Liu, Dipak D. Dey and Kun Chen (2016) Canonical variate regression. 
#'        Biostatistics, doi: 10.1093/biostatistics/kxw001.
#' @seealso \code{\link{CVR}}, \code{\link{CCA}}, \code{\link{CCA.permute}}.   
#' @import PMA
#' @export 
SparseCCA <- function(X1, X2, rank = 2, penaltyx1 = NULL, penaltyx2 = NULL, nperms = 25, ifplot = 0){
    if (nrow(X1) != nrow(X2)) stop("X1 and X2 should have the same row number!")
    if (ncol(X1) == 1 | ncol(X2) == 1) stop("X1 and X2 should have at least 2 columns!")
    if (is.null(penaltyx1) | is.null(penaltyx2)){
      penaltyx1 <- seq(0.1, 0.7, len = 20)
      penaltyx2 <- seq(0.1, 0.7, len = 20)
    }
    if (length(penaltyx1) > 1){
      perm.out <- CCA.permute(X1, X2, "standard", "standard", penaltyx1, penaltyx2, nperms = nperms, trace = F)  
      if (ifplot == 1 ) {
        print(c(perm.out$bestpenaltyx, perm.out$bestpenaltyz))
        plot(perm.out)
      }
      SparseCCA_X12 <- CCA(X1, X2, typex = "standard", typez = "standard", K = rank,
                    penaltyx = perm.out$bestpenaltyx, penaltyz = perm.out$bestpenaltyz,
                    v = perm.out$v.init, trace = FALSE)
    } else {
      SparseCCA_X12 <- CCA(X1, X2, typex = "standard", typez = "standard", K = rank, 
                      penaltyx = penaltyx1, penaltyz = penaltyx2, trace = FALSE)
    }
  D1 <- diag(diag(t(X1 %*% SparseCCA_X12$u) %*% X1 %*% SparseCCA_X12$u)^(-0.5), rank, rank)
  D2 <- diag(diag(t(X2 %*% SparseCCA_X12$v) %*% X2 %*% SparseCCA_X12$v)^(-0.5), rank, rank)
  W1 <- SparseCCA_X12$u %*% D1
  W2 <- SparseCCA_X12$v %*% D2 ## to make X1 %*% W1 orthogonal
  return(list(W1 = W1, W2 = W2))
}



#' @title  Generate simulation data.
#' 
#' @description Generate two sets of covariates and an univariate response driven by several latent factors.
#' 
#' @usage  SimulateCVR(family = c("gaussian", "binomial", "poisson", "cox"),  
#'          n = 100, rank = 4, p1 = 50, p2 = 70, pnz = 10, sigmax = 0.2,   
#'          sigmay = 0.5, beta = c(2, 1, 0, 0), 
#'          param.cox = list(weibull_scale = 200, weibull_shape = 20, 
#'          event.rate = 0.6, tie = F), standardization = TRUE)
#' 
#' @param  family  Type of response. \code{"gaussian"} for continuous response, \code{"binomial"}  
#'                  for binary response, and \code{"poisson"} for Poisson response. The default is \code{"gaussian"}.                             
#' @param  n       Number of rows. The default is 100.          
#' @param  rank   Number of latent factors generating the covariates.  The default is 4.         
#' @param  p1     Number of variables in X1. The default is 50.          
#' @param  p2     Number of variables in X2. The default is 70.          
#' @param  pnz    Number of variables in X1 and X2 related to the signal. The default is 10.           
#' @param  sigmax Standard deviation of normal noise in X1 and X2. The default is 0.2.          
#' @param  sigmay Standard deviation of normal noise in Y. Only used when the response is Gaussian. The default is 0.5.                   
#' @param  beta   Numeric vector, the coefficients used to generate respose from the latent factors. The default is c(2, 1, 0, 0).
#' @param  param.cox parameters for generating survival outcome
#' @param  standardization Logical. If TRUE, standardize X1 and X2 before output. The default is TRUE.
#' 
#' @details The latent factors in U are randomly generated normal vectors, 
#' 
#'          \eqn{X_1 = U*V_1 + \sigma_x*E_1, X_2 = U*V_2 + \sigma_x*E_2, E_1, E_2}{%
#'          X1 = U*V1 + sigmax*E1, X2 = U*V2 + sigmax*E2, E1, E2} are N(0,1) noise matrices.
#'          
#'          The nonzero entries of \eqn{V_1}{%
#'          V1} and \eqn{V_2}{%
#'          V2} are generated from Uniform([-1,-0.5]U[0.5,1]).
#'                   
#'          For Gaussian response, 
#'          
#'          \eqn{y = U*\beta + \sigma_y*e_y, e_y}{%
#'          y = U*\beta + sigmay*ey, ey} is N(0,1) noise vector,
#'          
#'          for binary response, 
#'          
#'          \eqn{y \sim rbinom(n, 1, 1/(1 + \exp(-U*\beta)))},
#'          
#'          and for Poisson response,
#'          
#'          \eqn{y \sim rpois(n, \exp(U*\beta))}.
#'          
#'          See the reference for more details.
#'          
#' @return \item{X1, X2}{The two sets of covariates with dimensions  n*p1 and n*p2 respectively.} 
#' @return \item{y}{The  response vector with length   n.}
#' @return \item{U}{The true latent factor matrix with dimension  n*rank.}
#' @return \item{beta}{The coefficients used to generate response from \code{U}. The length is rank.}
#' @return \item{V1, V2}{The true loading matrices for X1 and X2 with dimensions  p1*rank and p2*rank. The first \code{pnz} rows are nonzero.}
#' @author Chongliang Luo, Kun Chen.
#' @references Chongliang Luo, Jin Liu, Dipak D. Dey and Kun Chen (2016) Canonical variate regression. 
#'               Biostatistics, doi: 10.1093/biostatistics/kxw001.
#'              Jendoubi 2019 BMC. A whitening approach to probabilistic canonical correlation analysis for omics data integration
#' @seealso \code{\link{CVR}}, \code{\link{cvrsolver}}.
#' @export
SimulateCVR <- function(family = c("gaussian", "binomial", "poisson", "cox"), 
                        n = 100, rank = 4, rank.common=2,
                    p1 = 50, p2 = 70, pnz = 10, sigmax = 0.2, sigmay = 0.5, 
                    # beta = c(2, 1, 0, 0), 
                    beta = c(1,0, 1,0, 1,0), # shared Us, ind Us 
                    cc = NULL, # c(0.9, 0.5),  # correlation
                    param.cox =list(weibull_scale = 200, weibull_shape = 20, event.rate = 0.6, tie = F), 
                    standardization = TRUE){
  family <- match.arg(family)
  ## U = svd(U)$u*sqrt(n); #  n by r
  
  V1 <- matrix(0, p1, rank);
  V11 <- matrix(runif(rank * pnz) * (rbinom(rank * pnz, 1, 0.5) * 2-1), ncol = rank);
  V2 <- matrix(0, p2, rank);
  V21 <- matrix(runif(rank * pnz) * (rbinom(rank * pnz, 1, 0.5) * 2-1), ncol = rank);
  V1[1:pnz, ] <- svd(V11)$u;
  V2[1:pnz, ] <- svd(V21)$u; 
  
  if(!is.null(cc)){  # specified canonical cors, generate by shared and individual U's ()
    # [U0:shared, U1:X1 specific, U2:X2 specific]
    rank.ind = rank - rank.common
    U0 <- matrix(rnorm(2 * n * rank.common), 2 * n, rank.common);  
    U1 <- matrix(rnorm(2 * n * rank.ind), 2 * n, rank.ind) 
    U2 <- matrix(rnorm(2 * n * rank.ind), 2 * n, rank.ind) 
    U = cbind(U0,U1,U2)   
    x1 = cbind(U0, U1) %*% diag(sqrt(c(cc, 1-cc))) %*% t(V1)  
    x2 = cbind(U0, U2) %*% diag(sqrt(c(cc, 1-cc))) %*% t(V2)  
    # boxplot(x1)
    # tt=cancor(x1, x2)$cor 
  }else{
    U <- matrix(rnorm(2 * n * rank), 2 * n, rank);
    x1 <- U %*% t(V1); 
    x2 <- U %*% t(V2);
  }
  
  lp <- U %*% beta;
  e1 <- matrix(rnorm(2 * n * p1), 2 * n, p1) * sigmax;
  e2 <- matrix(rnorm(2 * n * p2), 2 * n, p2) * sigmax; 
  
  X1 <- x1[1:n, ] + e1[1:n, ];
  X2 <- x2[1:n, ] + e2[1:n, ];
  X1test <- x1[(n + 1):(2 * n), ] + e1[(n + 1):(2 * n), ];
  X2test <- x2[(n + 1):(2 * n), ] + e2[(n + 1):(2 * n), ];
  Utest = U[(n + 1):(2 * n),]
  U <- U[1:n, ]
  
  # generate outcome
  if (family == "gaussian") {
    e  <- rnorm(2 * n) * sigmay;
    y  <- lp[1:n] + e[1:n];
    ytest  <- lp[(n + 1):(2 * n)] + e[(n + 1):(2 * n)];    
    if (standardization) {
      V1 <- diag(apply(X1, 2, sd)) %*% V1 / sqrt(n - 1)
      V2 <- diag(apply(X2, 2, sd)) %*% V2 / sqrt(n - 1); 
      y <- scale(y); 
      X1 <- scale(X1); 
      X2 <- scale(X2);#/sqrt(n-1)
      ytest <- scale(ytest); 
      X1test <- scale(X1test); 
      X2test <- scale(X2test);#/sqrt(n-1)
    }
    snr.x1 <- sd(x1) / sd(e1)
    snr.x2 <- sd(x2) / sd(e2)
    snr.y <- sd(lp) / sd(e);
    snr <- c(snr.x1, snr.x2, snr.y)
  } else if (family == "binomial") {
    px <- 1 / (1 + exp(-lp));
    y <- rbinom(n, 1, px[1:n])
    ytest <- rbinom(n, 1, px[(n + 1):(2 * n)])
    if (standardization) {
      V1 <- diag(apply(X1, 2, sd)) %*% V1 / sqrt(n - 1); 
      V2 <- diag(apply(X2, 2, sd)) %*% V2 / sqrt(n - 1); 
      X1 <- scale(X1); 
      X2 <- scale(X2);
      X1test <- scale(X1test); 
      X2test <- scale(X2test);
    }
    snr.x1 <- sd(x1) / sd(e1)
    snr.x2 <- sd(x2) / sd(e2)
    snr <- c(snr.x1, snr.x2, 0)
  } else if (family == "poisson") {
    px <- exp(lp);
    y <- rpois(n, px[1:n])
    ytest <- rpois(n, px[(n + 1):(2 * n)])
    if (standardization) {
      V1 <- diag(apply(X1, 2, sd)) %*% V1 / sqrt(n - 1); 
      V2 <- diag(apply(X2, 2, sd)) %*% V2 / sqrt(n - 1); 
      X1 <- scale(X1); 
      X2 <- scale(X2);
      X1test <- scale(X1test); 
      X2test <- scale(X2test);
    }
    snr.x1 <- sd(x1) / sd(e1)
    snr.x2 <- sd(x2) / sd(e2)
    snr <- c(snr.x1, snr.x2, 0)
  } else if (family == "cox") {
    dd = generate_Cox_data(param = list(N = 2*n,
                                       px = rank,
                                       beta = beta,
                                       beta0 = 0,
                                       S0inv = S0inv_weibull,
                                       weibull_scale = param.cox$weibull_scale,
                                       weibull_shape = param.cox$weibull_shape,
                                       censor = T,
                                       event.rate = param.cox$event.rate,
                                       tie = param.cox$tie), X = rbind(U,Utest), verbose=F)
    y = cbind(dd$all_data$t_event, dd$all_data$ind_event)[1:n,]
    ytest = cbind(dd$all_data$t_event, dd$all_data$ind_event)[(n + 1):(2 * n),]
     
    if (standardization) {
      V1 <- diag(apply(X1, 2, sd)) %*% V1 / sqrt(n - 1); 
      V2 <- diag(apply(X2, 2, sd)) %*% V2 / sqrt(n - 1); 
      X1 <- scale(X1); 
      X2 <- scale(X2);
      X1test <- scale(X1test); 
      X2test <- scale(X2test);
    }
    
    idx = order(dd$all_data$t_event[1:n])
    y = y[idx,]
    X1 = X1[idx,]
    X2 = X2[idx,]
    U = U[idx,]
    
    idx = order(dd$all_data$t_event[(n + 1):(2 * n)])
    ytest = ytest[idx,]
    X1test = X1test[idx,]
    X2test = X2test[idx,]
    Utest = Utest[idx,]
  }

  
  ## X1test = X1test, X2test = X2test, ytest = as.matrix(ytest), x1 = x1, x2 = x2, lp = lp, snr = snr
  return(list(X1 = X1, X2 = X2, y = as.matrix(y), U = U, beta = beta, V1 = V1, V2 = V2,
              X1test = X1test, X2test = X2test, ytest = as.matrix(ytest), Utest=Utest))
}


 
## calculate area under ROC curve (AUC) and deviance in binomial case
CalcROC <- function(y, px, ifplot = 0){ 
  ## predict binomial response with given fitted prob px, 
  ## compute AUC of ROC / (neg) deviance (use it instead of AUC, for small sample size)
  ## input: 
  ## y: real response
  ## px: fitted prob seq
  ## output: 
  ## list(spec = spec, sens = sens, auc = auc, dev = dev)
  n <- length(px);
  np <- sum(y == 1);
  nn <- n - np;
  sens <- c();
  spec <- c();
  
  if (np == 0) {
    auc <- sum(px < 0.5) / nn;
  } else if (nn == 0) {
    auc <- sum(px > 0.5) / np;
  } else { 
    for (i in 1:n) {
      sens[i] <- (sum(px > px[i] & y == 1) + 0.5 * sum(px == px[i] & y == 1)) / np;
      spec[i] <- (sum(px < px[i] & y == 0) + 0.5 * sum(px == px[i] & y == 0)) / nn;
    }
  }
  if (ifplot == 1) plot(1 - spec, sens);
  ##auc = trapz(sort(1-spec),sort(sens));
  xx <- sort(1 - spec, index.return = TRUE)
  yy <- c(0, sens[xx$ix], 1)
  xx <- c(0, xx$x, 1)
  auc <- sum((xx[-1] - xx[-(n + 2)]) * (yy[-1] + yy[-(n + 2)])) / 2
  
  dev <- -2 * sum(y * log(px + 0.00001) + (1 - y) * log(1 - px + 0.0001))
  return(list(spec = spec, sens = sens, auc = auc, dev = dev))
}


## calculate the error of testing data (mse, auc or deviance, or Concordance) using fitted W1, W2, alpha, beta
PredCVR <- function(Ytest, X1test, X2test, eventtest, W1, W2, alpha = 0, beta = 0, Y, X1, X2, event, family = "gaussian",
                      refit = TRUE, combine = TRUE, type.measure = NULL){
  ## compute the pred mse/AUC(deviance) of testing data
  X1testW1 <- X1test %*% W1;
  X2testW2 <- X2test %*% W2;
  if(is.null(type.measure)){
    if (family == "gaussian") {
      type.measure <- "mse";
    } else if (family == "binomial") {
      if(length(Ytest) > 10 & is.null(type.measure)) type.measure <- "auc";
      if(length(Ytest) < 11 & is.null(type.measure)) type.measure <- "deviance";   
    } else if (family == "poisson") {
      type.measure <- "deviance";   
    } else if (family == "cox"){
      type.measure = 'UnoC'
    }
  }
  
  if (refit) { 
    X1W1 <- X1 %*% W1;  
    X2W2 <- X2 %*% W2;
    XW <- cbind(X1W1, X2W2)
    XtestW <- cbind(X1testW1, X2testW2)    
    if (family == "gaussian") {   
      ## mse is the average of two predictions from X1testW1 and X2testW2
      if (combine) {
        ab <- coef(lm(Y ~ XW));  
        ab[is.na(ab)] <- 0
        mse <- sum((Ytest - ab[1] - XtestW %*% ab[-1])^2) / length(Ytest)      
      } else {
        ab1 <- coef(lm(Y ~ X1W1))  
        ab2 <- coef(lm(Y ~ X2W2))
        ab1[is.na(ab1)] <- 0; 
        ab2[is.na(ab2)] <- 0;
        mse <- (sum((Ytest - ab1[1] - X1testW1 %*% ab1[-1])^2)
              + sum((Ytest - ab2[1] - X2testW2 %*% ab2[-1])^2)) / length(Ytest) / 2
      }
      if (is.na(mse)) mse <- sum(Ytest^2) / length(Ytest)  # for empty W1 and W2 
      return(mse)
    } else if (family == "binomial") {  
      if (sum(abs(W1)) == 0 | sum(abs(W2)) == 0) {   ## avoid fitting glm with zero predictors 
        if (type.measure == "auc") {
          return(0)
        } else if (type.measure == "deviance") {
          dev <- -2 * length(Ytest) * log(0.00001)
          return(dev)
        }
      } else {
        if (combine) {
          refit12 <- glm(Y ~ XW, family = "binomial")
          ab <- coef(refit12)
          ab[is.na(ab)] <- 0
          px <- 1 / (1 + exp(-XtestW %*% ab[-1] - ab[1]))
          if (type.measure == "auc") {
            auc <- CalcROC(Ytest, px)$auc;
            return(auc)
          } else if (type.measure == "deviance") {
            dev <- CalcROC(Ytest, px)$dev;
            return(dev)
          }
        } else {
            refit1 <- glm(Y ~ X1W1, family = "binomial")
            refit2 <- glm(Y ~ X2W2, family = "binomial")
            ab1 <- coef(refit1) 
            ab2 <- coef(refit2)
            ab1[is.na(ab1)] <- 0 
            ab2[is.na(ab2)] <- 0
            px1 <- 1 / (1 + exp(-X1testW1 %*% ab1[-1] - ab1[1]))
            px2 <- 1 / (1 + exp(-X2testW2 %*% ab2[-1] - ab2[1]))
            if (type.measure == "auc") {
              auc <- (CalcROC(Ytest, px1)$auc + CalcROC(Ytest, px2)$auc) / 2
              return(auc)
            } else if (type.measure == "deviance") {
              dev <- (CalcROC(Ytest, px1)$dev + CalcROC(Ytest, px2)$dev) / 2
              return(dev)
            }
          }
        }      
      } else if (family == "poisson") {  
          if (sum(abs(W1)) == 0 | sum(abs(W2)) == 0) {   
          dev <- -2 * length(Ytest) * log(0.000001)
          return(dev)
        } else {
            if (combine) {
              refit12 <- glm(Y ~ XW, family = "poisson")
              ab <- coef(refit12)
              ab[is.na(ab)] <- 0
              px <-  exp(XtestW %*% ab[-1] + ab[1])
              dev <- -2 * sum(Ytest * log(px) - px)
              return(dev)       
            } else {
              refit1 <- glm(Y ~ X1W1, family = "poisson")
              refit2 <- glm(Y ~ X2W2, family = "poisson")
              ab1 <- coef(refit1); 
              ab2 <- coef(refit2)
              ab1[is.na(ab1)] <- 0; 
              ab2[is.na(ab2)] <- 0
              px1 <- exp(X1testW1 %*% ab1[-1] + ab1[1])
              px2 <- exp(X2testW2 %*% ab2[-1] + ab2[1])
              dev <- -sum(Ytest * log(px1) - px1) - sum(Ytest * log(px2) - px2)
              return(dev)
            }
          }      
      }else if (family == "cox") {  
        if (sum(abs(W1)) == 0 | sum(abs(W2)) == 0) {   
          Cstat <- 0.5
          return(0.5)
        } else {
          if (combine) {
            refit12 <- coxph(Surv(Y, event) ~ .,data=data.frame(Y=Y,event=event,XW=XW))
            ab <- coef(refit12)
            ab[is.na(ab)] <- 0
            lpnew = XtestW%*%ab
            # Censoring-adjusted C-statistic by Uno et al. 2011
            if(type.measure=='UnoC'){
              Cstat <- survAUC::UnoC(Surv(Y, event), Surv(Ytest, eventtest), lpnew) 
            }else if(type.measure=='C'){
              Cstat <- survival::concordance(refit12, newdata=data.frame(Y=Ytest,event=eventtest,XW=XtestW))$concordance 
            }
            return(Cstat)       
          } else {
            refit1 <- coxph(Surv(Y, event) ~ X1W1 )
            refit2 <- coxph(Surv(Y, event) ~ X2W2 )
            ab1 <- coef(refit1); 
            ab2 <- coef(refit2)
            ab1[is.na(ab1)] <- 0; 
            ab2[is.na(ab2)] <- 0
            lpnew1 = X1testW1%*%ab1
            lpnew2 = X2testW2 %*% ab2
            if(type.measure=='UnoC'){
              Cstat1 <- survAUC::UnoC(Surv(Y, event), Surv(Ytest, eventtest), lpnew1)
              Cstat2 <- survAUC::UnoC(Surv(Y, event), Surv(Ytest, eventtest), lpnew2)
            }else if(type.measure=='C'){
              Cstat1 <- survival::concordance(refit1, newdata=data.frame(Y=Ytest,event=eventtest,XW=X1testW1))$concordance
              Cstat1 <- survival::concordance(refit2, newdata=data.frame(Y=Ytest,event=eventtest,XW=X2testW2))$concordance
            }
            return((Cstat1+Cstat1)/2)
          }
        }      
      }
      } else {
          if (family == "gaussian") {   
          lp1 <- alpha + X1testW1 %*% beta
          lp2 <- alpha + X2testW2 %*% beta
          mse <- sum((Ytest - lp1)^2 + (Ytest - lp2)^2) / length(Ytest) / 2
          return(mse)
        } else if (family == "binomial") {   
            if (sum(abs(W1)) == 0 | sum(abs(W2)) == 0) {   
              if (type.measure == "auc") {
              return(0)
            } else if (type.measure == "deviance") {
              dev <- -2 * length(Ytest) * log(0.00001)
              return(dev)
            }
          } else {
              lp1 <- alpha + X1testW1 %*% beta
              lp2 <- alpha + X2testW2 %*% beta
              px1 <- 1 / (1 + exp(-lp1)) 
              px2 <- 1 / (1 + exp(-lp2))
              if (type.measure == "auc") {
                auc <- (CalcROC(Ytest, px1)$auc + CalcROC(Ytest, px2)$auc) / 2
                return(auc)
              } else if (type.measure == "deviance") {
                dev <- (CalcROC(Ytest, px1)$dev + CalcROC(Ytest, px2)$dev) / 2
                return(dev)
              }
            }
          } else if (family == "poisson") {   
              if (sum(abs(W1)) == 0 | sum(abs(W2)) == 0) {   
                dev <- -2 * length(Ytest) * log(0.000001)
                return(dev)
              } else {
                lp1 <- alpha + X1testW1 %*% beta
                lp2 <- alpha + X2testW2 %*% beta
                px1 <- exp(lp1) 
                px2 <- exp(lp2)
                dev <- -sum(Ytest * log(px1) - px1) - sum(Ytest * log(px2) - px2)
                return(dev)
              }
          } else if (family == "cox") {   
            if (sum(abs(W1)) == 0 | sum(abs(W2)) == 0) {   
              Cstat <- 0.5
              return(Cstat)
            } else {
              lp1 <- X1testW1 %*% beta
              lp2 <- X2testW2 %*% beta
              if(type.measure=='UnoC'){
                Cstat1 <- survAUC::UnoC(Surv(Y, event), Surv(Ytest, eventtest), lp1)
                Cstat2 <- survAUC::UnoC(Surv(Y, event), Surv(Ytest, eventtest), lp2)
              }else if(type.measure=='C'){
                Cstat1 <- survival::concordance(coxph(Surv(Ytest,eventtest)~lp1), newdata=data.frame(Y=Ytest,event=eventtest,lp1=lp1))$concordance
                Cstat1 <- survival::concordance(coxph(Surv(Ytest,eventtest)~lp2), newdata=data.frame(Y=Ytest,event=eventtest,lp2=lp2))$concordance
              }
              return((Cstat1+Cstat1)/2)
            }
          }
  }
}

## calculate solution path of CVR (along the sequence of lambda)
CVRPath <- function(Y, Xlist, event=NULL, rank = 2, eta = 0.5, Lamseq, 
                     family = c("gaussian", "binomial", "poisson", "cox"),
                     Wini = NULL, penalty = c("GL1", "L1"), opts, warm = T){
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  Y <- as.matrix(Y)
  if(family=='cox') event=as.matrix(event)
  K <- length(Xlist)
  
  if (dim(Lamseq)[2] != K) 
    stop("The length of Lambda is not the same to K!");  
  nlam <- dim(Lamseq)[1]
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (is.null(Wini)) 
    Wini <- SparseCCA(Xlist[[1]], Xlist[[2]], rank, 0.7, 0.7)
  
  p <- c()
  WPath <- list()
  for (k in 1:K) {
    p[k] <- dim(Xlist[[k]])[2]
    WPath[[k]] <- array(0, c(nlam, p[k], rank));
  }
  betaPath <-  matrix(0, nlam, rank)
  alphaPath <- rep(0, nlam) 
  iterPath <- rep(0, nlam)
  fitwork <- list()
  # warm <- TRUE
  Wwork <- Wini
  
  for (i in 1:nlam){
    if (i == 1) {
      optswork <- opts;
      optswork$maxIters <- 2 * opts$maxIters;
    } else {
      optswork <- opts;
      if (warm) 
        Wwork <- fitwork$W; ## warm start
    }   
    
    fitwork <- cvrsolver(Y, Xlist, event=event, rank, eta, as.vector(Lamseq[i, ]), family, Wwork, penalty, optswork);
    Wzero <- 0
    for (k in 1:K) 
      Wzero <- Wzero + sum(fitwork$W[[k]]!=0);
    
    if  (Wzero == 0) 
      break;
    
    for (k in 1:K) 
      WPath[[k]][i, , ] <- fitwork$W[[k]];
    
    betaPath[i, ] <- fitwork$beta;
    alphaPath[i] <- fitwork$alpha;
    iterPath[i] <- fitwork$iter
  }
  return(list(WPath = WPath, betaPath = betaPath, alphaPath = alphaPath, iterPath = iterPath))
}

## cross validate to select lambda, with fixed rank and eta
TuneCVR <- function(Y, Xlist, event=NULL, rank, eta, Lamseq = NULL, 
                   family = c("gaussian", "binomial", "poisson", "cox"), 
                   Wini = NULL, penalty = c("GL1", "L1"), 
                   nfold = 10, foldid = NULL, type.measure = NULL, opts, warm=T){    
  ##  cross-validate to select lam's
  ##  input: Lamseq: nlam by K, K=2
  ##  output: list(Lamseq = Lamseq, cverror = pred, cvm = mpred, msparse = msparse, 
  ##           Lamhat = Lamhat, W1trace = W1trace, W2trace = W2trace, cvr.fit = fit,  
  ##           alpha = alpha.refit, beta = beta.refit, type.measure = type.measure)
  ##          cvr.fit:  Refit all the data using the selected parameters, 
  ##                    see the output of cvrsolver(). Also output the coefficients of regressing                      
  ##                    the response to the combined fitted canonical variates (Y ~ cbind(X1W1, X2W2)). 
  ##  
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  Y <- as.matrix(Y)
  if(family=='cox') event=as.matrix(event)
  n <- dim(Y)[1]
  K <- length(Xlist)
  X1 <- as.matrix(Xlist[[1]])
  X2 <- as.matrix(Xlist[[2]])
  p1 <- ncol(X1) 
  p2 <- ncol(X2)

  if (is.null(Wini)) 
    Wini <- SparseCCA(X1, X2, rank, 0.7, 0.7)
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
  opts$nrank <- rank
  opts$family <- family
  opts$penalty <- penalty
  
  Lamseq <- as.matrix(Lamseq)
  if (dim(Lamseq)[2] != K) {
    lamseq <- 10^(seq(-2, 1.3, len = 50))
    Lamseq <- cbind(lamseq, lamseq * p2 / p1)
  }
  nlam <- dim(Lamseq)[1];
  
  pred <- matrix(0, nlam, nfold);
  itertrace <- matrix(0, nlam, nfold)
  W1trace <- array(0, c(nlam, nfold, p1, rank))
  W2trace <- array(0, c(nlam, nfold, p2, rank))
  
  sparse <- matrix(0, nlam, nfold)    # monitor sparsity
  if (is.null(foldid)) {  
    tmp <- ceiling(c(1:n) / (n / nfold));
    foldid <- sample(tmp, n)   
  } 
  
  for (i in 1:nfold) {
    X1train <- Xlist[[1]][foldid != i, ];
    X2train <- Xlist[[2]][foldid != i, ];
    Xtrain  <- list(X1 = X1train, X2 = X2train);
    Ytrain  <- as.matrix(Y[foldid != i,  ]);
    X1test  <- Xlist[[1]][foldid == i, ];
    X2test  <- Xlist[[2]][foldid == i, ];
    Ytest   <- as.matrix(Y[foldid == i, ]);
    eventtrain = eventtest = matrix(NA)
    if(family=='cox'){
      eventtrain = event[foldid != i,  ]
      eventtest = event[foldid == i,  ]
    }
    if (is.null(opts$W)) 
      opts$W = SparseCCA(X1train, X2train, rank, 0.7, 0.7)
         
    obj <- CVRPath(Ytrain, Xtrain, eventtrain, rank, eta, Lamseq, family, opts$W, penalty, opts);  
    W1trace[, i, , ] <- obj$WPath[[1]]
    W2trace[, i, , ] <- obj$WPath[[2]]  
    itertrace[, i] <- obj$iterPath  
    for (j in 1:nlam){
      W1 <- obj$WPath[[1]][j, , ]
      W2 <- obj$WPath[[2]][j, , ]
      ## sparse = proportion of nonzero elements among W1 and W2
      sparse[j, i] <- (sum(W1 != 0) + sum(W2 != 0)) / (length(W1) + length(W2))  
      alpha <- as.numeric(obj$alphaPath[j])
      beta <- obj$betaPath[j, ]
      pred[j, i] <- PredCVR(Ytest, X1test, X2test, eventtest, W1, W2, alpha, beta,           
                          Ytrain, X1train, X2train, eventtrain, family = family, 
                          refit = TRUE, type.measure = type.measure)    
    }
  }
  
  mpred <- apply(pred, 1, mean)
  msparse <- apply(sparse, 1, mean)
  spthresh <- opts$spthresh
  if (type.measure == "mse" | type.measure == "deviance") {
    ## search the best lam only among sparse models
    ilam <- which.min(mpred[msparse <= spthresh]) + min(which(msparse <= spthresh)) - 1
    ifold <- which.min(pred[ilam, ])
    mmpred <- min(mpred[msparse <= spthresh])
  } else if (type.measure %in% c("auc", "UnoC", "C")) {
    ilam <- which.max(mpred[msparse <= spthresh]) + min(which(msparse <= spthresh)) - 1
    ifold <- which.max(pred[ilam, ])
    mmpred <- max(mpred[msparse <= spthresh])
  }  
  Lamhat <- Lamseq[ilam, ]    
  ## refit all the data with Lamhat   
  if (warm) {
    Wini <- list(W1 = as.matrix(W1trace[max(ilam - 1, 1), ifold, , ]), 
                 W2 = as.matrix(W2trace[max(ilam - 1, 1), ifold, , ]))
  }

  fit <- cvrsolver(Y, Xlist, event, rank, eta, Lamhat, family, Wini, penalty, opts)

  XW <- cbind(X1 %*% fit$W[[1]], X2 %*% fit$W[[2]])  
  if (family == "gaussian") {   
    ab <- coef(lm(Y ~ XW)); 
    ab[is.na(ab)] <- 0
    alpha.refit <- ab[1]
    beta.refit <- ab[-1]
  } else if (family == "cox") {   
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
    
  refit <- list(alpha = alpha.refit, beta = beta.refit, W = fit$W)
  cv.out <- list(Lamseq = Lamseq, eta = eta, rank = rank, 
                 cverror = pred, cvm = mpred, mmpred = mmpred, msparse = msparse, 
                 Lamhat = Lamhat, W1trace = W1trace, W2trace = W2trace, cvr.fit = fit,  
                 refit = refit, type.measure = type.measure, spthresh = spthresh)
  class(cv.out) <- "TuneCVR"
  return(cv.out) 
}

## plot mean cv error and sparsity of the solution path from TuneCVR objects
plot.TuneCVR <- function(x) {
  par(mfrow = c(2, 1))
  plot(log(x$Lamseq[, 1]), x$cvm, xlab="",
       ylab = x$type.measure, type = "l", col = "red",
       main = paste("CVR(eta = ", round(x$eta, 5), ", rank = ", x$rank, ")", sep=""))
  abline(v = log(x$Lamhat[1]))
  plot(log(x$Lamseq[, 1]), x$msparse, 
       xlab = expression(log ~ lambda),  
       ylab = "sparsity", type = "l")
  abline(h = x$spthresh)
}




#' @title Fit canonical variate regression with tuning parameters selected by cross validation.
#' 
#' @description This function fits the solution path of canonical variate regression, 
#'      with tuning parameters selected by  cross validation. The tuning parameters 
#'      include the rank, the  \eqn{\eta} and the \eqn{\lambda}.
#' @usage 
#'  CVR(Y, Xlist, event=NULL, rankseq = 2, neta = 10, etaseq = NULL, nlam = 50, 
#'      Lamseq = NULL, family = c("gaussian", "binomial", "poisson", "cox"),  
#'      Wini = NULL, penalty = c("GL1", "L1"), nfold = 10, foldid = NULL,  
#'      opts = list(), type.measure = NULL)
#'   
#' @param Y A univariate response variable.
#' @param Xlist A list of two covariate matrices as in \code{cvrsolver}.
#' @param event event indicator for survival outcome (i.e. \code{family}=\code{'cox'}). 
#' @param rankseq A sequence of candidate ranks. The default is a single value 2.
#' @param neta Number of \eqn{\eta} values. The default is 10.
#' @param etaseq A sequence of  length \code{neta} containing candidate \eqn{\eta} values between 0 and 1. 
#'          The default is 10^seq(-2, log10(0.9), length = neta).
#' @param nlam Number of  \eqn{\lambda} values. The default is 50.
#' @param Lamseq A matrix of   \eqn{\lambda}  values. The column number is the number of sets in \code{Xlist},
#'                and the row number is \code{nlam}. The default is 10^(seq(-2, 2, length = nlam)) for each column.
#' @param family Type of response as in \code{cvrsolver}. The default is \code{"gaussian"}.
#' @param Wini  A list of initial loading W's. The default is from the SparseCCA solution. See \code{SparseCCA}.
#' @param penalty Type of penalty on loading matrices W's as in \code{cvrsolver}. The default is \code{"GL1"}.              
#' @param nfold Number of folds in cross validation. The default is 10.
#' @param foldid Specifying training and testing sets in cross validation; random generated if not supplied. 
#'          It remains the same across different rank and \eqn{\eta}.
#' @param opts A list of options for controlling the algorithm. The default of \code{opts$spthresh} is 0.4, which means 
#'              we only search sparse models with at most 40\% nonzero entries in W1 and W2. See the other options 
#'              (\code{standardization}, \code{maxIters} and \code{tol}) in \code{cvrsolver}.                    
#' @param type.measure  Type of measurement used in cross validation. \code{"mse"} for Gaussian, \code{"auc"} for binomial,
#'                      and \code{"deviance"} for binomial and Poisson.
#'                        
#' @details  In this function, the rank, \eqn{\eta} and \eqn{\lambda} are tuned by  cross validation. CVR then is refitted with
#'            all  data using  the selected tuning parameters.  The \code{plot} function  shows the tuning of  \eqn{\lambda}, 
#'            with selected rank and \eqn{\eta}.
#'            
#' @return An object with S3 class "CVR" containing the following components
#' @return \item{cverror}{A matrix containing the CV errors. The number of rows is the length 
#'                    of \code{etaseq} and the number of columns is the length of \code{rankseq}.}  
#' @return \item{etahat}{Selected \eqn{\eta}.}
#' @return \item{rankhat}{Selected rank.}
#' @return \item{Lamhat}{Selected  \eqn{\lambda}'s.}
#' @return \item{Alphapath}{An array containing  the fitted paths of the intercept term \eqn{\alpha}.}
#' @return \item{Betapath}{An array containing the fitted paths of the regression coefficient \eqn{\beta}.}
#' @return \item{W1path, W2path}{Arrays containing the fitted paths of  W1 and W2.}
#' @return \item{foldid}{\code{foldid} used in cross validation.} 
#' @return \item{cvout}{Cross validation results using selected \eqn{\eta} and rank.}
#' @return \item{solution}{A list including the solutions of \eqn{\alpha}, \eqn{\beta}, W1 and W2, by refitting all the data 
#'                        using selected tuning parameters.}
#'                        
#' @author Chongliang Luo, Kun Chen.
#' 
#' @references Chongliang Luo, Jin Liu, Dipak D. Dey and Kun Chen (2016) Canonical variate regression. 
#'                Biostatistics, doi: 10.1093/biostatistics/kxw001.
#' @seealso \code{\link{cvrsolver}}, \code{\link{SparseCCA}}, \code{\link{SimulateCVR}}.
#' @examples 
#' ############## Gaussian response ###################### 
#' set.seed(42)   
#' mydata <- SimulateCVR(family = "g", n = 100, rank = 4, p1 = 50, p2 = 70,  
#'                   pnz = 10, beta = c(2, 1, 0, 0))
#' X1 <- mydata$X1;
#' X2 <- mydata$X2
#' Xlist <- list(X1 = X1, X2 = X2); 
#' Y <- mydata$y
#' ## fix rank = 4, tune eta and lambda   
#' ##out_cvr <- CVR(Y, Xlist, rankseq = 4, neta = 5, nlam = 25,  
#' ##               family = "g", nfold = 5)
#' ## out_cvr$solution$W[[1]];  
#' ## out_cvr$solution$W[[2]];     
#' ### uncomment to see plots 
#' ## plot.CVR(out_cvr)
#' ## 
#' ## Distance of subspaces 
#' ##U <- mydata$U
#' ##Pj <- function(U) U %*% solve(t(U) %*% U, t(U)) 
#' ##sum((Pj(U) - (Pj(X1 %*% out_cvr$sol$W[[1]]) + Pj(X2 %*% out_cvr$sol$W[[2]]))/2)^2) 
#' ## Precision/Recall rate
#' ## the first 10 rows of the true W1 and W2 are set to be nonzero
#' ##W12 <- rbind(out_cvr$sol$W[[1]], out_cvr$sol$W[[1]])
#' ##W12norm <- apply(W12, 1, function(a)sqrt(sum(a^2)))
#' ##prec <- sum(W12norm[c(1:10, 51:60)] != 0)/sum(W12norm != 0); prec 
#' ##rec <- sum(W12norm[c(1:10, 51:60)] != 0)/20; rec
#' ## sequential SparseCCA, compare the Distance of subspaces and Prec/Rec
#' ##W12s <- SparseCCA(X1, X2, 4)
#' ## Distance larger than CVR's 
#' ##sum((Pj(U) - (Pj(X1 %*% W12s$W1) + Pj(X2 %*% W12s$W2))/2)^2) 
#' ##W12snorm <- apply(rbind(W12s$W1, W12s$W2), 1, function(a)sqrt(sum(a^2)))
#' ## compare Prec/Rec 
#' ##sum(W12snorm[c(1:10, 51:60)] != 0)/sum(W12snorm != 0);  
#' ##sum(W12snorm[c(1:10, 51:60)] != 0)/20; 
#' 
#' ############## binary response ########################
#' set.seed(12) 
#' mydata <- SimulateCVR(family = "binomial", n = 300, rank = 4, p1 = 50,  
#'                        p2 = 70, pnz = 10, beta = c(2, 1, 0, 0))
#' X1 <- mydata$X1; X2 <- mydata$X2
#' Xlist <- list(X1 = X1, X2 = X2); 
#' Y <- mydata$y
#' ## out_cvr <- CVR(Y, Xlist, 4, neta = 5, nlam=25, family = "b", nfold = 5)  
#' ## out_cvr$sol$W[[1]];  
#' ## out_cvr$sol$W[[2]];    
#' ## plot.CVR(out_cvr)
#' 
#' ############## Poisson response ######################
#' set.seed(34)
#' mydata <- SimulateCVR(family = "p", n = 100, rank = 4, p1 = 50,    
#'                        p2 = 70, pnz = 10, beta = c(0.2, 0.1, 0, 0))
#' X1 <- mydata$X1; X2 <- mydata$X2
#' Xlist <- list(X1 = X1, X2 = X2); 
#' Y <- mydata$y
#' ## etaseq <- 10^seq(-3, log10(0.95), len = 10)
#' ## out_cvr <- CVR(Y, Xlist, 4, neta = 5, nlam = 25, family = "p", nfold = 5)
#' ## out_cvr$sol$W[[1]];  
#' ## out_cvr$sol$W[[2]];    
#' ## plot.CVR(out_cvr)  
#' @export 
CVR <- function(Y, Xlist, event=NULL, rankseq = 2, neta = 10, etaseq = NULL, nlam = 50,  
                Lamseq = NULL, family = c("gaussian", "binomial", "poisson", "cox"), 
                Wini = NULL, penalty = c("GL1", "L1"), nfold = 10, foldid = NULL, 
                opts = list(), type.measure = NULL, warm=T){
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if(is.null(family)) 
    stop('Must specify family of the response!')  
  if(is.null(penalty)) {
    penalty = "GL1"; 
    cat('Use default: penalty = "GL1".\n')
  }
  
  Y <- as.matrix(Y)
  if(family=='cox'){
    if(is.null(event)) stop("need event indicator for family cox")
  }else{
    event=rep(0,n)
  }
  event <- as.matrix(event)
  X1 <- as.matrix(Xlist[[1]])
  X2 <- as.matrix(Xlist[[2]])
  Xlist <- list(X1 = X1, X2 = X2)
  n <- dim(Y)[1]
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  if(p1 == 1 | p2 == 1) 
    stop("X1 and X2 should have at least 2 columns!")
  if(nrow(X1) != nrow(X2)) 
    stop("X1 and X2 should have the same row number!")
  if(n != nrow(X1)) 
    stop("Y and X1, X2 should have the same row number!")
  
  opts$n  <- n
  opts$p1 <- p1
  opts$p2 <- p2
  
  mrank <- max(rankseq)
  if (is.null(Wini)) {
    Wini <- SparseCCA(X1, X2, mrank, 0.7, 0.7); 
    cat('Use default: initial W = SparseCCA(X1, X2, rank, 0.7, 0.7).\n')
  }
  if (is.null(opts$standardization)) {
    opts$standardization <- TRUE; 
    #cat('Use default: standardization = TRUE.\n')
  }  
  if (is.null(opts$maxIters)) {
    opts$maxIters <- 300; 
    #cat('Use default: maxIters = 300.\n')
  }
  if (is.null(opts$tol)) {
    opts$tol <- 0.01; 
    #cat('Use default: tol = 0.01.\n')
  }
  if (is.null(opts$spthresh)) {
    spthresh <- 0.4
    opts$spthresh = spthresh; 
    #cat('Use default: spthresh = 0.4.\n')
  }
   
  if (is.null(etaseq)) {
    if(is.null(neta)) 
      neta <- 10
    etaseq <- 10^seq(-2, log10(0.9), len = neta)
  }
  if (is.null(Lamseq)) {
    if (is.null(nlam)) 
      nlam <- 50
    lamseq <- 10^(seq(-2, 2, len = nlam))
    Lamseq <- cbind(lamseq, lamseq)
  } 
  
  if (is.null(foldid)) { 
    if(is.null(nfold))
      nfold <- 10
    tmp <- ceiling(c(1:n)/(n/nfold));
    foldid <- sample(tmp, n) 
  } else {
    if(nfold != length(unique(foldid))){
      #warning("nfold not equal to no. unique values in foldid")
      nfold <- length(unique(foldid)) 
    }
  }
  
  if(is.null(type.measure)){ 
    if (family == "gaussian") {
      type.measure <- "mse";
    } else if (family == "binomial") {
      if (n / nfold > 10 & is.null(type.measure)) 
        type.measure <- "auc";
      if (n / nfold < 11 & is.null(type.measure)) 
        type.measure <- "deviance";   
    } else if (family == "poisson") {
      type.measure = "deviance";   
    } else if (family == "cox") {
      type.measure = "UnoC";   
    }
  }
  
  Lr <- length(rankseq)
  Leta <- length(etaseq)

  pred <- matrix(0, Leta, Lr)
  Lamhat <- matrix(0, Leta, Lr)
  Alphapath  <- matrix(0, Leta, Lr)
  Betapath <- array(0, c(Leta, Lr, mrank))
  W1path <- array(0, c(Leta, Lr, p1, mrank))
  W2path <- array(0, c(Leta, Lr, p2, mrank))
  
  
  for (ir in 1:Lr) {
    for (ieta in 1:Leta) {
      Wini0 = list(W1=as.matrix(Wini$W1[,1:rankseq[ir]]), W2=as.matrix(Wini$W2[,1:rankseq[ir]]))
      obj_cv <- TuneCVR(Y, Xlist, event, rankseq[ir], etaseq[ieta], Lamseq, family, Wini0, penalty, 
                      nfold = nfold, foldid = foldid, type.measure = type.measure, opts,  warm=warm)      
      pred[ieta, ir] <- obj_cv$mmpred  
      Lamhat[ieta, ir] <- obj_cv$Lamhat[1]
      Alphapath[ieta, ir]  <- obj_cv$cvr.fit$alpha
      Betapath[ieta, ir, 1:rankseq[ir]] <- obj_cv$cvr.fit$beta
      W1path[ieta, ir, , 1:rankseq[ir]] <- obj_cv$cvr.fit$W[[1]]
      W2path[ieta, ir, , 1:rankseq[ir]] <- obj_cv$cvr.fit$W[[2]]
    }
  }

  if (type.measure == "mse" | type.measure == "deviance") {
    ind <- apply(pred, 2, which.min)
    ind2 <- which.min(apply(pred, 2, min))
    mmpred <- min(pred)
  } else if (type.measure %in% c("auc", "UnoC", "C")) {
    ind <- apply(pred, 2, which.max)
    ind2 <- which.max(apply(pred, 2, max))
    mmpred <- max(pred)
  }
  
  etahat <- etaseq[ind[ind2]]
  rankhat <- rankseq[ind2]
  Wini0 = list(W1=as.matrix(Wini$W1[,1:rankhat]), W2=as.matrix(Wini$W2[,1:rankhat]))
  cvout <- TuneCVR(Y, Xlist, event, rankhat, etahat, Lamseq, family, Wini0,  
                  penalty, nfold, foldid = foldid, type.measure, opts, warm=warm)
  
  cvrout <- list(cverror = pred, etahat = etahat, rankhat = rankhat, Lamhat = Lamhat,  
                 Alphapath = Alphapath, Betapath = Betapath, W1path = W1path, 
                 W2path = W2path, foldid = foldid, cvout = cvout,  
                 solution = append(cvout$refit, list(rankseq = rankseq, etaseq = etaseq, Lamseq = Lamseq)))
  class(cvrout) <-  "CVR" 
  return(cvrout)
}


#' @title Plot a CVR object.
#' 
#' @description Plot the tuning of CVR
#' 
#' @usage 
#'  \method{plot}{CVR}(x, ...)           
#' 
#' @param x A CVR object.
#' @param ... Other graphical parameters used in plot.
#' 
#' @details  The first plot is mean cv error vs log(\eqn{\lambda}). The type of mean cv error is
#'           decided by \code{type.measure} (see parameters of \code{CVR}). The selected  \eqn{\lambda}  is marked 
#'           by a vertical line in the plot. The second plot is sparsity vs log(\eqn{\lambda}). 
#'           Sparsity is the proportion of non-zero elements in fitted W1 and W2. 
#'           The threshold is marked by a horizontal line. 
#'           Press ENTER to see the second plot, which shows the tuning of \eqn{\eta}.
#' @export 
plot.CVR <- function(x, ...) {
  plot.TuneCVR(x$cvout)
  cat ("Press [enter] to continue")
  line <- readline()
  par(mfrow = c(1, 1))
  plot(x$cverror[, which(x$solution$rankseq == x$rankhat)] ~ log10(x$solution$etaseq), 
       xlab = "log10(etaseq)", ylab = x$cvout$type.measure, 
       main = "CVR(tune eta)")
}



