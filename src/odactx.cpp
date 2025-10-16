#include <RcppArmadillo.h> 
//#include <Rcpp.h>
#include <R.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;
using namespace std;
using namespace arma;
//[[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
Rcpp::List rcpp_coxph_logL_xt(arma::mat mydata, vec beta, bool derivative){
  // names(mydata)[1:4] = c('id', 'time1', 'time2', 'event')
  vec id=mydata.col(0), time1=mydata.col(1), time2=mydata.col(2), event=mydata.col(3);
  mat dd = mydata.rows(find(event==1));
  vec time2d = dd.col(2), eXb;
  int nd = dd.n_rows, p = dd.n_cols-4, i, j;
  mat Grad(nd, p, fill::zeros), Hess(nd, pow(p,2), fill::zeros), SumeXbX(1,p);
  vec ddSumeXb(nd);
  double SumeXb;
  
  for(i=0; i<nd; i++){
    // mat X = mydata.submat(find(time1<time2d(i) && time2>=time2d(i)), span(4,p-1));  // risk set
    mat X0 = mydata.rows(find(time1<time2d(i) && time2>=time2d(i)));  
    mat X = X0.cols(4,p+4-1);
    eXb = exp(X*beta);
    SumeXb = sum(eXb); 
    ddSumeXb(i) = SumeXb;
    if(derivative==true){ 
      SumeXbX = eXb.t()*X;  
      vec mm(pow(p,2), fill::zeros);
      for (j=0; j<X.n_rows; j++) { 
        mm = mm + eXb(j) * vectorise(X.row(j).t() * X.row(j));          // XXeXb
      }
      Grad.row(i) = SumeXbX/SumeXb;  
      Hess.row(i) = -(SumeXb*mm.t() - vectorise(SumeXbX.t() * SumeXbX).t())/pow(SumeXb,2);
    }
  }
  
  mat X = dd.cols(4,p+4-1);
  vec Xb = X * beta, grad(p), hess(pow(p,2)); 
  if(derivative==true){
    grad = sum(X - Grad, 0).t();    // mean
    hess = sum(Hess, 0).t();        // mean // reshape(mean(Hess, 0), p, p);
  }  
  // cout << i << ", ckpt=9" << endl;
  Rcpp::List out;
  out["logL"] = -sum(Xb - log(ddSumeXb));  // mean
  out["grad"] = -grad;
  out["hess"] = -hess;
  
  return(out);
}