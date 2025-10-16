#include <RcppArmadillo.h>
#include <R.h>
#include <stdio.h>
#include <math.h>

using namespace Rcpp;
using namespace std;
using namespace arma;
using namespace sugar;
 
// row-wise matrix norm, faster than R version
// [[Rcpp::export]]
arma::sp_mat norm21sp_rcpp (arma::sp_mat A){
 arma::sp_mat A21 = sqrt(sum(square(A),1));
 return(A21);
}

// [[Rcpp::export]]
arma::mat norm21_rcpp (arma::mat A){
 arma::mat A21 = sqrt(sum(square(A),1));
 return(A21);
}
 

// [[Rcpp::export]]
arma::sp_mat MGlasso_ss0_Rcpp (arma::sp_mat XR, arma::sp_mat XX, double lam, 
                              arma::sp_mat A0, double conv, int maxiter) {
  // min |Y-XA|^2 + lam*|A| but use summary stat XY and XX 
  // subgradient method, vector version of the soft-thresholding rule: Chen&Huang2012, eq 13
  // XR = XY - XX*A1;  // pxq
  int  p=XX.n_cols, iter=0, j;   // q=XY.n_cols, i, 
  double  diff=10*conv, l2A1;
  arma::sp_mat XRj; 
  arma::sp_mat A1=A0; 
  arma::sp_mat t1;  // rowvec
  // if (lam.size() == 1) {lam = as_scalar(lam)*ones(p);}
  // cout << "lam=" << lam  << endl;
  while ((diff > conv) & (iter < maxiter)) {
    // cout << "iter=" << iter  << endl; 
    A0 = A1;
    for (j = 0; j < p; j++) { 
      XRj = XR.row(j) + XX(j,j)*A1.row(j);  // 1xq
      // cout << "XRj=" << XRj  << endl; 
      t1 = XRj / XX(j,j) * std::max(0.0, 1-lam/pow(accu(square(XRj)),0.5)/2);  
      // cout << "t1=" << t1  << endl; 
      // Ajdiff = t1 - A1.row(j);
      A1.row(j) = t1;
      XR = XR - XX.col(j)*(t1 - A0.row(j));
    }
    l2A1 = accu(square(A1));
    if (l2A1 <= 0.001) {
      iter = maxiter;
    } else {
      diff = pow(accu(square(A0 - A1))/l2A1,0.5);
      iter = iter + 1;
    }
  } 
  return(A1);
}
 
// this version returns dense mat
// [[Rcpp::export]]
arma::mat MGlasso_ss1_Rcpp (arma::mat XR, arma::mat XX, arma::vec lam, 
                             arma::mat A0, double conv, int maxiter) {
 // min |Y-XA|^2 + lam*|A| but use summary stat XY and XX 
 // subgradient method, vector version of the soft-thresholding rule: Chen&Huang2012, eq 13
 // XR = XY - XX*A1;  // pxq
 // lam is a vector with possible adaptive weights
 int  p=XX.n_cols, iter=0, j;   // q=XY.n_cols, i, 
 double  diff=10*conv, l2A1;
 arma::mat XRj; 
 arma::mat A1=A0; 
 arma::mat t1;  // rowvec
 // if (lam.size() == 1) {lam = as_scalar(lam)*ones(p);}
 // cout << "lam=" << lam  << endl;
 while ((diff > conv) & (iter < maxiter)) {
   // cout << "iter=" << iter  << endl; 
   A0 = A1;
   for (j = 0; j < p; j++) { 
     XRj = XR.row(j) + XX(j,j)*A1.row(j);  // 1xq
     // cout << "XRj=" << XRj  << endl; 
     t1 = XRj / XX(j,j) * std::max(0.0, 1-lam(j)/pow(accu(square(XRj)),0.5)/2);  
     // cout << "t1=" << t1  << endl; 
     // Ajdiff = t1 - A1.row(j);
     A1.row(j) = t1;
     XR = XR - XX.col(j)*(t1 - A0.row(j));
   }
   l2A1 = accu(square(A1));
   if (l2A1 <= 0.001) {
     iter = maxiter;
   } else {
     diff = pow(accu(square(A0 - A1))/l2A1,0.5);
     iter = iter + 1;
   }
 } 
 return(A1);
}

 
// Rcpp don't take list of sparse matrix as input from R, 
// so hard coded as XX1 ... XX6 (assume q=6), change if you have different q
// somehow using sparse mat is not faster as arma claimed, but save RAM. 
// [[Rcpp::export]]
Rcpp::List srrr_ss0_Rcpp (arma::sp_mat XYlist,   // arma::cube XXlist0,
                          arma::sp_mat XX1,
                          arma::sp_mat XX2,
                          arma::sp_mat XX3,
                          arma::sp_mat XX4,
                          arma::sp_mat XX5, // arma::sp_mat XX6,   // the same numbe as q=ncol(XYlist)
                          arma::sp_mat XX,    // sum(XXk)
                           int n, double Yss, // String method, 
                           arma::sp_mat A0, arma::sp_mat V0,
                           int nrank, double lambda,
                           double conv, int maxiter,
                           double inner_conv, int inner_iter, // arma::vec WA,
                           int mglasso_impute,
                           int procrustes_impute) {
   // Y = bdiag(Y1, ..., Yq)
   // min |Y-XAV'|^2 + lamA*|A|, s.t. V'V = I_nrank 
   int p=XYlist.n_rows, q=XYlist.n_cols, k;     
   sp_mat A1=A0, C0=sp_mat(A0*V0.t()), C1=C0;  // 
   // sp_mat XX(p,p);  
   sp_mat XYd(p,q);   // XY, =zeros
   arma::field<sp_mat> XXlist(q);
   
   XXlist(0) = XX1;
   XXlist(1) = XX2;
   XXlist(2) = XX3;
   XXlist(3) = XX4;
   XXlist(4) = XX5;
   // XXlist(5) = XX6;   // change if you have different q
   
   // // convert to field of sparse mat
   // for(k=0; k<q; k++){
   //   XXlist(k) = XXlist0.slice(k);
   //   XX = XX + XXlist(k);
   // }
 
   // double xrank = accu(svd(XX) > 0.0001);
   int iter = 0, xrank = min(n, p);  
   bool conv_flag;
   double l2C1, sse, df, dfu0, dfv0, BIC, BICP, AIC, GCV, GIC;
   sp_mat XR, V1=V0, W;
   arma::mat u, v;  
   arma::vec s, diff(maxiter+1); 
   diff.fill(2*conv);   
   Rcpp::List out;   
   
   while ((iter < maxiter) & (diff(iter) > conv)) {
     V0 = V1;
     A0 = A1;
     C0 = C1;
     
     //// fix A, update V by procrustes svd 
     if(procrustes_impute==1){  // impute Y
       for(k=0; k<q; k++){ 
         XYd.col(k) = XYlist.col(k) - XXlist(k)*C1.col(k); // pxq
       }
       // XY = XX*C1 + XYd;     // X'Y with off-diagonal Y imputed as X*C
       W = trans(XX*C1 + XYd) * A1;      // q r
     }else{                     // ignore off-diagonal Y
       W = XYlist.t() * A1;     // qxr
     } 
     svds(u, s, v, W, nrank);      //qXq qXr rXr    
     u = u.cols(0, nrank-1);
     V1 = u * v.t();
     C1 = A1*V1.t();
     // cout << "V1=" << V1 << endl;  //
     
     //// fix V, update A by mglasso
     // subgradient method, vector version of the soft-thresholding rule: Chen&Huang2012, eq 13
     if(mglasso_impute==1){    // impute off-diagonal Y
       for(k=0; k<q; k++){ 
         XYd.col(k) = XYlist.col(k) - XXlist(k)*C1.col(k);
       }
       // XY = XX*C1 + XYd;       
       // XR = XY*V1 - XX*A0;
       XR = XYd*V1;             // pxr
     }else{                     // ignore off-diagonal Y
       XR = XYlist*V1;
       for(k=0; k<q; k++){ 
         XR = XR - XXlist(k)*A1*trans(V1.row(k))*V1.row(k);
       }
     }
     // int it =  out_MGlasso["iter"]; 
     A1 = MGlasso_ss0_Rcpp(XR, XX, lambda, A1, inner_conv, inner_iter); 
     // cout << "A1.n_nonzero=" << A1.n_nonzero << endl;
     // if(A1.n_nonzero<3 | accu(abs(A1))<0.01) break;  // all entry are 0
     C1 = A1*V1.t();
     l2C1 = accu(square(C1));
     
     if (l2C1 <= 0.001) {
       diff(iter) = 0;
     } else {
       iter = iter + 1;
       diff(iter) = pow(accu(square(C0 - C1))/l2C1,0.5);
     }
     // cout << "diff=" << diff(iter)  << endl;
   }
   
   diff = diff.subvec(0, iter); 
   sse = Yss;
   for(k=0; k<q; k++){
     sse = sse + as_scalar(arma::trans(C1.col(k))*XXlist(k)*C1.col(k)-2*arma::trans(XYlist.col(k))*C1.col(k)); 
   }
   dfu0 = A1.n_nonzero;     // accu(A1 != 0);
   dfv0 = V1.n_nonzero;     // accu(V1 != 0);
   df = dfu0 * xrank/p + dfv0 - nrank*nrank;   
   double logqn = std::log(static_cast <double>(q * n));
   double logsse = std::log(sse);
   BIC = logsse + df*logqn/(q*n);
   BICP = logsse + 2*df*logqn/(q*n);
   AIC = logsse + 2*df/(q*n);
   GIC = logsse + df*std::log(logqn)*std::log(static_cast <double>(p*q))/(q*n);
   double dfqn2 = pow((1 - df/(q*n)), 2);
   GCV = sse/(q*n*dfqn2);
   if (diff(iter) <= conv) {
     conv_flag = true;
   } else {
     conv_flag = false;
   }
   
   // out["XXlist"] = XXlist;
   // out["XX"] = XX;
   out["diff"] = diff;
   out["iter"] = iter;
   out["BIC"]  = BIC;
   out["BICP"] = BICP;
   out["AIC"]  = AIC;
   out["GCV"]  = GCV;
   out["GIC"]  = GIC;
   out["sse"]  = sse;
   out["df"]   = df;
   out["conv_flag"] = conv_flag;
   out["A"]  = A1;
   out["V"]  = V1;
   // out["C"]  = C1;
   return(out);
}



// This version takes XXlist dense array as input
// [[Rcpp::export]]
Rcpp::List srrr_ss1_Rcpp (arma::mat XYlist, 
                          arma::cube XXlist, 
                          arma::mat XX,       // sum_XXk 
                          arma::vec popu_wt,  // weight for each popu
                          arma::vec nn,   // size of each popu, to calc IC 
                          arma::vec Yss,  
                          arma::mat A0, arma::mat V0,
                          String method, double wgamma, arma::vec WA, // adaptive lasso penalty 
                          int nrank, double lambda,
                          double conv, int maxiter,
                          double inner_conv, int inner_iter, 
                          int mglasso_impute,
                          int procrustes_impute) {
  // Y = bdiag(Y1, ..., Yq)
  // XYlist, Yss include the weights
  // min Sum_k |Yk*w_k-XkAv_k|^2 + lamA*WA*|A|, s.t. V'V = I_nrank
  // w_k = sqrt(popu_wt_k), B=AV'diag{W_k^-1} 
  int p=XYlist.n_rows, q=XYlist.n_cols, k, n=accu(nn);     
  arma::mat XY, XYd = arma::zeros(p,q), C0 = A0*V0.t(), XR;  
  // XX = arma::zeros(p,p), 
  // for(k=0; k<q; k++){ 
  //   XX = XX + XXlist.slice(k);
  // }
  // double xrank = accu(svd(XX) > 0.0001);
  int xrank = min(n, p);
  int iter = 0;
  double dfu0, dfv0;
  bool conv_flag;
  double l2C1, df, sse, BIC, BICP, AIC, GCV, GIC;
  arma::vec SSE, s, diff(maxiter+1); 
  diff.fill(2*conv);
  arma::mat V1;
  arma::mat A1;
  arma::mat C1; 
  arma::mat W;
  arma::mat u;
  arma::mat v;
  arma::mat residual;
  Rcpp::List out;

  if(WA.is_empty()) {
    WA=ones(p); 
    if((method=="adglasso")) {  // adaptive lasso penalty 
      for(int i=0;i<p;i++){
       WA(i) = pow(sqrt(accu(square(A0.row(i))))+0.001,-wgamma);
      }
      WA = WA / accu(WA) * p; 
    }
  }
  
  V1 = V0;
  A1 = A0;
  C1 = C0;
  
  while ((iter < maxiter) & (diff(iter) > conv)) {
    V0 = V1;
    A0 = A1;
    C0 = C1;
    
    //// update V by procrustes svd, with imputed XY
    if(procrustes_impute==1){
      for(k=0; k<q; k++){
        // XYd.col(k) = XYlist(k) - XXlist(k)*C1.col(k);
        XYd.col(k) = XYlist.col(k) - XXlist.slice(k)*C1.col(k);
      }
      XY = XX*C1 + XYd;     // X'Y with off-diagonal Y imputed as X*C
      W = XY.t() * A1;      //q r
    }else{
      W = XYlist.t() * A1; 
    }
    svd(u, s, v, W);      //qXq qXr rXr
    u = u.cols(0, nrank-1);
    V1 = u * v.t();
    C1 = A0*V1.t();
    
    //// update A by mglasso
    if(mglasso_impute==1){    // impute off-diagonal Y
      for(k=0; k<q; k++){ 
        XYd.col(k) = XYlist.col(k) - XXlist.slice(k)*C1.col(k);
      } 
      XR = XYd*V1;             // pxr
    }else{                     // ignore off-diagonal Y
      XR = XYlist*V1;
      for(k=0; k<q; k++){ 
        XR = XR - XXlist.slice(k)*A1*trans(V1.row(k))*V1.row(k);
      }
    } 
    A1 = MGlasso_ss1_Rcpp(XR, XX, lambda*WA, A1, inner_conv, inner_iter); 
    C1 = A1*V1.t();
    
    l2C1 = accu(square(C1));
    if (l2C1 == 0) {
      diff(iter) = 0;
    } else {
      iter = iter + 1;
      diff(iter) = pow(accu(square(C0 - C1))/l2C1,0.5);
    }
  }
  
  diff = diff.subvec(0, iter); 
  // residual = Y - X * C1;
  // sse = accu(square(residual));
  // sse = Yss + trace(C1.t()*XX*C1-2*YX*C1);
  SSE = Yss; // SSE for each popu
  for(k=0; k<q; k++){
    // weighted SSE
    SSE(k) = SSE(k) + as_scalar(arma::trans(C1.col(k))*XXlist.slice(k)*C1.col(k)-2*arma::trans(XYlist.col(k))*C1.col(k) ); 
    C1.col(k) = C1.col(k) / sqrt(popu_wt(k));   // scale back C = AV'/sqrt(popu_wt)
    // SSE(k) = SSE(k)/popu_wt(k) + as_scalar(arma::trans(C1.col(k))*XXlist.slice(k)*C1.col(k)-2*arma::trans(XYlist.col(k))*C1.col(k)/sqrt(popu_wt(k))); 
  }
  sse = accu(SSE);
  dfu0 = accu(A1 != 0);
  dfv0 = accu(V1 != 0);
  df = dfu0 * xrank/p + dfv0 - nrank*nrank;   
  double logqn = std::log(static_cast <double>(q * n));
  double logsse = std::log(sse); 
  BIC = logsse + df*logqn/(q*n);  // n instead of q*n?
  BICP = logsse + 2*df*logqn/(q*n);
  AIC = logsse + 2*df/(q*n);
  GIC = logsse + df*std::log(logqn)*std::log(static_cast <double>(p*q))/(q*n);
  double dfqn2 = pow((1 - df/(q*n)), 2);
  GCV = sse/(q*n*dfqn2);
  // IC for each popu: d.f. = total d.f. / q
  arma::vec logSSE = log(SSE);
  arma::vec lognn = log(nn);
  arma::vec dfqnn2 = pow((1 - df/q/nn), 2);
  arma::vec BIC_k = logSSE + df/q*lognn/(nn);
  arma::vec BICP_k = logSSE + 2*df/q*lognn/(nn);
  arma::vec AIC_k = logSSE + 2*df/q/(nn);
  arma::vec GIC_k = logSSE + df/q*log(lognn)*std::log(p)/(nn);
  arma::vec GCV_k = SSE/(nn%dfqnn2);
  double loglogp = std::log(std::log(static_cast <double>(p)));
  arma::vec mBIC_k = logSSE + df/q*loglogp*lognn/(nn); // modified BIC, Wang, Li, and Leng (2009)
  
  if (diff(iter) <= conv) {
    conv_flag = true;
  } else {
    conv_flag = false;
  }
  
  out["diff"] = diff;
  out["iter"] = iter;
  out["WA_lasso"] = WA;
  out["BIC"]  = BIC;
  out["BICP"] = BICP;
  out["AIC"]  = AIC;
  out["GCV"]  = GCV;
  out["GIC"]  = GIC;
  out["sse"]  = sse;
  out["df"]   = df;
  
  out["BIC_k"]  = BIC_k;
  out["BICP_k"] = BICP_k;
  out["AIC_k"] = AIC_k;
  out["GIC_k"]  = GIC_k;
  out["GCV_k"] = GCV_k;
  out["mBIC_k"] = mBIC_k;
  out["SSE"]  = SSE;
  
  out["conv_flag"] = conv_flag;
  out["A"]  = A1;
  out["V"]  = V1;
  out["C"]  = C1; // scale back C = AV'/sqrt(popu_wt)
  return(out);
}





// Rcpp::List srrr_ss0_path_Rcpp (arma::mat XYlist,    //  matrix of p x q
//                           arma::cube XXlist,        //  array of p x p x q
//                           int n, double Yss,
//                           int nrank, int batch_size, 
//                           arma::vec lamA,         // int nlam, double lam_max, double lam_min,
//                           double conv, int maxiter, double inner_conv, int inner_iter,
//                           int mglasso_impute, int procrustes_impute) {
//     // String method,       // c("glasso", "adglasso"),
//     // String ic.type,      //  = c("BIC", "BICP", "AIC", "GCV", "GIC"),
//     int p = XYlist.n_rows, q = XYlist.n_cols, nlam=lamA.n_elem, l, k, check, tune_counter;
//     arma::mat XYd(p,q), V0(q, nrank, fill::zeros);
//     arma::mat A0;
//     // arma::sp_mat A0(p, nrank);
//     // arma::field<sp_mat> Apath(nlam), XXlist_active(q);  // , p, nrank, fill::zeros
//     arma::cube Apath(nlam, p, nrank, fill::zeros);
//     arma::cube Vpath(nlam, q, nrank, fill::zeros);
//     arma::mat ICpath(nlam, 5, fill::zeros), norm21_A0;  // , norm21_XYd
//     uvec kkt_idx;
//     arma::uvec minid, nz, nz_, active, active_;
//     IntegerVector tt = seq_len(p)-1, tt1, tt0;
//     Rcpp::List fit;
//     // ic.idx = which(c("BIC", "BICP", "AIC", "GCV", "GIC")==ic.type)
// 
//     // convert to field of sparse mat
//     // arma::sp_mat XX(p,p);
//     // arma::field<sp_mat> XXlist(q);
//     // for(k=0; k<q; k++){
//     //   XXlist(k) = XXlist0.slice(k);
//     //   XX = XX + XXlist(k);
//     // }
//     arma::mat XX = sum(XXlist, 2);  
// 
//     // sol path: lamA from large to small
//     for (l=0; l<nlam; l++) {      // start from null model
//       // A0 = Apath.row(std::max(l-1,0));
//       // V0 = Vpath.row(std::max(l-1,0));
//       norm21_A0 = norm21_rcpp(A0);
//       nz = find(norm21_A0!=0);   // non-zero idx
//       nz_ = find(norm21_A0==0);
//       if(nz.n_elem >0){
//         XYd.fill(0);
//         for(k=0; k<q; k++){
//           uvec uveck(1,fill::value(k));
//           XYd.submat(nz_,uveck) = XYlist.submat(nz_,uveck) - XXlist.slice(k).submat(nz_,nz)*A0.rows(nz)*V0.row(k); //
//         }
//       }else{
//         XYd = XYlist;
//       }
//       active = unique(join_cols(nz, sort_index(norm21_rcpp(XYd), "descend").subvec(0,batch_size-1)));
//       // active_ = span(p).shed_rows(active);
//       tt1 = as<IntegerVector>(wrap(active));
//       tt0 = sugar::sort(setdiff(tt, tt1));
//       active_ = as<uvec>(wrap(tt0))
//       // cat("\n l =", l, ", active=", length(active), ' ...')
//       check = 0;
//       while(check==0){
//         arma::cube XXlist_active(active.n_elem, active.n_elem, q);
//         for(k=1; k<q; k++) XXlist_active.slice(k)=XXlist.slice(k).submat(active,active);
//         // XXlist.for_each( [&](arma::sp_mat& X) { X.submat(active,active); } ); // lambda function
//         fit = srrr_ss1_Rcpp_inner(XYlist.rows(active), XXlist_active, n, Yss,
//                             A0.rows(active), V0, nrank, lamA(l),
//                             conv, maxiter, inner_conv, inner_iter,
//                             mglasso_impute, procrustes_impute);
//         A0.rows(active) = fit["A"];
//         V0 = fit["V"];
//         // KTT check
//         XYd.fill(zeros);
//         for(k=0; k<q; k++){
//           XYd.submat(active_,k) = XYlist.submat(active_,k) - XXlist.slice(k).submat(active_,active)*A0.rows(active)*V0.row(k);
//         }
//         // norm21_XYd = norm21_rcpp(XYd);
//         kkt_idx = find(norm21_rcpp(XYd) > lamA(l));
//         if(any(kkt_idx)){  // any(norm21_XYd > lamA(l))
//           active = join_cols(active, kkt_idx);
//           // active = join_cols(active, sort_index(norm21_XYd, 'descend').subvec(0,batch_size-1));
//           // cat('add=', length(kkt_idx), '...')
//         } else check=1;
//       }
// 
//       Apath(l) = A0;
//       Vpath.row(l) = V0;
//       ICpath(l) = fit["BIC"];
//       // ICpath(l, 1) = fit$BICP
//       // ICpath(l, 1) = fit$AIC
//       // ICpath(l, 1) = fit$GCV
//       // ICpath(l, 1) = fit$GIC
// 
//       // early stop: if IC increases in previous 5 lam's, stop tuning
//       if(l>=1 & ICpath(l) - ICpath(l-1) > 0) tune_counter = tune_counter+1;
//       else tune_counter = 0;
//       if(tune_counter==5) break;
//     }
// 
//     // colnames(ICpath) = c('BIC', 'BICP', 'AIC', 'GCV', 'GIC')
//     minid = sort_index(ICpath)(0); // which_min(ICpath);
//     A = Apath(minid);
//     V = Vpath.row(minid);
//     // coefMat = A * V.t();
// 
//     // out['XYlist'] = XYlist;
//     // out['XXlist'] = XXlist;
//     out["rank"] = nrank;
//     out["Apath"] = Apath;
//     out["Vpath"] = Vpath;
//     out["icpath"] = ICpath;
//     out["lambda"] = lamA;
//     out["minid"] = minid;
//     // out["coef"] = coefMat;
//     out["A"] = A;
//     out["V"] = V;
//     return(out);
// }


// this is a inner func to be called in srrr_sum_path
// XXlist: input within Rcpp, a field of sparse mat
// XX: input within Rcpp,  sparse mat
// // [[Rcpp::export]]
// Rcpp::List srrr_ss0_Rcpp_inner (arma::mat XYlist,
//                           arma::field<sp_mat> XXlist,
//                           arma::sp_mat XX,
//                           int n, double Yss, // String method,
//                           arma::sp_mat A0, arma::mat V0,
//                           int nrank, double lambda,
//                           double conv, int maxiter,
//                           double inner_conv, int inner_iter, // arma::vec WA,
//                           int mglasso_impute,
//                           int procrustes_impute) {
//   // Y = bdiag(Y1, ..., Yq)
//   // min |Y-XAV'|^2 + lamA*|A|, s.t. V'V = I_nrank
//   int p=XYlist.n_rows, q=XYlist.n_cols, k;
//   sp_mat A1=A0, C0=sp_mat(A0*V0.t()), C1=C0;  // XX(p,p),
//   mat XYd=zeros(p,q);   // XY,
//   // arma::field<sp_mat> XXlist(q);
//   // // convert to field of sparse mat
//   // for(k=0; k<q; k++){
//   //   XXlist(k) = XXlist0.slice(k);
//   //   XX = XX + XXlist(k);
//   // }
// 
//   // double xrank = accu(svd(XX) > 0.0001);
//   int iter = 0, xrank = min(n, p);
//   bool conv_flag;
//   double l2C1, sse, df, dfu0, dfv0, BIC, BICP, AIC, GCV, GIC;
//   arma::mat XR, V1=V0, W, u, v;
//   arma::vec s, diff(maxiter+1);
//   diff.fill(2*conv);
//   Rcpp::List out;
// 
//   while ((iter < maxiter) & (diff(iter) > conv)) {
//     V0 = V1;
//     A0 = A1;
//     C0 = C1;
// 
//     //// fix A, update V by procrustes svd
//     if(procrustes_impute==1){  // impute Y
//       for(k=0; k<q; k++){
//         XYd.col(k) = XYlist.col(k) - XXlist(k)*C1.col(k); // pxq
//       }
//       // XY = XX*C1 + XYd;     // X'Y with off-diagonal Y imputed as X*C
//       W = trans(XX*C1 + XYd) * A1;      // q r
//     }else{                     // ignore off-diagonal Y
//       W = XYlist.t() * A1;     // qxr
//     }
//     svd(u, s, v, W);      //qXq qXr rXr
//     u = u.cols(0, nrank-1);
//     V1 = u * v.t();
//     C1 = A1*V1.t();
//     // cout << "V1=" << V1 << endl;  //
// 
//     //// fix V, update A by mglasso
//     // subgradient method, vector version of the soft-thresholding rule: Chen&Huang2012, eq 13
//     if(mglasso_impute==1){    // impute off-diagonal Y
//       for(k=0; k<q; k++){
//         XYd.col(k) = XYlist.col(k) - XXlist(k)*C1.col(k);
//       }
//       // XY = XX*C1 + XYd;
//       // XR = XY*V1 - XX*A0;
//       XR = XYd*V1;             // pxr
//     }else{                     // ignore off-diagonal Y
//       XR = XYlist*V1;
//       for(k=0; k<q; k++){
//         XR = XR - XXlist(k)*A1*trans(V1.row(k))*V1.row(k);
//       }
//     }
//     // int it =  out_MGlasso["iter"];
//     A1 = MGlasso_ss0_Rcpp(XR, XX, lambda, A1, inner_conv, inner_iter);
//     // cout << "A1.n_nonzero=" << A1.n_nonzero << endl;
//     // if(A1.n_nonzero<3 | accu(abs(A1))<0.01) break;  // all entry are 0
//     C1 = A1*V1.t();
//     l2C1 = accu(square(C1));
// 
//     if (l2C1 <= 0.001) {
//       diff(iter) = 0;
//     } else {
//       iter = iter + 1;
//       diff(iter) = pow(accu(square(C0 - C1))/l2C1,0.5);
//     }
//     // cout << "diff=" << diff(iter)  << endl;
//   }
// 
//   diff = diff.subvec(0, iter);
//   sse = Yss;
//   for(k=0; k<q; k++){
//     sse = sse + as_scalar(arma::trans(C1.col(k))*XXlist(k)*C1.col(k)-2*arma::trans(XYlist.col(k))*C1.col(k));
//   }
//   dfu0 = A1.n_nonzero; // accu(A1 != 0);
//   dfv0 = accu(V1 != 0);
//   df = dfu0 * xrank/p + dfv0 - nrank*nrank;
//   double logqn = std::log(static_cast <double>(q * n));
//   double logsse = std::log(sse);
//   BIC = logsse + df*logqn/(q*n);
//   BICP = logsse + 2*df*logqn/(q*n);
//   AIC = logsse + 2*df/(q*n);
//   GIC = logsse + df*std::log(logqn)*std::log(static_cast <double>(p*q))/(q*n);
//   double dfqn2 = pow((1 - df/(q*n)), 2);
//   GCV = sse/(q*n*dfqn2);
//   if (diff(iter) <= conv) {
//     conv_flag = true;
//   } else {
//     conv_flag = false;
//   }
// 
//   // out["XXlist"] = XXlist;
//   // out["XX"] = XX;
//   out["diff"] = diff;
//   out["iter"] = iter;
//   out["BIC"]  = BIC;
//   out["BICP"] = BICP;
//   out["AIC"]  = AIC;
//   out["GCV"]  = GCV;
//   out["GIC"]  = GIC;
//   out["sse"]  = sse;
//   out["df"]   = df;
//   out["conv_flag"] = conv_flag;
//   out["A"]  = A1;
//   out["V"]  = V1;
//   // out["C"]  = C1;
//   return(out);
// }
// 












// Row-sparse Reduced-Rank Regression using summary statistics
//
// @param YX summary stat matrix: t(Y)*X
// @param XX summary stat matrix: t(X)*X
// @param n  sample size
// @param Yss summary stat scalar: sum square Y (F norm), for AIC/BIC etc
// @param method method
// @param A0 initial value
// @param V0 initial value
// @param nrank rank
// @param lambda tuning parameter
// @param conv conv
// @param maxiter maxiter
// @param inner_conv inner conv
// @param inner_iter inner iter
// @param WA weights
// @return estimation results
// // [[Rcpp::export]]
// Rcpp::List srrr_ss_Rcpp (arma::mat YX, arma::mat XX, 
//                          int n, double Yss,  
//                          String method, arma::mat A0, arma::mat V0,
//                          int nrank, double lambda,
//                          double conv, int maxiter,
//                          double inner_conv, int inner_iter,
//                          arma::vec WA) {
//   // min |Y-XAV|^2 + lamA*|A|, s.t. V'V = I_nrank
//   
//   // int n = Y.n_rows, p = X.n_cols, q = Y.n_cols;
//   // arma::mat YX = Y.t() * X;
//   //vec lamA = as_scalar(lambda)*ones(p);
//   
//   int p = YX.n_cols, q = YX.n_rows;
//   double xrank = accu(svd(XX) > 0.0001);
//   int iter=0;
//   double dfu0, dfv0;
//   bool conv_flag;
//   double l2C1, sse, df, BIC, BICP, AIC, GCV, GIC;
//   arma::vec s;
//   arma::vec diff(maxiter+1);
//   diff.fill(2*conv);
//   arma::mat mat1 = eye(q, q);   ///, iniU, iniD, iniV, iniC;
//   arma::mat V1;
//   arma::mat A1;
//   arma::mat C1;
//   arma::mat C0;
//   arma::mat W;
//   arma::mat u;
//   arma::mat v;
//   arma::mat residual;
//   
//   Rcpp::List ini;
//   Rcpp::List out_MGlasso;
//   Rcpp::List out;
//   
//   //ini = rrr_cpp(Y, X, nrank, false, mat1, true, true);
//   //arma::mat iniU = ini["U"];
//   //arma::vec inid = ini["D"];
//   //arma::mat iniV = ini["V"];
//   //arma::mat iniC = ini["C"];
//   //arma::mat iniD = diagmat(inid);
//   
//   //if (V0.is_empty() || A0.is_empty()) {
//   //  V1 = iniV;  //q r
//   //  A1 = iniU * iniD;  //p r
//   //  C1 = iniC; //p q
//   //} else {
//   V1 = V0;
//   A1 = A0;
//   C1 = A0*V0.t();
//   //}
//   
//   //  if (WA.is_empty()) {
//   //    if((method=="glasso")) {
//   //        WA = ones(p);
//   //    } else if((method=="adglasso")) {
//   //      A1 = iniU * iniD;  //p r
//   //      //vec A1norm(p);
//   //      for(int i=0;i<p;i++){
//   //        WA(i) = pow(sqrt(accu(square(A1.row(i)))),-wgamma);
//   //      }
//   //      //WA = pow(A1norm, -wgamma);
//   //    }
//   //  }
//   
//   while ((iter < maxiter) & (diff(iter) > conv)) {
//     //while (iter < 10) {
//     V0 = V1;
//     A0 = A1;
//     C0 = C1;
//     
//     arma::mat XYV0 = YX.t()*V0;
//     out_MGlasso = MGlasso_ss_Rcpp(XYV0, XX, lambda*WA, A0, inner_conv, inner_iter);
//     arma::mat MGlassoB = out_MGlasso["B"];
//     A1 = MGlassoB;   //p r
//     W = YX * A1;      //q r
//     svd(u, s, v, W);  //qXq qXr rXr
//     u = u.cols(0, nrank-1);
//     V1 = u * v.t();
//     C1 = A1*V1.t();
//     l2C1 = accu(square(C1));
//     if (l2C1 == 0) {
//       diff(iter) = 0;
//     } else {
//       iter = iter + 1;
//       diff(iter) = pow(accu(square(C0 - C1))/l2C1,0.5);
//     }
//     
//   }
//   
//   diff = diff.subvec(0, iter);
//   // residual = Y - X * C1;
//   // sse = accu(square(residual));
//   sse = Yss + trace(C1.t()*XX*C1-2*YX*C1);
//   dfu0 = accu(A1 != 0);
//   dfv0 = accu(V1 != 0);
//   df = dfu0 * xrank/p + dfv0 - nrank*nrank;
//   double logqn = std::log(static_cast <double>(q * n));
//   double logsse = std::log(sse);
//   BIC = logsse + df*logqn/(q*n);
//   BICP = logsse + 2*df*logqn/(q*n);
//   AIC = logsse + 2*df/(q*n);
//   GIC = logsse + df*std::log(logqn)*std::log(static_cast <double>(p*q))/(q*n);
//   double dfqn2 = pow((1 - df/(q*n)), 2);
//   GCV = sse/(q*n*dfqn2);
//   if (diff(iter) <= conv) {
//     conv_flag = true;
//   } else {
//     conv_flag = false;
//   }
//   
//   out["diff"] = diff;
//   out["iter"] = iter;
//   out["BIC"]  = BIC;
//   out["BICP"] = BICP;
//   out["AIC"]  = AIC;
//   out["GCV"]  = GCV;
//   out["GIC"]  = GIC;
//   out["sse"]  = sse;
//   out["df"]   = df;
//   out["conv_flag"] = conv_flag;
//   out["A"]  = A1;
//   out["V"]  = V1;
//   out["C"]  = C1;
//   return(out);
// }







// // [[Rcpp::export]]
// Rcpp::List MGlasso_Rcpp (arma::mat Y, arma::mat X, arma::vec lam, arma::mat B0, double conv, int maxiter) {
//   // min |Y-XB|^2 + lam*|B|
//   int  p=X.n_cols, iter=0, j;  // n=Y.n_rows, q=Y.n_cols,
//   double  diff=10*conv, l2B1, sse;
//   arma::rowvec sh;
//   //arma::mat mat1=eye(p,p);
//   arma::mat B1;
//   arma::mat res1;
//   arma::mat res1j;
//   arma::mat XRj;
//   Rcpp::List out;
//   if (lam.size() == 1) {lam = as_scalar(lam)*ones(p);}
//   sh = sum(square(X), 0);
//   //if (B0.is_finite()) {
//   B1 = B0;
//   //} else {
//   //  Rcpp::List ini = rrr_cpp(Y, X, 1, false, mat1, true, true);
//   //  arma::mat iniC = ini["C_ls"];
//   //  B1 = iniC;
//   //}
//   res1 = Y - X * B1;
//   while ((diff > conv) & (iter < maxiter)) {
//     // cout << "iter=" << iter  << endl; 
//     B0 = B1;
//     for (j = 0; j < p; j++) {
//       res1j = res1 +  X.col(j)* B1.row(j); //n q
//       XRj =   trans(X.col(j)) * res1j;    //1 q
//       // cout << "XRj=" << XRj  << endl; 
//       arma::rowvec t1=XRj/as_scalar(sh(j))*max(0.0,1-lam(j)/pow(accu(square(XRj)),0.5));
//       // cout << "t1=" << t1  << endl; 
//       B1.row(j) = t1;
//       res1 = res1j - X.col(j)* B1.row(j);
//     }
//     l2B1 = accu(square(B1));
//     if (l2B1 == 0) {
//       iter = maxiter;
//     } else {
//       diff = pow(accu(square(B0 - B1))/l2B1,0.5);
//       iter = iter + 1;
//     }
//   }
//   sse = accu(square(Y - X * B0));
//   out["B"] = B1;
//   out["sse"] = sse;
//   out["iter"] = iter;
//   return(out);
// }



// // [[Rcpp::export]]
// Rcpp::List msrrr_ss_Rcpp (arma::mat XYlist,     // pxq
//                           arma::mat XXlist,     // (qp)xp
//                           int n, double Yss,   // int p, int q,  String method, 
//                           arma::mat B0, // int nrank, 
//                           double lambda1,  // rank reduction
//                           double lambda2,  // sparsity
//                           int xrank, 
//                           double conv, 
//                           int maxiter,
//                           double h_eta) {
//   // mix-type sparse rrr
//   // Y = bdiag(Y1, ..., Yq)
//   // min |Y-XB|^2 + lam1*sum(svd(B)$d) + lam2*|B|_21 
//   
//   int p=XYlist.n_rows, q=XYlist.n_cols, k;   
//   // arma::mat XX = arma::zeros(p,p), XY, XYd = arma::zeros(p,q); 
//   // for(k=0; k<q; k++){ 
//   //   XX = XX + XXlist.rows(k*p,(k+1)*p-1);
//   // }
//   // double xrank = accu(svd(XX) > 0.0001);
//   int iter=0;
//   double dfu0, dfv0, eta=1, nrank;
//   bool conv_flag;
//   double Phi=Yss/n, l2B1, sse, df, BIC, BICP, AIC, GCV, GIC;
//   arma::vec d, bgnorm, thresh;  // s, 
//   arma::vec diff(maxiter+1);
//   diff.fill(2*conv);
//   arma::mat mat1 = eye(q, q);   
//   arma::mat C0=B0, G0=arma::zeros(p,q);
//   arma::mat C1=C0, B1=B0, G1=G0;
//   arma::mat Xdel, BG, Bm , u, v; 
//   Rcpp::List out;
//   
//   while ((iter < maxiter) & (diff(iter) > conv)) {
//     B0 = B1;
//     C0 = C1;
//     G0 = G1;
//     
//     // C-step
//     BG = B1 + G1/eta;
//     bgnorm = sqrt(sum(square(BG), 1)); 
//     thresh = (max(1 - lambda2 / eta / bgnorm, zeros(p)));
//     // C1 = thresh % BG;
//     for(k=0; k<q; k++){
//       C1.col(k) = thresh % BG.col(k); 
//     }
//     // cout << "C1: " << iter  << endl;  //
//     
//     // B-step
//     Xdel = zeros(p,q);  // XYlist - reshape(XXlist.t() * vec(B1), p, q);
//     for(k=0; k<q; k++){
//       Xdel.col(k) = XXlist.rows(k*p,(k+1)*p-1)*B1.col(k); 
//     }
//     Xdel = XYlist - Xdel;
//     // cout << "Xdel: " << Xdel.rows(0,4)   << endl;  //
//     Bm = (B1+eta*C1-G1)/(1+eta) + Xdel/sqrt(Phi);   //   *diag(1/sqrt(Phi1))
//     svd(u,d,v, Bm);
//     nrank = sum(d>lambda1/(1+eta)); 
//     B1 = u.cols(0,nrank-1) * diagmat(d.subvec(0,nrank-1)) * arma::trans(v.cols(0,nrank-1));
//     // cout << "B1: " << B1.rows(0,4)  << endl;  //
//     
//     // Phi-step: for now assume common dispersion (sigma^2) for all Y
//     sse = Yss;
//     for(k=0; k<q; k++){
//       sse = sse + as_scalar(arma::trans(B1.col(k))*XXlist.rows(k*p,(k+1)*p-1)*B1.col(k)-2*arma::trans(XYlist.col(k))*B1.col(k)); 
//     }
//     Phi = sse / n;
//     cout << "Phi: " << Phi  << endl;  //
//     
//     // Gam step
//     G1 = G1 + eta*(B1-C1);
//     eta = eta*h_eta;
//     // cout << "G1: " << iter  << endl;  //
//     
//     // conv
//     l2B1 = accu(square(B1));
//     if (l2B1 == 0) {
//       diff(iter) = 0;
//     } else {
//       iter = iter + 1;
//       diff(iter) = pow(accu(square(B1 - B0))/l2B1,0.5);
//     }
//   }
//   
//   
//   diff = diff.subvec(0, iter);
//   // sse = Yss;
//   // for(k=0; k<q; k++){
//   //   sse = sse + as_scalar(arma::trans(B1.col(k))*XXlist.rows(k*p,(k+1)*p-1)*B1.col(k)-2*arma::trans(XYlist.col(k))*B1.col(k)); 
//   // }
//   dfu0 = accu(B1 != 0)/q*nrank;
//   dfv0 = nrank*q; // accu(V1 != 0);
//   df = dfu0 * xrank/p + dfv0 - nrank*nrank;   
//   double logqn = std::log(static_cast <double>(q * n));
//   double logsse = std::log(sse);
//   BIC = logsse + df*logqn/(q*n);
//   BICP = logsse + 2*df*logqn/(q*n);
//   AIC = logsse + 2*df/(q*n);
//   GIC = logsse + df*std::log(logqn)*std::log(static_cast <double>(p*q))/(q*n);
//   double dfqn2 = pow((1 - df/(q*n)), 2);
//   GCV = sse/(q*n*dfqn2);
//   if (diff(iter) <= conv) {
//     conv_flag = true;
//   } else {
//     conv_flag = false;
//   }
//   
//   out["diff"] = diff;
//   out["iter"] = iter;
//   out["BIC"]  = BIC;
//   out["BICP"] = BICP;
//   out["AIC"]  = AIC;
//   out["GCV"]  = GCV;
//   out["GIC"]  = GIC;
//   out["sse"]  = sse;
//   out["df"]   = df;
//   out["conv_flag"] = conv_flag;
//   // out["A"]  = A1;
//   // out["V"]  = V1;
//   out["B"]  = B1;
//   return(out);
// }
// 


// Row-sparse Reduced-Rank Regression
//
// @param Y response matrix
// @param X covariate matrix
// @param method method
// @param A0 initial value
// @param V0 initial value
// @param nrank rank
// @param lambda tuning parameter
// @param conv conv
// @param maxiter maxiter
// @param inner_conv inner conv
// @param inner_iter inner iter
// @param WA weights
// @return estimation results
// // [[Rcpp::export]]
// Rcpp::List srrr_Rcpp (arma::mat Y, arma::mat X, String method, arma::mat A0, arma::mat V0,
//                       int nrank, double lambda,
//                       double conv, int maxiter,
//                       double inner_conv, int inner_iter,
//                       arma::vec WA) {
//   // min |Y-XAV|^2 + lamA*|A|, s.t. V'V = I_nrank
//   
//   int n = Y.n_rows, p = X.n_cols, q = Y.n_cols;
//   arma::mat YX = Y.t() * X;
//   //vec lamA = as_scalar(lambda)*ones(p);
//   
//   double xrank = accu(svd(X) > 0.01);
//   int iter=0;
//   double dfu0, dfv0;
//   bool conv_flag;
//   double l2C1, sse, df, BIC, BICP, AIC, GCV, GIC;
//   arma::vec s;
//   arma::vec diff(maxiter+1);
//   diff.fill(2*conv);
//   arma::mat mat1 = eye(q, q);   ///, iniU, iniD, iniV, iniC;
//   arma::mat V1;
//   arma::mat A1;
//   arma::mat C1;
//   arma::mat C0;
//   arma::mat W;
//   arma::mat u;
//   arma::mat v;
//   arma::mat residual;
//   
//   Rcpp::List ini;
//   Rcpp::List out_MGlasso;
//   Rcpp::List out;
//   
//   //ini = rrr_cpp(Y, X, nrank, false, mat1, true, true);
//   //arma::mat iniU = ini["U"];
//   //arma::vec inid = ini["D"];
//   //arma::mat iniV = ini["V"];
//   //arma::mat iniC = ini["C"];
//   //arma::mat iniD = diagmat(inid);
//   
//   //if (V0.is_empty() || A0.is_empty()) {
//   //  V1 = iniV;  //q r
//   //  A1 = iniU * iniD;  //p r
//   //  C1 = iniC; //p q
//   //} else {
//   V1 = V0;
//   A1 = A0;
//   C1 = A0*V0.t();
//   //}
//   
//   //  if (WA.is_empty()) {
//   //    if((method=="glasso")) {
//   //        WA = ones(p);
//   //    } else if((method=="adglasso")) {
//   //      A1 = iniU * iniD;  //p r
//   //      //vec A1norm(p);
//   //      for(int i=0;i<p;i++){
//   //        WA(i) = pow(sqrt(accu(square(A1.row(i)))),-wgamma);
//   //      }
//   //      //WA = pow(A1norm, -wgamma);
//   //    }
//   //  }
//   
//   while ((iter < maxiter) & (diff(iter) > conv)) {
//     //while (iter < 10) {
//     V0 = V1;
//     A0 = A1;
//     C0 = C1;
//     
//     arma::mat YV0 = Y*V0;
//     out_MGlasso = MGlasso_Rcpp(YV0, X, lambda*WA, A0, inner_conv, inner_iter);
//     arma::mat MGlassoB = out_MGlasso["B"];
//     A1 = MGlassoB;   //p r
//     W = YX * A1;      //q r
//     svd(u, s, v, W);  //qXq qXr rXr
//     u = u.cols(0, nrank-1);
//     V1 = u * v.t();
//     C1 = A1*V1.t();
//     l2C1 = accu(square(C1));
//     if (l2C1 == 0) {
//       diff(iter) = 0;
//     } else {
//       iter = iter + 1;
//       diff(iter) = pow(accu(square(C0 - C1))/l2C1,0.5);
//     }
//     
//   }
//   
//   diff = diff.subvec(0, iter);
//   residual = Y - X * C1;
//   sse = accu(square(residual));
//   dfu0 = accu(A1 != 0);
//   dfv0 = accu(V1 != 0);
//   df = dfu0 * xrank/p + dfv0 - nrank*nrank;
//   double logqn = std::log(static_cast <double>(q * n));
//   double logsse = std::log(sse);
//   BIC = logsse + df*logqn/(q*n);
//   BICP = logsse + 2*df*logqn/(q*n);
//   AIC = logsse + 2*df/(q*n);
//   GIC = logsse + df*std::log(logqn)*std::log(static_cast <double>(p*q))/(q*n);
//   double dfqn2 = pow((1 - df/(q*n)), 2);
//   GCV = sse/(q*n*dfqn2);
//   if (diff(iter) <= conv) {
//     conv_flag = true;
//   } else {
//     conv_flag = false;
//   }
//   
//   out["diff"] = diff;
//   out["iter"] = iter;
//   out["BIC"]  = BIC;
//   out["BICP"] = BICP;
//   out["AIC"]  = AIC;
//   out["GCV"]  = GCV;
//   out["GIC"]  = GIC;
//   out["sse"]  = sse;
//   out["df"]   = df;
//   out["conv_flag"] = conv_flag;
//   out["A"]  = A1;
//   out["V"]  = V1;
//   out["C"]  = C1;
//   return(out);
// }