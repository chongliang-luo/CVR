#include <RcppArmadillo.h> 
//#include <Rcpp.h>
#include <R.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;
using namespace std;
using namespace arma;
//[[Rcpp::plugins("cpp11")]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::mat list1(Rcpp::List Xlist) {  // crash!
  arma::mat X = as<mat>(Xlist[1]); 
  return(X );
}

// [[Rcpp::export]]
arma::mat array1(arma::cube Xcube) {  //  
  arma::mat X = Xcube.subcube(0, 0, 0, 2, 2, 0); 
  return(X );
}


// [[Rcpp::export]]
arma::mat X2(arma::mat X) { 
  arma::mat X1 = X.rows(0,3).cols(0,3);
  return(X1  * 2);
}

// [[Rcpp::export]]
arma::mat MM(arma::mat A, arma::mat B) {   // similar time cost
  return(A*B);
}


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
Rcpp::List fcpp(arma::mat A0, arma::cube Q 
                ){ 
  // , arma::sp_mat A0sp 
  // arma::cube XXlist,  
  // arma::mat XX, 
  // arma::sp_mat XXsp, 
  //arma::vec x, arma::sp_mat A, arma::mat A0, 
  
  Rcpp::List out;
  // arma::sp_mat A0(20, 4);
  // arma::cube Apath(40, 3, 2); // , fill::zeros
  // Apath.slice(0) = sp_mat(Apath.slice(0));  # cube can't be sparse?
  // Apath.slice(1) = sp_mat(Apath.slice(1));
  // arma::mat XX = as<mat>(XXlist[0] + XXlist[1]);
  // Rcpp::NumericMatrix XX0 = XXlist[0];
  // Rcpp::NumericMatrix XX1 = XXlist[1];
  // sp_mat XX(3000,3000);  
  // int q=XXlist.n_slices;
  // arma::field<sp_mat> XXlist0(q);
  // for(int k=0; k<q; k++){ 
  //   XXlist0(k) =  XXlist.slice(k);
  //   XX = XX + XXlist0(k); // as<arma::sp_mat>(XXlist.slice(k));
  // }
  // sp_mat XA = XXsp*A0sp;
  // mat U;
  // vec s;
  // mat V;
  // svds(u,s,v,XA)
  
  mat norm21_A0 = norm21_rcpp(A0);
  uvec nz = find(norm21_A0!=0);
  // sp_mat A00 = A0.submat(2,3,4,5); 
  IntegerVector tt = seq_len( A0.n_rows)-1;
  // cout << "tt" << tt << endl;
  // cout << "nz" << nz << endl;
  IntegerVector tt2 = as<IntegerVector>(wrap(nz));
  // cout << "tt2" << tt2 << endl;
  IntegerVector nz_ = setdiff(tt, tt2);
  uvec nz0 = as<uvec>(wrap(nz_));
  // cout << "nz0" << nz0 << endl;
  mat A00 = A0.rows(nz0);
  // sp_mat A01 = A0(span(0,12), span::all);
  // arma::mat XYd(A0.n_rows, A0.n_cols, fill::zeros)
  // arma::mat XYd_nz_ = XYd.submat(nz_,1);
  
  mat sumQ0 = sum(Q, 0);
  mat sumQ1 = sum(Q, 1);
  mat sumQ2 = sum(Q, 2);
  uvec id0 = { 0,2},id1 = { 0,1},id2 = {1};
  int k=1;
  // uvec tt1(1,fill::value(k));
  A0.submat(id0, uvec tt1(1,fill::value(k))).fill(0.1);
  
  mat subQ = Q.slice(1).submat(id0, id1);
  cube Qtube(id0.n_elem, id1.n_elem, id2.n_elem);
  for(int k=0; k<id2.n_elem; k++) Qtube.slice(k) = Q.slice(k).submat(id0,id1);
  // cube Qtube = Q.each_slice( [&](mat& X){ X.submat({ 0,2},{ 0,1}); } );
  // out["order"] = sort_index(x);
  // arma::uvec tt = sort_index(x, "descend");
  // out["which_max"] = tt(0);
  // // out["norm21"] = sqrt(sum(square(A),1));
  // arma::sp_mat A2 = sqrt(sum(square(A),1)); 
  // out["norm21"] = A2;
  // // out["sp_cube"] = Apath;
  // 
  // // field with sparse mat
  // arma::field<arma::sp_mat> Apath(2);
  // // Apath.for_each( [&](arma::sp_mat& X) { X.set_size(40, 3); } );
  // 
  // Apath(0) = sp_mat(A0); 
  // Apath(1) = sp_mat(A0);   
  // 
  // out["sp_field"] = Apath;
  // sp_mat x = XXlist0(0).submat(0,0,400,3);
  // arma::field XXlist1=XXlist0.for_each( [&](arma::sp_mat& X) { X=X.submat(295,295,303,303); } );
  // out["sp_field"] = XXlist0;
  // out["sp_field_sum"] = XX;
  // out["sp_field_dim"] = x.n_cols+x.n_rows; // sum(x, 1);
  // out["XYd_nz_"] = XYd_nz_;
  // out["norm21_A0"] = norm21_A0;
  // out["nz"] = nz;
  // out["nz2"] = nz2;
  // out["tt"] = tt;
  // out["A01"] = A01;
  out["A00"] = A00;
  out["A0new"] = A0;
  out["sumQ0"] = sumQ0;
  out["sumQ1"] = sumQ1;
  out["sumQ2"] = sumQ2;
  out["subQ"] = subQ;
  out["Qtube"] = Qtube;
  return(out);
}

// Apath = array(A, c(40,4,2))
// tt = fcpp(x=c(2,5,3,1,9), A, A0) 
// tt = fcpp(x=c(2,5,3,1,9), A, A0, XXlist0) 
// tt = fcpp(XXlist, A0)
// tt$sp_field[[1]][295:303,295:303]
// tt$ sp_field_sum [295:303,295:303]
// Q = array(1:6, c(3,2,2))
  // tt = fcpp(A0, Q)
  
  
// [[Rcpp::export]]
arma::sp_mat ff (arma::sp_mat A0, arma::mat V0){
  
  sp_mat C0=sp_mat(A0*V0.t());
  
  return(C0);
}

// [[Rcpp::export]]
arma::sp_mat  spmat_test(arma::cube  M, arma::mat B, int sp) {
  cube Ms=M;
  int p = M.n_rows, q = M.n_slices, k;
  arma::sp_mat out;
  if(sp==1){
    sp_mat sMs(p,p);
     
    for(k=0; k<q; k++){
      M.slice(k) = sp_mat(M.slice(k));
      sMs = sMs + M.slice(k);
    }
    // sMs(99,101) = 2;
    // sMs(101,101) = 0;
    out = sMs*B;
    // sp_mat Ms = sp_mat(M);
    // sp_mat Bs = sp_mat(B);
    // out = arma::spsolve(Ms, B);     // ARMA_USE_SUPERLU must be enabled in config.hpp ...?
    // out = Ms*B;                     // matrix multiply 3x faster if M is 10% dense

    // mat U;
    // vec s;
    // svds(U, s, out, Ms, 10);           // 40x faster
    
    // out = Ms.cols(0,9);
    // out.rows(0,1).zeros(); 
  }else{
    // sp_mat Ms = sp_mat(M);
    // out = solve(M, B); 
    // mat out = M*B;
    mat sMs=zeros(p,p);
    
    for(k=0; k<q; k++){
      // M.slice(k) = sp_mat(M.slice(k));
      sMs = sMs + M.slice(k);
    }
    // sMs(99,101) = 2;
    // sMs(101,101) = 0;
    out = sMs*B;
    // mat U;
    // vec s;
    // svd(U, s, out, M);
  }
  
  return(out);
}
