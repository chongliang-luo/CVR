// the negative loglikelihood function based on the broslow's formula


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include "utils.h"
 
// [[Rcpp::export]]
arma::mat rcpp_coxph_deriv_lp(double beta,
                              arma::vec time,  // sorted
                              arma::vec event,
                              arma::vec z)  // lp vector the same len as time
{
  // sort event and z based on time
  // arma::uvec s_time_ind = arma::sort_index(time_);
  // arma::vec time = time_.elem(s_time_ind);
  // arma::vec event = event_.elem(s_time_ind);
  // arma::vec z = z_.elem(s_time_ind); 
  
  arma::vec exp_z_beta = arma::exp(z * beta);
  arma::vec h0_denom = Intsurv::cum_sum(exp_z_beta, true);
  arma::vec cumsum_h0_denom_inv = Intsurv::cum_sum(1/h0_denom, false);
  arma::vec grad = beta*(event - exp_z_beta % cumsum_h0_denom_inv);
  arma::vec cumsum_h0_denom_inv_2 = Intsurv::cum_sum(1/square(h0_denom), false); 
  arma::vec hess = (exp_z_beta % cumsum_h0_denom_inv - arma::square(exp_z_beta) % cumsum_h0_denom_inv_2) * beta * beta;
  // arma::mat deriv = join_cols(arma::conv_to<arma::vec>::from(s_time_ind), grad, hess);
  arma::mat deriv = join_cols(grad, hess);
  
  return deriv;
}


