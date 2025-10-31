#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda){
  // Your function code goes here
  if (a > lambda) {
    return a - lambda;
  } else if (a < -lambda) {
    return a + lambda;
  } else {
    return 0;
  }
}

// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Your function code goes here
  
  // Get rows of Xtilde and calculate residuals
  int n = Xtilde.n_rows;
  arma::colvec r = Ytilde - Xtilde * beta;
  
  // Calculate objective function
  double fobj = (1.0 / (2.0 * n)) * arma::dot(r, r) + lambda * arma::sum(arma::abs(beta));
  
  return fobj;
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  // Your function code goes here
  int n = Xtilde.n_rows;
  int p = Xtilde.n_cols;
  
  // Check for starting point beta_start
  // If none supplied, initialize with vector of zeroes
  arma::colvec beta;
  if (beta_start.is_empty()){
    beta = arma::zeros<arma::colvec>(p);
  } else{
    beta = beta_start;
  }
  
  // Calculate residuals
  arma::colvec r = Ytilde - Xtilde * beta;
  arma::rowvec denom = arma::sum(arma::square(Xtilde), 0) / n;
  
  // Set error for while loop
  double err = 1000.00;
  double obj_old = lasso_c(Xtilde, Ytilde, beta, lambda);
  
  while (err> eps){
    for(int j = 0; j < p; ++j){
      arma::colvec xj = Xtilde.col(j);
      double bj_old = beta(j);
      
      // Calculate rho_j
      double rho_j = arma::dot(xj, r + xj * bj_old) / n;
      
      // Use soft thresholding to calculating now beta_j
      double bj_new = soft_c(rho_j, lambda) / denom(j);
      
      // Recalculate residual
      if (bj_new != bj_old){
        double delta = bj_new - bj_old;
        r -= xj * delta;
        beta(j) = bj_new;
      }
    }
    
    // Calculate new fmin using objective function
    double fmin = lasso_c(Xtilde, Ytilde, beta, lambda);
    
    //Recalculate error
    err = obj_old - fmin;
    obj_old = fmin;
  }
  return beta;
}  

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Your function code goes here
  
  // Get necessary sizes of Xtilde and lambda_seq
  int p = Xtilde.n_cols;
  int l = lambda_seq.n_elem;
  
  // Initialize empty beta matrix and fmin vector
  arma::mat beta_mat(p, l, arma::fill::zeros);
  arma::vec fmin_vec = arma::zeros<arma::vec>(l);
  
  // Initialize beta_start
  arma::colvec beta_start = arma::zeros<arma::colvec>(p);
  
  // Apply fitlassostandardized_c going from largest to smallest lambda
  // Utilize warm starts method
  for(int i = 0; i < l; ++i){
    double lambda_i= lambda_seq(i);
    arma::colvec beta_i = fitLASSOstandardized_c(Xtilde, Ytilde, lambda_i, beta_start, eps);
    beta_mat.col(i) = beta_i;
    beta_start = beta_i;
  }
  
  return(beta_mat);
}