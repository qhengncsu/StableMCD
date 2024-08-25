#include "trimmed_glm.h"
using namespace Rcpp;
using namespace std;
using namespace arma;
extern boost::random::rand48 rEngine;

inline double logit_inverse(double x) { return 1.0/(1.0 + std::exp(-x)); }
inline double log_inverse(double x) { return std::exp(x); }

// [[Rcpp::export]]
NumericVector irls_glm(const arma::mat& X, const arma::vec& y, const arma::vec& weights,
                       const std::string& family, const double tol = 1e-8, const int max_iter = 25){
  int n = X.n_rows;
  int p = X.n_cols;
  
  arma::vec beta = arma::zeros(p);
  arma::vec beta_new = arma::zeros(p); 
  arma::vec eta = X * beta;
  arma::vec mu, z, w;
  
  double (*link_inverse)(double);
  
  if (family == "binomial") {
    link_inverse = logit_inverse;
  } else if (family == "poisson") {
    link_inverse = log_inverse;
  } else {
    Rcpp::stop("Unsupported family. Use 'binomial' or 'poisson'.");
  }
  arma::mat X_w;
  arma::vec z_w;
  for (int iter = 0; iter < max_iter; ++iter) {
    mu = arma::vec(n);
    for (int i = 0; i < n; ++i) {
      mu(i) = link_inverse(eta(i));
    }
    if (family == "binomial") {
      w = weights % mu % (1.0 - mu);
      z = eta + (y - mu) / (mu % (1.0 - mu));
    } else { // poisson
      w = weights % mu;
      z = eta + (y - mu) / mu;
    }
    
    X_w = diagmat(arma::sqrt(w)) *X;
    z_w = z % arma::sqrt(w);
    bool solving_success = arma::solve(beta_new, X_w.t() * X_w, X_w.t() * z_w);
    if (!solving_success) {
      Rcpp::warning("Failed to solve the system. The model may be overfitted or the data may be ill-conditioned.");
      break;
    }
    
    if (arma::norm(beta_new - beta) < tol) {
      beta = beta_new;
      break;
    }
    
    beta = beta_new;
    eta = X * beta;
  }
  
  return wrap(beta);
}

// [[Rcpp::export]]
NumericVector var_stab_res(const arma::mat& X, const arma::vec& y, const arma::vec& weights,
                           const std::string& family, const arma::vec& beta){
  int p = X.n_cols;
  if (family == "binomial") {
    arma::vec y_transformed = sqrt(weights)%asin(sqrt(y));
    arma::vec phat = 1/(1+exp(-X*beta(span(1,p))-beta(0)));
    arma::vec residual = y_transformed - sqrt(weights)%asin(sqrt(phat)) + (0.125*(1-2*phat))/sqrt(weights%phat%(1-phat));
    return wrap(residual);
  }else if (family == "poisson"){
    arma::vec y_transformed = sqrt(y);
    arma::vec lambdahat = exp(X*beta(span(1,p))+beta(0));
    arma::vec residual = y_transformed - sqrt(lambdahat) + 1/(8*sqrt(lambdahat));
    return wrap(residual);
  }else {
    Rcpp::stop("Unsupported family. Use 'binomial' or 'poisson'.");
  }
}

// [[Rcpp::export]]
List trimmed_glm(const arma::mat& X, const arma::colvec& y, const arma::colvec& alphas, 
                 const std::string& family, const arma::colvec& weights){
  int n = X.n_rows;
  int p = X.n_cols + 1;
  arma::mat Xy = arma::join_horiz(X,y);
  //copy the contents of Xy to a new double vector
  //traverse row by row!
  int h = 0;
  arma::uvec init_index = arma::zeros<arma::uvec>(floor(0.5*n));
  arma::uvec order = arma::zeros<arma::uvec>(n);
  arma::colvec depths = proj_depth(Xy, Xy);
  order = arma::sort_index(depths, "descend");
  init_index = order.head(floor(0.5*n));
  arma::mat X_with_intercept = arma::join_horiz(arma::ones(n), X);
  arma::colvec beta_init;
  if (family == "binomial" || family == "poisson") {
    beta_init = irls_glm(X_with_intercept.rows(init_index),y.elem(init_index),weights.elem(init_index),family);
  }else {
    Rcpp::stop("Unsupported family. Use 'binomial' or 'poisson'.");
  }
  arma::colvec beta = beta_init;
  arma::colvec beta_old = beta_init;
  arma::colvec res = arma::zeros<arma::colvec>(n);
  Rcpp::List betas(alphas.n_elem);
  Rcpp::List indexes(alphas.n_elem);
  for (int i = 0; i < alphas.n_elem; i++){
    double alpha = alphas(i);
    h = floor(alpha * n);
    arma::uvec index = arma::zeros<arma::uvec>(h);
    int j = 0;
    double rel_change = arma::datum::inf;
    beta = beta_init;
    beta_old = beta;
    while (j <= 100 && rel_change > 1e-6) {
      if (family == "binomial" || family == "poisson") {
        res = var_stab_res(X,y,weights,family,beta);
      }else {
        Rcpp::stop("Unsupported family. Use 'binomial' or 'poisson'.");
      }
      order = arma::sort_index(abs(res), "ascend");
      index = order.head(h);
      if (family == "binomial" || family == "poisson") {
        beta_init = irls_glm(X_with_intercept.rows(index),y.elem(index),weights.elem(index),family);
      }else {
        Rcpp::stop("Unsupported family. Use 'binomial' or 'poisson'.");
      }
      rel_change = arma::norm(beta - beta_old, 2) / (arma::norm(beta_old, 2) + 1);
      beta_old = beta;
      j++;
    }
    betas[i] = beta;
    indexes[i] = index;
  }
  return List::create(Named("betas") = betas,
                      Named("indexes") = indexes);
}

// [[Rcpp::export]]
List bootstrap_glm(const arma::mat& X, const arma::colvec& y, const arma::colvec& alphas, 
                   const std::string& family, const arma::colvec& weights, int B = 50){
  int n = X.n_rows;
  int p = X.n_cols;
  int nalpha = alphas.n_elem;
  arma::mat X1 = arma::mat(X);
  arma::mat X2 = arma::mat(X);
  arma::colvec y1 = y;
  arma::colvec y2 = y;
  arma::colvec insta_medians = arma::zeros<arma::colvec>(nalpha);
  arma::colvec insta_sds = arma::zeros<arma::colvec>(nalpha);
  arma::mat instas = arma::mat(nalpha, B);
  std::uniform_int_distribution<int> dis(0, n-1);
  int index1, index2;
  arma::colvec beta1,beta2;
  arma::colvec res1,res2;
  arma::uvec order1 = arma::zeros<arma::uvec>(n);
  arma::uvec order2 = arma::zeros<arma::uvec>(n);
  arma::ivec is_outlier1, is_outlier2;
  int h;
  for (int b = 0; b < B; b++){
    for (int i = 0; i<n; i++){
      index1 = dis(rEngine);
      index2 = dis(rEngine);
      X1.row(i) = X.row(index1);
      X2.row(i) = X.row(index2);
      y1(i) = y(index1);
      y2(i) = y(index2);
    }
    List L1 = trimmed_glm(X1,y1,alphas,family,weights);
    List L2 = trimmed_glm(X2,y2,alphas,family,weights);
    List betas1 = L1["betas"];
    List betas2 = L2["betas"];
    for(int i = 0; i<nalpha; i++){
      double alpha = alphas(i);
      h = floor(alpha * n);
      if (h <  floor(0.5 * n) || h >= n){
        throw std::runtime_error("Error： invalid h, h = " + std::to_string(h));
      }
      is_outlier1 = arma::ones<arma::ivec>(n);
      is_outlier2 = arma::ones<arma::ivec>(n);
      beta1 = as<arma::vec>(wrap(betas1[i]));
      beta2 = as<arma::vec>(wrap(betas2[i]));
      if (family == "binomial" || family == "poisson") {
        res1 = var_stab_res(X,y,weights,family,beta1);
        res2 = var_stab_res(X,y,weights,family,beta2);
      }else {
        Rcpp::stop("Unsupported family. Use 'binomial' or 'poisson'.");
      }
      order1 = arma::sort_index(abs(res1), "ascend");
      order2 = arma::sort_index(abs(res2), "ascend");
      is_outlier1.elem(order1.head(h)) = arma::zeros<arma::ivec>(h);
      is_outlier2.elem(order2.head(h)) = arma::zeros<arma::ivec>(h);
      double insta = get_insta(is_outlier1, is_outlier2, h);
      instas(i, b) = insta;
    }
    Rcout << "Bootstrap pair " << b+1 << " completed!" << std::endl;
  }
  insta_medians = arma::median(instas,1);
  insta_sds = arma::stddev(instas, 0, 1);
  double best_alpha = alphas(insta_medians.index_min());
  return List::create(Named("best_alpha") = best_alpha,
                      Named("insta_medians") = insta_medians,
                      Named("insta_sds") = insta_sds,
                      Named("alphas") = alphas);
}
