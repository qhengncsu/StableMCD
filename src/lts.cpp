#include "lts.h"
using namespace Rcpp;
using namespace std;
extern boost::random::rand48 rEngine;



// [[Rcpp::export]]
List lts(const arma::mat& X, const arma::colvec& y, const arma::colvec& alphas){
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
  arma::colvec beta_init = arma::solve(X_with_intercept.rows(init_index), y.elem(init_index));
  arma::colvec beta = beta_init;
  arma::colvec beta_old = beta_init;
  arma::colvec likelihood = arma::zeros<arma::colvec>(n);
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
      likelihood = -arma::square(y - X_with_intercept * beta);
      order = arma::sort_index(likelihood, "descend");
      index = order.head(h);
      beta = arma::solve(X_with_intercept.rows(index), y.elem(index));
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

double Choose(int n, int k){
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;
  double result = n;
  for( int i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}

double get_insta(const arma::ivec& is_outlier1, const arma::ivec& is_outlier2, int h){
  double n = is_outlier1.n_elem;
  double pX = sum(abs(is_outlier1-is_outlier2))/n;
  double c = (Choose(h,2)+Choose(n-h,2))/Choose(n,2);
  return pX*(1-pX)/(c*(1-c))-1;
}

// [[Rcpp::export]]
List bootstrap_lts(const arma::mat& X, const arma::colvec& y, const arma::colvec& alphas, int B = 50){
  int n = X.n_rows;
  int p = X.n_cols;
  int nalpha = alphas.n_elem;
  arma::mat X1 = arma::mat(X);
  arma::mat X2 = arma::mat(X);
  arma::colvec y1 = y;
  arma::colvec y2 = y;
  arma::colvec insta_means = arma::zeros<arma::colvec>(nalpha);
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
    List L1 = lts(X1,y1,alphas);
    List L2 = lts(X2,y2,alphas);
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
      res1 = abs(y - X * beta1(arma::span(1,p)) - beta1[0]);
      res2 = abs(y - X * beta2(arma::span(1,p)) - beta2[0]);
      order1 = arma::sort_index(res1, "ascend");
      order2 = arma::sort_index(res2, "ascend");
      is_outlier1.elem(order1.head(h)) = arma::zeros<arma::ivec>(h);
      is_outlier2.elem(order2.head(h)) = arma::zeros<arma::ivec>(h);
      double insta = get_insta(is_outlier1, is_outlier2, h);
      instas(i, b) = insta;
    }
    Rcout << "Bootstrap pair " << b+1 << " completed!" << std::endl;
  }
  insta_means = arma::mean(instas,1);
  insta_sds = arma::stddev(instas, 0, 1);
  double best_alpha = alphas(insta_means.index_min());
  return List::create(Named("best_alpha") = best_alpha,
                      Named("insta_means") = insta_medians,
                      Named("insta_sds") = insta_sds,
                      Named("alphas") = alphas);
}