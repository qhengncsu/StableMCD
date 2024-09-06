#include "mcd.h"
using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::export]]
List concentration(const arma::mat& X, int h, arma::uvec index, bool verbose=true){
  int n = X.n_rows;
  double old_sum_distances = numeric_limits<double>::infinity();
  if(index.n_elem != h){
    Rcpp::stop("Initial index set length does not equal specified subset size h.");
  }
  int iter = 0;
  arma::uvec new_index = index;
  arma::vec mu;
  arma::mat SigmaInv;
  arma::mat X_subset;
  arma::mat X_centered;
  arma::vec MDs(n);
  arma::uvec order;
  arma::vec row;
  while(iter < 100){
    X_subset = X.rows(new_index);
    mu = arma::mean(X_subset,0).as_col();
    SigmaInv = inv(cov(X_subset,1));
    X_centered = X.each_row() - mu.as_row();
    for(int i=0;i<n;i++){
      row = X_centered.row(i).as_col();
      MDs(i) = dot(row,SigmaInv*row);
    }
    order = arma::sort_index(MDs);
    new_index = order.head(h);
    double sum_distances = sum(MDs.elem(new_index));
    if(verbose){
      Rcout << "Iteration " << iter+1 << " log determinant: "<<real(-log_det(SigmaInv))<<endl;
    }
    if(std::abs(sum_distances-old_sum_distances) < 1e-4){
      break;
    }
    old_sum_distances = sum_distances;
    iter += 1;
  }
  return List::create(Named("index") = new_index,
                      Named("muhat") = mu,
                      Named("Sigmahat") = cov(X_subset,1));
}

// [[Rcpp::export]]
List mcd(const arma::mat& X, const arma::colvec& alphas, bool csteps = true, int direction_style=3){
  int n = X.n_rows;
  int p = X.n_cols;
  int h = 0;
  arma::colvec depths = proj_depth(X, X, direction_style);
  arma::uvec order = arma::sort_index(depths, "descend");
  arma::uvec index;
  Rcpp::List indexes(alphas.n_elem);
  Rcpp::List muhats(alphas.n_elem);
  Rcpp::List Sigmahats(alphas.n_elem);
  for (int i = 0; i < alphas.n_elem; i++){
    double alpha = alphas(i);
    h = floor(alpha * n);
    index = order.head(h);
    if(csteps){
      List cstep_result = concentration(X,h,index);
      indexes[i] = cstep_result["index"];
      muhats[i] = cstep_result["muhat"];
      Sigmahats[i] = cstep_result["Sigmahat"];
    }else{
      indexes[i] = index;
      muhats[i] = arma::mean(X.rows(index),0);
      Sigmahats[i] = cov(X.rows(index),1);
    }
  }
  return List::create(Named("indexes") = indexes,
                      Named("muhats") = muhats,
                      Named("Sigmahats") = Sigmahats);
}

// [[Rcpp::export]]
List bootstrapMcd(const arma::mat& X, const arma::colvec& alphas, bool csteps = true, int direction_style=3, int B=100){
  int n = X.n_rows;
  int p = X.n_cols;
  int nalpha = alphas.n_elem;
  arma::mat X1 = arma::mat(X);
  arma::mat X2 = arma::mat(X);
  arma::colvec insta_means = arma::zeros<arma::colvec>(nalpha);
  arma::colvec insta_sds = arma::zeros<arma::colvec>(nalpha);
  arma::mat instas = arma::mat(nalpha, B);
  std::uniform_int_distribution<int> dis(0, n-1);
  int idx1, idx2;
  arma::uvec order1 = arma::zeros<arma::uvec>(n);
  arma::uvec order2 = arma::zeros<arma::uvec>(n);
  arma::uvec order1_X = arma::zeros<arma::uvec>(n);
  arma::uvec order2_X = arma::zeros<arma::uvec>(n);
  arma::ivec is_outlier1, is_outlier2;
  int h;
  std::random_device rd;
  std::mt19937 gen(rd());
  vec depths = proj_depth(X, X, direction_style);
  Rcout << depths(0)<<endl;
  vec depths1(n), depths2(n);
  uvec index1, index2;
  for (int b = 0; b < B; b++){
    for (int i = 0; i<n; i++){
      idx1 = dis(gen);
      idx2 = dis(gen);
      while (idx1==idx2){
        idx2 = dis(gen);
      }
      X1.row(i) = X.row(idx1);
      X2.row(i) = X.row(idx2);
      depths1(i) = depths(idx1);
      depths2(i) = depths(idx2);
    }
    //depths1 = proj_depth(X1, X1, direction_style);
    //depths2 = proj_depth(X2, X2, direction_style);
    order1 = arma::sort_index(depths1, "descend");
    order2 = arma::sort_index(depths2, "descend");
    for(int i = 0; i<nalpha; i++){
      double alpha = alphas(i);
      h = floor(alpha * n);
      if (h <  floor(0.5 * n) || h >= n){
        throw std::runtime_error("Error： invalid h, h = " + std::to_string(h));
      }
      index1 = order1.head(h);
      index2 = order2.head(h);
      if(csteps){
        List cstep_result1 = concentration(X1,h,index1,false);
        List cstep_result2 = concentration(X2,h,index2,false);
        index1 = as<arma::uvec>(wrap(cstep_result1["index"]));
        index2 = as<arma::uvec>(wrap(cstep_result2["index"]));
      }
      depths1 = proj_depth(X,X1.rows(index1),1);
      depths2 = proj_depth(X,X2.rows(index2),1);
      order1_X = arma::sort_index(depths1, "descend");
      order2_X = arma::sort_index(depths2, "descend");
      is_outlier1 = arma::ones<arma::ivec>(n);
      is_outlier2 = arma::ones<arma::ivec>(n);
      is_outlier1.elem(order1_X.head(h)) = arma::zeros<arma::ivec>(h);
      is_outlier2.elem(order2_X.head(h)) = arma::zeros<arma::ivec>(h);
      double insta = get_insta(is_outlier1, is_outlier2, h);
      instas(i, b) = insta;
    }
    if(b%10 == 9){
      Rcout << "Bootstrap pair " << b+1 << " completed!" << std::endl;
    }
  }
  for(int i = 0; i<nalpha; i++){
    insta_means(i) = arma::mean(instas.row(i).as_col());
  }
  insta_sds = arma::stddev(instas, 0, 1);
  double best_alpha = alphas(insta_means.index_min());
  return List::create(Named("best_alpha") = best_alpha,
                      Named("insta_means") = insta_means,
                      Named("insta_sds") = insta_sds,
                      Named("alphas") = alphas,
                      Named("instas") = instas);
}
