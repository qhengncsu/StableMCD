#include <RcppArmadillo.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>

Rcpp::NumericVector proj_depth(const arma::mat& X, const arma::mat& data, int style);
double trimean(arma::vec x);

// extern "C" { 
void ProjectionDepth(double *points, double *objects, int *numObjects,
                     int *dimension, int *cardinalities, int *numClasses,
                     double *directions, double *projections, int *k,
                     int *newDirs, int *seed, double *depths, int *dirStyle); 

double Choose(int n, int k);
double get_insta(const arma::ivec& is_outlier1, const arma::ivec& is_outlier2, int h);