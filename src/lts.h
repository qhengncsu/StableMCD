#include <RcppArmadillo.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>
Rcpp::List lts(const arma::mat& X, const arma::colvec& y, const arma::colvec& alphas);
Rcpp::NumericVector proj_depth(const arma::mat& X, const arma::mat& data, int style);

// extern "C" { 
void ProjectionDepth(double *points, double *objects, int *numObjects,
                            int *dimension, int *cardinalities, int *numClasses,
                            double *directions, double *projections, int *k,
                            int *newDirs, int *seed, double *depths, int *dirStyle); 
// }