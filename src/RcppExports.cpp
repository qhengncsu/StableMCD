// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// proj_depth
Rcpp::NumericVector proj_depth(const arma::mat& X, const arma::mat& data);
RcppExport SEXP _StableMCD_proj_depth(SEXP XSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(proj_depth(X, data));
    return rcpp_result_gen;
END_RCPP
}
// lts
List lts(const arma::mat& X, const arma::colvec& y, const arma::colvec& alphas);
RcppExport SEXP _StableMCD_lts(SEXP XSEXP, SEXP ySEXP, SEXP alphasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type alphas(alphasSEXP);
    rcpp_result_gen = Rcpp::wrap(lts(X, y, alphas));
    return rcpp_result_gen;
END_RCPP
}
// trimean
double trimean(arma::vec x);
RcppExport SEXP _StableMCD_trimean(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(trimean(x));
    return rcpp_result_gen;
END_RCPP
}
// bootstrap_lts
List bootstrap_lts(const arma::mat& X, const arma::colvec& y, const arma::colvec& alphas, int B);
RcppExport SEXP _StableMCD_bootstrap_lts(SEXP XSEXP, SEXP ySEXP, SEXP alphasSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap_lts(X, y, alphas, B));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _StableMCD_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// irls_glm
NumericVector irls_glm(const arma::mat& X, const arma::vec& y, const arma::vec& weights, const std::string& family, const double tol, const int max_iter);
RcppExport SEXP _StableMCD_irls_glm(SEXP XSEXP, SEXP ySEXP, SEXP weightsSEXP, SEXP familySEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family(familySEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(irls_glm(X, y, weights, family, tol, max_iter));
    return rcpp_result_gen;
END_RCPP
}
// var_stab_res
NumericVector var_stab_res(const arma::mat& X, const arma::vec& y, const arma::vec& weights, const std::string& family, const arma::vec& beta);
RcppExport SEXP _StableMCD_var_stab_res(SEXP XSEXP, SEXP ySEXP, SEXP weightsSEXP, SEXP familySEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family(familySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(var_stab_res(X, y, weights, family, beta));
    return rcpp_result_gen;
END_RCPP
}
// trimmed_glm
List trimmed_glm(const arma::mat& X, const arma::colvec& y, const arma::colvec& alphas, const std::string& family, const arma::colvec& weights);
RcppExport SEXP _StableMCD_trimmed_glm(SEXP XSEXP, SEXP ySEXP, SEXP alphasSEXP, SEXP familySEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family(familySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(trimmed_glm(X, y, alphas, family, weights));
    return rcpp_result_gen;
END_RCPP
}
// bootstrap_glm
List bootstrap_glm(const arma::mat& X, const arma::colvec& y, const arma::colvec& alphas, const std::string& family, const arma::colvec& weights, int B);
RcppExport SEXP _StableMCD_bootstrap_glm(SEXP XSEXP, SEXP ySEXP, SEXP alphasSEXP, SEXP familySEXP, SEXP weightsSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family(familySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap_glm(X, y, alphas, family, weights, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StableMCD_proj_depth", (DL_FUNC) &_StableMCD_proj_depth, 2},
    {"_StableMCD_lts", (DL_FUNC) &_StableMCD_lts, 3},
    {"_StableMCD_trimean", (DL_FUNC) &_StableMCD_trimean, 1},
    {"_StableMCD_bootstrap_lts", (DL_FUNC) &_StableMCD_bootstrap_lts, 4},
    {"_StableMCD_rcpp_hello_world", (DL_FUNC) &_StableMCD_rcpp_hello_world, 0},
    {"_StableMCD_irls_glm", (DL_FUNC) &_StableMCD_irls_glm, 6},
    {"_StableMCD_var_stab_res", (DL_FUNC) &_StableMCD_var_stab_res, 5},
    {"_StableMCD_trimmed_glm", (DL_FUNC) &_StableMCD_trimmed_glm, 5},
    {"_StableMCD_bootstrap_glm", (DL_FUNC) &_StableMCD_bootstrap_glm, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_StableMCD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
