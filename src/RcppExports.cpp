// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// HierarCluster_paris
NumericMatrix HierarCluster_paris(Eigen::SparseMatrix<double> m);
RcppExport SEXP _HGC_core_HierarCluster_paris(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(HierarCluster_paris(m));
    return rcpp_result_gen;
END_RCPP
}
// HierarCluster_paris_time
NumericMatrix HierarCluster_paris_time(Eigen::SparseMatrix<double> m);
RcppExport SEXP _HGC_core_HierarCluster_paris_time(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(HierarCluster_paris_time(m));
    return rcpp_result_gen;
END_RCPP
}
// get_leaves
IntegerVector get_leaves(NumericMatrix hglink);
RcppExport SEXP _HGC_core_get_leaves(SEXP hglinkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type hglink(hglinkSEXP);
    rcpp_result_gen = Rcpp::wrap(get_leaves(hglink));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HGC_core_HierarCluster_paris", (DL_FUNC) &_HGC_core_HierarCluster_paris, 1},
    {"_HGC_core_HierarCluster_paris_time", (DL_FUNC) &_HGC_core_HierarCluster_paris_time, 1},
    {"_HGC_core_get_leaves", (DL_FUNC) &_HGC_core_get_leaves, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_HGC_core(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
