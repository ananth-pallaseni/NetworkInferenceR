// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// infer_mi_network
DataFrame infer_mi_network(const DataFrame& df);
RcppExport SEXP _NetworkInferenceR_infer_mi_network(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const DataFrame& >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(infer_mi_network(df));
    return rcpp_result_gen;
END_RCPP
}
// infer_puc_network
DataFrame infer_puc_network(const DataFrame& df);
RcppExport SEXP _NetworkInferenceR_infer_puc_network(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const DataFrame& >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(infer_puc_network(df));
    return rcpp_result_gen;
END_RCPP
}
// infer_pid_network
DataFrame infer_pid_network(const DataFrame& df);
RcppExport SEXP _NetworkInferenceR_infer_pid_network(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const DataFrame& >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(infer_pid_network(df));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NetworkInferenceR_infer_mi_network", (DL_FUNC) &_NetworkInferenceR_infer_mi_network, 1},
    {"_NetworkInferenceR_infer_puc_network", (DL_FUNC) &_NetworkInferenceR_infer_puc_network, 1},
    {"_NetworkInferenceR_infer_pid_network", (DL_FUNC) &_NetworkInferenceR_infer_pid_network, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_NetworkInferenceR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
