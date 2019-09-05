// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// MCMCclip
Rcpp::NumericVector MCMCclip(int Nclip, int Ntotal, Rcpp::NumericVector clipPrior, int iter, int burnIn, int thin, unsigned int seed);
RcppExport SEXP _devMCpbt_MCMCclip(SEXP NclipSEXP, SEXP NtotalSEXP, SEXP clipPriorSEXP, SEXP iterSEXP, SEXP burnInSEXP, SEXP thinSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Nclip(NclipSEXP);
    Rcpp::traits::input_parameter< int >::type Ntotal(NtotalSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type clipPrior(clipPriorSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type burnIn(burnInSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(MCMCclip(Nclip, Ntotal, clipPrior, iter, burnIn, thin, seed));
    return rcpp_result_gen;
END_RCPP
}
// MCpbt
Rcpp::NumericVector MCpbt(int Nclip, int Ntotal, Rcpp::NumericVector clipPrior, int iter, int burnIn, int thin, unsigned int seed, bool clippedBool);
RcppExport SEXP _devMCpbt_MCpbt(SEXP NclipSEXP, SEXP NtotalSEXP, SEXP clipPriorSEXP, SEXP iterSEXP, SEXP burnInSEXP, SEXP thinSEXP, SEXP seedSEXP, SEXP clippedBoolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Nclip(NclipSEXP);
    Rcpp::traits::input_parameter< int >::type Ntotal(NtotalSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type clipPrior(clipPriorSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type burnIn(burnInSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type clippedBool(clippedBoolSEXP);
    rcpp_result_gen = Rcpp::wrap(MCpbt(Nclip, Ntotal, clipPrior, iter, burnIn, thin, seed, clippedBool));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_devMCpbt_MCMCclip", (DL_FUNC) &_devMCpbt_MCMCclip, 7},
    {"_devMCpbt_MCpbt", (DL_FUNC) &_devMCpbt_MCpbt, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_devMCpbt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}