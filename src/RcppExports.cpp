// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// MCpbt
Rcpp::NumericVector MCpbt(int iter, int burnIn, int thin, unsigned int seed, int Nclip, int Nunclip, Rcpp::NumericVector clipPrior, bool clippedBool, Rcpp::NumericVector piTotPrior, Rcpp::NumericVector ohnc, Rcpp::NumericVector piTotInitial, Rcpp::NumericVector oUTInitial, Rcpp::NumericVector groups, Rcpp::List values, Rcpp::List pi_VInitial, Rcpp::List pi_Vohnc, Rcpp::List pi_Vprior, Rcpp::NumericMatrix v_ut, Rcpp::NumericVector initZ, Rcpp::NumericVector t);
RcppExport SEXP _devMCpbt_MCpbt(SEXP iterSEXP, SEXP burnInSEXP, SEXP thinSEXP, SEXP seedSEXP, SEXP NclipSEXP, SEXP NunclipSEXP, SEXP clipPriorSEXP, SEXP clippedBoolSEXP, SEXP piTotPriorSEXP, SEXP ohncSEXP, SEXP piTotInitialSEXP, SEXP oUTInitialSEXP, SEXP groupsSEXP, SEXP valuesSEXP, SEXP pi_VInitialSEXP, SEXP pi_VohncSEXP, SEXP pi_VpriorSEXP, SEXP v_utSEXP, SEXP initZSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type burnIn(burnInSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type Nclip(NclipSEXP);
    Rcpp::traits::input_parameter< int >::type Nunclip(NunclipSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type clipPrior(clipPriorSEXP);
    Rcpp::traits::input_parameter< bool >::type clippedBool(clippedBoolSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type piTotPrior(piTotPriorSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ohnc(ohncSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type piTotInitial(piTotInitialSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type oUTInitial(oUTInitialSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pi_VInitial(pi_VInitialSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pi_Vohnc(pi_VohncSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pi_Vprior(pi_VpriorSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type v_ut(v_utSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type initZ(initZSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(MCpbt(iter, burnIn, thin, seed, Nclip, Nunclip, clipPrior, clippedBool, piTotPrior, ohnc, piTotInitial, oUTInitial, groups, values, pi_VInitial, pi_Vohnc, pi_Vprior, v_ut, initZ, t));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_devMCpbt_MCpbt", (DL_FUNC) &_devMCpbt_MCpbt, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_devMCpbt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
