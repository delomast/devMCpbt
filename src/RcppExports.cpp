// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// MCpbt
Rcpp::List MCpbt(int iter, int burnIn, int thin, unsigned int seed, Rcpp::NumericVector piTotPrior, Rcpp::NumericVector ohnc, Rcpp::NumericVector piTotInitial, Rcpp::NumericVector oUTInitial, Rcpp::NumericVector groups, int nPBT, Rcpp::NumericVector GSI_values, Rcpp::NumericVector gsiUT, Rcpp::NumericMatrix pi_gsiInitial, Rcpp::NumericMatrix prior_pi_gsi, Rcpp::NumericMatrix ohnc_gsi, Rcpp::List values, Rcpp::List pi_VInitial, Rcpp::List pi_Vohnc, Rcpp::List pi_Vprior, Rcpp::NumericMatrix v_ut, Rcpp::NumericVector initZ, Rcpp::NumericVector t, Rcpp::List valuesOth, Rcpp::List pi_VInitialOth, Rcpp::List pi_VohncOth, Rcpp::List pi_VpriorOth, Rcpp::NumericMatrix v_utOth);
RcppExport SEXP _devMCpbt_MCpbt(SEXP iterSEXP, SEXP burnInSEXP, SEXP thinSEXP, SEXP seedSEXP, SEXP piTotPriorSEXP, SEXP ohncSEXP, SEXP piTotInitialSEXP, SEXP oUTInitialSEXP, SEXP groupsSEXP, SEXP nPBTSEXP, SEXP GSI_valuesSEXP, SEXP gsiUTSEXP, SEXP pi_gsiInitialSEXP, SEXP prior_pi_gsiSEXP, SEXP ohnc_gsiSEXP, SEXP valuesSEXP, SEXP pi_VInitialSEXP, SEXP pi_VohncSEXP, SEXP pi_VpriorSEXP, SEXP v_utSEXP, SEXP initZSEXP, SEXP tSEXP, SEXP valuesOthSEXP, SEXP pi_VInitialOthSEXP, SEXP pi_VohncOthSEXP, SEXP pi_VpriorOthSEXP, SEXP v_utOthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type burnIn(burnInSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type piTotPrior(piTotPriorSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ohnc(ohncSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type piTotInitial(piTotInitialSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type oUTInitial(oUTInitialSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< int >::type nPBT(nPBTSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type GSI_values(GSI_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gsiUT(gsiUTSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pi_gsiInitial(pi_gsiInitialSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type prior_pi_gsi(prior_pi_gsiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ohnc_gsi(ohnc_gsiSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pi_VInitial(pi_VInitialSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pi_Vohnc(pi_VohncSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pi_Vprior(pi_VpriorSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type v_ut(v_utSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type initZ(initZSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type valuesOth(valuesOthSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pi_VInitialOth(pi_VInitialOthSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pi_VohncOth(pi_VohncOthSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pi_VpriorOth(pi_VpriorOthSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type v_utOth(v_utOthSEXP);
    rcpp_result_gen = Rcpp::wrap(MCpbt(iter, burnIn, thin, seed, piTotPrior, ohnc, piTotInitial, oUTInitial, groups, nPBT, GSI_values, gsiUT, pi_gsiInitial, prior_pi_gsi, ohnc_gsi, values, pi_VInitial, pi_Vohnc, pi_Vprior, v_ut, initZ, t, valuesOth, pi_VInitialOth, pi_VohncOth, pi_VpriorOth, v_utOth));
    return rcpp_result_gen;
END_RCPP
}
// propClipped
Rcpp::NumericVector propClipped(int Nclip, int Nunclip, int NumResults, unsigned int seed, Rcpp::NumericVector clipPrior);
RcppExport SEXP _devMCpbt_propClipped(SEXP NclipSEXP, SEXP NunclipSEXP, SEXP NumResultsSEXP, SEXP seedSEXP, SEXP clipPriorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Nclip(NclipSEXP);
    Rcpp::traits::input_parameter< int >::type Nunclip(NunclipSEXP);
    Rcpp::traits::input_parameter< int >::type NumResults(NumResultsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type clipPrior(clipPriorSEXP);
    rcpp_result_gen = Rcpp::wrap(propClipped(Nclip, Nunclip, NumResults, seed, clipPrior));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_devMCpbt_MCpbt", (DL_FUNC) &_devMCpbt_MCpbt, 27},
    {"_devMCpbt_propClipped", (DL_FUNC) &_devMCpbt_propClipped, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_devMCpbt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
