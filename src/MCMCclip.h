// This is the header file for MCMCclip function
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#ifndef MCMCCLIP_H
#define MCMCCLIP_H

#include <Rcpp.h>
#include <random>
Rcpp::NumericVector MCMCclip(int Nclip, int Nunclip, Rcpp::NumericVector clipPrior, int NumResults, std::mt19937 * rNum);

#endif
