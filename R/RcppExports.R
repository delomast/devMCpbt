# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

MCMCclip <- function(Nclip, Ntotal, clipPrior, iter, burnIn, thin, seed) {
    .Call(`_devMCpbt_MCMCclip`, Nclip, Ntotal, clipPrior, iter, burnIn, thin, seed)
}

MCpbt <- function(Nclip, Ntotal, clipPrior, iter, burnIn, thin, seed, clippedBool) {
    .Call(`_devMCpbt_MCpbt`, Nclip, Ntotal, clipPrior, iter, burnIn, thin, seed, clippedBool)
}
