#' This function uses the input created by \code{prepMultStrata} to estimate the
#' proportion of fish that are clipped for each strata.
#' 
#' A wrapper for \code{MCpbt}
#' 
#' 

estimStrataMCpbt <- function(prepData, iter, burnIn, thin, seed = NA){
	
	# get seed
	if(is.na(seed)){
		seed <- ceiling(as.numeric(format(Sys.time(), "%S")) * as.numeric(format(Sys.time(), "%j")) * as.numeric(format(Sys.time(), "%M")))
	}
	
	
	# get samples
	strataEstimates <- list()
	for(i in 1:length(prepData)){
		cStrata <- prepData[[i]]
		strataEstimates[[i]] <- MCpbt(iter = iter, burnIn = burnIn , thin = thin, seed = seed, #overall parameters for the chain
					     piTotPrior = cStrata$piTotPrior, ohnc = cStrata$ohnc, piTotInitial = cStrata$piTotInitial, #piTotal parameters
					     oUTInitial = cStrata$oUTInitial, groups = cStrata$groups,
					     nPBT = cStrata$nPBT, GSI_values = cStrata$GSI_values, gsiUT = cStrata$gsiUT, #pi_gsi parameters
					     pi_gsiInitial = cStrata$pi_gsiInitial, prior_pi_gsi = cStrata$prior_pi_gsi,
					     ohnc_gsi = cStrata$ohnc_gsi,
					     values = cStrata$values, pi_VInitial = cStrata$pi_VInitial, pi_Vohnc = cStrata$pi_Vohnc, pi_Vprior = cStrata$pi_Vprior, #pi_V parameters
					     v_ut = cStrata$v_ut, 
					     initZ = cStrata$initZ, t = cStrata$t, #z parameters
					     valuesOth = cStrata$valuesOth, pi_VInitialOth = cStrata$pi_VInitialOth, pi_VohncOth = cStrata$pi_VohncOth, #pi_VOth parameters
					     pi_VpriorOth = cStrata$pi_VpriorOth, 
                 v_utOth = cStrata$v_utOth
					     )
	
	}
	
	return(strataEstimates)
	
}
