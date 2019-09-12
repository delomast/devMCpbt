#' This function prepares the data to run multiple strata
#' 
#' The model is actually run by a separate function in order to allow
#' a user to edit the priors and initial values before running.
#' 
#' 

prepMultStrata <- function(trapData, tags, GSIcol, PBTcol, strataCol, variableCols = c(), variableColsOth = c(), adFinCol, AI = TRUE, 
									 GSIgroups = NA,
									 variableValues = NA, variableValuesOth = NA, verbose = TRUE){
	
	#turn adFinCol into boolean if necessary
	if(!is.logical(trapData[,adFinCol])){
		nonValid <- sum(!is.na(trapData[,adFinCol]) & !(trapData[,adFinCol] %in% c("AD", "AI")))
		if(nonValid > 0){
			errMessage <- paste(nonValid, "observations that are not valid options for", adFinCol,
				"\nthe adFinCol must either be a logical variable, with TRUE for ad-intact,", 
				"or be a character variable with values of AD and AI for ad-clipped and ad-intact, respectively.", 
				"\n Missing data should have values of NA.")
			stop(errMessage)
		}
		trapData[,adFinCol] <- trapData[,adFinCol] == "AI"
		trapData[,adFinCol] <- as.logical(trapData[,adFinCol])
	}
	
	### first filter data to either Ad-intact or Ad-clipped
	
	if(AI) {
		trapData <- trapData[!is.na(trapData[,adFinCol]) & trapData[,adFinCol] & !is.na(trapData[,strataCol]),]
	} else if (!AI){
		trapData <- trapData[!is.na(trapData[,adFinCol]) & !trapData[,adFinCol] & !is.na(trapData[,strataCol]),]
	} else {
		stop("Non-boolean input for AI in prepOneStrataAI")
	}

	numObs <- nrow(trapData)
	if(verbose & AI) cat("\nFound", numObs, "ad-intact observations.\n")
	if(verbose & !AI) cat("\nFound", numObs, "ad-clipped observations.\n")
	trapData <- trapData[!is.na(trapData[,PBTcol]) & !is.na(trapData[,GSIcol]),]
	if(verbose) cat("\nDiscarding", numObs - nrow(trapData), "observations that were not", "attempted to be PBT and GSI assigned\n")
	
	#use all GSI groups present in the data if not specified
	if(is.na(GSIgroups)) GSIgroups <- sort(unique(trapData[,GSIcol])) # these are the different GSI groups in the dataset

	#use all variable values present in the data if not specified
	if (is.na(variableValues)){
		values <- list() # this is list of the categories in each variable
		for(v in names_variables){
			values[[v]] <- sort(unique(trapData[,v])) # these are the categories in var
		}
	} else {
		values <- variableValues
		names(values) <- names_variables
	}
	
	#use all "other" variable values present in the data if not specified
	if (is.na(variableValuesOth)){
		valuesOth <- list() # this is list of the categories in each variable
		for(v in names_variablesOth){
			valuesOth[[v]] <- sort(unique(trapData[,v])) # these are the categories in var
		}
	} else {
		valuesOth <- variableValuesOth
		names(valuesOth) <- names_variablesOth
	}
	
	# now prep each strata
	allStrata <- list() # list containing inputs for each strata
	for(s in sort(unique(trapData[,strataCol]))){
		strataData <- trapData[trapData[,strataCol] == s,] #select one strata
		allStrata[[s]] <- prepOneStrata(trapData = strataData, tags = tags, GSIcol = GSIcol, PBTcol = PBTcol, 
							variableCols = variableCols, variableColsOth = variableColsOth, adFinCol = adFinCol, AI = AI, 
							 verbose = FALSE, GSIgroups = GSIgroups,
							 variableValues = values, variableValuesOth = valuesOth)
	}
	
	return(allStrata)
}