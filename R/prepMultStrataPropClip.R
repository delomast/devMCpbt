#' This function prepares the data to run multiple strata for estiamting the proportion clipped
#' 
prepMultStrataPropClip <- function(trapData, strataCol, adFinCol, verbose = TRUE){
	
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
	
	allStrata <- list() # list containing inputs for each strata
	trapData <- trapData[!is.na(trapData[,adFinCol]) & !is.na(trapData[,strataCol]),]
	if(verbose) cat("Using", nrow(trapData), "observations with adFin status and an assigned strata")
	for(s in sort(unique(trapData[,strataCol]))){
		strataData <- trapData[trapData[,strataCol] == s,] #select one strata
		allStrata[[s]] <- c(sum(!trapData[,adFinCol], trapData[,adFinCol])) # clipped, unclipped
	}
}
