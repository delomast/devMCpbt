#' This function takes data in a standard table format and prepares it for one strata for MCpbt
#' 
#' @param trapData dataframe with data for fish sampled from one strata - trap data for dam escapement
#' @param tags dataframe with tag rates for the groups to consider as present in that strata. This should
#'   NOT include the Unassigned group
#' @param GSIcol column name of column containing GSI assignments. If NA, the function gives all fish the same gsi assignment,
#'   which is equivalent to not having GSI information.
#' @param PBTcol column name of column containing PBT assignments
#' @param variableCols vector of column names of columns containing the variables to estimate composition and 
#'   to inform the group selection for PBT untagged fish
#' @param variableColsOth vector of column names of columns containing the variables to estimate composition but 
#'   NOT inform the group selection for PBT untagged fish
#' @param adFinCol column name of column containing adipose fin status - TRUE being intact FALSE being clipped, NA missing
#' @param physTagCol column name of column containing physTag status
#' @param AI TRUE to analyse ad-intact fish, FALSE to analyze ad-clipped fish
#' @param verbose TRUE to print some messages, FALSE to not
#' @param GSIgroups These are the values for all the GSI groups that you expect to be present in the population. If NA, then
#'   the values are taken to be all the unique values in the GSI column in the dataset.
#' @param variableValues a list, in the order of \code{variableCols} with entries having the values expected for each variable. This
#'   is helpful to make sure the same variable values are estiamted in each strata even if one value is not observed in all strata. If
#'   NA, it uses the values present in the dataset.
#' @param variableValuesOth Same as variableValues, but for variableColsOth
#' 

prepOneStrataAI <- function(trapData, tags, GSIcol, PBTcol, variableCols = c(), variableColsOth = c(), adFinCol, AI = TRUE, 
									 verbose = TRUE, GSIgroups = NA,
									 variableValues = NA, variableValuesOth = NA){
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
		trapData <- trapData[!is.na(trapData[,adFinCol]) & trapData[,adFinCol],]
	} else if (!AI){
		trapData <- trapData[!is.na(trapData[,adFinCol]) & !trapData[,adFinCol],]
	} else {
		stop("Non-boolean input for AI in prepOneStrataAI")
	}
	numObs <- nrow(trapData)
	if(verbose & AI) cat("\nFound", numObs, "ad-intact observations.\n")
	if(verbose & !AI) cat("\nFound", numObs, "ad-clipped observations.\n")
	trapData <- trapData[!is.na(trapData[,PBTcol]) & !is.na(trapData[,GSIcol]),]
	if(verbose) cat("\nDiscarding", numObs - nrow(trapData), "observations that were not attempted to be PBT and GSI assigned.\n")
	
	
	### now calculate input values
	if(is.na(GSIgroups)) GSIgroups <- sort(unique(trapData[,GSIcol])) # these are the different GSI groups in the dataset
	nGSI <- length(GSIgroups) #number of wild GSI groups
	nPBT <- nrow(tags) #number of PBT groups

	groups <- c(tags[,1], GSIgroups) #this is list of all the PBT groups and the wild
	ohnc <- rep(0, nrow(tags)) # count how many assigned to each pbt group
	for(i in 1:nrow(tags)){
		ohnc[i] <- sum(trapData[,PBTcol] == tags[i,1])
	}
	ohnc <- c(ohnc, rep(0, nGSI)) # add wild with 0
	t <- as.numeric(c(tags[,2], rep(0, nGSI))) # tag rates, with wild as 0

	#GSI observations for untagged fish
	untagBool <- trapData[,PBTcol] == "Unassigned" # boolean variable to select PBT-untagged fish
	gsiUT <- trapData[untagBool, GSIcol]
	
	#GSI observations for PBT-tagged fish
	ohnc_GSI <- matrix(0, nrow = nrow(tags), ncol = nGSI)
	for(i in 1:nrow(tags)){
		for(j in 1:length(GSIgroups)){
			ohnc_GSI[i,j] <- sum(trapData[,PBTcol] == tags[i,1] & trapData[,GSIcol] == GSIgroups[j])
		}
	}
	
	### set up categorical variables
	names_variables <- variableCols
	v_ut <- list() # this is list of the observations for categorical variables for untagged fish
	for(v in names_variables){
		v_ut[[v]] <- trapData[untagBool, v] # these are the values for Var for untagged fish
	}
	if (is.na(variableValues)){
		values <- list() # this is list of the categories in each variable
		for(v in names_variables){
			values[[v]] <- sort(unique(trapData[,v])) # these are the categories in var
		}
	} else {
		values <- variableValues
		names(values) <- names_variables
	}
	
	V_ohnc <- list() # this is list of matrices of counts of values of each variable in the ohnc
	for(v in names_variables){
		V1_ohnc <- matrix(0, nrow = length(groups), ncol = length(values[[v]])) # this is a matrix of counts of each Var1 value in the ohnc
		for(i in 1:nPBT){
			tempData <- trapData[trapData[,PBTcol] == groups[i],]
			for(j in 1:length(values[[v]])){
					V1_ohnc[i,j] <- sum(tempData[,v] == values[[v]][j])
			}
		}
		V_ohnc[[v]] <- V1_ohnc
	}
	
	### set up "other" categorical variables
	names_variablesOth <- variableColsOth
	v_utOth <- list() # this is list of the observations for categorical variables for untagged fish
	for(v in names_variablesOth){
		v_utOth[[v]] <- trapData[untagBool, v] # these are the values for Var for untagged fish
	}
	if (is.na(variableValuesOth)){
		valuesOth <- list() # this is list of the categories in each variable
		for(v in names_variablesOth){
			valuesOth[[v]] <- sort(unique(trapData[,v])) # these are the categories in var
		}
	} else {
		valuesOth <- variableValuesOth
		names(valuesOth) <- names_variablesOth
	}
	
	V_ohncOth <- list() # this is list of matrices of counts of values of each variable in the ohnc
	for(v in names_variablesOth){
		V1_ohnc <- matrix(0, nrow = length(groups), ncol = length(valuesOth[[v]])) # this is a matrix of counts of each Var1 value in the ohnc
		for(i in 1:nPBT){
			tempData <- trapData[trapData[,PBTcol] == groups[i],]
			for(j in 1:length(valuesOth[[v]])){
					V1_ohnc[i,j] <- sum(tempData[,v] == valuesOth[[v]][j])
			}
		}
		V_ohncOth[[v]] <- V1_ohnc
	}
	
	### set up priors - these are default, user can modify them as they see fit before running the model
	# these are alphas for Dirichlet priors
	prior_piTot <- rep(1, length(groups))
	#for categorical variables
	prior_piV <- list()
	for(v in names_variables){
		prior_piV[[v]] <- matrix(1, nrow = length(groups), ncol = length(values[[v]])) # here just using a uniform for all
	}
	#for "other" categorical variables
	prior_piVOth <- list()
	for(v in names_variablesOth){
		prior_piVOth[[v]] <- matrix(1, nrow = length(groups), ncol = length(valuesOth[[v]])) # here just using a uniform for all
	}
	
	#prior for GSI composition of PBT groups
	prior_piGSI <- matrix(1, nrow = nrow(tags), ncol = nGSI)
	
	
	### initial values - these are default, user can modify them as they see fit before running the model
	piTot <- rep(1/length(ohnc), length(ohnc)) #all equal
	z <- sample(groups, sum(trapData[,PBTcol] == "Unassigned"), replace = TRUE) #all random
	oUT <- rep(0, length(groups))
	for(i in 1:length(groups)){
		oUT[i] <- sum(z == groups[i])
	}
	
	pi_V <- list()
	for(v in names_variables){
		pi_V[[v]] <- matrix(1/length(values[[v]]), nrow = length(groups), ncol = length(values[[v]])) #all equal
	}
	
	pi_VOth<- list()
	for(v in names_variablesOth){
		pi_VOth[[v]] <- matrix(1/length(valuesOth[[v]]), nrow = length(groups), ncol = length(valuesOth[[v]])) #all equal
	}
	
	#proportions of groups that GSI assign to each GSI category
	pi_gsi <- matrix(0, nrow = length(groups), ncol = nGSI)
	for(i in 1:nrow(tags)){
		pi_gsi[i,] <- 1/nGSI #all equal
	}
	currentCol <- 1
	for(i in (nrow(tags)+1):length(groups)){ #setting GSI groups to be exactly as observed
		pi_gsi[i,currentCol] <- 1 #others were initialized at 0, so only set this one
		currentCol <- currentCol + 1
	}
	
	#convert some variables to ints and matrices instead of characters and lists, respectively
	groupsInt <- 1:length(groups) #groups as ints instead of characters
	# make key to convert between ints and original values
	groupsKey <- data.frame(groupName = groups,
									groupID = groupsInt,
									stringsAsFactors = FALSE)
	zInt <- z # change z's accordingly
	for(i in groupsInt){
		zInt[z == groups[i]] <- i
	}
	zInt <- as.numeric(zInt) #turn types into ints
	
	# GSI groups
	#make key
	GSIkey <- data.frame(GSIvalue = GSIgroups,
									GSIid = 1:nGSI,
									stringsAsFactors = FALSE)
	GSIgroups <- 1:nGSI
	# gsiUT <- trapData[untagBool, GSIcol]
	old <- gsiUT
	for(i in 1:nGSI){
		gsiUT[old == GSIkey[i,1]] <- GSIkey[i,2]
	}
	rm(old)
	
	
	# categorical variables
	variKey <- list()
	for(v in names_variables){
		# make key to convert between ints and original values
		variKey[[v]] <- data.frame(category = values[[v]],
									categoryID = 1:length(values[[v]]),
									stringsAsFactors = FALSE)
		#convert values to ints in input
		values[[v]] <- 1:length(values[[v]])
		v_ut[[v]][is.na(v_ut[[v]])] <- -9
		old <- v_ut[[v]]
		for(i in 1:nrow(variKey[[v]])){
			v_ut[[v]][old == variKey[[v]][i,1]] <- variKey[[v]][i,2]
		}
		v_ut[[v]] <- as.numeric(v_ut[[v]])
		rm(old)
	}

	if(length(v_ut) == 0){
		v_utMat <- matrix(0,0,0)
	} else {
		v_utMat <- matrix(nrow = length(v_ut[[1]]), ncol = 0)
		for(i in v_ut){
			v_utMat <- cbind(v_utMat, i)
		}
		colnames(v_utMat) <- names_variables
	}
	
	# "other" categorical variables
	variKeyOth <- list()
	for(v in names_variablesOth){
		# make key to convert between ints and original values
		variKeyOth[[v]] <- data.frame(category = valuesOth[[v]],
									categoryID = 1:length(valuesOth[[v]]),
									stringsAsFactors = FALSE)
		#convert values to ints in input
		valuesOth[[v]] <- 1:length(valuesOth[[v]])
		v_utOth[[v]][is.na(v_utOth[[v]])] <- -9
		old <- v_utOth[[v]]
		for(i in 1:nrow(variKeyOth[[v]])){
			v_utOth[[v]][old == variKeyOth[[v]][i,1]] <- variKeyOth[[v]][i,2]
		}
		v_utOth[[v]] <- as.numeric(v_utOth[[v]])
		rm(old)
	}
	
	if(length(v_utOth) == 0){
		v_utMatOth <- matrix(0,0,0)
	} else {
		v_utMatOth <- matrix(nrow = length(v_utOth[[1]]), ncol = 0)
		for(i in v_utOth){
			v_utMatOth <- cbind(v_utMatOth, i)
		}
		colnames(v_utMatOth) <- names_variablesOth
	}
	

	# now package all values and return
	return(list( #first all input values to MCpbt in order
		piTotPrior = prior_piTot, ohnc = ohnc, piTotInitial = piTot, oUTInitial = oUT, groups = groupsInt, 
		nPBT = nPBT, GSI_values = GSIgroups, gsiUT = gsiUT,
		pi_gsiInitial = pi_gsi, prior_pi_gsi = prior_piGSI, ohnc_gsi = ohnc_GSI, 
		values = values, pi_VInitial = pi_V, pi_Vohnc = V_ohnc, pi_Vprior = prior_piV,
		v_ut = v_utMat, initZ = zInt, t = t, valuesOth = valuesOth, pi_VInitialOth = pi_VOth, 
		pi_VohncOth = V_ohncOth, pi_VpriorOth = prior_piVOth,
		v_utOth = v_utMatOth,
		#then all keys
		groupsKey = groupsKey, GSIkey = GSIkey, variKey = variKey, variKeyOth = variKeyOth
	))
}