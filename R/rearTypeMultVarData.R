#' H HNC W data simulation function
#'
#'
#' @param varArray list of arrays with first dimension for each PBT groups (in the same order
#'   as relSizePBTgroups) plus wild as the last entry of the first dimension
#'   second dimension for each strata, and third dimension for each category in the variable.
#'   Values are the proportion of the strata and PBT group that belong to each category of the variable

rearTypeMultVarData <- function(sampRates = .1, censusSizes = 3000, relSizePBTgroups = 1, tagRates = .8, obsTagRates = .8, physTagRates = 0,
				    true_clipped = .5, true_noclip_H = .2, true_wild = .3, varArray = NA){

	# check tag rates all between 0 and 1
	if(sum(tagRates > 1) > 0 || sum(tagRates < 0) > 0){
		stop("tagRates cannot be greater than 1 or less than 0")
	}
	if(sum(obsTagRates > 1) > 0 || sum(obsTagRates < 0) > 0){
		stop("obsTagRates cannot be greater than 1 or less than 0")
	}
	if(sum(physTagRates > 1) > 0 || sum(physTagRates < 0) > 0){
		stop("physTagRates cannot be greater than 1 or less than 0")
	}
	#sample rates check
	if(sum(sampRates > 1) > 0 || sum(sampRates < 0) > 0){
		stop("sampRates cannot be greater than 1 or less than 0")
	}
	# true props
	propCheck <- mapply(function(x,y,z){
			return(!isTRUE(all.equal(sum(x,y,z), 1))) #floating point arithmetic error avoidance
		}, true_clipped, true_noclip_H, true_wild)
	if(sum(propCheck) > 0){
		stop("true_clipped, true_noclip_H, and true_wild must sum to 1 for all strata")
	}

	# check varArray
	if(!is.list(varArray)){
		stop("No list of varArrays provided")
	}
	for(i in 1:length(varArray)){
		if(length(dim(varArray[[i]])) != 3){
			stop("wrong dimensions for one varArray provided")
		}
		if(dim(varArray[[i]])[1] != length(relSizePBTgroups) + 1 || dim(varArray[[i]])[2] != length(sampRates)){
			stop("Dimensions one and two of varArrays must be number of PBT groups + 1 and number of strata, respectively")
		}
	}


	#count number of strata
	nStrata <- length(sampRates)

	#empty dataframe for this simulation
	simData <- data.frame(NULL, stringsAsFactors = FALSE)
	#simulate data
	for(i in 1:nStrata){
		# choose number sampled (trapped)
		nSampled <- rbinom(1, censusSizes[i], sampRates[i])
		if(nSampled < 1){ #prevent strata with no trapped fish - if there were no trapped fish, the user woudl combine them with others
			nSampled <- 1
		}
		#choose numbers clipped, noclip_H, and wild
			## returns 1col matrix, turn into vector
		nRear <- as.vector(rmultinom(1, nSampled, c(true_clipped[i], true_noclip_H[i], true_wild[i])))
		#choose number of noclip_H that belong to each PBT group
			## returns 1col matrix, turn into vector
		PBTgroupNums <- as.vector(rmultinom(1, nRear[2], relSizePBTgroups))
		#determine number PBT tagged, physTagged, and untagged for each group
		#this mapply returns a matrix with columns as the PBT groups and rows
		#  as numBoth, numPBTonly, numPhysOnly, numUntag
		tagUntagNums <- mapply(function(x,y,z){
			# calculate probabilities of all four groups, assuming that being
			#   physically tagged and PBT tagged are independent
			#prob pbt and phys
			pBoth <- y*z
			#prob phys only
			pPhysOnly <- (1-y)*z
			#prob pbt only
			pPBTonly <- y*(1-z)
			#prob untagged
			pUntag <- 1 - pPBTonly - pPhysOnly - pBoth
			return(rmultinom(1, x, c(pBoth, pPBTonly, pPhysOnly, pUntag)))
			}, PBTgroupNums, tagRates, physTagRates)
		# build datasets
		pbtNames <- paste0("pbtGroup", 1:length(tagRates))

		#initiate dataset with clipped and wild fish
		#columns needed are clip, phystag, pbt group
		strataData <- data.frame(clip = c(rep("AD", nRear[1]), rep("AI", nRear[3])),
							physTag = c(rep("notag", (nRear[1] + nRear[3]))),
							pbtGroup = c(rep("Unassigned", (nRear[1] + nRear[3]))),
								stringsAsFactors = FALSE)
		# add variable columns
		for(j in 1:length(varArray)){
			strataData <- cbind(strataData,
							c(rep(NA, nRear[1]),
				    sample(1:dim(varArray[[j]])[3], nRear[3], replace = TRUE, prob = varArray[[j]][dim(varArray[[j]])[1],i,]))
						)
		}
		# name variable columns
		colnames(strataData)[(ncol(strataData) - length(varArray) + 1):ncol(strataData)] <- paste0("Var", 1:length(varArray))

		#now add various types of HNC fish
				#this mapply returns a matrix with columns as the PBT groups and rows
				#  as numBoth, numPBTonly, numPhysOnly, numUntag
				# tagUntagNums
		#first untagged
		for(c in which(tagUntagNums[4,] > 0)){
			numAdd <- tagUntagNums[4,c]
			tempData <- data.frame(clip = rep("AI", numAdd),
									physTag = rep("notag", numAdd),
									pbtGroup = rep("Unassigned", numAdd),
										stringsAsFactors = FALSE
							)
			# add variable columns
			for(j in 1:length(varArray)){
				tempData <- cbind(tempData,
								sample(1:dim(varArray[[j]])[3], numAdd, replace = TRUE, prob = varArray[[j]][c,i,])
							)
			}
			# name variable columns
			colnames(tempData) <- colnames(strataData)
			# add to data for the strata
			strataData <- rbind(strataData, tempData)
		}
		#Then PhysOnly
		for(c in which(tagUntagNums[3,] > 0)){
			numAdd <- tagUntagNums[3,c]
			tempData <- data.frame(clip = rep("AI", numAdd),
									physTag = rep("tag", numAdd),
									pbtGroup = rep("Unassigned", numAdd),
										stringsAsFactors = FALSE
							)
			# add variable columns
			for(j in 1:length(varArray)){
				tempData <- cbind(tempData,
								sample(1:dim(varArray[[j]])[3], numAdd, replace = TRUE, prob = varArray[[j]][c,i,])
							)
			}
			# name variable columns
			colnames(tempData) <- colnames(strataData)
			# add to data for the strata
			strataData <- rbind(strataData, tempData)
		}
		#pbt only
		#only loop through groups with one or more tagged fish present
		for(c in which(tagUntagNums[2,] > 0)){
			numAdd <- tagUntagNums[2,c]
			tempData <- data.frame(clip = rep("AI", numAdd),
								physTag = rep("notag", numAdd),
								pbtGroup = rep(pbtNames[c], numAdd), #pbtNames is same order as columns
									stringsAsFactors = FALSE)
			# add variable columns
			for(j in 1:length(varArray)){
				tempData <- cbind(tempData,
								sample(1:dim(varArray[[j]])[3], numAdd, replace = TRUE, prob = varArray[[j]][c,i,])
							)
			}
			# name variable columns
			colnames(tempData) <- colnames(strataData)
			# add to data for the strata
			strataData <- rbind(strataData, tempData)
		}
		#both
		#only loop through groups with one or more tagged fish present
		for(c in which(tagUntagNums[1,] > 0)){
			numAdd <- tagUntagNums[1,c]
			tempData <- data.frame(clip = rep("AI", numAdd),
								physTag = rep("tag", numAdd),
								pbtGroup = rep(pbtNames[c], numAdd), #pbtNames is same order as columns
									stringsAsFactors = FALSE)
			# add variable columns
			for(j in 1:length(varArray)){
				tempData <- cbind(tempData,
								sample(1:dim(varArray[[j]])[3], numAdd, replace = TRUE, prob = varArray[[j]][c,i,])
							)
			}
			# name variable columns
			colnames(tempData) <- colnames(strataData)
			# add to data for the strata
			strataData <- rbind(strataData, tempData)
		}
		#add simulated data for the strata to the main data frame for this simulations
		simData <- rbind(simData, cbind(i, strataData)) #columns are: strata, clip, physTag, pbtGroup
	}
	colnames(simData)[1:4] <- c("WeekNumber", "AdClip", "PhysTag", "GenParentHatchery")
	#generate window counts and tag rates for convenience
	windowData <- data.frame(week = 1:nStrata,
						count = censusSizes,
						collapse = 1:nStrata, stringsAsFactors = FALSE) #data is simulated by strata
	tagRateInput <- data.frame(pbtGroups = pbtNames, tagRates = obsTagRates, stringsAsFactors = FALSE) #give the function observed tag rates

	#output all three as a list
	return(list(trapData = simData, windowData = windowData, tagRates = tagRateInput))

}
