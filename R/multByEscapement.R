#' Multiplies escapement estimates by composition estimates
#' 
#' 


multByEscapement <- function(	input, clipInput, aiRes, clipRes, escapementByStrata, verbose = TRUE){
	#### Need to:
	# make clip and unclip optional
	# have optional summary of posterior means, medians, and CI's written to file
	
	
	##  note that clipInput and clipRes are optional
	
	## check inputs
	numStrata <- length(input)
	## strata names the same and in order
	## escapement same length as others
	
	## identify AI
	AI <- input[[1]]$AI
	if(verbose && AI) cat("\nDecomposing ad-intact fish")
	if(verbose && !AI) cat("\nDecomposing ad-clipped fish")
	
	strataEstimates <- list()
	#for each strata
	for(i in 1:length(input)){
		strataName <- input[[i]]$strataName
		#clipped and unclipped breakdown
		clipUnclip <- escapementByStrata[[i]]*clipRes[[i]] ## calculate number clipped
		clipUnclip <- cbind(clipUnclip, escapementByStrata[[i]] - clipUnclip) ##remainder are unclipped
		colnames(clipUnclip) <- c("Number_clipped", "Number_unclipped")
		#select clipped or unclipped total to use downstream
		if (AI) {
			total <- clipUnclip[,2]
		} else{
			total <- clipUnclip[,1] 
		}
		
		#save current strata to make access easier
		cStrataInput <- input[[i]]
		cStrata <- aiRes[[i]]
		
		# breakdown by piTot
		piTotNumbers <- cStrata$piTot * total
		colnames(piTotNumbers) <- cStrataInput$groupsKey[match(cStrataInput$groups, cStrataInput$groupsKey[,2]),1]
		
		# piGSI - this doesn't seem that useful, but why not calculate it all the same
		cur <- 1 # to iterate appropriately across columns
		end <- length(cStrataInput$GSI_values)
		GSInumbers <- matrix(0, nrow = nrow(cStrata$piGSI), ncol = ncol(cStrata$piGSI))
		colnames(GSInumbers) <- paste0("temp", 1:ncol(cStrata$piGSI))
		for(j in 1:ncol(piTotNumbers)){
			GSInumbers[,cur:end] <- cStrata$piGSI[,cur:end]*piTotNumbers[,j]
			# assign column names
			colnames(GSInumbers)[cur:end] <- paste0(cStrataInput$groupsKey[cStrataInput$groupsKey[,2] == cStrataInput$groups[j],1],
																 "_",
																 cStrataInput$GSIkey[match(cStrataInput$GSI_values, cStrataInput$GSIkey[,2]),1])
			cur <- end + 1
			end <- end + length(cStrataInput$GSI_values)
		}
		rm(cur)
		rm(end)
		
		# piV
		piVnumbers <- list()
		pos <- 1
		for(v in names(cStrataInput$values)){
			props <- cStrata$piV[[pos]]
			cur <- 1 # to iterate appropriately across columns
			end <- length(cStrataInput$values[[pos]])
			piVcounts <- matrix(0, nrow = nrow(props), ncol = ncol(props))
			colnames(piVcounts) <- paste0("temp", 1:ncol(props))
			for(j in 1:ncol(piTotNumbers)){
				piVcounts[,cur:end] <- props[,cur:end]*piTotNumbers[,j]
				# assign column names
				colnames(piVcounts)[cur:end] <- paste0(cStrataInput$groupsKey[cStrataInput$groupsKey[,2] == cStrataInput$groups[j],1],
																	 "_",
																	 cStrataInput$variKey[[pos]][match(cStrataInput$values[[pos]], cStrataInput$variKey[[pos]][,2]),1])
				cur <- end + 1
				end <- end + length(cStrataInput$values[[pos]])
			}
			rm(cur)
			rm(end)
			piVnumbers[[v]] <- piVcounts
			pos <- pos + 1
		}
		rm(pos)
		
		# piVOth
		piVnumbersOth <- list()
		pos <- 1
		for(v in names(cStrataInput$valuesOth)){
			props <- cStrata$piVOth[[pos]]
			cur <- 1 # to iterate appropriately across columns
			end <- length(cStrataInput$valuesOth[[pos]])
			piVcounts <- matrix(0, nrow = nrow(props), ncol = ncol(props))
			colnames(piVcounts) <- paste0("temp", 1:ncol(props))
			for(j in 1:ncol(piTotNumbers)){
				piVcounts[,cur:end] <- props[,cur:end]*piTotNumbers[,j]
				# assign column names
				colnames(piVcounts)[cur:end] <- paste0(cStrataInput$groupsKey[cStrataInput$groupsKey[,2] == cStrataInput$groups[j],1],
																	 "_",
																	 cStrataInput$variKeyOth[[pos]][match(cStrataInput$valuesOth[[pos]], cStrataInput$variKeyOth[[pos]][,2]),1])
				cur <- end + 1
				end <- end + length(cStrataInput$valuesOth[[pos]])
			}
			rm(cur)
			rm(end)
			piVnumbersOth[[v]] <- piVcounts
			pos <- pos + 1
		}
		rm(pos)
		
		strataEstimates[[strataName]] <- list(clipEstim = clipUnclip, 
														  piTotEstim = piTotNumbers,
														  GSIestim = GSInumbers,
														  variableEstim = piVnumbers,
														  variableEstimOth = piVnumbersOth,
														  strataName = strataName
														  )
	}
	
	# now sum all strata together
	
	
		# clipEstim = clipUnclip, 
	totClipUnclip <- matrix(0, nrow = nrow(strataEstimates[[1]]$clipEstim), ncol = 2)
	for(i in 1:numStrata){
		totClipUnclip[,1] <- totClipUnclip[,1] + strataEstimates[[i]]$clipEstim[,1]
		totClipUnclip[,2] <- totClipUnclip[,2] + strataEstimates[[i]]$clipEstim[,2]
	}
	
	## general strategy: 
	## get all column names across all strata
	## build matrix
	## sum up by mathcing column names across and within strata
	##		the within option allows for changing the key after running
	##		to easily pool groups
	
		# piTotEstim = piTotNumbers,
	## get all column names across all strata
	columns <- c()
	for(i in 1:numStrata){
		columns <- c(columns, colnames(strataEstimates[[i]]$piTotEstim))
	}
	columns <- unique(columns)
	## build matrix
	totPiTotEstim <- matrix(0, nrow = nrow(strataEstimates[[1]]$piTotEstim), ncol = length(columns))
	colnames(totPiTotEstim) <- columns
	## sum up by mathcing column names across and within strata
	for(c in columns){
		for(i in 1:numStrata){
			tempEstim <- strataEstimates[[i]]$piTotEstim
			for(j in which(colnames(tempEstim) == c)){
				totPiTotEstim[,c] <- totPiTotEstim[,c] + tempEstim[,j]
			}
		}
	}
	rm(columns)
	
		# GSIestim = GSInumbers,
	columns <- c()
	for(i in 1:numStrata){
		columns <- c(columns, colnames(strataEstimates[[i]]$GSIestim))
	}
	columns <- unique(columns)
	totGSIestim <- matrix(0, nrow = nrow(strataEstimates[[1]]$GSIestim), ncol = length(columns))
	colnames(totGSIestim) <- columns
	for(c in columns){
		for(i in 1:numStrata){
			tempEstim <- strataEstimates[[i]]$GSIestim
			for(j in which(colnames(tempEstim) == c)){
				totGSIestim[,c] <- totGSIestim[,c] + tempEstim[,j]
			}
		}
	}
	rm(columns)
	
		# variableEstim = piVnumbers,
	totpiVnumbers <- list()
	for(v in names(input[[1]]$values)){
		columns <- c()
		for(i in 1:numStrata){
			columns <- c(columns, colnames(strataEstimates[[i]]$variableEstim[[v]]))
		}
		columns <- unique(columns)
		totVarestim <- matrix(0, nrow = nrow(strataEstimates[[1]]$variableEstim[[v]]), ncol = length(columns))
		colnames(totVarestim) <- columns
		for(c in columns){
			for(i in 1:numStrata){
				tempEstim <- strataEstimates[[i]]$variableEstim[[v]]
				for(j in which(colnames(tempEstim) == c)){
					totVarestim[,c] <- totVarestim[,c] + tempEstim[,j]
				}
			}
		}
		rm(columns)
		totpiVnumbers[[v]] <- totVarestim
	}
		
		# variableEstimOth = piVnumbersOth,
	totpiVnumbersOth <- list()
	for(v in names(input[[1]]$valuesOth)){
		columns <- c()
		for(i in 1:numStrata){
			columns <- c(columns, colnames(strataEstimates[[i]]$variableEstimOth[[v]]))
		}
		columns <- unique(columns)
		totVarestimOth <- matrix(0, nrow = nrow(strataEstimates[[1]]$variableEstimOth[[v]]), ncol = length(columns))
		colnames(totVarestimOth) <- columns
		for(c in columns){
			for(i in 1:numStrata){
				tempEstim <- strataEstimates[[i]]$variableEstimOth[[v]]
				for(j in which(colnames(tempEstim) == c)){
					totVarestimOth[,c] <- totVarestimOth[,c] + tempEstim[,j]
				}
			}
		}
		rm(columns)
		totpiVnumbersOth[[v]] <- totVarestimOth
	}
		# strataName = strataName
	strataNames <- c()
	for(i in 1:numStrata){
		strataNames <- c(strataNames, strataEstimates[[i]]$strataName)
	}
	
	
	# to return, per strata estimates, total estimates, and strata names
	
	return(list(totClipUnclip = totClipUnclip,
					totPiTotEstim = totPiTotEstim,
					totGSIestim = totGSIestim,
					totpiVnumbers = totpiVnumbers,
					totpiVnumbersOth = totpiVnumbersOth,
					strataNames = strataNames,
					strataEstimates = strataEstimates
					))
}
