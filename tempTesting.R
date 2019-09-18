### testing script for devMCpbt
# not built as part of package

##########################################
### after changes to make GSI groups separate
###############

pbtGSImat <- matrix(c(.1, .8, .1,.8, .1, .1,.1, .1, .8), nrow = 3, ncol = 3, byrow = TRUE)
varMat <- list(
	matrix(c(rep(c(.1,.9), 3), rep(c(.9,.1), 3)), nrow = 6, ncol = 2, byrow = TRUE),
	matrix(c(rep(c(.4,.6), 3), rep(c(.6,.4), 3)), nrow = 6, ncol = 2, byrow = TRUE)
)

## two variables
testData <- generatePBTGSIdata(sampRate = .18, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = varMat)

## one variable
testData <- generatePBTGSIdata(sampRate = .18, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = varMat[1])

## no variables
testData <- generatePBTGSIdata(sampRate = .18, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = NA)

trapData <- testData[[1]]
tags <- testData[[2]]

# trapData <- trapData[,!grepl("Var", colnames(trapData))]
# trapData$GSI <- 1
head(trapData)


#data
groups <- c(tags[,1], paste0(sort(unique(trapData$GSI)))) #this is list of all the PBT groups and the wild
nGSI <- length(sort(unique(trapData$GSI)))
ohnc <- table(trapData$GenParentHatchery) # count how many of each pbt group
ohnc <- c(ohnc[tags[,1]], rep(0, nGSI)) # order counts the same as "groups" and add wild with 0
names(ohnc) <- groups # adding names, but not used by algorithm
t <- c(tags[,2], rep(0, nGSI)) # tag rates, with wild as 0
nUT <- sum(trapData$GenParentHatchery == "Unassigned") #this is the total number of untagged fish

#set up GSI observations for untagged fish
gsiUT <- trapData[trapData$GenParentHatchery == "Unassigned", "GSI"]

gsiValues <- sort(unique(trapData$GSI)) # these are the different GSI groups in the dataset

#set up GSI for tagged fish
ohnc_GSI <- matrix(0, nrow = nrow(tags), ncol = nGSI)
for(i in 1:nrow(tags)){
	for(j in 1:length(gsiValues)){
		ohnc_GSI[i,j] <- sum(trapData$GenParentHatchery == tags[i,1] & trapData$GSI == gsiValues[j])
	}
	
}

### set up categorical variables
names_variables <- colnames(trapData)[grepl("^Var", colnames(trapData))]
v_ut <- list() # this is list of the values for categorical variables for untagged fish
for(v in names_variables){
	v_ut[[v]] <- trapData[trapData$GenParentHatchery == "Unassigned", v] # these are the values for Var for untagged fish
}

values <- list() # this is list of the categories in each variable
for(v in names_variables){
	values[[v]] <- sort(unique(trapData[,v])) # these are the categories in var
}

V_ohnc <- list() # this is list of matrices of counts of values of each variable in the ohnc
for(v in names_variables){
	V1_ohnc <- matrix(0, nrow = length(groups), ncol = length(values[[v]])) # this is a matrix of counts of each Var1 value in the ohnc
	for(i in 1:length(groups)){
		tempData <- trapData[trapData$GenParentHatchery == groups[i],]
		for(j in 1:length(values[[v]])){
				V1_ohnc[i,j] <- sum(tempData[,v] == values[[v]][j])
		}
	}
	V_ohnc[[v]] <- V1_ohnc
}

# set up priors
# these are alphas for a Dirichlet prior for piTot
prior_piTot <- rep(1, length(groups))
#for all variables
prior_piV <- list()
for(v in names_variables){
	prior_piV[[v]] <- matrix(1, nrow = length(groups), ncol = length(values[[v]])) # here just using a uniform for all
}

#prior for GSI composition of PBT groups
prior_piGSI <- matrix(1, nrow = nrow(tags), ncol = nGSI)

#initial values
piTot <- rep(1/length(ohnc), length(ohnc))
z <- sample(groups, nUT, replace = TRUE)
oUT <- table(z)
oUT <- oUT[groups]

pi_V <- list()
for(v in names_variables){
	pi_V[[v]] <- matrix(1/length(values[[v]]), nrow = length(groups), ncol = length(values[[v]]))
}

#proportions of groups that GSI assign to each GSI category
pi_gsi <- matrix(0, nrow = length(groups), ncol = nGSI)
for(i in 1:nrow(tags)){
	pi_gsi[i,] <- 1/nGSI
}
currentCol <- 1
for(i in (nrow(tags)+1):length(groups)){ #setting GSI groups to be exactly as observed
	pi_gsi[i,currentCol] <- 1 #others were initialized at 0, so only set this one
	currentCol <- currentCol + 1
}

# number of iterations to run the sampler
reps <- 2000
# reps <- 10000

## other variables to pass
groupsInt <- 1:length(groups) #groups as ints instead of characters
nPBT <- nrow(tags)
zInt <- z
for(i in 1:length(groups)){
	zInt[zInt == groups[i]] <- i + length(groups)
}
zInt <- as.numeric(zInt)
zInt <- zInt - length(groups)

if(length(v_ut) == 0){
	v_utMat <- matrix(0,0,0)
} else {
	v_utMat <- matrix(nrow = length(v_ut[[1]]), ncol = 0)
	for(i in v_ut){
		v_utMat <- cbind(v_utMat, i)
	}
}


result <- MCpbt(iter = 2000, burnIn = 0 , thin = 1, seed = 7, #overall parameters for the chain
				     Nclip = NA, Nunclip = NA, clipPrior = NA, clippedBool = FALSE, #prop clipped parameters
				     piTotPrior = prior_piTot, ohnc = ohnc, piTotInitial = piTot, #piTotal parameters
				     oUTInitial = oUT, groups = groupsInt,
				     nPBT = nPBT, GSI_values = gsiValues, gsiUT = gsiUT, #pi_gsi parameters
				     pi_gsiInitial = pi_gsi, prior_pi_gsi = prior_piGSI,
				     ohnc_gsi = ohnc_GSI,
				     values = values, pi_VInitial = pi_V, pi_Vohnc = V_ohnc, pi_Vprior = prior_piV, #pi_V parameters
				     v_ut = v_utMat, 
				     initZ = zInt, t = t, #z parameters
         	   valuesOth = list(), pi_VInitialOth = list(), pi_VohncOth = list(), pi_VpriorOth = list(), #pi_VOth parameters
                 v_utOth = matrix(0,0,0)
				     )

result <- MCpbt(iter = 2000, burnIn = 0 , thin = 1, seed = 7, #overall parameters for the chain
				     Nclip = NA, Nunclip = NA, clipPrior = NA, clippedBool = FALSE, #prop clipped parameters
				     piTotPrior = prior_piTot, ohnc = ohnc, piTotInitial = piTot, #piTotal parameters
				     oUTInitial = oUT, groups = groupsInt,
				     nPBT = nPBT, GSI_values = gsiValues, gsiUT = gsiUT, #pi_gsi parameters
				     pi_gsiInitial = pi_gsi, prior_pi_gsi = prior_piGSI,
				     ohnc_gsi = ohnc_GSI,
				     values = list(), pi_VInitial = list(), pi_Vohnc = list(), pi_Vprior = list(), #pi_V parameters
				     v_ut = matrix(0,0,0), 
				     initZ = zInt, t = t, #z parameters
         	   valuesOth = values, pi_VInitialOth = pi_V, pi_VohncOth = V_ohnc, pi_VpriorOth = prior_piV, #pi_VOth parameters
                 v_utOth = v_utMat
				     )

str(result)
head(result$piTot)
unique(rowSums(result$piTot))
head(result$z)
head(result$piGSI)

unique(rowSums(result$piGSI[,1:3]))
unique(rowSums(result$piGSI[,4:6]))
unique(rowSums(result$piGSI[,7:9]))


head(result$piV[[1]])
unique(rowSums(result$piV[[1]][,1:2]))
unique(rowSums(result$piV[[1]][,3:4]))
unique(rowSums(result$piV[[1]][,5:6]))
unique(rowSums(result$piV[[1]][,7:8]))
unique(rowSums(result$piV[[1]][,9:10]))
unique(rowSums(result$piV[[1]][,11:12]))

head(result$piV[[2]])
unique(rowSums(result$piV[[2]][,1:2]))
unique(rowSums(result$piV[[2]][,3:4]))
unique(rowSums(result$piV[[2]][,5:6]))
unique(rowSums(result$piV[[2]][,7:8]))
unique(rowSums(result$piV[[2]][,9:10]))
unique(rowSums(result$piV[[2]][,11:12]))



#Checking results
plot(result$piTot[,5])

apply(result$piTot[500:reps,], 2, quantile, c(.05, .5, .95))
#           [,1]       [,2]       [,3]      [,4]      [,5]
# 5%  0.04396798 0.07923265 0.09881588 0.1770856 0.3221508
# 50% 0.05965463 0.10103489 0.12120468 0.2052440 0.3555022
# 95% 0.08031701 0.12395459 0.14673945 0.2344450 0.3887464
#          [,6]
# 5%  0.1290980
# 50% 0.1549492
# 95% 0.1837211

#           [,1]       [,2]       [,3]      [,4]      [,5]
# 5%  0.04280399 0.07906851 0.09858023 0.1781737 0.3212372
# 50% 0.05958460 0.10030789 0.11996547 0.2066282 0.3562417
# 95% 0.08098295 0.12320489 0.14647693 0.2360139 0.3895964
#          [,6]
# 5%  0.1297746
# 50% 0.1549512
# 95% 0.1827319



apply(result$piGSI[500:reps,], 2, quantile, c(.05, .5, .95))

apply(result$piV[[1]][500:reps,], 2, quantile, c(.05, .5, .95))

apply(result$piV[[2]][500:reps,], 2, quantile, c(.05, .5, .95))


pbtGSImat
varMat



## two variables
generatePBTGSIdata(sampRate = .18, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = varMat)

.3/6
.7/4




quantile(r_piTot[2,500:reps]/r_piTot[1,500:reps], c(.05, .5, .95))
quantile(r_piTot[3,500:reps]/r_piTot[1,500:reps], c(.05, .5, .95))

plot(r_pi_V[[1]][4,2,])

apply(r_pi_V[[1]][,1,500:reps], 1, quantile, c(.05, .5, .95))
apply(r_pi_V[[1]][,2,500:reps], 1, quantile, c(.05, .5, .95))

apply(r_pi_V[[2]][,1,500:reps], 1, quantile, c(.05, .5, .95))
apply(r_pi_V[[2]][,2,500:reps], 1, quantile, c(.05, .5, .95))


table(trapData$GenParentHatchery, trapData$Var1)
apply(r_z[v_ut[[1]] == 1,500:reps][1:10,], 1, table)
apply(r_z[v_ut[[1]] == 2,500:reps][1:10,], 1, table)

hist(r_piTot[4,500:reps])



###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
#### look at piTot CIs

countright <- 0
countwrong <- 0

for(r in 1:100){
	## no variables
	testData <- generatePBTGSIdata(sampRate = .18, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
											 obsTagRates = c(.8, .85,.9), physTagRates = 0,
					    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = NA)
	
	trapData <- testData[[1]]
	tags <- testData[[2]]
	
	u1 <- .3/6
	u2 <- .7/4
	truePiTot <- c(u1, u1*2, u1*3, u2, u2*2, u2)
	
	
	
	#data
	groups <- c(tags[,1], paste0(sort(unique(trapData$GSI)))) #this is list of all the PBT groups and the wild
	nGSI <- length(sort(unique(trapData$GSI)))
	ohnc <- table(trapData$GenParentHatchery) # count how many of each pbt group
	ohnc <- c(ohnc[tags[,1]], rep(0, nGSI)) # order counts the same as "groups" and add wild with 0
	names(ohnc) <- groups # adding names, but not used by algorithm
	t <- c(tags[,2], rep(0, nGSI)) # tag rates, with wild as 0
	nUT <- sum(trapData$GenParentHatchery == "Unassigned") #this is the total number of untagged fish
	
	#set up GSI observations for untagged fish
	gsiUT <- trapData[trapData$GenParentHatchery == "Unassigned", "GSI"]
	
	gsiValues <- sort(unique(trapData$GSI)) # these are the different GSI groups in the dataset
	
	#set up GSI for tagged fish
	ohnc_GSI <- matrix(0, nrow = nrow(tags), ncol = nGSI)
	for(i in 1:nrow(tags)){
		for(j in 1:length(gsiValues)){
			ohnc_GSI[i,j] <- sum(trapData$GenParentHatchery == tags[i,1] & trapData$GSI == gsiValues[j])
		}
		
	}
	
	### set up categorical variables
	names_variables <- colnames(trapData)[grepl("^Var", colnames(trapData))]
	v_ut <- list() # this is list of the values for categorical variables for untagged fish
	for(v in names_variables){
		v_ut[[v]] <- trapData[trapData$GenParentHatchery == "Unassigned", v] # these are the values for Var for untagged fish
	}
	
	values <- list() # this is list of the categories in each variable
	for(v in names_variables){
		values[[v]] <- sort(unique(trapData[,v])) # these are the categories in var
	}
	
	V_ohnc <- list() # this is list of matrices of counts of values of each variable in the ohnc
	for(v in names_variables){
		V1_ohnc <- matrix(0, nrow = length(groups), ncol = length(values[[v]])) # this is a matrix of counts of each Var1 value in the ohnc
		for(i in 1:length(groups)){
			tempData <- trapData[trapData$GenParentHatchery == groups[i],]
			for(j in 1:length(values[[v]])){
					V1_ohnc[i,j] <- sum(tempData[,v] == values[[v]][j])
			}
		}
		V_ohnc[[v]] <- V1_ohnc
	}
	
	# set up priors
	# these are alphas for a Dirichlet prior for piTot
	prior_piTot <- rep(1, length(groups))
	#for all variables
	prior_piV <- list()
	for(v in names_variables){
		prior_piV[[v]] <- matrix(1, nrow = length(groups), ncol = length(values[[v]])) # here just using a uniform for all
	}
	
	#prior for GSI composition of PBT groups
	prior_piGSI <- matrix(1, nrow = nrow(tags), ncol = nGSI)
	
	#initial values
	piTot <- rep(1/length(ohnc), length(ohnc))
	z <- sample(groups, nUT, replace = TRUE)
	oUT <- table(z)
	oUT <- oUT[groups]
	
	pi_V <- list()
	for(v in names_variables){
		pi_V[[v]] <- matrix(1/length(values[[v]]), nrow = length(groups), ncol = length(values[[v]]))
	}
	
	#proportions of groups that GSI assign to each GSI category
	pi_gsi <- matrix(0, nrow = length(groups), ncol = nGSI)
	for(i in 1:nrow(tags)){
		pi_gsi[i,] <- 1/nGSI
	}
	currentCol <- 1
	for(i in (nrow(tags)+1):length(groups)){ #setting GSI groups to be exactly as observed
		pi_gsi[i,currentCol] <- 1 #others were initialized at 0, so only set this one
		currentCol <- currentCol + 1
	}
	
	# number of iterations to run the sampler
	reps <- 2000
	# reps <- 10000
	
	## other variables to pass
	groupsInt <- 1:length(groups) #groups as ints instead of characters
	nPBT <- nrow(tags)
	zInt <- z
	for(i in 1:length(groups)){
		zInt[zInt == groups[i]] <- i + length(groups)
	}
	zInt <- as.numeric(zInt)
	zInt <- zInt - length(groups)
	
	if(length(v_ut) == 0){
		v_utMat <- matrix(0,0,0)
	} else {
		v_utMat <- matrix(nrow = length(v_ut[[1]]), ncol = 0)
		for(i in v_ut){
			v_utMat <- cbind(v_utMat, i)
		}
	}
	
	
	result <- MCpbt(iter = 2000, burnIn = 0 , thin = 1, seed = 7, #overall parameters for the chain
					     Nclip = NA, Nunclip = NA, clipPrior = NA, clippedBool = FALSE, #prop clipped parameters
					     piTotPrior = prior_piTot, ohnc = ohnc, piTotInitial = piTot, #piTotal parameters
					     oUTInitial = oUT, groups = groupsInt,
					     nPBT = nPBT, GSI_values = gsiValues, gsiUT = gsiUT, #pi_gsi parameters
					     pi_gsiInitial = pi_gsi, prior_pi_gsi = prior_piGSI,
					     ohnc_gsi = ohnc_GSI,
					     values = values, pi_VInitial = pi_V, pi_Vohnc = V_ohnc, pi_Vprior = prior_piV, #pi_V parameters
					     v_ut = v_utMat, 
					     initZ = zInt, t = t #z parameters
					     )
	
	
	cis <- apply(result$piTot[500:reps,], 2, quantile, c(.05, .5, .95))
	
	for(i in 1:6){
		if(cis[1,i] <= truePiTot[i] && cis[3,i] >= truePiTot[i]){
			countright <- countright + 1
		} else {
			countwrong <- countwrong + 1
		}
	}
}

# this is for a 90% CI
countright / (countright + countwrong)
# [1] 0.9016667
countwrong / (countright + countwrong)
# [1] 0.09833333

# BEAUTIFUL!!!

#################################################################################################################
#################################################################################################################
## testing data prep and wrapper functions

pbtGSImat <- matrix(c(.1, .8, .1,.8, .1, .1,.1, .1, .8), nrow = 3, ncol = 3, byrow = TRUE)
varMat <- list(
	matrix(c(rep(c(.1,.9), 3), rep(c(.9,.1), 3)), nrow = 6, ncol = 2, byrow = TRUE),
	matrix(c(rep(c(.4,.6), 3), rep(c(.6,.4), 3)), nrow = 6, ncol = 2, byrow = TRUE)
)

## two variables
testData2 <- generatePBTGSIdata(sampRate = .2, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = varMat)
testData2[[1]][,"AdClip"] <- testData2[[1]][,"AdClip"] == "AI"
nrow(testData2[[1]])

input <- prepOneStrataAI(testData2[[1]], testData2[[2]], "GSI", "GenParentHatchery", variableCols = c("Var1"), 
								 variableColsOth = c("Var2"), adFinCol = "AdClip", AI = TRUE, verbose = TRUE, GSIgroups = NA,
									 variableValues = NA, variableValuesOth = NA)

head(testData2[[1]])

#with multiple strata
nStrata <- 2
multStratData <- data.frame()
for(i in 1:nStrata){
	tempData <- generatePBTGSIdata(sampRate = .2, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = .4, true_noclip_H = .2, true_wild = .4, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = varMat)
	tags <- tempData[[2]]
	tempData <- tempData[[1]]
	tempData$Strata <- i
	multStratData <- rbind(multStratData, tempData)
	
}
tags
head(multStratData)

clipInput <- prepStrataPropClip(multStratData, "Strata", "AdClip")

input <- prepStrata(multStratData, tags, "GSI", "GenParentHatchery", "Strata", variableCols = c("Var1", "Var2"), variableColsOth = c(), "AdClip",
								AI = TRUE, GSIgroups = NA,
									 variableValues = NA, variableValuesOth = NA, verbose = TRUE)
	
str(input)
input[[1]]
input[[2]]

names(input[[1]])


aiRes <- estimStrataMCpbt(input, 1000, 10, 1)
clipRes <- estimStrataPropClip(clipInput, nrow(aiRes[[1]]$piTot), seed = NA, priors = NA, verbose = TRUE)



str(aiRes)
str(clipRes)

# escapement estimates for each strata
escapementByStrata <- list(rep(3000, nrow(aiRes[[1]]$piTot)),
									rep(3000, nrow(aiRes[[1]]$piTot))
									)


escapeResults <- multByEscapement(input, aiRes, clipRes, escapementByStrata, verbose = TRUE)
# escapeResults <- multByEscapement(input, aiRes, NA, escapementByStrata, verbose = TRUE)

str(escapeResults)

escapeResults$totPiTotEstim
rowSums(escapeResults$totPiTotEstim)
unique(rowSums(escapeResults$totClipUnclip))
unique(rowSums(escapeResults$totPiTotEstim) - escapeResults$totClipUnclip[,2])

head(escapeResults$totClipUnclip)

str(escapeResults$strataEstimates)
str(escapeResults$strataEstimates[[1]]$piTotEstim)
head(escapeResults$strataEstimates[[1]]$piTotEstim)


generatePBTGSIdata(sampRate = .2, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = .4, true_noclip_H = .2, true_wild = .4, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat)


######################################################################################################
######################################################################################################
######################################################################################################

setwd("S:\\Eagle Fish Genetics Lab\\Tom\\Scobi_Deux\\OtsBON_SCOBI2018")


##Now load the three data sources into R

tags <- read.csv("PBT tag rates.csv", stringsAsFactors = FALSE)
tags <- tags[tags[,1] != "Unassigned",] #remove Unassigned from tag rate dataframe
trap <- read.csv("2018RearHNCscobi_TD.csv", stringsAsFactors = FALSE)
window <- read.csv("HNC_clip_window.csv", stringsAsFactors = FALSE)

trap$WeekNumber <- as.character(trap$WeekNumber)

#just showing what has been loaded into R
head(tags)
head(trap)
head(window)

# note the AI=TRUE option b/c these are ad-intact fish we are analyzing
mainInput <- prepStrata(trapData = trap, tags = tags, GSIcol = "GenStock", PBTcol = "GenParentHatchery", strataCol = "WeekNumber", adFinCol = "AdClip", AI = TRUE)

str(tags)

propEstimates <- estimStrataMCpbt(mainInput, iter = 11000, burnIn = 1000, thin = 1, seed = 7)






## one variable
testData <- generatePBTGSIdata(sampRate = .18, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = 0, true_noclip_H = 0, true_wild = 1, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = varMat[1])

testData[[1]]$Strata <- 1
mainInput <- prepStrata(trapData = testData[[1]], tags = testData[[2]], GSIcol = "GSI", PBTcol = "GenParentHatchery", strataCol = "Strata", 
								adFinCol = "AdClip", AI = TRUE)


propEstimates <- estimStrataMCpbt(mainInput, iter = 1000, burnIn = 100, thin = 1, seed = 7)

countEstimates <- multByEscapement(mainInput, mainRes = propEstimates, popSizes = popSizeEstimates,
writeSummary = TRUE)
