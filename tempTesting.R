### testing script for devMCpbt
# not built as part of package



MCpbt(1, 1, 5, 1)
MCpbt(100, 1000, 7, 10)

rec <- c()
for(i in 1:10000){
	rec <- c(rec,MCpbt(c(1,1,1,1), i)[1])
}
MCpbt(c(1,1,1,1,1), 2)
hist(rec)
summary(rec)

library(MCMCpack)
summary(rdirichlet(10000,c(1,1,1,1))[,1])
rdirichlet(10,c(1,1,1,1,1))
hist(rdirichlet(10000,c(1,1,1,1))[,1])

# storage calculations
iter = 150
burnin = 30
thin = 2
count=0
for(i in 0:(iter - 1)){
	if(i< burnin){next}
	if((i - (burnin -1)) %% thin == 0){
		count = count + 1
	}
}
count
(iter - burnin) %/% thin

# testing clip routine
hist(MCMCclip(100, 500, c(1,1), 15000, 30, 2, 7))

system.time(MCMCclip(100, 500, c(1,1), 15000, 30, 2, 7))
MCMCclip(10, 50, c(1,1), 15000, 30, 2, 7)

#### make some data
subCatArray <- array(dim = c(4,3,3))
### variable 1
subCatArray[1:3,,1] <- .1
subCatArray[1:3,,2] <- .8
subCatArray[1:3,,3] <- .1
subCatArray[4,,1] <- .4
subCatArray[4,,2] <- .2
subCatArray[4,,3] <- .4

subCatArray2 <- array(dim = c(4,3,3))
### variable 2
subCatArray2[1:3,,1] <- .7
subCatArray2[1:3,,2] <- .2
subCatArray2[1:3,,3] <- .1
subCatArray2[4,,1] <- .2
subCatArray2[4,,2] <- .2
subCatArray2[4,,3] <- .6

subCats <- list(subCatArray, subCatArray2)

testData <- rearTypeMultVarData(sampRates = c(.1, .1, .5), censusSizes = c(3000, 2500, 500), relSizePBTgroups = c(1,2,3), tagRates = c(.8, .9, .95), 
						 obsTagRates = c(.8, .9, .95), 
					physTagRates = c(0,0,0), true_clipped = c(0, 0, 0), true_noclip_H = c(.3, .3, .3), true_wild = c(.7, .7, .7),
					varArray = subCats)
trapData <- testData[[1]]
window <- testData[[2]]
tags <- testData[[3]]

#data
groups <- c(tags[,1], "wild") #this is list of all the PBT groups and the wild
ohnc <- table(trapData$GenParentHatchery) # count how many of each pbt group
ohnc <- c(ohnc[tags[,1]],0) # order counts the same as "groups" and add wild with 0
names(ohnc) <- groups # adding names, but not used by algorithm
t <- c(tags[,2], 0) # tag rates, with wild as 0
nUT <- sum(trapData$GenParentHatchery == "Unassigned") #this is the total number of untagged fish

### set up categorical variables
names_variables <- c("Var1", "Var2")
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

#initial values
piTot <- rep(1/length(ohnc), length(ohnc))
z <- sample(groups, nUT, replace = TRUE)
oUT <- table(z)
oUT <- oUT[groups]

#need to turn into a list, one for each var
pi_V <- list()
for(v in names_variables){
	pi_V[[v]] <- matrix(1/length(values[[v]]), nrow = length(groups), ncol = length(values[[v]]))
}

groupsInt <- 1:length(groups) #groups as ints instead of characters
v_utMat <- matrix(nrow = length(v_ut[[1]]), ncol = 0)
for(i in v_ut){
	v_utMat <- cbind(v_utMat, i)
}
zInt <- z
for(i in 1:length(groups)){
	zInt[zInt == groups[i]] <- i
}
zInt <- as.numeric(zInt)


# number of iterations to run the sampler
reps <- 5000
burn <- 500
thin <- 5
system.time(
result <- MCpbt(iter = reps, burnIn = burn, thin = thin, seed = 7, #overall parameters for the chain
                          Nclip = 1000, Nunclip = 756, clipPrior = c(1,1), clippedBool = TRUE, #prop clipped parameters
                          piTotPrior = prior_piTot, ohnc = ohnc, piTotInitial = piTot, #piTotal parameters
                          oUTInitial = oUT, groups = groupsInt,
                          values = values, pi_VInitial = pi_V, pi_Vohnc = V_ohnc, pi_Vprior = prior_piV, #pi_V parameters
                          v_ut = v_utMat, 
                          initZ = zInt, t = t #z parameters
                          )
)

apply(result, 2, quantile, c(.05, .5, .95))

plot(result[,4])

nrow(result)

##########################################
### after changes to make GSI groups separate
###############

pbtGSImat <- matrix(c(.1, .8, .1,.8, .1, .1,.1, .1, .8), nrow = 3, ncol = 3, byrow = TRUE)
varMat <- list(
	matrix(c(rep(c(.1,.9), 3), rep(c(.9,.1), 3)), nrow = 6, ncol = 2, byrow = TRUE),
	matrix(c(rep(c(.4,.6), 3), rep(c(.6,.4), 3)), nrow = 6, ncol = 2, byrow = TRUE)
)

testData <- generatePBTGSIdata(sampRate = .18, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = varMat)
testData <- generatePBTGSIdata(sampRate = .18, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = varMat[1])
testData <- generatePBTGSIdata(sampRate = .18, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = NA)

trapData <- testData[[1]]
tags <- testData[[2]]



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



#Checking results
plot(result[,5])

apply(result[500:reps,], 2, quantile, c(.05, .5, .95))




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

