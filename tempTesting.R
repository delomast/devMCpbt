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
