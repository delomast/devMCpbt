#include <Rcpp.h>
#include <random>
#include <vector>
#include "misc_math.h"
#include "MCMCclip.h"

using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

//This is the main function to run the MCMC composition model
//
//right now it is in development and so the input and output chagne based on checking
//
//groups is numeric vector of all the groups ie 1,2,3,...in the same order as all other group based variables
//
//values is a list of NumericVectors giving the different categories for each variable
//pi_VInitial is a list of NumericMatrices giving the composition of each variable(matrix) by group(rows) and values (column)
//pi_Vohnc is a list of NumericMatrices giving the counts of ohnc each variable(matrix) by group(rows) and values (column)
//  It includes the wild, with counts of zero
//pi_Vprior is a list of NumericMatrices giving the alphas for a Dirichlet prior of each variable(matrix) by group(rows) and values (column)
//v_ut is a matrix with a row for each untagged individual, a column for each variable, and the data being the values for each variable.
//  -9 is missing
//t is tag rates in same order as groups
// [[Rcpp::export]]
Rcpp::NumericVector MCpbt(int iter, int burnIn, int thin, unsigned int seed, //overall parameters for the chain
                          int Nclip, int Nunclip, Rcpp::NumericVector clipPrior, bool clippedBool, //prop clipped parameters
                          Rcpp::NumericVector piTotPrior, Rcpp::NumericVector ohnc, Rcpp::NumericVector piTotInitial, //piTotal parameters
                          Rcpp::NumericVector oUTInitial, Rcpp::NumericVector groups,
                          Rcpp::List values, Rcpp::List pi_VInitial, Rcpp::List pi_Vohnc, Rcpp::List pi_Vprior, //pi_V parameters
                          Rcpp::NumericMatrix v_ut, 
                          Rcpp::NumericVector initZ, Rcpp::NumericVector t //z parameters
                          )
{
	// initiate random number generator
	mt19937 rg (seed);
	mt19937 * rgPoint = &rg; //pointer to random number generator to pass to subfunctions
	

	//allocate result storage
	int NumResults = (iter - burnIn) / thin; //number of results to store
	Rcpp::NumericVector r_Propclip(NumResults, NA_REAL); //prop clipped
	
	////////////////////////////////
	//get prop clipped estimates
	////////////////////////////////
	if (clippedBool){
		r_Propclip = MCMCclip(Nclip, Nunclip, clipPrior, NumResults, rgPoint);
	}
	
	////////////////////////////////
	//get other estimates
	///////////////////////////////
	
	//set quantities frequently used
	int nGroups = groups.size(); //number of groups in piTot, including wild/unassigned
	vector <double> priorPlusOhnc (nGroups); //creating prior plus observed PBT assigned
	vector <double> untagRates (nGroups);
	for(int i=0; i < nGroups; i++){
		priorPlusOhnc[i] = piTotPrior[i] + ohnc[i];
		untagRates[i] = 1 - t[i];
	}
	int tempCol; //stores position of vector to take a value from
	
	
	//define temporary variables used in calculations
	vector <double> piTot_tempAlphas (nGroups, 0);
	vector <double> tempZprobs (nGroups, 0);
	
	//define main variables and set initial values
	vector <double> piTot (nGroups, 0); //proportions of each group in the population
	vector <int> oUT (nGroups,0); //number of untagged fish in each group
	for(int i=0; i < nGroups; i++){
		piTot[i] = piTotInitial[i];
		oUT[i] = oUTInitial[i];
		
	}
	vector <int> z; //declare vector of z's
	for(int i=0; i < initZ.size(); i++){
		z.push_back(initZ[i]); //initialize with given values
	}
	int nZ = z.size(); //number of untagged fish
	/*now for the variables
	 we need:
	 	a structure to be able to loop through the variables - 0:nVar
	 	then, for each variable - set up as vectors for each need, with element 0 being var1, 1 being var2, etc...
		 	a list of values - valuesC
			a list of proportions for each group (a matrix of proportions within groups) - pi_V
	 		a temporary list of counts, reset for each group when calculating - could allocate this on the fly, but that's a lot of allocating

	  		a prior for each group - combination of prior and observed (ohnc) values
	 Instead of a specific structure to loop through, just keep the order the same for all
	 */ 
	int nVar = values.size(); //number of variables
	vector <vector <int>> valuesC (nVar); //possible values of variables represented as (positive) ints, (with -9 in data indicates missing)
	vector <vector <vector <double>>> pi_V (nVar); //a list each variable, with proportions of each group for that variable
	vector <vector <int>> tempCounts (nVar); //temporary counts for each variable - used in calculations

	vector <vector <vector <double>>> priorPlusOhncV (nVar); //creating prior plus observed PBT assigned

	for(int i=0; i<nVar; i++){
		//add values to valuesC
		Rcpp::NumericVector tempV = values[i]; //pull NumericVector out of the list of the possible values for a variable
		valuesC[i].assign(tempV.size(), 0);
		for(int j=0; j<tempV.size(); j++){
			valuesC[i][j] = tempV[j];
		}
		//add values to temporary counts
		tempCounts[i].assign(tempV.size(), 0); //don't really need to do this, but good practice in case later changes to code
		//add initial values to pi_V and priorPlusOhncV
		Rcpp::NumericMatrix tempVmat = pi_VInitial[i]; //this is matrix of initial values for pi_V of that variable
		Rcpp::NumericMatrix tempVohnc = pi_Vohnc[i]; //this is matrix of ohnc counts for pi_V of that variable
		Rcpp::NumericMatrix tempVprior = pi_Vprior[i];//this is matrix of alphas for Dirichlet prior for pi_V of that variable
		for(int j=0; j < nGroups; j++){ //for each row - group
			vector <double> tempVec (tempV.size()); //make empty vector to accept values for one group
			vector <double> tempVecPrior (tempV.size()); //make empty vector to accept values for one group
			for(int k=0; k < tempV.size(); k++){//for each column -variable value
				tempVec[k] = tempVmat(j,k); //take value from initial values matrix and assign to vector
				tempVecPrior[k] = tempVohnc(j,k) + tempVprior(j,k); //take value from initial values matrix and assign to vector
			}
			pi_V[i].push_back(tempVec); //add initial values for variable i and group j
			priorPlusOhncV[i].push_back(tempVec); //add prior + ohnc values for variable i and group j
		}
	}
	

	//cycle through iterations
	for (int r=0; r < iter; r++){
		//cycle through thinning reps between recording values
		for (int t=0; t < thin; t++){
			
			// sample from pi_tot
			for(int i=0; i < nGroups; i++){ //calculating the alphas of the Dirichlet
				piTot_tempAlphas[i] = priorPlusOhnc[i] + oUT[i];
			}
			piTot = randDirich(piTot_tempAlphas, rgPoint);
			
			// sample from pi_v's
			for(int v=0; v < nVar; v++){ //cycle through variables
				for(int g=0; g < nGroups; g++){ //cycle through groups
					tempCounts[v].assign(valuesC[v].size(),0); //zero the vector
					for(int i=0; i < valuesC[v].size(); i++){ //cycle through categories to count
						//now count all in that group with that value for that variable
						for(int j=0; j < nZ; j++){ //for each untagged individual
							if(z[j] == groups[g] && v_ut(j,v) == valuesC[v][i]) //v_ut is kept as a NumericMatrix b/c no manipulation needed. row is indiv col is variable
								tempCounts[v][i] += 1; //this is utVcounts equivalent in the R script
						}
					}
					//now tempCounts has the counts for that group and variable
					//so now add to ohnc counts, prior, and then sample from Dirichlet
					for(int i=0; i < valuesC[v].size(); i++){
						tempCounts[v][i] += priorPlusOhncV[v][g][i];
					}
					pi_V[v][g] = randDirich(tempCounts, rgPoint);
				}
			}
			
			// sample from z's
			for(int i=0; i < nZ; i++){ //for each z
				//first assign just the piTot and tag rates to the probs
				for(int g=0; g < nGroups; g++){
					tempZprobs[g] = piTot[g] * untagRates[g];
				}
				//now refine probs using pi_V and the observed variables
				for(int v=0; v < nVar; v++){
					if(v_ut(i,v) == -9) continue; //if value is missing, skip
					for(int g=0; g < nGroups; g++){ //for each group
						//find position of that value - they should be in order, so if speed is a big issue, can replace this just using the value as the position
						for(int k=0; k < valuesC[v].size(); k++){
							if(valuesC[v][k] == v_ut(i,v)){
								tempCol = k;
								break;
							}
						}
						// find probability of that value for that group and multiply
						tempZprobs[g] *= pi_V[v][g][k];
					}
				}
				// now sample z for that individual
				/* left off here
				 * Maybe need to define a function that is equivalent of R's sample
				 * This could help streamline the code
				 */
			}
			
			
#sample from z's
	for(i in 1:length(z)){
		tempZprobs <- piTot * (1-t) # these are probs for each group not taking into account observed variables
		for(v in names_variables){
			if(is.na(v_ut[[v]][i])){ # if value is missing, ignore this variable for this individual
				next
			}
			tempCol <- which(values[[v]] == v_ut[[v]][i]) #this is the column of pi_V that has the proportions for the category observed in that individual
			tempZprobs <- tempZprobs*pi_V[[v]][,tempCol]
		}
		z[i] <- sample(groups, 1, prob = tempZprobs) #tempZprobs is normalized by R
	}
			
			
			// calculate oUT for ease of later steps
			
		}
		//record values
		
		//check user interrupt every 500 recordings
		if (r % 500 == 0){
			Rcpp::checkUserInterrupt();
		}
	}
	
	
	// int dims = alphas.size();
	// Rcpp::NumericVector t(dims);
	// t = randDirich(alphas, rgPoint); //random beta computed from two random gammas

	//organize output to return to R
	
	Rcpp::NumericVector t(1);
	t[0] = 3;
	return t;
}
