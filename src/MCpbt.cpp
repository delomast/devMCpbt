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
// [[Rcpp::export]]
Rcpp::NumericVector MCpbt(int iter, int burnIn, int thin, unsigned int seed, //overall parameters for the chain
                          int Nclip, int Nunclip, Rcpp::NumericVector clipPrior, bool clippedBool, //prop clipped parameters
                          Rcpp::NumericVector piTotPrior, Rcpp::NumericVector ohnc, Rcpp::NumericVector piTotInitial, //piTotal parameters
                          Rcpp::NumericVector oUTInitial
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
	int nGroups = ohnc.size(); //number of groups in piTot, including wild/unassigned
	vector <double> priorPlusOhnc (nGroups); //creating prior plus observed PBT assigned
	for(int i=0; i < nGroups; i++){
		priorPlusOhnc[i] = piTotPrior[i] + ohnc[i];
	}
	
	//define temporary variables used in calculations
	vector <double> piTot_tempAlphas (nGroups, 0);
	
	//define main variables and set initial values
	vector <double> piTot (nGroups, 0); //proportions of each group in the population
	vector <int> oUT (nGroups,0); //number of untagged fish in each group
	for(int i=0; i < nGroups; i++){
		piTot[i] = piTotInitial[i];
		oUT[i] = oUTInitial[i];
		
	}
		/*now for the variables
		 we need:
		 	a structure to be able to loop through the variables
		 	then, for each variable:
		 	a list of values
			a list of proportions for each group, or a matrix
		 	
		 */ 
	vector <vector> 
	
	
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
			
			
		for(v in names_variables){
			for(i in 1:nrow(pi_V[[v]])){
				utVcounts <- rep(0, length(values))
				tempData <- v_ut[[v]][z == groups[i]]
				for(j in 1:length(values[[v]])){  #//can't use table in case any is zero
					utVcounts[j] <- sum(tempData == values[[v]][j])
				}
				#//conditional is dirichlet with prior alphas and observations added up
				tempAlphas <- prior_piV[[v]][i,] + V_ohnc[[v]][i,] + utVcounts
				pi_V[[v]][i,] <- rdirichlet(1,tempAlphas)
			}
		}
			
			// sample from z's
			// calculate oUT for ease of later steps
			
		}
		//record values
		
	}
	
	
	// int dims = alphas.size();
	// Rcpp::NumericVector t(dims);
	// t = randDirich(alphas, rgPoint); //random beta computed from two random gammas

	//organize output to return to R
	
	Rcpp::NumericVector t(1);
	t[0] = 3;
	return t;
}
