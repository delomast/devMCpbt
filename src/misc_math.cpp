#include <Rcpp.h>
#include <random>
#include <vector>
using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// this function generates random numbers from a beta distribution
//not exported
double randBeta(double alpha, double beta, mt19937 * rNum) {

	//input checking
	if (alpha + beta == 0 || alpha < 0 || beta < 0){
		Rcpp::stop("invalid alpha or beta to randBeta");
	}
	
	// initiate random gammas
	gamma_distribution<double> rGamma1 (alpha,1); //alpha of beta distribution
	gamma_distribution<double> rGamma2 (beta,1); //beta of beta distribution

	double x1;
	double x2;
	x1 = rGamma1(*rNum);
	x2 = rGamma2(*rNum);

	return x1 / (x1 + x2); //random beta computed from two random gammas
}

// this function generates random numbers from a Dirichlet distribution
//not exported
vector<double> randDirich(vector<double> alphas, mt19937 * rNum) {

	int dim = alphas.size(); //get number of dimensions of requested Dirichlet
	
	//input checking
	double sumA = 0;
	for(int i=0; i < dim; i++){
		if(alphas[i] < 0){
			Rcpp::stop("invalid alphas to randDirichlet");
		}
		sumA += alphas[i];
	}
	if (sumA == 0){
		Rcpp::stop("invalid alphas to randDirichlet");
	}
	
	// allocate storage of random gammas
	vector<double> rGammaValues (dim);
	double sum = 0; //sum for normalizing
	//sample gammas
	for (int i = 0; i < dim; i++){
		gamma_distribution<double> rGamma1 (alphas[i],1);
		rGammaValues[i] = rGamma1(*rNum);
		sum += rGammaValues[i]; //calculate sum
	}
	//normalize
	for (int i = 0; i < dim; i++){
		rGammaValues[i] = rGammaValues[i] / sum;
	}

	return rGammaValues; //random Dirichlet computed from random gammas
}
