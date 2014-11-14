/*
 * functions.c
 *
 *  Created on: Jan 27, 2014
 *      Author: MAPF
 */

#include<math.h>

#include"lnLfunctions.h"

// Define functions

#define ZERO 1e-4
#define BIG 1e100

/****************** ******************/
//Harmonic number Hn calculated "crudely" approximating 1 + 1/2 + ...+ 1/k
// TODO Maybe there is already implemented harmonic number
double Hn (int k){
	double sum = 0.0;
	int i;
	for(i=1; i<=k; i++){
		sum += 1.0/i;
	}

	return(sum);
}

double log_NegativeBinomial(int k, double mu, double shape){

	if(k==0 && mu==0) { return(0);}

	//MyImplementation
	// transform mu to probability
	double p = (shape/(shape+mu));

	// dnbinom according to R documentation
	// Obs: tgamma(k+1)= k!
	double res = (tgamma(k + shape)/(tgamma(shape) * tgamma(k + 1))) * pow(p,shape) * pow((1-p),k);

	return(log(res));
}

// j and n are int,
double XP (double j, double n, double S){
	return(n/(j*(n-j))*(1-exp(-S *(1-j/n))) /(1-exp(-S)));
}

// j and n are int,
double sumXP (double n, double S){
	unsigned int i;
	double sumXP = 0.0;

	for (i=1; i<n; i++){
		sumXP += XP(i, n, S);
	}

	return sumXP;

}

//EPs, ESs and EDs
double EPs(double Ls, double r, double theta, int n){
	return(Ls*r*theta*Hn(n-1));
}

double ESs(double Ls, double theta, double n){
	return(Ls*theta*(n/(n-1)));
}

double EDs(double Ls, double lambda, double theta, int n){
	return(Ls*(lambda + theta/n));
}

/**************************************************/

// EPn_ana, ESn_ana, EDn_ana
double EPn_ana (double Ln, double theta, double r, int n, double S){

	double epn = 0.0;

	if (fabs(S) < 0){ S = ZERO; }

	//calculate intrasum:
	double sum = sumXP(n, S);

	epn = Ln*r*sum*theta;

	return(epn);
}

double ESn_ana (double Ln, double theta, int n, double S){

	double esn;

	if(fabs(S) < ZERO) {S = ZERO;} // if fabs(S) < Zero we force S to the small number Zero
	esn = Ln*theta*(XP(1,n,S)+XP((n-1),n,S));
	return(esn);

}

double EDn_ana (double Ln, double lambda, double theta, double r, int n, double S){
	double edn = 0.0;

	if(fabs(S) < ZERO){S = ZERO;}
	edn = Ln*(lambda*(S/(1-exp(-S)))+ 2*(1-exp(-S/(2*n)))/(1-exp(-S)) * theta*r);
	return(edn);

}

/*************************************************/

double phiS (double smean, double smax, double shapeBeta, double S){

	double res = 0.0;
	res = exp((smax-S)*(shapeBeta)/smean)*pow((smax-S),(shapeBeta-1))*pow((-smean/shapeBeta),(-shapeBeta))/tgamma(shapeBeta);

	return(res);

}

double expS (double Sb, double S){
	double res = 0.0;
	res = Sb * exp(-Sb * S);

	return(res);
}

/**************************************************/
// EPn_ana, ESn_ana, EDn_ana OVER phi(S)

double EPn_OverPhiS_adapt(unsigned dim, const double *x, void *fparams){

	//params *par = *(struct params **)fparams;

	double smean = ((struct params_s *)fparams)->smean;
	double smax = ((struct params_s *)fparams)->smax;
	double shapeBeta = ((struct params_s *)fparams)->shapeBeta;

	double Ln = ((struct params_s *)fparams)->Ln;
	double theta = ((struct params_s *)fparams)->theta;
	double r = ((struct params_s *)fparams)->r;
	double n = ((struct params_s *)fparams)->n;

	unsigned i;

    double prod = 1.0;
    for (i = 0; i < dim; ++i)
	  prod *= phiS(smean, smax, shapeBeta, x[i]) * EPn_ana(Ln, theta, r, n, x[i]);

    return prod;
}

double ESn_OverPhiS_adapt(unsigned dim, const double *x, void *fparams){

	double smean = ((struct params_s *)fparams)->smean;
	double smax = ((struct params_s *)fparams)->smax;
	double shapeBeta = ((struct params_s *)fparams)->shapeBeta;

	double Ln = ((struct params_s *)fparams)->Ln;
	double theta = ((struct params_s *)fparams)->theta;
	double n = ((struct params_s *)fparams)->n;

    unsigned int i;

	double prod = 1.0;

    for (i = 0; i < dim; ++i){
    	prod *= phiS(smean, smax, shapeBeta, x[i]) * ESn_ana(Ln, theta, n, x[i]);
    }
    return prod;

}

double EDn_OverPhiS_adapt(unsigned dim, const double *x, void *fparams){

	double smean = ((struct params_s *)fparams)->smean;
	double smax = ((struct params_s *)fparams)->smax;
	double shapeBeta = ((struct params_s *)fparams)->shapeBeta;

	double Ln = ((struct params_s *)fparams)->Ln;
	double lambda = ((struct params_s *)fparams)->lambda;
	double theta = ((struct params_s *)fparams)->theta;
	double r = ((struct params_s *)fparams)->r;
	double n = ((struct params_s *)fparams)->n;

    unsigned int i;

	double prod = 1.0;
    for (i = 0; i < dim; ++i){
    	prod *= phiS(smean, smax, shapeBeta, x[i]) * EDn_ana(Ln, lambda, theta, r, n, x[i]);
    }

    return prod;

}

/**************************************************/
// EPn_ana, ESn_ana, EDn_ana OVER Exp(S)

double EPn_OverExpS_adapt(unsigned dim, const double *x, void *fparams){

	double Sb = ((struct params_s *)fparams)->Sb;

	double Ln = ((struct params_s *)fparams)->Ln;
	double theta = ((struct params_s *)fparams)->theta;
	double r = ((struct params_s *)fparams)->r;
	double n = ((struct params_s *)fparams)->n;

	unsigned i;

    double prod = 1.0;
    for (i = 0; i < dim; ++i)
	  prod *= expS(Sb, x[i]) * EPn_ana(Ln, theta, r, n, x[i]);

    return prod;
}

double ESn_OverExpS_adapt(unsigned dim, const double *x, void *fparams){

	double Sb = ((struct params_s *)fparams)->Sb;

	double Ln = ((struct params_s *)fparams)->Ln;
	double theta = ((struct params_s *)fparams)->theta;
	double n = ((struct params_s *)fparams)->n;

    unsigned int i;

	double prod = 1.0;

    for (i = 0; i < dim; ++i){
    	prod *= expS(Sb, x[i]) * ESn_ana(Ln, theta, n, x[i]);
    }
    return prod;

}

double EDn_OverExpS_adapt(unsigned dim, const double *x, void *fparams){

	double Sb = ((struct params_s *)fparams)->Sb;

	double Ln = ((struct params_s *)fparams)->Ln;
	double lambda = ((struct params_s *)fparams)->lambda;
	double theta = ((struct params_s *)fparams)->theta;
	double r = ((struct params_s *)fparams)->r;
	double n = ((struct params_s *)fparams)->n;

    unsigned int i;

	double prod = 1.0;
    for (i = 0; i < dim; ++i){
    	prod *= expS(Sb, x[i]) * EDn_ana(Ln, lambda, theta, r, n, x[i]);
    }

    return prod;

}
