/*
 * likelihood.c
 *
 *  Created on: Jan 29, 2014
 *      Author: MAPF
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cubature.h"
#include <ctype.h>

#include "lnLfunctions.h"

#if defined(PCUBATURE)
#  define cubature pcubature
#else
#  define cubature hcubature
#endif

int *f_integrand = NULL;
unsigned integrand_functions = 0;

int integrands_phi(unsigned dim, const double *x, void *data_,
	   unsigned fdim, double *retval)
{
     double val;
     unsigned j;

     for (j = 0; j < fdim; ++j) {
     switch (f_integrand[j]) {
		case 0: // EPn_OverPhiS_adapt
			val = EPn_OverPhiS_adapt(dim, x, data_);
			break;
		case 1: // ESn_OverPhiS_adapt
			val = ESn_OverPhiS_adapt(dim, x, data_);
			break;
		case 2: // EDn_OverPhiS_adapt
			val = EDn_OverPhiS_adapt(dim, x, data_);
			break;
		default:
	      fprintf(stderr, "unknown integrand %d\n", f_integrand[j]);
	      exit(EXIT_FAILURE);
     }
     retval[j] = val;
     }
     return 0;
}

int integrands_exp(unsigned dim, const double *x, void *data_,
	   unsigned fdim, double *retval)
{
     double val;
     unsigned j;

     for (j = 0; j < fdim; ++j) {
     switch (f_integrand[j]) {
		case 0: // EPn_OverExpS_adapt
			val = EPn_OverExpS_adapt(dim, x, data_);
			break;
		case 1: // ESn_OverExpS_adapt
			val = ESn_OverExpS_adapt(dim, x, data_);
			break;
		case 2: // EDn_OverExpS_adapt
			val = EDn_OverExpS_adapt(dim, x, data_);
			break;
		default:
	      fprintf(stderr, "unknown integrand %d\n", f_integrand[j]);
	      exit(EXIT_FAILURE);
     }
     retval[j] = val;
     }
     return 0;
}

// Obs: Ls, Ln and n are constants not parameters
double lnL_MKS_PhiS_PhiMu(char model
		,double Ln,double Ls, int n
		,int obsPs, int obsPn, int obsSn, int obsSs, int obsDn, int obsDs //THESE PARAMETERS CAN BE IN ARRAY FORM
		,double Mutheta, double lambda, double r, double smean, double smax, double Pb, double Sb, double shapeBeta, double shapeMu
		,double cubTol, int cubMaxEval, double cubZero, double cubSInf
		){

    double *xmin, *xmax;
    double *val, *valSum;
    double tol, *err;
    unsigned i, dim, maxEval;

	dim = 1;
	tol = cubTol;
	maxEval = cubMaxEval;

	integrand_functions = 3;
	f_integrand = (int *) malloc(sizeof(int) * integrand_functions);
	for(i=0; i<integrand_functions;i++){
		f_integrand[i]=i;
	}
	/*****************************/

	 err = (double *) malloc(sizeof(double) * integrand_functions);

     xmin = (double *) malloc(dim * sizeof(double));
     xmax = (double *) malloc(dim * sizeof(double));

     val = (double *) malloc(sizeof(double) * integrand_functions);
     valSum = (double *) malloc(sizeof(double) * integrand_functions);

     params fparams;

     fparams.Ls = Ls;
     fparams.Ln = Ln;
     fparams.theta = Mutheta;
     fparams.lambda = lambda;
     fparams.r = r;
     fparams.n = n;
     fparams.smean = smean;
     fparams.smax = smax;
     fparams.Pb = Pb;
     fparams.Sb = Sb;
     fparams.shapeBeta = shapeBeta;
     fparams.shapeMu = shapeMu;

	 /////////////////////////////////////////////////

	 //NEGATIVE PART - Integral over Phi
	 for (i = 0; i < dim; ++i) {
		 xmin[i] = -cubSInf;
		 if(smax > cubZero){ //In this case, the calculation will be performed 2 times
			 xmax[i] = -cubZero;
		 }
		 else
			 xmax[i] = smax;
     }
	 ///////////////////////

     cubature(integrand_functions, integrands_phi, &fparams,
	      dim, xmin, xmax,
	      maxEval, 0, tol, ERROR_INDIVIDUAL, val, err);

     //sum of val:
     for(i=0;i<integrand_functions;i++){
    	 valSum[i] = (1 - Pb) * val[i];
     }

     //POSITIVE PART
     if(model == 'A'){
     // model_A : Integral over Phi
         if(smax > cubZero){
        	 for (i = 0; i < dim; ++i) {
            	 xmin[i] = cubZero;
            	 xmax[i] = smax;
        	 }

             cubature(integrand_functions, integrands_phi, &fparams,
        	      dim, xmin, xmax,
        	      maxEval, 0, tol, ERROR_INDIVIDUAL, val, err);

             //sum of val:
             for(i=0;i<integrand_functions;i++){
            	 valSum[i] += val[i];
             }

         }
     }

     else if (model == 'B'){
         //model_B : Constant
         // Sum of val
         valSum[0] += Pb * EPn_ana(Ln, Mutheta, r, n, Sb);
         valSum[1] += Pb * ESn_ana(Ln,Mutheta,n,Sb);
         valSum[2] += Pb * EDn_ana(Ln,lambda,Mutheta,r,n,Sb);
     }
     else if (model == 'C'){
         //POSITIVE PART - Integral over Exp
    	 for (i = 0; i < dim; ++i) {
    		 xmin[i] = cubZero;
    		 xmax[i] = cubSInf;
    	 }

    	  cubature(integrand_functions, integrands_exp, &fparams,
    		  dim, xmin, xmax,
    		  maxEval, 0, tol, ERROR_INDIVIDUAL, val, err);

    	  //sum of val:
    	  for(i=0;i<integrand_functions;i++){
    		 valSum[i] += Pb * val[i];
    	  }
     }

     //NOW THE LIKELIHOOD CALCULATION
     double lnL = 0.0;
     lnL = log_NegativeBinomial(obsPs,EPs(Ls,r,Mutheta,n),shapeMu)
		+ log_NegativeBinomial(obsSs,ESs(Ls,Mutheta,n),shapeMu)
		+ log_NegativeBinomial(obsDs,EDs(Ls,lambda,Mutheta,n),shapeMu)
		+ log_NegativeBinomial(obsPn,valSum[0],shapeMu)
		+ log_NegativeBinomial(obsSn,valSum[1],shapeMu)
		+ log_NegativeBinomial(obsDn,valSum[2],shapeMu)
		;

     return(lnL);

}
