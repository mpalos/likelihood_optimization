/*
 * lnLfunctions.h
 *
 *  Created on: Jan 29, 2014
 *      Author: MAPF
 */

#ifndef LNLFUNCTIONS_H_
#define LNLFUNCTIONS_H_

/******************************************************************************************/
/* CUSTOM TYPES */
/******************************************************************************************/
typedef struct params_s {
	char model;
	char method;
	double Ls;
	double Ln;
	double theta;
	double lambda;
	double r;
	double n;
	double smean;
	double smax;
	double Pb;
	double Sb;
	double shapeBeta;
	double shapeMu;
	int nObs;
	int* PSobs;
	int* SSobs;
	int* DSobs;
	int* PNobs;
	int* SNobs;
	int* DNobs;
} params;

/******************************************************************************************/
/* FUNCTION PROTOTYPES */
/******************************************************************************************/
double log_NegativeBinomial(int k,double mu, double shape);
double EPs(double Ls, double r, double theta, int n);
double ESs(double Ls, double theta, double n);
double EDs(double Ls, double lambda, double theta, int n);
double EPn_ana (double Ln, double theta, double r, int n, double S);
double ESn_ana (double Ln, double theta, int n, double S);
double EDn_ana (double Ln, double lambda, double theta, double r, int n, double S);
double EPn_OverPhiS_adapt(unsigned dim, const double *x, void *params);
double ESn_OverPhiS_adapt(unsigned dim, const double *x, void *params);
double EDn_OverPhiS_adapt(unsigned dim, const double *x, void *params);
double EPn_OverExpS_adapt(unsigned dim, const double *x, void *fparams);
double ESn_OverExpS_adapt(unsigned dim, const double *x, void *fparams);
double EDn_OverExpS_adapt(unsigned dim, const double *x, void *fparams);

double lnL_MKS_PhiS_PhiMu(char model
		,double Ln, double Ls, int n
		,int obsPs,int obsPn, int obsSn, int obsSs, int obsDn, int obsDs //THESE PARAMETERS CAN BE IN ARRAY FORM
		,double Mutheta, double lambda, double r, double smean, double smax, double Pb, double Sb
		,double shapeBeta, double shapeMu
		,double cubTol, int cubMaxEval, double cubZero, double cubSInf);

/******************************************************************************************/

#endif /* LNLFUNCTIONS_H_ */
