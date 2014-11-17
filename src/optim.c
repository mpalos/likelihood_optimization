/*
 * TODO Include verbose levels
 * TODO Create new method parameters:
 * 	- bfgs: MAXITER, MAXSAME
 * 	- simplex: MAXITER,
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
# include <unistd.h>

#include <gsl/gsl_vector.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>

#include <time.h>
#include <sys/timeb.h>
#include <sys/types.h>

#include "lnLfunctions.h"

/******************************************************************************************/
/* FUNCTION PROTOTYPES */
/******************************************************************************************/
void readStartPoints(char *paramsfile, int startPointIn, char model);
void readParams(char *paramsfile, int nParamIn, char model, char method);
void readObsData(char *paramsfile);
void readConst(char *paramsfile, int nConstIn);
void freeAll();
/******************************************************************************************/

#define handle_error(msg) \
	do {perror(msg); exit(EXIT_FAILURE);} while (0)

// Global variables:
int nObs = 0;
int *PSobs = NULL;
int *SSobs = NULL;
int *DSobs = NULL;
int *PNobs = NULL;
int *SNobs = NULL;
int *DNobs = NULL;
double *x_init = NULL;
double step_size, tol, epsabs;
unsigned int maxSame = 10, maxIter = 1500;
double Ls = 3000, Ln = 7000, n = 24;
double cubTol, cubZero, cubSInf;
int cubMaxEval;
long nData = 0;
int startPoint = -1;

// Global variables:
double *x_params = NULL;
int nParam = -1;
int nRerun = 0;

double lnL_f(const gsl_vector * x, void *params) {

	double lnL = 0.0;
	int j;
	const char model = ((struct params_s *) params)->model;

	double theta;
	double lambda;
	double r;
	double smean;
	double smax;
	double Pb;
	double Sb;
	double shapeBeta;
	double shapeMu;

	// PARAMETERS
	if(model == 'A'){
		theta 			= exp(gsl_vector_get(x, 0));
		lambda 			= exp(gsl_vector_get(x, 1));
		r 				= exp(gsl_vector_get(x, 2));
		smean 			= gsl_vector_get(x, 3);
		smax 			= gsl_vector_get(x, 4);
		Pb				= 0;
		shapeBeta 		= exp(gsl_vector_get(x, 7));
		shapeMu 		= exp(gsl_vector_get(x, 8));
	}

	else if(model == 'B' || model == 'C'){
		theta 			= exp(gsl_vector_get(x, 0));
		lambda 			= exp(gsl_vector_get(x, 1));
		r 				= exp(gsl_vector_get(x, 2));
		smean 			= gsl_vector_get(x, 3);
		smax		= -cubZero;
		Pb			= 1 / (1 + exp(-1*gsl_vector_get (x, 5)));
		Sb			= exp(gsl_vector_get (x, 6));
		shapeBeta 		= exp(gsl_vector_get(x, 7));
		shapeMu 		= exp(gsl_vector_get(x, 8));
	}


	for (j = 0; j < ((struct params_s *) params)->nObs; j++) {
		lnL += lnL_MKS_PhiS_PhiMu(model,
				((struct params_s *) params)->Ln,
				((struct params_s *) params)->Ls,
				((struct params_s *) params)->n,
				((struct params_s *) params)->PSobs[j],
				((struct params_s *) params)->PNobs[j],
				((struct params_s *) params)->SNobs[j],
				((struct params_s *) params)->SSobs[j],
				((struct params_s *) params)->DNobs[j],
				((struct params_s *) params)->DSobs[j],
				theta, lambda, r, smean, smax, Pb, Sb, shapeBeta, shapeMu,
				cubTol, cubMaxEval, cubZero, cubSInf);
	}

	const double y = lnL * (-1);

	return y;
}

int lnL_f1(const gsl_vector * x, void *params, gsl_vector * f) {

	gsl_vector_set(f, 0, lnL_f(x, params));

	return GSL_SUCCESS;

}

int lnL_df(const gsl_vector *x, void *params, gsl_vector *df) {

	//HERE IS THE DEFERENTIAL USING FINITE DIFFERENCES (multifit_fdfsolver)
	gsl_multifit_function_fdf fdf;
	fdf.f = &lnL_f1;
	fdf.df = NULL;
	fdf.fdf = NULL;
	fdf.n = 1;
	fdf.p = 9;
	fdf.params = params;

	const char model = ((struct params_s *) params)->model;

	gsl_matrix *J = gsl_matrix_alloc(1, 9);
	J->size1 = 1;
	J->size2 = 9;

	gsl_vector *f = gsl_vector_alloc(1);

	gsl_vector_set(f, 0, lnL_f(x, params));

	int status = gsl_multifit_fdfsolver_dif_df(x, &fdf, f, J);

	if (status) {
		printf("ERROR #%d", status);
		return status;
	}

	gsl_vector_set(df, 0, gsl_matrix_get(J, 0, 0));
	gsl_vector_set(df, 1, gsl_matrix_get(J, 0, 1));
	gsl_vector_set(df, 2, gsl_matrix_get(J, 0, 2));
	gsl_vector_set(df, 3, gsl_matrix_get(J, 0, 3));
	gsl_vector_set(df, 4, gsl_matrix_get(J, 0, 4));
	if(model == 'B' || model == 'C'){
		gsl_vector_set(df, 5, gsl_matrix_get(J, 0, 5));
		gsl_vector_set(df, 6, gsl_matrix_get(J, 0, 6));
	}
	gsl_vector_set(df, 7, gsl_matrix_get(J, 0, 7));
	gsl_vector_set(df, 8, gsl_matrix_get(J, 0, 8));

	return GSL_SUCCESS;

}

int lnL_fdf(const gsl_vector * x, void *params, double * f, gsl_vector *df) {

	*f = lnL_f(x, params);
	lnL_df(x, params, df);
	return GSL_SUCCESS;
}

int optim_lnL_derivatives(const gsl_multimin_fdfminimizer_type *T, void *params,
		gsl_vector *x) {

	unsigned int iter = 0;
	unsigned int countSame = 0;
	int status;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *prev = gsl_vector_alloc(9);

	gsl_multimin_function_fdf my_func;
	my_func.n = 9;
	my_func.f = lnL_f;
	my_func.df = lnL_df;
	my_func.fdf = lnL_fdf;
	my_func.params = params;

	const char model = ((struct params_s *) params)->model;

	s = gsl_multimin_fdfminimizer_alloc(T, 9);

	gsl_multimin_fdfminimizer_set(s, &my_func, x, step_size, tol); //ORIGINAL

	//PRINT HEADER IF NECESSARY
//	printf ("startPoint;nData;iter;"
//			"theta;lambda;r;smean;smax;shapeBeta;shapeMu;"
//			"f;grad;status\n"
//	 );

	double mod_grad = 0.0;
	int n_dim = 2;

	//PRINT START POINT
	mod_grad = pow(gsl_vector_get(s->gradient, 0), n_dim)
			+ pow(gsl_vector_get(s->gradient, 1), n_dim)
			+ pow(gsl_vector_get(s->gradient, 2), n_dim)
			+ pow(gsl_vector_get(s->gradient, 3), n_dim)
			+ pow(gsl_vector_get(s->gradient, 4), n_dim)
			+ pow(gsl_vector_get(s->gradient, 7), n_dim)
			+ pow(gsl_vector_get(s->gradient, 8), n_dim);

	  if(model == 'B' || model == 'C'){
		  mod_grad += pow(gsl_vector_get(s->gradient, 5), n_dim)
					+ pow(gsl_vector_get(s->gradient, 6), n_dim);
	  }

	  printf("%d;",startPoint);
	  printf("%ld;",nData);
	  printf("%d;",iter);
	  printf("%.20f;",exp(gsl_vector_get(s->x, 0)));
	  printf("%.20f;",exp(gsl_vector_get(s->x, 1)));
	  printf("%.20f;",exp(gsl_vector_get(s->x, 2)));
	  printf("%.20f;",gsl_vector_get(s->x, 3));
	  printf("%.20f;",gsl_vector_get(s->x, 4));
	  if(model == 'B' || model == 'C'){
		  printf("%.20f;",1 / (1 + exp(-1*gsl_vector_get (s->x, 5))));
		  printf("%.20f;",exp(gsl_vector_get (s->x, 6)));
	  }
	  printf("%.20f;",exp(gsl_vector_get(s->x, 7)));
	  printf("%.20f;",exp(gsl_vector_get(s->x, 8)));
	  printf("%.20f;", s->f);
	  printf("%.20f;",pow(mod_grad, 1.0 / n_dim));
	  printf("-2\n");

	gsl_vector_memcpy(prev, x);

	do {
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);

		if (status) {
			// IF IS NECESSARY TO PRINT THE ERRORS ...
			//if (status==27){
			//	printf("ERROR #%d: iteration is not making progress towards solution\n",status);
			//}

			  printf("%d;",startPoint);
			  printf("%ld;",nData);
			  printf("%d;",iter);
			  printf("%.20f;",exp(gsl_vector_get(s->x, 0)));
			  printf("%.20f;",exp(gsl_vector_get(s->x, 1)));
			  printf("%.20f;",exp(gsl_vector_get(s->x, 2)));
			  printf("%.20f;",gsl_vector_get(s->x, 3));
			  printf("%.20f;",gsl_vector_get(s->x, 4));
			  if(model == 'B' || model == 'C'){
				  printf("%.20f;",1 / (1 + exp(-1*gsl_vector_get (s->x, 5))));
				  printf("%.20f;",exp(gsl_vector_get (s->x, 6)));
			  }
			  printf("%.20f;",exp(gsl_vector_get(s->x, 7)));
			  printf("%.20f;",exp(gsl_vector_get(s->x, 8)));
			  printf("%.20f;", s->f);
			  printf("%.20f;",pow(mod_grad, 1.0 / n_dim));
			  printf("%d\n",status);

			break;
		}

		status = gsl_multimin_test_gradient(s->gradient, epsabs);
		// IF IS NECESSARY TO PRINT THE MIN ...
		//if (status == GSL_SUCCESS){
		//	printf ("Minimum found at:\n");
		//}

		mod_grad = pow(gsl_vector_get(s->gradient, 0), n_dim)
				+ pow(gsl_vector_get(s->gradient, 1), n_dim)
				+ pow(gsl_vector_get(s->gradient, 2), n_dim)
				+ pow(gsl_vector_get(s->gradient, 3), n_dim)
				+ pow(gsl_vector_get(s->gradient, 4), n_dim)
				+ pow(gsl_vector_get(s->gradient, 7), n_dim)
				+ pow(gsl_vector_get(s->gradient, 8), n_dim);

		  if(model == 'B' || model == 'C'){
			  mod_grad += pow(gsl_vector_get(s->gradient, 5), n_dim)
						+ pow(gsl_vector_get(s->gradient, 6), n_dim);
		  }

		  printf("%d;",startPoint);
		  printf("%ld;",nData);
		  printf("%d;",iter);
		  printf("%.20f;",exp(gsl_vector_get(s->x, 0)));
		  printf("%.20f;",exp(gsl_vector_get(s->x, 1)));
		  printf("%.20f;",exp(gsl_vector_get(s->x, 2)));
		  printf("%.20f;",gsl_vector_get(s->x, 3));
		  printf("%.20f;",gsl_vector_get(s->x, 4));
		  if(model == 'B' || model == 'C'){
			  printf("%.20f;",1 / (1 + exp(-1*gsl_vector_get (s->x, 5))));
			  printf("%.20f;",exp(gsl_vector_get (s->x, 6)));
		  }
		  printf("%.20f;",exp(gsl_vector_get(s->x, 7)));
		  printf("%.20f;",exp(gsl_vector_get(s->x, 8)));
		  printf("%.20f;", s->f);
		  printf("%.20f;",pow(mod_grad, 1.0 / n_dim));
		  printf("%d\n",status);

		  //TEST IF THE ITERATION STUCKS ON SAME ANSWER AFTER maxSame TIMES
		  if(gsl_vector_equal(prev,s->x)){ countSame++; }
		  else { countSame =0;}

		  if(countSame >= maxSame){
			  printf("%d;",startPoint);
			  printf("%ld;",nData);
			  printf("%d;",iter);
			  printf("%.20f;",exp(gsl_vector_get(s->x, 0)));
			  printf("%.20f;",exp(gsl_vector_get(s->x, 1)));
			  printf("%.20f;",exp(gsl_vector_get(s->x, 2)));
			  printf("%.20f;",gsl_vector_get(s->x, 3));
			  printf("%.20f;",gsl_vector_get(s->x, 4));
			  if(model == 'B' || model == 'C'){
				  printf("%.20f;",1 / (1 + exp(-1*gsl_vector_get (s->x, 5))));
				  printf("%.20f;",exp(gsl_vector_get (s->x, 6)));
			  }
			  printf("%.20f;",exp(gsl_vector_get(s->x, 7)));
			  printf("%.20f;",exp(gsl_vector_get(s->x, 8)));
			  printf("%.20f;", s->f);
			  printf("%.20f;",pow(mod_grad, 1.0 / n_dim));
			  printf("27\n");

			  break;
			}

		gsl_vector_memcpy(prev, s->x);

	} while (status == GSL_CONTINUE && iter < maxIter);

	//VERIFY IF AFTER maxIter ITERATIONS IT COULDN'T FIND SOLUTION
	if (iter >= maxIter) {
		  printf("%d;",startPoint);
		  printf("%ld;",nData);
		  printf("%d;",iter);
		  printf("%.20f;",exp(gsl_vector_get(s->x, 0)));
		  printf("%.20f;",exp(gsl_vector_get(s->x, 1)));
		  printf("%.20f;",exp(gsl_vector_get(s->x, 2)));
		  printf("%.20f;",gsl_vector_get(s->x, 3));
		  printf("%.20f;",gsl_vector_get(s->x, 4));
		  if(model == 'B' || model == 'C'){
			  printf("%.20f;",1 / (1 + exp(-1*gsl_vector_get (s->x, 5))));
			  printf("%.20f;",exp(gsl_vector_get (s->x, 6)));
		  }
		  printf("%.20f;",exp(gsl_vector_get(s->x, 7)));
		  printf("%.20f;",exp(gsl_vector_get(s->x, 8)));
		  printf("%.20f;", s->f);
		  printf("%.20f;",pow(mod_grad, 1.0 / n_dim));
		  printf("27\n");
	}

	gsl_multimin_fdfminimizer_free(s);

	return 0;

}

int optim_lnL_simplex(const gsl_multimin_fminimizer_type *T, void *params, gsl_vector *x){

	  unsigned int n = 9;

	  gsl_multimin_fminimizer *s = NULL;
	  gsl_vector *ss;
	  gsl_multimin_function minex_func = {&lnL_f, n, params};

	  unsigned int iter = 0;
	  int status;
	  double size = 0.0;
	  const char model = ((struct params_s *) params)->model;

	  ss = gsl_vector_alloc (9);

	  //Use x_params variable
	  gsl_vector_set(ss, 0, log(1 + (x_params[0] / exp(gsl_vector_get(x,0)))));
	  gsl_vector_set(ss, 1, log(1 + (x_params[1] / exp(gsl_vector_get(x,1)))));
	  gsl_vector_set(ss, 2, log(1 + (x_params[2]   / exp(gsl_vector_get(x,2)))));
	  gsl_vector_set(ss, 3, x_params[3]);
	  gsl_vector_set(ss, 4, x_params[4]);

	  if (model == 'B' || model == 'C'){
		  gsl_vector_set(x,5,1 / (1 + exp(-1*gsl_vector_get (x, 5)))); // Back to the Pb original value
	  	  gsl_vector_set(ss, 5, log((gsl_vector_get(x,5) + x_params[5])/(1-(gsl_vector_get(x,5) + x_params[5]))) - log(gsl_vector_get(x,5)/(1-gsl_vector_get(x,5))));
	  	  gsl_vector_set(x,5,log(gsl_vector_get(x,5)/(1-gsl_vector_get(x,5)))); // Now the Pb transformation

	  	  gsl_vector_set(ss, 6, x_params[6]);
	  }

	  gsl_vector_set(ss, 7, log(1 + (x_params[7]  / exp(gsl_vector_get(x,7)))));
	  gsl_vector_set(ss, 8, log(1 + (x_params[8]  / exp(gsl_vector_get(x,8)))));


	  s = gsl_multimin_fminimizer_alloc (T, 9);
	  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	  //PRINT START POINT
	  size = gsl_multimin_fminimizer_size (s);
	  printf("%d;",startPoint);
	  printf("%ld;",nData);
	  printf("%d;",iter);
	  printf("%.20f;",exp(gsl_vector_get (s->x, 0)));
	  printf("%.20f;",exp(gsl_vector_get (s->x, 1)));
	  printf("%.20f;",exp(gsl_vector_get (s->x, 2)));
	  printf("%.20f;",gsl_vector_get (s->x, 3));
	  printf("%.20f;",gsl_vector_get (s->x, 4));
	  if(model == 'B' || model == 'C'){
		  printf("%.20f;",1 / (1 + exp(-1*gsl_vector_get (s->x, 5))));
		  printf("%.20f;",exp(gsl_vector_get (s->x, 6)));
	  }
	  printf("%.20f;",exp(gsl_vector_get (s->x, 7)));
	  printf("%.20f;",exp(gsl_vector_get (s->x, 8)));
	  printf("%.20f;",lnL_f(x,params));
	  printf("%.20f;",size);
	  printf("-2\n");

	  gsl_set_error_handler_off();

	  	unsigned int var = 0;
	  	unsigned int iterPerRun = 0;
		for (var = 0; var <= nRerun; var++) {
		  	iterPerRun = 0; //set iterPerRun
			do
			  {
				iter++;
				iterPerRun++;
				status = gsl_multimin_fminimizer_iterate(s);

				if (status){ //print iteration for error status
					if(var!=nRerun){
						status = status + (100 * (nRerun - var)); // To indicate intermediate execs
					}

						printf("%d;",startPoint);
						printf("%ld;",nData);
						printf("%d;",iter);
						printf("%.20f;",exp(gsl_vector_get (s->x, 0)));
						printf("%.20f;",exp(gsl_vector_get (s->x, 1)));
						printf("%.20f;",exp(gsl_vector_get (s->x, 2)));
						printf("%.20f;",gsl_vector_get (s->x, 3));
						printf("%.20f;",gsl_vector_get (s->x, 4));
						  if(model == 'B' || model == 'C'){
							  printf("%.20f;",1 / (1 + exp(-1*gsl_vector_get (s->x, 5))));
							  printf("%.20f;",exp(gsl_vector_get (s->x, 6)));
						  }
						printf("%.20f;",exp(gsl_vector_get (s->x, 7)));
						printf("%.20f;",exp(gsl_vector_get (s->x, 8)));
						printf("%.20f;",s->fval);
						printf("%.20f;",size);
						printf("%d\n",status);

						if(var!=nRerun){ //Reset rerun
							gsl_multimin_fminimizer_set (s, &minex_func, s->x, ss); //set the simplex again!!!
							//iterPerRun = 0; //set iterPerRun
						}

						break;
				}

				size = gsl_multimin_fminimizer_size (s);
				status = gsl_multimin_test_size (size, epsabs);

				if(status == 0 && var!=nRerun){ //Force to status if do not reach the # of re-executions
					status = status + (100 * (nRerun - var)); // To indicate intermediate execs
					gsl_multimin_fminimizer_set (s, &minex_func, s->x, ss); //set the simplex again!!!
					//iterPerRun = 0; //set iterPerRun
				}

				printf("%d;",startPoint);
				printf("%ld;",nData);
				printf("%d;",iter);
				printf("%.20f;",exp(gsl_vector_get (s->x, 0)));
				printf("%.20f;",exp(gsl_vector_get (s->x, 1)));
				printf("%.20f;",exp(gsl_vector_get (s->x, 2)));
				printf("%.20f;",gsl_vector_get (s->x, 3));
				printf("%.20f;",gsl_vector_get (s->x, 4));
				  if(model == 'B' || model == 'C'){
					  printf("%.20f;",1 / (1 + exp(-1*gsl_vector_get (s->x, 5))));
					  printf("%.20f;",exp(gsl_vector_get (s->x, 6)));
				  }
				printf("%.20f;",exp(gsl_vector_get (s->x, 7)));
				printf("%.20f;",exp(gsl_vector_get (s->x, 8)));
				printf("%.20f;",s->fval);
				printf("%.20f;",size);
				printf("%d\n",status);


			  }
			while (status == GSL_CONTINUE && iterPerRun < maxIter);

			//PRINT STATUS IF IT HAVEN'T FOUND THE SOLUTION AFTER maxIter ITERATIONS PER RUN
			status = 27 + (100 * (nRerun - var));
			if(iterPerRun >= maxIter){
				printf("%d;",startPoint);
				printf("%ld;",nData);
				printf("%d;",iter);
				printf("%.20f;",exp(gsl_vector_get (s->x, 0)));
				printf("%.20f;",exp(gsl_vector_get (s->x, 1)));
				printf("%.20f;",exp(gsl_vector_get (s->x, 2)));
				printf("%.20f;",gsl_vector_get (s->x, 3));
				printf("%.20f;",gsl_vector_get (s->x, 4));
				  if(model == 'B' || model == 'C'){
					  printf("%.20f;",1 / (1 + exp(-1*gsl_vector_get (s->x, 5))));
					  printf("%.20f;",exp(gsl_vector_get (s->x, 6)));
				  }
				printf("%.20f;",exp(gsl_vector_get (s->x, 7)));
				printf("%.20f;",exp(gsl_vector_get (s->x, 8)));
				printf("%.20f;",s->fval);
				printf("%.20f;",size);
				printf("%d\n",status);
			}
		}

	  gsl_vector_free(ss);
	  gsl_multimin_fminimizer_free (s);

	  return status;
}

int main(int argc, char **argv) {

	char *obsFile, *startFile, *paramFile, *constFile;

	// Parameters ID:
	int nStart = -1;
	int nParam = -1;
	int nConst = -1;

	//Default Parameters (cubature)
	cubTol = 1e-5;
	cubMaxEval = 1e4;
	cubZero = 1e-4;
	cubSInf = 500;

	char method = 'B';
	char model = 'A';

	//NEW: Use getopt()
	//TODO make some options mandatory
	int flagArgs = 0;
	int opt;
	while ((opt = getopt(argc, argv, "o:p:s:c:t:m:z:i:d:h:")) != -1) {
		switch (opt) {
		//File Parameters
		case 'd':
			if(strcmp(optarg,"model_A") == 0){ model = 'A';	}
			else if(strcmp(optarg,"model_B") == 0){	model = 'B'; }
			else if(strcmp(optarg,"model_C") == 0){	model = 'C'; }
			// TODO error if not provide the correct model
			break;
		case 'h':
			if(strcmp(optarg,"bfgs") == 0){ method = 'B';	}
			else if(strcmp(optarg,"simplex") == 0){	method = 'S'; }
			else if(strcmp(optarg,"conj_pr") == 0){ method = 'P';}
			else if(strcmp(optarg,"conj_fr") == 0){ method = 'F';}
			else if(strcmp(optarg,"steep_desc") == 0){ method = 'D';}
			// TODO Include MCMC
			// TODO error if not provide the correct method
			break;
		case 'o': //
			obsFile = optarg;
			break;
		case 's': //
			startFile = optarg;

			if (optind < argc && *argv[optind] != '-') {
				nStart = atoi(argv[optind]);
				optind++;
			} else {
				fprintf(stderr,
						"\n-s option requires TWO arguments <startPointFile> "
								"<startPointID>\n\n");
				//usage();
				flagArgs = 1;
			}

			break;
		case 'p': //
			paramFile = optarg;

			if (optind < argc && *argv[optind] != '-') {
				nParam = atoi(argv[optind]);
				optind++;
			} else {
				fprintf(stderr,
						"\n-p option requires TWO arguments <paramFile> "
								"<paramID>\n\n");
				//usage();
				flagArgs = 1;
			}

			break;
		case 'c': //
			constFile = optarg;

			if (optind < argc && *argv[optind] != '-') {
				nConst = atoi(argv[optind]);
				optind++;
			} else {
				fprintf(stderr,
						"\n-c option requires TWO arguments <constFile> "
								"<constID>\n\n");
				//usage();
				flagArgs = 1;
			}

			break;
			//Cubature parameters:
		case 't': //
			cubTol = atof(optarg);
			break;
		case 'm': //
			cubMaxEval = atoi(optarg);
			break;
		case 'z': //
			cubZero = atof(optarg);
			break;
		case 'i': //
			cubSInf = atof(optarg);
			break;
		default:
			fprintf(stderr,
					"Usage: %s -d model(model_A, model_B, model_C) -h method(bfgs, simplex) -o obsFile -s startPointFile startPointID -p paramFile paramID -c constFile constID"
							" [-t tol(cubature)] [-m maxEval (cubature)] [-z zero(cubature)] [-i SInfVal(cubature)]\n",
					argv[0]);
			return EXIT_FAILURE;
		}
	}

	//Verify if generate any errors reading the arguments
	if (flagArgs == 1) {
		fprintf(stderr,
				"Usage: %s -d model(model_A, model_B, model_C) -h method(bfgs, simplex) -o obsFile -s startPointFile startPointID -p paramFile paramID -c constFile constID"
						" [-t tol(cubature)] [-m maxEval (cubature)] [-z zero(cubature)] [-i SInfVal(cubature)]\n",
				argv[0]);
		return EXIT_FAILURE;
	}

	//READ FILES AS ARGUMENT
//    if (argc <= 7) {
//	  fprintf(stderr, "Usage: %s obsFile startPointFile startPointID paramFile paramID constFile constID"
//			  "[tol(cubature)] [maxEval (cubature)] [zero(cubature)] [SInfVal(cubature)]\n",
//		  argv[0]);
//	  return EXIT_FAILURE;
//    }
//
//    obsFile = argv[1];
//    startFile = argv[2];
//    paramFile = argv[4];
//    constFile = argv[6];

	readObsData(obsFile);
	readStartPoints(startFile, nStart,model);
	readParams(paramFile, nParam,model,method);

	readConst(constFile, nConst);

	//cubature parameters:
//	cubTol = argc > 8 ? atof(argv[8]) : 1e-5;
//	cubMaxEval = argc > 9 ? atoi(argv[9]) : 1e4;
//	cubZero = argc > 10 ? atof(argv[10]) : 1e-4;
//	cubSInf = argc > 11 ? atof(argv[11]) : 500;

	params par;
	par.Ls = Ls;
	par.Ln = Ln;
	par.n = n;
	par.nObs = nObs;
	par.PSobs = PSobs;
	par.SSobs = SSobs;
	par.DSobs = DSobs;
	par.PNobs = PNobs;
	par.SNobs = SNobs;
	par.DNobs = DNobs;
	par.method = method;
	par.model = model;

	gsl_vector_view x = gsl_vector_view_array(x_init, 9);

	// Define optimization algorithm/ run the optimizer

	const gsl_multimin_fdfminimizer_type *Tfdfmin;
	const gsl_multimin_fminimizer_type *Tfmin;

	if(method == 'B'){
		Tfdfmin = gsl_multimin_fdfminimizer_vector_bfgs2;
		optim_lnL_derivatives(Tfdfmin, &par, &x.vector);
	}

	else if (method == 'S'){
		Tfmin = gsl_multimin_fminimizer_nmsimplex2;
		optim_lnL_simplex(Tfmin, &par, &x.vector);
	}

	else if (method == 'P'){
		Tfdfmin = gsl_multimin_fdfminimizer_conjugate_pr;
		optim_lnL_derivatives(Tfdfmin, &par, &x.vector);
	}

	else if (method == 'F'){
		Tfdfmin = gsl_multimin_fdfminimizer_conjugate_fr;
		optim_lnL_derivatives(Tfdfmin, &par, &x.vector);
	}

	else if (method == 'D'){
		Tfdfmin = gsl_multimin_fdfminimizer_steepest_descent;
		optim_lnL_derivatives(Tfdfmin, &par, &x.vector);
	}


	return EXIT_SUCCESS;
}

void readObsData(char *paramsfile) {

	FILE *filein, *bufin;

	char *linebuf = NULL;
	size_t buflen = 0;

	int d; 		//temporary store data
	int pp; 	//temporary store data
	int line;
	int j;

	// Try to find file
	if (paramsfile == NULL ) {
		filein = stdin;
	} else {
		filein = fopen(paramsfile, "r");
		if (filein == NULL ) {
			handle_error("obs - fopen");
		}
	}

	// Read the file
	line = 0;
	// find # of lines
	while ((pp = getline(&linebuf, &buflen, filein)) != -1) {
		// Ignore empty line and line beginning with '#'
		if (pp <= 1 || linebuf[0] == '#') {
			continue;
		}
		line++;
	}

	fclose(filein);

	nObs = line;

	//Observation memory allocation
	PSobs = malloc(sizeof(int) * nObs);
	if (PSobs == NULL )
		handle_error("malloc: PSobs");
	SSobs = malloc(sizeof(int) * nObs);
	if (SSobs == NULL )
		handle_error("malloc: SSobs");
	DSobs = malloc(sizeof(int) * nObs);
	if (DSobs == NULL )
		handle_error("malloc: DSobs");
	PNobs = malloc(sizeof(int) * nObs);
	if (PNobs == NULL )
		handle_error("malloc: PNobs");
	SNobs = malloc(sizeof(int) * nObs);
	if (SNobs == NULL )
		handle_error("malloc: SNobs");
	DNobs = malloc(sizeof(int) * nObs);
	if (DNobs == NULL )
		handle_error("malloc: DNobs");

	// Read again the file to load the observation values
	// Try to find file
	if (paramsfile == NULL ) {
		filein = stdin;
	} else {
		filein = fopen(paramsfile, "r");
		if (filein == NULL ) {
			handle_error("obs - fopen");
		}
	}

	// Read the obs file
	line = 0;
	long temp[7];
	while ((pp = getline(&linebuf, &buflen, filein)) != -1) {
		// Ignore empty line and line beginning with '#'
		if (pp <= 1 || linebuf[0] == '#') {
			continue;
		}
		bufin = fmemopen(linebuf, buflen, "r");
		if (bufin == NULL )
			handle_error("fmemopen");
		for (j = 0; j < 7; j++) {
			if (fscanf(bufin, "%d", &d) != 1) {
				fprintf(stderr, "params file format error - Obs - line: %d\n",
						line);
				freeAll();
				exit(EXIT_FAILURE);
			}
			temp[j] = d;
		}
		fclose(bufin);

		nData = temp[0];
		PSobs[line] = temp[1];
		SSobs[line] = temp[2];
		DSobs[line] = temp[3];
		PNobs[line] = temp[4];
		SNobs[line] = temp[5];
		DNobs[line] = temp[6];

		line++;

	}
	fclose(filein);

}

void readStartPoints(char *paramsfile, int startPointIn, char model) {

	FILE *filein, *bufin;

	char *linebuf = NULL;
	size_t buflen = 0;
	int nStartModel = -1;

	if (model == 'A') {nStartModel = 8;}
	else if(model == 'B' || model == 'C') {nStartModel = 10;}

	double d; 	//temporary store data
	int pp; 	//temporary store data
	int j;

	//startpoint memory allocation
	x_init = malloc(sizeof(double) * 9); if (x_init == NULL) {handle_error("malloc: x_init");}

	// Try to find file
	if (paramsfile == NULL ) {
		filein = stdin;
	} else {
		filein = fopen(paramsfile, "r");
		if (filein == NULL ) {
			handle_error("startPoint - fopen");
		}
	}

	double temp[nStartModel];
	while ((pp = getline(&linebuf, &buflen, filein)) != -1) {
		// Ignore empty line and line beginning with '#'
		if (pp <= 1 || linebuf[0] == '#') {
			continue;
		}
		bufin = fmemopen(linebuf, buflen, "r");
		if (bufin == NULL )
			handle_error("fmemopen");
		for (j = 0; j < nStartModel; j++) {
			if (fscanf(bufin, "%lf", &d) != 1) {
				fprintf(stderr, "params file format error - x_init\n");
				freeAll();
				exit(EXIT_FAILURE);
			}
			temp[j] = d;
		}

		fclose(bufin);

		if (temp[0] == startPointIn) {
			startPoint = temp[0];

			if(model == 'A'){
				x_init[0] = log(temp[1]);				// theta
				x_init[1] = log(temp[2]);				// lambda
				x_init[2] = log(temp[3]);				// r
				x_init[3] = temp[4];					// smean
				x_init[4] = temp[5];					// smax
				x_init[7] = log(temp[6]);				// shapeBeta
				x_init[8] = log(temp[7]);				// shapeMu
			}

			else if(model == 'B' || model == 'C'){
				x_init[0] = log(temp[1]); 				// theta
				x_init[1] = log(temp[2]);				// lambda
				x_init[2] = log(temp[3]);				// r
				x_init[3] = temp[4];					// smean
				x_init[4] = temp[5];					// smax
				x_init[5] = log(temp[6]/(1-temp[6]));	// Pb
				x_init[6] = log(temp[7]);   			// Sb
				x_init[7] = log(temp[8]);				// shapeBeta
				x_init[8] = log(temp[9]);				// shapeMu
			}
			break;
		}

	}
	fclose(filein);

	//IF COULDN'T FIND STARTPOINT
	if (startPoint == -1) {
		fprintf(stderr, "couldn't find startpoint - %d", startPointIn);
		freeAll();
		exit(EXIT_FAILURE);
	}

}

void readParams(char *paramsfile, int nParamIn, char model, char method) {

	FILE *filein, *bufin;

	char *linebuf = NULL;
	size_t buflen = 0;
	int nParamModelMethod = -1;

	double d; 	//temporary store data
	int pp; 	//temporary store data
	int j;

	if(method == 'B' || method == 'P' || method == 'F' || method == 'D'){
		nParamModelMethod = 6;
	}

	else if(method == 'S'){
		if(model == 'A'){
			nParamModelMethod = 11;
		}
		else if(model == 'B' || model == 'C'){
			nParamModelMethod = 13;
		}
		//params memory allocation (for simplex method)
		x_params = malloc(sizeof(double) * 10); if (x_init == NULL) handle_error("malloc: x_params");
	}

	// Try to find file
	if (paramsfile == NULL ) {
		filein = stdin;
	} else {
		filein = fopen(paramsfile, "r");
		if (filein == NULL ) {
			handle_error("paramsfile - fopen");
		}
	}

	int nParam = -1;

	double temp[nParamModelMethod];
	while ((pp = getline(&linebuf, &buflen, filein)) != -1) {
		// Ignore empty line and line beginning with '#'
		if (pp <= 1 || linebuf[0] == '#') {
			continue;
		}
		bufin = fmemopen(linebuf, buflen, "r");
		if (bufin == NULL )
			handle_error("fmemopen");
		for (j = 0; j < nParamModelMethod; j++) {
			if (fscanf(bufin, "%lf", &d) != 1) {
				fprintf(stderr, "params file format error\n");
				freeAll();
				exit(EXIT_FAILURE);
			}
			temp[j] = d;
		}

		fclose(bufin);

		if (temp[0] == nParamIn) {
			nParam = temp[0];

			if(method == 'B' || method == 'P' || method == 'F' || method == 'D'){ // from BFGS method
				epsabs = temp[1];
				step_size = temp[2];
				tol = temp[3];
				maxIter = temp[4];
				maxSame = temp[5];
			}

			else if(method == 'S'){ // from Simplex method
				maxIter = temp[1];
				nRerun = temp[2];
				epsabs = temp[3];

				if(model == 'A'){
					x_params[0] = temp[4];	//theta
					x_params[1] = temp[5];	// lambda
					x_params[2] = temp[6];	// r
					x_params[3] = temp[7];	// smean
					x_params[4] = temp[8];	// smax
					x_params[7] = temp[9];	// shapeBeta
					x_params[8] = temp[10];	// shapeMu
				}

				else if(model == 'B' || model == 'C'){
					x_params[0] = temp[4];	// theta
					x_params[1] = temp[5];	// lambda
					x_params[2] = temp[6];	// r
					x_params[3] = temp[7];	// smean
					x_params[4] = temp[8];	// smax
					x_params[5] = temp[9];	// Pb
					x_params[6] = temp[10];	// Sb
					x_params[7] = temp[11];	// shapeBeta
					x_params[8] = temp[12];	// shapeMu
				}

			}

			break;
		}

	}
	fclose(filein);

	//IF COULDN'T FIND NPARAM
	if (nParam == -1) {
		fprintf(stderr, "couldn't find nParam - %d", nParamIn);
		freeAll();
		exit(EXIT_FAILURE);
	}

}

void readConst(char *paramsfile, int nConstIn) {

	FILE *filein, *bufin;

	char *linebuf = NULL;
	size_t buflen = 0;

	double d; 	//temporary store data
	int pp; 	//temporary store data
	int j;

	// Try to find file
	if (paramsfile == NULL ) {
		filein = stdin;
	} else {
		filein = fopen(paramsfile, "r");
		if (filein == NULL ) {
			handle_error("constFile - fopen");
		}
	}

	int nParam = -1;

	double temp[4];
	while ((pp = getline(&linebuf, &buflen, filein)) != -1) {
		// Ignore empty line and line beginning with '#'
		if (pp <= 1 || linebuf[0] == '#') {
			continue;
		}
		bufin = fmemopen(linebuf, buflen, "r");
		if (bufin == NULL )
			handle_error("fmemopen");
		for (j = 0; j < 4; j++) {
			if (fscanf(bufin, "%lf", &d) != 1) {
				fprintf(stderr, "params file format error\n");
				freeAll();
				exit(EXIT_FAILURE);
			}
			temp[j] = d;
		}

		fclose(bufin);

		if (temp[0] == nConstIn) {
			nParam = temp[0];
			Ls = temp[1];
			Ln = temp[2];
			n = temp[3];
			break;
		}

	}
	fclose(filein);

	//IF COULDN'T FIND NPARAM
	if (nParam == -1) {
		fprintf(stderr, "couldn't find nParam - %d", nConstIn);
		freeAll();
		exit(EXIT_FAILURE);
	}

}

void freeAll() {
	if (PSobs)
		free(PSobs);
	if (SSobs)
		free(SSobs);
	if (DSobs)
		free(DSobs);
	if (PNobs)
		free(PNobs);
	if (SNobs)
		free(SNobs);
	if (DNobs)
		free(DNobs);
}
