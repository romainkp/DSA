/** 
 * dsa_model.h provides a wrapper for the various methods to fit generalized
 * linear models. It makes calls to various model fitting libraries such as
 * NAG or LINPACK depending on how things were compiled. It tries to sensibly
 * wrap errors so that they can be dealt with trivially in the code. 
 *
 * Author: James Bullard
 * Date: 10/24/05
 *
 */

#ifndef __HAS_DSA_MODEL__
#define __HAS_DSA_MODEL__

#include <math.h>

#ifdef __USE_NAG__
#include <nag.h>
#include <nagg02.h>
#endif

#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>
#include <Rmath.h>

#define DSA_MODEL_NO_ERROR                  0
#define DSA_MODEL_ERROR_AT_BOUNDARY         1
#define DSA_MODEL_ERROR_RANK_CHANGE         2
#define DSA_MODEL_ERROR_NO_CONVERGENCE      3
#define DSA_MODEL_ERROR_UNKNOWN             10

/** global constants specifying the debug level. **/
int R_DSA_DEBUG_LEVEL;   

#define R_DSA_DEBUG_TRACE   1
#define R_DSA_DEBUG_LOW     2 
#define R_DSA_DEBUG_WARN    3 

void dsa_xtx (double* x, int n, int p, double* z);
void dsa_transpose (double* x, int n, int p, double* xt);

/**  prints the level of debugging which we are currently using. */
void dsa_print_debug_level();

void dsa_fit_lr (double* X, int n, int tdx, int ip, double* Y, 
		 double tol, double* coefficients, 
		 double* residuals, double* effects, int* rank, 
		 int* pivot, double* qraux, double* work, double* dev,
		 double* weights, int max_iter, double* new_coefficients, 
		 double* WX, double* y_shifted, int transpose, int* error);

void dsa_fit_mlr (double* X, int n, int tdx, int ip, double* Y, 
		  double* weights, double tol, double* coefficients, 
		  double* residuals, double* effects, int* rank, 
		  int* pivot, double* qraux, double* work, double* rss,
		  int transpose, int* fail);

void run_lr(double* dataframeX, int samplesize, int TOTALTERMS, int IP, double* dataframeY, 
	    double TOL, double* PARA, double* RES, int* rank, double* WT,  
	    double* DEV, int* fail);

void run_ols(double* dataframeX, int samplesize, int TOTALTERMS, 
	     int IP, double* dataframeY, double* WT, 
	     double TOL, double* PARA, double* RES, int* rank, double* RSS, 
	     int* fail);

#endif /* __HAS_DSA_MODEL__ */
