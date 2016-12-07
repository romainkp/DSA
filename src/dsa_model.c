

#include "dsa_model.h"

/* #include "dsa.h" */

/* taken right from family.c in the R code. */
double y_log_y(double y, double mu)
{
    return (y) ? (y * log(y/mu)) : 0;
}

/** functions independent of NAG **/
/*
void dsa_xtx (double* x, int n, int p, double* z) {
  char* transa = "T";
  char* transb = "N";
  double one   = 1.0; 
  double zero  = 0.0;
  
  F77_CALL(dgemm)(transa, transb, &p, &p, &n, &one,
		  x, &n, x, &n, &zero, z, &p);
}
*/

void dsa_transpose (double* x, int n, int p, double* xt) {
  int i;

  for (i = 0; i < n*p; i++)
    xt[i] = x[(i / p) + (i % p) * n];
}


#ifdef __USE_NAG__

void dsa_fit_mlr (double* X, int n, int tdx, int ip, double* Y, 
		  double* weights, double tol, double* coefficients, 
		  double* residuals, double* effects, int* rank, 
		  int* pivot, double* qraux, double* work, double* rss,
		  int transpose, int* fail) 
{
  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW)
    Rprintf("Entering dsa_fit_mlr (NAG) with: n: %d, p: %d, tdx: %d\n", 
	    n, ip, tdx); 
  
  NagError nag_fail;
  Nag_IncludeMean MEAN = Nag_MeanZero;
  Boolean SVD;
  INIT_FAIL(nag_fail);

  int j;

  Integer Irank;
  double df = 0.0;

  double* XT;

  double* se  = (double *) R_alloc(ip, sizeof(double)); 
  double* cov = (double *) R_alloc(ip*(ip+1)/2, sizeof(double));  
  double* h   = (double *) R_alloc(n, sizeof(double)); 
  double* q   = (double *) R_alloc(n*(ip+1), sizeof(double)); 
  double* pp  = (double *) R_alloc(2*ip + ip*ip, sizeof(double));  
  double* com_ar = (double *) R_alloc(5*(ip-1) + ip*ip, sizeof(double));
  
  Integer* sx = (Integer *) R_alloc(n, sizeof(Integer)); 
  memset(sx, 0, tdx*sizeof(Integer));
  for (j = 0; j < ip; j++) {
    sx[j] = 1; /** include the columns of the design matrix. **/
  }

  /** this does not have to happen like this, but for now it seems ok. **/
  if (transpose != 0) { 
    XT = (double *) R_alloc(n*tdx, sizeof(double));  
    dsa_transpose(X, n, tdx, XT);
  }
  else {
    XT = X;
  }
 
  g02dac(MEAN, n, XT, tdx, tdx, sx, ip, Y, weights, rss, &df, coefficients, se,
	 cov, residuals, h, q, (ip + 1), &SVD, &Irank, pp, tol, com_ar, &nag_fail);
 
  *rank = (int) Irank;

  if (nag_fail.code != NE_NOERROR) {
    if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW)
      Rprintf("nag failed with: %s\n", nag_fail.message);
  }
  
  if(SVD) {
    if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_WARN)
      Rprintf("\nModel ignored because the covariate matrix is not of full rank. ");
    *rss = R_PosInf;
  }
  
}

void dsa_fit_lr (double* X, int n, int tdx, int ip, double* Y, 
		 double tol, double* coefficients, 
		 double* residuals, double* effects, int* rank, 
		 int* pivot, double* qraux, double* work, double* dev,
		 double* weights, int max_iter, double* new_coefficients, 
		 double* WX, double* y_shifted, int transpose, int* error)

{
  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW)
    Rprintf("Entering dsa_fit_lr (NAG) with: N: %d, P: %d, tdx: %d\n", n, ip, tdx); 
  
  NagError fail;
  Nag_IncludeMean MEAN = Nag_MeanZero;
  Boolean SVD;

  INIT_FAIL(fail);
  Nag_Link LINK = Nag_Logistic;

  Integer r = 1; 
  int j; 

  double df   = 0.0;
  double* v   = (double *) R_alloc(n*(ip+6), sizeof(double)); 
  double* cov = (double *) R_alloc(ip*(ip+1)/2, sizeof(double)); 
  Integer* sx = (Integer *) R_alloc(tdx, sizeof(Integer)); 
  double* bt  = (double *) R_alloc(n, sizeof(double));
  
  memset(sx, 0, sizeof(Integer)*tdx);
  for (j = 0; j < ip; j++) {
    sx[j] = 1; /** include the columns of the design matrix. **/
  }
    
  for (j = 0; j < n; j++) {
    bt[j] = 1.0;
  }
  
  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_TRACE) {
    for (j = 0; j < tdx; j++) { 
      Rprintf("sx[%d]: %d\n", j, sx[j]);
    }
  }
  
  if (transpose != 0) { 
    if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW)
      Rprintf("transposing X.\n");
    dsa_transpose(X, n, tdx, WX);
  }
  else {
    memcpy(WX, X, sizeof(double)*n*tdx);
  }
  
  g02gbc(LINK, MEAN, n, WX, tdx, tdx, sx, ip, Y, bt, weights, (double *)0,
	 dev, &df, coefficients, &r, new_coefficients, cov, v, (ip+6), tol, max_iter, 
	 0, NULL, tol, &fail);

  *rank = (int) r;

  if (fail.code != NE_NOERROR) {
    if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_WARN)
      Rprintf("NAG failed with: %s\n", fail.message);
    
    if (fail.code == NE_VALUE_AT_BOUNDARY_B) {
      /*      *dev = R_PosInf; */
      *error = DSA_MODEL_ERROR_AT_BOUNDARY;
    }
    else if (fail.code == NE_RANK_CHANGED) {
      /*      *dev = R_PosInf; */
      *error = DSA_MODEL_ERROR_RANK_CHANGE;
    }
    else if (fail.code == NE_LSQ_ITER_NOT_CONV) {
      *error = DSA_MODEL_ERROR_NO_CONVERGENCE;
      /*      *dev = R_PosInf; */
    }
    else {
      Rprintf("NAG unhandled error: %s", fail.message);
      *error = DSA_MODEL_ERROR_UNKNOWN;
      /*      *dev = R_PosInf; */
    }
  }
  
}

#else /** __USE_NAG__ **/ 

void dsa_fit_mlr (double* X, int n, int tdx, int ip, double* Y, 
		  double* weights, double tol, double* coefficients, 
		  double* residuals, double* effects, int* rank, 
		  int* pivot, double* qraux, double* work, double* rss,
		  int transpose, int* fail)
{
  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW)
    Rprintf("Entering R_DSA_fit_mlr with: n: %d, ip: %d, tdx: %d, transpose: %d\n", n, ip, tdx, transpose);

  int k = 0, j = 0, i = 0, one = 1;
  double RSS = 0;
  double* XT; 
  int two = 2;
 
  *fail = DSA_MODEL_NO_ERROR;

  for (i = 0; i < ip; i++) {
    pivot[i] = i + 1;
  }
  memset(qraux, 0, sizeof(double)*ip);
  memset(effects, 0, sizeof(double)*n);
  memset(work, 0, sizeof(double)*2*ip);
  memset(coefficients, 0, sizeof(double)*tdx);
    
  /* we must allocate to transpose, the ensuing calculations are performed on the 
     transpose and then this is freed. */
  if (transpose == 0) {
    if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW)
      Rprintf("In R_DSA_fit_mlr, transposing.\n");
    
    XT = (double *) R_alloc(n*tdx, sizeof(double));  
    dsa_transpose(X, tdx, n, XT);
  }
  else {
    XT = X; 
  }

  /* apply weights to the transposed matrix. */
  for (j = 0; j < n; j++) {
    for (i = 0; i < tdx; i++) {
      XT[j + i*n] = XT[j + i*n] * sqrt(weights[j]);
    }
    Y[j] = Y[j] * sqrt(weights[j]);
  }

  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_TRACE)
    Rprintf("In R_DSA_fit_mlr, applied weights\n");

  F77_CALL(dqrls)(XT, &n, &ip, Y, &one, &tol, coefficients, residuals, 
		  effects, rank, pivot, qraux, work);

  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_TRACE)
    Rprintf("Fit model, calculating rss\n");

  if (*rank < ip) { 
    if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_WARN)
      Rprintf("Rank deficiency, setting RSS to infinity.");
    
    *rss = R_PosInf;
  }
  else {
    for(k = 0; k < n; k++) {
      RSS += R_pow_di(residuals[k], two);
    }
    *rss = RSS; 
  } 
  
  /* pivot the coefficients. */
  for (j = 0; j < ip; j++) {
    work[j] = coefficients[pivot[j] - 1]; 
  }
  for (j = 0; j < ip; j++) {
    coefficients[j] = work[j];
  }
  
  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW)
    Rprintf("Calculated rss of: %f\n", RSS);
  
}


/**
 * X:            the design matrix. (it will be transposed if the transpose != 0).
 * n:            the number of independent observations. 
 * tdx:          the number of columns in X
 * ip:           the number of independent variables, these will occurr in the first p - 1 columns. 
 *
 * weights:      an (n x 1) vector of weights which will be applied to the observations. 
                 
*/
void dsa_fit_lr (double* X, int n, int tdx, int ip, double* Y, 
		 double tol, double* coefficients, 
		 double* residuals, double* effects, int* rank, 
		 int* pivot, double* qraux, double* work, double* dev,
		 double* weights, int max_iter, double* new_coefficients, 
		 double* WX, double* y_shifted, int transpose, int* error)
{
  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW)
    Rprintf("Entering R_DSA_fit_glm_lr with: n: %d, ip: %d, tdx: %d\n", n, ip, tdx);
  
  int iter = 0, k = 0, j = 0, i = 0; 
  double sum = 0, theta_j = 0, pi_j = 0, var_pi_j = 0;
  double abs_diff = 0.0;
  double DEV = 0.0, DEV_OLD = 0.0;
  int one = 1;
  double sqrt_var_pi_j,sqrt_weight_j;

  *error = DSA_MODEL_NO_ERROR;
  
  if (tol == 0) { 
    if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_WARN)
      Rprintf("tolerance is precisely 0 resetting to something sensible.\n");
    tol = 1e-9;
  }

  /* zero the coefficients in case they haven't been zeroed. */
  memset(coefficients, 0.0, sizeof(double)*tdx); 

  if (transpose == 0) {
    dsa_transpose(X, tdx, n, WX);
    memcpy(X, WX, sizeof(double)*n*tdx); 
  }
    
  while (1) {
    if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW)
      Rprintf("Iteration: %d\n", iter);
    
    if (iter++ >= max_iter) {
      *error = DSA_MODEL_ERROR_NO_CONVERGENCE;
      DEV = R_PosInf;
      break;
    }
    
    DEV = 0.0;
    
    for (j = 0; j < n; j++) {
      theta_j = 0;
      
      for (i = 0; i < ip; i++) {
	if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_TRACE)
	  Rprintf("%f * %f\n", X[j + i*n], coefficients[i]);
	
 	theta_j += X[j + i*n] * coefficients[i];
      }

      pi_j = exp(theta_j)/(1 + exp(theta_j));
      var_pi_j = pi_j * (1 - pi_j); 

      /* decide if we are at a boundary or not. */
      /*
      if (pi_j > 1.0 - tol || pi_j < tol) {
	if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW) 
	  Rprintf("model at boundary; aborting.\n");
	*error = DSA_MODEL_ERROR_AT_BOUNDARY;
	DEV = R_PosInf;
	break;
      }
      */
      
      sqrt_var_pi_j = sqrt(var_pi_j);
      sqrt_weight_j = sqrt(weights[j]);

      /* avoid dividing by 0, this effectively ignores the row for this fit. */
      y_shifted[j] = (sqrt_var_pi_j != 0) ?
	sqrt_weight_j * (Y[j] - pi_j)/sqrt_var_pi_j :
	0;

      DEV += 2*weights[j]*(y_log_y(Y[j], pi_j) + y_log_y(1 - Y[j], 1 - pi_j));
      
      for (i = 0; i < ip; i++) {
	WX[j + i*n] = X[j + i*n] * sqrt_var_pi_j * sqrt_weight_j; 
      }
      
      if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_TRACE) {
	Rprintf("theta_j: %f\n", theta_j);
	Rprintf("pi_j: %f\n", pi_j);
	Rprintf("weight[%d] = %f\n", j, weights[j]);
	Rprintf("y_shifted[%d] = %f\n", j, y_shifted[j]);
      }
    }

    F77_CALL(dqrls)(WX, &n, &ip, y_shifted, &one, &tol, new_coefficients, residuals, 
		    effects, rank, pivot, qraux, work);

    if (*rank < ip) {
      *error = DSA_MODEL_ERROR_RANK_CHANGE;
      DEV = R_PosInf;
      break;
    }

    abs_diff = 0.0;

    for (j = 0; j < ip; j++) {
      /* converge based on the absolute difference in the coefficients. */
      /* abs_diff += fabs(new_coefficients[j]); */
    
      /* pivot the coefficients. */
      coefficients[j] += new_coefficients[pivot[j] - 1]; 

      /* zero the coefficients. */
      new_coefficients[pivot[j] - 1] = 0;

      /* reset the pivots. */
      pivot[j] = j + 1; 
    }
    
    if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW)
      Rprintf("deviance: %g\n", DEV);

    if (fabs(DEV_OLD - DEV)/(0.1 + fabs(DEV)) < tol)
      break;

    /** R uses the difference in deviance to decide when to exit. */
    DEV_OLD = DEV;
  }
  
  /* set some things on the way out. */
  *dev = DEV;

  /* transpose on the way out for cleanliness. */
  if (transpose == 0) {
    dsa_transpose(X, n, tdx, WX);
    memcpy(X, WX, sizeof(double)*n*tdx);
  }
}

#endif /* __USE_NAG__ */

/** hooks to call the functions from R. **/

void R_DSA_transpose (double* x, int* n, int* p, double* xt) {
  dsa_transpose(x, *n, *p, xt);
}

void R_DSA_fit_mlr(double* X, int* n, int* tdx, int* ip, double* Y, 
		   double* weights, double* tol, double* coefficients, 
		   double* residuals, double* effects, int* rank, 
		   int* pivot, double* qraux, double* work, double* rss,
		   int* transpose, int* fail) 
{  
  dsa_fit_mlr(X, *n, *tdx, *ip, Y, weights, *tol, coefficients, 
	      residuals, effects,
	      rank, pivot, qraux, work, rss, *transpose, fail);
}


void R_DSA_fit_glm_lr(double* X, int* n, int* tdx, int* ip, double* Y, 
		      double* tol, double* coefficients, 
		      double* residuals, double* effects, int* rank, 
		      int* pivot, double* qraux, double* work, double* dev,
		      double* weights, int* max_iter, double* new_coefficients, 
		      double* WX, double* y_shifted, int* transpose, int* error)
{
  dsa_fit_lr (X, *n, *tdx, *ip, Y, *tol, coefficients,
	      residuals, effects, rank, 
	      pivot, qraux, work, dev,
	      weights, *max_iter, new_coefficients, WX, y_shifted,
	      *transpose, error); 
}

void R_DSA_set_debug_level(int* level) 
{
  R_DSA_DEBUG_LEVEL = *level;
}

/*************************************** Interface to DSA ******************************************/

void run_ols(double* dataframeX, int samplesize, int TOTALTERMS, 
	     int IP, double* dataframeY, double* WT, 
	     double TOL, double* PARA, double* RES, int* rank, double* RSS, 
	     int* fail) 
{
  double* qraux = (double *) R_alloc(IP, sizeof(double)); 
  double* effects = (double *) R_alloc(samplesize, sizeof(double));
  int* pivot = (int *) R_alloc(IP, sizeof(int));
  double* work  = (double *) R_alloc(2*IP, sizeof(double));
  int i;
    
  for (i = 0; i < IP; i++) {
    pivot[i] = i + 1;
  }
	
  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW) 
    Rprintf("calling dsa_fit_mlr\n");
   
  dsa_fit_mlr(dataframeX, samplesize, TOTALTERMS, IP, dataframeY, WT, 
	      TOL, PARA,
	      RES, effects, rank, pivot, qraux, work, RSS, 0, fail);
   
  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW) 
    Rprintf("returning from dsa_fit_mlr\n");

  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW) {
    Rprintf("Rank: %d, Fail: %d RSS: %f \n", *rank, *fail, *RSS);
    for (i = 0; i < IP; i++)
      Rprintf("coefficient[%d] = %f\n", i, PARA[i]);
  }
}

void R_run_ols(double* dataframeX, int* samplesize, int* TOTALTERMS, int* IP, double* dataframeY, double* WT, 
	       double* TOL, double* PARA, double* RES) /* , int* rank, double* RSS, int* fail) */
{
  int rank    = 0;
  double RSS  = 0;
  int fail    = 0;

  run_ols(dataframeX, *samplesize, *TOTALTERMS, *IP, dataframeY, WT, *TOL, PARA, RES, &rank, &RSS, &fail); 
}

void run_lr(double* dataframeX, int samplesize, int TOTALTERMS, int IP, double* dataframeY, 
	    double TOL, double* PARA, double* RES, int* rank, double* WT,  
	    double* DEV, int* fail) 
{
  int i; 
  
  double* new_coefficients = (double *) R_alloc(TOTALTERMS, sizeof(double)); 
  double* WX = (double *) R_alloc((TOTALTERMS) * samplesize, sizeof(double)); 
  double* Z = (double *) R_alloc(samplesize, sizeof(double)); 
  double* qraux = (double *) R_alloc(IP, sizeof(double)); 
  double* effects = (double *) R_alloc(samplesize, sizeof(double));
  int* pivot = (int *) R_alloc(IP, sizeof(int));
  double* work  = (double *) R_alloc(2*(IP), sizeof(double));
  
  for (i = 0; i < IP; i++) {
    pivot[i] = i + 1;
  }
  
  dsa_fit_lr(dataframeX, samplesize, TOTALTERMS, IP, dataframeY, TOL, PARA, RES, effects, rank, 
	     pivot, qraux, work, DEV, WT, 25, new_coefficients, WX, Z, 0, fail);
  
  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW) 
    Rprintf("Returned from the dsa_fit_lr, ending rank: %d\n", *rank);
  	
  if (R_DSA_DEBUG_LEVEL <= R_DSA_DEBUG_LOW)
    for (i = 0; i < IP; i++)
      Rprintf("coefficient[%d] = %f\n", i, PARA[i]);
}

void R_run_lr(double* dataframeX, int* samplesize, int* TOTALTERMS, int* IP, double* dataframeY, 
	      double* TOL, double* PARA, double* RES, double* WT, int* rank, double* DEV, int* fail)
{
  run_lr(dataframeX, *samplesize, *TOTALTERMS, *IP, dataframeY, 
	 *TOL, PARA, RES, rank, WT, DEV, fail);
}
