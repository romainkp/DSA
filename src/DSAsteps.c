#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>

#include "dsa.h"
#include "dsa_model.h"
#include "PRIMES.h"

#define DOUBLE_EPS 2.2204460492503131e-16

void (*model_fitter)(
  int*, double*, double*, int*, int*, int*, double*, double*, double*, double*,
  int*, double*, int*, int*, double*
);

static const double THRESH = 30.;
static const double MTHRESH = -30.;
static const double INVEPS = 1/DOUBLE_EPS;
static R_INLINE double x_d_opx(double x) {return x/(1 + x);}

double __identity_link(double x)
{
  return x;
}

double __logit_link(double x) {
  return  (x_d_opx((x < MTHRESH) ? DOUBLE_EPS : ((x > THRESH) ? INVEPS : exp(x))));
}

/*Need to avoid duplication of this code from modelUtils in the future*/
double DSA_PACK_multinomial_ll(double* X, double* Y, int N, int P, int M, int tdx,
			       double* weights, double* betas) 
{
  int n, m, p;
  double denom, sm, eta, yy;
  double ll = 0;

  double* mus = (double *)R_alloc(M, sizeof(double));

  for (n = 0; n < N; n++) {
    if(abs(weights[n])>DSA_PACK_epsilonCompare) 
      { 
	denom = sm = yy = 0;
    
	for (m = 0; m < M; m++) {
	  eta = 0;
	  for (p = 0; p < P; p++) {
	    eta += (X[n+p*N] * betas[p + m*P]); /* we have to jump through the rows. */
	  }
	  /*mus[m] = exp(eta);*/
	  mus[m] = ((eta < MTHRESH) ? DOUBLE_EPS :
		    ((eta > THRESH) ? INVEPS : exp(eta)));
	  denom += mus[m];
	  yy += Y[m*N + n];
	}
  
	for (m = 0; m < M; m++) {
	  mus[m] = mus[m]/(1 + denom);
	  /*sm += mus[m];*/
	}

	for (m = 0; m < M; m++) {
	  ll += Y[m*N + n] * log(mus[m])*weights[n];
	}
	/*ll += (1 - yy) * log(1 - sm)*weights[n];*/
	ll += (1 - yy) * log(1/(1+denom))*weights[n];
      }
  }

  return -ll;
}

double DSA_PACK_l2_loss(double (*link)(double),double* X, double* Y, int N, int P, int M, 
			int tdx, double* glmWT, double* binWTdata, double* estimates)
{
  double loss = 0.0;
  double tmp;
  int i,j;

  for(i = 0; i < N; i++) {
    if(abs(glmWT[i])>DSA_PACK_epsilonCompare) 
      { 
	tmp = 0;
	for(j = 0; j < P; j++) {
	  tmp += estimates[j] * X[i+j*N];
	}
	tmp = R_pow_di((Y[i] - ((*link)(tmp))), 2);
	loss += ((glmWT[i]*binWTdata[i])*tmp); 
      } 
  }
  return loss;
}
  
/**
 * we could definitely think of better names here.
 */
void DSA_PACK_evaluate_loss(int lossFunction, double* X, double* Y, int N, int P, int M, 
			    int tdx, double* glmWT, double* binWTdata, double* estimates, 
			    double* loss) 
{
  /** L2(identity) == 1 */
  if (lossFunction == 1) {
    *loss = DSA_PACK_l2_loss(&(__identity_link), X, Y, N, P, M, tdx, glmWT, binWTdata, estimates);
  }
  /** L2(logit) == 1 */
  else if (lossFunction == 2) {
    *loss = DSA_PACK_l2_loss(&(__logit_link), X, Y, N, P, M, tdx, glmWT/* +4*P+N*(2*P+M+2) */, binWTdata, estimates);
  }
  else if (lossFunction == 3) {
    *loss = DSA_PACK_multinomial_ll(X, Y, N, P, M, tdx, glmWT, estimates);
  }
  else {
    error("\n Fatal error. Aborting C subroutine (Unknown loss function)."); 
  }
}


/* ------------------------------------------------------------------------
   Substitution Moves on Learning Set
   ------------------------------------------------------------------------ */

int DSA_PACK_Sstep(int *gCOUNT,int *gCONST,double *bestarss,int **bestmodels,int *tempmodel,
		   int *currentmodel,int *nextmodel,int *newterm,double *Xdata,double *tempXdata,
		   double *Ydata,double *tempYdata,double *tempterm,double *glmWT,double *PARA,
		   double *RES,int *ithterm,double *WTdata,int trainFlag,double *currentarss,int nforced,
		   double *binWTdata,int *modelfit_work_i,double *modelfit_work_d,int *lastXdesign)
{
  int IP,msize,samplesize,nvarX,maxsize,nfolds,maxsumofpow,orderint,i,j,pl,newtermindex,usetree;
  int jthvar,nglmcall,orderint_in_newterm,sumofpow_in_newterm,CENSOR,is_newterm_unique,binind;
  double oldrss,newrss;

  /* Initialization */  
  IP = gCOUNT[0]+1;
  msize = gCOUNT[0];
  nglmcall = gCOUNT[2];

  maxsize = gCONST[0];
  maxsumofpow = gCONST[2];
  nfolds = gCONST[3];
  nvarX = gCONST[5];
  CENSOR = gCONST[6];
  binind = gCONST[9];
  usetree = gCONST[11];

  if(trainFlag == 1)
    {
      samplesize = gCOUNT[4];
      orderint = gCOUNT[1];
    }
  else 
    {
      samplesize = gCONST[4];
      orderint = gCONST[8];
    }
  oldrss = currentarss[IP-1];
  newrss = currentarss[IP-1];

  /* Start all the sub moves and swap moves */
  for(newtermindex=nforced;newtermindex<msize;newtermindex++)
    for(jthvar=0;jthvar<nvarX;jthvar++) {
      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\n\nSUB/SWAP on term %i with var %i",newtermindex,jthvar);
      for(i=0;i<maxsize;i++)for(j=0;j<nvarX;j++)tempmodel[j*maxsize+i] = currentmodel[j*maxsize+i];
      for(j=0;j<nvarX;j++)newterm[j] = tempmodel[j*maxsize+newtermindex];

      /* TRY ADDING ONE TO P-VECTOR */
      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\n\tTry to add 1 to var %i",jthvar);
      newterm[jthvar] = newterm[jthvar] + 1;
	  
      DSA_PACK_checknewterm(&orderint_in_newterm,&sumofpow_in_newterm,&is_newterm_unique,newterm,nvarX,
			    tempmodel,ithterm,maxsize,newtermindex,0);
      /* Test if submove accepted */
      if(is_newterm_unique == 1 && sumofpow_in_newterm <= maxsumofpow && orderint_in_newterm <= orderint)
	{
	  if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)
	    Rprintf("\n\tadd 1 accepted");
		  
	  DSA_PACK_eval_model(maxsize,msize,tempmodel,newtermindex,newterm,nvarX,samplesize,
			      tempXdata,Xdata,tempterm,CENSOR,glmWT,WTdata,IP,tempYdata,Ydata,PARA,
			      RES,&nglmcall,&newrss,nextmodel,0,binind, binWTdata, 0,usetree,
			      modelfit_work_i,modelfit_work_d,lastXdesign);
	}
      else
	{
	  if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\n\tadd 1 rejected");
	  /* Test if swap moves allowed and possible */
	  if(is_newterm_unique == 1 && orderint_in_newterm > orderint)
	    {
	      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\n\t\tTry Swap moves");
	      /* TRY ALTERNATE-SUB MOVES */
	      for(pl=0;pl<nvarX;pl++)
		{
		  if(newterm[pl] != 0)
		    {
		      i=newterm[pl];
		      newterm[pl] = 0;
		      DSA_PACK_checknewterm(&orderint_in_newterm,&sumofpow_in_newterm,&is_newterm_unique,
					    newterm,nvarX,tempmodel,ithterm,maxsize,newtermindex,0);

		      /*Test if swap move accepted */
		      if(is_newterm_unique == 1 && sumofpow_in_newterm <= maxsumofpow && 
			 orderint_in_newterm <= orderint)
			{
			  if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\n\t\tSubswap accepted");
			  DSA_PACK_eval_model(maxsize,msize,tempmodel,newtermindex,newterm,
					      nvarX,samplesize,tempXdata,Xdata,tempterm,CENSOR,
					      glmWT,WTdata,IP,tempYdata,Ydata,PARA,
					      RES,&nglmcall,&newrss,nextmodel,0,binind,
					      binWTdata,0,usetree,modelfit_work_i,
					      modelfit_work_d,lastXdesign);
			}
		      else
			{
			  if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\n\t\tSubswap rejected");
			}
		      newterm[pl] = i;
		    }
		}
	    }
	} /* end of else do swaps */
      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\n\tTry to subs 1 to var %i",jthvar);
      for(i=0;i<maxsize;i++)for(j=0;j<nvarX;j++)tempmodel[j*maxsize+i] = currentmodel[j*maxsize+i];
      for(j=0;j<nvarX;j++)newterm[j] = tempmodel[j*maxsize+newtermindex];

      /* TRY SUBTRACTING ONE TO P-VECTOR */
      if(newterm[jthvar] > 0)
	{
	  if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\n\tPossible to subs");
	  newterm[jthvar] = newterm[jthvar] - 1;
	  DSA_PACK_checknewterm(&orderint_in_newterm,&sumofpow_in_newterm,&is_newterm_unique,newterm,
				nvarX,tempmodel,ithterm,maxsize,newtermindex,0);

	  /* Test if submove accepted */
	  if(is_newterm_unique == 1 && sumofpow_in_newterm <= maxsumofpow && 
	     orderint_in_newterm <= orderint && orderint_in_newterm>=1)
	    {
	      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\n\tsubs 1 accepted");
	      DSA_PACK_eval_model(maxsize,msize,tempmodel,newtermindex,newterm,nvarX,samplesize,
				  tempXdata,Xdata,tempterm,CENSOR,glmWT,WTdata,IP,tempYdata,Ydata,
				  PARA,RES,&nglmcall,&newrss,nextmodel,0,binind, binWTdata,0,usetree,
				  modelfit_work_i,modelfit_work_d,lastXdesign);

	    }
	  else
	    {
	      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\n\tsubs 1 rejected");
	    }
	}
    }   /* loops over i-th term and jth variable  */

  return(DSA_PACK_exportresult(&oldrss,newrss,msize,nvarX,currentmodel,nextmodel,maxsize,bestarss,
			       gCOUNT,IP,nglmcall,bestmodels,0,currentarss));
}       


/* ------------------------------------------------------------------------
   Addition Moves on Learning Set
   ------------------------------------------------------------------------ */

int DSA_PACK_Astep(int *gCOUNT,int *gCONST,double *bestarss,int **bestmodels,int *tempmodel,int *currentmodel,
		   int *nextmodel,int *newterm,double *Xdata,double *tempXdata,double *Ydata,
		   double *tempYdata,double *tempterm,double *glmWT,double *PARA,
		   double *RES,int *ithterm,double *WTdata,int trainFlag,double *currentarss,double *binWTdata,
		   int *modelfit_work_i,double *modelfit_work_d,int *lastXdesign)
{
  int IP,msize,nvarX,maxsize,nfolds,maxsumofpow,orderint,i,j,pl,newtermindex,jthvar,samplesize,nglmcall,CENSOR,orderint_in_newterm,sumofpow_in_newterm,is_newterm_unique,binind,usetree;
  double oldrss,newrss;

  /* Initialization */  
  IP = gCOUNT[0]+1;
  msize = gCOUNT[0];
  nglmcall = gCOUNT[2];

  maxsize = gCONST[0];
  maxsumofpow = gCONST[2];
  nfolds = gCONST[3];
  nvarX = gCONST[5];
  CENSOR = gCONST[6];
  binind = gCONST[9];
  usetree = gCONST[11];
  if(trainFlag == 1)
    {
      samplesize = gCOUNT[4];
      orderint = gCOUNT[1];
    }
  else 
    {
      samplesize = gCONST[4];
      orderint = gCONST[8];
    }
  oldrss = R_PosInf;
  newrss = R_PosInf;

  /* Start all the addsub moves and addswap moves */
  for(newtermindex=0;newtermindex<msize;newtermindex++)
    for(jthvar=0;jthvar<nvarX;jthvar++) {
      for(i=0;i<maxsize;i++)
	for(j=0;j<nvarX;j++)
	  tempmodel[j*maxsize+i] = currentmodel[j*maxsize+i];
      for(j=0;j<nvarX;j++)
	newterm[j] = tempmodel[j*maxsize+newtermindex];
	    
      /* TRY ADDING ONE TO P-VECTOR */
      newterm[jthvar] = newterm[jthvar] + 1;
      DSA_PACK_checknewterm(&orderint_in_newterm,&sumofpow_in_newterm,&is_newterm_unique,newterm,
			    nvarX,tempmodel,ithterm,maxsize,newtermindex,1);

      /* Test if addsub move accepted */
      if(is_newterm_unique == 1 && sumofpow_in_newterm <= maxsumofpow && orderint_in_newterm <= orderint)
	{
	  msize++;
	  IP++;
		    
	  DSA_PACK_eval_model(maxsize,msize,tempmodel,msize-1,newterm,nvarX,samplesize,tempXdata,
			      Xdata,tempterm,CENSOR,glmWT,WTdata,IP,tempYdata,Ydata,PARA,RES,
			      &nglmcall,&newrss,nextmodel,0,binind,binWTdata,0,usetree,
			      modelfit_work_i,modelfit_work_d,lastXdesign);
		    
	  IP--;
	  msize--;
	} 
      else
	{
	  /* Test if swap moves allowed and possible */
	  if(is_newterm_unique == 1 && orderint_in_newterm > orderint)
	    {
	      /* TRY ALTERNATE-SUBADD MOVES */
	      for(pl=0;pl<nvarX;pl++)
		{
		  if(newterm[pl] != 0)
		    {
		      i=newterm[pl];
		      newterm[pl] = 0;
		      DSA_PACK_checknewterm(&orderint_in_newterm,&sumofpow_in_newterm,&is_newterm_unique,
					    newterm,nvarX,tempmodel,ithterm,maxsize,newtermindex,1);
		      /*Test if addswap move accepted */
		      if(is_newterm_unique == 1 && sumofpow_in_newterm <= maxsumofpow && 
			 orderint_in_newterm <= orderint)
			{
			  msize++;
			  IP++;
						    
			  DSA_PACK_eval_model(maxsize,msize,tempmodel,msize-1,newterm,nvarX,
					      samplesize,tempXdata,Xdata,tempterm,CENSOR,glmWT,
					      WTdata,
					      IP,tempYdata,Ydata,PARA,RES,&nglmcall,
					      &newrss,nextmodel,0,binind, binWTdata,0,usetree,
					      modelfit_work_i,modelfit_work_d,lastXdesign);
			  IP--;
			  msize--;
			}
		      newterm[pl] = i;
		    }
		}
	    }
	} /* end of else do swaps */
	    
      for(i=0;i<maxsize;i++)for(j=0;j<nvarX;j++)tempmodel[j*maxsize+i] = currentmodel[j*maxsize+i];
      for(j=0;j<nvarX;j++)newterm[j] = tempmodel[j*maxsize+newtermindex];
	    
      /* TRY SUBTRACTING ONE TO P-VECTOR */
      if(newterm[jthvar] > 0)
	{
	  newterm[jthvar] = newterm[jthvar] - 1;
	  DSA_PACK_checknewterm(&orderint_in_newterm,&sumofpow_in_newterm,&is_newterm_unique,
				newterm,nvarX,tempmodel,ithterm,maxsize,newtermindex,1);
	  /* Test if addsub move accepted */
	  if(is_newterm_unique == 1 && sumofpow_in_newterm <= maxsumofpow && orderint_in_newterm <= orderint && 
	     orderint_in_newterm>=1)
	    {
	      msize++;
	      IP++;
	      DSA_PACK_eval_model(maxsize,msize,tempmodel,msize-1,newterm,nvarX,samplesize,tempXdata,
				  Xdata,tempterm,CENSOR,glmWT,WTdata,IP,tempYdata,Ydata,PARA,RES,
				  &nglmcall,&newrss,nextmodel,0,binind, binWTdata,0,usetree,
				  modelfit_work_i,modelfit_work_d,lastXdesign);
	      IP--;
	      msize--;
	    } 
	}
    }   /* loops over i-th term and jth variable */
    
  /* TRY ADDING A MAIN TERM (UNIT VECTORS) */
  for(jthvar=0;jthvar<nvarX;jthvar++)
    {
      for(i=0;i<maxsize;i++)for(j=0;j<nvarX;j++)tempmodel[j*maxsize+i] = currentmodel[j*maxsize+i];
      memset(newterm, 0, sizeof(int)*nvarX);
      newterm[jthvar] = 1;
      DSA_PACK_checknewterm(&orderint_in_newterm,&sumofpow_in_newterm,&is_newterm_unique,newterm,nvarX,
			    tempmodel,ithterm,maxsize,newtermindex,1);
	    
      /* Test if addmain move accepted */
      if(is_newterm_unique == 1)
	{
	  msize++;
	  IP++;
	  DSA_PACK_eval_model(maxsize,msize,tempmodel,msize-1,newterm,nvarX,samplesize,tempXdata,Xdata,
			      tempterm,CENSOR,glmWT,WTdata,IP,tempYdata,Ydata,PARA,RES,
			      &nglmcall,&newrss,nextmodel,0,binind, binWTdata,0,usetree,
			      modelfit_work_i,modelfit_work_d,lastXdesign);
	  IP--;
	  msize--;
	}
    }

  msize++;
  IP++;
  return(DSA_PACK_exportresult(&oldrss,newrss,msize,nvarX,currentmodel,nextmodel,maxsize,bestarss,gCOUNT,IP,
			       nglmcall,bestmodels,0,currentarss));
}       /* end of addition routine */      


/* ------------------------------------------------------------------------
   Deletion Moves on Learning Set
   ------------------------------------------------------------------------ */

int DSA_PACK_Dstep(int *gCOUNT,int *gCONST,double *bestarss,int **bestmodels,int *tempmodel,int *currentmodel,
		   int *nextmodel,int *newterm,double *Xdata,double *tempXdata,double *Ydata,double *tempYdata,
		   double *tempterm,double *glmWT,double *PARA,
		   double *RES,double *WTdata,int trainFlag,double *currentarss,int nforced,double *binWTdata,
		   int *modelfit_work_i,double *modelfit_work_d,int *lastXdesign)
{
  int IP,msize1,msize2,nvarX,maxsize,nfolds,maxsumofpow,orderint,newtermindex,samplesize,nglmcall,CENSOR,binind,usetree;
  double EPSILON,oldrss, newrss;

  /* Initialization */  
  IP = gCOUNT[0]+1;
  msize1 = gCOUNT[0];
  nglmcall = gCOUNT[2];

  maxsize = gCONST[0];
  maxsumofpow = gCONST[2];
  nfolds = gCONST[3];
  nvarX = gCONST[5];
  CENSOR = gCONST[6];
  binind = gCONST[9];
  usetree = gCONST[11];
  if(trainFlag == 1)
    {
      samplesize = gCOUNT[4];
      orderint = gCOUNT[1];
    }
  else 
    {
      samplesize = gCONST[4];
      orderint = gCONST[8];
    }

  IP = IP - 1;
  msize2 = msize1 - 1;
  oldrss = bestarss[IP-1];
  newrss = bestarss[IP-1];
  EPSILON = 0.01;

  /* Start all the del moves */
  for(newtermindex=nforced;newtermindex<msize1;newtermindex++)
    {
      memset(tempmodel, 0, sizeof(int)*maxsize*nvarX);
      DSA_PACK_eval_model(maxsize,msize2,tempmodel,newtermindex,currentmodel,nvarX,samplesize,
			  tempXdata,Xdata,tempterm,CENSOR,glmWT,WTdata,IP,tempYdata,Ydata,PARA,RES,
			  &nglmcall,&newrss,nextmodel,1,binind, binWTdata,0,usetree,
			  modelfit_work_i,modelfit_work_d,lastXdesign);
    }
  IP++; /* term not yet deleted */

  if((oldrss-newrss)>DSA_PACK_epsilonCompare)
    {
      msize1--;
      IP--;
    }
  return(DSA_PACK_exportresult(&oldrss,newrss,msize1,nvarX,currentmodel,nextmodel,maxsize,bestarss,gCOUNT,IP,
			       nglmcall,bestmodels,1,currentarss));
}       /* end of deletion routine */      


void DSA_PACK_eval_model(int maxsize,int msize,int *tempmodel,int newterm_index,
			 int *changeindexset,int nvarX, int samplesize,double *tempXdata,double *Xdata,
			 double *tempterm, int CENSOR,double *glmWT,double *WTdata, int IP,
			 double *tempYdata, double *Ydata,double *PARA,double *RES,
			 int *nglmcall,double *newrss,int *nextmodel,int inddel, int binind, double* binWTdata,
			 int always_doit,int usetree,int *modelfit_work_i,double *modelfit_work_d,int *lastXdesign) /*always_doit means fit the model even if in the tree, e.g. to fit in on the validation sets*/
{
  int i, j, k, TDQ, IRANK = IP, NnonNA, MAXITER, PRINTITER,ii;
  double DEV = 0.0, IDF, TOL, temprss,*WTnoNA;

  int fail = DSA_MODEL_NO_ERROR;
  int M;
  int ms = maxsize + 1;
  double tol = 1e-8;
  int maxiter = 25;
  int mtype = min(binind + 1,3);

  /*M:*/
  if(binind>1)
    M=binind;
  else 
    M=1;
 
  if(inddel==0) {
    for(i=0;i<nvarX;i++)
      tempmodel[i*maxsize+newterm_index] = changeindexset[i];
  }
  else {
    k = 0;
    for(i=0;i<maxsize;i++) {
      if(i != newterm_index)  {
	for(j=0;j<nvarX;j++)
	  tempmodel[j*maxsize+k] = changeindexset[j*maxsize+i];
	k++;
      }
    }
  }

  p_node_t* temp_node = NULL;  
  unsigned long int term_composite = 1;

  if(usetree==1)
    {
      for (j = 0; j < msize; j++) {
	term_composite = 1;
	for (k = 0; k < nvarX; k++) {
	  /*term_composite = (term_composite * pow(PRIMES[k], tempmodel[k*maxsize + j]));*/
	  term_composite = (term_composite * DSA_PACK_ulongpow(PRIMES[k], tempmodel[k*maxsize + j]));
	}
	MODEL_TO_FIT[j] = term_composite;
      }
    
      if (msize > 0) 
	/** sort the input vector. **/
	/*R_isort(MODEL_TO_FIT, msize);*/
	DSA_PACK_ulongsort(MODEL_TO_FIT, msize);
      temp_node = get_node(MODEL_TO_FIT, msize, TREE);
    }
    

  if (always_doit == 1 || temp_node == NULL || temp_node->rss < 0) {
    /* Print model tried on stdout */
    if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh) {
      Rprintf("\nModel tried %i:",*nglmcall+1);
      Rprintf("\n");
      for(i=0;i<msize;i++) {
	Rprintf("%d \t ",i);
	for(j=0;j<nvarX;j++)
	  if(tempmodel[j*maxsize+i] > 0)Rprintf("%d^%d ",j,tempmodel[j*maxsize+i]);
	Rprintf("\n");
      }
    } 
    /* Form tempXdata,tempYdata and glmWT */
    DSA_PACK_get_tempXdata(samplesize,tempXdata,msize,maxsize,tempmodel,
			   nvarX,Xdata,tempterm,lastXdesign);

    model_fitter = (void(*)(int*, double*, double*, int*, int*, int*, double*, double*, double*, double*, int*, double*, int*, int*, double*)) R_GetCCallable("modelUtils", "fit_model");

    (*model_fitter)(&mtype, tempXdata, tempYdata, &samplesize, &IRANK, &M,
		    glmWT, PARA, &DEV, RES, &fail, &tol, &maxiter,modelfit_work_i,modelfit_work_d);
    NnonNA=*modelfit_work_i;

    /** error checking. **/
    if (fail != DSA_MODEL_NO_ERROR) {
      DEV = R_PosInf;
	    
      if(fail == DSA_MODEL_ERROR_AT_BOUNDARY) {
	if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh) {
	  Rprintf("\nModel ignored (a fitted value is at a boundary) - rss set to infinity for:\n");
	}
      }
      else if (fail == DSA_MODEL_ERROR_RANK_CHANGE) {
	if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh) {
	  Rprintf("\nModel ignored because not of full rank.\n");
	}
      }
      else if (fail == DSA_MODEL_ERROR_NO_CONVERGENCE) {
	if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh) {
	  Rprintf("\nModel ignored (no convergence) - rss set to infinity for:\n");
	}
      }
      else {
	Rprintf("Witnessed unknown error, ignoring model. fail code is: %d\n", fail);
	error("\n Fatal error. Aborting C subroutine."); 
      }
    }
    else {
      if (mtype == 1)WTnoNA=modelfit_work_d + 3*IRANK+samplesize*(IRANK+M+1);
      else if(mtype == 2)WTnoNA=modelfit_work_d + 4*IRANK+samplesize*(2*IRANK+M+2);
      else if(mtype == 3)WTnoNA=modelfit_work_d + 2*M*IRANK*(M*IRANK+1) + samplesize*(IRANK+2*M+M*M);
      DSA_PACK_evaluate_loss(mtype, tempXdata,tempYdata, samplesize, IRANK, M, ms,
			     WTnoNA /*glmWT*/ , binWTdata, PARA, &DEV);
    }
	
    (*nglmcall)++;
    temprss = DEV/NnonNA;

    if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh) {
      Rprintf("\nCoefficients for model with ARSS of %lf based on NnonNA=%i obs: ", temprss, NnonNA);
      if(DEV==R_PosInf)
	Rprintf("not printed because the deviance is R_PosInf.");
      else if(binind>1)for(ii = 0; ii < binind; ii++)
	{
	  Rprintf("\n");
	  DSA_PACK_printmatrix(PARA+ii*IP,IP,1);
	}
      else DSA_PACK_printmatrix(PARA,IP,1);
      Rprintf("\n");
    }
  
    if(usetree==1)
      {
	/** add it to the tree. **/
	if (msize > 0 && always_doit == 0) {
	  if (temp_node == NULL)
	    add_node(MODEL_TO_FIT, msize, temprss, TREE);
	  else
	    /* this is the case when you have node, but haven't fit that model. **/
	    temp_node->rss = temprss; 
	}
      }
  }
  else {
    temprss = temp_node->rss;
  }
  

  /* Compare loss of new model to current model. */ 
  if((*newrss - temprss) > DSA_PACK_epsilonCompare) {
    (*newrss) = temprss;
    for(i=0;i<maxsize;i++)
      for(j=0;j<nvarX;j++)nextmodel[j*maxsize+i] = tempmodel[j*maxsize+i];
  }
}


void DSA_PACK_checknewterm(int *orderint_in_newterm,int *sumofpow_in_newterm,int *is_newterm_unique,
			   int *newterm,int nvarX,int *tempmodel,int *ithterm,int maxsize,
			   int newterm_index,int addmove_notsubmove)
{
  *orderint_in_newterm = DSA_PACK_num_nonzero(newterm,nvarX);
  *sumofpow_in_newterm = DSA_PACK_sumofints(newterm,nvarX);

  /* check whether the new term is unique */
  *is_newterm_unique = DSA_PACK_istermunique(newterm,tempmodel,ithterm,maxsize,nvarX,newterm_index,
					     addmove_notsubmove,*orderint_in_newterm,*sumofpow_in_newterm);
}

int DSA_PACK_istermunique(int *newterm,int *tempmodel,int *ithterm,int maxsize,int nvarX,
			  int newterm_index,int addmove_notsubmove,int orderint_in_newterm, 
			  int sumofpow_in_newterm)
{
  int j,ithtermindex,orderint_in_ithterm,sumofpow_in_ithterm,ithterm_equal_newterm,k,uniqueind;

  uniqueind = 1;
  for(ithtermindex=0;ithtermindex<maxsize;ithtermindex++)
    {
      if((addmove_notsubmove == 1) || (addmove_notsubmove == 0 && ithtermindex != newterm_index))
	{	 
	  for(j=0;j<nvarX;j++)ithterm[j] = tempmodel[j*maxsize+ithtermindex];
	  orderint_in_ithterm = DSA_PACK_num_nonzero(ithterm,nvarX);
	  sumofpow_in_ithterm = DSA_PACK_sumofints(ithterm,nvarX);
	  if(orderint_in_newterm == orderint_in_ithterm && sumofpow_in_newterm == sumofpow_in_ithterm)
	    {
	      ithterm_equal_newterm=1;
	      for(j=0;j<nvarX;j++)
		{
		  if(newterm[j] == ithterm[j])k=1;
		  else k=0;
		  ithterm_equal_newterm*=k;
		}
	      if(ithterm_equal_newterm == 1)uniqueind=0;
	    }
	  if(uniqueind==0)break;
	}
    }
  return(uniqueind);
}


int DSA_PACK_exportresult(double *oldrss,double newrss,int msize,int nvarX,int *currentmodel,
			  int *nextmodel,int maxsize,double *bestarss,int *gCOUNT,int IP,
			  int nglmcall,int **bestmodels,int delind,double *currentarss)
{
  int sflag;
  int i,j;

  sflag=0;
  if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh) {
    Rprintf("\noldrss:%e - newrss:%e - diff: %e - eps: %e",*oldrss,newrss,*oldrss - newrss,
	    DSA_PACK_epsilonCompare);
  }


  if(newrss==R_PosInf)
    {
      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)
	Rprintf("\nNone of the candidate models in this step can be fitted.\n");
    }
  else
    {
      if((*oldrss - newrss)>DSA_PACK_epsilonCompare)
	{
	  sflag = 1;
	  (*oldrss) = newrss;
	  for(i=0;i<maxsize;i++)for(j=0;j<nvarX;j++)currentmodel[j*maxsize+i] = nextmodel[j*maxsize+i];

	  currentarss[IP-1]=*oldrss;
	  if(delind==1)
	    {
	      bestarss[IP-1]=currentarss[IP-1];
	      memset(bestmodels[IP-2], 0, sizeof(int)*maxsize*nvarX);
	      for(i=0;i<maxsize;i++)
		for(j=0;j<nvarX;j++)
		  bestmodels[IP-2][j*maxsize+i] = currentmodel[j*maxsize+i];
	      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)
		{
		  Rprintf("\nModel accepted %i - bestarss: %lf\n",nglmcall,bestarss[IP-1]);
		  
		  
		  for(i=0;i<msize;i++)
		    {
		      Rprintf("%d \t ",i);
		      for(j=0;j<nvarX;j++)
			if(bestmodels[IP-2][j*maxsize+i] > 0)
			  Rprintf("%d^%d ",j,bestmodels[IP-2][j*maxsize+i]);
		      Rprintf("\n");
		    }
		}
	    }
	}

      gCOUNT[0] = msize;
      gCOUNT[2] = nglmcall;

      /* Update Best Index Set */
      if((delind==0) && ( (bestarss[IP-1]-currentarss[IP-1])>DSA_PACK_epsilonCompare))
	{
	  bestarss[IP-1]=currentarss[IP-1];
	  memset(bestmodels[IP-2], 0, sizeof(int)*maxsize*nvarX);
	  for(i=0;i<maxsize;i++)for(j=0;j<nvarX;j++)bestmodels[IP-2][j*maxsize+i] = currentmodel[j*maxsize+i];
	  if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)
	    {
	      Rprintf("\nModel accepted %i - bestarss: %lf\n",nglmcall,bestarss[IP-1]);
	      for(i=0;i<msize;i++)
		{
		  Rprintf("%d \t ",i);
		  for(j=0;j<nvarX;j++)if(bestmodels[IP-2][j*maxsize+i] > 0)
		    Rprintf("%d^%d ",j,bestmodels[IP-2][j*maxsize+i]);
		  Rprintf("\n");
		}
	    }
	}
    }
  return(sflag);
}

void  DSA_PACK_get_tempXdata(int samplesize,double *tempXdata,int msize,int maxsize,int *tempmodel, int nvarX,double *Xdata, double *tempterm,int *lastXdesign)
{
  int size,sumpower,singlepj,d,m,onevarmissing,updateX;
  double carrot;

  for(size=0;size<msize;size++)
    {
      updateX=DSA_PACK_update_design_matrix(tempmodel,lastXdesign,nvarX,maxsize,size);
/*       if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)if(updateX==1)Rprintf("\nWould need to update term %li\n",size+1); */

      if(updateX==1) 
	{ 
	  sumpower=DSA_PACK_rowsumints(tempmodel,maxsize,nvarX,size);

	  /*Just a main term */
	  if(sumpower==1)
	    {  
	      for(d=0;d<nvarX;d++)
		{
		  if(tempmodel[d*maxsize+size] == 1)
		    {
		      memcpy(tempXdata+(size+1)*samplesize,Xdata+d*samplesize,sizeof(double)*samplesize);
		      break;
		    }
		}
	    }

	  /*Interactions or main term to the power>1*/
	  if(sumpower > 1)
	    {
	      memcpy(tempterm,tempXdata,sizeof(double)*samplesize); /*set to 1.0*/
	      for(d=0;d<nvarX;d++)
		{
		  singlepj = tempmodel[d*maxsize+size];
		  if(singlepj > 1)
		    {
		      for(m=0;m<samplesize;m++)
			{
			  carrot = R_pow_di(Xdata[m+d*samplesize],singlepj);
			  tempterm[m] = tempterm[m] * carrot;
			}
		    }
		  if(singlepj == 1)for(m=0;m<samplesize;m++)tempterm[m] = tempterm[m] * Xdata[m+d*samplesize];
		}
	      memcpy(tempXdata+(size+1)*samplesize,tempterm,sizeof(double)*samplesize);
	    }
	} 
    }
}

void  DSA_PACK_get_Ytempdata_glmWT(int binind, int samplesize,double *tempYdata,double *Ydata,double *tempXdata,int CENSOR,double *glmWT,double *WTdata)
{
  int m,d;

  if(binind>1)memcpy(tempYdata,Ydata,sizeof(double)*samplesize*binind);
  else memcpy(tempYdata,Ydata,sizeof(double)*samplesize);
  if(CENSOR!=0)memcpy(glmWT,WTdata,sizeof(double)*samplesize);

  for(m=0;m<samplesize;m++)
    {
      *(tempXdata+m)=1.0;
      if(CENSOR==0)*(glmWT+m)=1.0;
    }

}

int  DSA_PACK_update_design_matrix(int *tempmodel,int *lastXdesign,int nvarX,int maxsize,int rowtocompare)
{
  int i,res=0;
  
  for(i=0;i<nvarX;i++)if(*(tempmodel+rowtocompare+maxsize*i)!=*(lastXdesign+rowtocompare+maxsize*i))
    {
      res=1;
      *(lastXdesign+rowtocompare+maxsize*i)=*(tempmodel+rowtocompare+maxsize*i);
    }

  return(res);
}

