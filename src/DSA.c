/*****************************************************************************
THIS FILE CONTAINS THE C FUNCTION CALLED IN R

   - runs DSA on each training set
   - uses empirical RSS to form models
   - uses cross-validation to pick size and order of interaction of final model
   - runs swap-sub. moves if complexity measure is exceeded
   - incorporates rules on what variables cannot interact
   - fits weighted regression models
   - automatically fills cross-validated risks grid for maxsize and orderint

  
--Counters and Constants
gCOUNT:
0 number of terms in the last model accepted - intercept excluded (msize)
1 current maximum order of interaction allowed (orderint)
2 number of calls to glm routine (nglmcall)
3 fold number (fold)
4 size of training set (ntrain)
5 size of validation set (nval)

gCONST:
0 max number of terms allowed - intercept excluded (maxsize)
1 max interaction order (maxorderint)
2 max sum of powers in each term (maxsumofpow)
3 number of data splits (nfolds)
4 total number of observations in learning set (nlearn)
5 number of candidate variables (nvarX)
6 indicator of weight use (CENSOR)
7 optimal size chosen via cross-validation (bestsize)
8 optimal order of interaction chosen via cross-validation (bestorderint)
9 equal to 0 if linear model, 1 if logistic and >1 if multinomial with the 
  value indicating the number of categories MINUS 1 (binind)
10 indicator of CVDSA (vs no CVDSA) use
11 indicator of TREE use
12 nsplits
13 Dmove - indicator that Deletion moves are allowed
14 Smove - indicator that Substitution moves are allowed
***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>
#include "dsa.h"
#include "dsa_tree.h"


/* define imports from modelUtils */
int* (*model_allocate_i)(int,int,int,int);
double* (*model_allocate_d)(int,int,int,int);

void DSA_PACK_Rentry(int *gCONST,double *Xlearn,double *Ylearn,double *CVrisks,double *PARA,
		     int *vsplits,double *WTtrainvallearn,int *nforced,int *forcedterms,
		     int *maxntrain, int *maxnval,int *bestmodelsspace,double *binWTtrainvallearn)
{
  int i,j,k,maxsize,orderint,maxsumofpow,maxorderint,bestsize,bestorderint,nfolds,fold,totalglmcalls,nlearn,nvarX,CENSOR,**bestmodels,ntrain,nval,binind,*trainindex,*valindex,*currentmodel,*tempmodel,*nextmodel,*newterm,*ithterm,*gCOUNT,CVDSA,usetree,nsplits,mtype,M,*lastXdesign,*modelfit_work_i; 

  double grandmin,meanCVrisk,*CVrisksfold,*WK,*glmWT,*H,*RES,*Q,*PP,*bestarss,*currentarss,*Ytrain,*Xtrain,*Yval,*Xval,*tempXval,*tempYval,*temptermval,*temptermlearntrain,*tempXlearntrain,*tempYlearntrain,*WTlearntrain,*WTval,*binWTlearntrain,*binWTval,*modelfit_work_d;

  /* INITIALIZATION */
  maxsize = gCONST[0] ;        /* maxsize: max number of terms allowed - intercept excluded */
  maxorderint = gCONST[1] ;    /* maxorderint: max interaction order */
  maxsumofpow = gCONST[2] ;    /* maxsumofpow: max sum of powers in each term */
  nfolds = gCONST[3] ;         /* nfolds: number of data splits */
  nlearn = gCONST[4] ;         /* nlearn: total number of observations in learning set */
  nvarX = gCONST[5] ;          /* nvarX: number of candidate variables*/
  CENSOR = gCONST[6] ;         /* CENSOR: indicator of weight use */
                               /* gCONST[7] = bestsize: optimal size chosen via cross-validation */
                               /* gCONST[8] = bestorderint: optimal order of interaction chosen via cross-validation */
  binind = gCONST[9] ;         /* binind: indicator of the outcome being continuous, binary or categorical>2*/
  CVDSA = gCONST[10] ;         /* CVDSA: indicator of CVDSA (or no CV) */
  usetree = gCONST[11] ;       /* indicator of TREE use */
  nsplits = gCONST[12] ;       /* nsplits from R */
                               /* gCONST[13] = indicator that Deletion moves are allowed */
                               /* gCONST[14] = indicator that Substitution moves are allowed */

  totalglmcalls = 0;
  mtype=min(binind + 1,3);
  if(binind>1)
    M=binind;
  else 
    M=1;

  if(usetree==1)MODEL_TO_FIT = (unsigned long int* )R_alloc(maxsize, sizeof(unsigned long int));
  /* MEMORY ALLOCATION */ /* NEED TO MAKE SURE TO INITIALIZE IF NECESSARY */
  trainindex=(int *)R_alloc(*maxntrain,sizeof(int));
  valindex=(int *)R_alloc(*maxnval,sizeof(int));
  CVrisksfold=(double *)R_alloc((maxsize+1),sizeof(double));
  glmWT=(double *)R_alloc(nlearn,sizeof(double));
  RES=(double *)R_alloc(nlearn,sizeof(double));

  bestarss=(double *)R_alloc(maxsize+1,sizeof(double));
  currentarss=(double *)R_alloc(maxsize+1,sizeof(double));
  if(binind>1)Ytrain=(double *)R_alloc(*maxntrain*binind,sizeof(double));
  else Ytrain=(double *)R_alloc(*maxntrain,sizeof(double));
  Xtrain=(double *)R_alloc(*maxntrain*nvarX,sizeof(double));
  memset(Xtrain, 0, sizeof(double)*(*maxntrain)*nvarX); /** JHB valgrind. **/

  if(binind>1)Yval=(double *)R_alloc(*maxnval*binind,sizeof(double));
  else Yval=(double *)R_alloc(*maxnval,sizeof(double));

  Xval=(double *)R_alloc(*maxnval*nvarX,sizeof(double));
  memset(Xval, 0, sizeof(double)*(*maxnval) * nvarX); /** JHB valgrind. **/

  tempXval=(double *)R_alloc(*maxnval*(maxsize+1),sizeof(double)); 
  memset(tempXval, 0, sizeof(double)*(*maxnval) * (maxsize+1)); /** JHB valgrind. **/

  if(binind>1)tempYval=(double *)R_alloc(*maxnval*binind,sizeof(double));
  else tempYval=(double *)R_alloc(*maxnval,sizeof(double));

  temptermval=(double *)R_alloc(*maxnval,sizeof(double));
  temptermlearntrain=(double *)R_alloc(nlearn,sizeof(double));

  tempXlearntrain=(double *)R_alloc(nlearn*(maxsize+1),sizeof(double));
  memset(tempXlearntrain, 0, sizeof(double)*nlearn*(maxsize+1)); /** JHB valgrind. **/

  if(binind>1)tempYlearntrain=(double *)R_alloc(nlearn*binind,sizeof(double));
  else tempYlearntrain=(double *)R_alloc(nlearn,sizeof(double));
  bestmodels=(int **)R_alloc(maxsize,sizeof(int *));
  for(i=0;i<maxsize;i++)*(bestmodels+i) = bestmodelsspace+i*(maxsize*nvarX);
  currentmodel=(int *)R_alloc(maxsize*nvarX,sizeof(int));
  tempmodel=(int *)R_alloc(maxsize*nvarX,sizeof(int));
  lastXdesign=(int *)R_alloc(maxsize*nvarX,sizeof(int));
  nextmodel=(int *)R_alloc(maxsize*nvarX,sizeof(int));
  newterm=(int *)R_alloc(nvarX,sizeof(int));
  ithterm=(int *)R_alloc(nvarX,sizeof(int));
  WTlearntrain=(double *)R_alloc(nlearn,sizeof(double));
  WTval=(double *)R_alloc(*maxnval,sizeof(double));
  binWTlearntrain=(double *)R_alloc(nlearn,sizeof(double));
  binWTval=(double *)R_alloc(*maxnval,sizeof(double));
  gCOUNT=(int *)R_alloc(6,sizeof(int));

  model_allocate_d = (double*(*)(int,int,int,int)) R_GetCCallable("modelUtils", "model_allocate_double");
  model_allocate_i = (int*(*)(int,int,int,int)) R_GetCCallable("modelUtils", "model_allocate_int");
  modelfit_work_d=(*model_allocate_d)(mtype,nlearn,maxsize+1,M);
  modelfit_work_i=(*model_allocate_i)(mtype,nlearn,maxsize+1,M);

  if(CVDSA==1)
    {
      /* START DATA-ADAPTIVE SELECTION OF BEST SIZE AND ORDER OF INTERACTION */
      for(fold = 0; fold < nfolds; fold++)
	{
	  /* ALLOCATE MEMORY FOR A TREE FOR THIS FOLD ONLY */
	  if(usetree==1)
	    {
	      TREE = make_new_node(1, -1); /** the root node. **/ 
	      if(DSA_PACK_userMlevel >= DSA_PACK_Mbase)Rprintf("\nTree is initialized for training set %i: name = %lu, rss = %g \n", fold+1, TREE->name, TREE->rss);    
	    }

	  /* COMPUTE FOR THE CORRESPONDING SPLIT THE SIZE OF THE TRAINING AND VALIDATION SET */
	  nval=0;
	  for(i=0;i<nlearn;i++)if(*(vsplits+fold+i*nfolds)==1)nval++;
	  ntrain=nlearn-nval;

	  /*IDENTIFIES WHICH OBSERVATION IS IN THE TRAIN AND VAL SETS FOR THE CURRENT SPLIT OF THE DATA*/
	  memset(trainindex, 0, sizeof(int)*ntrain);
	  memset(valindex, 0, sizeof(int)*nval);
	  j=0;
	  k=0;
	  for(i=0; i<nlearn; i++)
	    {
	      if(vsplits[i*nfolds+fold] == 0)
		{
		  trainindex[j] = i;
		  j++;
		}
	      else
		{
		  valindex[k] = i;
		  k++;
		}
	    }

	  /* CREATE TWO MATRICES CONTAINING THE TRAINING DATA: Y AND X */
	  for(i = 0; i < ntrain; i++)
	    {
	      if(binind>1)for(j = 0; j < binind; j++)Ytrain[j*ntrain+i] = Ylearn[j*nlearn+trainindex[i]];
	      else Ytrain[i] = Ylearn[trainindex[i]];

	      WTlearntrain[i] = WTtrainvallearn[trainindex[i]*(nfolds+1) + fold];
	      binWTlearntrain[i] = binWTtrainvallearn[trainindex[i]];
	      for(j = 0; j < nvarX; j++)Xtrain[i+j*ntrain] = Xlearn[trainindex[i]+ j*nlearn];
	    }
	  /* CREATE TWO MATRICES CONTAINING THE VALIDATION DATA: Y AND X */
	  for(i = 0; i < nval; i++)
	    {
	      if(binind>1)for(j = 0; j < binind; j++)Yval[j*nval+i] = Ylearn[j*nlearn+valindex[i]];
	      else Yval[i] = Ylearn[valindex[i]];

	      WTval[i] = WTtrainvallearn[valindex[i]*(nfolds+1) + fold];
	      binWTval[i] = binWTtrainvallearn[valindex[i]];
	      for(j=0; j<nvarX; j++)Xval[i+j*nval] = Xlearn[valindex[i] + j*nlearn];
	    }

	  for(orderint = 1; orderint <= maxorderint; orderint++)
	    {
	      if(DSA_PACK_userMlevel >= DSA_PACK_Mbase)
		Rprintf("\n****** Model selection among models of maximum interaction order = %i based on split V%i ******",orderint,fold+1);

	      gCOUNT[0] = 0;              /* msize: number of terms in the current model - intercept excluded */
	      gCOUNT[1] = orderint;       /* orderint: current maximum order of interaction allowed */
	      gCOUNT[2] = 0;              /* nglmcall: number of calls to glm routine */
	      gCOUNT[3] = fold;           /* fold: fold number (fold) */ 
	      gCOUNT[4] = ntrain;         /* ntrain: size of training set */
	      gCOUNT[5] = nval;           /* nval: size of validation set */

	      if(DSA_PACK_userMlevel>=DSA_PACK_Mlow)
		{
		  Rprintf("\n");
		  Rprintf("\tMaximal number of terms in each model considered - intercept excluded = %d \n",
			  maxsize);
		  Rprintf("\tMaximal sum of powers in each term of each model considered = %d \n",maxsumofpow);
		  Rprintf("\tSize of the validation set (nval): %d \n",nval);
		  Rprintf("\tSize of the training set (ntrain): %d \n",ntrain);
		  Rprintf("\n");
		}

	      /* RUN THE DSA ALGORITHM */
	      if(DSA_PACK_userMlevel>=DSA_PACK_Mbase)
		Rprintf("Selecting models on the training set based on average risks.\n");

	      DSA_PACK_startAlgo(glmWT,PARA,RES,gCOUNT,gCONST,bestarss,Ytrain,
				 Xtrain,bestmodels,currentmodel,tempmodel,nextmodel,newterm,ithterm,
				 temptermlearntrain,tempXlearntrain,tempYlearntrain,
				 WTlearntrain,1,currentarss,
				 nforced,forcedterms,binWTlearntrain,modelfit_work_i,modelfit_work_d,
				 lastXdesign);

	      if(DSA_PACK_userMlevel>=DSA_PACK_Mlow)
		Rprintf("The largest model considered is of size %i.\n",gCOUNT[0]);
	  
	      totalglmcalls = totalglmcalls + gCOUNT[2];
	      if(DSA_PACK_userMlevel>=DSA_PACK_Mlow)
		Rprintf("The number of times glm was called on this training set is %d.\n",gCOUNT[2]);
	    
	      /* COMPUTE THE CROSS-VALIDATED RISKS */
	      if(DSA_PACK_userMlevel>=DSA_PACK_Mlow)
		Rprintf("Computing the average cross-validated risks of the models selected.\n");
	      gCOUNT[2] = 0; /* reset nglmcall */

	      DSA_PACK_getaverageCVrisks(gCOUNT,gCONST,Xtrain,Ytrain,Xval,Yval,bestmodels,tempXlearntrain,
					 tempYlearntrain,tempXval,tempYval,
					 temptermlearntrain,temptermval,
					 PARA,RES,glmWT,CVrisksfold,WTlearntrain,WTval,*nforced,
					 tempmodel,binWTlearntrain,binWTval,modelfit_work_i,modelfit_work_d,
					 lastXdesign);
				     
	      totalglmcalls = totalglmcalls + gCOUNT[2];
	      if(DSA_PACK_userMlevel>=DSA_PACK_Mlow)
		Rprintf("The number of times glm was called on this validation set is %d.\n",gCOUNT[2]);

	      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)
		Rprintf("The average cross-validated risks by split (maxsize by nfolds) are: \n");

	      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)
		DSA_PACK_printmatrix(CVrisksfold,maxsize+1,1);

	      /* COLLECT CV RISKS ACROSS VALIDATION SETS FOR EACH SIZE */
	      for(i=0;i<maxsize+1;i++)
		CVrisks[i*maxorderint+(orderint-1)] += CVrisksfold[i];

	    }/*end of for(orderint... ) loop*/

	  /* DELETE THE TREE FOR THIS FOLD */
	  if(usetree==1)
	    {
	      if(DSA_PACK_userMlevel >= DSA_PACK_Mbase)Rprintf("\nTree is going to be deleted for training set %i",fold+1);
	      delete_tree(TREE);
	      Free(TREE);
	      TREE = NULL;
	      if(DSA_PACK_userMlevel >= DSA_PACK_Mbase)Rprintf("\nTree is deleted for training set %i",fold+1);
	    }
	}/*end of for(fold =...)  loop*/

      /*Take average of all CVrisks by number of data splits*/
      for(orderint = 1; orderint <= maxorderint; orderint++)for(i=0;i<maxsize+1;i++)
	CVrisks[i*maxorderint+(orderint-1)]/=nfolds;

      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)
	Rprintf("\nThe cross-validated risks averaged over each split (maxorderint by maxsize) are: \n");
      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)DSA_PACK_printmatrix(CVrisks,maxorderint,maxsize+1);

      /* PICK BEST CV RISK TO PICK SELECT BEST SIZE AND ORDER OF INTERACTION */
      grandmin = DSA_PACK_getmin(CVrisks,maxorderint*(maxsize+1));
      bestsize=bestorderint=-1;
      for(i=0;i<maxorderint;i++)
	for(j=0;j<maxsize+1;j++)
	  if(CVrisks[j*maxorderint+i]!=R_PosInf)
	    if(abs(CVrisks[j*maxorderint+i] - grandmin) < DSA_PACK_epsilonCompare)
	      {
		bestsize = j ;
		bestorderint = i + 1;
		i=maxorderint+1;
		break;
	      }
      if(DSA_PACK_userMlevel>=DSA_PACK_Mlow)
	Rprintf("\nMinimum average cross-validated risk of %lf is reached for a model of size %i with a maximum order of interaction of %i.\n",grandmin,bestsize,bestorderint);
      gCONST[7]=bestsize;
      gCONST[8]=bestorderint;
    } /* end if(CVDSA==1)*/
  else
    {
      /* set bestsize to 0 in order not to get bestmodel=interecpt model in DSA.R (value corresponding with gCONST[7] not used when CVDSA=0 so it does not matter) */
      gCONST[7]=0;
      gCONST[8]=maxorderint;      
    }

  if(nsplits==1)
    {
      /* GET FINAL BEST MODEL BY RUNNING DSA ON LEARNING SET */
      if(DSA_PACK_userMlevel>=DSA_PACK_Mbase && CVDSA==1)
	Rprintf("\n****** Selecting the best model of size = %i and maximum order of interaction = %i on the learning set based on average residuals sum of squares ******\n",gCONST[7]+1,gCONST[8]);
      gCOUNT[2] = 0;   /* reset nglmcall */  
      for(i = 0; i < nlearn; i++)
	{    
	  WTlearntrain[i] = WTtrainvallearn[i*(nfolds+1) + nfolds];
	  binWTlearntrain[i] = binWTtrainvallearn[i];
	}

      /* ALLOCATE MEMORY FOR A TREE FOR THE LEARNING SET ONLY */
      if(usetree==1)
	{
	  TREE = make_new_node(1, -1); /** the root node. **/ 
	  if(DSA_PACK_userMlevel >= DSA_PACK_Mbase)Rprintf("\nTree is initialized for the learning set: name = %lu, rss = %g \n", TREE->name, TREE->rss);
	}
  
      DSA_PACK_startAlgo(glmWT,PARA,RES,gCOUNT,gCONST,bestarss,Ylearn,Xlearn,
			 bestmodels,currentmodel,tempmodel,nextmodel,newterm,ithterm,temptermlearntrain,
			 tempXlearntrain,tempYlearntrain,WTlearntrain,0,currentarss
			 ,nforced,forcedterms,
			 binWTlearntrain,modelfit_work_i,modelfit_work_d,lastXdesign);

      /* DELETE THE TREE FOR THE LEARNING SET */
      if(usetree==1)
	{
	  if(DSA_PACK_userMlevel >= DSA_PACK_Mbase)Rprintf("\nTree is going to be deleted for the learning set");
	  delete_tree(TREE);
	  Free(TREE);
	  TREE = NULL;
	  if(DSA_PACK_userMlevel >= DSA_PACK_Mbase)Rprintf("\nTree is deleted for the learning set");
	}

      /*Set to NAs all elements of bestmodelsspace corresponding with model sizes that were not considered in startAlgo due to a premature stop of the DSA (before reaching maxsize)*/
      for(i=gCOUNT[0]+1; i<maxsize; i++)for(k=0;k<maxsize;k++)for(j=0;j<nvarX;j++)*(*(bestmodels+i-1)+k+j*maxsize)=NA_INTEGER;

      totalglmcalls = totalglmcalls + gCOUNT[2];

      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)
	Rprintf("\nThe number of times glm was called on the learning set is %d.\n",gCOUNT[2]);
      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)
	Rprintf("The total number of calls to glm is %d.\n",totalglmcalls);
    }
  
  if(DSA_PACK_userMlevel>=DSA_PACK_Mlow)Rprintf("\nLeaving the C subroutine.\n");
}

void DSA_PACK_setMessagelevel(int *newMlevel)
{
  DSA_PACK_userMlevel=*newMlevel;
}

void DSA_PACK_getMessageLevel(int* currentLevel) 
{
  *currentLevel = DSA_PACK_userMlevel; 
}

void DSA_PACK_setEpsilonCompare(double* epsilon)
{
  DSA_PACK_epsilonCompare = *epsilon;
}

void DSA_PACK_getPrimes(int *prime)
{
  extern int PRIMES[];
  int i;

  for(i=0;i<1000;i++)*(prime+i)=*(PRIMES+i);
}

void DSA_PACK_DeleteTree()
{
  extern p_node_t* TREE;

  delete_tree(TREE);
}
