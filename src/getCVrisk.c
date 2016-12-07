#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include "dsa.h"


void DSA_PACK_getaverageCVrisks(int *gCOUNT,int *gCONST,double *Xtrain,double *Ytrain,
				double *Xval,double *Yval,int **bestmodels,double *tempXtrain,
				double *tempYtrain,double *tempXval,
				double *tempYval,double *temptermtrain,double *temptermval,double *PARA,
				double *RES, double *glmWT,
				double *CVrisksfold,double *WTtrain,double *WTval,int nforced,
				int *tempmodel,double *binWTtrain,double *binWTval,
				int *modelfit_work_i,double *modelfit_work_d,int *lastXdesign)
{
  int IP,msize,maxsize,ntrain,nval,nvarX,fold,kk,i,nglmcall,CENSOR,orderint,maxorderint,NnonNA,binind,*Pmodel,mtype,M,ii,usetree,m,d;
  double temprss;

  /* INITIALIZATION */
  msize = gCOUNT[0];
  orderint = gCOUNT[1];
  nglmcall = gCOUNT[2];
  fold = gCOUNT[3];
  ntrain = gCOUNT[4];
  nval = gCOUNT[5];

  maxsize = gCONST[0];
  nvarX = gCONST[5];
  CENSOR = gCONST[6];
  maxorderint = gCONST[1];
  binind = gCONST[9];
  usetree = gCONST[11];
  mtype = min(binind + 1,3);
  if(binind>1)M=binind;
  else M=1;
  
  if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)
    {
      Rprintf("orderint: %ld and maxorderint: %ld\n",orderint,maxorderint);
      Rprintf("msize is: %li\n",(msize+1));
    }

  /* GET CV RISK FOR EACH MODEL SIZE */
  for(kk=0;kk<nforced;kk++)CVrisksfold[kk]=R_PosInf;
  for(kk=msize+1;kk<maxsize+1;kk++)CVrisksfold[kk]=R_PosInf;

 
  for(kk=nforced;kk<(msize+1);kk++) {
    /*FIT MODELS ON TRAINING SET*/
    IP = kk + 1;
    temprss=R_PosInf;
    if(kk==0)
      {
	Pmodel=tempmodel;
	memset(Pmodel, 0, sizeof(int)*maxsize*nvarX);
      }
    else Pmodel=bestmodels[kk-1];
    
    DSA_PACK_get_Ytempdata_glmWT(binind,ntrain,tempYtrain,Ytrain,tempXtrain,CENSOR,glmWT,WTtrain);
    memset(lastXdesign, 0, sizeof(int)*maxsize*nvarX);
    DSA_PACK_eval_model(maxsize,kk,Pmodel, -1, Pmodel, nvarX, ntrain,
			tempXtrain, Xtrain, temptermtrain, CENSOR, glmWT, WTtrain, IP,
			tempYtrain, Ytrain, PARA, RES,
			&nglmcall, &temprss, Pmodel, 1, binind, binWTtrain, 1 /* the always doit flag. */,usetree,
			modelfit_work_i,modelfit_work_d,lastXdesign);
    /*Get DATA FROM VALIDATION SET*/
    DSA_PACK_get_Ytempdata_glmWT(binind,nval,tempYval,Yval,tempXval,CENSOR,glmWT,WTval);
    memset(lastXdesign, 0, sizeof(int)*maxsize*nvarX);
    DSA_PACK_get_tempXdata(nval,tempXval,kk,maxsize,Pmodel,nvarX,Xval,temptermval,lastXdesign);

    /*GET NnonNA for validation set*/
    NnonNA=nval;
    for(m=0;m<nval;m++)
      { 
 	if(ISNA(*(tempYval+m))!=0)*(glmWT+m)=0; 
 	for(d=1;d<kk+1;d++)
	  if(ISNA(*(tempXval+m+d*nval))!=0)*(glmWT+m)=0; 
	if(abs(*(glmWT+m))<DSA_PACK_epsilonCompare)
	  NnonNA--;
      } 
    /* GET CV ON VALIDATION SET */
    if(NnonNA==0 || temprss==R_PosInf)error("\nThe DSA was interrupted because the data splits cannot be used to compute the cross-validated risk for one of the models selected (there is no complete observation in one of the test sets or non-convergence of the model fitting procedure on one of the training sets). Try another data split.\n");

    DSA_PACK_evaluate_loss(mtype, tempXval,tempYval, nval, (kk + 1) , 
			   M, maxsize+1, glmWT, binWTval, PARA, &temprss);
      
    temprss/=NnonNA;
    if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)
      {
	Rprintf("\nCoef on training:\n");
	if(binind>1)for(ii = 0; ii < binind; ii++)
	  {
	    DSA_PACK_printmatrix(PARA+ii*(kk+1),kk+1,1);
	    Rprintf("\n");
	  }
	else DSA_PACK_printmatrix(PARA,kk+1,1);
	Rprintf("\nNnonNA: %li\n",NnonNA);
      }


    if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)Rprintf("\nFor k=%ld, cvrisk=%lf \n",kk,temprss);
    CVrisksfold[kk] = temprss;

  }
  
  gCOUNT[2] = nglmcall;
}
