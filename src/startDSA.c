/*****************************************************************************
- runs on learning or training set 
- calls the deletion, substitution, and addition routines
- fits weighted regressions

***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>

/**
 #include <nag.h>
 #include <nagg02.h>
**/

#include <R.h>
#include <R_ext/Utils.h>
#include "dsa.h"

int DSA_PACK_userMlevel = -1;
double DSA_PACK_epsilonCompare = MY_ZERO;
unsigned long int* MODEL_TO_FIT = NULL;
p_node_t* TREE = NULL;
int R_DSA_DEBUG_LEVEL = 10;

void DSA_PACK_startAlgo(double *glmWT,double *PARA,
			double *RES,int *gCOUNT,int *gCONST,double *bestarss,
			double *Ydata,double *Xdata,int **bestmodels,int *currentmodel,
			int *tempmodel,int *nextmodel,int *newterm,int *ithterm,double *tempterm,
			double *tempXdata,double *tempYdata,double *WTdata,int trainFlag,
			double *currentarss,int *nforced, int *forcedterms,double *binWTdata,
			int *modelfit_work_i,double *modelfit_work_d,int *lastXdesign)
{
  int IP,msize,maxsize,nvarX,samplesize,nval,nexttoadd,i,j,bestsize,nglmcall,ii,CENSOR,usetree;
  int NnonNA,deletion,subsmade,extend,FLAG,binind,Dmove,Smove;
  double oldarss,newarss,temparss,keepBeta;

  /** JHB
      
   Nag_IncludeMean MEAN=Nag_MeanZero;
   Nag_Link LINK=Nag_Logistic;
   Boolean SVD; 
   NagError fail; 

  **/
  
  /* INITIALIZATION */
  IP = 0;
  msize = gCOUNT[0];
  nglmcall = gCOUNT[2];
  nval = gCOUNT[5];

  maxsize = gCONST[0];
  bestsize = gCONST[7];
  nvarX = gCONST[5];
  CENSOR = gCONST[6];
  binind= gCONST[9];
  usetree = gCONST[11];

  if(trainFlag == 1)
    samplesize = gCOUNT[4];
  else 
    samplesize = gCONST[4];

  Dmove = gCONST[13];  
  Smove = gCONST[14];  

  for(i=0;i<(maxsize+1);i++)bestarss[i]=currentarss[i]=R_PosInf;
  DSA_PACK_get_Ytempdata_glmWT(binind,samplesize,tempYdata,Ydata,tempXdata,CENSOR,glmWT,WTdata);
  memset(lastXdesign, 0, sizeof(int)*maxsize*nvarX);

  if(*nforced==0)
    {
      /* INTERCEPT ONLY MODEL */
      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)Rprintf("\nFit intercept model");
      IP = 1;
      msize=0;
      memset(tempmodel, 0, sizeof(int)*maxsize*nvarX);
      oldarss=R_PosInf;
      
      DSA_PACK_eval_model(maxsize, msize, nextmodel, -1, tempmodel, nvarX, samplesize,
			  tempXdata, Xdata, tempterm, CENSOR, glmWT, WTdata, IP,tempYdata, Ydata, PARA,RES,
			  &nglmcall, &oldarss, currentmodel, 1, binind, binWTdata,0,usetree,
			  modelfit_work_i,modelfit_work_d,lastXdesign);
      bestarss[IP-1] = oldarss;
      currentarss[IP-1] = oldarss;
    }
  else
    {
      /* FORCED TERMS */
      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)Rprintf("\nFit forced terms only");
      IP = *nforced+1;
      msize=*nforced;
      memset(tempmodel, 0, sizeof(int)*maxsize*nvarX);
      for(ii=0;ii<*nforced;ii++)for(j=0;j<nvarX;j++)*(tempmodel+ii+j*maxsize)=*(forcedterms+ii+j**nforced);
      oldarss=R_PosInf;

      DSA_PACK_eval_model(maxsize,msize, nextmodel, -1, tempmodel, nvarX, samplesize,
			  tempXdata, Xdata, tempterm, CENSOR, glmWT, WTdata, IP,tempYdata, Ydata, PARA,RES,
			  &nglmcall, &oldarss, currentmodel, 1, binind, binWTdata,0,usetree,
			  modelfit_work_i,modelfit_work_d,lastXdesign);

      if(oldarss==R_PosInf)error("\nThe DSA was interrupted because the base model or model selected cannot be fitted on one of the training sets. Try another data split or change the base model.\n");

      bestarss[IP-1] = oldarss;
      currentarss[IP-1] = oldarss;
      for(ii=0;ii<IP-1;ii++)
	{
	  bestarss[ii]=R_PosInf;
	  currentarss[ii]=R_PosInf;
	}
    }

  if(*nforced==0)
    {
      /* SELECT BEST ONE-TERM  MODEL */
      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)Rprintf("\nFit best one-term model");
      newarss=R_PosInf;
      nexttoadd=-1;
      keepBeta=R_PosInf;
      IP=2;
      msize=1;
      for(j=0;j<nvarX;j++)
	{
	  memset(tempmodel, 0, sizeof(int)*maxsize*nvarX);
	  *(tempmodel+j*maxsize)=1;
	  temparss=R_PosInf;

	  DSA_PACK_eval_model(maxsize,msize, nextmodel, -1, tempmodel, nvarX, samplesize,
			      tempXdata, Xdata, tempterm, CENSOR, glmWT, WTdata, IP,tempYdata, Ydata, PARA,RES,
			      &nglmcall, &temparss, currentmodel, 1, binind, binWTdata,0,usetree,
			      modelfit_work_i,modelfit_work_d,lastXdesign);

	  if((newarss-temparss)>DSA_PACK_epsilonCompare)
	    {
	      newarss=temparss;
	      nexttoadd=j;
	      keepBeta = PARA[1];
	    }
	}
      if(nexttoadd==-1)error("\nThe DSA was interrupted because the base model or model selected cannot be fitted on one of the training sets. Try another data split or change the base model.\n");
      oldarss=newarss;
    }


  /* START DSA - INITIALIZE BESTARSS, BESTINDEXSET, AND CURRINDEXSET */
  if(*nforced==0)
    {
      bestarss[IP-1] = oldarss;
      currentarss[IP-1] = oldarss;
      memset(currentmodel, 0, sizeof(int)*maxsize*nvarX);
      currentmodel[nexttoadd*maxsize+0] = 1;
      memset(bestmodels[0], 0, sizeof(int)*maxsize*nvarX);
      bestmodels[0][nexttoadd*maxsize+0] = 1;
      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)Rprintf("\nFirst best model: ARSS=%lf",oldarss);
    }
  else
    {
      memset(bestmodels[IP-2], 0, sizeof(int)*maxsize*nvarX);
      for(ii=0;ii<*nforced;ii++)
	for(j=0;j<nvarX;j++)
	  *(bestmodels[IP-2]+ii+j*maxsize)=*(currentmodel+ii+j*maxsize);
      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)Rprintf("\nFirst model: ARSS=%lf",oldarss);
    }

  if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)Rprintf("\n");
  for(i=0;i<msize;i++)
    {
      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)Rprintf("%d \t ",i);
      for(j=0;j<nvarX;j++)
	if(currentmodel[j*maxsize+i] > 0)
	  if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)Rprintf("%d^%d ",j,currentmodel[j*maxsize+i]);
      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)Rprintf("\n");
    }
  if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)Rprintf("\n");

  gCOUNT[0]=msize;
  gCOUNT[2]=nglmcall;
  FLAG = 1;

  if(gCOUNT[0]<maxsize) /*No model selection if the maximum size is already reached before starting, that means the user just cares about CV risks computation*/
    {
      if(DSA_PACK_userMlevel>=DSA_PACK_Mlow)Rprintf("\n");
      if(DSA_PACK_userMlevel>=DSA_PACK_Mlow)Rprintf("Start! \n");
      if(DSA_PACK_userMlevel>=DSA_PACK_Mlow)Rprintf("Running algorithm ..... \n");
      if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)Rprintf("\nNumber of terms in model: %i",gCOUNT[0]);
      while(FLAG==1 && gCOUNT[0]<=maxsize && (gCOUNT[0]+1)<samplesize)
	{
	  deletion=0;
	  subsmade=0;
	  extend=1;

	  if(gCOUNT[0]>max(1,*nforced) && Dmove==1)
	    {
	      /* D step */  
	      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\nDeletion step\n");
	      deletion = DSA_PACK_Dstep(gCOUNT,gCONST,bestarss,bestmodels,tempmodel,currentmodel,
					nextmodel,newterm,Xdata,tempXdata,Ydata,tempYdata,tempterm,
					glmWT,PARA,RES,WTdata,trainFlag,currentarss,*nforced,binWTdata
					,modelfit_work_i,modelfit_work_d,lastXdesign);

	      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)
		Rprintf("\nNumber of terms in model: %i, deletion=%i",gCOUNT[0],deletion);
	      R_CheckUserInterrupt();
	    }
	  if(bestarss[gCOUNT[0]] < DSA_PACK_epsilonCompare)
	    {
	      FLAG=0;
	      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)
		Rprintf("\nStop DSA after del. because the last model accepted fits the data very well (loss of 0).\n");
	    }
	  if(deletion==0 && FLAG==1)
	    {
	      if(Smove==1)
		{
		  /* S step */
		  if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\nSubstitution step\n");
		  subsmade = DSA_PACK_Sstep(gCOUNT,gCONST,bestarss,bestmodels,tempmodel,currentmodel,
					    nextmodel,newterm,Xdata,tempXdata,Ydata,tempYdata,tempterm,
					    glmWT,PARA,RES,ithterm,WTdata,trainFlag,currentarss,
					    *nforced,binWTdata,modelfit_work_i,modelfit_work_d,lastXdesign);
		  if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)
		    Rprintf("\nNumber of terms in model: %i, subsmade=%i",gCOUNT[0],subsmade);
		  R_CheckUserInterrupt();
		}
	      if(bestarss[gCOUNT[0]] < DSA_PACK_epsilonCompare)
		{
		  FLAG=0;
		  if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)
		    Rprintf("\nStop DSA after sub. because the last model accepted fits the data very well (loss of 0).\n");
		}
	      if(subsmade==0 && FLAG==1)
		{
		  if(gCOUNT[0]<maxsize)
		    {
		      /* A step */
		      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)Rprintf("\nAddition step\n");
		      extend = DSA_PACK_Astep(gCOUNT,gCONST,bestarss,bestmodels,tempmodel,currentmodel,
					      nextmodel,newterm,Xdata,tempXdata,Ydata,tempYdata,tempterm,
					      glmWT,PARA,RES,ithterm,WTdata,trainFlag,currentarss,binWTdata,
					      modelfit_work_i,modelfit_work_d,lastXdesign);	 
		      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)
			Rprintf("\nNumber of terms in model: %i, extend=%i",gCOUNT[0],extend);
		      R_CheckUserInterrupt();
		    }
		  else if(gCOUNT[0] == maxsize)
		    {
		      extend = 0;
		      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)
			Rprintf("\nNot allowed to add anymnore.");
		    }
		  if(extend == 0)
		    {
		      FLAG=0;
		      if(DSA_PACK_userMlevel>=DSA_PACK_Mhigh)
			Rprintf("\nStop DSA after add. because the last model accepted fits the data very well (loss of 0).\n");
		    }
		}
	    }
	}    /* end of while loop */
    } /* end of if(gCOUNT[0]<maxsize) */
  if(DSA_PACK_userMlevel>=DSA_PACK_Mmedium)
    Rprintf("\nAfter DSA in startalgo: globalcount[1]=msize: %i\n\n",gCOUNT[0]);

  if(trainFlag == 0)
    {
      /* OPTIMAL MODEL IS PRINTED TO FORMATOUT */
      IP = bestsize + 1;
      oldarss=R_PosInf;
      
      if(bestsize>0) {
	DSA_PACK_eval_model(maxsize, bestsize, tempmodel, -1, bestmodels[bestsize-1], nvarX, samplesize,
			    tempXdata, Xdata, tempterm, CENSOR, glmWT, WTdata, IP,tempYdata, Ydata, PARA, RES,
			    &nglmcall, &oldarss, currentmodel, 1, binind, binWTdata, 1,usetree,
			    modelfit_work_i,modelfit_work_d,lastXdesign);

      }
      else {
	memset(nextmodel, 0, sizeof(int)*maxsize*nvarX);
	DSA_PACK_eval_model(maxsize, bestsize, tempmodel, -1, nextmodel, nvarX, samplesize,
			    tempXdata, Xdata, tempterm, CENSOR, glmWT, WTdata, IP,tempYdata, Ydata, PARA, RES,
			    &nglmcall, &oldarss, currentmodel, 1, binind, binWTdata, 1,usetree,
			    modelfit_work_i,modelfit_work_d,lastXdesign);

      }
      if(DSA_PACK_userMlevel>=DSA_PACK_Mlow)
	{
	  Rprintf("OPTIMAL MODEL: \n");
	  Rprintf("Final ARSS: %lf\n",oldarss);
	  Rprintf("terms in model w intercept %d \n",IP);
	  Rprintf("terms in model wo intercept %d \n",bestsize);
	  Rprintf("Estimates of coefficients: \n");
	  if(binind>1)for(ii = 0; ii < binind; ii++)
	    {
	      DSA_PACK_printmatrix(PARA+ii*IP,IP,1);
	      Rprintf("\n");
	    }
	  else DSA_PACK_printmatrix(PARA,IP,1);
	  Rprintf("\n");
	  Rprintf("bestmodels for size %d \n",bestsize);
	  if(bestsize>0)
	    {
	      DSA_PACK_printmatrixl(bestmodels[bestsize-1],maxsize,nvarX);
	      Rprintf("\n");
	      Rprintf("Terms in Final Model: \n");
	      for(ii=0;ii<bestsize;ii++)
		{
		  Rprintf("%d \t ",ii);
		  for(i=0;i<nvarX;i++)
		    {
		      j = bestmodels[bestsize-1][i*maxsize+ii];
		      if(j > 0)Rprintf("%d^%d ",i,j);
		    }
		  Rprintf("\n");
		}
	      Rprintf("\n");
	    }
	}


    } /* end of printing if this is DSA on learning set */ 
 

}
