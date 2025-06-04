#ifndef __HAS_dsa__ 
#define __HAS_dsa__ 

#include "dsa_tree.h"

#define DSA_PACK_Mbase 0
#define DSA_PACK_Mlow 1
#define DSA_PACK_Mmedium 2
#define DSA_PACK_Mhigh 3
#define abs(a) (((a)>0)?(a):(-(a)))
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)>(b))?(b):(a))

/** 
 * some global constants which the user sets at package load time. 
 */
extern int DSA_PACK_userMlevel;
extern double DSA_PACK_epsilonCompare;

#ifndef DOUBLE_EPS
# define DOUBLE_EPS 2.2204460492503131e-16
#endif

/** snagged from R. **/
static const double MY_ZERO = (DOUBLE_EPS*100);

/** 
    we are going to make a couple of global variables here
    - an array to hold the composite numbers. 
    - The tree (which needs to be dealt with split wise. 
      so probably needs to go in DSA.c
#define __USE_TREE__
*/

extern unsigned long int* MODEL_TO_FIT;
extern p_node_t* TREE; 

void DSA_PACK_printmatrix(double *x, int nrow, int ncol);
void DSA_PACK_printmatrixl(int *x, int nrow, int ncol);
void DSA_PACK_printmatrix_to_file(double *x, int nrow, int ncol, char* f);

double DSA_PACK_getmin(double *x,int n);
int DSA_PACK_rowsumints(int *x,int nrow,int ncol,int whichrow);
int DSA_PACK_num_nonzero(int *x,int n);
int DSA_PACK_sumofints(int *x,int n);


void DSA_PACK_startAlgo(double *glmWT,double *PARA,
			double *RES,int *gCOUNT,int *gCONST,double *bestarss,
			double *Ydata,double *Xdata,int **bestmodels,int *currentmodel,
			int *tempmodel,int *nextmodel,int *newterm,int *ithterm,double *tempterm,
			double *tempXdata,double *tempYdata,double *WTdata,int trainFlag,
			double *currentarss,int *nforced, int *forcedterms,double *binWTdata,
			int *modelfit_work_i,double *modelfit_work_d,int *lastXdesign);		       

void DSA_PACK_getaverageCVrisks(int *gCOUNT,int *gCONST,double *Xtrain,double *Ytrain,
				double *Xval,double *Yval,int **bestmodels,double *tempXtrain,
				double *tempYtrain,double *tempXval,
				double *tempYval,double *temptermtrain,double *temptermval,double *PARA,
				double *RES, double *glmWT,
				double *CVrisksfold,double *WTtrain,double *WTval,int nforced,
				int *tempmodel,double *binWTtrain,double *binWTval,
				int *modelfit_work_i,double *modelfit_work_d,int *lastXdesign);

int DSA_PACK_Sstep(int *gCOUNT,int *gCONST,double *bestarss,int **bestmodels,int *tempmodel,
		   int *currentmodel,int *nextmodel,int *newterm,double *Xdata,double *tempXdata,
		   double *Ydata,double *tempYdata,double *tempterm,double *glmWT,double *PARA,
		   double *RES,int *ithterm,double *WTdata,int trainFlag,double *currentarss,int nforced,
		   double *binWTdata,int *modelfit_work_i,double *modelfit_work_d,int *lastXdesign);

int DSA_PACK_Astep(int *gCOUNT,int *gCONST,double *bestarss,int **bestmodels,int *tempmodel,int *currentmodel,
		   int *nextmodel,int *newterm,double *Xdata,double *tempXdata,double *Ydata,
		   double *tempYdata,
		   double *tempterm,double *glmWT,double *PARA,
		   double *RES,int *ithterm,double *WTdata,int trainFlag,double *currentarss,double *binWTdata,
		   int *modelfit_work_i,double *modelfit_work_d,int *lastXdesign);

int DSA_PACK_Dstep(int *gCOUNT,int *gCONST,double *bestarss,int **bestmodels,int *tempmodel,int *currentmodel,
		   int *nextmodel,int *newterm,double *Xdata,double *tempXdata,double *Ydata,
		   double *tempYdata,
		   double *tempterm,double *glmWT,double *PARA,
		   double *RES,double *WTdata,int trainFlag,double *currentarss,int nforced,double *binWTdata,
		   int *modelfit_work_i,double *modelfit_work_d,int *lastXdesign);

void  DSA_PACK_get_tempXdata(int samplesize,double *tempXdata,int msize,int maxsize,
			     int *tempmodel, int nvarX,double *Xdata, double *tempterm,
			     int *lastXdesign);
void  DSA_PACK_get_Ytempdata_glmWT(int binind, int samplesize,double *tempYdata,double *Ydata,double *tempXdata,int CENSOR,double *glmWT,double *WTdata);

int  DSA_PACK_update_design_matrix(int *tempmodel,int *lastXdesign,int nvarX,int maxsize,int rowtocompare);

/** temporary declaration. **/
void DSA_PACK_eval_model(int maxsize,int msize,int *tempmodel,int newterm_index,
			 int *changeindexset,int nvarX, int samplesize,double *tempXdata,double *Xdata,
			 double *tempterm, int CENSOR,double *glmWT,double *WTdata, int IP,
			 double *tempYdata, double *Ydata,double *PARA,double *RES,
			 int *nglmcall,double *newrss,int *nextmodel,int inddel, int binind, double* binWTdata, 
			 int always_doit, int usetree,int *modelfit_work_i,double *modelfit_work_d,int *lastXdesign);



void DSA_PACK_checknewterm(int *orderint_in_newterm,int *sumofpow_in_newterm,int *is_newterm_unique,
			   int *newterm,int nvarX,int *tempmodel,int *ithterm,int maxsize,int newterm_index,
			   int addmove_notsubmove);

int DSA_PACK_istermunique(int *newterm,int *tempmodel,int *ithterm,int maxsize,int nvarX,int newterm_index,
			  int addmove_notsubmove,int orderint_in_newterm, int sumofpow_in_newterm);

int DSA_PACK_exportresult(double *oldrss,double newrss,int msize,int nvarX,int *currentmodel,int *nextmodel,
			  int maxsize,double *bestarss,int *gCOUNT,int IP,int nglmcall,int **bestmodels,
			  int delind,double *currentarss);

void DSA_PACK_evaluate_loss(int lossFunction, double* X, double* Y, int N, int P, int M, 
			    int tdx, double* glmWT, double* binWTdata, double* estimates, 
			    double* loss);

unsigned long int DSA_PACK_ulongpow(int x, int p);
void DSA_PACK_ulongsort(unsigned long int *x, int n);

#endif 
