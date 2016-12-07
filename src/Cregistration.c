#include <Rdefines.h>
#include <R_ext/Rdynload.h>

void DSA_PACK_Rentry(int *gCONST,double *Xlearn,double *Ylearn,double *CVrisks,double *PARA,
		     int *vsplits,double *WTtrainvallearn,int *nforced,int *forcedterms,
		     int *maxntrain, int *maxnval,int *bestmodelsspace,double *binWTtrainvallearn);

void DSA_PACK_setMessagelevel(int *newMlevel);

void DSA_PACK_getPrimes(int *prime);

void DSA_PACK_DeleteTree();

static R_NativePrimitiveArgType DSA_PACK_Rentry_t[13]= {INTSXP,REALSXP,REALSXP,REALSXP,
							REALSXP,INTSXP,REALSXP,INTSXP,INTSXP,
							INTSXP,INTSXP,INTSXP,REALSXP};

static R_NativePrimitiveArgType DSA_PACK_setMessagelevel_t[1]={INTSXP};

static R_NativePrimitiveArgType DSA_PACK_getPrimes_t[1]={INTSXP};

static R_NativePrimitiveArgType DSA_PACK_DeleteTree_t[0]={};

static const  R_CMethodDef cMethods[] = {
      {"DSA_PACK_Rentry", (DL_FUNC) &DSA_PACK_Rentry, 13, DSA_PACK_Rentry_t},
      {"DSA_PACK_setMessagelevel", (DL_FUNC) &DSA_PACK_setMessagelevel, 1, DSA_PACK_setMessagelevel_t},
      {"DSA_PACK_getPrimes", (DL_FUNC) &DSA_PACK_getPrimes, 1, DSA_PACK_getPrimes_t},
      {"DSA_PACK_DeleteTree", (DL_FUNC) &DSA_PACK_DeleteTree, 0, DSA_PACK_DeleteTree_t},
      {NULL, NULL, 0}
    };

 void R_init_DSA(DllInfo *info)
 {
   R_registerRoutines(info,cMethods,NULL,NULL,NULL);
 }


