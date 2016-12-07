#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <R.h>
#include <R_ext/Utils.h>
#include "dsa.h"




void DSA_PACK_printmatrix_to_file(double *x, int nrow, int ncol, char* f)
{
  FILE* fp = fopen(f, "w");

  int i,j;
  for (i=0;i<nrow;i++)
    {
      for(j=0;j<ncol;j++)
	fprintf(fp, "%lf ", x[j*nrow+i]); 
      fprintf(fp, "\n");
    }

  fclose(fp);
}


void DSA_PACK_printmatrix(double *x, int nrow, int ncol)
{
  int i,j;
  for (i=0;i<nrow;i++)
    {
      for(j=0;j<ncol;j++)Rprintf("%lf ",x[j*nrow+i]); 
      Rprintf("\n");
    }
}

void DSA_PACK_printmatrixl(int *x, int nrow, int ncol)
{
  int i,j;
  for (i=0;i<nrow;i++)
    {
      Rprintf("\n");
      for(j=0;j<ncol;j++)Rprintf("%ld ",x[j*nrow+i]); 
    }
}

double DSA_PACK_getmin(double *x,int n)
{
  double mymin;
  int i;

  mymin = *x;
  for(i=1;i<n;i++)if( (x[i]<mymin)  & (x[i]!=R_PosInf))mymin=x[i];
  return(mymin);
}

int DSA_PACK_rowsumints(int *x,int nrow,int ncol,int whichrow)
{
  int j,m;

  m=0;
  for(j=0;j<ncol;j++)m += x[j*nrow+whichrow];
  return(m);
}

int DSA_PACK_num_nonzero(int *x,int n)
{
  int i,count;

  count=0;
  for(i=0;i<n;i++)if(x[i]!=0)count++;
  return(count);
}

int DSA_PACK_sumofints(int *x,int n)
{
  int i,m;

  for(i=0,m=0;i<n;i++,m+=*x,x++);
  return(m);
}

/*This is not the most efficient way of computing the power but good enough for now*/
R_INLINE unsigned long int DSA_PACK_ulongpow(int x, int p)
{
  int i;
  unsigned long int res=1;

  for(i=0;i<p;i++)res*=(unsigned long int)(x);
  return(res);
}

/*Lifted from sort.c in R source code and adapted from R_isort - will not work with NAs*/
void DSA_PACK_ulongsort(unsigned long int *x, int n)
{
  unsigned long int v;

  int i, j, h;

  for (h = 1; h <= n / 9; h = 3 * h + 1);
  for (; h > 0; h /= 3)
    for (i = h; i < n; i++) {
      v = x[i];
      j = i;
      while (j >= h && ((x[j - h]<v)?(-1):(1)) > 0)
	{ x[j] = x[j - h]; j -= h; }
      x[j] = v;
    }

}
