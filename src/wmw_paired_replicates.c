#include <float.h> //DBL_EPSILON
#include <R_ext/Lapack.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#define RUNIF runif
//#include <R_ext/Applic.h>
//#include <R_ext/BLAS.h>
//#include <R_ext/RS.h> //R definitions for 'extending' R, registering functions,...

#define PRINTF Rprintf
#define MAX(A,B)    ((A) > (B) ? (A) : (B))
#define MIN(A,B)    ((A) < (B) ? (A) : (B))

/* ********************************************************************** */


/* ---------------------------------------------------------------------- */



/* ********************************************************************** */

static const double* sortdata = NULL; /* used in the quicksort algorithm */

/* ---------------------------------------------------------------------- */

static
int compare(const void* a, const void* b)
/* Helper function for sort. Previously, this was a nested function under
 * sort, which is not allowed under ANSI C.
 */
{ const int i1 = *(const int*)a;
  const int i2 = *(const int*)b;
  const double term1 = sortdata[i1];
  const double term2 = sortdata[i2];
  if (term1 < term2) return -1;
  if (term1 > term2) return +1;
  return 0;
}


static void sort(int n, const double data[], int index[])
/* Sets up an index table given the data, such that data[index[]] is in
 * increasing order. Sorting is done on the indices; the array data
 * is unchanged.
 */
{ int i;
  sortdata = data;
  for (i = 0; i < n; i++) index[i] = i;
  qsort(index, n, sizeof(int), compare);
}

static void getrank (int n, double data[], double rank[])
/* Calculates the ranks of the elements in the array data. Two elements with
 * the same value get the same rank, equal to the average of the ranks had the
 * elements different values. The ranks are returned as a newly allocated
 * array that should be freed by the calling routine. If getrank fails due to
 * a memory allocation error, it returns NULL.
 */
{ int i;
  int* index;
  index = malloc(n*sizeof(int));
  if (!index)
  { free(rank);
    return;
  }
  /* Call sort to get an index table */
  sort (n, data, index);
  /* Build a rank table */
  for (i = 0; i < n; i++) rank[index[i]] = i;
  /* Fix for equal ranks */
  i = 0;
  while (i < n)
  { int m;
    double value = data[index[i]];
    int j = i + 1;
    while (j < n && data[index[j]] == value) j++;
    m = j - i; /* number of equal ranks found */
    value = rank[index[i]] + (m-1)/2.;
    for (j = i; j < i + m; j++) rank[index[j]] = value;
    i += m;
  }
  free (index);
}

static void SampleNoReplace(int k, int n, int *y, int *x)
{
    int i, j;
    int n0=n;
#ifndef SCYTHE_COMPILE_DIRECT    
    GetRNGstate();    
#endif
    for (i = 0; i < n; i++) x[i] = i;
    for (i = 0; i < k; i++) {
    	j = (int)((double)n * RUNIF(0.0,1.0));
 	    y[i] = x[j] + 1;
 	    x[j] = x[--n];
    }
    // make the first n-k element of x complement y
    for (i = 0; i < n0-k; i++) x[i] = x[i]+1;
#ifndef SCYTHE_COMPILE_DIRECT    
    PutRNGstate();    
#endif
}




// xy is passed in so that there is no need to allocate memory for getting rank
// m=n
double compute_wmw_paired_replicates_stat(double *X0, double *Y0, double *X, double *Y, int * lenY, int m) {

    int i,j,k;
        
	int cnt=0;
    for (i = 0; i < m; i++) {
        k=lenY[i];
    	int* index = (int*)malloc(1*sizeof(int));	
    	int* perm =  (int*)malloc((1+k)*sizeof(int));	
        SampleNoReplace(1, 1+k, index, perm);//index is the new X1, perm[1:n] is the new Z  
        //PRINTF("%i %i", k, index[0]); PRINTF("\n");
        
        for (j=1; j<index[0]; j++) {Y[cnt]=Y0[cnt]; cnt++;}
        if (index[0]==k+1) X[i]=X0[i]; else {
            X[i]=Y0[cnt];
            Y[cnt]=X0[i];
            cnt++;
            for (j=index[0]+1; j<=k; j++) {Y[cnt]=Y0[cnt]; cnt++;}
        }
        
        free(index); free(perm);                
    }
    //PRINTF("%f %f", stat, mu); PRINTF("\n");

    return 0;
}

// remember only use sexp to return value
SEXP wmw_paired_replicates(SEXP _X, SEXP _Y, SEXP _lenY, SEXP _corr, SEXP _method, SEXP _mc_rep){
    
    int m=length(_X);
    double *X0=REAL(_X), *Y0=REAL(_Y);
    int *lenY=INTEGER(_lenY);  
    double Umw;
    int i,b;
//    int corr=asInteger(_corr);// 0,1,-1,2. correction
//    int method=asInteger(_method);// 0,1,2,3,4

	int n=0;
    for (i = 0; i < m; i++) n+=lenY[i];
            
    double *X = malloc(m*sizeof(double));
    double *Y = malloc(n*sizeof(double));
    double *xy = malloc((m+n)*sizeof(double));
    double *r = malloc((m+n)*sizeof(double));
    
    int mc_rep=asInteger(_mc_rep); // 0: exact perm, 1: z only, 1e4: mc
    SEXP _ans=PROTECT(allocVector(REALSXP, mc_rep));
    double *ans=REAL(_ans);    

    // Monte Carlo
    for (b=0; b<mc_rep; b++) {                
        ans[b]=compute_wmw_paired_replicates_stat(X0, Y0, X, Y, lenY, m);
        //for (i=0; i<n; i++) PRINTF("%f ", Y[i]); PRINTF("\n");
        //for (i=0; i<m; i++) PRINTF("%f ", X[i]); PRINTF("\n");
        for (i=0; i<m; i++) xy[i]=X[i];
        for (i=0; i<n; i++) xy[i+m]=Y[i];
        getrank (m+n, xy, r);
        for (i = 0; i < m+n; i++) r[i]+=1; // 0 based in C, 1 based in R
        Umw=0;
        for (i = 0; i < m; i++) Umw+=r[i];
        Umw=Umw - m*(m+1)/2.0;
        ans[b]=Umw;
    }   
    //for (b=0; b<mc_rep; b++) PRINTF("%f ", ans[b]); PRINTF("\n");
            
    free(X); free(Y); free(xy); free(r);

    UNPROTECT(1);
    return _ans;
}



