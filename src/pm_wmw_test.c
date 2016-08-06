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
// Z is just a shorter name for Yprime 
// X and Y form Up
// X1 and Z form Umw, it may or may not be the same as X
double compute_pm_wmw_Z(double * X, double * Y, double * xy, double * X1, double * Z, double * x1z, int m, int n, int corr, int method, double adj) {

    int i,j;
    int N=m+n;
    double alpha=1.0*n/N;
    double Up, Umw, tmp, mean_x, mean_y, z, var, ratio;
    double var_U_p_0, var_U_mw_0, cov_U_mw_p_0, cov_G_F;
    double *r=malloc((m+m)*sizeof(double));
    double *r2=malloc(N*sizeof(double));
    double *Fy_X = malloc(m*sizeof(double));
    double *Fx_Y = malloc(n*sizeof(double));
    
    for (i = 0; i < m; i++) xy[i]=X[i];
    for (j = 0; j < m; j++) xy[m+j]=Y[j];    
    for (i = 0; i < m; i++) x1z[i]=X1[i];
    for (j = 0; j < n; j++) x1z[m+j]=Z[j];    

    getrank (m+m, xy, r);
    getrank (m+n, x1z, r2);
    for (i = 0; i < m+m; i++) r[i]+=1; // 0 based in C, 1 based in R
    for (i = 0; i < N; i++)   r2[i]+=1; // 0 based in C, 1 based in R
    
    Up=0;
    for (j = 0; j < m; j++) Up+=r[m+j];    
    Up=(Up - m*(m+1.0)/2)/m/m;
    if(corr!=0) {
        if(corr==2) corr=Up>0.5?1:-1;        
        Up=Up-corr*0.5/m/m;
    }
    //PRINTF("%f ", Up); PRINTF("\n");
    
    Umw=0;
    for (j = 0; j < n; j++) Umw+=r2[m+j];    
    Umw=(Umw - n*(n+1.0)/2)/m/n;
    if(corr!=0) {
        if(corr==2) corr=Umw>0.5?1:-1;        
        Umw=Umw-corr*0.5/m/n;
    }
    //PRINTF("%f ", Umw); PRINTF("\n");
    
    for (i = 0; i < m; i++) {
        Fy_X[i]=0;
        for (j = 0; j < n; j++) Fy_X[i]+= Y[j]<X[i]?1:(Y[j]==X[i]?0.5:0);
        Fy_X[i]/=n;
    }
    for (j = 0; j < n; j++) {
        Fx_Y[j]=0;
        for (i = 0; i < m; i++) Fx_Y[j]+= X[i]<Y[j]?1:(X[i]==Y[j]?0.5:0);
        Fx_Y[j]/=m;
    }
    
    mean_x=0; for (i = 0; i < m; i++) mean_x+=Fy_X[i]; mean_x/=m;
    mean_y=1-mean_x; 
    tmp=0; for (i = 0; i < m; i++) tmp+=Fy_X[i] * Fx_Y[i]; tmp/=m;
    cov_G_F=(tmp - mean_x * mean_y)*m/(m-1);
    
    var_U_p_0=1./6-2*cov_G_F;
    var_U_mw_0=1./(12*alpha);
    cov_U_mw_p_0=1./12-cov_G_F;

    ratio=var_U_p_0/var_U_mw_0;
    z=sqrt(m*1.)*(Up-0.5 + ratio*(Umw-0.5));
    var=var_U_p_0 + 2*cov_U_mw_p_0*ratio + var_U_mw_0*ratio*ratio;
    //PRINTF("%f %f %f %f %f \n", Up, Umw, ratio, var, z);
    z=z/sqrt(var);
  

    free(Fy_X);
    free(Fx_Y);
    free(r);
    free(r2);
    
    return z;
}

// remember only use sexp to return value
SEXP pm_wmw_test(SEXP _X, SEXP _Y, SEXP _Z, SEXP _corr, SEXP _method, SEXP _mc_rep){
    
    int m=length(_X);
    int n=length(_Z);
    int N=m+n;
    int i,j,b;
    int corr=asInteger(_corr);// 0,1,-1,2
    int method=asInteger(_method);// 0,1,2,3,4
    double *X0=REAL(_X), *Y0=REAL(_Y), *Z0=REAL(_Z);  
    double ind;
      
    double *X = malloc(m*sizeof(double));
    double *Y = malloc(m*sizeof(double));
    double *X1 = malloc(m*sizeof(double));
    double *Z = malloc(n*sizeof(double));
    double *xy = malloc((m+m)*sizeof(double));
    double *x1z = malloc(N*sizeof(double));
    double *xz0 = malloc(N*sizeof(double));
    double *unique = malloc(N*sizeof(double));
    int *nties = malloc(N*sizeof(int));
        
    int mc_rep=asInteger(_mc_rep); // 0: exact perm, 1: z only, 1e4: mc
    SEXP _ans=PROTECT(allocVector(REALSXP, mc_rep));
    double *ans=REAL(_ans);
    
	int* index = (int*)calloc(m, sizeof(int));	
	int* perm = (int*)calloc(N, sizeof(int));	
        
    // save a copy in xy0
    for (i = 0; i < m; i++) xz0[i]=X0[i];
    for (j = 0; j < n; j++) xz0[m+j]=Z0[j];    
    
//    // NTIES <- table(r)
//    for (i = 0; i < N; i++) nties[i]=1;
//    int n_unique=0;
//    int flag;
//    for (i = 0; i < N; i++) {
//        flag=0;
//        for (j= 0; j < n_unique; j++) {
//            if (xy0[i]==unique[j]) {
//            //PRINTF("inside\n");
//                nties[j]++;
//                flag=1;
//                break;
//            }
//        }    
//        if(flag==0) unique[n_unique++]=xy0[i];
//    }
//    //for (i = 0; i < n_unique; i++) PRINTF("%f ", unique[i]); PRINTF("\n");
//    //for (i = 0; i < n_unique; i++) PRINTF("%i ", nties[i]); PRINTF("\n");
//    // sum(NTIES^3-NTIES)/(12*m*n*N*(N - 1))
    double adj=0;
//    for (i = 0; i < n_unique; i++) {
//        if(nties[i]>1) adj+= (pow(nties[i],3)-nties[i]);
//    }
//    adj/=(12.*m*n*N*(N-1));


    if(mc_rep==1) {
        ans[0]=compute_pm_wmw_Z(X0, Y0, xy, X0, Z0, x1z, m, n, corr, method, adj); 
    } else {
        if (mc_rep>1) {
            // Monte Carlo
            for (b=0; b<mc_rep; b++) {                
                // resample X and Y
                for (i = 0; i < m; i++) {
                    ind=RUNIF(0.0,1.0);
                    X[i]=ind< 0.5?X0[i]:Y0[i];
                    Y[i]=ind>=0.5?X0[i]:Y0[i];
                }
                // resample X and Z
                SampleNoReplace(m, N, index, perm);//index is the new X1, perm[1:n] is the new Z
                //for (i = 0; i < m; i++) PRINTF("%i ", index[i]); PRINTF("\n");
                for (i = 0; i < m; i++) X1[i]=xz0[index[i]-1];
                for (j = 0; j < n; j++) Z[j] =xz0[perm[j]-1];
                ans[b]=compute_pm_wmw_Z(X, Y, xy, X1, Z, x1z, m, n, corr, method, adj);        
            }    
        } else {
            // exact   
            // there is no exact for pm_wmw
        }
        //for (b=0; b<mc_rep; b++) PRINTF("%f ", ans[b]); PRINTF("\n");
    }
            
    free(X); free(Y); free(xy); free(x1z);
    free(unique); free(nties);

    UNPROTECT(1);
    return _ans;
}



