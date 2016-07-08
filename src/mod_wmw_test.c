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


void sort(int n, const double data[], int index[])
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


void SampleNoReplace(int k, int n, int *y, int *x)
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

double compute_z(double * X, double * Y, double * xy, int m, int n, int corr, int method, double adj) {

    int i,j;
    int N=m+n;
    double U, a, b, mean_x, mean_y, z; // a1, b1;
    double var_wmw, var_exact, var_large, var_3, var, var_nsm3;
    double *r=malloc(N*sizeof(double));
    double *Fy_X = malloc(m*sizeof(double));
    double *Fx_Y = malloc(n*sizeof(double));
    
    for (i = 0; i < m; i++) xy[i]=X[i];
    for (j = 0; j < n; j++) xy[m+j]=Y[j];    

    getrank (N, xy, r);
    for (i = 0; i < N; i++) r[i]+=1; // 0 based in C, 1 based in R
    
    U=0;
    for (j = 0; j < n; j++) U+=r[m+j];    
    U=(U - n*(n+1.0)/2)/m/n;
    if(corr!=0) {
        if(corr==2) corr=U>0.5?1:-1;        
        U=U-corr*0.5/m/n;
    }
    //PRINTF("%f ", U); PRINTF("\n");
    
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
    
    a=0; for (i = 0; i < m; i++) a+=pow(Fy_X[i],2); a/=m;
    mean_x=0; for (i = 0; i < m; i++) mean_x+=Fy_X[i]; mean_x/=m;
    a=(a-pow(mean_x,2))*m/(m-1);

    b=0; for (j = 0; j < n; j++) b+=pow(Fx_Y[j],2); b/=n;
    mean_y=1-mean_x;
    b=(b-pow(mean_y,2))*n/(n-1);
    
    var_wmw=1.0/12/m + 1.0/12/n + 1.0/12/m/n - adj;
    var_exact=a/m*(1-1.0/n) + b/n*(1-1.0/m) + mean_x/m * mean_y/n;
    var_nsm3= a/m*(1-1.0/m) + b/n*(1-1.0/n) + mean_x/m * mean_y/n;
    var_large=a/m + b/n;
    var_3=var_wmw<var_exact?var_wmw:var_exact;
    
    // attempt to address numerical instability
    //a1=0; for (i = 0; i < m; i++) a1+=pow(Fy_X[i],2);     
    //b1=0; for (j = 0; j < n; j++) b1+=pow(Fx_Y[j],2); 
    //var_nsm3=(a1*n*n + b1*m*m)/m/m/n/n - (m-(2*m+1)*mean_x+(N+1)*mean_x*mean_x)/m/n; 

    if(method==0) var=var_wmw; else if (method==1) var=var_exact; else if (method==2) var=var_large; else if (method==3) var=var_3; else if (method==4) var=var_nsm3; else var=0;//last one cannot happen
    z=(U-0.5)/sqrt(var);
    //PRINTF("%f %f %f %f %f \n", U, var_wmw, var_large, var_exact, z);
    //PRINTF("%f %f %f %f %f %f\n", a1, b1, a1/m/m + b1/n/n, - (m-(2*m+1)*mean_x+(N+1)*mean_x*mean_x)/m/n, mean_x, z);
    
    free(Fy_X);
    free(Fx_Y);
    free(r);
    
    return z;
}


SEXP mod_wmw_test(SEXP _X, SEXP _Y, SEXP _corr, SEXP _method, SEXP _mc_rep, SEXP _comb){
    
    int m=length(_X);
    int n=length(_Y);
    int N=m+n;
    int corr=asInteger(_corr);// 0,1,-1,2
    int method=asInteger(_method);// 0,1,2,3
    double *X0=REAL(_X), *Y0=REAL(_Y);
    int *comb=INTEGER(_comb);
      
    double *xy0 = malloc(N*sizeof(double));
    double *X = malloc(m*sizeof(double));
    double *Y = malloc(n*sizeof(double));
    double *xy = malloc(N*sizeof(double));
    double *unique = malloc(N*sizeof(double));
    int *nties = malloc(N*sizeof(int));

    int mc_rep=asInteger(_mc_rep); // 0: exact perm, 1: z only, 1e4: mc
    int nperm=length(_comb)/m; // number of permutation
    SEXP _ans=PROTECT(allocVector(REALSXP, mc_rep==0?nperm:mc_rep));
    double *ans=REAL(_ans);
    
    int i,j,b;
	int* index = (int*)calloc(m, sizeof(int));	
	int* perm = (int*)calloc(N, sizeof(int));	
    int c,index_cur,perm_cur;
        
    // save a copy in xy0
    for (i = 0; i < m; i++) xy0[i]=X0[i];
    for (j = 0; j < n; j++) xy0[m+j]=Y0[j];    
    
    // NTIES <- table(r)
    for (i = 0; i < N; i++) nties[i]=1;
    int n_unique=0;
    int flag;
    for (i = 0; i < N; i++) {
        flag=0;
        for (j= 0; j < n_unique; j++) {
            //PRINTF("%f %f", xy0[i], unique[j]); PRINTF("\n");
            if (xy0[i]==unique[j]) {
            //PRINTF("inside\n");
                nties[j]++;
                flag=1;
                break;
            }
        }    
        if(flag==0) unique[n_unique++]=xy0[i];
    }
    //for (i = 0; i < n_unique; i++) PRINTF("%f ", unique[i]); PRINTF("\n");
    //for (i = 0; i < n_unique; i++) PRINTF("%i ", nties[i]); PRINTF("\n");
    // sum(NTIES^3-NTIES)/(12*m*n*N*(N - 1))
    double adj=0;
    for (i = 0; i < n_unique; i++) {
        if(nties[i]>1) adj+= (pow(nties[i],3)-nties[i]);
    }
    adj/=(12.*m*n*N*(N-1));
    

    if(mc_rep==1) {
        ans[0]=compute_z(X0, Y0, xy, m, n, corr, method, adj); 
    } else {
        if (length(_comb)==1) {
            // Monte Carlo
            for (b=0; b<mc_rep; b++) {
                SampleNoReplace(m, N, index, perm);//index is the new X, perm[1:n] is the new Y
                //for (i = 0; i < m; i++) PRINTF("%i ", index[i]); PRINTF("\n");
                for (i = 0; i < m; i++) X[i]=xy0[index[i]-1];
                for (j = 0; j < n; j++) Y[j]=xy0[perm[j]-1];
                ans[b]=compute_z(X, Y, xy, m, n, corr, method, adj);        
            }    
        } else {
            // exact   
         	c=0;
            for (b=0; b<nperm; b++) {
                //PRINTF("%i ", b); 
                for (i = 0; i < m; i++) index[i]=comb[c++];        
                index_cur=0; perm_cur=0;
                for (i=1; i<=N; i++) {
                    if (i==index[index_cur]) {index_cur++; continue;} else perm[perm_cur++]=i;
                }  
                //for (i = 0; i < m; i++) PRINTF("%i ", index[i]); PRINTF("\n");
                //for (j = 0; j < n; j++) PRINTF("%i ", perm [j]); PRINTF("\n");
                
                for (i = 0; i < m; i++) X[i]=xy0[index[i]-1];
                for (j = 0; j < n; j++) Y[j]=xy0[perm[j]-1];
                ans[b]=compute_z(X, Y, xy, m, n, corr, method, adj);        
            }    
        }
    }
                
    free(X); free(Y); free(xy); free(xy0); 
    free(index); free(perm); 
    free(unique); free(nties);
    
    UNPROTECT(1);
    return _ans;
}



