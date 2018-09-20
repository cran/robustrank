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
// X and Y form Up. Xprime and Yprime form Umw. Xprime is allowed to be the same as X
double _mw_mw_2_perm (double * X, double * Y, double * xy, double * Xprime, double * Yprime, double * xyprime, int m, int n, int l, int corr, double adj_ties) {

    int i,j;
    double alpha=1.0*n/(m+n);
    double beta=1.0*l/(m+l);
    double Up, Umw, tmp, mean_x, mean_y, z, var, ratio;
    double var_U_p_0, var_U_mw_0;
    double *Fy_X = malloc(m*sizeof(double));
    double *Fx_Y = malloc(m*sizeof(double));
    
    double corr1;
    
    
    for (i = 0; i < m; i++) xy[i]=X[i];
    for (j = 0; j < m; j++) xy[m+j]=Y[j];    
    double *r=malloc((m+m)*sizeof(double));
    getrank (m+m, xy, r);
    for (i = 0; i < m+m; i++) r[i]+=1; // 0 based in C, 1 based in R
    Up=0;
    for (j = 0; j < m; j++) Up+=r[m+j];    
    Up=(Up - m*(m+1.0)/2)/m/m;
    //PRINTF("Up %f corr %d", Up, corr); PRINTF("\n");
    // corr is 0, 1, -1 or 2
    if(corr!=0) {
        if(corr==2) {
            if(Up==0.5) corr1=0; else corr1=Up>0.5?1:-1; 
        } else corr1=corr;
        Up=Up-corr1*0.5/m/m;
    }
    //PRINTF("Up %f corr %d", Up, corr); PRINTF("\n");
    
    //PRINTF("l %d n %d", l, n); PRINTF("\n");
    for (i = 0; i < l; i++) xyprime[i]=Xprime[i];
    for (j = 0; j < n; j++) xyprime[l+j]=Yprime[j];    
    double *r2=malloc((n+l)*sizeof(double));
    getrank (l+n, xyprime, r2);
    for (i = 0; i < l+n; i++)   r2[i]+=1; // 0 based in C, 1 based in R
    Umw=0;
    for (j = 0; j < n; j++) Umw+=r2[l+j];    
    //for (i = 0; i < l+n; i++) PRINTF("%f ", xyprime[i]); PRINTF("\n");
    //for (i = 0; i < l+n; i++) PRINTF("%f ", r2[i]); PRINTF("\n");
    //PRINTF("Umw before centering scaling %f ", Umw); PRINTF("\n");
    Umw=(Umw - n*(n+1.0)/2)/l/n;
    //PRINTF("Umw %f corr %d", Umw, corr); PRINTF("\n");
    if(corr!=0) {
        if(corr==2) {
            if(Umw==0.5) corr1=0; else corr1=Umw>0.5?1:-1;        
        } else corr1=corr;
        Umw=Umw-corr1*0.5/l/n;
    }
    //PRINTF("%f ", Umw); PRINTF("\n");
    
    for (i = 0; i < m; i++) {
        Fy_X[i]=0;
        for (j = 0; j < m; j++) Fy_X[i]+= Y[j]<X[i]?1:(Y[j]==X[i]?0.5:0);
        Fy_X[i]/=m;
    }
    for (j = 0; j < m; j++) {
        Fx_Y[j]=0;
        for (i = 0; i < m; i++) Fx_Y[j]+= X[i]<Y[j]?1:(X[i]==Y[j]?0.5:0);
        Fx_Y[j]/=m;
    }
    
    mean_x=0; for (i = 0; i < m; i++) mean_x+=Fy_X[i]; mean_x/=m;
    mean_y=1-mean_x; 
    tmp=0; for (i = 0; i < m; i++) tmp+=Fy_X[i] * Fx_Y[i]; tmp/=m;
    double cov_G_F=(tmp - mean_x * mean_y)*m/(m-1); // there are better estimates of this in pair_wmw_test.c, but for perm test, this should be okay
    
    var_U_p_0=1./6-2*cov_G_F;
    var_U_mw_0=(1/alpha+1/beta-2)/12.;

    ratio=var_U_p_0/var_U_mw_0;
    z=sqrt(m*1.)*(Up-0.5 + ratio*(Umw-0.5));
    var=var_U_p_0 + var_U_mw_0*ratio*ratio;
    //PRINTF("Up %f Umw %f var_U_p_0 %f var_U_mw_0 %f z %f \n", Up, Umw, var_U_p_0, var_U_mw_0, z);
    z=z/sqrt(var);
  

    free(Fy_X);
    free(Fx_Y);
    free(r);
    free(r2);
    
    return z;
}


SEXP mw_mw_2_perm (SEXP _X, SEXP _Y, SEXP _Xprime, SEXP _Yprime, 
     SEXP _corr, 
     SEXP _mc_rep, SEXP _xy_idx, SEXP _xyprime_idx){
    
    double *X0=REAL(_X), *Y0=REAL(_Y), *Xprime0=REAL(_Xprime), *Yprime0=REAL(_Yprime);  
    int corr=asInteger(_corr);// continuity correction: 0,1,-1,2
    int mc_rep=asInteger(_mc_rep); // 0: exact perm, 1: z only, 1e4: mc
    int *xy_idx=INTEGER(_xy_idx);
    int *xyprime_idx=INTEGER(_xyprime_idx);

    int m=length(_X);
    int n=length(_Yprime);
    int l=length(_Xprime);
    int ln=l+n;

    int i,j;
    double ind;

    // save a copy of data in xyprime0 for resampling
    double *xyprime0 = malloc(ln*sizeof(double));
    for (i = 0; i < l; i++) xyprime0[i]  =Xprime0[i];
    for (j = 0; j < n; j++) xyprime0[l+j]=Yprime0[j];        
        
    // for resampling use
    double *X = malloc(m*sizeof(double));
    double *Y = malloc(m*sizeof(double));
    double *xy = malloc((m+m)*sizeof(double));    
    double *Xprime = malloc(l*sizeof(double));
    double *Yprime = malloc(n*sizeof(double));
    double *xyprime = malloc(ln*sizeof(double));
    
    int nperm_1=length(_xy_idx)/m;
    int nperm_2=length(_xyprime_idx)/ln; // number of permutation
    // mc_rep 0 means we are only doing           original stat
    // mc_rep -1 means we are doing permutation + original
    // mc_rep other numbers are MC replicates   + original
    SEXP _ans=PROTECT(allocVector(REALSXP, 1+(mc_rep==-1?(nperm_1*nperm_2):mc_rep) )); // +1 b/c always compute stat for the original data
    double *ans=REAL(_ans);
    
//    double *unique = malloc(N*sizeof(double));
//    int *nties = malloc(N*sizeof(int));      
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
    double adj_ties=0;
//    for (i = 0; i < n_unique; i++) {
//        if(nties[i]>1) adj_ties+= (pow(nties[i],3)-nties[i]);
//    }
//    adj_ties/=(12.*m*n*N*(N-1));


    ans[0]=_mw_mw_2_perm(X0, Y0, xy, Xprime0, Yprime0, xyprime, m, n, l, corr, adj_ties); 
    
    if(mc_rep>0) {
        // Monte Carlo
    	int* index = (int*)calloc(l, sizeof(int));	
    	int* perm = (int*)calloc(ln, sizeof(int));	
            
        for (int b=0; b<mc_rep; b++) {                
            // resample X and Y
            for (i = 0; i < m; i++) {
                ind=RUNIF(0.0,1.0);
                X[i]=ind< 0.5?X0[i]:Y0[i];
                Y[i]=ind>=0.5?X0[i]:Y0[i];
            }
            // resample Xprime and Yprime
            SampleNoReplace(l, ln, index, perm);//index is the new Xprime, perm[1:n] is the new Yprime
            //for (i = 0; i < m; i++) PRINTF("%i ", index[i]); PRINTF("\n");
            for (i = 0; i < l; i++) Xprime[i]=xyprime0[index[i]-1];
            for (j = 0; j < n; j++) Yprime[j] =xyprime0[perm[j]-1];
            ans[b+1]=_mw_mw_2_perm(X, Y, xy, Xprime, Yprime, xyprime, m, n, l, corr, adj_ties); // +1 b/c first one is the original stat
        }    
        
        free(index); free(perm); 
    
    } else if(mc_rep<0) {
        // exact, aka permutation
     	int c1=0; // index of xy_idx 
        for (int b1=0; b1<nperm_1; b1++) {
            
            // resample paired X and Y
            for (i = 0; i < m; i++) {
                ind=xy_idx[c1++];   
                X[i]=ind==0?X0[i]:Y0[i];
                Y[i]=ind==1?X0[i]:Y0[i];
            }            
            
            //PRINTF("X : "); for (i = 0; i < m; i++) PRINTF("%f ", X[i]); PRINTF("\n");
            //PRINTF("Y : "); for (i = 0; i < m; i++) PRINTF("%f ", Y[i]); PRINTF("\n");
                
            int c2=0; // index of xyprime_idx
            int c_xprime, c_yprime; // indexes
            for (int b2=0; b2<nperm_2; b2++) {
                
                // resample unpaired Xprime and Yprime
                c_xprime=0; c_yprime=0; // reset indexes
                for (i = 0; i < l+n; i++) {
                    if (xyprime_idx[c2++]==1) Xprime[c_xprime++]=xyprime0[i]; else Yprime[c_yprime++]=xyprime0[i]; 
                }            
                
                //PRINTF("X': "); for (i = 0; i < l; i++) PRINTF("%f ", Xprime[i]); PRINTF("\n");
                //PRINTF("Y': "); for (i = 0; i < n; i++) PRINTF("%f ", Yprime[i]); PRINTF("\n");
                ans[1+nperm_2*b1+b2]=_mw_mw_2_perm(X, Y, xy, Xprime, Yprime, xyprime, m, n, l, corr, adj_ties); // +1 b/c first one is the original stat            
            }    
        }    

        //for (b=0; b<mc_rep; b++) PRINTF("%f ", ans[b]); PRINTF("\n");
    }
            
    free(X); free(Y); free(xy); free(xyprime); free(Xprime); free(Yprime); free(xyprime0); 
    //free(unique); free(nties); 

    UNPROTECT(1);
    return _ans;
}

