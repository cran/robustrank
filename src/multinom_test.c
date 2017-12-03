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
double compute_multinom_stat(int * lenY, int m) {

    int i,k;
    double stat,mu;
        
	stat=0; mu=0;
    for (i = 0; i < m; i++) {
        k=lenY[i];
    	int* index = (int*)malloc(1*sizeof(int));	
    	int* perm =  (int*)malloc((1+k)*sizeof(int));	
        SampleNoReplace(1, 1+k, index, perm);//index is the new X1, perm[1:n] is the new Z  
        //PRINTF("%i %i", k, index[0]); PRINTF("\n");
        stat+=index[0]-1;      
        mu+=1*(k+1-1)/2.0; // it has to be 2.0!
        free(index); free(perm);                
    }
    //PRINTF("%f %f", stat, mu); PRINTF("\n");

    return stat-mu;
}

// remember only use sexp to return value
SEXP multinom_test(SEXP _X, SEXP _Y, SEXP _lenY, SEXP _corr, SEXP _method, SEXP _mc_rep){
    
    int m=length(_X);
//    double *X0=REAL(_X), *Y0=REAL(_Y);
    int *lenY=INTEGER(_lenY);  
//    int corr=asInteger(_corr);// 0,1,-1,2. correction
//    int method=asInteger(_method);// 0,1,2,3,4
            
    int mc_rep=asInteger(_mc_rep); // 0: exact perm, 1: z only, 1e4: mc
    SEXP _ans=PROTECT(allocVector(REALSXP, mc_rep));
    double *ans=REAL(_ans);    

    // Monte Carlo
    int b;
    for (b=0; b<mc_rep; b++) {                
        ans[b]=compute_multinom_stat(lenY, m);
    }   
    //for (b=0; b<mc_rep; b++) PRINTF("%f ", ans[b]); PRINTF("\n");
            

    UNPROTECT(1);
    return _ans;
}



