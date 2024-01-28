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




// xy is passed in so that there is no need to allocate memory for getting rank
// m=n
double compute_pair_wmw_Z(double * X, double * Y, double * xy, int m, int n, int corr, int method, double adj, int return_var) {

    int i,j,k;
    int N=m+n;
    double U, a, b, c, mean_x, mean_y, z;
    double tmp_1=0, tmp_2=0, tmp_3=0, tmp_4=0, tmp_5=0, tmp_6=0, theta=0, tao=0; 
    double A, B, C_1, v_f, A_0, B_0, C_1_0, v_f_0;
    double C_2_cov, C_2_cov_0, C_2_cov_1, C_2_var, C_2_var_0;
    double var_exact, var_exact_0, var_exact_1, var_exact_2, var_exact_3, var_large, var_large_0, var;
    
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
    double corr1;
    if(corr!=0) {
        if(corr==2) {
            if(U==0.5) corr1=0; else corr1=U>0.5?1:-1;        
        } else corr1=corr;
        U=U-corr1*0.5/m/n;
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
    a=(a-pow(mean_x,2))*m/(m-1); // *m/(m-1) is finite sample adj.

    b=0; for (j = 0; j < n; j++) b+=pow(Fx_Y[j],2); b/=n;
    mean_y=1-mean_x; 
    b=(b-pow(mean_y,2))*n/(n-1);
    
    // - 2*cov(Fy(X), Fx(Y))
    c=0; for (i = 0; i < m; i++) c+=Fy_X[i] * Fx_Y[i]; c/=m;
    c=(c - mean_x * mean_y)*m/(m-1);
    var_large  =    a + b     - 2*c;
    var_large_0=1./12 + 1./12 - 2*c;
    
    //tmp.1=mean(outer(1:m,1:n, function(i,j) ifelse(i!=j,Y[i]>X[i] & Y[j]>X[i],NA)), na.rm=TRUE) 
    //tmp.2=mean(outer(1:m,1:n, function(i,j) ifelse(i!=j,Y[i]>X[i] & Y[i]>X[j],NA)), na.rm=TRUE)
    //tmp.3=mean(outer(1:m,1:n, function(i,j) ifelse(i!=j,Y[j]>X[i] & Y[i]>X[j],NA)), na.rm=TRUE)
    //theta=mean(outer(1:m,1:n, function(i,j) ifelse(i!=j,Y[j]>X[i],NA)), na.rm=TRUE) # note that theta!=FxY here b/c i==j not counted
    for(i = 0;i < m;i++){
    for(j = 0;j < m;j++){
          if(i==j) continue;
          tmp_1+=(Y[i]>X[i])*(Y[j]>X[i]);
          tmp_2+=(Y[i]>X[i])*(Y[i]>X[j]);
          tmp_3+=(Y[j]>X[i])*(Y[i]>X[j]);
          theta+=(Y[j]>X[i]);
    }}    
    tmp_1/=m*(m-1);
    tmp_2/=m*(m-1);
    tmp_3/=m*(m-1);
    theta/=m*(m-1);
    
    //tmp.4=outer.3.mean(1:m,1:m,1:m, function(i,j,k) ifelse(i==j|i==k|j==k,NA,(Y[j]>X[i])*(Y[i]>X[k])) )    
    //tmp.5=outer.3.mean(1:m,1:m,1:m, function(i,j,k) ifelse(i==j|i==k|j==k,NA,(Y[j]>X[i])*(Y[k]>X[i])) )
    //tmp.6=outer.3.mean(1:m,1:m,1:m, function(i,j,k) ifelse(i==j|i==k|j==k,NA,(Y[j]>X[i])*(Y[j]>X[k])) )
    for(i = 0;i < m;i++){
    for(j = 0;j < m;j++){
          if(i==j) continue;
          for(k = 0;k < m;k++){
              if (k==i || k==j) continue;
              tmp_4 +=(Y[j]>X[i])*(Y[i]>X[k]);
              tmp_5 +=(Y[j]>X[i])*(Y[k]>X[i]);
              tmp_6 +=(Y[j]>X[i])*(Y[j]>X[k]);
          }
    }}    
    tmp_4/=m*(m-1)*(m-2);
    tmp_5/=m*(m-1)*(m-2);
    tmp_6/=m*(m-1)*(m-2);
    
    //C.2.cov=    m*(m-1)*(m-2)* (tmp.4 - theta.sq) * 2 
    //C.2.var=    m*(m-1)*(m-2)* (tmp.5 - theta.sq + tmp.6 - theta.sq) *finite.adj 
    C_2_cov_0=m*(m-1)*(m-2)* (tmp_4 - 1.0/4) * 2 ;
    C_2_cov  =m*(m-1)*(m-2)* (tmp_4 - theta*theta) * 2 ;
    C_2_var=  m*(m-1)*(m-2)* (tmp_5 - theta*theta + tmp_6 - theta*theta);
    C_2_var_0=m*(m-1)*(m-2)* (1./12 + 1./12);
    //PRINTF("%f %f %f \n", tmp_4, theta, C_2_cov);
        
    //tao=mean(Y>X)
    for(i = 0;i < m;i++) tao+=Y[i]>X[i]; 
    tao/=m;
    A  =m*        tao*(1-tao);
    B  =m*(m-1)*  (tmp_1 - tao*theta + tmp_2 - tao*theta) ;
    C_1=m*(m-1)*  (theta*(1-theta)  + (tmp_3 - theta*theta));
    v_f=A + 2*B + C_1;

    //var.exact.0=(C.2.var.0 + C.2.cov + v.f(theta,tao))/m^3 
    //var.exact  =(C.2.var   + C.2.cov + v.f(theta,tao))/m^3 
    var_exact  =(C_2_var   + C_2_cov + v_f)/pow(m,3);    
    var_exact_0=(C_2_var_0 + C_2_cov_0+v_f)/pow(m,3);
    var_exact_1=(C_2_var_0 + C_2_cov + v_f)/pow(m,3);
    C_2_cov_1=m*(m-1)*(m-2)* (tmp_4 - (theta*theta-var_exact_1/m)) * 2; // a less biased estimate of theta^2, a one-step correction
    var_exact_2=(C_2_var_0 + C_2_cov_1+v_f)/pow(m,3);

    A_0  =m*        1./2*(1-1./2);
    B_0  =m*(m-1)*  (tmp_1 - tao*theta + tmp_2 - tao*theta) ;
    C_1_0=m*(m-1)*  (1./2*(1-1./2)  + (tmp_3 - (theta*theta-var_exact_1/m)));
    v_f_0=A_0 + 2*B_0 + C_1_0;    
    var_exact_3=(C_2_var_0 + C_2_cov_1+v_f_0)/pow(m,3);
    //PRINTF("%f %f %f %f \n", C_2_var_0, C_2_cov, v_f, C_2_var);

    //PRINTF("%f %f %f %f %f \n", U, var_large_0, var_large, var_exact_0, var_exact);
    switch(method){
        case -1: var=var_large_0; break;
        case -2: var=var_large; break;
        case -3: var=var_exact; break;
        case 0: var=var_exact_0; break;
        case 1: var=var_exact_1; break;
        case 2: var=var_exact_2; break;
        case 3: var=var_exact_3; break;
        default: var=0;
    }
    //if(method==0) var=var_large_0; else if (method==1) var=var_large; else if (method==2) var=var_exact_0; else if (method==3) var=var_exact; else if (method==4) var=var_exact_2; else if (method==5) var=var_exact_3; else var=0;//last one cannot happen
    if(return_var) z=var; else z=sqrt(m*1.)*(U-0.5)/sqrt(var); 
    //PRINTF("%f %f %f %f %f %f\n", a1, b1, a1/m/m + b1/n/n, - (m-(2*m+1)*mean_x+(N+1)*mean_x*mean_x)/m/n, mean_x, z);
  

    free(Fy_X);
    free(Fx_Y);
    free(r);
    
    return z;
}

// m and n should be the same. This code should be updated at some point.
SEXP pair_wmw_test(SEXP _X, SEXP _Y, SEXP _corr, SEXP _method, SEXP _mc_rep, SEXP _comb){
    
    int m=length(_X);
    int n=length(_Y);
    int N=m+n;
    int i,b;
    int corr=asInteger(_corr);// 0,1,-1,2
    int method=asInteger(_method);// 0,1,2,3,4
    double *X0=REAL(_X), *Y0=REAL(_Y);  
    int *comb=INTEGER(_comb);
    double ind;
      
    double *X = malloc(m*sizeof(double));
    double *Y = malloc(n*sizeof(double));
    double *xy = malloc(N*sizeof(double));
    double *unique = malloc(N*sizeof(double));
    int *nties = malloc(N*sizeof(int));
        
    int mc_rep=asInteger(_mc_rep); // 0: exact perm, 1: z only, 1e4: mc
    int nperm=length(_comb)/m; // number of permutation
    SEXP _ans=PROTECT(allocVector(REALSXP, mc_rep==0?nperm:mc_rep));
    double *ans=REAL(_ans);
    
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
        ans[0]=compute_pair_wmw_Z(X0, Y0, xy, m, n, corr, method, adj, 0); 
    } else {
        if (mc_rep>1) {
            // Monte Carlo
            for (b=0; b<mc_rep; b++) {                
                for (i = 0; i < m; i++) {
                    ind=RUNIF(0.0,1.0);
                    X[i]=ind< 0.5?X0[i]:Y0[i];
                    Y[i]=ind>=0.5?X0[i]:Y0[i];
                }
                ans[b]=compute_pair_wmw_Z(X, Y, xy, m, n, corr, method, adj, 0);        
            }    
        } else {
            // exact   
         	int c=0;
            for (b=0; b<nperm; b++) {
                //PRINTF("%i ", b); 
                for (i = 0; i < m; i++) {
                    ind=comb[c++];   
                    X[i]=ind==0?X0[i]:Y0[i];
                    Y[i]=ind==1?X0[i]:Y0[i];
                }
                ans[b]=compute_pair_wmw_Z(X, Y, xy, m, n, corr, method, adj, 0);        
            }    
        }
        //for (b=0; b<mc_rep; b++) PRINTF("%f ", ans[b]); PRINTF("\n");
    }
            
    free(X); free(Y); free(xy); 
    free(unique); free(nties);

    UNPROTECT(1);
    return _ans;
}


// return just var
SEXP pair_wmw_var(SEXP _X, SEXP _Y, SEXP _method){
    
    int m=length(_X);
    int n=length(_Y);
    int method=asInteger(_method);// 0,1,2,3,4
    double *X=REAL(_X), *Y=REAL(_Y);  
      
    SEXP _ans=PROTECT(allocVector(REALSXP, 1));
    double *ans=REAL(_ans);
    
    double *xy = malloc((m+n)*sizeof(double));        
    ans[0]=compute_pair_wmw_Z(X, Y, xy, m, n, 0, method, 0, 1); 
    free(xy); 
    //free(unique); free(nties);

    UNPROTECT(1);
    return _ans;
}

