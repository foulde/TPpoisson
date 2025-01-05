/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
  int n= *la; 
  for (int i =0 ; i <*la ; i++){
    eigval[i] = 4*pow(sin ((i+1)*M_PI/(2.0*(n+1))),2); 
  } 
}

double eigmax_poisson1D(int *la){
  double * eigval  = malloc((*la) *sizeof(double )); 
  eig_poisson1D(eigval , la);
  double max =eigval[0]; 
  for(int i =0 ; i<*la ; i++){
    max = fmax(eigval[i] , max); 
  }

  // return 0;
  free(eigval);
  return max ;
}

double eigmin_poisson1D(int *la){

  double * eigval  = malloc((*la) *sizeof(double )); 
  eig_poisson1D(eigval , la);

  double min =eigval[0]; 
  for(int i =0 ; i<*la ; i++){
    min = fmin(eigval[i] , min); 
  }

  free(eigval);
  return min ;

}

double richardson_alpha_opt(int *la){

  double max =  eigmax_poisson1D( la); 
  double min =  eigmin_poisson1D( la); 
  
  return 2.0/(max+min);

  // return 0;
  // return 2.0/(eigmax_poisson1D(la) + eigminpoisson1D(la));
}









// void richardson_alpha(double *AB, double RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku/*1*/, int*kl/*1*/, double *tol/*tolerance */, int *maxit /*max iteration */, double *resvec/*norm of the */, int *nbite){
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit , double *resvec, int *nbite){
  *nbite = 0;
  double * residus = malloc((*la) * sizeof(double)); 
  // double *temp = malloc(*la sizeof(double)); 
  cblas_dcopy(*la, residus,1,RHS ,1 );  /*on stocke b dans le residus de la fonction
  et après on soustrera AX*/

  double nresidus = cblas_dnrm2(*la,residus , 1); 
  cblas_dgbmv(CblasColMajor , CblasNoTrans ,*la , *la , *kl , *ku, -1.0 , AB , *lab, X ,1 ,1.0 ,residus,1); 

  double norm_res = cblas_dnrm2(*la , residus , 1); 
  double res = norm_res/nresidus; 

  while (res> *tol  && *nbite<*maxit)
  {
    cblas_daxpy(*la , *alpha_rich , residus ,1,X  ,1); 



    ///////////////////////////////////////////////////////////

    cblas_dcopy(*la, RHS, 1, residus, 1);

    // cblas_dcopy(*la , residus,1,RHS,1);  /*on stocke b dans le residus de la fonction
  // et après on soustrera AX*/

    nresidus = cblas_dnrm2(*la ,residus , 1); 
    cblas_dgbmv(CblasColMajor , CblasNoTrans ,*la , *la , *kl , *ku, -1.0 , AB , *lab, X ,1 ,1.0 ,residus,1); 

    norm_res = cblas_dnrm2(*la ,residus , 1); 
    res = norm_res/nresidus; 

    resvec[*nbite] = res; /*we store the error in here to later plot the history */

    ///////////////////////////////////////////////////////////


    (*nbite) ++; 

  }

  free (residus); 
     
}


//   // /*xkp1= xk + alpha(b-Ax)      */ 
//   // /*X -> X + alpha*(b-A*X)      */ 
//   // *nbite =0; 
//   // for (int iteration =0 ; iteration < *maxit ; i++){
//   //   /* appel function b-A*X  */
//   //   /* appel function RHS-A*X  */
    
    
//   //   for (int i =0 ; i <*la ; i++){
//   //     for(int j =0 ; j<){

//   //     }
//   //     X[i] = 

//   //   }

//   //   *nbite ++; 

//   // }



// // xkp1 = xk




















void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

