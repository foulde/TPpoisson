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
    res = norm_res/nresidus; /*we calculate the relative error */

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

  printf("this is kv %d \n" , *kv); 

  /* we set MB to 0 */
    for(int j = 0; j < (*la)*(*lab); j++) {
        MB[j] = 0.0;
    }

    /*get the diag elements and invert them  in MB*/
    for(int i = 0; i < *la; i++) {
        double diag = AB[*kv + i * (*lab)];/*temp element*/
        MB[*kv + i * (*lab)] = 1.0 / diag;  
    }


}








double* invert_lower_tri(const double *MB, int la, int ldMB, int kl)
{
    /*We will allocate a dense NxN array "Minv" in col-major.*/
    /*So Minv has size la * la. The (i,j) element is Minv[i + j*la].*/
    double *Minv = (double*)calloc(la * la, sizeof(double));
    if (!Minv) {
        fprintf(stderr, "Memory allocation failed for Minv\n");
        return NULL;
    }

    /* For each column j in [0..la-1], solve M * z = e_j and store z in col j of Minv.*/
    for (int j = 0; j < la; j++)
    {
        
        double *rhs = (double*)calloc(la, sizeof(double));
        if (!rhs) {
            fprintf(stderr, "Memory allocation failed for rhs\n");
            free(Minv);
            return NULL;
        }
        rhs[j] = 1.0;  


        cblas_dtbsv(CblasColMajor,
                    CblasLower,
                    CblasNoTrans,
                    CblasNonUnit,
                    la, kl,
                    MB, ldMB,
                    rhs, 1);

        for (int i = 0; i < la; i++)
        {
            Minv[i + j*la] = rhs[i];
        }

        free(rhs);
    }

    // Now Minv holds M^-1 in col-major.
    // You can either return it or do something else with it here.
    return Minv;
}










void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){


    int n    = *la;       // matrix size
    int ldA  = *lab;      // leading dimension in band form
    int diag = *kv;       // row index of main diagonal in AB, typically 1

    // 1) Zero out MB
    for (int i = 0; i < n * ldA; i++) {
        MB[i] = 0.0;
    }

    // 2) Copy the main diagonal
    for (int j = 0; j < n; j++) {
        MB[diag + j*ldA] = AB[diag + j*ldA];
    }

    // 3) Copy the subdiagonal
    //    If A is truly tridiagonal, we have subdiagonal data in AB[diag+1 + j*ldA]
    //    for j=0..(n-2). We keep the same sign as in A (which might already be "−E").
    if (*kl >= 1) {
        for (int j = 0; j < n - 1; j++) {
            MB[diag+1 + j*ldA] = AB[diag+1 + j*ldA];
        }
    }


    invert_lower_tri(MB, *la, *lab, *kl);


}















// // void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite);



// void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){


//   *nbite = 0;
//   double * residus = malloc((*la) * sizeof(double)); 
//   double *temp  = malloc((*la) * sizeof(double)); 

//   // double *temp = malloc(*la sizeof(double)); 
//   cblas_dcopy(*la, residus,1,RHS ,1 );  /*on stocke b dans le residus de la fonction
//   et après on soustrera AX*/

//   double nresidus = cblas_dnrm2(*la,residus , 1); 
//   cblas_dgbmv(CblasColMajor , CblasNoTrans ,*la , *la , *kl , *ku, -1.0 , AB , *lab, X ,1 ,1.0 ,residus,1); 

//   double norm_res = cblas_dnrm2(*la , residus , 1); 
//   double res = norm_res/nresidus; 

//   while (res> *tol  && *nbite<*maxit)
//   {
//     // cblas_daxpy(*la , *alpha_rich , residus ,1,X  ,1);  
//   /* matrix computation of temp = MB * residus */
//     cblas_dgbmv(CblasColMajor, CblasNoTrans,
//                     *la, *la, *kl, *ku,
//                      1.0, MB, *lab, residus, 1,
//                      0.0, temp, 1);
    


//     // X = X + temp 
//     cblas_daxpy(*la, 1.0, temp, 1, X, 1);

//     // cblas_dcopy(*la, RHS, 1, residus, 1);
//     // /*we compute residus = b - A*x for the next iteration*/


//     // nresidus = cblas_dnrm2(*la ,residus , 1); 

//     // cblas_dgbmv(CblasColMajor , CblasNoTrans ,*la , *la , *kl , *ku, -1.0 , AB , *lab, X ,1 ,1.0 ,residus,1); 





//     ///////////////////////////////////////////////////////////

//     cblas_dcopy(*la, RHS, 1, residus, 1);



//     nresidus = cblas_dnrm2(*la ,residus , 1); 
//     cblas_dgbmv(CblasColMajor , CblasNoTrans ,*la , *la , *kl , *ku, -1.0 , AB , *lab, X ,1 ,1.0 ,residus,1); 

//     norm_res = cblas_dnrm2(*la ,residus , 1); 
//     res = norm_res/nresidus; /*we calculate the relative error */

//     resvec[*nbite] = res; /*we store the error in here to later plot the history */

//     ///////////////////////////////////////////////////////////


//     (*nbite) ++; 

//   }

//   free (residus); 
     


// }






void richardson_MB(double *AB, double *RHS, double *X, double *MB,
                   int *lab, int *la, int *ku, int *kl,
                   double *tol, int *maxit, double *resvec, int *nbite)
{
    *nbite = 0;

    double *residus = malloc((*la) * sizeof(double));
    double *temp    = malloc((*la) * sizeof(double));

    // Initialize residus = b
    // (Previously the arguments to cblas_dcopy were reversed.)
    cblas_dcopy(*la, RHS, 1, residus, 1);

    // Compute the norm of b to use for the relative residual denominator
    double nresidus = cblas_dnrm2(*la, residus, 1);

    // If b is the zero vector, then nresidus=0 => subsequent division would be invalid.
    // We'll handle that by skipping the relative measure or setting res=0 if b=0.
    double norm_res = 0.0;
    double res      = 0.0;

    if (nresidus > 1.0e-30)
    {
        // residus = b - A*x
        cblas_dgbmv(CblasColMajor, CblasNoTrans,
                    *la, *la, *kl, *ku,
                    -1.0, AB, *lab, X, 1,
                     1.0, residus, 1);

        norm_res = cblas_dnrm2(*la, residus, 1);
        res      = norm_res / nresidus;
    }
    else
    {
        // If b is effectively zero, we can treat the residual ratio as the absolute residual
        // or consider that the system might be trivial.
        cblas_dgbmv(CblasColMajor, CblasNoTrans,
                    *la, *la, *kl, *ku,
                    -1.0, AB, *lab, X, 1,
                     1.0, residus, 1);
        norm_res = cblas_dnrm2(*la, residus, 1);
        res      = norm_res; 
    }

    // Main iteration loop
    while ((res > *tol) && (*nbite < *maxit))
    {
        // temp = MB * residus
        cblas_dgbmv(CblasColMajor, CblasNoTrans,
                    *la, *la, *kl, *ku,
                     1.0, MB, *lab, residus, 1,
                     0.0, temp, 1);

        // X = X + temp
        cblas_daxpy(*la, 1.0, temp, 1, X, 1);

        // Recompute residual: residus = b
        cblas_dcopy(*la, RHS, 1, residus, 1);

        nresidus = cblas_dnrm2(*la, residus, 1);

        // residus = b - A*x
        cblas_dgbmv(CblasColMajor, CblasNoTrans,
                    *la, *la, *kl, *ku,
                    -1.0, AB, *lab, X, 1,
                     1.0, residus, 1);

        norm_res = cblas_dnrm2(*la, residus, 1);

        // Safely compute the relative residual
        if (nresidus > 1.0e-30)
            res = norm_res / nresidus;
        else
            res = norm_res;  // fallback if b is zero

        // Store the residual history
        resvec[*nbite] = res;

        (*nbite)++;
    }

    free(residus);
    free(temp);
}