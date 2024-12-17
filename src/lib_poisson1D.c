/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"


/*la =3 */
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){

  
  if ((*kv)!=0){
  for (int i =0 ; i< (*la) ; i++ ){
    AB[i*(*lab)]  = 0;  
    printf("this is the value of i : %d\n" ,i);

  }
  }

  for (int i =0; i< (*la) ; i++ ){
    AB[i*(*lab)+1] =-1; 
    printf("this is the value of i : %d\n" ,i);

  }
  AB[(*kv)] =0 ; 
  // AB[(*la) ] =0 ; 

  for (int i =0 ; i< (*la) ; i++ ){
    AB[i*(*lab)+2] =2; 
    printf("this is the value of i : %d\n" ,i);

  }

  for (int i =0 /*derniere ligne remplie de -1*/ ; i< (*la)  ; i++ ){
    AB[i*(*lab)+3] =-1; 
    printf("this is the value of i : %d\n" ,i);
  }

  // AB[(*la)*(*kv)+3*(*la)] =0; 
  int last_element = (*lab ) * (*la) -1; 
  AB[last_element] =0; 

}




void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){


  
  if ((*kv)!=0){
  for (int i =0 ; i< (*la) ; i++ ){
    AB[i*(*lab)]  = 0;  
    printf("this is the value of i : %d\n" ,i);

  }
  }

  for (int i =0; i< (*la) ; i++ ){
    AB[i*(*lab)+1] =0; 
    printf("this is the value of i : %d\n" ,i);

  }
  AB[(*kv)] =0 ; 
  // AB[(*la) ] =0 ; 

  for (int i =0 ; i< (*la) ; i++ ){
    AB[i*(*lab)+2] =1; 
    printf("this is the value of i : %d\n" ,i);

  }

  for (int i =0 /*derniere ligne remplie de -1*/ ; i< (*la)  ; i++ ){
    AB[i*(*lab)+3] =0; 
    printf("this is the value of i : %d\n" ,i);
  }

  // AB[(*la)*(*kv)+3*(*la)] =0; 
  int last_element = (*lab ) * (*la) -1; 
  AB[last_element] =0; 



}/*make identity matrix */






void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  /*this function sets up the vector b */
  /*in this specific case there are no source terms so we should get 0 everywhere 
  except at the begining where we will have BC0 and at the end where we will have BC1 */
  /*considering we are dealing with Dirichlet boundary conditions*/

  for (int i = 0; i < *la; i++) {/* initializing everything to 0 might be uncessary */
    RHS[i] = 0; 
  }

  RHS[0] = *BC0;         // Left boundary condition BC0
  RHS[*la - 1] = *BC1;   // Right boundary condition BC1




}  /*Ax=b */

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){

  /*the analytical solution is */
  /*EX_SOL = BC0 +  (BC1 - BC0) *X /L  */
  double L = X[(*la)-1] - X[0]; 

  for (int i =0 ; i < (*la) ;i++){
    EX_SOL[i] = (*BC0)  +  ( (*BC1) - (*BC0) )* (X[i] - X[0])/L; 
  }

  /*EX_SOL[i] = BC0 + (X[i] - X[0]) / (X[la-1] - X[0]) * (BC1 - BC0)*/
    
}  




void set_grid_points_1D(double* x, int* la){
  
  /*this function role is to initialize the vector x as a uniformly spaced grid */

  /*the formula for computing the space between each point of the grid is */

  /*h=  L /(N-1) */
  /*here we assume L is 1 based on exercice 1 set up */

  double L = 1; 

  double h = L /((*la) -1); 

  for (int i =0 ; i < *la ; i ++){
    x[i] = h*i; 
  }
}

double relative_forward_error(double* x, double* y, int* la){
}

int indexABCol(int i, int j, int *lab){
  return 0;
}
/**/
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}/*to implement this one you need to implement all the ones before */
