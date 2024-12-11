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

  // AB[(*kv)*(*la)] =0 ; /* car */

  // for (int i =(*la)+1 ; i< (*la)*(*kv)+(*la) ; i++ ){

  for (int i =0; i< (*la) ; i++ ){
    AB[i*(*lab)+1] =-1; 
    printf("this is the value of i : %d\n" ,i);
    // printf("this is the value of i : %d" ,i);
    // printf("this is the value of izzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz : %d" ,i);

  }
  AB[(*kv)] =0 ; 
  // AB[(*la) ] =0 ; 

  // for (int i =(*la)*(*kv)+(*la) ; i< (*la)*(*kv)+2*(*la) ; i++ ){
  for (int i =0 ; i< (*la) ; i++ ){
    AB[i*(*lab)+2] =2; 
    printf("this is the value of i : %d\n" ,i);

  }



  // for (int i =(*la)*(*kv)+2*(*la) /*derniere ligne remplie de -1*/ ; i< (*la)*(*kv)+3*(*la)-1  ; i++ ){
  for (int i =0 /*derniere ligne remplie de -1*/ ; i< (*la)  ; i++ ){
    AB[i*(*lab)+3] =-1; 
    printf("this is the value of i : %d\n" ,i);

  }
  // AB[(*la)*(*kv)+3*(*la)] =0; 
  int last_element = (*lab ) * (*la) -1; 
  AB[last_element] =0; 

}




void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){


}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
}  

void set_grid_points_1D(double* x, int* la){
}

double relative_forward_error(double* x, double* y, int* la){
}

int indexABCol(int i, int j, int *lab){
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}
