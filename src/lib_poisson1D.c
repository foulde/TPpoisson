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
  /*in this case x is the exact solution and y is the prediction */
  /*we set up : */
  /*Relative_forward_error = ||x -xpred || / ||x||*/


  double sumofdif =0; /* we calculate ||x -xpred||*/
  double sumofx =0; /* we calculate ||x||*/
  
  // double sumofx

  for(int i =0  ; i < *la ; i++){
    sumofdif += (x[i]-y[i])*(x[i]-y[i]); 
    sumofx += x[i]*x[i]; 
  }

  sumofdif = sqrt(sumofdif);
  sumofx = sqrt(sumofx);

  if (sumofx==0){
    printf("The vector x provided is null \n"); 
    return -1;
  }


  printf("this is the diff error %f\n", sumofdif);
  return sumofdif/sumofx; /*error calcul*/

}

int indexABCol(int i, int j, int *lab){

  int index=0; 
  /*The role of this function is too return the position of the value of index i,j */
  /*in a Colmajor matrix */

  /*we have to account for multiple cases */

  /*CASE1 i,j is not on the tridiagonal*/

  if (abs(i-j) > 1){
    return 0; /*we chose 0 as in AB the value at position 0 is also 0 so looking at 
    this place for value is the same as looking at other place not on the triadiagonals */
  }

  /*CASE2 i,j is on the tridiagonal*/

  /*CASE UPER diagonal */
  if (j == i + 1)        /*we multiply by 3 but *lab works too */
    // index = 3 * j - 1;
    index = (*lab) * j - 1;

  /*CASE main diagonal */
  else if(j==i){
    // index = 3*j
    index = (*lab) * j ;
    
  } 

  else if (j== i-1){
    // index = 3 * j + 1;
    index = (*lab) * j + 1; 

  }


  return index;
}



int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info) {
  double factor;
  *info = 0;

  // Validate that the matrix is tridiagonal
  if ((*kl)!= 1 ||  (*ku) != 1) {
      *info = -1;
      printf("\nthis matrix is not tridiag\n");
      return -1;
  }

  // LU factorization loop
  for (int i = 0; i < *n - 1; i++) {

    // Get diagonal and subdiagonal indices using indexABCol
    int diag_index = indexABCol(i, i, lab);           // Main diagonal A[i, i]
    int subdiag_index = indexABCol(i + 1, i, lab);    // Subdiagonal A[i + 1, i]
    int next_diag_index = indexABCol(i + 1, i + 1, lab); // Next diagonal A[i + 1, i + 1]

    // Check for zero pivot
    if (AB[diag_index] == 0.0) {
        *info = i + 1; // Report the problematic pivot
        return -1;
        }

    //   factor computation
    factor = AB[subdiag_index] / AB[diag_index];

    // we put  the factor in the  subdiagonal 
    AB[subdiag_index] = factor;

      // Update the diagonal of the next row
      int superdiag_index = indexABCol(i, i + 1, lab); // Superdiagonal A[i, i+1]
      AB[next_diag_index] -= factor * AB[superdiag_index];
  }

  // Initialize the pivot array
  for (int i = 0; i < *n; i++) {
      ipiv[i] = i + 1;
  }

  return 0; 
}

/*en methode de validation on peut comparer avec l'erreur avant et quand on descend en dessous de 1e-6
on valide */