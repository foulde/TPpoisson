/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include <time.h>


#define TRF 0
#define TRI 1
#define SV 2

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info = 1;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;

  double relres;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  // kv=0;
  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  double *Id; 
  Id = (double *) malloc(sizeof(double)*lab*la);
  set_GB_operator_colMajor_poisson1D_Id(Id, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(Id, &lab, &la, "ID.dat");

  
  // write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB.dat");

  printf("Solution with LAPACK\n");
  ipiv = (int *) calloc(la, sizeof(int));

  /* LU Factorization */
  if (IMPLEM == TRF) {
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  }

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  if (IMPLEM == TRI) {
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    // if (info!=0){printf("\n INFO dgbtridiag = %d\n",info);}


  }

  if (IMPLEM == TRI || IMPLEM == TRF){
    /* Solution (Triangular) */
    if (info==0){
      dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }else{
      printf("\n INFO = %d\n",info);
    }
  }

  /* It can also be solved with dgbsv */
  if (IMPLEM == SV) { /*exercice 5*/
    // TODO : use dgbsv
    // dgbsv_
    // dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);



    clock_t start, end; 
    double cpu_time_used;

    start = clock();
    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    /*return value in info*/
    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time used by DGBSV: %f seconds\n", cpu_time_used);


    /* we could plot the difference of time between dgbsv and the combination of dgbtrs_
    and dgbtrf_    */


  }

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  relres = relative_forward_error(RHS, EX_SOL, &la);
  
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
