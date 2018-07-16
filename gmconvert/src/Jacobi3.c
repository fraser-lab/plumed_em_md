/*

<Jacobi3.c>

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>



/*** FUNCTIONS (GLOBAL) ***/
int Jacobi_Wilkinson3();
void Sort_Eigen_Value3();
void print_matrix3();

/*** FUNCTIONS (LOCAL) ***/
static void find_max_abs3();
static void prod_matrix3();
static void equal_matrix3();




void Sort_Eigen_Value3(E,ind,SortType)
 double E[3][3];  /* Diagonal matrix for eigen values (input) */
 int    ind[3];   /* Index for sorted array (to be calculated, not necessary to be initialized) */
 char   SortType; /* 'D'ecreasing, 'I'ncreasing */
{
 int i;
 double ev[3];

 for (i=0;i<3;++i) ev[i] = E[i][i];

 if (SortType == 'I'){
        if ((ev[0]<=ev[1])&&(ev[1]<=ev[2])) {ind[0] = 0; ind[1] = 1; ind[2] = 2;} 
   else if ((ev[0]<=ev[2])&&(ev[2]<=ev[1])) {ind[0] = 0; ind[1] = 2; ind[2] = 1;} 
   else if ((ev[1]<=ev[0])&&(ev[0]<=ev[2])) {ind[0] = 1; ind[1] = 0; ind[2] = 2;} 
   else if ((ev[1]<=ev[2])&&(ev[2]<=ev[0])) {ind[0] = 1; ind[1] = 2; ind[2] = 0;} 
   else if ((ev[2]<=ev[0])&&(ev[0]<=ev[1])) {ind[0] = 2; ind[1] = 0; ind[2] = 1;} 
   else if ((ev[2]<=ev[1])&&(ev[1]<=ev[0])) {ind[0] = 2; ind[1] = 1; ind[2] = 0;} 
 }
 else if (SortType == 'D'){
        if ((ev[0]>=ev[1])&&(ev[1]>=ev[2])) {ind[0] = 0; ind[1] = 1; ind[2] = 2;} 
   else if ((ev[0]>=ev[2])&&(ev[2]>=ev[1])) {ind[0] = 0; ind[1] = 2; ind[2] = 1;} 
   else if ((ev[1]>=ev[0])&&(ev[0]>=ev[2])) {ind[0] = 1; ind[1] = 0; ind[2] = 2;} 
   else if ((ev[1]>=ev[2])&&(ev[2]>=ev[0])) {ind[0] = 1; ind[1] = 2; ind[2] = 0;} 
   else if ((ev[2]>=ev[0])&&(ev[0]>=ev[1])) {ind[0] = 2; ind[1] = 0; ind[2] = 1;} 
   else if ((ev[2]>=ev[1])&&(ev[1]>=ev[0])) {ind[0] = 2; ind[1] = 1; ind[2] = 0;} 
 }

} /* end of Sort_Eigen_Value3() */




int Jacobi_Wilkinson3(A,U)
 double A[3][3]; /* Input matrix whose eigen value/vector will be calculated. 
      After calculation, it becomes diagonal matrix of eigen values. */
 double U[3][3]; /* Output matrix corresponding to eigen vectors. 
   Eigen vector corresponds to each "column" vector of matrix U.
   For example k-th eigen vector 
   evec_k[0] =  U[0][k]; 
   evec_k[1] =  U[1][k]; 
   evec_k[2] =  U[2][k]; 
   */ 
{ 
 double R[3][3],TR[3][3],BUFF[3][3];
 double max,co,si;
 double wa,sa,r,sqrt2;
 int i,j,mi,mj,c;

 sqrt2 = sqrt(2.0);

   /* --------SETTING U[][] = I --------*/
   
   for (i=0;i<3;++i){
     for (j=0;j<3;++j){ if (i==j) U[i][j] = 1.0;  else U[i][j] = 0.0;}
   } 
    find_max_abs3(&mi,&mj,&max,A);
    c = 0;
    
    while ((max>0.00000001)&&(c<100)){
      /*--------- SET of cos sin -----------*/ 
      
      wa = (A[mi][mi] + A[mj][mj])/2;  
      sa = (A[mi][mi] - A[mj][mj])/2;  
      r = sqrt(sa*sa + A[mi][mj]*A[mi][mj]);
      co =  sqrt(1.0+sa/r)/sqrt2;  si =  A[mi][mj]/2.0/r/co; 
      if (sa>0.0) 
       { co =  sqrt(1.0+sa/r)/sqrt2;  si =  A[mi][mj]/2.0/r/co; }
      else 
       { co =   sqrt(1.0-sa/r)/sqrt2; si = - A[mi][mj]/2.0/r/co; }
      
      /* -------- SET of Rot Matrix R and TR----- */  
      for (i=0;i<3;++i){
	for (j=0;j<3;++j){ if (i==j)  R[i][j] = 1.0;  else R[i][j] =0.0;}
       }
 
       R[mi][mi] =   co; R[mi][mj] = -si;
       R[mj][mi] =   si; R[mj][mj] = co;

      for (i=0;i<3;++i){
	for (j=0;j<3;++j){TR[j][i] = R[i][j];}
      }

      prod_matrix3(BUFF,A,R);
      prod_matrix3(A,TR,BUFF);
      prod_matrix3(BUFF,U,R);
      equal_matrix3(U,BUFF);

      find_max_abs3(&mi,&mj,&max,A);
      ++c;
 
     } /* end of while */

  return(1);
}/* end of Jacobi_Wilkinson() */




void find_max_abs3(mi,mj,max,A) /* Max_Abs element A[mi][mj] = max */
 int *mi,*mj;
 double *max,A[3][3];
{ int i,j,p,q;

  p = 0; q = 0;
  *max =0.0;
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      if (i!=j)
       if (fabs(A[i][j])> *max)
	      { *max = fabs(A[i][j]);     p = i; q= j;}  
    }
  }
   *mi =p; *mj = q;
  
} /* end of find_max_abs3() */

                                                                                                                        
                                                                                                                        
void prod_matrix3(C,A,B) /* C = A * B */
 double C[3][3],A[3][3],B[3][3];
{ int i,j,k;
  double sum;
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){ 
       sum =0;
       for (k=0;k<3;++k){ sum = sum + A[i][k] * B[k][j];}
       C[i][j] = sum;
    }
  }
 } /* end of prod_matrix3() */

                                                                                                                        
void equal_matrix3(A,B)    /* A = B */
 double A[3][3],B[3][3];
{ int i,j;
  for (i=0;i<3;++i) {
    for (j=0;j<3;++j) { A[i][j] = B[i][j];}
   }
}/* end of equal_matrix3() */


void print_matrix3(A,comment)
 double A[3][3];
 char *comment;
{ int i,j;
  printf("\n");
  printf("#%s\n",comment);
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){ 
      printf(" [%d%d]:%+8.3lf",i,j,A[i][j]);
      if (j==2)  printf("\n");
    }
  }
  printf("\n");
}/* end of print_matrix3() */
