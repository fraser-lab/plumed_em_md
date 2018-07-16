/*
 
   <Matrix3D.c>

   Various functions for dealing with 3x3 matrix.

   The function caclcularing inverse matrix using "Cramer's rule"
   is included.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Jacobi3.h" 

/** FUNCTIONS (GLOBAL) **/
double Dot_Prod3D();
void Cross_Prod3D();
void Add_Vec3D();
void Sub_Vec3D();
void Equal_Vec3D();
void Multiply_Scalar_Vec3D();
double Norm_Vec3D();
double Length_Vec3D();
double Normalize_Vec3D();
double Length_Normalize_Vec3D();
void Print_Vec3D();
int Cal_Inverse_Matrix3D_by_Cramer_Rule();
int Cal_Inverse_Matrix3D_by_Cramer_Rule_Symmetric();
void Set_Unit_Matrix3D();
void Equal_Matrix3D();
void Add_Matrix3D();
void Sub_Matrix3D();
void Multiply_Scalar_Matrix3D();
void Multiply_Scalar_Matrix3D_Self();
void Multiply_Scalar_Matrix3D_Symmetric(); 
void Multiply_Matrix3D();
void Multiply_Matrix3D_Vec3D();
void Multiply_Transpose_Matrix3D();
void Transform_RAtR_Matrix3D();
void Transform_RAtQ_Matrix3D();
void Multiply_Matrix_By_Vector3D();
void Transform_Rx_plus_tvec_3D();
void Transform_Rx_minus_x0_plus_tvec_3D();
void Transform_Rmat_around_tvec_3D();
void Print_Matrix3D();
void Transpose_Matrix3D();
void Make_Matrix3D_By_Prod_Vector3D();
double Determinant_Matrix3D();
double Quadratic_Form_3D();
void Cal_EigenVectors_For_Matrix3D();

/** FUNCTIONS (LOCAL) **/

double Dot_Prod3D(a,b)
 double a[3],b[3];
{
   return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
} /* end of Dot_Prod3D() */



void Cross_Prod3D(c,a,b) /* c = a x b */
 double c[3],a[3],b[3];
{
 c[0] = a[1]*b[2] - b[1]*a[2];
 c[1] = a[2]*b[0] - b[2]*a[0];
 c[2] = a[0]*b[1] - b[0]*a[1];
} /* end of Cross_Prod3D() */


void Add_Vec3D(c,a,b) /* c = a + b */
 double c[3], a[3], b[3];
{
 c[0] = a[0] + b[0];
 c[1] = a[1] + b[1];
 c[2] = a[2] + b[2];
} /* end of Sub_Vec3D() */


void Sub_Vec3D(c,a,b) /* c = a - b */
 double c[3], a[3], b[3];
{
 c[0] = a[0] - b[0];
 c[1] = a[1] - b[1];
 c[2] = a[2] - b[2];
} /* end of Sub_Vec3D() */


void Equal_Vec3D(a,b) /* a := b */
 double  a[3], b[3];
{
 a[0] =  b[0];
 a[1] =  b[1];
 a[2] =  b[2];
} /* end of Equal_Vec3D() */



void Multiply_Scalar_Vec3D(c,t,a) /* c = t * a */
 double c[3],t,a[3];
{
 c[0] = t * a[0];
 c[1] = t * a[1];
 c[2] = t * a[2];
} /* end of  Multiply_Scalar_Vec3D()  */


double Norm_Vec3D(a)
 double a[3];
{
 return(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
} /* end of Norm_Vec3D() */


double Length_Vec3D(a)
 double a[3];
{
 double len;
 len = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
 if (len>0.0) len = sqrt(len);
 return(len);
} /* end of Length_Vec3D() */


double Normalize_Vec3D(a)
 double a[3];
{
 /* return length of vector a */
 double len;

 len = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
 if (len>0.0){ 
   len = sqrt(len);
   a[0] /= len;  a[1] /= len;  a[2] /= len;
 }

 return(len);
} /* end of Normalize_Vec3D() */



double Length_Normalize_Vec3D(a,len_new)
 double a[3];
 double len_new;
{
 /* return length of vector a */
 double len,cons;
 len = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
 if (len>0.0){ 
   len = sqrt(len);
   cons = len_new/len; 
   a[0] *= cons;  a[1] *= cons;  a[2] *= cons; }
 return(len);

} /* end of Length_Normalize_Vec3D() */


void Print_Vec3D(a,title)
 double a[3];
 char *title;
{
 printf("#>%s %lf %lf %lf\n",title,a[0],a[1],a[2]);
} /* end of Print_Vec3D() */






void Set_Unit_Matrix3D(A) /* A := I */
 double A[3][3];
{
  A[0][0] = 1.0; A[0][1] = 0.0; A[0][2] = 0.0;
  A[1][0] = 0.0; A[1][1] = 1.0; A[1][2] = 0.0;
  A[2][0] = 0.0; A[2][1] = 0.0; A[2][2] = 1.0;

} /* end of Set_Unit_Matrix3D() */




void Equal_Matrix3D(A,B) /* A := B */
 double A[3][3],B[3][3];
{
 int i,j;
 for (i=0;i<3;++i){
   for (j=0;j<3;++j){A[i][j] = B[i][j];}
 }
} /* end of Equal_Matrix3D() */



void Multiply_Matrix3D(C,A,B)    /* C = A*B */
 double C[3][3],A[3][3],B[3][3];
{
 int i,j,k;

 for (i=0;i<3;++i){
  for (j=0;j<3;++j){
    C[i][j] = 0.0;
    for (k=0;k<3;++k){C[i][j] += A[i][k]*B[k][j];}
  }
 }

} /* end of Multiply_Matrix3D() */


void Multiply_Matrix3D_Vec3D(y,A,x)    /* y = A*x */
 double y[3],A[3][3],x[3];
{
 int i,k;

 for (i=0;i<3;++i){
   y[i] = 0.0;
   for (k=0;k<3;++k){y[i] += A[i][k]*x[k];}
 }

} /* end of Multiply_Matrix3D_Vec3D() */





void Multiply_Transpose_Matrix3D(C,A,B)    /* C = A*trans(B) */
 double C[3][3],A[3][3],B[3][3];
{
 int i,j,k;

 for (i=0;i<3;++i){
   for (j=0;j<3;++j){
     C[i][j] = 0.0;
     for (k=0;k<3;++k){C[i][j] += A[i][k]*B[j][k];}
   }
 }

} /* end of Multiply_Transpose_Matrix3D() */




void Transform_RAtR_Matrix3D(C,A,R)    /* C = R*A*trans(R) */
 double C[3][3],A[3][3],R[3][3];
{
 int i,j,k;
 double buff[3][3];

 /* (1) buff = R * A */
 for (i=0;i<3;++i){
  for (j=0;j<3;++j){
    buff[i][j] = 0.0;
    for (k=0;k<3;++k){buff[i][j] += R[i][k]*A[k][j];}
  }
 }

 /* (2) C = buff * tR **/
 for (i=0;i<3;++i){
  for (j=0;j<3;++j){
    C[i][j] = 0.0;
    for (k=0;k<3;++k){C[i][j] += buff[i][k]*R[j][k];}
  }
 }
} /* end of Transform_RAtR_Matrix3D() */



void Transform_RAtQ_Matrix3D(C,A,R,Q)    /* C = R*A*trans(Q) */
 double C[3][3],A[3][3],R[3][3],Q[3][3];
{
 int i,j,k;
 double buff[3][3];

 /* (1) buff = R * A */
 for (i=0;i<3;++i){
  for (j=0;j<3;++j){
    buff[i][j] = 0.0;
    for (k=0;k<3;++k){buff[i][j] += R[i][k]*A[k][j];}
  }
 }
/* (2) C = buff * tQ **/
 for (i=0;i<3;++i){
   for (j=0;j<3;++j){
     C[i][j] = 0.0;
     for (k=0;k<3;++k){C[i][j] += buff[i][k]*Q[j][k];}
   }
 }

} /* end of Transform_RAtQ_Matrix3D() */



void Multiply_Matrix_By_Vector3D(y,R,x)    /* y = R*x  */
 double y[3];
 double R[3][3];
 double x[3]; 
{
 int i,j;
 for (i=0;i<3;++i)
 { y[i] = 0.0;
   for (j=0;j<3;++j) y[i] += R[i][j]*x[j]; }

} /* end Multiply_Matrix_By_Vector3D()  */





void Transform_Rx_plus_tvec_3D(xT,x,R,t)    /* xT = R*x + t */
 double xT[3];
 double x[3]; 
 double R[3][3];
 double t[3];
{
 int i,j;

 for (i=0;i<3;++i){
   xT[i] = t[i];
   for (j=0;j<3;++j){xT[i] += R[i][j]*x[j];}
 }

} /* end of Transform_Rx_plus_tvec_3D() */



void Transform_Rx_minus_x0_plus_tvec_3D(xT,x,R,t,x0)    /* xT = R*(x-x0) + x0 + t */
 double xT[3];
 double x[3]; 
 double R[3][3];
 double t[3];
 double x0[3]; 
{
 int i,j;

 for (i=0;i<3;++i){
   xT[i] = x0[i]+t[i];
   for (j=0;j<3;++j) xT[i] += R[i][j]*(x[j]-x0[j]);
 }

} /* end of Transform_Rx_minus_x0_plus_tvec_3D() */


void Transform_Rmat_around_tvec_3D(xT,x,R,t)    /* xT = R*(x-t) + t */
 double xT[3];
 double x[3]; 
 double R[3][3];
 double t[3];
{
 int i,j;

 for (i=0;i<3;++i)
 {
   xT[i] = t[i];
   for (j=0;j<3;++j) xT[i] += R[i][j]*(x[j]-t[j]);
 }

} /* end of Transform_Rmat_around_tvec_3D() */




void Add_Matrix3D(C,A,B)   /* C = A + B */
 double C[3][3],A[3][3],B[3][3];
{
 int i,j;
 for (i=0;i<3;++i)
  for (j=0;j<3;++j) C[i][j] = A[i][j] + B[i][j];

} /* end of Add_Matrix3D() */


void Sub_Matrix3D(C,A,B)  /* C = A - B */ 
 double C[3][3],A[3][3],B[3][3];
{
 int i,j;
 for (i=0;i<3;++i)
  for (j=0;j<3;++j) C[i][j] = A[i][j] - B[i][j];

} /* end of Sub_Matrix3D() */


void Multiply_Scalar_Matrix3D(C,t,A)  /* C = t * A */
 double C[3][3],t, A[3][3];
{
 int i,j;
 for (i=0;i<3;++i){
   for (j=0;j<3;++j){C[i][j] = t * A[i][j];}
 }
} /* end of Multiply_Scalar_Matrix3D() */


void Multiply_Scalar_Matrix3D_Self(A,t)  /* A = t * A */
 double A[3][3],t;
{
 int i,j;
 for (i=0;i<3;++i){
   for (j=0;j<3;++j){A[i][j] = t * A[i][j];}
 }
} /* end of Multiply_Scalar_Matrix3D_Self() */

                                                                                         

void Multiply_Scalar_Matrix3D_Symmetric(C,t,A)  /* C = t * A */
 double C[3][3],t, A[3][3];
{
 C[0][0] = t * A[0][0];
 C[1][1] = t * A[1][1];
 C[2][2] = t * A[2][2];
 C[0][1] = C[1][0] = t * A[0][1];
 C[0][2] = C[2][0] = t * A[0][2];
 C[1][2] = C[2][1] = t * A[1][2];

} /* end of Multiply_Scalar_Matrix3D_Symmetric() */
                                                                                         


void Transpose_Matrix3D(tA,A) /* tA = transpose[A] */
 double tA[3][3],A[3][3];
{
 int i,j;

 for (i=0;i<3;++i){
   for (j=0;j<3;++j) {tA[i][j] = A[j][i];}
 }
} /* end of Transpose_Matrix3D() */


void Make_Matrix3D_By_Prod_Vector3D(A,b) /* A = b * transpose[b] */
 double A[3][3],b[3];
{
 int i,j;

 for (i=0;i<3;++i){
   for (j=0;j<3;++j){A[i][j] = b[i]*b[j];}
 }
} /* end of Make_Matrix3D_By_Prod_Vector3D() */






void Print_Matrix3D(A,title)
 double A[3][3];
 char *title;
{
 printf("#>%s\n",title);
 printf("#%lf %lf %lf\n",A[0][0],A[0][1],A[0][2]);
 printf("#%lf %lf %lf\n",A[1][0],A[1][1],A[1][2]);
 printf("#%lf %lf %lf\n",A[2][0],A[2][1],A[2][2]);
 printf("\n");

} /* end of Print_Matrix3D() */




double Determinant_Matrix3D(A)
 double A[3][3];
{
 double det;
 
 det =  A[0][0]*A[1][1]*A[2][2] 
      + A[0][2]*A[1][0]*A[2][1] 
      + A[0][1]*A[1][2]*A[2][0]
      - A[0][2]*A[1][1]*A[2][0]
      - A[0][0]*A[1][2]*A[2][1]
      - A[0][1]*A[1][0]*A[2][2];
 return(det);

} /* end of Determinant_Matrix3D() */


int Cal_Inverse_Matrix3D_by_Cramer_Rule(InvA,A,Det)
 double InvA[3][3];
 double A[3][3];
 double *Det;   /* Determinant of matrix A (to be calculated)*/
{
 /*
 
 <Cramer's Rule>

 Inv[A] = 1/|A| transpose[Delta_ij] 

 Delta_ij = (-1)**(i+j) * |Aij|

 Aij = matrix removing i-th row and j-th column.

 */ 
 int i,j;
 double det;


 det =  A[0][0]*A[1][1]*A[2][2] + A[0][2]*A[1][0]*A[2][1]
      + A[0][1]*A[1][2]*A[2][0] - A[0][2]*A[1][1]*A[2][0]
      - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2];

 InvA[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1]);
 InvA[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1]);
 InvA[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1]);

 InvA[1][0] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0]);
 InvA[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0]);
 InvA[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0]);

 InvA[2][0] =  (A[1][0]*A[2][1] - A[1][1]*A[2][0]);
 InvA[2][1] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0]);
 InvA[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0]);

 for (i=0;i<3;++i){
   for (j=0;j<3;++j){ 
     InvA[i][j] /= det;
   }
 }

 *Det = det;

 if (fabs(det)<0.00000000001){ return(0); }
 return(1); 

 /*
 printf("#detA = %lf %lf\n",det,*Det); 
 Print_Matrix3D(A,"A");
 Print_Matrix3D(InvA,"InvA");
 */

} /* end of  Cal_Inverse_Matrix3D_by_Cramer_Rule() */





int Cal_Inverse_Matrix3D_by_Cramer_Rule_Symmetric(InvA,A,Det)
 double InvA[3][3];
 double A[3][3];
 double *Det;   /* Determinant of matrix A */
{
 /*
 
   <Cramer's Rule>
 
   Inv[A] = 1/|A| transpose[Delta_ij]

   Delta_ij = (-1)**(i+j) * |Aij|
 
    Aij = matrix removing i-th row and j-th column.
 
 */
 double det;

 det =  A[0][0]*A[1][1]*A[2][2]
      + A[0][2]*A[1][0]*A[2][1]
      + A[0][1]*A[1][2]*A[2][0]
      - A[0][2]*A[1][1]*A[2][0]
      - A[0][0]*A[1][2]*A[2][1]
      - A[0][1]*A[1][0]*A[2][2];

 InvA[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1])/det;
 InvA[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1])/det;
 InvA[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1])/det;
 InvA[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0])/det;
 InvA[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0])/det;
 InvA[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0])/det;

 InvA[1][0] = InvA[0][1];
 InvA[2][0] = InvA[0][2];
 InvA[2][1] = InvA[1][2];

 *Det = det;
 if (fabs(det)<0.00000000001){ return(0); }
 return(1); 

} /* end of  Cal_Inverse_Matrix3D_by_Cramer_Rule_Symmetric() */


double Quadratic_Form_3D(A,v) /* transpose[v] * A * v */
 double A[3][3];  /* matrix */
 double v[3];     /* vector */
{
 double qform;
                                                                                         
 qform = A[0][0]*v[0]*v[0] + A[1][1]*v[1]*v[1] + A[2][2]*v[2]*v[2];
 qform += 2.0*(A[0][1]*v[0]*v[1] + A[0][2]*v[0]*v[2] + A[1][2]*v[1]*v[2]);
 return(qform);
                                                                                         
} /* end of Quadratic_Form_3D() */



void Cal_EigenVectors_For_Matrix3D(CovM,PCvar,PCaxis)
 double CovM[3][3];
 double PCvar[3];     /* eigen values  0:largest, 1:middle, 3:smallest. */
 double PCaxis[3][3]; /* eigen vectors [1st,2nd,3rd][x,y,z]*/
{
  double A[3][3],U[3][3],XxY[3];
  int ind[3],i,j;

/*
  printf("#Cal_EigenVectors_For_Matrix3D(CovM,PCvar,PCaxis)\n");
*/
 
  Equal_Matrix3D(A,CovM);
  Jacobi_Wilkinson3(A,U);
/*
  Print_Matrix3D(CovM,"CovM");
  Print_Matrix3D(A,"A");
  Print_Matrix3D(U,"U");
*/ 
 
/*
  int Jacobi_Wilkinson3(A,U)
    double A[3][3];  Input matrix whose eigen value/vector will be calculated.
                     After calculation, it becomes diagonal matrix of eigen values. 
    double U[3][3];  Output matrix corresponding to eigen vectors.
      Eigen vector corresponds to each "column" vector of matrix U.
       For example k-th eigen vector
        evec_k[0] =  U[0][k];
        evec_k[1] =  U[1][k];
        evec_k[2] =  U[2][k];
*/

  Sort_Eigen_Value3(A,ind,'D');
  for (i=0;i<3;++i){
    PCvar[i] = A[ind[i]][ind[i]];
    for (j=0;j<3;++j){
      PCaxis[i][j] = U[j][ind[i]]; 
    }
 }

  /*** Correlation for 'right-handed ***/
  Cross_Prod3D(XxY,PCaxis[0],PCaxis[1]);
  if (Dot_Prod3D(PCaxis[2],XxY)<0.0){
    for (i=0;i<3;++i){ PCaxis[2][i] = - PCaxis[2][i];}
  }

} /* end of Cal_EigenVectors_For_Matrix3D() */

