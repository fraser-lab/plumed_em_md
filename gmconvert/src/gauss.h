/*
 
 <gauss.h>

 structure for Gaussian distribution 

 */

#define one_over_2pi_32  0.0634936359342410
#define two_pi_32       15.7496099457224190


struct GAUSS3D{
 int num;
 double M[3];         /* Mean Position */
 double CovM[3][3];   /* Covariance Matrix (Sigma) */
 double iCovM[3][3];  /* Inverse of Covariance Matrix      */
 double det;          /* Determinant of Covariance Matrix  */
 double Cons;         /* 1/{2pi)**3/2 * sqrt(det)} */
 double Weight;       /* Weight  */
 double PCvar[3];     /* Variant for PC axis for CovM.      0:largest, 1:middle, 3:smallest. */
 double PCaxis[3][3]; /* Principal Component Axes for CovM. [1st,2nd,3rd][x,y,z]*/
};
