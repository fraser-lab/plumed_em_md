/*

< qRMS.c >

 RMS calculation based using quaternion-rotation

 based on 
  Charles F.F. Karney "Quaternions in molecular modeling"
  E-print: arXiv:physics/0506177

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include "pdbstruct.h"

struct POS_ARRAY{
  int N;
  int Nmalloc;
  double **pos; /* [N][3] (malloc later) */
};

/*** FUNCTIONS (GLOBAL) ***/
double Calculate_CRMS_Bwn_Two_ATOMs();
double Calculate_CRMS_Bwn_Two_ATOMs_Using_CA_With_Same_Rnum();
double Rotate_ATOMs();

/*** FUNCTIONS (LOCAL) ***/
static void Malloc_POS_ARRAY();
static void Free_POS_ARRAY();
static void Cal_Optimal_Rmatrix_Quaternion();
static void Cal_Rmatrix_From_Quaternion();
static double Cal_RMS();
static void Jacobi_Wilkinson4();
static void Set_G_to_Zero();
static void Rotation();
static void find_max_abs4();
static void Find_Minimum_Eigen_Value4();
static void equal_matrix4();
static void prod_matrix4();
static void Mult_Mat_Vec3();
static void print_matrix4();



void Malloc_POS_ARRAY(PosArray,N)
  struct POS_ARRAY *PosArray;
  int N;
{
 int i;
 PosArray->N = N;
 PosArray->Nmalloc = N;
 PosArray->pos = (double **)malloc(sizeof(double *) * N); 
 for (i=0;i<N;++i) PosArray->pos[i] = (double *)malloc(sizeof(double)*3);
} /* end of Malloc_POS_ARRAY() */


void Free_POS_ARRAY(PosArray)
  struct POS_ARRAY *PosArray;
{
 int i;
 for (i=0;i<PosArray->Nmalloc;++i) free(PosArray->pos[i]);
 free(PosArray->pos);
} /* end of Free_POS_ARRAY() */



double Calculate_CRMS_Bwn_Two_ATOMs(atmheadA,atmheadB,gA,gB,Rmat)
 struct ATOM  *atmheadA,*atmheadB; 
 double gA[3],gB[3]; /* Center position */
 double Rmat[3][3];  /* Rotation Matrix */
{ 
  int j,n,NatomA, NatomB;
  double rms;
  struct ATOM *an,*bn;
  struct POS_ARRAY mch_posA,mch_posB;  /* matched atoms */
 
  printf("#Calculate_CRMS_Bwn_Two_ATOMs(atmheadA,atmheadB,gA,gB,Rmat)\n");
  NatomA = 0;
  an = atmheadA;
  while (an->next!=NULL){
    an = an->next;
    NatomA += 1; }
 
  NatomB = 0;
  bn = atmheadB;
  while (bn->next!=NULL){
    bn = bn->next;
    NatomB += 1; }

  if (NatomA != NatomB){
    printf("#ERROR:NatomA(%d) is not equal to NatomB(%d).\n",NatomA,NatomB);
    exit(1);
    gA[0] = gA[1] = gA[2] = gB[0] = gB[1] = gA[2] = 0.0;
    Rmat[0][0] = 1.0; Rmat[0][1] = 0.0; Rmat[0][2] = 0.0; 
    Rmat[1][0] = 0.0; Rmat[1][1] = 1.0; Rmat[1][2] = 0.0; 
    Rmat[2][0] = 0.0; Rmat[2][1] = 0.0; Rmat[2][2] = 1.0; 
    return(0.0);
  }

  Malloc_POS_ARRAY(&mch_posA,NatomA);
  Malloc_POS_ARRAY(&mch_posB,NatomB);

  an = atmheadA;
  bn = atmheadB;
  n = 0;
  while ((an->next!=NULL)&&(bn->next!=NULL)){
    an = an->next;
    bn = bn->next;
    for (j=0;j<3;++j){
      mch_posA.pos[n][j] = an->Pos[j]; 
      mch_posB.pos[n][j] = bn->Pos[j]; 
    }
    printf("#%d A %lf %lf %lf B %lf %lf %lf\n",n,
      mch_posA.pos[n][0], mch_posA.pos[n][1], mch_posA.pos[n][2],
      mch_posB.pos[n][0], mch_posB.pos[n][1], mch_posB.pos[n][2]);
    n += 1;
  }

  Set_G_to_Zero(&mch_posA,gA);
  Set_G_to_Zero(&mch_posB,gB);
  printf("#n %d gA %lf %lf %lf gB %lf %lf %lf\n",n,
gA[0],gA[1],gA[2], gB[0],gB[1],gB[2]);
  Cal_Optimal_Rmatrix_Quaternion(&mch_posA,&mch_posB,Rmat);
  Rotation(&mch_posA,Rmat);
  rms = (float)Cal_RMS(&mch_posA,&mch_posB);
  Free_POS_ARRAY(&mch_posB);
  Free_POS_ARRAY(&mch_posA);
  return(rms);
} /* end of Calculate_CRMS_Bwn_Two_ATOMs() */






double Calculate_CRMS_Bwn_Two_ATOMs_Using_CA_With_Same_Rnum(atmheadA,atmheadB,gA,gB,Rmat)
 struct ATOM  *atmheadA,*atmheadB; 
 double gA[3],gB[3]; /* Center position */
 double Rmat[3][3];  /* Rotation Matrix */
{ 
  int j,hit,NatomA, Natom;
  double rms;
  struct ATOM *an,*bn;
  struct POS_ARRAY mch_posA,mch_posB;  /* matched atoms */
 
  printf("#Calculate_CRMS_Bwn_Two_ATOMs_Using_CA_With_Same_Resi_Rnum(atmheadA,atmheadB,gA,gB,Rmat)\n");

 /** [1] Malloc mch_posA and mch_posB **/
  NatomA = 0;
  an = atmheadA;
  while (an->next!=NULL){
    an = an->next;
    if (strcmp(an->Atom," CA ")==0) {NatomA += 1;}
  }

  Malloc_POS_ARRAY(&mch_posA,NatomA);
  Malloc_POS_ARRAY(&mch_posB,NatomA);

 /** [2] Make mch_posA and mch_posB for CA atoms with the same Rnum**/

  Natom = 0;
  an = atmheadA;
  while (an->next!=NULL){
    an = an->next;
    if (strcmp(an->Atom," CA ")==0){
      bn = atmheadB;
      hit  = 0;
      while ((hit==0)&&(bn->next!=NULL)){
        bn = bn->next;
        if ((strcmp(bn->Atom," CA ")==0)&&(strcmp(an->Rnum,bn->Rnum)==0)){
           hit = 1;
           for (j=0;j<3;++j){
             mch_posA.pos[Natom][j] = an->Pos[j];
             mch_posB.pos[Natom][j] = bn->Pos[j];
           }
           Natom += 1;
        }
        
      } 
    }
  }

  mch_posA.N = Natom;
  mch_posB.N = Natom;

  if (Natom == 0){
    printf("#ERROR: two molecules have no common CA atoms with the same Rnum.\n");
    exit(1);
    gA[0] = gA[1] = gA[2] = gB[0] = gB[1] = gA[2] = 0.0;
    Rmat[0][0] = 1.0; Rmat[0][1] = 0.0; Rmat[0][2] = 0.0; 
    Rmat[1][0] = 0.0; Rmat[1][1] = 1.0; Rmat[1][2] = 0.0; 
    Rmat[2][0] = 0.0; Rmat[2][1] = 0.0; Rmat[2][2] = 1.0; 
    return(0.0);
  }

 /** [3] Calculate RMSD for mch_posA and mch_posB **/
  Set_G_to_Zero(&mch_posA,gA);
  Set_G_to_Zero(&mch_posB,gB);
  printf("#Natom %d gA %lf %lf %lf gB %lf %lf %lf\n",Natom, gA[0],gA[1],gA[2], gB[0],gB[1],gB[2]);
  Cal_Optimal_Rmatrix_Quaternion(&mch_posA,&mch_posB,Rmat);
  Rotation(&mch_posA,Rmat);
  rms = (float)Cal_RMS(&mch_posA,&mch_posB);
  Free_POS_ARRAY(&mch_posB);
  Free_POS_ARRAY(&mch_posA);
  return(rms);
} /* end of Calculate_CRMS_Bwn_Two_ATOMs_Using_CA_With_Same_Rnum() */









double Rotate_ATOMs(atmhead,gA,gB,Rmat)
 struct ATOM *atmhead;
 double gA[3];
 double gB[3];
 double Rmat[3][3];
 /**
  return [rmsd between original and rotated molecules].
  **/ 
{
 int j,Natom;
 double newpos[3],rpos[3];
 double rmsd;
 struct ATOM *an;

 an = atmhead;
 rmsd = 0.0;
 Natom = 0;
 while (an->next != NULL){
   an = an->next;
   for (j=0;j<3;++j) an->Pos[j] = an->Pos[j] - gA[j];

   Mult_Mat_Vec3(rpos,Rmat,an->Pos);
   for (j=0;j<3;++j){  
     newpos[j] = rpos[j] + gB[j];
     rmsd += (newpos[j]-an->Pos[j])*(newpos[j]-an->Pos[j]);
     an->Pos[j] = newpos[j];
   }
   Natom += 1;
 }
 return(sqrt(rmsd/Natom));
 
} /* end of Rotate_ATOMs() */




void Set_G_to_Zero(parray,G)
 struct POS_ARRAY *parray;
 double G[3];
{
 int i,j;
 double sum[3];
 
 printf("#Set_G_to_Zero(parray N %d Nmalloc %d)\n",parray->N,parray->Nmalloc);
 sum[0] = sum[1] = sum[2] = 0.0;
 for (i=0;i<parray->N;++i){ 
   for (j=0;j<3;++j){sum[j] += parray->pos[i][j];}
  } 

 for (j=0;j<3;++j){G[j] = sum[j]/parray->N;}

 for (i=0;i<parray->N;++i){
  for (j=0;j<3;++j){parray->pos[i][j] = parray->pos[i][j] - G[j];}
 }
} /* end of Set_G_to_Zero() */


void Rotation(parray,mat)
 struct POS_ARRAY *parray; 
 double mat[3][3];
{
 double p[3],q[3];
 int i;

 for (i=0;i<parray->N;++i){ 
   p[0] = parray->pos[i][0];  
   p[1] = parray->pos[i][1];  
   p[2] = parray->pos[i][2];  
   Mult_Mat_Vec3(q,mat,p);
   parray->pos[i][0] = q[0]; 
   parray->pos[i][1] = q[1]; 
   parray->pos[i][2] = q[2];
 }

} /* end of Rotation() */




double Cal_RMS(pA,pB)
 struct POS_ARRAY *pA,*pB;
{ int i;
  double dx,dy,dz;
  double RM;
  double rms,dd,d;

 RM = 0.0; 
 for (i=0;i<pA->N;++i) 
 { dx = pA->pos[i][0] - pB->pos[i][0]; 
   dy = pA->pos[i][1] - pB->pos[i][1]; 
   dz = pA->pos[i][2] - pB->pos[i][2]; 
   dd = (dx*dx)+(dy*dy)+(dz*dz);
   if (dd>0.0) d = sqrt(dd); else d = 0.0;
   RM += dd;
  }

 if (RM>0.0) rms = sqrt(RM/pA->N); else rms = 0.0;
 return(rms);
} /* end of Cal_RMS() */






void Cal_Optimal_Rmatrix_Quaternion(pA,pB,R)
 struct POS_ARRAY *pA,*pB; 
 double R[3][3];
{
 int i,j,k;
 double B[4][4],E[4][4],V[4][4];
 double a[3],b[3],aa[3],bb[3];
 double min_eval,min_evec[4];

/*
 printf("#int Cal_Optimal_Rmatrix_Quaternion(pA,pB,R)\n");
*/

 /**** (1) Make Matrix B  *******/ 
 /*
   a = yk+xk
   b = yk-xk
   Ak = |0.0 -b0 -b1  -b2|
        |b0  0.0 -a2   a1|
        |b1  a2  0.0  -a0|
        |b2 -a1  a0   0.0|
  
   Bk = tra[Ak] * Ak 
      = |0.0  b0  b1    b2| |0.0 -b0 -b1  -b2|
        |-b0  0.0 a2   -a1| |b0  0.0 -a2   a1|
        |-b1 -a2  0.0   a0| |b1  a2  0.0  -a0|
        |-b2  a1  -a0  0.0| |b2 -a1  a0   0.0|

      = |b0b0+b1*b1+b2*b2 a2b1-a1b2        -a2b0+a0b2       a1b0-a0b1       |
        |a2b1-a1b2        b0b0+a1*a1+b2*b2 b0b1-a0a1        b0b2-a0a2       |
        |-a2b0+a0b2       b0b1-a0a1        a0a0+b1*b1+a2*a2 b1b2-a1a2       |
        | a1b0-a0b1       b0b2-a0a2        b1b2-a1a2        a0a0+a1*a1+b2*b2|

   B = 1/N * \sum_{k}^{N} Bk

*/ 

 for (i=0;i<4;++i)
  for (j=0;j<4;++j) B[i][j] = 0.0;


 for (k=0;k<pA->N;++k){
 
  for (i=0;i<3;++i) {
   a[i] = pB->pos[k][i] + pA->pos[k][i];  
   b[i] = pB->pos[k][i] - pA->pos[k][i]; 
   aa[i] = a[i]*a[i];
   bb[i] = b[i]*b[i]; }  

  B[0][0] += bb[0]+bb[1]+bb[2];
  B[0][1] +=  a[2]*b[1] - a[1]*b[2];
  B[0][2] += -a[2]*b[0] + a[0]*b[2];
  B[0][3] +=  a[1]*b[0] - a[0]*b[1];
  
  B[1][1] += bb[0]+aa[1]+aa[2];
  B[1][2] +=  b[0]*b[1] - a[0]*a[1];
  B[1][3] +=  b[0]*b[2] - a[0]*a[2];
 
  B[2][2] += aa[0]+bb[1]+aa[2];
  B[2][3] +=  b[1]*b[2] - a[1]*a[2];
  B[3][3] += aa[0]+aa[1]+bb[2];

 } /* k */

  for (i=0;i<4;++i)
   for (j=i;j<4;++j) 
    { B[i][j] /= pA->N;
      B[j][i] = B[i][j]; }

 /**** (2) Calculate Eigen Value/Vector of Matrix B  *******/ 
 equal_matrix4(E,B); 
 Jacobi_Wilkinson4(E,V);
 Find_Minimum_Eigen_Value4(E,V,&min_eval,min_evec);
 
 /**** (3) Calculate Rmat from the eigen vector quaternion ****/
 Cal_Rmatrix_From_Quaternion(R,min_evec);

 /* print_matrix3(R,"Rmat"); */

} /* end of Cal_Optimal_Rmatrix_Quaternion() */







void Find_Minimum_Eigen_Value4(E,V,min_eval,min_evec)
 double E[4][4]; /* Diagonal matrix for eigen values */
 double V[4][4]; /* Matrix for eigen vectors */
 double *min_eval;   /* minimum_evalue to be calulated */
 double min_evec[4];   /* evector for minimum_evalue to be calulated */
{
 int i,min_i;
 double ev[4];

 min_i = 0;
 for (i=0;i<4;++i) ev[i] = E[i][i];
      if ((ev[0]<=ev[1])&&(ev[0]<=ev[2])&&(ev[0]<=ev[3])) min_i = 0;
 else if ((ev[1]<=ev[0])&&(ev[1]<=ev[2])&&(ev[1]<=ev[3])) min_i = 1;
 else if ((ev[2]<=ev[0])&&(ev[2]<=ev[1])&&(ev[2]<=ev[3])) min_i = 2;
 else if ((ev[3]<=ev[0])&&(ev[3]<=ev[1])&&(ev[3]<=ev[2])) min_i = 3;

 *min_eval = ev[min_i];
 for (i=0;i<4;++i) min_evec[i] = V[i][min_i];
} /* end of Find_Minimum_Eigen_Value4() */



void Cal_Rmatrix_From_Quaternion(R,q)
 double R[3][3],q[4];
{
 int i,j;
 double Q[4][4];

 for (i=0;i<4;++i) 
  for (j=i;j<4;++j) Q[i][j] = Q[j][i] = q[i]*q[j]; 
 
 R[0][0] = 2.0*(Q[0][0]+Q[1][1])-1.0;
 R[0][1] = 2.0*(Q[1][2]-Q[0][3]);
 R[0][2] = 2.0*(Q[1][3]+Q[0][2]);
 
 R[1][0] = 2.0*(Q[1][2]+Q[0][3]);
 R[1][1] = 2.0*(Q[0][0]+Q[2][2])-1.0;
 R[1][2] = 2.0*(Q[2][3]-Q[0][1]);

 R[2][0] = 2.0*(Q[1][3]-Q[0][2]);
 R[2][1] = 2.0*(Q[2][3]+Q[0][1]);
 R[2][2] = 2.0*(Q[0][0]+Q[3][3])-1.0;

 /*
 R[0][0] = 1.0-2.0*(Q[2][2]+Q[3][3]);
 R[0][1] = 2.0*(Q[1][2]+Q[0][3]);
 R[0][2] = 2.0*(Q[1][3]-Q[0][2]);
 
 R[1][0] = 2.0*(Q[1][2]-Q[0][3]);
 R[1][1] = 1.0-2.0*(Q[1][1]+Q[3][3]);
 R[1][2] = 2.0*(Q[2][3]+Q[0][1]);

 R[2][0] = 2.0*(Q[1][3]+Q[0][2]);
 R[2][1] = 2.0*(Q[2][3]-Q[0][1]);
 R[2][2] = 1.0-2.0*(Q[1][1]+Q[2][2]);
 */

} /* end of Cal_Rmatrix_From_Quaternion() */



void Jacobi_Wilkinson4(A,U)
 double A[4][4]; /* Input matrix whose eigen value/vector will be calculated. 
      After calculation, it becomes diagonal matrix of eigen values. */
 double U[4][4]; /* Output matrix corresponding to eigen vectors. 
   Eigen vector corresponds to each "column" vector of matrix U.
   For example k-th eigen vector 
   evec_k[0] =  U[0][k]; 
   evec_k[1] =  U[1][k]; 
   evec_k[2] =  U[2][k]; 
   evec_k[3] =  U[3][k]; 
   */ 
{ 
 double R[4][4],TR[4][4],BUFF[4][4];
 double max,co,si;
 double wa,sa,r,sqrt2;
 int i,j,mi,mj,c;

 sqrt2 = sqrt(2.0);

   /* --------SETTING U[][] = I --------*/
   
   for (i=0;i<4;++i)
     for (j=0;j<4;++j) if (i==j) U[i][j] = 1.0;  else U[i][j] = 0.0;
    
    find_max_abs4(&mi,&mj,&max,A);
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
      for (i=0;i<4;++i)
	for (j=0;j<4;++j) if (i==j)  R[i][j] = 1.0;  else R[i][j] =0.0;
     
       R[mi][mi] =   co; R[mi][mj] = -si;
       R[mj][mi] =   si; R[mj][mj] = co;

      for (i=0;i<4;++i)
	for (j=0;j<4;++j) TR[j][i] = R[i][j]; 

       prod_matrix4(BUFF,A,R);
       prod_matrix4(A,TR,BUFF);
       prod_matrix4(BUFF,U,R);
       equal_matrix4(U,BUFF);

       find_max_abs4(&mi,&mj,&max,A);
       ++c;
 
      } /* end of while */

}/* end of Jacobi_Wilkinson() */




void find_max_abs4(mi,mj,max,A) /* Max_Abs element A[mi][mj] = max */
 int *mi,*mj;
 double *max,A[4][4];
{ int i,j,p,q;

 *max =0.0;
  p = q = 0;
  for (i=0;i<4;++i)
    for (j=0;j<4;++j)
      if (i!=j)
       if (fabs(A[i][j])> *max)
	      { *max = fabs(A[i][j]);     p = i; q= j;}  
   *mi =p; *mj = q;
} /* end of find_max_abs4() */


void equal_matrix4(A,B)    /* A = B */
 double A[4][4],B[4][4];
{ int i,j;
  for (i=0;i<4;++i){
    for (j=0;j<4;++j)  A[i][j] = B[i][j];
   }
}/* end of equal_matrix4() */


void prod_matrix4(C,A,B) /* C = A * B */
 double C[4][4],A[4][4],B[4][4];
{ int i,j,k;
  double sum;
  for (i=0;i<4;++i){
    for (j=0;j<4;++j){ 
       sum =0;
       for (k=0;k<4;++k)  sum = sum + A[i][k] * B[k][j];
       C[i][j] = sum; }
   }
 } /* end of prod_matrix4() */


void Mult_Mat_Vec3(y,mat,x)
 double y[3],mat[3][3],x[3];  /* This function is only for 3-D */
{ 
  y[0] = mat[0][0]*x[0] + mat[0][1]*x[1] + mat[0][2] * x[2];
  y[1] = mat[1][0]*x[0] + mat[1][1]*x[1] + mat[1][2] * x[2];
  y[2] = mat[2][0]*x[0] + mat[2][1]*x[1] + mat[2][2] * x[2]; 
}


void print_matrix4(A,comment)
 double A[4][4];
 char *comment;
{ int i,j;
  printf("\n");
  printf("#%s\n",comment);
  for (i=0;i<4;++i){
    for (j=0;j<4;++j){ 
       printf(" %8.3lf",A[i][j]);
       if (j==3)  printf(";\n"); else printf(",");
    }
  }
  printf("\n");
}/* end of print_matrix4() */


