/*
 
 <Ellipsoid.c>

 functions for output Ellipsoidal shape in VRML.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include "globalvar.h"
#include "gauss.h"
#include "Matrix3D.h"
#include "Jacobi3.h"

/*** FUNCTIONS (GLOBAL) ***/
void write__GAUSS3Ds_in_ellipsoid_VRML();
void Setup_PCvar_and_PCaxis();

/*** FUNCTIONS (LOCAL) ***/
static void Set_One_Gdf_for_GMM();
static double Set_QformThre_from_UpperCumuChiSquare();
static void Get_RotAxis_and_RotAngle_from_RotMatrix();
static void value_to_RGB();
static void find_min_max_Weight();

void write__GAUSS3Ds_in_ellipsoid_VRML(ofname,Ngauss, G, cover_ratio, color_type, RGBT)
  char   *ofname;
  int    Ngauss;
  struct GAUSS3D *G; /* array of GAUSS3D [0..Ngauss-1]*/
  double cover_ratio;  /* cover volume of ellipsoid for total gaussian density */
  char   color_type;   /* 'W':by weight, 'N':by gdf number, otherwise:given colorR,G,B,T */
  float RGBT[4];
{
  FILE *fpo;
  int g,c;
  double Rmat[3][3],Raxis[3],Rangle,qform_thre,dens_thre;
  int i,j;
  float rgbt[4],val;
  double minW,maxW;

  printf("#write__GAUSS3Ds_in_ellipsoid_VRML(color_type '%c' cover_ratio %lf RGBT %.2lf %.2lf %.2lf %.2lf) --> '%s'\n",color_type,cover_ratio,RGBT[0],RGBT[1],RGBT[2],RGBT[3],ofname);
 
 for (g=0;g<Ngauss;++g){
    Setup_PCvar_and_PCaxis(G[g].CovM,G[g].PCvar,G[g].PCaxis);
  }
 find_min_max_Weight(&minW,&maxW,Ngauss,G);

  /** make one gauss for GMM */
  qform_thre = 1.0;
  dens_thre  = 1.0;


  fpo = fopen(ofname,"w");
  fprintf(fpo,"#VRML V2.0 utf8\n");
  fprintf(fpo,"\n");
  fprintf(fpo,"\n");


  for (g=0;g<Ngauss;++g){
    qform_thre =  Set_QformThre_from_UpperCumuChiSquare(1.0-cover_ratio);

    for (i=0;i<3;++i){
      for (j=0;j<3;++j){
        Rmat[i][j] = G[g].PCaxis[j][i]; 
      }
    }
    Get_RotAxis_and_RotAngle_from_RotMatrix(Raxis,&Rangle,Rmat);

    fprintf(fpo,"Transform{\n");
    fprintf(fpo,"  translation %lf %lf %lf\n",G[g].M[0],G[g].M[1],G[g].M[2]);
    if (qform_thre>0){ 
      fprintf(fpo,"  scale  %lf %lf %lf\n",
       sqrt(G[g].PCvar[0]*qform_thre), sqrt(G[g].PCvar[1]*qform_thre), sqrt(G[g].PCvar[2]*qform_thre) );
    }
    else{
      fprintf(fpo,"  scale  1.0 1.0 1.0\n");
    }
    fprintf(fpo,"  rotation  %lf %lf %lf %lf\n",Raxis[0],Raxis[1],Raxis[2],Rangle);
    fprintf(fpo,"  children[\n");
    fprintf(fpo,"    Shape{\n");
    fprintf(fpo,"      appearance Appearance{\n");
    fprintf(fpo,"        material Material{\n");
  
  if (color_type=='W'){
      val = (G[g].Weight-minW)/(maxW-minW);
      val = (G[g].Weight-minW)/(maxW-minW);
      value_to_RGB(&rgbt[0], &rgbt[1], &rgbt[2], val, "BGR");
     /*
      printf("#Weight %lf minW %lf maxW %lf val %f rgbt %f %f %f %f\n",G[g].Weight,minW, maxW, val,rgbt[0],rgbt[1],rgbt[2],rgbt[3]);
      */
      rgbt[3] = RGBT[3];
    }
    else if (color_type=='N'){
      if (Ngauss>1){
        value_to_RGB(&rgbt[0], &rgbt[1], &rgbt[2], (float)g/(float)(Ngauss-1), "BGR");
      }
      else{
        value_to_RGB(&rgbt[0], &rgbt[1], &rgbt[2], 0.5, "BGR");
      } 
      rgbt[3] = RGBT[3];
    }
    else{
      for (c=0;c<4;++c){
        rgbt[c] = RGBT[c];
      }
    }

    fprintf(fpo,"          diffuseColor %f %f %f\n",rgbt[0], rgbt[1], rgbt[2]);
    fprintf(fpo,"          transparency %f\n",rgbt[3]);
    fprintf(fpo,"        }\n");
    fprintf(fpo,"      }\n");
    fprintf(fpo,"      geometry Sphere{\n");
    fprintf(fpo,"        radius 1.0\n");
    fprintf(fpo,"      }\n");
    fprintf(fpo,"    }\n");
    fprintf(fpo,"  ]\n");
    fprintf(fpo,"}\n");
  }
  fclose(fpo);

} /* end of write_GAUSS_MOLECULE_in_VRML_ellipsoid() */




void Set_One_Gdf_for_GMM(OneGdf, Ngauss, Garray)
  struct GAUSS3D *OneGdf; /* OneGdf for GMM (to be calculated) */ 
  int Ngauss;             /* number of gdf */
  struct GAUSS3D *Garray; /* GMM [0..Ngauss-1] */

{
  int i,j,g;
  double WeightAll;
/*

  Sigma = \sum[Wi*(Sigma_i + m_i*m_iT)] - m*mT


*/ 
  for (i=0;i<3;++i){
    OneGdf->M[i] = 0.0;
    for (j=0;j<3;++j){
      OneGdf->CovM[i][j] = 0.0;
    }
  }

  WeightAll = 0.0;
  for (g=0;g<Ngauss;++g){
    WeightAll += Garray[g].Weight; 
    for (i=0;i<3;++i){
      OneGdf->M[i] += Garray[g].Weight * Garray[g].M[i];
      for (j=0;j<3;++j){
        OneGdf->CovM[i][j] +=
         Garray[g].Weight * (Garray[g].CovM[i][j] + Garray[g].M[i]*Garray[g].M[j]);
      }
    }
  }
  for (i=0;i<3;++i){
    OneGdf->M[i] =  Garray[g].M[i]/WeightAll;
  }

  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      OneGdf->CovM[i][j] = OneGdf->CovM[i][j]/WeightAll - OneGdf->M[i]*OneGdf->M[j];
    }
  }

  Cal_Inverse_Matrix3D_by_Cramer_Rule(OneGdf->iCovM,OneGdf->CovM,&(OneGdf->det));
  OneGdf->det = Determinant_Matrix3D(OneGdf->CovM);
  OneGdf->Cons = 1.0/(pow(2*M_PI,1.5)*sqrt(OneGdf->det));


}/* end of Set_One_Gdf_for_GMM() */




void Setup_PCvar_and_PCaxis(CovM,PCvar,PCaxis)
  double CovM[3][3];   /* covariance matrix (input) */
  double PCvar[3];     /* variance for each PC axis (to be calculated) */
  double PCaxis[3][3]; /* PC axis [axis_num][x,y,z] (to be calculated) */
{
 int i,j,index[3];
 double A[3][3],U[3][3],PC0xPC1[3];
 double y[3];

 for (i=0;i<3;++i){
   y[i] = 0.0;
   for (j=0;j<3;++j){ A[i][j] = CovM[i][j];}
 }

 Jacobi_Wilkinson3(A,U);

 Sort_Eigen_Value3(A,index,'D'); /* sort by decreasing order of eigen values */



 for (i=0;i<3;++i){
   PCvar[i]  = A[index[i]][index[i]];
   for (j=0;j<3;++j){
     PCaxis[i][j] = U[j][index[i]];
   }
 }

  /* Check handedness, force the right-handed system */
  Cross_Prod3D(PC0xPC1,PCaxis[0],PCaxis[1]);
  if (Dot_Prod3D(PC0xPC1,PCaxis[2])<0.0){
    for (i=0;i<3;++i){ PCaxis[2][i] = -PCaxis[2][i];}
  }

} /* end of Setup_PCvar_and_PCaxis() */



double Set_QformThre_from_UpperCumuChiSquare(chi2)
  double chi2;
{
 double chi2_tab[50],uppcumu_tab[50];
 double eps, uppcumu_interpolate;
 int    Ntab,t; 
 /*
  printf("#Set_QformThre_from_UpperCumuChiSquare(chi2:%lf %e)\n",chi2,chi2);
 */ 
  eps = 0.0000000001; 
/*
   Chi2 vs UpperCumurativeDistributionFunction
   are taken from Excel Function "CHI2.INV.RT(value, degree_of_freedom)"
   with degree of freedom = 3;
 */
  chi2_tab[ 0] = 0.000001; uppcumu_tab[ 0] = 30.66484971;
  chi2_tab[ 1] = 0.00001;  uppcumu_tab[ 1] = 25.90174975;
  chi2_tab[ 2] = 0.0001;   uppcumu_tab[ 2] = 21.10751347;
  chi2_tab[ 3] = 0.001;    uppcumu_tab[ 3] = 16.2662362;
  chi2_tab[ 4] = 0.005;    uppcumu_tab[ 4] = 12.83815647;
  chi2_tab[ 5] = 0.01;     uppcumu_tab[ 5] = 11.34486673;
  chi2_tab[ 6] = 0.05;     uppcumu_tab[ 6] = 7.814727903;
  chi2_tab[ 7] = 0.1;      uppcumu_tab[ 7] = 6.251388631;
  chi2_tab[ 8] = 0.2;      uppcumu_tab[ 8] = 4.641627676;
  chi2_tab[ 9] = 0.3;      uppcumu_tab[ 9] = 3.664870783;
  chi2_tab[10] = 0.4;      uppcumu_tab[10] = 2.946166073;
  chi2_tab[11] = 0.5;      uppcumu_tab[11] = 2.365973884;
  chi2_tab[12] = 0.6;      uppcumu_tab[12] = 1.869168403;
  chi2_tab[13] = 0.7;      uppcumu_tab[13] = 1.423652243;
  chi2_tab[14] = 0.8;      uppcumu_tab[14] = 1.005174013;
  chi2_tab[15] = 0.9;      uppcumu_tab[15] = 0.584374374;
  chi2_tab[16] = 0.99;     uppcumu_tab[16] = 0.114831802;
  chi2_tab[17] = 0.999;    uppcumu_tab[17] = 0.024297586;
  chi2_tab[18] = 0.9999;   uppcumu_tab[18] = 0.005214832;
  chi2_tab[19] = 0.99999;  uppcumu_tab[19] = 0.001122583;
  Ntab = 20; 

  t = 0;
  while (t<Ntab){
    if (fabs(chi2 - chi2_tab[t])<eps){ 
      return(uppcumu_tab[t]);
     }
    t += 1;
  }

  t = 0;
  while (t<(Ntab-1)){
    if ((chi2_tab[t]<=chi2)&&(chi2<chi2_tab[t+1])){ 
      uppcumu_interpolate = 
       (uppcumu_tab[t+1] - uppcumu_tab[t])/(chi2_tab[t+1] - chi2_tab[t])
        * (chi2 - chi2_tab[t]) + uppcumu_tab[t];
    return(uppcumu_interpolate);  
    }
    t += 1;
  }

  printf("#ERROR(Set_QformThre_from_UpperCumuChiSquare()):Qform for chi2 %lf is not prepared.\n",chi2);
   exit(1);
} /* end of Set_QformThre_from_UpperCumuChiSquare() */







void Get_RotAxis_and_RotAngle_from_RotMatrix(axis,angle,Rmat)
  double axis[3]; /* rotational axis (to be calculated) */
  double *angle;  /* rotational angle [radian](to be calculated) */
  double Rmat[3][3];  /* rotational matrix (input) */
{
  double cos_t,sin_t;
  double sq_axis[3],max_sq_axis;
  int i,max_i,alpha,beta,gamma;
 
  cos_t = (Rmat[0][0]+Rmat[1][1]+Rmat[2][2]-1.0)/2.0;
  max_i = -1; max_sq_axis = 0.0;
  for (i=0;i<3;++i){
    sq_axis[i] = (Rmat[i][i]-cos_t)/(1.0 - cos_t);
    if ((i==0) || (sq_axis[i]>max_sq_axis)){
      max_i = i;
      max_sq_axis = sq_axis[i];
    }
  }  

  alpha = max_i;
  beta  = (alpha+1)%3;
  gamma = (alpha+2)%3;
  axis[alpha] = sqrt(sq_axis[alpha]);
  sin_t = (Rmat[gamma][beta] - Rmat[beta][gamma])/(2.0*axis[alpha]);
  *angle = acos(cos_t);
  if (sin_t<0.0){
    *angle = -(*angle);
  }
  axis[beta]  = (Rmat[alpha][gamma] - Rmat[gamma][alpha])/(2.0*sin_t);
  axis[gamma] = (Rmat[beta][alpha]  - Rmat[alpha][beta] )/(2.0*sin_t);

} /* end of Get_RotAxis_and_RotAngle_from_RotMatrix() */



void value_to_RGB(R, G, B, x, type)
  float *R,*G,*B;
  float x; /* input value 0<=x<=1 */
  char *type;  /* "RWB", "BGR", "GWB", "GWC", "BWR", "RGB", "BWG", "CWG", */
{
  if ((strcmp(type,"RWB")==0)||(strcmp(type,"BGR")==0)||
      (strcmp(type,"GWB")==0)|| (strcmp(type,"GWC")==0)){
    x = 1.0 - x;
  }

  if ((strcmp(type,"BWR")==0)||(strcmp(type,"RWB")==0)){
    if (x<=0.5){
      *R = 2.0*x;
      *G = 2.0*x;
      *B = 1.0;
    } 
    else if (x<=1.0){
      *R = 1.0;
      *G = 1.0-2.0*(x-0.5);
      *B = 1.0-2.0*(x-0.5);
    }
  }
  if ((strcmp(type,"BWG")==0)||(strcmp(type,"GWB")==0)){
    if (x<=0.5){
      *R = 2.0*x;
      *G = 2.0*x;
      *B = 1.0;
    }
    else if (x<=1.0){
      *R = 1.0-2.0*(x-0.5);
      *G = 1.0;
      *B = 1.0-2.0*(x-0.5);
    }
  }
  if ((strcmp(type,"CWR")==0)||(strcmp(type,"RWC")==0)){
    if (x<=0.5){
      *R = 2.0*x;
      *G = 1.0;
      *B = 1.0;
    }
    else if (x<=1.0){
      *R = 1.0;
      *G = 1.0-2.0*(x-0.5);
      *B = 1.0-2.0*(x-0.5);
    }
  }

  if ((strcmp(type,"RGB")==0)||(strcmp(type,"BGR")==0)){
    if (x<=0.25){
      *R = 1.0;
      *G = 4.0*x;
      *B = 0.0;
    }
    else if (x<=0.50){
      *R = 1.0 -4.0*(x-0.25);
      *G = 1.0;
      *B = 0.0;
    }
    else if (x<=0.75){
      *R = 0.0;
      *G = 1.0;
      *B = 4.0*(x-0.5);
    }
    else if (x<=1.0){
      *R = 0.0;
      *G = 1.0-4.0*(x-0.75);
      *B = 1.0;
    }
  }
} /* end of value_to_RGB() */



void find_min_max_Weight(minW,maxW,Ngauss,G)
  double *minW;
  double *maxW;
  int Ngauss;
  struct GAUSS3D *G; /* array .[0..Ngauss-1] */
{
  int g;

  *minW = 0.0; 
  *maxW = 0.0; 
  for (g=0;g<Ngauss;++g){
    if ((g==0)||(G[g].Weight< *minW)){
      *minW = G[g].Weight;
    }
    if ((g==0)||(G[g].Weight> *maxW)){
      *maxW = G[g].Weight;
    }
  }
} /* end of find_min_max_Weight() */

