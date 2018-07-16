/*

 <3DmapEM.c>
 for EM fitting for Gaussian Mixture Model for 3D Voxel Map Data

 LastModified: 2015/08/03 

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "globalvar.h"
#include "pdbstruct.h"
#include "PdbIO.h"
#include "Radius.h"
#include "gauss.h"
#include "Matrix3D.h"
#include "Voxel.h"
#include "MCubeFunc.h"
#include "PointEM.h"
#include "GaussIO.h"
#include "GridMap.h"

/*** Functions (GLOBAL) ***/
void Set_Less_Than_Threshold_Voxels_to_Zero();
void Cal_Mean_and_CovarMat_from_Voxel();
int  EM_optimize_GaussMix_for_3D_Map_SmallMemory();
int  EM_optimize_GaussMix_for_3D_Map();
double Corr_Coeff_Bwn_Two_Voxels();
double Corr_Coeff_Bwn_Voxel_and_GMM();
double Log_Likelihood_Of_GaussMix_For_3Dmap();
void Set_GaussMix_Density_To_3Dmap();

/*** Functions (LOCAL) ***/
static void change_gdf_to_unit_CovM_center_M_zero_Weight();

/*****************/
/*** FUNCTIONS ***/
/*****************/



void Set_Less_Than_Threshold_Voxels_to_Zero(vox,threshold)
 struct VOXEL *vox;
 float threshold;  /* normally set to 0.0 */
{
 int x,y,z,Nback,Nfore,Nvoxel;
 double Vvoxel;

 /*
 Modified By T.K. (2015/08/03) 
 */
 Nfore = 0;
 Nback = 0;
 Nvoxel = vox->N[0]*vox->N[1]*vox->N[2];
 Vvoxel = vox->grid_width *  vox->grid_width * vox->grid_width;

 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){ 
       if (vox->dat[x][y][z] < threshold){ 
        vox->dat[x][y][z] = 0.0; 
        Nback += 1;
       }
       else{
         Nfore += 1;
       }
     }
   } 
 }
 printf("#Set_Less_Than_Threshold_Voxels_to_Zero(threshold %f Nzero %d/%d %.2lf %%)\n",threshold,
   Nback,Nvoxel,100.0*(double)Nback/(double)Nvoxel);
 printf("THRESHOLD         %lf\n",threshold); 
 printf("NVOXEL_FOREGROUND %d\n",Nfore); 
 printf("NVOXEL_BACKGROUND %d\n",Nback); 
 printf("VOLUME_FOREGROUND %lf\n",Vvoxel * Nfore); 
 printf("VOLUME_BACKGROUND %lf\n",Vvoxel * Nback); 
  
 
} /* end of Set_Less_Than_Threshold_Voxels_to_Zero() */


                                                                                                               
void Cal_Mean_and_CovarMat_from_Voxel(vox,M,Cov,Min,Max)
 struct VOXEL *vox;
 double  M[3];
 double Cov[3][3];
 double Min[3],Max[3];
{
 int x,y,z,i,j;
 double sumDensity,X[3],SD,mag;
 
 /** Cal Mean Vector **/
 sumDensity = 0.0;
 M[0] = M[1] = M[2] = 0.0;
 for (x=0;x<vox->N[0];++x){  
   X[0] = vox->OrigPos[0] + vox->grid_width *x;
   for (y=0;y<vox->N[1];++y){ 
     X[1] = vox->OrigPos[1] + vox->grid_width *y;
     for (z=0;z<vox->N[2];++z) { 
       if (vox->dat[x][y][z] > 0.0){
         X[2] = vox->OrigPos[2] + vox->grid_width *z;
         sumDensity  += vox->dat[x][y][z];
         for (i=0;i<3;++i){
            M[i] += vox->dat[x][y][z] * X[i];
         }
       }
     }
   }
 }

 printf("#sumDensity %lf\n",sumDensity);

 for (i=0;i<3;++i){ M[i] /= sumDensity;}

 /** Cal Covariance Matrix **/

  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      Cov[i][j] = 0.0;
    }
  }

  for (x=0;x<vox->N[0];++x){  
    X[0] = vox->OrigPos[0] + vox->grid_width *x;
    for (y=0;y<vox->N[1];++y){ 
      X[1] = vox->OrigPos[1] + vox->grid_width *y;
      for (z=0;z<vox->N[2];++z){ 
        if (vox->dat[x][y][z] > 0.0){
          X[2] = vox->OrigPos[2] + vox->grid_width *z;
          for (i=0;i<3;++i){ 
            for (j=0;j<3;++j){ 
              Cov[i][j] += vox->dat[x][y][z]*(X[i]-M[i])*(X[j]-M[j]);
            }
          }
        }
      }
    }
  }


  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      Cov[i][j] /= sumDensity;
    }
  }


  mag = 2.0;
  for (i=0;i<3;++i){
    SD = sqrt(Cov[i][i]);
    Min[i] = M[i] - mag * SD;
    Max[i] = M[i] + mag * SD;
  }

} /* end of Cal_Mean_and_CovarMat_from_Voxel() */





int EM_optimize_GaussMix_for_3D_Map_SmallMemory(Nrepeat,Ngauss,Garray,vox,logLike_final)
 int Nrepeat;
 int Ngauss;
 struct GAUSS3D *Garray;
 struct VOXEL  *vox;
 double *logLike_final;
{
 int r,g,h,i,j,x,y,z;
 struct VOXEL Hvox;
 double sumG,val,Hg_sum_a,logLike,sumDensity;
 float X[3];

 printf("#EM_optimize_GaussMix_for_3D_Map_SmallMemory(Nrepeat,Ngauss,Garray,vox,logLike_final)\n");

 /*** [1] Malloc matrix Hvox ***/
 Malloc_Voxel(&Hvox,vox->N[0],vox->N[1],vox->N[2]);
 sumDensity = 0.0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){ sumDensity  += vox->dat[x][y][z];}
   }
 }
 /*** [3] Iteration(r) for EM  ***/
 
 r = 0;
 while (r<Nrepeat){
  for (g=0;g<Ngauss;++g){ 
    /** (1) Calculate Hvox **/ 
    for (x=0;x<vox->N[0];++x){  
      X[0] = vox->OrigPos[0] + vox->grid_width *x;
      for (y=0;y<vox->N[1];++y){ 
        X[1] = vox->OrigPos[1] + vox->grid_width *y;
        for (z=0;z<vox->N[2];++z){ 
          X[2] = vox->OrigPos[2] + vox->grid_width *z;
          sumG = 0.0;
          for (h=0;h<Ngauss;++h){ 
            val = Garray[h].Weight * Gaussian3D_Value_Table(&(Garray[h]),X);
            sumG += val;
            if (h==g){Hvox.dat[x][y][z] = val;}
          } 
          Hvox.dat[x][y][z] /= sumG;
       } 
     } 
   } 

   /** (2) Calculate Weight[g] **/ 
   Hg_sum_a =  0.0;
   for (x=0;x<vox->N[0];++x){
     for (y=0;y<vox->N[1];++y){
       for (z=0;z<vox->N[2];++z){Hg_sum_a += vox->dat[x][y][z] * Hvox.dat[x][y][z];}
     }
   }

   Garray[g].Weight = Hg_sum_a/sumDensity;

   for (i=0;i<3;++i) Garray[g].M[i] = 0.0;
    
    /* (3) Calculate Mean M[] */
    for (x=0;x<vox->N[0];++x){
      X[0] = vox->OrigPos[0] + vox->grid_width *x;
      for (y=0;y<vox->N[1];++y){ 
        X[1] = vox->OrigPos[1] + vox->grid_width *y;
        for (z=0;z<vox->N[2];++z){ 
          X[2] = vox->OrigPos[2] + vox->grid_width *z;
          for (i=0;i<3;++i) { Garray[g].M[i] += vox->dat[x][y][z] * Hvox.dat[x][y][z] * X[i];}
        }
      }
    } 

   for (i=0;i<3;++i){Garray[g].M[i]  /=  Hg_sum_a;}

    /* (4) Caluclate Covariance Matrix CovM[] */
    for (i=0;i<3;++i){ 
     for (j=0;j<3;++j){Garray[g].CovM[i][j] = 0.0;}
    }

   for (x=0;x<vox->N[0];++x){
     X[0] = vox->OrigPos[0] + vox->grid_width *x;
     for (y=0;y<vox->N[1];++y){ 
       X[1] = vox->OrigPos[1] + vox->grid_width *y;
       for (z=0;z<vox->N[2];++z){ 
         X[2] = vox->OrigPos[2] + vox->grid_width *z;
         for (i=0;i<3;++i){ 
           for (j=0;j<3;++j){ 
             Garray[g].CovM[i][j] +=  
             vox->dat[x][y][z] * Hvox.dat[x][y][z]  * (X[i]- Garray[g].M[i]) * (X[j] - Garray[g].M[j]);
           }
         }
       }
      }
    }

   for (i=0;i<3;++i){ 
     for (j=0;j<3;++j){ Garray[g].CovM[i][j] /= Hg_sum_a;}
   }
   
  /* (5) Set up iCovM */
   Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));
   Garray[g].Cons = 1.0/(two_pi_32*sqrt(Garray[g].det));

  } /* g */

  logLike = Log_Likelihood_Of_GaussMix_For_3Dmap(Ngauss,Garray,vox);

  printf("Step %d logLike %lf\n",r,logLike);
  *logLike_final = logLike;
  ++r;
 } /* r */

 Free_Voxel(&Hvox,vox->N[0],vox->N[1],vox->N[2]);
 return(1);
} /* end of EM_optimize_GaussMix_for_3D_Map_SmallMemory() */








int EM_optimize_GaussMix_for_3D_Map(Nrepeat,Ngauss,Garray,vox,logLike_final,ConvThre)
 int Nrepeat;
 int Ngauss;
 struct GAUSS3D *Garray;
 struct VOXEL  *vox;
 double *logLike_final;
 double ConvThre;  /* Threshold value for logLike convergence */
{
 struct FLOAT4DMAP Ppost; /* prob(g/x,y,z) [x][y][z][g] */
 /* struct DOUBLE4DMAP Ppost; */ /* prob(g/x,y,z) [x][y][z][g] */
 float *X,*Y,*Z;          /* X,Y,Z coordinates [x], [y], [z] */
 double sumDensity;             /* sum of density vox->dat[][][] */
 double Sxyz_vox_hpost;    /* sum of density x hpos */
 int r,g,i,j,x,y,z;
 double sumWG,logLike,logLikePre,DiffLogLike,minVar;
 float Pos[3];

 minVar = (0.5*vox->grid_width)*(0.5*vox->grid_width);
 
 printf("#EM_optimize_GaussMix_for_3D_Map(N %d %d %d Ngauss %d Nrep %d ConvThre %lf grid_width %lf minVar %lf)\n",
  vox->N[0],vox->N[1],vox->N[2],Ngauss,Nrepeat,ConvThre,vox->grid_width,minVar);

 if (PAR.ologfile[0]!='\0'){
   fprintf(PAR.fp_olog,"#EM_optimize_GaussMix_for_3D_Map(N %d %d %d Ngauss %d Nrep %d ConvThre %lf)\n",
    vox->N[0],vox->N[1],vox->N[2],Ngauss,Nrepeat,ConvThre);
 }
  
 logLike = Log_Likelihood_Of_GaussMix_For_3Dmap(Ngauss,Garray,vox);
 logLikePre = logLike - 10000.0*ConvThre;

 if (PAR.ologfile[0]!='\0'){
   fprintf(PAR.fp_olog,"#Initial_logLike %lf %e\n",logLike,logLike);
   fprintf(PAR.fp_olog,"#[Step] [logLike]  [Diff] [Voxel_sizeX]:[Y]:[Z]\n");
 }




 /*** [1] Malloc matrix Ppost[][][][] and X[],Y[],Z[] ***/
 X = (float *)malloc(sizeof(float)*vox->N[0]);
 Y = (float *)malloc(sizeof(float)*vox->N[1]);
 Z = (float *)malloc(sizeof(float)*vox->N[2]);

 Malloc_FLOAT4DMAP(&Ppost,vox->N[0],vox->N[1],vox->N[2],Ngauss);
/* 
 Malloc_DOUBLE4DMAP(&Ppost,vox->N[0],vox->N[1],vox->N[2],Ngauss);
 */

 /*** [2] Calculate sumDensity and X[],Y[],Z[] ***/
 sumDensity = 0.0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){
       if (vox->dat[x][y][z]>0.0){
         sumDensity  += vox->dat[x][y][z];
       }
     }
   }
 }

 for (x=0;x<vox->N[0];++x){ X[x] = vox->OrigPos[0] + vox->grid_width * x;}
 for (y=0;y<vox->N[1];++y){ Y[y] = vox->OrigPos[1] + vox->grid_width * y;}
 for (z=0;z<vox->N[2];++z){ Z[z] = vox->OrigPos[2] + vox->grid_width * z;}

 /*** [3] Iteration(r) for EM  ***/
 
 r = 0;

 do {
   
   /** [E-step] Calculate Ppost **/ 
   logLike = 0.0;
   for (x=0;x<vox->N[0];++x){
     Pos[0] = X[x];
     for (y=0;y<vox->N[1];++y){
       Pos[1] = Y[y];
       for (z=0;z<vox->N[2];++z){
         Pos[2] = Z[z];
         if (vox->dat[x][y][z]>0.0){
           sumWG = 0.0;
           for (g=0;g<Ngauss;++g){  
             if (Garray[g].Weight > 0.0){ 
               Ppost.dat[x][y][z][g] = Garray[g].Weight * Gaussian3D_Value_Table(&(Garray[g]),Pos);
               sumWG += Ppost.dat[x][y][z][g];
              }
              else {
                Ppost.dat[x][y][z][g] = 0.0;
              }
           } 

           if (sumWG > 0.0){ 
             for (g=0;g<Ngauss;++g){
                Ppost.dat[x][y][z][g] /= sumWG;
             }
             logLike += vox->dat[x][y][z] * log(sumWG);
           }
         }
       }
     }
   } 
 
   DiffLogLike = logLike - logLikePre; 
   logLikePre = logLike;

   printf("#step %d %lf %lf %d:%d:%d\n",r,logLike,DiffLogLike,vox->N[0],vox->N[1],vox->N[2]);
   if (PAR.ologfile[0]!='\0'){
     fprintf(PAR.fp_olog,"%d %lf  %lf %d:%d:%d\n",r,logLike,DiffLogLike,vox->N[0],vox->N[1],vox->N[2]);
   }
 
   /** [M-step] Calculate Weight[],M[], Sig[][] by g-loop  **/ 
   for (g=0;g<Ngauss;++g){
    if (Garray[g].Weight > 0.0){
      /** (2-0) Initialize **/ 
      Sxyz_vox_hpost = 0.0; 
      for (i=0;i<3;++i){ 
        Garray[g].M[i] = 0.0;
        for (j=0;j<3;++j){Garray[g].CovM[i][j] = 0.0;}
      }
  
      /** (2-1) Summation for x,y,z to get Sxyz_vox_hpost, M[],CovM **/ 
      for (x=0;x<vox->N[0];++x){
        Pos[0] = X[x];
        for (y=0;y<vox->N[1];++y){
          Pos[1] = Y[y];
          for (z=0;z<vox->N[2];++z){ 
            Pos[2] = Z[z];
            if (vox->dat[x][y][z]>0.0){   
              Sxyz_vox_hpost  += vox->dat[x][y][z] * Ppost.dat[x][y][z][g]; 
          
              for (i=0;i<3;++i){   
                Garray[g].M[i] += vox->dat[x][y][z] * Ppost.dat[x][y][z][g] * Pos[i]; 
                for (j=i;j<3;++j){
                  Garray[g].CovM[i][j] += vox->dat[x][y][z] * Ppost.dat[x][y][z][g] * Pos[i]*Pos[j]; 
                }
              }
            }
          }
        }
      }

        Garray[g].Weight = Sxyz_vox_hpost/sumDensity;


        for (i=0;i<3;++i){Garray[g].M[i]  /=  Sxyz_vox_hpost;}
      
        for (i=0;i<3;++i){ 
          for (j=i;j<3;++j){ 
            Garray[g].CovM[i][j] /= Sxyz_vox_hpost;
            Garray[g].CovM[i][j] -= Garray[g].M[i]*Garray[g].M[j];
          } 
        }

        Garray[g].CovM[1][0] = Garray[g].CovM[0][1]; 
        Garray[g].CovM[2][0] = Garray[g].CovM[0][2]; 
        Garray[g].CovM[2][1] = Garray[g].CovM[1][2]; 

       /** Enforce Minimal Variance (grid_width/2)^2 for CovM[0][0], CovM[1][1], CovM[2][2].**/
       for (i=0;i<3;++i){
         if (Garray[g].CovM[i][i]<minVar){ Garray[g].CovM[i][i]=minVar;}
       }
/*
       if ((Garray[g].CovM[0][0]<minVar) || (Garray[g].CovM[1][1]<minVar) || (Garray[g].CovM[2][2]<minVar)){ 
         change_gdf_to_unit_CovM_center_M_zero_Weight(Garray,g,vox);
       }
*/
       /* NaN GDF check */
        if (Check_NaN_Inf_Gaussian3D(&Garray[g])==1){
          printf("#WARNING:Non-normal value is observed in %d-th GMM parameters.\n",g+1);
          show_Gaussian3D(&Garray[g]);
          printf("#Forced into unit CovM[][] with the center M[]  with the zero weight.\n");
          change_gdf_to_unit_CovM_center_M_zero_Weight(Garray,g,vox);
        } 
 
   
      /* (2-4) Set up iCovM */
      Cal_Inverse_Matrix3D_by_Cramer_Rule_Symmetric(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));
      Garray[g].Cons = 1.0/(two_pi_32*sqrt(Garray[g].det));
    } 
  } /* g */

  /*
  printf("#loglike %e\n",Log_Likelihood_Of_GaussMix_For_3Dmap(Ngauss,Garray,vox));
  */
 
  ++r;

 } while ((r<Nrepeat)&&(DiffLogLike>ConvThre));

  printf("#DiffLogLike %lf %e\n",DiffLogLike,DiffLogLike);


 *logLike_final = logLike;

 logLike = Log_Likelihood_Of_GaussMix_For_3Dmap(Ngauss,Garray,vox);
 printf("#logLike_final %lf %e\n",logLike,logLike);
 if (PAR.ologfile[0]!='\0'){
   fprintf(PAR.fp_olog,"#logLike_final %lf %e\n",logLike,logLike);
 }
 /*** [4] Free  Variables (Ppost,fixg,X,Y,Z) ***/
 Free_FLOAT4DMAP(&Ppost); 
 /* Free_DOUBLE4DMAP(&Ppost);  */
 free(X); free(Y); free(Z);
 return(1);
} /* end of EM_optimize_GaussMix_for_3D_Map() */














void change_gdf_to_unit_CovM_center_M_zero_Weight(Garray,g,vox)
  struct GAUSS3D *Garray;
  int g;
  struct VOXEL *vox;
{
  int i,j;

  Garray[g].Weight = 0.0;
  for (i=0;i<3;++i){
     Garray[g].M[i] = vox->OrigPos[i] + vox->grid_width * vox->N[i]/2;
     for (j=0;j<3;++j){ 
       Garray[g].CovM[i][j] = 0.0;
       if (i==j) {Garray[g].CovM[i][j] = 1.0;}
     }
  } 
  Cal_Inverse_Matrix3D_by_Cramer_Rule_Symmetric(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));
} /* end of change_to_unit_CovM_center_M_zero_Weight() */




double Corr_Coeff_Bwn_Two_Voxels(voxA,voxB)
 struct VOXEL  *voxA,*voxB;
{
 int x,y,z;
 double Sa,Sb,Saa,Sbb,Sab;
 double Ma,Mb,Va,Vb,Vab,CC,N;
 
 Sa = Sb = Saa = Sbb = Sab = 0.0;
 
 for (x=0;x<voxA->N[0];++x){
   for (y=0;y<voxA->N[1];++y){
     for (z=0;z<voxA->N[2];++z){
       Sa  += voxA->dat[x][y][z];
       Sb  += voxB->dat[x][y][z];
       Saa += voxA->dat[x][y][z] * voxA->dat[x][y][z];
       Sbb += voxB->dat[x][y][z] * voxB->dat[x][y][z];
      Sab += voxA->dat[x][y][z] * voxB->dat[x][y][z];
     }
   }
 }

 N = voxA->N[0] * voxA->N[1] * voxA->N[2];
 Ma  = Sa/N;
 Mb  = Sb/N;
 Va  = Saa/N - Ma*Ma;
 Vb  = Sbb/N - Mb*Mb;
 Vab = Sab/N - Ma*Mb;
 if ((Va>0.0)&&(Vb>0.0)) CC = Vab/sqrt(Va)/sqrt(Vb);
  else CC = 0.0;

 return(CC);
} /* end of Corr_Coeff_Bwn_Two_Voxels() */




double Corr_Coeff_Bwn_Voxel_and_GMM(voxA,Ngauss,Garray)
 struct VOXEL  *voxA;
 int Ngauss;
 struct GAUSS3D *Garray;
{
 int x,y,z,g;
 double Sa,Sb,Saa,Sbb,Sab,fb;
 double Ma,Mb,Va,Vb,Vab,CC,N;
 float Pos[3]; 
 Sa = Sb = Saa = Sbb = Sab = 0.0;
 
 printf("#Corr_Coeff_Voxel_and_GMM(voxA,Ngauss:%d)\n",Ngauss);

 for (x=0;x<voxA->N[0];++x){
   Pos[0] = voxA->OrigPos[0] + voxA->grid_width * x;
   for (y=0;y<voxA->N[1];++y){
     Pos[1] = voxA->OrigPos[1] + voxA->grid_width * y;
     for (z=0;z<voxA->N[2];++z){
       Sa  += voxA->dat[x][y][z];
       Pos[2] = voxA->OrigPos[2] + voxA->grid_width * z;
       fb = 0.0;
       for (g=0;g<Ngauss;++g){ 
         fb += Garray[g].Weight*Gaussian3D_Value(&(Garray[g]),Pos);
         /* printf(" %d %e weight %lf fb %e\n",g, Gaussian3D_Value(&(Garray[g]),Pos),Garray[g].Weight,fb); */
       }
/*
       printf("%d %d %d %f %f %f f %e\n",x,y,z,Pos[0],Pos[1],Pos[2],fb);
*/

       Sb  += fb;
       Saa += voxA->dat[x][y][z] * voxA->dat[x][y][z];
       Sbb += fb*fb;
      Sab  += voxA->dat[x][y][z] * fb;
     }
   }
 }

 N = voxA->N[0] * voxA->N[1] * voxA->N[2];
 Ma  = Sa/N;
 Mb  = Sb/N;
 Va  = Saa/N - Ma*Ma;
 Vb  = Sbb/N - Mb*Mb;
 Vab = Sab/N - Ma*Mb;
 if ((Va>0.0)&&(Vb>0.0)) {CC = Vab/sqrt(Va)/sqrt(Vb);} else {CC = 0.0;}

 return(CC);
} /* end of Corr_Coeff_Bwn_Voxel_and_GMM() */





double Log_Likelihood_Of_GaussMix_For_3Dmap(Ngauss, Garray,vox)
 int Ngauss;
 struct GAUSS3D *Garray;
 struct VOXEL  *vox;
{
 int x,y,z,g;
 double PrX_G,logLike;
 float X[3];

 logLike = 0.0;
 
 for (x=0;x<vox->N[0];++x){  
   X[0] = vox->OrigPos[0] + vox->grid_width *x;
   for (y=0;y<vox->N[1];++y){ 
    X[1] = vox->OrigPos[1] + vox->grid_width *y;
    for (z=0;z<vox->N[2];++z) { 
     if (vox->dat[x][y][z]>0.0){
      X[2] = vox->OrigPos[2] + vox->grid_width *z;
      PrX_G = 0.0;
      for (g=0;g<Ngauss;++g) 
       if (Garray[g].Weight>0.0) 
       PrX_G += Garray[g].Weight * Gaussian3D_Value(&(Garray[g]),X);
      if (PrX_G>0.0) logLike += vox->dat[x][y][z] * log(PrX_G); 
     }
   } /* z */
  } /* y */
 } /* x */

 return(logLike);

} /* end of Log_Likelihood_Of_Gaussians_For_3Dmap() */


void Set_GaussMix_Density_To_3Dmap(Ngauss, Garray,vox)
 int Ngauss;
 struct GAUSS3D *Garray;
 struct VOXEL  *vox;
{
 int x,y,z,g;
 double PrX_G;
 float X[3];

 /** Cal Gaussian Values **/
 for (x=0;x<vox->N[0];++x){  
   X[0] = vox->OrigPos[0] + vox->grid_width *x;
   for (y=0;y<vox->N[1];++y){ 
    X[1] = vox->OrigPos[1] + vox->grid_width *y;
    for (z=0;z<vox->N[2];++z){  
      X[2] = vox->OrigPos[2] + vox->grid_width *z;
      PrX_G = 0.0;
      for (g=0;g<Ngauss;++g){ 
       if (Garray[g].Weight>0.0) PrX_G += Garray[g].Weight * Gaussian3D_Value(&(Garray[g]),X);
     } 
      vox->dat[x][y][z] = PrX_G; 
    } 
   } 
  } 

} /* end of Set_GaussMix_Density_To_3Dmap() */
