/*

 <Gmm3DmapEM.c>
 for EM fitting for Gaussian Mixture Model for 3D Voxel Map Data

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
int  EM_optimize_GaussMix_for_Gaussian_3D_Map();
double Corr_Coeff_Bwn_VoxelGMM_and_GMM();
int  GaussMix_by_One_Voxel_One_Gdf_Assignment();
void cal_sumDensity_and_Nvox_posi_dens();



/*** Functions (LOCAL) ***/
static double Estep_of_GMM_For_Gaussian_3Dmap();
static double Likelihood_for_Observed_Isotropic_GDF();





int EM_optimize_GaussMix_for_Gaussian_3D_Map(Nrepeat,Ngauss,Garray,vox,VarType,Cww2var,resolution,logLike_final,ConvThre)
 int Nrepeat;
 int Ngauss;
 struct GAUSS3D *Garray;
 struct VOXEL  *vox;
 char   VarType;      /* G: Var = Cww2var * grid_width*grid_width R: var = (resolution/2)^2 */
 double Cww2var;      /* Constant bwn variance and grid_width. Variance = Cww2var * grid_width * grid_width */
 double resolution;   /* Var = (resolution/2)^2 */
 double *logLike_final;
 double ConvThre;  /* Threshold value for logLike convergence */
{
 /* struct FLOAT4DMAP PostProb; */ /* prob(g/x,y,z) [x][y][z][g] */
 struct DOUBLE4DMAP PostProb;  /* prob(g/x,y,z) [x][y][z][g] */
 float *X,*Y,*Z;          /* X,Y,Z coordinates [x], [y], [z] */
 double sumDensity;             /* sum of density vax->dat[][][] */
 double sum_w_post;    /* sum of density x hpos */
 int Nvoxel,r,g,i,j,x,y,z;
 double logLike,logLikePre,DiffLogLike,logLikeTemp,w_post,Variance;
 float Pos[3];

 if ((VarType=='R') && (resolution>0.0)){
   Variance = (resolution/2.0)*(resolution/2.0);
 }
 else{ 
   Variance = Cww2var * vox->grid_width * vox->grid_width;
 }

 printf("#EM_optimize_GaussMix_for_Gaussian_3D_Map(N %d %d %d Ngauss %d Nrep %d ConvThre %lf)\n", vox->N[0],vox->N[1],vox->N[2],Ngauss,Nrepeat,ConvThre);
 
 printf("#grid_width %lf VarType:%c Cww2var:%lf resolution:%lf Variance %lf)\n",vox->grid_width,VarType,Cww2var,resolution,Variance);

 if (PAR.ologfile[0]!='\0'){
   fprintf(PAR.fp_olog,"#EM_optimize_GaussMix_for_Gaussian_3D_Map(N %d %d %d Ngauss %d Nrep %d ConvThre %lf)\n",
    vox->N[0],vox->N[1],vox->N[2],Ngauss,Nrepeat,ConvThre);
 }
  

 /*** [1] Malloc matrix PostProb[][][][] and X[],Y[],Z[] ***/
 X = (float *)malloc(sizeof(float)*vox->N[0]);
 Y = (float *)malloc(sizeof(float)*vox->N[1]);
 Z = (float *)malloc(sizeof(float)*vox->N[2]);
/*
 Malloc_FLOAT4DMAP(&PostProb,vox->N[0],vox->N[1],vox->N[2],Ngauss);
*/
 Malloc_DOUBLE4DMAP(&PostProb,vox->N[0],vox->N[1],vox->N[2],Ngauss);
 /*** [2] Calculate Nvoxel,sumDensity and X[],Y[],Z[] ***/
 cal_sumDensity_and_Nvox_posi_dens(vox,&sumDensity,&Nvoxel);
 
 printf("#Nvoxel %d (%d) sumDensity %lf\n",Nvoxel,vox->N[0]*vox->N[1]*vox->N[2],sumDensity);

 for (x=0;x<vox->N[0];++x){ X[x] = vox->OrigPos[0] + vox->grid_width * x;}
 for (y=0;y<vox->N[1];++y){ Y[y] = vox->OrigPos[1] + vox->grid_width * y;}
 for (z=0;z<vox->N[2];++z){ Z[z] = vox->OrigPos[2] + vox->grid_width * z;}

 logLike = Estep_of_GMM_For_Gaussian_3Dmap(Ngauss,Garray,vox,Nvoxel,sumDensity,X,Y,Z,Variance,'-',&PostProb); 

 logLikePre = logLike - 10000.0*ConvThre;
 printf("#Initial_logLike %lf %e\n",logLike,logLike);
 
 if (PAR.ologfile[0]!='\0'){
   fprintf(PAR.fp_olog,"#Initial_logLike %lf %e\n",logLike,logLike);
   fprintf(PAR.fp_olog,"#[Step] [logLike]  [Diff] [Voxel_sizeX]:[Y]:[Z]\n");
 }

 /*** [3] Iteration(r) for EM  ***/
 
 r = 0;

 do {
   
   /** [E-step] Calculate PostProb  **/ 
   logLike = Estep_of_GMM_For_Gaussian_3Dmap(Ngauss,Garray,vox,Nvoxel,sumDensity,X,Y,Z,Variance,'U',&PostProb); 

   DiffLogLike = logLike - logLikePre; 
   logLikePre = logLike;
   logLikeTemp = Estep_of_GMM_For_Gaussian_3Dmap(Ngauss,Garray,vox,Nvoxel,sumDensity,X,Y,Z,Variance,'-',&PostProb); 

   printf("#step %d %lf %lf logLike %lf\n",r,logLike,logLikePre,logLikeTemp);
   if (PAR.ologfile[0]!='\0'){
     fprintf(PAR.fp_olog,"%d %lf  %lf %d:%d:%d\n",r,logLike,DiffLogLike,vox->N[0],vox->N[1],vox->N[2]);
   }
 
   /** [M-step] Calculate Weight[],M[], CovM[][] by 'g'-loop  **/ 
   for (g=0;g<Ngauss;++g){
      /** (2-0) Initialize **/ 
      for (i=0;i<3;++i){ 
        Garray[g].M[i] = 0.0;
        for (j=0;j<3;++j){Garray[g].CovM[i][j] = 0.0;}
      }
      sum_w_post = 0.0; 
      for (x=0;x<vox->N[0];++x){
        Pos[0] = X[x];
        for (y=0;y<vox->N[1];++y){
          Pos[1] = Y[y];
          for (z=0;z<vox->N[2];++z){ 
            Pos[2] = Z[z];
            if (vox->dat[x][y][z]>0.0){ 
              /*  
              w_post = vox->dat[x][y][z] * PostProb.dat[x][y][z][g] /sumDensity; 
              */
              w_post = vox->dat[x][y][z] * PostProb.dat[x][y][z][g];
/*
              printf("#PostProb.dat[%d/%d][%d/%d][%d/%d][%d] %lf vox->dat %e\n",x,vox->N[0],y,vox->N[1],z,vox->N[2],g,PostProb.dat[x][y][z][g],vox->dat[x][y][z]); 
*/
              sum_w_post  += w_post;
              for (i=0;i<3;++i){   
                Garray[g].M[i] += w_post * Pos[i]; 
                for (j=i;j<3;++j){
                  if (i==j){
                    Garray[g].CovM[i][j] += w_post * (Variance + Pos[i]*Pos[j]);
                  }
                  else{
                    Garray[g].CovM[i][j] += w_post * Pos[i]*Pos[j]; 
                  }
                }
              }
            }
          }
        }
      }

      Garray[g].Weight = sum_w_post/sumDensity;
 /*
      printf("#Garray[%d].Weight %lf sum_w_post %lf sumDensity %lf \n",g,Garray[g].Weight,sum_w_post,sumDensity); 
 */
      for (i=0;i<3;++i){Garray[g].M[i]  /=  sum_w_post;}
      
      for (i=0;i<3;++i){ 
        for (j=i;j<3;++j){ 
          Garray[g].CovM[i][j] /= sum_w_post;
          Garray[g].CovM[i][j] -=  Garray[g].M[i]*Garray[g].M[j];
        } 
      }

      Garray[g].CovM[1][0] = Garray[g].CovM[0][1]; 
      Garray[g].CovM[2][0] = Garray[g].CovM[0][2]; 
      Garray[g].CovM[2][1] = Garray[g].CovM[1][2]; 

      /* Set up iCovM */
      Cal_Inverse_Matrix3D_by_Cramer_Rule_Symmetric(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));
      Garray[g].Cons = 1.0/(two_pi_32*sqrt(Garray[g].det));

  } /* g */

  ++r;

 } while ((r<Nrepeat)&&(DiffLogLike>ConvThre));

 printf("#DiffLogLike %lf %e\n",DiffLogLike,DiffLogLike);

 *logLike_final = logLike;

 printf("#Converged after %d steps (logLike_fin:%e)\n",r,*logLike_final);


 logLike = Estep_of_GMM_For_Gaussian_3Dmap(Ngauss,Garray,vox,Nvoxel,sumDensity,X,Y,Z,Variance,'-',&PostProb); 
 printf("#logLike_final %lf %e\n",logLike,logLike);
 if (PAR.ologfile[0]!='\0'){
   fprintf(PAR.fp_olog,"#logLike_final %lf %e\n",logLike,logLike);
 }

 /*** [4] Free  Variables (PostProb,fixg,X,Y,Z) ***/
 /* Free_FLOAT4DMAP(&PostProb);  */
 Free_DOUBLE4DMAP(&PostProb);  
 free(X); free(Y); free(Z);
 return(1);

} /* end of EM_optimize_GaussMix_for_Gaussian_3D_Map() */




double Estep_of_GMM_For_Gaussian_3Dmap(Ngauss, Garray,vox,Nvoxel,sumDensity,X,Y,Z,Variance,UpdatePostProb,PostProb)
 int Ngauss;
 struct GAUSS3D *Garray;
 struct VOXEL  *vox;
 int    Nvoxel;
 double sumDensity;
 float  *X,*Y,*Z;      /* X,Y,Z coordinates [x], [y], [z] */
 double Variance; 
 char   UpdatePostProb;       /* 'U':pdate Post Prob */
 /* struct FLOAT4DMAP *PostProb; */ /* posterior probability for GDF g for (x,y,z) */
 struct DOUBLE4DMAP *PostProb; /* posterior probability for GDF g for (x,y,z) */
{
 int x,y,z,g;
 double sum_postprob,postprob,logLike;
 float Pos[3];

 logLike = 0.0;

 for (x=0;x<vox->N[0];++x){
   Pos[0] = X[x];
   for (y=0;y<vox->N[1];++y){
     Pos[1] = Y[y];
     for (z=0;z<vox->N[2];++z) {
       Pos[2] = Z[z];
       if (vox->dat[x][y][z]>0.0){
          sum_postprob = 0.0;
          postprob = 0.0;
          for (g=0;g<Ngauss;++g){
            if (Garray[g].Weight>0.0){
         /* 
              wLike = Garray[g].Weight * Likelihood_for_Observed_Isotropic_GDF(&(Garray[g]),Pos,Sigma); 
              postprob = pow(wLike,Nvoxel*vox->dat[x][y][z]/sumDensity);
         */
              postprob  = Garray[g].Weight * Likelihood_for_Observed_Isotropic_GDF(&(Garray[g]),Pos,Variance); 
              if (UpdatePostProb=='U'){ PostProb->dat[x][y][z][g] = postprob;}
              sum_postprob += postprob;
            }
            else {
              if (UpdatePostProb=='U'){ PostProb->dat[x][y][z][g] = 0.0;}
            }
          }

          if (sum_postprob>0.0){ 
             if (UpdatePostProb=='U'){ 
               for (g=0;g<Ngauss;++g){
                 PostProb->dat[x][y][z][g] /= sum_postprob;
                 if (PostProb->dat[x][y][z][g]>1.0){
                   PostProb->dat[x][y][z][g] = 1.0;
/*
                   printf("#woops PostProb[%d][%d][%d][%d] %e sum_postprob %e\n",x,y,z,g,PostProb->dat[x][y][z][g],sum_postprob);
 */
                 }
/*
                if (PostProb->dat[x][y][z][g]<0.99){
                   printf("#woops PostProb[%d][%d][%d][%d] %e sum_postprob %e\n",x,y,z,g,PostProb->dat[x][y][z][g],sum_postprob);
                }
*/
               }
             }
            logLike += log(sum_postprob);
          }
          else{
             if (UpdatePostProb=='U'){ 
               for (g=0;g<Ngauss;++g){
                 PostProb->dat[x][y][z][g] = 1.0/(double)Ngauss;
               }
             }
          }
        }
     }  
   } 
 } 

 return(logLike);

} /* end of Estep_Of_GMM_For_Gaussian_3Dmap() */











double Likelihood_for_Observed_Isotropic_GDF(G,floatX,variance)
 struct GAUSS3D *G;
 float  floatX[3];  /* Center for the observed GDF */
 double variance;   /* variance for the observed gdf */
{
 /*
   L = 1/(2pi)^(3/2)/|CovM|^(1/2) *
       exp[ -1/2*trace(iCovM*variance*I) - 1/2 (X-M)CovM(X-M)^T ]
     = 1/(2pi)^(3/2)/|CovM|^(1/2) *
       exp[ -1/2*variance*trace(iCovM) - 1/2 (X-M)CovM(X-M)^T ]
     = 1/(2pi)^(3/2)/|CovM|^(1/2) *
       exp[ -1/2*variance*(iCovM[0][0] + iCovM[1][1] + iCovM[2][2]) - 1/2 (X-M)CovM(X-M)^T ]
 */
 double D[3],xSx,val,trace;

 D[0] = floatX[0] - G->M[0];
 D[1] = floatX[1] - G->M[1];
 D[2] = floatX[2] - G->M[2];

 /*
     xSx += D[0]*(G->iCovM[0][0]*D[0] + G->iCovM[0][1]*D[1] + G->iCovM[0][2]*D[2]);
     xSx += D[1]*(G->iCovM[1][0]*D[0] + G->iCovM[1][1]*D[1] + G->iCovM[1][2]*D[2]);
     xSx += D[2]*(G->iCovM[2][0]*D[0] + G->iCovM[2][1]*D[1] + G->iCovM[2][2]*D[2]);
 */
 xSx = 0.0;
 xSx += G->iCovM[0][0]*D[0]*D[0] + G->iCovM[1][1]*D[1]*D[1] + G->iCovM[2][2]*D[2]*D[2];
 xSx +=(G->iCovM[0][1]*D[0]*D[1] + G->iCovM[0][2]*D[0]*D[2] + G->iCovM[1][2]*D[1]*D[2]) * 2.0;

 trace = G->iCovM[0][0] + G->iCovM[1][1] + G->iCovM[2][2];
 val = G->Cons * exp(-0.5*variance*trace-0.5*xSx);

 return(val);
}/* end of Likelihood_for_Observed_Isotropic_GDF() */








double Corr_Coeff_Bwn_VoxelGMM_and_GMM(vox,Ngauss,Garray,VarType,Cww2var,resolution)
 struct VOXEL  *vox;
 int Ngauss;
 struct GAUSS3D *Garray;
 char   VarType;      /* G: Var = Cww2var * grid_width*grid_width R: var = (resolution/2)^2 */
 double Cww2var;      /* Constant bwn variance and grid_width. Variance = Cww2var * grid_width * grid_width */
 double resolution;   /* Var = (resolution/2)^2 */
{
 int i,j,x,y,z,p,q,r,g,h,xx,yy,zz;
 double CC,Variance,varGG,varVV,varGV,sumDensity,Wv,Wn; 
 struct GAUSS3D Gv,Gn;
 int Nvoxel,Ngrid_near;
 
 Ngrid_near = 3;

 if ((VarType=='R') && (resolution>0.0)){
   Variance = (resolution/2.0)*(resolution/2.0);
 }
 else{ 
   Variance = Cww2var * vox->grid_width * vox->grid_width;
 }

 printf("#Corr_Coeff_Bwn_VoxelGMM_and_GMM(vox,Ngauss:%d VarType:%c Cww2var:%lf resolution:%lf Variance:%lf)\n",Ngauss,VarType,Cww2var,resolution,Variance);

 /*  [1] calculate varGG */
 varGG = 0.0;
 for (g=0;g<Ngauss;++g){
   for (h=0;h<Ngauss;++h){
     varGG += Garray[g].Weight * Garray[h].Weight *
               Overlap_Integral_Bwn_Two_GAUSS3Ds(&(Garray[g]),&(Garray[h]));
   }
 }

 /*  [2] Set up GMM for one-grid */

 for (i=0;i<3;++i){
   Gv.CovM[i][i] = Variance;
   for (j=0;j<3;++j){
     if (i!=j){Gv.CovM[i][j] = 0.0;}
   } 
 }

 Cal_Inverse_Matrix3D_by_Cramer_Rule(Gv.iCovM, Gv.CovM, &(Gv.det));
 Gv.Cons = 1.0/(pow(2*M_PI,1.5)*sqrt(Gv.det));
 Copy_GAUSS3D(&Gn,&Gv);

 /*  [3] calculate sumDensity */
 cal_sumDensity_and_Nvox_posi_dens(vox,&sumDensity,&Nvoxel);

 /*  [3] calculate varVG */
 varGV = 0.0;
 for (x=0;x<vox->N[0];++x){
   Gv.M[0] = vox->OrigPos[0] + vox->grid_width * x;
   for (y=0;y<vox->N[1];++y){
     Gv.M[1] = vox->OrigPos[1] + vox->grid_width * y;
     for (z=0;z<vox->N[2];++z){
       if (vox->dat[x][y][z]>0.0){
         Gv.M[2] = vox->OrigPos[2] + vox->grid_width * z;
         Wv = vox->dat[x][y][z]/sumDensity;
         for (g=0;g<Ngauss;++g){
           varGV += Garray[g].Weight * Wv * Overlap_Integral_Bwn_Two_GAUSS3Ds(&(Garray[g]),&Gv);
         }
       }
     }
   }
 }

 /*  [3] calculate varVV */
 varVV = 0.0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){
       if (vox->dat[x][y][z]>0.0){
         Wv = vox->dat[x][y][z]/sumDensity;
         Gv.M[0] = vox->OrigPos[0] + vox->grid_width * x;
         Gv.M[1] = vox->OrigPos[1] + vox->grid_width * y;
         Gv.M[2] = vox->OrigPos[2] + vox->grid_width * z;
         for (p=-Ngrid_near;p<=Ngrid_near;++p){
           xx = x + p;
           for (q=-Ngrid_near;q<=Ngrid_near;++q){
             yy = y + q;
             for (r=-Ngrid_near;r<=Ngrid_near;++r){
               zz = z + r;
               if ((xx>=0)&&(xx<vox->N[0])&&(yy>=0)&&(yy<vox->N[1])&&(zz>=0)&&(zz<vox->N[2])){
                 if (vox->dat[xx][yy][zz]>0.0){
                    Gn.M[0] = vox->OrigPos[0] + vox->grid_width * xx;
                    Gn.M[1] = vox->OrigPos[1] + vox->grid_width * yy;
                    Gn.M[2] = vox->OrigPos[2] + vox->grid_width * zz;
                    Wn = vox->dat[xx][yy][zz]/sumDensity;
                    varVV += Wv*Wn * Overlap_Integral_Bwn_Two_GAUSS3Ds(&Gv,&Gn);
                 }
               }
             }
           }
         }
       } 
     }
   }
 }

/*
 varVV = 0.0;
 NNvoxel = 0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){
       Gv.M[0] = vox->OrigPos[0] + vox->grid_width * x;
       Gv.M[1] = vox->OrigPos[1] + vox->grid_width * y;
       Gv.M[2] = vox->OrigPos[2] + vox->grid_width * z;
       Wv = vox->dat[x][y][z]/sumDensity;
       if (vox->dat[x][y][z]>0.0){
         for (xx=0;xx<vox->N[0];++xx){
           for (yy=0;yy<vox->N[1];++yy){
             for (zz=0;zz<vox->N[2];++zz){
               if (vox->dat[xx][yy][zz]>0.0){
                  Gn.M[0] = vox->OrigPos[0] + vox->grid_width * xx;
                  Gn.M[1] = vox->OrigPos[1] + vox->grid_width * yy;
                  Gn.M[2] = vox->OrigPos[2] + vox->grid_width * zz;
                  Wn = vox->dat[xx][yy][zz]/sumDensity;
                  varVV += Wv*Wn * Overlap_Integral_Bwn_Two_GAUSS3Ds(&Gv,&Gn);
               }
             }
           }
         }
       } 
     }
   }
 }
*/

 printf("#sumDensity %lf varGG %e varVV %e varGV %e Nvoxel %d\n",sumDensity,varGG,varVV,varGV,Nvoxel);

 CC = 0.0;
 if ((varGV>0.0) && (varGG>0.0)){ CC = varGV/sqrt(varGG*varVV); }
 return(CC);

} /* end of Corr_Coeff_Bwn_VoxelGMM_and_GMM() */




void cal_sumDensity_and_Nvox_posi_dens(vox,sumDensity,Nvox_posi_dens)
 struct VOXEL  *vox;
 double *sumDensity;
 int *Nvox_posi_dens; /* Number of voxels with positive density */
{
 int x,y,z;

 *Nvox_posi_dens = 0;
 *sumDensity = 0.0;

 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){
       if (vox->dat[x][y][z]>0.0){
         *Nvox_posi_dens += 1;
         *sumDensity  += vox->dat[x][y][z];
       }
     }
   }
 }

} /* end of cal_sumDensity_and_Nvox_posi_dens() */



int GaussMix_by_One_Voxel_One_Gdf_Assignment(Ngauss,Garray,vox,VarType,Cww2var,resolution)
 int    Ngauss;
 struct GAUSS3D *Garray;
 struct VOXEL  *vox;
 char   VarType;      /* G: Var = Cww2var * grid_width*grid_width R: var = (resolution/2)^2 */
 double Cww2var;      /* Constant bwn variance and grid_width. Variance = Cww2var * grid_width * grid_width */
 double resolution;   /* Var = (resolution/2)^2 */
{
  double sumDensity,Variance;
  int  n,i,j,x,y,z,Nvox_posi_dens;
  struct GAUSS3D Gv;  

  printf("#GaussMix_by_One_Voxel_One_Gdf_Assignment(Ngauss:%d VarType:%c Cww2var:%lf resolution:%lf)\n",Ngauss,VarType,Cww2var,resolution);
  /** Initialization **/
  if ((VarType=='R') && (resolution>0.0)){
    Variance = (resolution/2.0)*(resolution/2.0);
  }
  else{ 
    Variance = Cww2var * vox->grid_width * vox->grid_width;
  }

  for (i=0;i<3;++i){
    Gv.CovM[i][i] = Variance;
    for (j=0;j<3;++j){
      if (i!=j){Gv.CovM[i][j] = 0.0;}
    } 
  }
  Cal_Inverse_Matrix3D_by_Cramer_Rule(Gv.iCovM, Gv.CovM, &(Gv.det));
  Gv.Cons = 1/(pow(2*M_PI,1.5)*sqrt(Gv.det));

  cal_sumDensity_and_Nvox_posi_dens(vox,&sumDensity,&Nvox_posi_dens);

 
  /** Assign gdf to each voxel **/
  n = 0;
  for (x=0;x<vox->N[0];++x){
    for (y=0;y<vox->N[1];++y){
      for (z=0;z<vox->N[2];++z){
        if (vox->dat[x][y][z]>0.0){
          Garray[n].num = n;
          Copy_GAUSS3D(&(Garray[n]),&Gv);
          Garray[n].M[0] = vox->OrigPos[0] + vox->grid_width * x;
          Garray[n].M[1] = vox->OrigPos[1] + vox->grid_width * y;
          Garray[n].M[2] = vox->OrigPos[2] + vox->grid_width * z;
          Garray[n].Weight = vox->dat[x][y][z]/sumDensity;
          n += 1;
        }
      }
    }
  }
 return(1);

} /* end of GaussMix_by_One_Voxel_One_Gdf_Assignment() */
