/*

 <3DmapKmean.c>
  K-means Clustering for Point Data

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "pdbstruct.h"
#include "PdbIO.h"
#include "Radius.h"
#include "gauss.h"
#include "GaussIO.h"
#include "Matrix3D.h"
#include "Voxel.h"
#include "MCubeFunc.h"
#include "PointEM.h"
#include "3DmapEM.h"


/** Functions(global) **/
void K_Means_Clustering_for_3D_Map_Multiple_Start();
double K_Means_Clustering_for_3D_Map();
void Initialize_Gauss_From_Randomly_Chosen_Voxels();

/** Functions(local) **/


/*****************/
/*** FUNCTIONS ***/
/*****************/


void K_Means_Clustering_for_3D_Map_Multiple_Start(Ngauss,Garray,vox,ObjFunc_fin,Ntry,VarType,Cgg2var,resolution)
 int Ngauss;
 struct GAUSS3D *Garray;
 struct VOXEL  *vox;
 double *ObjFunc_fin; /* Final Objective Function */
 int  Ntry;
 char   VarType;      /* G: Var = Cgg2var * grid_width*grid_width R: var = (resolution/2)^2 */
 double Cgg2var;      /* Constant bwn variance and grid_width. Variance = Cgg2var * grid_width * grid_width */
 double resolution;   /* Var = (resolution/2)^2 */
{
 struct GAUSS3D *GarrayTry;
 int t,t_min;
 double ObjFunc_min,ObjFunc;
 
 printf("#K_Means_Clustering_for_3D_Map_Multiple_Start(%d %d %d Ngauss %d Ntry %d)\n",
   vox->N[0],vox->N[1],vox->N[2],Ngauss,Ntry);
 
 /* (1) Malloc Gauss Try */
 GarrayTry = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss);

 /* (2) try K_Means Ntry times */
 t_min = -1;
 ObjFunc_min = -1.0;
 for (t=0;t<Ntry;++t){
   Initialize_Gauss_From_Randomly_Chosen_Voxels(vox,Ngauss,GarrayTry);
   ObjFunc = K_Means_Clustering_for_3D_Map(Ngauss,GarrayTry,vox);
   printf("#try %d ObjFunc %lf\n",t,ObjFunc);
   if ((t==0)||(ObjFunc<ObjFunc_min)){ 
     ObjFunc_min  = ObjFunc;
     t_min = t;
     Copy_GAUSS3D_Array(Ngauss,Garray,GarrayTry);
   } 
 }
 
 printf("#t_min %d ObjFunc_min %f\n",t_min,ObjFunc_min);

 /* (3) Free Gauss Try */
 free(GarrayTry);



} /* end of K_Means_Clustering_for_3D_Map_Multiple_Start() */





 
double K_Means_Clustering_for_3D_Map(Ngauss,Garray,vox,VarType,Cgg2var,resolution)
 int Ngauss;
 struct GAUSS3D *Garray;
 struct VOXEL  *vox;
 char   VarType;      /* G: Var = Cgg2var * grid_width*grid_width R: var = (resolution/2)^2 */
 double Cgg2var;      /* Constant bwn variance and grid_width. Variance = Cgg2var * grid_width * grid_width */
 double resolution;   /* Var = (resolution/2)^2 */
{
 int x,y,z,i,j,g,r;
 int ***Class,***ClassPre;  /* class number for each grid [x][y][z] */
 int Nrepeat,Nclass_change;
 double DD,minDD,Dpos[3];
 double SumDensAll,*SumDensClass; /*  [Ngauss] */
 double *X,*Y,*Z; /* X,Y,Z coordinates [x], [y], [z] */
 double tole,minVar;
 double ObjFunc_fin; /* Final Objective Function */
 double Variance;
 
 /*
 printf("#K_Means_Clustering_for_3D_Map(Ngauss %d)\n",Ngauss);
 */

 if ((VarType=='R') && (resolution>0.0)){
   Variance = (resolution/2.0)*(resolution/2.0);
 }
 else{
   Variance = Cgg2var * vox->grid_width * vox->grid_width;
 }



 minVar = (0.5*vox->grid_width)*(0.5*vox->grid_width);

 Nrepeat = 1000;
 tole = 0.000001;
 /* [1] Malloc Class[][][], ClassPre[][][], SumDensClass[], X[],Y[],Z[] */
 SumDensClass = (double *)malloc(sizeof(double)*Ngauss);
 
 X = (double *)malloc(sizeof(double)*vox->N[0]);
 Y = (double *)malloc(sizeof(double)*vox->N[1]);
 Z = (double *)malloc(sizeof(double)*vox->N[2]);
  
 Class    = (int ***)malloc(sizeof(int **)*vox->N[0]);
 ClassPre = (int ***)malloc(sizeof(int **)*vox->N[0]);
 for (x=0;x<vox->N[0];++x){
   Class[x]    = (int **)malloc(sizeof(int *)*vox->N[1]);
   ClassPre[x] = (int **)malloc(sizeof(int *)*vox->N[1]);
   for (y=0;y<vox->N[1];++y){ 
     Class[x][y]    = (int *)malloc(sizeof(int)*vox->N[2]);
     ClassPre[x][y] = (int *)malloc(sizeof(int)*vox->N[2]);
    }
 }

 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){Class[x][y][z] = ClassPre[x][y][z] = 0;}
   }
 }

 for (x=0;x<vox->N[0];++x){X[x] = vox->OrigPos[0] + vox->grid_width * x;}
 for (y=0;y<vox->N[1];++y){Y[y] = vox->OrigPos[1] + vox->grid_width * y;}
 for (z=0;z<vox->N[2];++z){Z[z] = vox->OrigPos[2] + vox->grid_width * z;}

 /*** [2] Iteration(r) for K-means ***/
  
 r = 0;
 minDD = -1.0;
 do{
     Nclass_change = 0;
     SumDensAll = 0.0;
     for (g=0;g<Ngauss;++g){ SumDensClass[g] = 0.0; } 

     /** (i) Calculate Class[x][y][z] **/
     for (x=0;x<vox->N[0];++x){
       for (y=0;y<vox->N[1];++y){
         for (z=0;z<vox->N[2];++z){
           if (vox->dat[x][y][z]>0.0){
             for (g=0;g<Ngauss;++g){
               DD =  (X[x]-Garray[g].M[0])*(X[x]-Garray[g].M[0])
                    +(Y[y]-Garray[g].M[1])*(Y[y]-Garray[g].M[1])
                    +(Z[z]-Garray[g].M[2])*(Z[z]-Garray[g].M[2]);
               if ((DD<minDD)||(g==0)) {minDD = DD; Class[x][y][z] = g;}
             } 
            if (Class[x][y][z] != ClassPre[x][y][z]) ++Nclass_change;
            SumDensClass[Class[x][y][z]] += vox->dat[x][y][z];
          } 
        } 
      } 
    } 
  
 
     /** (ii) Output WARNING if SumDensClass[g]==0.0 **/
     for (g=0;g<Ngauss;++g){ 
       if (SumDensClass[g]==0.0){
         printf("#WARNING gcluster %d is SumDensClass=%lf (M %f %f %f)\n",g,SumDensClass[g],Garray[g].M[0],Garray[g].M[1],Garray[g].M[2]); 
         /*
         printf("%f..%f   %f..%f   %f..%f\n",
            vox->OrigPos[0], vox->OrigPos[0] + vox->grid_width * vox->N[0],
            vox->OrigPos[1], vox->OrigPos[1] + vox->grid_width * vox->N[1],
            vox->OrigPos[2], vox->OrigPos[2] + vox->grid_width * vox->N[2]);
         */
       } 
     } 


 
     /** (iii) Update Mean[g] (only if SumDensClass[g]>0.0) **/
     for (g=0;g<Ngauss;++g){ 
       if (SumDensClass[g]>0.0){
         for (i=0;i<3;++i){ Garray[g].M[i] = 0.0;}
       }
     } 

     for (x=0;x<vox->N[0];++x){
       for (y=0;y<vox->N[1];++y){
         for (z=0;z<vox->N[2];++z){
           if (vox->dat[x][y][z]>0.0){
             g = Class[x][y][z];
             Garray[g].M[0] += X[x] * vox->dat[x][y][z];
             Garray[g].M[1] += Y[y] * vox->dat[x][y][z];
             Garray[g].M[2] += Z[z] * vox->dat[x][y][z];
             ClassPre[x][y][z] = Class[x][y][z];  
           }
        }
      }
     } 
   
     for (g=0;g<Ngauss;++g){ 
       if (SumDensClass[g]>0.0){
         for (i=0;i<3;++i) { Garray[g].M[i] /= SumDensClass[g]; }
         SumDensAll += SumDensClass[g];
       }
     } 
   
   /*
   printf("r %d Nclass_change %d\n",r,Nclass_change);
   */
   ++r; 


 } while ((r<Nrepeat)&&(Nclass_change>0));


 /** [3] Caluclate Final Objective Function **/
 ObjFunc_fin = 0.0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){
       g = Class[x][y][z];
       DD =  (X[x]-Garray[g].M[0])*(X[x]-Garray[g].M[0])
            +(Y[y]-Garray[g].M[1])*(Y[y]-Garray[g].M[1])
            +(Z[z]-Garray[g].M[2])*(Z[z]-Garray[g].M[2]);
 
       ObjFunc_fin += vox->dat[x][y][z] * DD;
    } 
   }
 }

 /** [4] Caluclate Covar Matrix (CovM) for each cluster **/
 
 for (g=0;g<Ngauss;++g){ 
   SumDensClass[g] = 0.0;
   for (i=0;i<3;++i){
     for (j=0;j<3;++j){Garray[g].CovM[i][j] = 0.0;}
   }
 }
  for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
    for (z=0;z<vox->N[2];++z){
      if (vox->dat[x][y][z]>0.0){
        g = Class[x][y][z];
        SumDensClass[g] += vox->dat[x][y][z];
        Dpos[0] = X[x] - Garray[g].M[0];
        Dpos[1] = Y[y] - Garray[g].M[1];
        Dpos[2] = Z[z] - Garray[g].M[2];
        for (i=0;i<3;++i){
          for (j=0;j<3;++j){ 
            Garray[g].CovM[i][j] += vox->dat[x][y][z]*Dpos[i]*Dpos[j]; 
            /*
            if (i==j){Garray[g].CovM[i][j] += Variance;}
            Garray[g].CovM[i][j] += vox->dat[x][y][z]*Dpos[i]*Dpos[j]; 
            */
          }
        }
      }
    }
  }
 }

 for (g=0;g<Ngauss;++g){
   for (i=0;i<3;++i){
     for (j=0;j<3;++j){
      Garray[g].CovM[i][j] /= SumDensClass[g];
     }
   } 
   /** Enforce minimum Variance (1/2*grid_width)^2 for CovM[0][0],CovM[1][1],CovM[2][2] **/
   for (i=0;i<3;++i){
     if (Garray[g].CovM[i][i]<minVar){Garray[g].CovM[i][i]=minVar;}
   } 
   Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));

   /* if rank[CovM] < 3, then add tole into diagonal component, and set Weight=0.0. */
   if (fabs(Garray[g].det)<tole){
      for (i=0;i<3;++i){Garray[g].CovM[i][i] += tole;}
      Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));
      Garray[g].Weight = 0.0;
      Garray[g].Cons = 1/(pow(2*M_PI,1.5)*sqrt(Garray[g].det));
   }
   else{ 
     Garray[g].Weight = SumDensClass[g]/SumDensAll; 
     Garray[g].Cons = 1/(pow(2*M_PI,1.5)*sqrt(Garray[g].det));
   }
 }
 
 /** [4] Free Class[][][] and ClassPre[][][][] **/

 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){ 
     free(Class[x][y]);
     free(ClassPre[x][y]);
   }
 }

 for (x=0;x<vox->N[0];++x){ 
   free(Class[x]);
   free(ClassPre[x]);
 }

 free(Class); free(ClassPre);
 free(SumDensClass);
 return(ObjFunc_fin);

} /* end of K_Means_Clustering_for_3D_Map() */






void Initialize_Gauss_From_Randomly_Chosen_Voxels(vox,Ngauss,Garray)
  struct VOXEL  *vox;
  int Ngauss;
  struct GAUSS3D *Garray;
{
  struct VOXEL cvox;  /* voxels for cumulative density map */
  double  sum_density,cum_prob,randnum;
  int x,y,z,i,g;
  char hit;

  for (g=0;g<Ngauss;++g){
    Garray[g].num    = g;
    Garray[g].Weight = 0.0;
  }
 
  Malloc_Voxel(&cvox,vox->N[0],vox->N[1],vox->N[2]);

  /*[1] cal sum_density for vox */
  sum_density = 0.0;
  for (x=0;x<vox->N[0];++x){
    for (y=0;y<vox->N[1];++y){
      for (z=0;z<vox->N[2];++z){
        if (vox->dat[x][y][z]>0.0){
          sum_density += vox->dat[x][y][z];
        }
      }
    } 
  } 

  /*[2] set up cumulative density for cvox */
  cum_prob = 0.0;
  for (x=0;x<vox->N[0];++x){
    for (y=0;y<vox->N[1];++y){
      for (z=0;z<vox->N[2];++z){
        if (vox->dat[x][y][z]>0.0){
          cum_prob += (double)vox->dat[x][y][z]/sum_density; 
        }
        cvox.dat[x][y][z] = (float)cum_prob;
      }
    } 
  } 

  /* [3] Random chosen Ngauss points from the cvox voxels */
  for (i=0;i<Ngauss;++i){
    randnum = (double)random()/(double)RAND_MAX;
    hit = 0;
    for (x=0;x<vox->N[0];++x){
      for (y=0;y<vox->N[1];++y){
        for (z=0;z<vox->N[2];++z){
          if (randnum <= cvox.dat[x][y][z]){
             Garray[i].M[0] = vox->OrigPos[0] + vox->grid_width * x;
             Garray[i].M[1] = vox->OrigPos[1] + vox->grid_width * y;
             Garray[i].M[2] = vox->OrigPos[2] + vox->grid_width * z;
             /* printf("#gnum %d randnum %lf cvox %f xyz %d %d %d vox %f\n",i,randnum,cvox.dat[x][y][z],x,y,z,vox->dat[x][y][z]); */
             x = vox->N[0] + 1;
             y = vox->N[1] + 1;
             z = vox->N[2] + 1;
             hit = 1;
          }
        }
      }
    }
    if (hit==0){
      printf("#ERROR(Initialize_Gauss_From_Randomly_Chosen_Voxels()): NoHit for  gnum %d randnum %lf\n",i,randnum);
      exit(1);
    }
  }

  Free_Voxel(&cvox);
} /* end of Initialize_Gauss_From_Randomly_Chosen_Voxels() */



