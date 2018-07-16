/*

 <AtomKmean.c>
 K-means Clustering for Atom Data (PDB atom data)

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
#include "Voxel.h"
#include "MCubeFunc.h"
#include "PointEM.h"
#include "Matrix3D.h"
#include "GaussIO.h"

/*** Functions (GLOBAL) ***/
void K_Means_Clustering_For_Atoms();
void K_Means_Clustering_For_Atoms_Multiple_Start();

/*** Functions (LOCAL) ***/
static void update_Class();
static void update_M_Weight_CovM();

/*****************/
/*** FUNCTIONS ***/
/*****************/

void K_Means_Clustering_For_Atoms_Multiple_Start(Ngauss,Garray,Ahead,ObjFunc_fin,Ntry)
 int Ngauss;
 struct GAUSS3D *Garray;
 struct ATOM *Ahead;
 double *ObjFunc_fin; /* Final Objective Function */
 int    Ntry;
{
 struct GAUSS3D **GarrayTry;
 int t,g,t_min;
 double ObjFunc,ObjFunc_min; 

 printf("#K_Means_Clustering_For_Atoms_Multiple_Start(Ngauss:%d)\n",Ngauss);

 /* (1) Malloc Gauss Try */
 GarrayTry = (struct GAUSS3D **)malloc(sizeof(struct GAUSS3D*)*Ntry);
 for (t=0;t<Ntry;++t){
   GarrayTry[t] = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss);
 }

 /* (2) try K_Means Ntry times */
 t_min = -1;
 ObjFunc_min = -1.0;
 for (t=0;t<Ntry;++t){
   Initialize_Gauss_From_Randomly_Chosen_Atoms(Ahead,Ngauss,GarrayTry[t]);
   K_Means_Clustering_For_Atoms(Ngauss,GarrayTry[t],Ahead,&ObjFunc);
   if ((t==0)||(ObjFunc<ObjFunc_min)){ 
     ObjFunc_min  = ObjFunc;
     t_min = t;
   }
 } 

 /* printf("#t_min %d ObjFunc_min %f\n",t_min,ObjFunc_min); */

 /* (3) Copy the best gaussians */
 for (g=0;g<Ngauss;++g){
   Copy_GAUSS3D(&(Garray[g]),&(GarrayTry[t_min][g]));
 } 
 
  *ObjFunc_fin = ObjFunc_min; 

 /* (4) Free Gauss Try */
 for (t=0;t<Ntry;++t) free(GarrayTry[t]);
 free(GarrayTry);

} /* end of K_Means_Clustering_For_Atoms_Multiple_Start() */




void K_Means_Clustering_For_Atoms(Ngauss,Garray,Ahead,ObjFunc_fin)
 int Ngauss;
 struct GAUSS3D *Garray;
 struct ATOM *Ahead;
 double *ObjFunc_fin; /* Final Objective Function */
{
  int Nrepeat,r,a,Natom,Nclass_change; 
  int  *Class,*ClassPre; /* [0..Natom-1] */
  struct ATOM *an;
  double ObjFunc,sumWeight;

  /* printf("#K_Means_Clustering_For_Atoms()\n"); */
  Nrepeat = 1000; 

  /*** [1] Count Natom and malloc array class ***/
  Natom = 0;
  sumWeight = 0.0; 
  an = Ahead;
  while (an->next != NULL){
    an = an->next;
    an->num = Natom;
    ++Natom;
    sumWeight += an->Weight;
  }

  Class    = (int *)malloc(sizeof(int)*Natom);
  ClassPre = (int *)malloc(sizeof(int)*Natom);
  for (a=0;a<Natom;++a){ Class[a] = ClassPre[a] = 0;}

  /*** [2] Iteration(r) for K-means ***/
 
  r = 0;
  do{
   /* (2-1) Calculate Class[a] */ 
   update_Class(Ngauss, Garray, Ahead,Class,ClassPre,&Nclass_change,&ObjFunc);

   /* (2-2) Update M[], Weight */
   if (Nclass_change>0){
     update_M_Weight_CovM(Ngauss, Garray,Ahead,Class,sumWeight,'-');
   }

   for (a=0;a<Natom;++a){ClassPre[a] = Class[a];}

  }while ((r<Nrepeat)&&(Nclass_change>0));


  /* [3] Finishing. Update M[], Weight,CovM[][]  */
  update_M_Weight_CovM(Ngauss, Garray,Ahead,Class,sumWeight,'C');
  *ObjFunc_fin = ObjFunc; 

  free(Class);

} /* end of K_Means_Clustering_For_Atoms() */





void update_Class(Ngauss, Garray, Ahead,Class,ClassPre,Nclass_change,ObjFunc)
 int Ngauss;
 struct GAUSS3D *Garray;
 struct ATOM *Ahead;
 int  *Class;    /* [0..Natom-1] */
 int  *ClassPre; /* [0..Natom-1] */
 int  *Nclass_change;
 double *ObjFunc;
{
  struct ATOM *an;
  int g;
  double DD,minDD;
 
  *Nclass_change = 0;
  an = Ahead;
  DD = 0.0;
  minDD = -1.0;
  *ObjFunc = 0.0;
  while (an->next != NULL){ 
    an = an->next;
   
    for (g=0;g<Ngauss;++g){ 
      DD =  (an->Pos[0] - Garray[g].M[0])*(an->Pos[0] - Garray[g].M[0])
           +(an->Pos[1] - Garray[g].M[1])*(an->Pos[1] - Garray[g].M[1])
           +(an->Pos[2] - Garray[g].M[2])*(an->Pos[2] - Garray[g].M[2]);
      if ((DD<minDD)||(g==0)) { minDD = DD; Class[an->num] = g; }
    } 
 
    if (Class[an->num] != ClassPre[an->num]){ *Nclass_change = *Nclass_change + 1;} 
    *ObjFunc += DD;  
  }

} /* end of update_Class() */







void update_M_Weight_CovM(Ngauss, Garray, Ahead, Class,sumWeight,UpdateCovMtype)
 int Ngauss;
 struct GAUSS3D *Garray;
 struct ATOM *Ahead;
 int  *Class; /* [0..Natom-1] */
 double sumWeight;
 char   UpdateCovMtype;  /* 'C':update covariance matrix */ 
{
 int g,i,j;
 double sumWeight_g;
 struct ATOM *an;

 for (g=0;g<Ngauss;++g){ 
   sumWeight_g = 0.0; 
   for (i=0;i<3;++i){
     Garray[g].M[i] = 0.0;
     for (j=0;j<3;++j){
       Garray[g].CovM[i][j] = 0.0;
     }
   }

    /* update Garray[g].M[] and Garray[g].Weight */
    an = Ahead;
    while (an->next != NULL){ 
      an = an->next; 
      if (Class[an->num]==g){ 
        for (i=0;i<3;++i){Garray[g].M[i]  +=  an->Weight*an->Pos[i];}
        sumWeight_g += an->Weight;
      }   
    }

    if (sumWeight_g >0.0){ 
      for (i=0;i<3;++i) {Garray[g].M[i]  /=  sumWeight_g;} 
      Garray[g].Weight = sumWeight_g/sumWeight;  
    }
    else {Garray[g].Weight = 0.0;}  


    /* update Garray[g].CovM[] */
    if (UpdateCovMtype=='C'){
      an = Ahead;
      while (an->next != NULL){ 
        an = an->next; 
        if (Class[an->num]==g){ 
          for (i=0;i<3;++i){
            for (j=0;j<3;++j){
              if (i==j){ 
                Garray[g].CovM[i][j]  +=  an->Weight * (double)an->variance;
              }
              Garray[g].CovM[i][j]  +=  an->Weight * (an->Pos[i] - Garray[g].M[i])*(an->Pos[j] - Garray[g].M[j]);
            }
          }
        }   
      }

      if (sumWeight_g >0.0){ 
       for (i=0;i<3;++i) {
         for (j=0;j<3;++j) {
           Garray[g].CovM[i][j]  /=  sumWeight_g;
         } 
       }
       Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].iCovM, Garray[g].CovM,&(Garray[g].det));
       Garray[g].Cons = 1/(pow(2*M_PI,1.5)*sqrt(Garray[g].det));
      }
    }

  }  /* g loop */

} /* end of update_M_Weight_CovM() */


