/*

 <GmmAtomEM.c>

 for EM fit of Gaussian Mixture Model for for Gaussian Atom Observations 

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
#include "GaussIO.h"
#include "BasicIO.h"
#include "GridMap.h"

/*** Functions (GLOBAL) ***/
void EM_optimize_GaussMix_For_Gaussian_Atoms();
void GaussMix_by_One_Atom_One_Gdf_Assignment();
float CorrCoeff_Bwn_Atoms_and_GMM();


/*** Functions (LOCAL) ***/
static double Estep_of_GMM_For_Gaussian_Atoms();
static double Likelihood_for_Observed_Isotropic_GDF();
static double sumWeight_of_Atoms();
static void assign_gdf_by_one_atom();
/*
static double SqDistanceDoubleDouble();
static double SqDistanceFloatDouble();
static double SqDistanceFloatFloat();
*/




void EM_optimize_GaussMix_For_Gaussian_Atoms(Nrepeat,Ngauss,Garray,Ahead,logLike_fin,ConvThre)
 int Nrepeat;
 int Ngauss;
 struct GAUSS3D *Garray;
 struct ATOM *Ahead;
 double *logLike_fin; /* Final Log Likelihood */
 double ConvThre;  /* Threshold value for logLike convergence */
{
  int r,g,a,i,j;
  int Natom;
  struct ATOM *an;
  double logLike,logLikePre,DiffLogLike,sumWeight,w_post,sum_w_post;
  /* struct FLOAT2DMAP PostProb; */ /* [a][g]:Posterior probability of gdf [g] for atom [a].*/
  struct DOUBLE2DMAP PostProb; /* [a][g]:Posterior probability of gdf [g] for atom [a].*/
 
  printf("#EM_optimize_GaussMix_For_Gaussian_Atoms()\n");
  /*** [1] Count Natom and malloc matrix H ***/
  Natom = Number_Of_Atom(Ahead);
  sumWeight = sumWeight_of_Atoms(Ahead);
  /* Malloc_FLOAT2DMAP(&PostProb,Natom,Ngauss); */
  Malloc_DOUBLE2DMAP(&PostProb,Natom,Ngauss); 
  
  for (a=0;a<Natom;++a){
    for (g=0;g<Ngauss;++g){ PostProb.dat[a][g] = 0.0;}
  }


 logLikePre = Estep_of_GMM_For_Gaussian_Atoms(Ngauss,Garray,Ahead,sumWeight,'-',&PostProb);
 printf("#logLikeInit %lf %e\n",logLikePre,logLikePre);

 /*** [2] Iteration(r) for EM  ***/
 
 r = 0;
 do{
  /** (E-step) Calculate H[g][a] (posterior probablility )**/ 
  logLike = Estep_of_GMM_For_Gaussian_Atoms(Ngauss,Garray,Ahead,sumWeight,'U',&PostProb);

  /** (M-step) Update Weight,Mean[],CovM[][] **/ 
  for (g=0;g<Ngauss;++g){ 
     /* initialize */
    sum_w_post = 0.0;
    for (i=0;i<3;++i){
      Garray[g].M[i] = 0.0;
      for (j=0;j<3;++j){
        Garray[g].CovM[i][j] = 0.0;
      }
    }

    /* summation for atoms */
    an = Ahead;
    while (an->next != NULL){
      an = an->next;
      w_post = an->Weight*PostProb.dat[an->num][g] /sumWeight;
      sum_w_post += w_post;
      for (i=0;i<3;++i){ 
        Garray[g].M[i]  +=  w_post * an->Pos[i];
        for (j=0;j<3;++j){
          if (i==j){
            Garray[g].CovM[i][j] +=  w_post * ((double)an->variance + an->Pos[i] * an->Pos[j]);
          }
          else{
            Garray[g].CovM[i][j] +=  w_post * an->Pos[i] * an->Pos[j];
          }
        }
      }
    }


     /**  Update Weight,M[],CovM[][] **/ 
    Garray[g].Weight = sum_w_post;

    for (i=0;i<3;++i){ 
      Garray[g].M[i]  /=  sum_w_post;
    }
    
    for (i=0;i<3;++i){ 
      for (j=0;j<3;++j){
        Garray[g].CovM[i][j] /=  sum_w_post;
        Garray[g].CovM[i][j] -=  Garray[g].M[i]*Garray[g].M[j];
      }
    } 

   Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));
   Garray[g].Cons = 1/(pow(2*M_PI,1.5)*sqrt(Garray[g].det));

  } /* g */

  /* Calculate Likelihood */
  logLike = Estep_of_GMM_For_Gaussian_Atoms(Ngauss,Garray,Ahead,sumWeight,'-',&PostProb);
 
  DiffLogLike = logLike - logLikePre;
  if (PAR.OutStdLog=='T'){
    printf("Step %d logLike %lf Diff %lf\n",r,logLike,DiffLogLike);
  }
  logLikePre = logLike;
  *logLike_fin = logLike;
  ++r;
 
 } while ((r<Nrepeat)&&(DiffLogLike>ConvThre));

 printf("#Converged after %d steps (logLike_fin:%e)\n",r,*logLike_fin);

 /* Free_FLOAT2DMAP(&PostProb,Natom,Ngauss); */
 Free_DOUBLE2DMAP(&PostProb,Natom,Ngauss); 

} /* end of EM_optimize_GaussMix_For_Gaussian_Atoms() */





double Estep_of_GMM_For_Gaussian_Atoms(Ngauss,Garray,Ahead,sumWeight,UpdatePostProb,PostProb)
  int Ngauss;
  struct GAUSS3D *Garray;
  struct ATOM *Ahead;
  double sumWeight;
  char   UpdatePostProb;       /* 'U':pdate Post Prob */
  /* struct FLOAT2DMAP *PostProb; */
  struct DOUBLE2DMAP *PostProb; 
{
  struct ATOM *an;
  int g;
  double logLike,sumPostProb,atom_weight,postprob;
 
  logLike = 0.0;
  an = Ahead;
  while (an->next != NULL){ 
    an = an->next;
    /* printf("#>>%d %s %s R %f W %f\n",an->num,an->Atom,an->Resi,an->R,an->Weight); */
    sumPostProb = 0.0;
    for (g=0;g<Ngauss;++g){
     if (Garray[g].Weight>0.0){ 
       atom_weight = an->Weight/sumWeight;
  /*
       wLike = Garray[g].Weight * Likelihood_for_Observed_Isotropic_GDF(&(Garray[g]),an->Pos,inv_sqrt3*(double)an->R); 
       postprob = pow(wLike,Natom*atom_weight);
 */
       postprob = Garray[g].Weight * Likelihood_for_Observed_Isotropic_GDF(&(Garray[g]),an->Pos,(double)an->variance);
/*
       printf("#an %d g %d Garray[g].Weight %lf atom_weight %lf likelihood %lf\n",an->num,g,Garray[g].Weight,atom_weight,Likelihood_for_Observed_Isotropic_GDF(&(Garray[g]),an->Pos,(double)an->variance));
 */
       if (UpdatePostProb=='U'){
         PostProb->dat[an->num][g] = postprob;
       }
/*
       printf("#an %d g %d atom_weight %lf wLike %e PostProb %e\n",an->num,g,atom_weight,wLike,PostProb->dat[an->num][g]);
 */
       sumPostProb += postprob;
      }
    } 
    if (sumPostProb>0.0){
      if (UpdatePostProb=='U'){
        for (g=0;g<Ngauss;++g){
          PostProb->dat[an->num][g] /= sumPostProb;
         }
      }
      logLike += log(sumPostProb); 
    }
  } 
  return(logLike);
} /* end of Estep_of_GMM_For_Gaussian_Atoms() */











double Likelihood_for_Observed_Isotropic_GDF(G,floatX,variance)
 struct GAUSS3D *G;
 float  floatX[3];  /* Center for the observed GDF */
 double variance;   /* Variance for the observed gdf */
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
 xSx = 0.0;
 /*
   xSx += D[0]*(G->iCovM[0][0]*D[0] + G->iCovM[0][1]*D[1] + G->iCovM[0][2]*D[2]);
   xSx += D[1]*(G->iCovM[1][0]*D[0] + G->iCovM[1][1]*D[1] + G->iCovM[1][2]*D[2]);
   xSx += D[2]*(G->iCovM[2][0]*D[0] + G->iCovM[2][1]*D[1] + G->iCovM[2][2]*D[2]);
 */
 xSx += G->iCovM[0][0]*D[0]*D[0] + G->iCovM[1][1]*D[1]*D[1] + G->iCovM[2][2]*D[2]*D[2];
 xSx +=(G->iCovM[0][1]*D[0]*D[1] + G->iCovM[0][2]*D[0]*D[2] + G->iCovM[1][2]*D[1]*D[2]) * 2.0;

 trace = G->iCovM[0][0] + G->iCovM[1][1] + G->iCovM[2][2];
 val = G->Cons * exp(-0.5*variance*trace-0.5*xSx);
 return(val);
}/* end of Likelihood_for_Observed_Isotropic_GDF() */



double sumWeight_of_Atoms(Ahead)
  struct ATOM *Ahead;
{
  struct ATOM *an;
  double sumWeight;

  sumWeight = 0.0;
  an = Ahead;
  while (an->next != NULL){ 
    an = an->next;
    sumWeight += an->Weight;
  }
  return(sumWeight);
} /* end of sumWeight_of_Atoms() */


void GaussMix_by_One_Atom_One_Gdf_Assignment(Ngauss,Garray,Ahead)
 int Ngauss;
 struct GAUSS3D *Garray;
 struct ATOM *Ahead;
{
  int n,Natom;
  struct ATOM *an;
  double sumWeight;
 
  Natom     =  Number_Of_Atom(Ahead);
  sumWeight =  sumWeight_of_Atoms(Ahead);
  printf("#GaussMix_by_One_Atom_One_Gdf_Assignment(Ngauss:%d, Natom:%d, sumWeight:%lf)\n",Ngauss,Natom,sumWeight);

  n = 0;
  an = Ahead;
  while (an->next != NULL){
    an = an->next;
    assign_gdf_by_one_atom(&(Garray[n]),an,sumWeight);
    printf("#Weight %lf Gweight:%lf M %lf %lf %lf\n",an->Weight,Garray[n].Weight,Garray[n].M[0],Garray[n].M[1],Garray[n].M[2]); 
    n += 1;
  }

} /* end of GaussMix_by_One_Atom_One_Gdf_Assignment() */



float CorrCoeff_Bwn_Atoms_and_GMM(HeadAtom,Ngauss, Garray)
  struct ATOM *HeadAtom;
  int Ngauss;
  struct GAUSS3D *Garray;
{
  float CC;
  double varGG,varAA,covarGA,sumWeight;
  int g,h;
  struct ATOM *an,*bn;
  struct GAUSS3D gA,gB;

  printf("#CorrCoeff_Bwn_Atoms_and_GMM(Ngauss %d)\n",Ngauss);

  /*  [1] calculate varGG */
  varGG = 0.0;
  for (g=0;g<Ngauss;++g){
    varGG += Garray[g].Weight * Garray[g].Weight *
              Overlap_Integral_Bwn_Two_GAUSS3Ds(&(Garray[g]),&(Garray[g]));
    for (h=g+1;h<Ngauss;++h){
      varGG += 2*Garray[g].Weight * Garray[h].Weight *
                Overlap_Integral_Bwn_Two_GAUSS3Ds(&(Garray[g]),&(Garray[h]));
    }
  }
  printf("#varGG %lf\n",varGG);
  /* [2] calculate varAA */
  sumWeight = sumWeight_of_Atoms(HeadAtom);

  varAA = 0.0;
  bn = NULL;
  an = HeadAtom;
  while (an->next != NULL){
    an = an->next;
    assign_gdf_by_one_atom(&gA,an,sumWeight);
    varAA +=gA.Weight * gA.Weight * Overlap_Integral_Bwn_Two_Isotropic_GAUSS3Ds(&(gA),&(gA)); 

    bn = an;
    while (bn->next != NULL){
      bn = bn->next;
      assign_gdf_by_one_atom(&gB,bn,sumWeight);
      varAA += 2.0*gA.Weight * gB.Weight * Overlap_Integral_Bwn_Two_Isotropic_GAUSS3Ds(&(gA),&(gB)); 
    }
  }
  printf("#varAA %lf\n",varAA);

  /* [3] calculate covarGA */
  covarGA = 0.0;
  for (g=0;g<Ngauss;++g){
    an = HeadAtom;
    while (an->next != NULL){
      an = an->next;
      assign_gdf_by_one_atom(&gA,an,sumWeight);
      covarGA += Garray[g].Weight * gA.Weight *
                 Overlap_Integral_Bwn_Two_GAUSS3Ds(&(Garray[g]),&gA);
   }
  }
  printf("#covarGA %lf\n",covarGA);

  /* [4] calculate CC */
  printf("#varGG %e varAA %e covarGA %e\n",varGG, varAA, covarGA);

  CC = 0.0;
  if ((varGG>0.0)&&(varAA>0.0)){
    CC = covarGA/sqrt(varGG)/sqrt(varAA);
  }

  return(CC);
} /* end of CorrCoeff_Bwn_Atoms_and_GMM() */



void assign_gdf_by_one_atom(gauss,atom,sumWeight)
  struct GAUSS3D *gauss;
  struct ATOM    *atom;
  double sumWeight;
{
  int i,j;

  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      gauss->CovM[i][j]  = 0.0;
      gauss->iCovM[i][j] = 0.0;
    }
    gauss->M[i]     = atom->Pos[i]; 
    gauss->CovM[i][i]  = atom->variance; 
    gauss->iCovM[i][i] = 1.0/atom->variance; 
  }   

  gauss->Weight = atom->Weight/sumWeight;
  gauss->det    = atom->variance * atom->variance * atom->variance;
  gauss->Cons   = 1.0/(pow(2*M_PI,3.0/2.0)*sqrt(gauss->det));

} /* end of assign_gdf_by_one_atom() */

/*
double SqDistanceDoubleDouble(dposA,dposB,dd)
  double dposA[3],dposB[3];
  double dd[3];
{
  int i;
  double DD;

  DD = 0.0;
  for (i=0;i<3;++i){
    dd[i] = dposA[i] - dposB[i];
    dd[i] = dd[i]*dd[i];
    DD += dd[i];
  }
  return(DD);
}


double SqDistanceFloatDouble(fposA,dposB,dd)
  float  fposA[3];
  double dposB[3];
  double dd[3];
{
  int i;
  double DD;

  DD = 0.0;
  for (i=0;i<3;++i){
    dd[i] = fposA[i] - dposB[i];
    dd[i] = dd[i]*dd[i];
    DD += dd[i];
  }
  return(DD);
}

double SqDistanceFloatFloat(fposA,fposB,dd)
  float  fposA[3];
  float  fposB[3];
  double dd[3];
{
  int i;
  double DD;

  DD = 0.0;
  for (i=0;i<3;++i){
    dd[i] = fposA[i] - fposB[i];
    dd[i] = dd[i]*dd[i];
    DD += dd[i];
  }
  return(DD);
}
*/
