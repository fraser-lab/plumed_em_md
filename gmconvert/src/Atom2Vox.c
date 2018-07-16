/**

 <Atom2Vox.c>

 Functions from atomic data to Voxel data

**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdbool.h>
#include "pdbstruct.h"
#include "PdbIO.h"
#include "Radius.h"
#include "gauss.h"
#include "Matrix3D.h"
#include "Voxel.h"
#include "MCubeFunc.h"
#include "GaussIO.h"

/* one_over_2pi_32  = 1.0/pow(2*M_PI,3/2) = pow(1/2/M_PI,3/2) */
#define one_over_2pi_32   0.0634936359342410
/* three_over_2pi_32  = pow(3/2/M_PI,3/2) */
#define three_over_2pi_32 0.3299226101861591


/** FUNCTIONS (GLOBAL) **/
void Set_Voxel_Value_By_Atoms();
void Get_Min_Max_From_Atoms();
void Malloc_Voxel_From_Min_Max();
void Malloc_Voxel_From_Atoms();

/** FUNCTIONS (LOCAL) **/
static double Gaussian3D_Isotropic();
static double Gaussian3D_Isotropic_PeakDensityOne();

void Set_Voxel_Value_By_Atoms(vox,HeadAtom,Atom2VoxType,sigma_isotropic,init)
 struct VOXEL *vox;
 struct ATOM *HeadAtom;
 char   Atom2VoxType;  /* 'I':isotropic gaussian, 'R':Rvdw-based gaussian */ 
 float  sigma_isotropic;
 bool   init;
{
/*
 [Atom2VoxType=='I'] 
   gauss(r) = 1/pow(2pi,3/2)/sigma^3 * exp[-r^2/sigma^2/2.0] 
 [Atom2VoxType=='R'] 
   Based on the idea of Laskowski (J.Mol.Graph., 13,323-330, (1995) ).  
   Sigma is decided so that gauss(Rvdw) is half of gauss(0.0).
   Therefore,
      exp[-Rvdw^2/sigma^2/2]  = 1/2
   ==>  sigma^2  = Rvdw^2/(2log2) 
        sigma    = Rvdw/sqrt(2log2) 
  
   To simplify the problem, we set one to peak density of gauss(r). Threfore,we use:
      gauss(r) = exp[-r^2/sigma^2/2.0] 
      sigma    = Rvdw/sqrt(2log2) 
   For this case, the threshold value for the iscontour should be '0.5'.
*/
 int i,x,y,z;
 float X[3];
 struct ATOM *an; 
 float min,max,Sigma,SDcutoff;
 int Nmin[3],Nmax[3];

 printf("#Set_Voxel_Value_By_Atoms(Atom2VoxType:%c sigma_isotropic:%lf)\n",Atom2VoxType,sigma_isotropic);
 SDcutoff = 3.0;

 /** Initialize **/
 if(init){
  for (x=0;x< vox->N[0];++x){
    for (y=0;y< vox->N[1];++y){
      for (z=0;z< vox->N[2];++z){vox->dat[x][y][z] = 0.0; }
    }
  }  
 }

 /** Add density for each atoms **/
 Sigma = sigma_isotropic; 
 an = HeadAtom;
 while (an->next != NULL){
  an = an->next;
  if (Atom2VoxType=='R'){ Sigma = an->R/sqrt(2*log(2.0)); }

  for (i=0;i<3;++i){ 
    min = an->Pos[i] - an->R;
    max = an->Pos[i] + an->R;
    Nmin[i] = (int)floor((min - SDcutoff*Sigma - vox->OrigPos[i])/vox->grid_width);
    Nmax[i] = (int)ceil( (max + SDcutoff*Sigma - vox->OrigPos[i])/vox->grid_width);
    if (Nmin[i]<0)        Nmin[i] = 0;
    if (Nmax[i]>=vox->N[i]) Nmax[i] = vox->N[i]-1; }

   for (x=Nmin[0];x<=Nmax[0];++x){
    X[0] = vox->OrigPos[0] + vox->grid_width * x;
     for (y=Nmin[1];y<=Nmax[1];++y){
      X[1] = vox->OrigPos[1] + vox->grid_width * y;
       for (z=Nmin[2];z<=Nmax[2];++z){
        X[2] = vox->OrigPos[2] + vox->grid_width * z;
        if (Atom2VoxType=='R'){ vox->dat[x][y][z] += (float)Gaussian3D_Isotropic_PeakDensityOne(X,an->Pos,Sigma); }
        else { vox->dat[x][y][z] += (float)Gaussian3D_Isotropic(X,an->Pos,Sigma); }

       } /* z */
     } /* y */
   } /* x */

 } /* an */

} /* end of Set_Voxel_Value_By_Atoms() */


void Get_Min_Max_From_Atoms(HeadAtom,Min,Max)
 struct ATOM *HeadAtom;
 double (*Min)[3], (*Max)[3];
{
 int i;
 struct ATOM *an; 
 float min[3],max[3];
 char init;

 /** (1) Find min[3] and max[3] **/
 an = HeadAtom; init = 1;
 while (an->next != NULL){
  an = an->next;

  for (i=0;i<3;++i){ 
    min[i] = an->Pos[i] - an->R;
    max[i] = an->Pos[i] + an->R;
  }
 
  if (init==1){ 
    for (i=0;i<3;++i){ (*Min)[i] = min[i]; (*Max)[i] = max[i]; }
    init = 0;
  }
  else{ 
    for (i=0;i<3;++i){ 
      if (min[i]<(*Min)[i]) (*Min)[i] = min[i]; 
      if (max[i]>(*Max)[i]) (*Max)[i] = max[i];
   }
 }
}

}

void Malloc_Voxel_From_Atoms(vox,HeadAtom,margin)
 struct VOXEL *vox;
 struct ATOM *HeadAtom;
 float  margin;
{
 int i;
 struct ATOM *an; 
 float min[3],max[3],Min[3],Max[3];
 char init;
 printf("#Malloc_Voxel_Value_From_Atoms()\n");

 /** (1) Find min[3] and max[3] **/
 an = HeadAtom; init = 1;
 while (an->next != NULL){
  an = an->next;

  for (i=0;i<3;++i){ 
    min[i] = an->Pos[i] - an->R;
    max[i] = an->Pos[i] + an->R;
  }
 
  if (init==1){ 
    for (i=0;i<3;++i){ Min[i] = min[i]; Max[i] = max[i]; }
    init = 0;
  }
  else{ 
    for (i=0;i<3;++i){ 
      if (min[i]<Min[i]) Min[i] = min[i]; 
      if (max[i]>Max[i]) Max[i] = max[i];
   }
 }


 } /* an */
 printf("[0]Min %f Max %f [1]Min %f Max %f [2]Min %f Max %f\n",
 Min[0], Max[0], Min[1], Max[1], Min[2], Max[2]);

 /** (2) Set up N[] **/
 
 for (i=0;i<3;++i){ 
   vox->OrigPos[i] = (int)floor((Min[i]-margin)/vox->grid_width)*vox->grid_width; 
   vox->N[i] = (int)ceil((Max[i]- vox->OrigPos[i] + margin)/vox->grid_width);
 }
 
 /** (3) Malloc Voxel **/
 printf("#OrigPos %f %f %f\n",
 vox->OrigPos[0], vox->OrigPos[1], vox->OrigPos[2]);

 Malloc_Voxel(vox,vox->N[0],vox->N[1],vox->N[2]);

} /* end of Malloc_Voxel_From_Atoms() */

void Malloc_Voxel_From_Min_Max(vox,Min,Max,margin)
 struct VOXEL *vox;
 double Min[3], Max[3]; 
 float  margin;
{
 int i;
 printf("#Malloc_Voxel_Value_From_Min_Max()\n");

 /** (1) Set up N[] **/
 
 for (i=0;i<3;++i){ 
   vox->OrigPos[i] = (int)floor((Min[i]-margin)/vox->grid_width)*vox->grid_width; 
   vox->N[i] = (int)ceil((Max[i]- vox->OrigPos[i] + margin)/vox->grid_width);
 }
 
 /** (2) Malloc Voxel **/
 printf("#OrigPos %f %f %f\n",
 vox->OrigPos[0], vox->OrigPos[1], vox->OrigPos[2]);

 Malloc_Voxel(vox,vox->N[0],vox->N[1],vox->N[2]);

} /* end of Malloc_Voxel_From_Min_Max() */


double Gaussian3D_Isotropic(Pos,Cen,Sigma)
 float  Pos[3]; /* Position */
 float  Cen[3]; /* Center   */
 float  Sigma;  /* Standard Deviation */
{
 double D[3],xSx,val;
 
 D[0] = Pos[0] - Cen[0];
 D[1] = Pos[1] - Cen[1];
 D[2] = Pos[2] - Cen[2];
 xSx = 0.0;
 xSx = (D[0]*D[0] + D[1]*D[1] + D[2]*D[2])/(Sigma*Sigma);
 val =  three_over_2pi_32/(Sigma*Sigma*Sigma)*exp(-0.5*xSx);
 return(val);

}/* end of Gaussian3D_Isotropic() */


double Gaussian3D_Isotropic_PeakDensityOne(Pos,Cen,Sigma)
 float  Pos[3]; /* Position */
 float  Cen[3]; /* Center   */
 float  Sigma;  /* Standard Deviation */
{
 double D[3],xSx,val;
 
 D[0] = Pos[0] - Cen[0];
 D[1] = Pos[1] - Cen[1];
 D[2] = Pos[2] - Cen[2];
 xSx = 0.0;
 xSx = (D[0]*D[0] + D[1]*D[1] + D[2]*D[2])/(Sigma*Sigma);
 val =  exp(-0.5*xSx);
 return(val);

}/* end of Gaussian3D_Isotropic_PeakDensityOne() */

