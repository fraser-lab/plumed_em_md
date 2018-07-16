/*

<MCubeVolSur.c>
 
 functions for calculating volumes or surfaces 
 for surface model generted by Marching Cubes.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include "globalvar.h" 
#include "Voxel.h" 
#include "mc_verface.h" 
#include "BasicIO.h" 
#include "MCubeFunc.h" 



/*** FUNCTIONS (GLOBAL) ***/
float Volume_Inside_Faces();
void Cal_Gcenter_Of_Vertexes();
float Area_Faces();
void MinMax_XYZ_of_MCVertex();
void setup_static_inout_variables_for_MCFaces();
void free_static_inout_variables_for_MCFaces();
void check_inside_MCFace_for_VOXEL();
void check_inside_MCFace_for_VOXEL_old();
void set_zero_for_higher_threshold_VOXEL();


/*** FUNCTIONS (LOCAL) ***/
static float VolumeTetrahedron();

static float **MCF_u;     /* [Nface][3]    :vector a->b */
static float **MCF_v;     /* [Nface][3]    :vector a->c */
static float **MCF_w;     /* [Nface][3]    : a->b cross_prod a->c */
static float ***MCF_imat; /* [Nface][2][2] : 
  inverse matrix for 
     |u^2 uv| 
     |uv u^2| 
 */



float Volume_Inside_Faces(grid_width,Vhead,Fhead)
 float  grid_width;
 struct MC_VERTEX *Vhead;
 struct MC_FACE   *Fhead;
{
 struct MC_FACE   *fn;
 float G[3],V,v;

 Cal_Gcenter_Of_Vertexes(G,Vhead);
 
 V = 0.0;
 fn = Fhead;
 while (fn->next != NULL){
  fn = fn->next;
  v = VolumeTetrahedron(fn->a->pos,fn->b->pos,fn->c->pos,G);
  /* printf("v %f\n",v); */
  V += v;
 }

 V *= grid_width*grid_width*grid_width;

 return(V);

} /* end of Volume_Inside_Faces() */




void Cal_Gcenter_Of_Vertexes(G,Vhead)
 float G[3]; /* Gcenter to be calculated. */
 /* CAUTION: unit of G[3] is "grid_width".
             The real cooditnate of Gcenter should be 
               (grid_width*G[0], grid_width*G[1], grid_width*G[2] ).
 */
 struct MC_VERTEX *Vhead;
{
 struct MC_VERTEX *vn;
 int Natom,i;

 G[0] = G[1] = G[2] = 0.0; Natom = 0;
 vn = Vhead;
 while (vn->next != NULL){
  vn = vn->next;
  ++Natom;
  for (i=0;i<3;++i) G[i] += vn->pos[i];
  }

  for (i=0;i<3;++i){ G[i] /= Natom;}

} /*  end of Cal_Gcenter_Of_Vertexes() */



float Area_Faces(grid_width,Fhead)
 float grid_width;
 struct MC_FACE   *Fhead;
{
 struct MC_FACE   *fn;
 int i;
 float X[3],Y[3],Z[3],S,DD,s;
 S = 0.0;

 fn = Fhead;
 while (fn->next != NULL){
   fn = fn->next;
   for (i=0;i<3;++i){ 
     X[i] = fn->a->pos[i] - fn->c->pos[i];
     Y[i] = fn->b->pos[i] - fn->c->pos[i];
   }
 
   Z[0] = X[1]*Y[2] - X[2]*Y[1];
   Z[1] = X[2]*Y[0] - X[0]*Y[2];
   Z[2] = X[0]*Y[1] - X[1]*Y[0];

   DD = Z[0]*Z[0] + Z[1]*Z[1] + Z[2]*Z[2];
   if (DD>0.0) s = 0.5 * sqrt(DD); else s = 0.0;
   S += s;
   /* printf("s %f\n",s);  */
 }

 S *= grid_width*grid_width;
 return(S);

} /* end of Area_Faces() */


float VolumeTetrahedron(a,b,c,x)
 float a[3],b[3],c[3],x[3];
 /*
 The order of a,b,c must be defined by counter-clockwise
 from the outside.

 if x is inside,  V is positive.
 if x is outside, V is negative.
 */
{
 float A[3],B[3],C[3],V;
 int i;

 for (i=0;i<3;++i){ 
   A[i] = a[i] - x[i];
   B[i] = b[i] - x[i];
   C[i] = c[i] - x[i];
 }

  V =   A[0]*B[1]*C[2] + A[1]*B[2]*C[0] + A[2]*B[0]*C[1]
      - A[2]*B[1]*C[0] - A[1]*B[0]*C[2] - A[0]*B[2]*C[1];

 return(V/6.0);

} /* end of Volume_Tetrahedron() */



void MinMax_XYZ_of_MCVertex(Min,Max,Vhead)
  float Min[3],Max[3];    /* (to be calculated) */
  struct MC_VERTEX *Vhead;
{
  struct MC_VERTEX *vn;
  unsigned char init;
  int i;
  vn = Vhead;
  init = 1;

  while (vn->next != NULL){
    vn = vn->next;
    for (i=0;i<3;++i){
      if ((init==1)||(vn->pos[i]<Min[i])) { Min[i] = vn->pos[i];}
      if ((init==1)||(vn->pos[i]>Max[i])) { Max[i] = vn->pos[i];}
    }
    init = 0;
  }

} /* end of MinMax_XYZ_of_MCVertex() */




void setup_static_inout_variables_for_MCFaces(Fhead)
  struct MC_FACE *Fhead;
{
  int Nface,i,j;
  struct MC_FACE *f;
  float uu,vv,uv,det;
  /** (1) Malloc inout_variables **/
  Nface = Number_Of_MC_FACE(Fhead);
  MCF_u    = (float **)malloc(sizeof(float *)*Nface);
  MCF_v    = (float **)malloc(sizeof(float *)*Nface);
  MCF_w    = (float **)malloc(sizeof(float *)*Nface);
  MCF_imat = (float ***)malloc(sizeof(float **)*Nface);

  for (i=0;i<Nface;++i){
    MCF_u[i] = (float *)malloc(sizeof(float)*3);
    MCF_v[i] = (float *)malloc(sizeof(float)*3);
    MCF_w[i] = (float *)malloc(sizeof(float)*3);
    MCF_imat[i] = (float **)malloc(sizeof(float *)*2);
    for (j=0;j<2;++j){
      MCF_imat[i][j] = (float *)malloc(sizeof(float)*2);
    }
  }

  /** (2) Setup inout_variables **/
  f = Fhead;
  while (f->next != NULL){
    f = f->next;
    for (i=0;i<3;++i){
      MCF_u[f->num][i] = f->b->pos[i] - f->a->pos[i];
      MCF_v[f->num][i] = f->c->pos[i] - f->a->pos[i];
    }
    Cross_Prod_F(MCF_w[f->num],MCF_u[f->num],MCF_v[f->num]);
    uu = Dot_Prod_F(MCF_u[f->num],MCF_u[f->num]);
    vv = Dot_Prod_F(MCF_v[f->num],MCF_v[f->num]);
    uv = Dot_Prod_F(MCF_u[f->num],MCF_v[f->num]);
    det = uu*vv - uv*uv;
    MCF_imat[f->num][0][0] =  vv/det;
    MCF_imat[f->num][0][1] = -uv/det;
    MCF_imat[f->num][1][0] = -uv/det;
    MCF_imat[f->num][1][1] =  uu/det;
  }

} /* end of setup_static_inout_variables_for_MCFaces() */


void free_static_inout_variables_for_MCFaces(Fhead)
  struct MC_FACE *Fhead;
{
  int Nface,i,j;

  Nface = Number_Of_MC_FACE(Fhead);

  for (i=0;i<Nface;++i){
    free(MCF_u[i]);
    free(MCF_v[i]);
    free(MCF_w[i]);
    for (j=0;j<2;++j){ free(MCF_imat[i][j]); }
    free(MCF_imat[i]);      
  }

  free(MCF_u);
  free(MCF_v);
  free(MCF_w);
  free(MCF_imat);

} /* end of free_static_inout_variables_for_MCFaces() */




void check_inside_MCFace_for_VOXEL_old(vox,Fhead,val_increment)
  struct VOXEL *vox;
  struct MC_FACE *Fhead;
  float  val_increment;  /* increment value for inside voxel */
{
  int i,j,k,m,Ncross;
  float d[3],p[3],pa[3],ax[3],wpa,wd,t,uax,vax,beta,gamma;
  struct MC_FACE *f;

  for (i=0;i<3;++i){  d[i] = -0.5 + (float)rand()/RAND_MAX; }
  d[0] = 1.0;
  d[1] = 0.0;
  d[2] = 0.0;

  for (i=0;i<vox->N[0];++i){
    p[0] = vox->OrigPos[0] + vox->grid_width * i;
    printf("[%d/%d]\n",i,vox->N[0]);
    for (j=0;j<vox->N[1];++j){
      p[1] = vox->OrigPos[1] + vox->grid_width * j;
      for (k=0;k<vox->N[2];++k){
        p[2] = vox->OrigPos[2] + vox->grid_width * k;
        Ncross = 0;
        f = Fhead;
        while (f->next != NULL){
          f = f->next;
          for (m=0;m<3;++m){ pa[m] = f->a->pos[m] - p[m];}
          wpa = Dot_Prod_F(MCF_w[f->num],pa);
          wd  = Dot_Prod_F(MCF_w[f->num],d);
          t = wpa/wd;
/*
          printf("#[%2d %2d %2d][%f %f %f] w %f %f %f wpa %f wd %f t %f\n",i,j,k,p[0],p[1],p[2], MCF_w[f->num][0], MCF_w[f->num][1], MCF_w[f->num][2], wpa,wd,t);
*/
          if (t>=0.0){
            for (m=0;m<3;++m){ ax[m] = t * d[m] + p[m] - f->a->pos[m];}
            uax = Dot_Prod_F(MCF_u[f->num],ax);
            vax = Dot_Prod_F(MCF_v[f->num],ax);
            beta  = MCF_imat[f->num][0][0] * uax +  MCF_imat[f->num][0][1] * vax;
            gamma = MCF_imat[f->num][1][0] * uax +  MCF_imat[f->num][1][1] * vax;
            if ((0.0<=beta)&&(beta<=1.0)&&(0.0<=gamma)&&(gamma<=1.0)&&((beta+gamma)<=1.0)){
              /* printf("beta %f gamma %f beta+gamma %f\n",beta,gamma,beta+gamma); */
              Ncross += 1;
            }
          }
        }

         /* printf("#[%2d %2d %2d][%f %f %f] Ncross %d\n",i,j,k,p[0],p[1],p[2], Ncross); */
        /* if (Ncross>0) {printf("#[%2d %2d %2d] Ncross %d\n",i,j,k,Ncross);  } */
        if ((Ncross%2)==1) { vox->dat[i][j][k] += val_increment;} 

      }
    }
  }

} /* end of check_inside_MCFace_for_VOXEL_old() */





void check_inside_MCFace_for_VOXEL(vox,Fhead,val_increment)
  struct VOXEL *vox;
  struct MC_FACE *Fhead;
  float  val_increment;  /* increment value for inside voxel */
{
  int i,j,k,m,ti,Ncross;
  float d[3],p[3],pa[3],ax[3],wpa,wd,t,uax,vax,beta,gamma;
  struct MC_FACE *f;
  int Ntri_intersect;
  struct MC_FACE *TRI_INTERSECT[100];

  d[0] = 0.0; d[1] = 0.0; d[2] = 1.0;

  for (i=0;i<vox->N[0];++i){
    p[0] = vox->OrigPos[0] + vox->grid_width * i;
    /* printf("[%d/%d]\n",i,vox->N[0]); */
    for (j=0;j<vox->N[1];++j){
      p[1] = vox->OrigPos[1] + vox->grid_width * j;

      p[2] = vox->OrigPos[2];
      /* [1] Check intersect triangle through p = (i,j,0) */
      Ntri_intersect = 0;
      f = Fhead;
      while (f->next != NULL){
        f = f->next;
        for (m=0;m<3;++m){ pa[m] = f->a->pos[m] - p[m];}
        wpa = Dot_Prod_F(MCF_w[f->num],pa);
        wd  = Dot_Prod_F(MCF_w[f->num],d);
        t = wpa/wd;
        for (m=0;m<3;++m){ ax[m] = t * d[m] + p[m] - f->a->pos[m];}
        uax = Dot_Prod_F(MCF_u[f->num],ax);
        vax = Dot_Prod_F(MCF_v[f->num],ax);
        beta  = MCF_imat[f->num][0][0] * uax +  MCF_imat[f->num][0][1] * vax;
        gamma = MCF_imat[f->num][1][0] * uax +  MCF_imat[f->num][1][1] * vax;
        if ((0.0<=beta)&&(beta<=1.0)&&(0.0<=gamma)&&(gamma<=1.0)&&((beta+gamma)<=1.0)){
          TRI_INTERSECT[Ntri_intersect] = f; 
          Ntri_intersect += 1;
        }
      }

      /* printf("#[%d %d *] Ntri_intersect %d\n",i,j,Ntri_intersect); */

      /** [2] Check intersect triangle only for in TRI_INTESECT[] */
      for (k=0;k<vox->N[2];++k){
        p[2] = vox->OrigPos[2] + vox->grid_width * k;
        Ncross = 0;
        for (ti=0;ti<Ntri_intersect;++ti){
          f = TRI_INTERSECT[ti]; 
          for (m=0;m<3;++m){ pa[m] = f->a->pos[m] - p[m];}
          wpa = Dot_Prod_F(MCF_w[f->num],pa);
          wd  = Dot_Prod_F(MCF_w[f->num],d);
          t = wpa/wd;
          /* printf("#t %f\n",t); */
          if (t>=0.0){
            Ncross += 1;
          }
        }

        if ((Ncross%2)==1) { vox->dat[i][j][k] += val_increment;} 
      }

    }
  } 

} /* end of check_inside_MCFace_for_VOXEL() */








void set_zero_for_higher_threshold_VOXEL(vox,threshold)
  struct VOXEL *vox;
  float threshold;
{
  int i,j,k,Nzero;
  Nzero = 0;
  for (i=0;i<vox->N[0];++i){
    for (j=0;j<vox->N[1];++j){
      for (k=0;k<vox->N[2];++k){
        if (vox->dat[i][j][k]>=threshold){
          vox->dat[i][j][k] = 0.0;
          Nzero += 1;
        }
      }
    }
  }
 printf("#set_zero_for_higher_threshold_VOXEL(threshold:%lf Nzero:%d)\n",threshold,Nzero);

} /* end of set_zero_for_higher_threshold_VOXEL() */

