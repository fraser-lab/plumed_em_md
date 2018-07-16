/*

<MCubeFunc.c>
 
 functions for calculating surface data from voxel data
 using Marching-Cube algorithm of the tetrahedron

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

/** Connection Matrix **/
struct CMATRIX{
 int N;
 char **m;
};


/*** FUNCTIONS (GLOBAL) ***/
void Marching_Cube_Tetrahedral();
void Cal_Gcenter_Of_Vertexes();
void Free_MC_Vertex_Stack();
void Free_MC_Face_Stack();
struct MC_VERTEX *Add_MC_Vertex_Stack();
struct MC_FACE *Add_MC_Face_Stack();
struct MC_FACE *Add_MC_Face_Stack_Simple();

void Sub_Vec_F();
void Sub_Vec_FI();
void Cross_Prod_F();
float Dot_Prod_F();
void Make_Mid_Point();
struct MC_VERTEX *Add_Pnt_to_Ver_Map();
struct MC_VERTEX *Find_Ver_Map();
void Malloc_CMATRIX();
void Free_CMATRIX();
void Malloc_MIDMAP_YZ();
void Free_MIDMAP_YZ();
void Reset_MIDMAP_YZ();
int Number_Of_MC_VERTEX();
int Number_Of_MC_FACE();
void Free_MC_Edge_Stack();
struct MC_EDGE *Add_MC_Edge_Stack();
int  check_MC_Vertex_neighbors();
int Check_Overlap_MC_Edge();
int common_xyz_in_two_vertexes();


/*** FUNCTIONS (LOCAL) ***/
static int  march_each_teterahedron();


void Marching_Cube_Tetrahedral(vox,Vhead,Fhead,Ehead,thre,MidPntType)
 struct VOXEL  *vox;
 struct MC_VERTEX *Vhead;
 struct MC_FACE   *Fhead;
 struct MC_EDGE   *Ehead;
 float thre;
 char  MidPntType; /* 'H'alf-point, 'S'imple interpolation, 'C'omplicate Interpolation */
{
 struct MIDMAP_YZ Map0,Map1,*Curr,*Post;
 int x,y,z;
 int a[3],b[3],c[3],d[3],e[3],f[3],g[3],h[3]; /* Eight vertex of cube */
 
 /*
  To decide a new mid-point is already constructed or not,
  two maps (map[y][z]) are used : Map0 and Map1.
  if ((x%2)==0)  {Cur = Map0; Post = Map1;}
     (otherwise) {Cur = Map1; Post = Map0;}
 */

 printf("#Marching_Cube_Tetrahedral(%f %e,midpnttype:%c)\n",thre,thre,MidPntType);
 printf("MC_THRESHOLD_DENSITY %e\n",thre);

 Free_MC_Vertex_Stack(Vhead);
 Free_MC_Face_Stack(Fhead);
 
 Vhead->num = -1; Vhead->next = NULL;
 Ehead->next = NULL; 
 /* Fhead->num = -1; */
 Fhead->next = NULL;
 Malloc_MIDMAP_YZ(&Map0,vox->N[1], vox->N[2]);
 Reset_MIDMAP_YZ(&Map0);
 Malloc_MIDMAP_YZ(&Map1,vox->N[1], vox->N[2]);
 Reset_MIDMAP_YZ(&Map1);

 /* printf("#grid %f X %d Y %d Z %d minX %f minY %f minZ %f \n",
   vox->grid_width,vox->N[0],vox->N[1],vox->N[2],vox->minX,vox->minY,vox->minZ);  */

 /*** Marching Cube ***/
 for (x=0;x<(vox->N[0]-1);++x){ 
    if ((x%2)==0) {Curr = &Map0; Post = &Map1; Reset_MIDMAP_YZ(&Map1);}
             else {Curr = &Map1; Post = &Map0; Reset_MIDMAP_YZ(&Map0);}

   printf("#%d ",x); 
   if (Vhead->next!=NULL) printf(" %d ",Vhead->next->num);
   fflush(stdout);
   if (((x+1)%5)==0) {printf("\n"); fflush(stdout);} 
   if (x==(vox->N[0]-2)) printf("\n"); 
  
   for (y=0;y<(vox->N[1]-1);++y){
     for (z=0;z<(vox->N[2]-1);++z){ 
      /* Eight vertices of the focused cube */ 
       a[0] = x;   a[1] = y;   a[2] = z;   
       b[0] = x+1; b[1] = y;   b[2] = z;   
       c[0] = x+1; c[1] = y+1; c[2] = z;   
       d[0] = x;   d[1] = y+1; d[2] = z;   
       e[0] = x;   e[1] = y;   e[2] = z+1; 
       f[0] = x+1; f[1] = y;   f[2] = z+1; 
       g[0] = x+1; g[1] = y+1; g[2] = z+1; 
       h[0] = x;   h[1] = y+1; h[2] = z+1; 
   
       if (((x+y+z)%2)==0){
         march_each_teterahedron(x,a,b,c,f,vox,Vhead,Fhead,Ehead,Curr,Post,thre,MidPntType);
         march_each_teterahedron(x,f,g,h,c,vox,Vhead,Fhead,Ehead,Curr,Post,thre,MidPntType);
         march_each_teterahedron(x,f,e,h,a,vox,Vhead,Fhead,Ehead,Curr,Post,thre,MidPntType);
         march_each_teterahedron(x,a,c,d,h,vox,Vhead,Fhead,Ehead,Curr,Post,thre,MidPntType);
         march_each_teterahedron(x,f,c,h,a,vox,Vhead,Fhead,Ehead,Curr,Post,thre,MidPntType);
        }
       else{
         march_each_teterahedron(x,e,f,g,b,vox,Vhead,Fhead,Ehead,Curr,Post,thre,MidPntType);
         march_each_teterahedron(x,b,c,d,g,vox,Vhead,Fhead,Ehead,Curr,Post,thre,MidPntType);
         march_each_teterahedron(x,e,g,h,d,vox,Vhead,Fhead,Ehead,Curr,Post,thre,MidPntType);
         march_each_teterahedron(x,a,b,d,e,vox,Vhead,Fhead,Ehead,Curr,Post,thre,MidPntType);
         march_each_teterahedron(x,e,g,b,d,vox,Vhead,Fhead,Ehead,Curr,Post,thre,MidPntType);
       }
     } /* z */
   } /* y */
 } /* x */

 /* printf("\n#done\n"); */

 Free_MIDMAP_YZ(&Map1);
 Free_MIDMAP_YZ(&Map0);
 
} /* end of Marching_Cube_Tetrahedral() */




int march_each_teterahedron(cur_x,a,b,c,d,vox,Vhead,Fhead,Ehead,Curr,Post,thre,MidPntType)
 int cur_x;  /* current x */
 int a[3],b[3],c[3],d[3];
 struct VOXEL *vox; 
 struct MC_VERTEX *Vhead;
 struct MC_FACE   *Fhead;
 struct MC_EDGE   *Ehead;
 struct MIDMAP_YZ *Curr,*Post;
 float thre;
 char  MidPntType; /* 'H'alf-point, 'S'imple interpolation, 'C'omplicate Interpolation */
{
 int Nin,Nout;
 int IN[4][3],OUT[4][3];   /* Lattice Number for INside and OUTside vertex */
 struct MC_VERTEX *pn,*qn,*rn,*sn;
 char A,B,C,D;  /* if >thre 1, otherwise 0 for a,b,c,d points */

 /*
     a
    /|\
   / | \
  /  |  \
 b---|---d 
  \  |  / 
   \ | /
    \|/
     c
*/

 
 if (vox->dat[a[0]][a[1]][a[2]] > thre) A = 1; else A = 0;
 if (vox->dat[b[0]][b[1]][b[2]] > thre) B = 1; else B = 0;
 if (vox->dat[c[0]][c[1]][c[2]] > thre) C = 1; else C = 0;
 if (vox->dat[d[0]][d[1]][d[2]] > thre) D = 1; else D = 0;

 /* printf("#A %d B %d C %d D %d\n",A,B,C,D); */

 Nin = Nout = 0;
 if (A==1) {IN[Nin][0]  =a[0]; IN[Nin][1]  =a[1]; IN[Nin][2]  =a[2];  ++Nin;} 
      else {OUT[Nout][0]=a[0]; OUT[Nout][1]=a[1]; OUT[Nout][2]=a[2]; ++Nout;} 
 if (B==1) {IN[Nin][0]  =b[0]; IN[Nin][1]  =b[1]; IN[Nin][2]  =b[2];  ++Nin;} 
      else {OUT[Nout][0]=b[0]; OUT[Nout][1]=b[1]; OUT[Nout][2]=b[2]; ++Nout;} 
 if (C==1) {IN[Nin][0]  =c[0]; IN[Nin][1]  =c[1]; IN[Nin][2]  =c[2];  ++Nin;} 
      else {OUT[Nout][0]=c[0]; OUT[Nout][1]=c[1]; OUT[Nout][2]=c[2]; ++Nout;} 
 if (D==1) {IN[Nin][0]  =d[0]; IN[Nin][1]  =d[1]; IN[Nin][2]  =d[2];  ++Nin;} 
      else {OUT[Nout][0]=d[0]; OUT[Nout][1]=d[1]; OUT[Nout][2]=d[2]; ++Nout;} 

 
 /* Completely inside or completely outside */
 if ((Nin==0)||(Nin==4)) return(0);

 if (Nin==1){
 /*   i0  
      /|\
     p | r 
    /  q  \
  o0---|---o2 
    \  |  / 
     \ | /
      \|/
       o1   */

  pn=Find_Ver_Map(cur_x,Curr,Post,IN[0],OUT[0]); 
  if(pn==NULL) {pn=Add_Pnt_to_Ver_Map(cur_x,IN[0],OUT[0],Vhead,vox,Curr,Post,thre,MidPntType); } 
  
  qn=Find_Ver_Map(cur_x,Curr,Post,IN[0],OUT[1]); 
  if(qn==NULL) {qn=Add_Pnt_to_Ver_Map(cur_x,IN[0],OUT[1],Vhead,vox,Curr,Post,thre,MidPntType); } 
 
  rn=Find_Ver_Map(cur_x,Curr,Post,IN[0],OUT[2]); 
  if(rn==NULL) {rn=Add_Pnt_to_Ver_Map(cur_x,IN[0],OUT[2],Vhead,vox,Curr,Post,thre,MidPntType);} 


  if ((check_MC_Vertex_neighbors(pn,qn)==0)&&(common_xyz_in_two_vertexes(pn,qn)==1)){Add_MC_Edge_Stack(Ehead,pn,qn);}
  if ((check_MC_Vertex_neighbors(qn,rn)==0)&&(common_xyz_in_two_vertexes(qn,rn)==1)){Add_MC_Edge_Stack(Ehead,qn,rn);}
  if ((check_MC_Vertex_neighbors(rn,pn)==0)&&(common_xyz_in_two_vertexes(rn,pn)==1)){Add_MC_Edge_Stack(Ehead,rn,pn);}

  Add_MC_Face_Stack(Fhead,pn,qn,rn,IN[0]);
 }

 else if (Nin==2){
 /*   o0  
      /|\
     p | \
    /  |  \
  i0--q|---o1 
    \  r  / 
     \ | s
      \|/
       i1   */

   pn=Find_Ver_Map(cur_x,Curr,Post,IN[0],OUT[0]); 
   if(pn==NULL) {pn=Add_Pnt_to_Ver_Map(cur_x,IN[0],OUT[0],Vhead,vox,Curr,Post,thre,MidPntType); }
   qn=Find_Ver_Map(cur_x,Curr,Post,IN[0],OUT[1]); 
   if(qn==NULL) {qn=Add_Pnt_to_Ver_Map(cur_x,IN[0],OUT[1],Vhead,vox,Curr,Post,thre,MidPntType); } 
   rn=Find_Ver_Map(cur_x,Curr,Post,IN[1],OUT[0]); 
   if(rn==NULL) {rn=Add_Pnt_to_Ver_Map(cur_x,IN[1],OUT[0],Vhead,vox,Curr,Post,thre,MidPntType); } 
   sn=Find_Ver_Map(cur_x,Curr,Post,IN[1],OUT[1]); 
   if(sn==NULL) {sn=Add_Pnt_to_Ver_Map(cur_x,IN[1],OUT[1],Vhead,vox,Curr,Post,thre,MidPntType); }
  
   if ((check_MC_Vertex_neighbors(pn,qn)==0)&&(common_xyz_in_two_vertexes(pn,qn)==1)){Add_MC_Edge_Stack(Ehead,pn,qn);}
   if ((check_MC_Vertex_neighbors(qn,rn)==0)&&(common_xyz_in_two_vertexes(qn,rn)==1)){Add_MC_Edge_Stack(Ehead,qn,rn);}
   if ((check_MC_Vertex_neighbors(rn,pn)==0)&&(common_xyz_in_two_vertexes(rn,pn)==1)){Add_MC_Edge_Stack(Ehead,rn,pn);}

   if ((check_MC_Vertex_neighbors(rn,qn)==0)&&(common_xyz_in_two_vertexes(rn,qn)==1)){Add_MC_Edge_Stack(Ehead,rn,qn);}
   if ((check_MC_Vertex_neighbors(qn,sn)==0)&&(common_xyz_in_two_vertexes(qn,sn)==1)){Add_MC_Edge_Stack(Ehead,qn,sn);}
   if ((check_MC_Vertex_neighbors(sn,rn)==0)&&(common_xyz_in_two_vertexes(sn,rn)==1)){Add_MC_Edge_Stack(Ehead,sn,rn);}

   Add_MC_Face_Stack(Fhead,pn,qn,rn,IN[0]);
   Add_MC_Face_Stack(Fhead,qn,rn,sn,IN[0]);
 }

 else if (Nin==3){
 /*   i0  
      /|\
     / | \ 
    /  |  \
   i1--|---i2 
    \  p  / 
     q | r
      \|/
       o0   */
   pn=Find_Ver_Map(cur_x,Curr,Post,IN[0],OUT[0]); 
   if(pn==NULL) {pn=Add_Pnt_to_Ver_Map(cur_x,IN[0],OUT[0],Vhead,vox,Curr,Post,thre,MidPntType); }
   qn=Find_Ver_Map(cur_x,Curr,Post,IN[1],OUT[0]); 
   if(qn==NULL) {qn=Add_Pnt_to_Ver_Map(cur_x,IN[1],OUT[0],Vhead,vox,Curr,Post,thre,MidPntType); } 
   rn=Find_Ver_Map(cur_x,Curr,Post,IN[2],OUT[0]); 
   if(rn==NULL) {rn=Add_Pnt_to_Ver_Map(cur_x,IN[2],OUT[0],Vhead,vox,Curr,Post,thre,MidPntType); } 

   if ((check_MC_Vertex_neighbors(pn,qn)==0)&&(common_xyz_in_two_vertexes(pn,qn)==1)){Add_MC_Edge_Stack(Ehead,pn,qn);}
   if ((check_MC_Vertex_neighbors(qn,rn)==0)&&(common_xyz_in_two_vertexes(qn,rn)==1)){Add_MC_Edge_Stack(Ehead,qn,rn);}
   if ((check_MC_Vertex_neighbors(rn,pn)==0)&&(common_xyz_in_two_vertexes(rn,pn)==1)){Add_MC_Edge_Stack(Ehead,rn,pn);}

   Add_MC_Face_Stack(Fhead,pn,qn,rn,IN[0]);
 }
 return(1);

} /* end of march_each_teterahedron() */ 




void Make_Mid_Point(m,a,b,vox,thre,MidPntType)
 float m[3];     /* mid point (to be calculated) */
 int a[3],b[3];  /* start and end points in lattice */
 struct VOXEL *vox;
 float thre;
 char  MidPntType; /* 'H'alf-point, 'S'imple interpolation, 'C'omplicate Interpolation */
{
 float Va,Vb,Vc,Vd,r,rsimp;
 float ab[3];
 double A,B,C,det,eps; 
 int diff[3],Ndiff;
 int c[3],d[3];

 eps = 0.0001;

 /***********************/
 /*** [1] DECIDE 'r'  ***/
 /***********************/
 
 /** (H) HalfPoint Interpolation  **/
 if (MidPntType == 'H') r = 0.5;
 else{ 
  Va = vox->dat[a[0]][a[1]][a[2]];
  Vb = vox->dat[b[0]][b[1]][b[2]];
  rsimp = (thre - Va)/(Vb - Va); 
 } 
 
 /** (S) Simple Interpolation  **/
 if (MidPntType == 'S'){ r = rsimp; }

 /** (C) Complicated Interpolation  **/
 if (MidPntType == 'C'){

   diff[0] = b[0] - a[0];
   diff[1] = b[1] - a[1];
   diff[2] = b[2] - a[2];
   
   Ndiff = 0;
   if (diff[0]!=0) ++Ndiff;
   if (diff[1]!=0) ++Ndiff;
   if (diff[2]!=0) ++Ndiff;
  
   /* for non-diagonal */
   if (Ndiff==1) r = rsimp; 
   
   /* for diagonal */
   if (Ndiff==2){
    c[0] = a[0]; c[1] = a[1]; c[2] = a[2];
    d[0] = a[0]; d[1] = a[1]; d[2] = a[2];
    
    if (diff[0]==0) { c[1] += diff[1]; d[2] += diff[2]; } 
    if (diff[1]==0) { c[0] += diff[0]; d[2] += diff[2]; } 
    if (diff[2]==0) { c[0] += diff[0]; d[1] += diff[1]; } 
  
    Vc = vox->dat[c[0]][c[1]][c[2]];
    Vd = vox->dat[d[0]][d[1]][d[2]];
  
    A = Va + Vb - Vc - Vd;
    B = Vc - 2 *Va + Vd;
    C = Va - thre;
    
    if (A==0.0) r = rsimp;
    else{
      det = (B*B-4*A*C);
      if (det>0.0) det = sqrt(det);
      r = (-B + det)/2.0/A; 
      if ((r<0.0)||(r>1.0)) r = (-B-det)/2.0/A;
     } 
    }
 
 } /* 'C' */

 /*************************************/
 /*** [2] Decide midpoint using 'r' ***/
 /*************************************/

 if (r<0.0){ 
   if ((r+eps) <0.0) {printf("#ERROR:r %f is outof range\n",r); exit(1);}
   else r = 0.0;
 }

 if (r>1.0){ 
   if ((r-eps) >1.0) {printf("#ERROR:r %f is outof range\n",r); exit(1);}
   else r = 1.0;
 }

 ab[0] = (float)b[0] - (float)a[0];
 ab[1] = (float)b[1] - (float)a[1];
 ab[2] = (float)b[2] - (float)a[2];
 
 m[0] = a[0] + r * ab[0];
 m[1] = a[1] + r * ab[1];
 m[2] = a[2] + r * ab[2];

} /* end of Make_Mid_Point() */









void Malloc_MIDMAP_YZ(M,Y,Z)
 struct MIDMAP_YZ *M;
 int Y,Z;
{
 int y;
 double Mbyte;

 M->Y = Y; M->Z = Z;

 Mbyte = (double)sizeof(struct MIDPNT)*Y*Z /1024.0/1024.0; 
 printf("#Malloc_MIDMAP_YZ(%d %d) %.2lf Mbyte\n",Y,Z,Mbyte);
 
 M->map = (struct MIDPNT **)malloc(sizeof(struct MIDPNT *) * Y); 
 for (y=0;y<Y;++y)
   M->map[y] = (struct MIDPNT *)malloc(sizeof(struct MIDPNT) * Z); 
 
} /* end of Malloc_MIDMAP_YZ() */





void Reset_MIDMAP_YZ(M)
 struct MIDMAP_YZ *M;
{
 int y,z,i,j,k;
 
  for (y=0;y<M->Y;++y){
    for (z=0;z<M->Z;++z){
      for (i=0;i<3;++i){ 
        for (j=0;j<3;++j){
          for (k=0;k<3;++k){ M->map[y][z].ver[i][j][k] = NULL;}
        }
      }
    }
  }
} /* end of  Reset_MIDMAP_YZ() */




void Free_MIDMAP_YZ(M)
 struct MIDMAP_YZ *M;
{
  int y;

 for (y=0;y<M->Y;++y){free(M->map[y]);}
 free(M->map);
 
} /* end of Free_MIDMAP_YZ() */





struct MC_VERTEX  *Add_Pnt_to_Ver_Map(cur_x,I,J,Vhead,vox,Curr,Post,thre,MidPntType)
 int    cur_x;
 int    I[3],J[3];
 struct MC_VERTEX *Vhead;
 struct VOXEL  *vox;
 struct MIDMAP_YZ *Curr,*Post;
 float  thre;
 char   MidPntType;
{
  struct MC_VERTEX *vn;
  float P[3]; 
  int i,d[3];

  Make_Mid_Point(P,I,J,vox,thre,MidPntType); 
  vn = Add_MC_Vertex_Stack(Vhead,P);

  for (i=0;i<3;++i) d[i] = J[i]-I[i]+1;
       if (I[0]==cur_x) Curr->map[I[1]][I[2]].ver[d[0]][d[1]][d[2]] = vn;
  else if (I[0]> cur_x) Post->map[I[1]][I[2]].ver[d[0]][d[1]][d[2]] = vn;
  
  for (i=0;i<3;++i) d[i] = I[i]-J[i]+1;
       if (J[0]==cur_x) Curr->map[J[1]][J[2]].ver[d[0]][d[1]][d[2]] = vn;
  else if (J[0]> cur_x) Post->map[J[1]][J[2]].ver[d[0]][d[1]][d[2]] = vn;

  return(vn);

} /* end of Add_Pnt_to_Ver_Map() */




struct MC_VERTEX *Find_Ver_Map(cur_x,Curr,Post,I,J)
 int cur_x;
 struct MIDMAP_YZ *Curr,*Post;
 int I[3],J[3]; 
{
 int i,d[3];

 for (i=0;i<3;++i) d[i] = J[i]-I[i]+1;
      if (I[0]==cur_x) return(Curr->map[I[1]][I[2]].ver[d[0]][d[1]][d[2]]); 
 else if (I[0]> cur_x) return(Post->map[I[1]][I[2]].ver[d[0]][d[1]][d[2]]); 
 else return(NULL); 

} /* end of Find_Ver_Map() */




struct MC_VERTEX *Add_MC_Vertex_Stack(Vhead,posf)
 struct MC_VERTEX *Vhead;
 float posf[3];
{
 struct MC_VERTEX *newn;
 int i;

 newn = (struct MC_VERTEX*)malloc(sizeof(struct MC_VERTEX));
 newn->pos[0] = posf[0]; 
 newn->pos[1] = posf[1]; 
 newn->pos[2] = posf[2]; 

 newn->next = Vhead->next;

 if (Vhead->next != NULL) { newn->num = newn->next->num + 1; }
 else newn->num = 0;

 Vhead->next = newn;
 for (i=0;i<MAX_MC_VERTEX_NEIGHBOR;++i) newn->neighbors[i] = NULL;

 return(newn);

} /* end of Add_MC_Vertex_Stack() */ 



void Free_MC_Vertex_Stack(Vhead)
 struct MC_VERTEX *Vhead;
{
 struct MC_VERTEX *cur,*next;

 cur = Vhead;
 next = cur->next;
 while (next != NULL){ 
   cur = next;
   next = cur->next;
   free(cur);
 }
 Vhead->next = NULL;

} /* end of Free_MC_Vertex_Stack() */ 





struct MC_FACE *Add_MC_Face_Stack(Fhead,an,bn,cn,in)
 struct MC_FACE *Fhead;
 struct MC_VERTEX *an,*bn,*cn;
 int in[3];       /* Inside Vertex */
{
 struct MC_FACE *newn;
 float norm[3],AB[3],AC[3],cross[3]; 
 int ok,i;

 /* Malloc New Face Node */
 newn = (struct MC_FACE*)malloc(sizeof(struct MC_FACE));

 /*** Set Vertex a,b,c by proper order ***/
 Sub_Vec_FI(norm,an->pos,in);
 Sub_Vec_F(AB,bn->pos,an->pos);
 Sub_Vec_F(AC,cn->pos,an->pos);
 Cross_Prod_F(cross,AB,AC);
 if (Dot_Prod_F(norm,cross)>0.0) { newn->a = an; newn->b = bn; newn->c = cn;}
 else { newn->a = an; newn->b = cn; newn->c = bn;}

 i = 0; ok = 0;
 while ((i<MAX_MC_VERTEX_NEIGHBOR)&&(ok==0)){
   if (an->neighbors[i]==NULL) {an->neighbors[i] = bn; ok = 1;}
   i += 1; } 

 i = 0; ok = 0;
 while ((i<MAX_MC_VERTEX_NEIGHBOR)&&(ok==0)){
   if (bn->neighbors[i]==NULL) {bn->neighbors[i] = an; ok = 1;}
   i += 1; } 


 i = 0; ok = 0;
 while ((i<MAX_MC_VERTEX_NEIGHBOR)&&(ok==0)){
   if (bn->neighbors[i]==NULL) {bn->neighbors[i] = cn; ok = 1;}
   i += 1; } 

 i = 0; ok = 0;
 while ((i<MAX_MC_VERTEX_NEIGHBOR)&&(ok==0)){
   if (cn->neighbors[i]==NULL) {cn->neighbors[i] = bn; ok = 1;}
   i += 1; } 

 i = 0; ok = 0;
 while ((i<MAX_MC_VERTEX_NEIGHBOR)&&(ok==0)){
   if (cn->neighbors[i]==NULL) {cn->neighbors[i] = an; ok = 1;}
   i += 1; } 

 i = 0; ok = 0;
 while ((i<MAX_MC_VERTEX_NEIGHBOR)&&(ok==0)){
   if (an->neighbors[i]==NULL) {an->neighbors[i] = cn; ok = 1;}
   i += 1; } 

 newn->next = Fhead->next;
 
 if (Fhead->next != NULL) { newn->num = newn->next->num + 1; }
 else newn->num = 0;

 Fhead->next = newn;
  
 return(newn);

} /* end of Add_MC_Face_Stack() */ 



struct MC_FACE *Add_MC_Face_Stack_Simple(Fhead,an,bn,cn)
 struct MC_FACE *Fhead;
 struct MC_VERTEX *an,*bn,*cn;
{
 struct MC_FACE *newn;

 /* Malloc New Face Node */
 newn = (struct MC_FACE*)malloc(sizeof(struct MC_FACE));

 newn->a = an;
 newn->b = bn;
 newn->c = cn;

 newn->next = Fhead->next;

 if (Fhead->next != NULL) { newn->num = newn->next->num + 1; }
 else newn->num = 0;

 Fhead->next = newn;
 /* printf("#Add_MC_Face_Stack_Simple(%d)\n",newn->num); */
  
 return(newn);

} /* end of Add_MC_Face_Stack_Simple() */ 







void Free_MC_Face_Stack(Fhead)
 struct MC_FACE *Fhead;
{
 struct MC_FACE *cur,*next;

 cur = Fhead;
 next = cur->next;

 while (next != NULL){ 
   cur = next;
   next = cur->next;
   free(cur);
}
 Fhead->next = NULL;

} /* end of Free_MC_Face_Stack() */ 



struct MC_EDGE *Add_MC_Edge_Stack(Ehead,an,bn)
 struct MC_EDGE *Ehead;
 struct MC_VERTEX *an,*bn;
{
 struct MC_EDGE *newn;


 /* Malloc New Edge Node */
 newn = (struct MC_EDGE*)malloc(sizeof(struct MC_EDGE));
 if (an->num < bn->num){
   newn->a = an; newn->b = bn;
 } 
 else{
   newn->a = bn; newn->b = an;
 }
 newn->next = Ehead->next;
 Ehead->next = newn;
 return(newn);
} /* end of Add_MC_Edge_Stack() */ 


int  check_MC_Vertex_neighbors(an,bn)
 struct MC_VERTEX *an,*bn;
{
 int i;

 if (an->neighbors[0] == NULL) return(0);

 for (i=0;i<MAX_MC_VERTEX_NEIGHBOR;++i){
   if ((an->neighbors[i] != NULL) && (an->neighbors[i]->num == bn->num)) return(1);
 }
 return(0);

} /* end of check_MC_Vertex_neighbors() */






void Free_MC_Edge_Stack(Ehead)
 struct MC_EDGE *Ehead;
{
 struct MC_EDGE *cur,*next;

 cur = Ehead;
 next = cur->next;

 while (next != NULL){ 
   cur = next;
   next = cur->next;
   free(cur);
 }
 Ehead->next = NULL;

} /* end of Free_MC_Edge_Stack() */ 



int Check_Overlap_MC_Edge(Ehead,an,bn)
 struct MC_EDGE *Ehead;
 struct MC_VERTEX *an,*bn;
{
 struct MC_EDGE *e;
 int anumQ, bnumQ;
 
 if (an->num < bn->num){ anumQ = an->num; bnumQ = bn->num;}
                  else { anumQ = bn->num; bnumQ = an->num;}
 e = Ehead;
 while (e->next != NULL){ 
   e = e->next;
   if ((e->a->num == anumQ) && (e->b->num == bnumQ)) return(1);
 }
 return(0);
} /* end of Check_Overlap_MC_Edge() */ 




int common_xyz_in_two_vertexes(an,bn)
 struct MC_VERTEX *an,*bn;
{
  if ((an->pos[0]==bn->pos[0])||(an->pos[1]==bn->pos[1])||(an->pos[2]==bn->pos[2])) return(1);
  else return(0);
} /* end of Check_Overlap_MC_Edge() */ 







void Sub_Vec_F(c,a,b)
 float c[3],a[3],b[3];
{ c[0] = a[0]-b[0];
  c[1] = a[1]-b[1];
  c[2] = a[2]-b[2]; }

void Sub_Vec_FI(c,a,b)
 float c[3],a[3];
 int b[3];
{ c[0] = a[0]-(float)b[0];
  c[1] = a[1]-(float)b[1];
  c[2] = a[2]-(float)b[2]; }


void Cross_Prod_F(c,a,b)
 float c[3],a[3],b[3];
{ c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0]; }

float Dot_Prod_F(a,b)
 float a[3],b[3];
{ return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]); }





void Malloc_CMATRIX(M,Nmalloc)
 struct CMATRIX *M;
 int Nmalloc;
{
 int i;

 M->N = Nmalloc;
 M->m = (char **)malloc(sizeof(char*)*M->N);
 for (i=0;i<M->N;++i) M->m[i] = (char *)malloc(sizeof(char)*M->N);
} /* end of Malloc_CMATRIX() */


void Free_CMATRIX(M)
 struct CMATRIX *M;
{
 int i;
 for (i=0;i<M->N;++i) free(M->m[i]);
 free(M->m);

} /* end of Free_CMATRIX() */



int Number_Of_MC_VERTEX(Vhead)
 struct MC_VERTEX *Vhead;
{ int Nver;
  struct MC_VERTEX *v;
 Nver = 0;
 v = Vhead;
 while (v->next != NULL) {v = v->next; ++Nver;}
 return(Nver);
} /* end of Number_Of_MC_VERTEX() */


int Number_Of_MC_FACE(Fhead)
 struct MC_FACE *Fhead;
{ int Nface,Nerror,a,b,c;
  struct MC_FACE *f;
 Nface = 0;
 Nerror = 0;
 f = Fhead;
 while (f->next != NULL){
   f = f->next; 
   if ((f->a==NULL)||(f->b==NULL) || (f->c==NULL)){ Nerror += 1; }
   a = f->a->num;
   b = f->b->num;
   c = f->c->num;
   /* printf("%d %d %d\n",a,b,c); */
   ++Nface;
 }
 if (Nerror>0){printf("#Nerror_face %d\n",Nerror);}

 return(Nface);
} /* end of Number_Of_MC_FACE() */

