/*
                                                                                
 <GridMap.c>
                                                                                
  functions for dealing with  struct CHAR3DMAP.
                                                                                
*/
                                                                                
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "globalvar.h"
#include "GridMap.h"

/*** FUNCIONS (GLOBAL) ***/
int Malloc_CHAR2DMAP();
void Initialize_CHAR2DMAP();
void Free_CHAR2DMAP();

int Malloc_CHAR3DMAP();
void Initialize_CHAR3DMAP();
void Free_CHAR3DMAP();

int Malloc_CHAR4DMAP();
void Initialize_CHAR4DMAP();
void Free_CHAR4DMAP();

int  Malloc_FLOAT1DMAP();
void Initialize_FLOAT1DMAP();
void Free_FLOAT1DMAP();

int  Malloc_FLOAT2DMAP();
void Initialize_FLOAT2DMAP();
void Free_FLOAT2DMAP();

int  Malloc_FLOAT3DMAP();
void Initialize_FLOAT3DMAP();
void Free_FLOAT3DMAP();

int Malloc_FLOAT4DMAP();
void Initialize_FLOAT4DMAP();
void Free_FLOAT4DMAP();

int  Malloc_DOUBLE1DMAP();
void Initialize_DOUBLE1DMAP();
void Free_DOUBLE1DMAP();

int  Malloc_DOUBLE2DMAP();
void Initialize_DOUBLE2DMAP();
void Free_DOUBLE2DMAP();

int  Malloc_DOUBLE3DMAP();
void Initialize_DOUBLE3DMAP();
void Free_DOUBLE3DMAP();

int Malloc_DOUBLE4DMAP();
void Initialize_DOUBLE4DMAP();
void Free_DOUBLE4DMAP();


/*****************/
/*** FUNCTIONS ***/
/*****************/ 
int Malloc_CHAR2DMAP(M,x,y)
 struct CHAR2DMAP *M;
 int x,y;
{
 int i;
 double byte,Mbyte;
                                                                                          
 byte = (double)x*y*(double)sizeof(char);
 Mbyte = byte/1024/1024;
 printf("#Malloc_CHAR2DMAP(%d %d) %.2lf byte -> %.2lf Mbyte\n",x,y,byte,Mbyte);
 if (Mbyte > MAX_MEMORY_MEGA) { 
   printf("#ERROR:Memory size (%.2lf MB) is larger than MAX_MEMORY_MEGA(%d MB).\n",Mbyte,MAX_MEMORY_MEGA); 
 exit(1); }
                                                                                          
 M->N[0] = x; M->N[1] = y;
                                                                                          
 M->dat = (unsigned char **)malloc(sizeof(unsigned char *) * x);
 for (i=0;i<x;++i){
  M->dat[i] = (unsigned char *)malloc(sizeof(unsigned char) * y);
 }
 
 Initialize_CHAR2DMAP(M);
 return(1);
} /* end of Malloc_CHAR2DMAP() */


void Initialize_CHAR2DMAP(M)
 struct CHAR2DMAP *M;
{
 int i,j;
 for (i=0;i<M->N[0];++i){
   for (j=0;j<M->N[1];++j){M->dat[i][j] = 0;}
 }
} /* end of Initialize_CHAR2DMAP() */
                                                                                          
                                                                                          
void Free_CHAR2DMAP(M)
 struct CHAR2DMAP *M;
{
 int i;
 printf("#Free_CHAR2DMAP(%d %d)\n",M->N[0],M->N[1]);
 for (i=M->N[0]-1;i>=0;--i) {free(M->dat[i]);}
 free(M->dat);
} /* end of Free_CHAR2DMAP() */



int Malloc_CHAR3DMAP(M,x,y,z)
 struct CHAR3DMAP *M;
 int x,y,z;
{
 int i,j;
 double byte,Mbyte;
                                                                                          
 byte = (double)x*y*z*(double)sizeof(char);
 Mbyte = byte/1024/1024;
 printf("#Malloc_CHAR3DMAP(%d %d %d) %.2lf byte -> %.2lf Mbyte\n",x,y,z,byte,Mbyte);
 if (Mbyte > MAX_MEMORY_MEGA) { 
   printf("#ERROR:Memory size (%.2lf MB) is larger than MAX_MEMORY_MEGA(%d MB).\n",Mbyte,MAX_MEMORY_MEGA); 
 exit(1); }
                                                                                          
 M->N[0] = x; M->N[1] = y; M->N[2] = z;
                                                                                          
 M->dat = (unsigned char ***)malloc(sizeof(unsigned char **) * x);
 for (i=0;i<x;++i){
   M->dat[i] = (unsigned char **)malloc(sizeof(unsigned char *) * y);
   for (j=0;j<y;++j){
     M->dat[i][j] = (unsigned char *)malloc(sizeof(unsigned char) * z);
   }
 }
                                                                                          
  Initialize_CHAR3DMAP(M);
                                                                                          
  return(1);
} /* end of Malloc_CHAR3DMAP() */


void Initialize_CHAR3DMAP(M)
 struct CHAR3DMAP *M;
{
 int i,j,k;
 for (i=0;i<M->N[0];++i)
  for (j=0;j<M->N[1];++j)
   for (k=0;k<M->N[2];++k) M->dat[i][j][k] = 0;
} /* end of Initialize_CHAR3DMAP() */
                                                                                          
                                                                                          
void Free_CHAR3DMAP(M)
 struct CHAR3DMAP *M;
{
 int i,j;
 printf("#Free_CHAR3DMAP(%d %d %d)\n",M->N[0],M->N[1],M->N[2]);
 for (i=M->N[0]-1;i>=0;--i){
  for (j=M->N[1]-1;j>=0;--j) free(M->dat[i][j]);
  free(M->dat[i]);
 }
 free(M->dat);
} /* end of Free_CHAR3DMAP() */




int Malloc_CHAR4DMAP(M,x,y,z,w)
 struct CHAR4DMAP *M;
 int x,y,z,w;
{
 int i,j,k;
 double byte,Mbyte;
                                                                                          
 byte = (double)x*y*z*(double)sizeof(char);
 Mbyte = byte/1024/1024;
 printf("#Malloc_CHAR4DMAP(%d %d %d %d) %.2lf byte -> %.2lf Mbyte\n",x,y,z,w,byte,Mbyte);
 if (Mbyte > MAX_MEMORY_MEGA) { 
   printf("#ERROR:Memory size (%.2lf MB) is larger than MAX_MEMORY_MEGA(%d MB).\n",Mbyte,MAX_MEMORY_MEGA); 
 exit(1); }
                                                                                          
 M->N[0] = x; M->N[1] = y; M->N[2] = z; M->N[3] = w;
                                                                                          
 M->dat = (unsigned char ****)malloc(sizeof(unsigned char ***) * x);
 for (i=0;i<x;++i){
  M->dat[i] = (unsigned char ***)malloc(sizeof(unsigned char **) * y);
  for (j=0;j<y;++j){
   M->dat[i][j] = (unsigned char **)malloc(sizeof(unsigned char*) * z);
   for (k=0;k<z;++k){
    M->dat[i][j][k] = (unsigned char *)malloc(sizeof(unsigned char) * w);
   } 
  } 
 }
  
  Initialize_CHAR4DMAP(M);
  return(1);
} /* end of Malloc_CHAR4DMAP() */



void Initialize_CHAR4DMAP(M)
 struct CHAR4DMAP *M;
{
 int i,j,k,m;
 for (i=0;i<M->N[0];++i){
  for (j=0;j<M->N[1];++j){
   for (k=0;k<M->N[2];++k){ 
    for (m=0;m<M->N[3];++m){M->dat[i][j][k][m] = 0;}
   }
  }
 }
} /* end of Initialize_CHAR4DMAP() */
                                                                                          
                                                                                          
void Free_CHAR4DMAP(M)
 struct CHAR4DMAP *M;
{
 int i,j,k;
 printf("#Free_CHAR4DMAP(%d %d %d %d)\n",M->N[0],M->N[1],M->N[2],M->N[3]);
 for (i=M->N[0]-1;i>=0;--i) {
  for (j=M->N[1]-1;j>=0;--j){ 
   for (k=M->N[2]-1;k>=0;--k){ free(M->dat[i][j][k]);}
   free(M->dat[i][j]);
  }
  free(M->dat[i]);
 }
 free(M->dat);
} /* end of Free_CHAR4DMAP() */


int Malloc_FLOAT1DMAP(M,N)
 struct FLOAT1DMAP *M;
 int N;
{
 double byte,Mbyte;
 byte = (double)N*(double)sizeof(float);
 Mbyte = byte/1024/1024;
 printf("#Malloc_FLOAT1DMAP(%d) %.2lf byte -> %.2lf Mbyte\n",N,byte,Mbyte);
 if (Mbyte > MAX_MEMORY_MEGA) { 
   printf("#ERROR:Memory size (%.2lf MB) is larger than MAX_MEMORY_MEGA(%d MB).\n",Mbyte,MAX_MEMORY_MEGA); 
 exit(1); }
 M->N = N; 
 M->dat = (float *)malloc(sizeof(float) * N);
 Initialize_FLOAT1DMAP(M);
 return(1);
} /* end of Malloc_FLOAT1DMAP() */



void Initialize_FLOAT1DMAP(M)
 struct FLOAT1DMAP *M;
{
 int i;
 for (i=0;i<M->N;++i) M->dat[i] = 0.0;
} /* end of Initialize_FLOAT1DMAP() */
                                                                                          
                                                                                          
void Free_FLOAT1DMAP(M)
 struct FLOAT1DMAP *M;
{
 printf("#Free_FLOAT1DMAP(N %d)\n",M->N);
 free(M->dat);
} /* end of Free_FLOAT1DMAP() */


int Malloc_FLOAT2DMAP(M,x,y)
 struct FLOAT2DMAP *M;
 int x,y;
{
 int i;
 double byte,Mbyte;
                                                                                          
 byte = (double)x*y*(double)sizeof(float);
 Mbyte = byte/1024/1024;
 printf("#Malloc_FLOAT2DMAP(%d %d) %.2lf byte -> %.2lf Mbyte\n",x,y,byte,Mbyte);
 M->N[0] = x; M->N[1] = y; 

 M->dat = (float **)malloc(sizeof(float *) * M->N[0]);
 for (i=0;i<M->N[0];++i){ M->dat[i] = (float *)malloc(sizeof(float) * M->N[1]); }
 Initialize_FLOAT2DMAP(M);
 return(1);
} /* end of Malloc_FLOAT2DMAP() */



void Initialize_FLOAT2DMAP(M)
 struct FLOAT2DMAP *M;
{
 int i,j;
 for (i=0;i<M->N[0];++i){
  for (j=0;j<M->N[1];++j){M->dat[i][j] = 0.0;}
 }
} /* end of Initialize_FLOAT2DMAP() */
                                                                                          
                                                                                          
void Free_FLOAT2DMAP(M)
 struct FLOAT2DMAP *M;
{
 int i;
 printf("#Free_FLOAT2DMAP(%d %d)\n",M->N[0],M->N[1]);
 for (i=M->N[0]-1;i>=0;--i) { free(M->dat[i]); }
 free(M->dat);
} /* end of Free_FLOAT2DMAP() */



int Malloc_FLOAT3DMAP(M,x,y,z)
 struct FLOAT3DMAP *M;
 int x,y,z;
{
 int i,j;
 double byte,Mbyte;
                                                                                          
 byte = (double)x*y*z*(double)sizeof(float);
 Mbyte = byte/1024/1024;
 printf("#Malloc_FLOAT3DMAP(%d %d %d) %.2lf byte -> %.2lf Mbyte\n",x,y,z,byte,Mbyte);
 if (Mbyte > MAX_MEMORY_MEGA) { 
   printf("#ERROR:Memory size (%.2lf MB) is larger than MAX_MEMORY_MEGA(%d MB).\n",Mbyte,MAX_MEMORY_MEGA); 
 exit(1); }
                                                                                          
 M->N[0] = x; M->N[1] = y; M->N[2] = z; 
 
 M->dat = (float ***)malloc(sizeof(float **) * x);
 for (i=0;i<x;++i){
  M->dat[i] = (float **)malloc(sizeof(float *) * y);
  for (j=0;j<y;++j){
   M->dat[i][j] = (float *)malloc(sizeof(float) * z);
  } 
 }
  
 Initialize_FLOAT3DMAP(M);
 return(1);
} /* end of Malloc_FLOAT3DMAP() */



void Initialize_FLOAT3DMAP(M)
 struct FLOAT3DMAP *M;
{
 int i,j,k;
 for (i=0;i<M->N[0];++i){
  for (j=0;j<M->N[1];++j){
   for (k=0;k<M->N[2];++k){M->dat[i][j][k] = 0.0;}
  }
 }
} /* end of Initialize_FLOAT3DMAP() */
                                                                                          
                                                                                          
void Free_FLOAT3DMAP(M)
 struct FLOAT3DMAP *M;
{
 int i,j;
 printf("#Free_FLOAT3DMAP(%d %d %d)\n",M->N[0],M->N[1],M->N[2]);
 for (i=M->N[0]-1;i>=0;--i) {
  for (j=M->N[1]-1;j>=0;--j){ free(M->dat[i][j]);}
 free(M->dat[i]);
 }
 free(M->dat);
} /* end of Free_FLOAT3DMAP() */








int Malloc_FLOAT4DMAP(M,x,y,z,w)
 struct FLOAT4DMAP *M;
 int x,y,z,w;
{
 int i,j,k;
 double byte,Mbyte;
                                                                                          
 byte = (double)x*y*z*w*(double)sizeof(float);
 Mbyte = byte/1024/1024;
 printf("#Malloc_FLOAT4DMAP(%d %d %d %d) %.2lf byte -> %.2lf Mbyte\n",x,y,z,w,byte,Mbyte);
 if (Mbyte > MAX_MEMORY_MEGA) { 
   printf("#ERROR:Memory size (%.2lf MB) is larger than MAX_MEMORY_MEGA(%d MB).\n",Mbyte,MAX_MEMORY_MEGA); 
 exit(1); }
                                                                                          
 M->N[0] = x; M->N[1] = y; M->N[2] = z; M->N[3] = w;
                                                                                          
 M->dat = (float ****)malloc(sizeof(float ***) * x);
 for (i=0;i<x;++i){
  M->dat[i] = (float ***)malloc(sizeof(float **) * y);
  for (j=0;j<y;++j){
   M->dat[i][j] = (float **)malloc(sizeof(float*) * z);
   for (k=0;k<z;++k){
    M->dat[i][j][k] = (float *)malloc(sizeof(float) * w);
   } 
  } 
 }
  
  Initialize_FLOAT4DMAP(M);
  return(1);
} /* end of Malloc_FLOAT4DMAP() */



void Initialize_FLOAT4DMAP(M)
 struct FLOAT4DMAP *M;
{
 int i,j,k,m;
 for (i=0;i<M->N[0];++i)
  for (j=0;j<M->N[1];++j)
   for (k=0;k<M->N[2];++k) 
    for (m=0;m<M->N[3];++m) M->dat[i][j][k][m] = 0.0;

} /* end of Initialize_FLOAT4DMAP() */
                                                                                          
                                                                                          
void Free_FLOAT4DMAP(M)
 struct FLOAT4DMAP *M;
{
 int i,j,k;
 printf("#Free_FLOAT4DMAP(%d %d %d %d)\n",M->N[0],M->N[1],M->N[2],M->N[3]);
 for (i=M->N[0]-1;i>=0;--i) {
  for (j=M->N[1]-1;j>=0;--j){ 
   for (k=M->N[2]-1;k>=0;--k){ free(M->dat[i][j][k]);}
   free(M->dat[i][j]);
  }
  free(M->dat[i]);
 }
 free(M->dat);
} /* end of Free_FLOAT4DMAP() */


int Malloc_DOUBLE1DMAP(M,N)
 struct DOUBLE1DMAP *M;
 int N;
{
 double byte,Mbyte;
 byte = (double)N*(double)sizeof(double);
 Mbyte = byte/1024/1024;
 printf("#Malloc_DOUBLE1DMAP(%d) %.2lf byte -> %.2lf Mbyte\n",N,byte,Mbyte);
 if (Mbyte > MAX_MEMORY_MEGA) { 
   printf("#ERROR:Memory size (%.2lf MB) is larger than MAX_MEMORY_MEGA(%d MB).\n",Mbyte,MAX_MEMORY_MEGA); 
 exit(1); }
 M->N = N; 
 M->dat = (double *)malloc(sizeof(double) * N);
 Initialize_DOUBLE1DMAP(M);
 return(1);
} /* end of Malloc_DOUBLE1DMAP() */



void Initialize_DOUBLE1DMAP(M)
 struct DOUBLE1DMAP *M;
{
 int i;
 for (i=0;i<M->N;++i) M->dat[i] = 0.0;
} /* end of Initialize_DOUBLE1DMAP() */
                                                                                          
                                                                                          
void Free_DOUBLE1DMAP(M)
 struct DOUBLE1DMAP *M;
{
 printf("#Free_DOUBLE1DMAP(N %d)\n",M->N);
 free(M->dat);
} /* end of Free_DOUBLE1DMAP() */


int Malloc_DOUBLE2DMAP(M,x,y)
 struct DOUBLE2DMAP *M;
 int x,y;
{
 int i;
 double byte,Mbyte;
                                                                                          
 byte = (double)x*y*(double)sizeof(double);
 Mbyte = byte/1024/1024;
 printf("#Malloc_DOUBLE2DMAP(%d %d) %.2lf byte -> %.2lf Mbyte\n",x,y,byte,Mbyte);
 if (Mbyte > MAX_MEMORY_MEGA) { 
   printf("#ERROR:Memory size (%.2lf MB) is larger than MAX_MEMORY_MEGA(%d MB).\n",Mbyte,MAX_MEMORY_MEGA); 
 exit(1); }
 M->N[0] = x; M->N[1] = y; 

 M->dat = (double **)malloc(sizeof(double *) * x);
 for (i=0;i<x;++i){ M->dat[i] = (double *)malloc(sizeof(double) * y); }
 Initialize_DOUBLE2DMAP(M);
 return(1);
} /* end of Malloc_DOUBLE2DMAP() */



void Initialize_DOUBLE2DMAP(M)
 struct DOUBLE2DMAP *M;
{
 int i,j;
 for (i=0;i<M->N[0];++i){
  for (j=0;j<M->N[1];++j){M->dat[i][j] = 0.0;}
 }
} /* end of Initialize_DOUBLE2DMAP() */
                                                                                          
                                                                                          
void Free_DOUBLE2DMAP(M)
 struct DOUBLE2DMAP *M;
{
 int i;
 printf("#Free_DOUBLE2DMAP(%d %d)\n",M->N[0],M->N[1]);
 for (i=M->N[0]-1;i>=0;--i) { free(M->dat[i]); }
 free(M->dat);
} /* end of Free_DOUBLE2DMAP() */



int Malloc_DOUBLE3DMAP(M,x,y,z)
 struct DOUBLE3DMAP *M;
 int x,y,z;
{
 int i,j;
 double byte,Mbyte;
                                                                                          
 byte = (double)x*y*z*(double)sizeof(double);
 Mbyte = byte/1024/1024;
 printf("#Malloc_DOUBLE3DMAP(%d %d %d) %.2lf byte -> %.2lf Mbyte\n",x,y,z,byte,Mbyte);
 if (Mbyte > MAX_MEMORY_MEGA) { 
   printf("#ERROR:Memory size (%.2lf MB) is larger than MAX_MEMORY_MEGA(%d MB).\n",Mbyte,MAX_MEMORY_MEGA); 
 exit(1); }
                                                                                          
 M->N[0] = x; M->N[1] = y; M->N[2] = z; 
 
 M->dat = (double ***)malloc(sizeof(double **) * x);
 for (i=0;i<x;++i){
  M->dat[i] = (double **)malloc(sizeof(double *) * y);
  for (j=0;j<y;++j){
   M->dat[i][j] = (double *)malloc(sizeof(double) * z);
  } 
 }
  
 Initialize_DOUBLE3DMAP(M);
 return(1);
} /* end of Malloc_DOUBLE3DMAP() */



void Initialize_DOUBLE3DMAP(M)
 struct DOUBLE3DMAP *M;
{
 int i,j,k;
 for (i=0;i<M->N[0];++i){
  for (j=0;j<M->N[1];++j){
   for (k=0;k<M->N[2];++k){ M->dat[i][j][k] = 0.0;}
  }
 }
} /* end of Initialize_DOUBLE3DMAP() */
                                                                                          
                                                                                          
void Free_DOUBLE3DMAP(M)
 struct DOUBLE3DMAP *M;
{
 int i,j;
 printf("#Free_DOUBLE3DMAP(%d %d %d)\n",M->N[0],M->N[1],M->N[2]);
 for (i=M->N[0]-1;i>=0;--i) {
  for (j=M->N[1]-1;j>=0;--j){ free(M->dat[i][j]);}
 free(M->dat[i]);
 }
 free(M->dat);
} /* end of Free_DOUBLE3DMAP() */








int Malloc_DOUBLE4DMAP(M,x,y,z,w)
 struct DOUBLE4DMAP *M;
 int x,y,z,w;
{
 int i,j,k;
 double byte,Mbyte;
                                                                                          
 byte = (double)x*y*z*w*(double)sizeof(double);
 Mbyte = byte/1024/1024;
 printf("#Malloc_DOUBLE4DMAP(%d %d %d %d) %.2lf byte -> %.2lf Mbyte\n",x,y,z,w,byte,Mbyte);
 if (Mbyte > MAX_MEMORY_MEGA) { 
   printf("#ERROR:Memory size (%.2lf MB) is larger than MAX_MEMORY_MEGA(%d MB).\n",Mbyte,MAX_MEMORY_MEGA); 
 exit(1); }
                                                                                          
 M->N[0] = x; M->N[1] = y; M->N[2] = z; M->N[3] = w;
                                                                                          
 M->dat = (double ****)malloc(sizeof(double ***) * x);
 for (i=0;i<x;++i){
  M->dat[i] = (double ***)malloc(sizeof(double **) * y);
  for (j=0;j<y;++j){
   M->dat[i][j] = (double **)malloc(sizeof(double*) * z);
   for (k=0;k<z;++k){
    M->dat[i][j][k] = (double *)malloc(sizeof(double) * w);
   } 
  } 
 }
  
  Initialize_DOUBLE4DMAP(M);
  return(1);
} /* end of Malloc_DOUBLE4DMAP() */



void Initialize_DOUBLE4DMAP(M)
 struct DOUBLE4DMAP *M;
{
 int i,j,k,m;
 for (i=0;i<M->N[0];++i){
  for (j=0;j<M->N[1];++j){
   for (k=0;k<M->N[2];++k){ 
    for (m=0;m<M->N[3];++m){ M->dat[i][j][k][m] = 0.0;}
   }
  } 
 }

} /* end of Initialize_DOUBLE4DMAP() */
                                                                                          
                                                                                          
void Free_DOUBLE4DMAP(M)
 struct DOUBLE4DMAP *M;
{
 int i,j,k;
 printf("#Free_DOUBLE4DMAP(%d %d %d %d)\n",M->N[0],M->N[1],M->N[2],M->N[3]);
 for (i=M->N[0]-1;i>=0;--i) {
  for (j=M->N[1]-1;j>=0;--j){ 
   for (k=M->N[2]-1;k>=0;--k){ free(M->dat[i][j][k]);}
   free(M->dat[i][j]);
  }
  free(M->dat[i]);
 }
 free(M->dat);
} /* end of Free_DOUBLE4DMAP() */


