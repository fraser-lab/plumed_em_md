/*

 <Voxel.c> 

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include "globalvar.h"
#include "Voxel.h" 
#include "mc_verface.h" 
#include "BasicIO.h" 
#include "HeapSort.h" 

#define MAX_LINE 10000

/** Functions(GLOBAL) **/
int  Malloc_Voxel();
void Read_Voxel_File();
void Write_Voxel_File();
void Write_Voxel_File_In_PDB();
void Free_Voxel();
void Initialize_Voxel();
void Copy_Voxel();
void Voxel_Statistics();
float Find_Threshold_Density_For_NumVoxel();
void Output_Voxel_Histgram();
void Make_Reduced_Size_Voxel();

/** Functions(LOCAL) **/



int Malloc_Voxel(vox,x,y,z)
 struct VOXEL *vox;
 int x,y,z;
{
 int i,j;
 double byte,Mbyte;

 byte = (double)x*y*z*(double)sizeof(float); 
 Mbyte = byte/1024/1024;
/*
 printf("#Malloc_Voxel(%d %d %d) %.2lf byte %.2lf Mbyte\n",x,y,z,byte,Mbyte);
*/
 if (Mbyte > MAX_MEMORY_MEGA){
   printf("#ERROR(Malloc_Voxel):Memory size (%d %d %d) is too large (%.2lf Mbyte > MAX_MEMORY_MEGA %d Mbyte).\n",x,y,z,Mbyte,MAX_MEMORY_MEGA);
   return(0); 
 }

 vox->dat = (float ***)malloc(sizeof(float **) * x);
 for (i=0;i<x;++i){
   vox->dat[i] = (float **)malloc(sizeof(float *) * y);
   for (j=0;j<y;++j){ 
     vox->dat[i][j] = (float *)malloc(sizeof(float) * z);
   }
 }

 vox->mal = 1;
 vox->N[0] = x; vox->N[1] = y; vox->N[2] = z;
 return(1);
} /* end of Malloc_Voxel() */




void Free_Voxel(vox)
 struct VOXEL *vox;
{
 int i,j;
 
 for (i=0;i<vox->N[0];++i){
   for (j=0;j<vox->N[1];++j){free(vox->dat[i][j]);}
   free(vox->dat[i]);
 }

 free(vox->dat);
 vox->mal = 0;

} /* end of Free_Voxel() */


void Initialize_Voxel(vox)
  struct VOXEL *vox;
{
  int i,j,k;
  for (i=0;i<vox->N[0];++i){
    for (j=0;j<vox->N[1];++j){
      for (k=0;k<vox->N[2];++k){ vox->dat[i][j][k] = 0.0; }
    }
  }
} /* end of Initialize_Voxel() */



void Copy_Voxel(nV,oV)
  struct VOXEL *nV,*oV;  /* nV := oV */ 
{
  int x,y,z,k;
  
  for (k=0;k<3;++k){ nV->N[k] = oV->N[k]; }
  
  nV->grid_width = oV->grid_width; 

  for (x=0;x<oV->N[0];++x){
    for (y=0;y<oV->N[1];++y){
      for (z=0;z<oV->N[2];++z){
        nV->dat[x][y][z] = oV->dat[x][y][z];
      }
    }
  }

} /* end of Copy_Voxel() */



void Read_Voxel_File(file,vox)
  char *file;
  struct VOXEL *vox;
{
 FILE *fp;
 char ln[MAX_LINE],buff1[64],buff2[64];
 int Len,K,x,y,z;
 float val;
                                                                                                                  
 /** [FILE FORMAT]
 Values of grid appear  by x-y-z forloop orders.
  [ for (x=0;x<X;++x)
     for (y=0;y<Y;++y)
      for (z=0;z<Z;++z) ]
                                                                                                                  
 [FILE EXAMPLE]
 GRIDLEN 0.500000
 X 23
 Y 28
 Z 22
 MINX 4.750000
 MINY 1.500000
 MINZ 4.750000
 0.000742
 0.000986
 0.001272
 0.001597
 0.001948
 0.002310
 :
 ***/
 printf("#Read_Voxel_File \"%s\"\n",file);
 fp = fopen(file,"r");
 if (fp==NULL) {printf("#ERROR:Can't open file \"%s\"\n",file); exit(1);}
                                                                                                                  
 vox->mal = 0;
 vox->N[0] = vox->N[1] = vox->N[2] = 0;
 vox->OrigPos[0] = vox->OrigPos[1] = vox->OrigPos[2] = 0.0;
 vox->grid_width = 1.0;
 x = y = z = 0;
 K = 0;
 while (feof(fp)==0){
   ln[0] = '\0';
   fgets(ln,MAX_LINE-1,fp);
   /* printf("%s",ln); */
   Len = strlen(ln);
   if (ln[Len-1]=='\n') {ln[Len-1] = '\0'; Len = strlen(ln);}
                                                                                                                  
   if ((ln[0]!='#')&&(Len>0)){
         if (ln[0]=='X')
     { sscanf(ln,"%s %s",buff1,buff2); vox->N[0] = atoi(buff2); }
    else if (ln[0]=='Y')
     { sscanf(ln,"%s %s",buff1,buff2); vox->N[1] = atoi(buff2); }
    else if (ln[0]=='Z')
     { sscanf(ln,"%s %s",buff1,buff2); vox->N[2] = atoi(buff2); }
    else if (strncmp(ln,"GRIDLEN",7)==0)
     { sscanf(ln,"%s %s",buff1,buff2); vox->grid_width = atof(buff2); }
    else if (strncmp(ln,"MINX",4)==0)
     { sscanf(ln,"%s %s",buff1,buff2); vox->OrigPos[0] = atof(buff2); }
    else if (strncmp(ln,"MINY",4)==0)
     { sscanf(ln,"%s %s",buff1,buff2); vox->OrigPos[1] = atof(buff2); }
    else if (strncmp(ln,"MINZ",4)==0)
     { sscanf(ln,"%s %s",buff1,buff2); vox->OrigPos[2] = atof(buff2); }
    else if ((vox->mal==1)&&(Len>0))
    { sscanf(ln,"%f",&val);
      ++z;
      if (z==vox->N[2]) { z = 0; ++y; }
      if (y==vox->N[1]) { y = 0; ++x; }
      if (x<vox->N[0]) { vox->dat[x][y][z] = val;  }
      }
  if ((vox->mal==0)&&(vox->N[0]*vox->N[1]*vox->N[2]>0))
      Malloc_Voxel(vox,vox->N[0],vox->N[1],vox->N[2]);
  }
                                                                                                                  
  } /* while */
                                                                                                                  
 fclose(fp);
 printf("MIN %f %f %f GRIDLEN %f\n",
 vox->OrigPos[0],vox->OrigPos[1],vox->OrigPos[2],vox->grid_width);
                                                                                                                  
} /* end of Read_Voxel_File() */








void Write_Voxel_File(vox,fname)
 struct VOXEL *vox;
 char *fname;
{
 FILE *fp; 
 int x,y,z;

 fp = fopen(fname,"w"); 
 if (fp==NULL) {printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1);}
 printf("#Write_Voxel_File(%d %d %d) -> \"%s\"\n",vox->N[0],vox->N[1],vox->N[2],fname);
 /*
 fprintf(fp,"#COMMAND %s\n",PAR.COMMAND);
 fprintf(fp,"#DATE    %s\n",Get_Date_String());
 */
 fprintf(fp,"X %d\n",vox->N[0]);
 fprintf(fp,"Y %d\n",vox->N[1]);
 fprintf(fp,"Z %d\n",vox->N[2]);
 fprintf(fp,"GRIDLEN %f\n",vox->grid_width);
 fprintf(fp,"MINX %f\n",vox->OrigPos[0]);
 fprintf(fp,"MINY %f\n",vox->OrigPos[1]);
 fprintf(fp,"MINZ %f\n",vox->OrigPos[2]);
 for (x=0;x<vox->N[0];++x){
  for (y=0;y<vox->N[1];++y){
   for (z=0;z<vox->N[2];++z){ fprintf(fp,"%.10f\n",vox->dat[x][y][z]); }
    fprintf(fp,"\n"); 
   } 
  }
 fclose(fp);

} /* end of Write_Voxel_File() */


void Write_Voxel_File_In_PDB(vox,fname)
 struct VOXEL *vox;
 char *fname;
{
 FILE *fp; 
 int x,y,z,Natm,Nres;

 fp = fopen(fname,"w"); 
 if (fp==NULL) {printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1);}
 printf("#Write_Voxel_File_In_PDB(%d %d %d) -> \"%s\"\n",vox->N[0],vox->N[1],vox->N[2],fname);
 /*
 fprintf(fp,"REMARK #COMMAND %s\n",PAR.COMMAND);
 fprintf(fp,"REMARK #DATE    %s\n",Get_Date_String());
 */
 fprintf(fp,"REMARK X %d\n",vox->N[0]);
 fprintf(fp,"REMARK Y %d\n",vox->N[1]);
 fprintf(fp,"REMARK Z %d\n",vox->N[2]);
 fprintf(fp,"REMARK GRIDLEN %f\n",vox->grid_width);
 fprintf(fp,"REMARK MINX %f\n",vox->OrigPos[0]);
 fprintf(fp,"REMARK MINY %f\n",vox->OrigPos[1]);
 fprintf(fp,"REMARK MINZ %f\n",vox->OrigPos[2]);

 Natm = Nres = 0;
 for (x=0;x<vox->N[0];++x){
  for (y=0;y<vox->N[1];++y){
   ++Nres;
   if (Nres>=1000)Nres = 0;
   for (z=0;z<vox->N[2];++z){
     if (vox->dat[x][y][z]>0.0){ 
     ++Natm; 
     if (Natm>=100000)Natm = 0;
     fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f %e\n",
       Natm,"GRD","GRD",' ',Nres, 
       vox->OrigPos[0]+vox->grid_width*x,
       vox->OrigPos[1]+vox->grid_width*y,
       vox->OrigPos[2]+vox->grid_width*z,vox->dat[x][y][z],vox->dat[x][y][z],vox->dat[x][y][z]);
     } 
   } 
  }
 }
 fclose(fp);

} /* end of Write_Voxel_File_In_PDB() */







void Voxel_Statistics(vox,MIN,MAX,AVE,SD)
 struct VOXEL *vox;
 float *MIN,*MAX,*AVE,*SD;
{
 int x,y,z;
 double Sum,SSum,Nall;
 float min,max,v,ave,sd;

 printf("#Voxel_Statistics()\n");
 Sum = SSum = Nall = 0.0;  
 min = max = vox->dat[0][0][0];

 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){
       v = vox->dat[x][y][z];
       Sum += v;
       SSum += v*v;
       if (v>max) max = v;  
       if (v<min) min = v;  
       ++Nall;
     }
   }
 }

 ave = (float)(Sum/Nall);
 sd = (float)(SSum/Nall - ave*ave);
 if (sd>0.0) sd = sqrt(sd);
 printf("#ave %f min %f max %f sd %f\n",ave,min,max,sd);
 *MIN = min;
 *MAX = max;
 *AVE = ave;
 *SD  = sd;

} /* end of Voxel_Statistics() */





float Find_Threshold_Density_For_NumVoxel(vox,Nvoxel_thre)
 struct VOXEL *vox;
 int  Nvoxel_thre;   /* number of voxel for threshold */
{
 int x,y,z,n,Nvoxel;
 float *array,thre;
 
 Nvoxel = vox->N[0] * vox->N[1] * vox->N[2];
 printf("#Find_Threshold_Density_For_NumVoxel(vox,Nvoxel_thre:%d Nvoxel %d)\n",Nvoxel_thre,Nvoxel);
 if (Nvoxel_thre>=Nvoxel){
   printf("#ERROR(Find_Threshold_Density_For_NumVoxel(): Nvoxel_thre[%d] is over Nvoxel[%d].\n",Nvoxel_thre,Nvoxel);
   exit(1);
 }

 array = (float *)malloc(sizeof(float)*Nvoxel);

 n = 0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){
       array[n] = vox->dat[x][y][z];
       n += 1;
     }
   }
 }

 HeapSort_NoIndex(Nvoxel,array,'R');
 thre = (array[Nvoxel_thre-1] + array[Nvoxel_thre])/2.0; 
 printf("#density[%d]:%e density[%d]:%e thre %e\n",Nvoxel_thre-1,array[Nvoxel_thre-1],Nvoxel,array[Nvoxel_thre],thre); 

 free(array);
 return(thre);

} /* end of Find_Threshold_Density_For_NumVoxel() */






void Output_Voxel_Histgram(ofname,vox,MIN,MAX,WIDTH)
 char   *ofname;
 struct VOXEL *vox;
 float  MIN,MAX,WIDTH;
{
 FILE *fp;
 long Nall,Nmax;
 int *hist; 
 int i,x,y,z,ind;
 double Mbyte;
 Nmax = (int)floor((MAX-MIN)/WIDTH); 
 
 printf("#Output_Voxel_Histgram(MIN %f MAX %f WIDTH %f Nmax %d)\n",MIN,MAX,WIDTH,(int)Nmax);
 Mbyte = (double)sizeof(int)*(double)Nmax/1024.0/1024.0;
 printf("#%.2lf Mbyte is necessary\n",Mbyte); 
 if (Mbyte > 512) { printf("#ERROR:Memory size is too large\n"); exit(1); } 
 
 hist = (int *)malloc(sizeof(int)*Nmax);
 for (i=0;i<Nmax;++i) hist[i] = 0; 
 
 /** Make Histgram **/
 Nall = 0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){
       ind = (int)floor((vox->dat[x][y][z]-MIN)/WIDTH);
       hist[ind] += 1;
       ++Nall; 
     }
   }
 }

 /** Output **/
 fp = fopen(ofname,"w");
 if (fp==NULL) {printf("#ERROR:Can't write to histfile \"%s\"\n",ofname); exit(1);}
 for (i=0;i<Nmax;++i){
   fprintf(fp,"%f %d %f\n",i*WIDTH+MIN,hist[i],(float)hist[i]/(float)Nall);
 }
 fclose(fp);
 free(hist);

} /* end of Output_Voxel_Histgram() */


void Make_Reduced_Size_Voxel(vH,vO,size_scale)
 struct VOXEL *vH;  /* half sized voxel (to be caluclated) */
 struct VOXEL *vO;  /* original voxl */
 int    size_scale; /* 2,3,4,.... */
{
  int x,y,z,p,q,r,xx,yy,zz,k;

 /*

  Just adding the densities for 2*2*2,or 3*3*3, or 4*4*4,... voxels.
  Not take an average (not dividing by 2*2*2,or 3*3*3, or 4*4*4,...).
 */

  for (k=0;k<3;++k){
    vH->N[k] = vO->N[k]/size_scale;
    vH->OrigPos[k] = vO->OrigPos[k];
  }
  vH->grid_width =size_scale*vO->grid_width;
  
  printf("#Make_Half_Size_Voxel(size_scale:%d) [%d %d %d]->[%d %d %d]\n",
   size_scale,vO->N[0], vO->N[1], vO->N[2], vH->N[0], vH->N[1], vH->N[2]);

  for (x=0;x<vH->N[0];++x){
    for (y=0;y<vH->N[1];++y){
      for (z=0;z<vH->N[2];++z){
        vH->dat[x][y][z] = 0.0;
        for (p=0;p<size_scale;++p){
          xx = size_scale*x + p;
          for (q=0;q<size_scale;++q){
            yy = size_scale*y + q;
            for (r=0;r<size_scale;++r){
              zz = size_scale*z + r;
              if ((xx <vO->N[0])&&(yy <vO->N[1])&&(zz <vO->N[2])){
                vH->dat[x][y][z] += vO->dat[xx][yy][zz];
              }
            }
          } 
        }
/*
        printf("%d %d %d %f-->%f\n",x,y,z,vO->dat[x][y][z],vH->dat[x][y][z]);
 */
      }
    }
  }

} /* end of Make_Reduced_Size_Voxel() */

