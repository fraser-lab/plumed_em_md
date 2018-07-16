/*

 <GaussIO.c>
 for Input/Output of Gaussian Mixture Model


*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "globalvar.h"
#include "gauss.h"
#include "BasicIO.h"
#include "TabExpMinX.h"
#include "Matrix3D.h"
#include "Voxel.h"

/*** Functions (GLOBAL) ***/
double Gaussian3D_Value();
double Gaussian3D_Value_Table();
int Write_Gaussian3D_File();
int Write_Gaussian3D_File_with_PCaxis();
int Write_Gaussian3D_Center_Points();
char  Read_Gaussian3D_File();
int   Number_Of_Gaussian3D_in_File();
double Corr_Coeff_Bwn_Two_GAUSS3D_Arrays();
double Overlap_Integral_Bwn_Two_GAUSS3D_Arrays();
double Overlap_Integral_Bwn_Two_GAUSS3Ds();
double Overlap_Integral_Bwn_Two_Isotropic_GAUSS3Ds();
void Malloc_Voxel_From_Gaussians();
void Set_Voxel_Value_By_Gaussians();
void Set_Voxel_Value_By_Gaussian();
void Transform_Gaussians_By_Rmat_Gorig_Gnew();
int Check_NaN_Inf_Gaussian3D();
void show_Gaussian3D();
void Initialize_Gauss_From_M_and_CovM();
void Delete_ZeroWeight_gdfs_of_GMM();
void Delete_Identical_gdfs_of_GMM();
void Copy_GAUSS3D();
void Copy_GAUSS3D_Array();


/*** Functions (LOCAL) ***/
static void Random_Point_from_GaussianFunction();
static void Setup_LowerTriangle_Matrix();

/*****************/
/*** FUNCTIONS ***/
/*****************/


double Gaussian3D_Value(G,floatX)
 struct GAUSS3D *G;
 float  floatX[3];
{
 double D[3],xSx,val;

 D[0] = floatX[0] - G->M[0];
 D[1] = floatX[1] - G->M[1];
 D[2] = floatX[2] - G->M[2];

 xSx = 0.0;
 /*
 xSx += D[0]*(G->iCovM[0][0]*D[0] + G->iCovM[0][1]*D[1] + G->iCovM[0][2]*D[2]);
 xSx += D[1]*(G->iCovM[1][0]*D[0] + G->iCovM[1][1]*D[1] + G->iCovM[1][2]*D[2]);
 xSx += D[2]*(G->iCovM[2][0]*D[0] + G->iCovM[2][1]*D[1] + G->iCovM[2][2]*D[2]);
 */
 xSx += G->iCovM[0][0]*D[0]*D[0] + G->iCovM[1][1]*D[1]*D[1]    + G->iCovM[2][2]*D[2]*D[2];
 xSx += 2.0*(G->iCovM[0][1]*D[0]*D[1] + G->iCovM[0][2]*D[0]*D[2] + G->iCovM[1][2]*D[1]*D[2]);
 val = G->Cons * exp(-0.5*xSx);
 return(val);
}/* end of Gaussian3D_Value() */



double Gaussian3D_Value_Table(G,floatX)
 struct GAUSS3D *G;
 float  floatX[3];
{
 double D[3],xSx,val;

 D[0] = floatX[0] - G->M[0];
 D[1] = floatX[1] - G->M[1];
 D[2] = floatX[2] - G->M[2];

 xSx = 0.0;
 /*
 xSx += D[0]*(G->iCovM[0][0]*D[0] + G->iCovM[0][1]*D[1] + G->iCovM[0][2]*D[2]);
 xSx += D[1]*(G->iCovM[1][0]*D[0] + G->iCovM[1][1]*D[1] + G->iCovM[1][2]*D[2]);
 xSx += D[2]*(G->iCovM[2][0]*D[0] + G->iCovM[2][1]*D[1] + G->iCovM[2][2]*D[2]);
 */
 xSx += G->iCovM[0][0]*D[0]*D[0] + G->iCovM[1][1]*D[1]*D[1]    + G->iCovM[2][2]*D[2]*D[2];
 xSx += 2.0*(G->iCovM[0][1]*D[0]*D[1] + G->iCovM[0][2]*D[0]*D[2] + G->iCovM[1][2]*D[1]*D[2]);
 if (isnan(D[0])!=0){printf("D[0] is nan!\n"); exit(1);}
 if (isnan(D[1])!=0){printf("D[1] is nan!\n"); exit(1);}
 if (isnan(D[2])!=0){printf("D[2] is nan!\n"); exit(1);}
 if (isnan(G->iCovM[0][0])!=0){printf("iCovM[0][0] is nan!\n"); exit(1);}
 if (isnan(G->iCovM[0][1])!=0){printf("iCovM[0][1] is nan!\n"); exit(1);}
 if (isnan(G->iCovM[0][2])!=0){printf("iCovM[0][2] is nan!\n"); exit(1);}
 if (isnan(G->iCovM[1][0])!=0){printf("iCovM[1][0] is nan!\n"); exit(1);}
 if (isnan(G->iCovM[1][1])!=0){printf("iCovM[1][1] is nan!\n"); exit(1);}
 if (isnan(G->iCovM[1][2])!=0){printf("iCovM[1][2] is nan!\n"); exit(1);}
 if (isnan(G->iCovM[2][0])!=0){printf("iCovM[2][0] is nan!\n"); exit(1);}
 if (isnan(G->iCovM[2][1])!=0){printf("iCovM[2][1] is nan!\n"); exit(1);}
 if (isnan(G->iCovM[2][2])!=0){printf("iCovM[2][2] is nan!\n"); exit(1);}

 if (isnan(xSx)!=0){
   printf("xSx is nan!\n"); 
   printf("floatX %f %f %f\n",floatX[0],floatX[1],floatX[2]); 
   printf("G->M   %lf %lf %lf\n",G->M[0],G->M[1],G->M[2]);
   printf("iCovM %lf %lf %lf\n", G->iCovM[0][0], G->iCovM[0][1], G->iCovM[0][2]);
   printf("iCovM %lf %lf %lf\n", G->iCovM[1][0], G->iCovM[1][1], G->iCovM[1][2]);
   printf("iCovM %lf %lf %lf\n", G->iCovM[2][0], G->iCovM[2][1], G->iCovM[2][2]);
   exit(1);}

 if (PAR.TabExpMinX == 'T'){
  val = G->Cons * Exp_minus_X_Table(0.5*xSx,&TABEXP);
 }
 else
 { val = G->Cons * exp(-0.5*xSx);}
 
 return(val);
}/* end of Gaussian3D_Value_Table() */







int Write_Gaussian3D_File_with_PCaxis(fname,file_mode,Ngauss,Garray,chain,comment)
 char *fname;
 char file_mode; /* 'a'ppend, 'w'rite */
 int Ngauss;
 struct GAUSS3D *Garray;
 char  chain; 
 char *comment;
{
 FILE *fp;
 int g,i; 

 if (file_mode=='a'){ fp = fopen(fname,"a"); }
 else{ fp = fopen(fname,"w"); }
 if (fp==NULL) {printf("#ERROR:Can't write to gaussfile \"%s\"\n",fname); return(0);} 
 printf("#Write_Gaussian3D_File_with_PCaxis()-->\"%s\"\n",fname);
 fprintf(fp,"HEADER 3D Gaussian Mixture Model\n");
 fprintf(fp,"REMARK COMMAND %s\n",PAR.COMMAND);
 Set_END_DATE();
 fprintf(fp,"REMARK START_DATE %s\n",PAR.START_DATE);
 fprintf(fp,"REMARK END_DATE   %s\n",PAR.END_DATE);
 fprintf(fp,"REMARK COMP_TIME_SEC  %lf %e\n",PAR.COMP_TIME_SEC,PAR.COMP_TIME_SEC);
 fprintf(fp,"REMARK FILENAME %s\n",fname);
 fprintf(fp,"REMARK NGAUSS %d\n",Ngauss);
 if (comment[0]!='\0') fprintf(fp,"REMARK COMMENT %s\n",comment);

 for (g=0;g<Ngauss;++g){
   fprintf(fp,"REMARK GAUSS%4d W %.10lf\n",g+1,Garray[g].Weight);
   fprintf(fp,"REMARK GAUSS%4d M %lf %lf %lf\n",g+1,Garray[g].M[0], Garray[g].M[1], Garray[g].M[2]);
   fprintf(fp,"REMARK GAUSS%4d CovM  xx %15.10lf xy %15.10lf xz %15.10lf\n",
        g+1,Garray[g].CovM[0][0],Garray[g].CovM[0][1],Garray[g].CovM[0][2]);
   fprintf(fp,"REMARK GAUSS%4d CovM  yy %15.10lf yz %15.10lf zz %15.10lf\n",
        g+1,Garray[g].CovM[1][1],Garray[g].CovM[1][2],Garray[g].CovM[2][2]);
   fprintf(fp,"REMARK GAUSS%4d PCvar   %lf %lf %lf\n",g+1,Garray[g].PCvar[0], Garray[g].PCvar[1], Garray[g].PCvar[2]);
   fprintf(fp,"REMARK GAUSS%4d PCaxis1 %lf %lf %lf\n",g+1,Garray[g].PCaxis[0][0], Garray[g].PCaxis[0][1], Garray[g].PCaxis[0][2]);
   fprintf(fp,"REMARK GAUSS%4d PCaxis2 %lf %lf %lf\n",g+1,Garray[g].PCaxis[1][0], Garray[g].PCaxis[1][1], Garray[g].PCaxis[1][2]);
   fprintf(fp,"REMARK GAUSS%4d PCaxis3 %lf %lf %lf\n",g+1,Garray[g].PCaxis[2][0], Garray[g].PCaxis[2][1], Garray[g].PCaxis[2][2]);
   fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n",
       7*g+1,"CEN","GAU",chain,g+1, Garray[g].M[0],Garray[g].M[1],Garray[g].M[2],Garray[g].Weight,Garray[g].Weight);
   for (i=0;i<3;++i){
     fprintf(fp,"HETATM%5d %3s%d %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n",
         7*g+2*i+2,"PC",i+1,"GAU",chain,g+1, 
         Garray[g].M[0] + sqrt(Garray[g].PCvar[i])*Garray[g].PCaxis[i][0],
         Garray[g].M[1] + sqrt(Garray[g].PCvar[i])*Garray[g].PCaxis[i][1],
         Garray[g].M[2] + sqrt(Garray[g].PCvar[i])*Garray[g].PCaxis[i][2], Garray[g].Weight,Garray[g].Weight);
     fprintf(fp,"HETATM%5d %3s%d %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n",
         7*g+2*i+3,"PC",i+1,"GAU",chain,g+1, 
         Garray[g].M[0] - sqrt(Garray[g].PCvar[i])*Garray[g].PCaxis[i][0],
         Garray[g].M[1] - sqrt(Garray[g].PCvar[i])*Garray[g].PCaxis[i][1],
         Garray[g].M[2] - sqrt(Garray[g].PCvar[i])*Garray[g].PCaxis[i][2], Garray[g].Weight,Garray[g].Weight);
   } 
 } 
 fprintf(fp,"TER\n");
 for (g=0;g<Ngauss;++g){
   fprintf(fp,"CONECT %4d %4d %4d %4d %4d\n", 7*g+1, 7*g+2, 7*g+3, 7*g+4, 7*g+5);
   fprintf(fp,"CONECT %4d %4d %4d\n", 7*g+1, 7*g+6, 7*g+7);
 }
 fclose(fp);
 return(1);
} /* end of Write_Gaussian3D_File_with_PCaxis() */









int Write_Gaussian3D_File(fname,file_mode,Ngauss,Garray,chain,comment)
 char *fname;
 char file_mode; /* 'a'ppend, 'w'rite */
 int Ngauss;
 struct GAUSS3D *Garray;
 char  chain; 
 char *comment;
{
 FILE *fp;
 int g; 

 if (file_mode=='a'){ fp = fopen(fname,"a"); }
 else{ fp = fopen(fname,"w"); }
 if (fp==NULL) {printf("#ERROR:Can't write to gaussfile \"%s\"\n",fname); return(0);} 
 printf("#Write_Gaussian3D_File()-->\"%s\"\n",fname);
 fprintf(fp,"HEADER 3D Gaussian Mixture Model\n");
 fprintf(fp,"REMARK COMMAND %s\n",PAR.COMMAND);
 Set_END_DATE();
 fprintf(fp,"REMARK START_DATE %s\n",PAR.START_DATE);
 fprintf(fp,"REMARK END_DATE   %s\n",PAR.END_DATE);
 fprintf(fp,"REMARK COMP_TIME_SEC  %lf %e\n",PAR.COMP_TIME_SEC,PAR.COMP_TIME_SEC);
 fprintf(fp,"REMARK FILENAME %s\n",fname);
 fprintf(fp,"REMARK NGAUSS %d\n",Ngauss);
 if (comment[0]!='\0') fprintf(fp,"REMARK COMMENT %s\n",comment);

 for (g=0;g<Ngauss;++g){
   fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n",
       g+1,"GAU","GAU",chain,g+1, Garray[g].M[0],Garray[g].M[1],Garray[g].M[2],Garray[g].Weight,Garray[g].Weight);
   fprintf(fp,"REMARK GAUSS%4d W %.10lf\n",g+1,Garray[g].Weight);
   fprintf(fp,"REMARK GAUSS%4d det  %.10lf\n",g+1,Garray[g].det);
   fprintf(fp,"REMARK GAUSS%4d Cons %.10lf\n",g+1,Garray[g].Cons);
   fprintf(fp,"REMARK GAUSS%4d M %lf %lf %lf\n",g+1,Garray[g].M[0], Garray[g].M[1], Garray[g].M[2]);
   fprintf(fp,"REMARK GAUSS%4d CovM  xx %15.10lf xy %15.10lf xz %15.10lf\n",
        g+1,Garray[g].CovM[0][0],Garray[g].CovM[0][1],Garray[g].CovM[0][2]);
   fprintf(fp,"REMARK GAUSS%4d CovM  yy %15.10lf yz %15.10lf zz %15.10lf\n",
        g+1,Garray[g].CovM[1][1],Garray[g].CovM[1][2],Garray[g].CovM[2][2]);
 } 
 fprintf(fp,"TER\n");
 fclose(fp);
 return(1);
} /* end of Write_Gaussian3D_File() */










int Write_Gaussian3D_Center_Points(fname,fmode,Ngauss,Garray,chain,offset_atomnum,comment)
 char *fname;
 char fmode;  /* 'a'ppend or 'w'rite */
 int Ngauss;
 struct GAUSS3D *Garray;
 char  chain; 
 int   offset_atomnum;
 char *comment;
{
 FILE *fp;
 int g;

 if (fmode=='a'){ fp = fopen(fname,"a");}
          else  { fp = fopen(fname,"w");}

 if (fp==NULL) {printf("#ERROR:Can't write to gaussfile \"%s\"\n",fname); return(0);}  

 printf("#Write_Gaussian3D_Center_Points()-->\"%s\"\n",fname);
 fprintf(fp,"HEADER 3D Gaussian Mixture Model\n");
 fprintf(fp,"REMARK NGAUSS %d\n",Ngauss);
 if (comment[0]!='\0') fprintf(fp,"REMARK COMMENT %s\n",comment);
 for (g=0;g<Ngauss;++g){
   fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n",
       g+1+offset_atomnum,"CEN","CEN",chain,g+1, Garray[g].M[0],Garray[g].M[1],Garray[g].M[2],Garray[g].Weight,Garray[g].Weight);
 } 
 fprintf(fp,"TER\n");
 fclose(fp);
 return(1);
} /* end of Write_Gaussian3D_Center_Points() */






char Read_Gaussian3D_File(fname,Ngauss_malloc,Ngauss,Garray)
 char *fname;
 int Ngauss_malloc;
 int *Ngauss;
 struct GAUSS3D *Garray;
{
 FILE *fp;
 int  g,gnum,gnum0,w,L,i;
 char line[512],buff[512],word[10][100],chain;
 int  Wsta[100],Wend[100],Nword; 
 double det,Weight_all;
 
 for (g=0;g<(*Ngauss);++g){
   Garray[g].det = 0.0;
   for (i=0;i<3;++i){
     Garray[g].M[i] = 0.0;
     Garray[g].CovM[i][i] = 0.0;
     Garray[g].iCovM[i][i] = 0.0;
   }
 }

 *Ngauss = 0; g = -1; gnum0 = -1;  chain = ' ';
 Weight_all = 0.0;
 fp = fopen(fname,"r");
 if (fp==NULL) {printf("#ERROR:Can't open gaussfile \"%s\"\n",fname); exit(1);}
 printf("#Read_Gaussian3D_File(\"%s\")\n",fname);
/*
>> EXAMPLE OF THE FILE FOR GAUSSIAN MIXTURE MODEL << 
HEADER 3D Gaussian Mixture Model
REMARK COMMAND gmconvert -ipdb pdb1mbd.ent -ng 2
REMARK START_DATE Feb 8,2013 12:0:31
REMARK END_DATE   Feb 8,2013 12:0:31
REMARK COMP_TIME_SEC  0.098650 9.865010e-02
REMARK FILENAME pdb1mbd_2.gmm
REMARK NGAUSS 2
HETATM    1  GAU GAU A   1      14.387  16.299  11.800 0.471 0.471
REMARK GAUSS   1 W 0.4707285444
REMARK GAUSS   1 M 14.386588 16.298907 11.799829
REMARK GAUSS   1 CovM  xx   56.9545342387 xy  -14.7197891094 xz   -8.8706625730
REMARK GAUSS   1 CovM  yy   50.4607933346 yz    1.0600698705 zz   32.0378248033
HETATM    2  GAU GAU A   2      16.337  24.890  -3.056 0.529 0.529
REMARK GAUSS   2 W 0.5292714556
REMARK GAUSS   2 M 16.336505 24.889884 -3.055859
REMARK GAUSS   2 CovM  xx  106.4774042405 xy  -19.7884844863 xz  -16.8298857592
REMARK GAUSS   2 CovM  yy   37.3625392323 yz    7.5580564129 zz   34.9972721500
TER
*/
 while (feof(fp)==0){
   line[0] = '\0';
   fgets(line,511,fp);

   L = strlen(line); 
   if ((L>0)&&(line[L-1]=='\n')) line[L-1] = '\0';

   if (strncmp(line,"HETATM",6) == 0) chain = line[21];
 
   if ((strncmp(line,"REMARK GAUSS_MOLECULE",21)!=0)&& (strncmp(line,"REMARK GAUSS-MOLECULE",21)!=0)&&
       (strncmp(line,"REMARK GAUSS",12)==0)){
     Get_Part_Of_Line(buff,line,12,512);
     Split_to_Words(buff,' ',&Nword,Wsta,Wend,100);
     for (w=0;w<Nword;++w) Get_Part_Of_Line(word[w],buff,Wsta[w],Wend[w]);
     gnum = atoi(word[0]); 
     if (gnum0!=gnum) { g += 1; *Ngauss += 1;}
     if (*Ngauss > Ngauss_malloc){
       printf("#ERROR:Ngauss %d is over Ngauss_malloc %d.\n",*Ngauss, Ngauss_malloc);
       exit(1);
     } 
     gnum0 = gnum;
     /* printf("g %d gnum %d '%s'\n",g,gnum,line); */
         if (strncmp(word[1],"W",1)==0){ 
      Garray[g].Weight = atof(word[2]);
      Weight_all += Garray[g].Weight;
    } 
    else if (strncmp(word[1],"M",1)==0){
      Garray[g].M[0] = atof(word[2]);
      Garray[g].M[1] = atof(word[3]);
      Garray[g].M[2] = atof(word[4]); 
    }
    else if ((strcmp(word[1],"CovM")==0)&&((strcmp(word[2],"00")==0)||(strcmp(word[2],"xx")==0))){
      Garray[g].CovM[0][0] = atof(word[3]);
      Garray[g].CovM[0][1] = Garray[g].CovM[1][0] = atof(word[5]);
      Garray[g].CovM[0][2] = Garray[g].CovM[2][0] = atof(word[7]);
    }
    else if ((strcmp(word[1],"CovM")==0)&&((strcmp(word[2],"11")==0)||(strcmp(word[2],"yy")==0))){
      Garray[g].CovM[1][1] = atof(word[3]);
      Garray[g].CovM[1][2] = Garray[g].CovM[2][1] = atof(word[5]);
      Garray[g].CovM[2][2] = atof(word[7]); 
    }

  } /* REMARK GAUSS */

 } /* while */

 fclose(fp);
 
 /*** Normalize Weights ***/

 /*
 for (g=0;g<(*Ngauss);++g){
   Garray[g].Weight  = Garray[g].Weight/Weight_all;
 }
  */

 /*** Calculate iCovM and det ***/
 for (g=0;g<(*Ngauss);++g){
   Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].iCovM, Garray[g].CovM, &det);
   Garray[g].det = Determinant_Matrix3D(Garray[g].CovM);
   Garray[g].Cons = 1.0/(pow(2*M_PI,1.5)*sqrt(Garray[g].det));
  /* printf("#Gauss[%d] det %lf\n",g,Garray[g].det);
    printf("#[%d/%d] det %lf\n",g, *Ngauss, Garray[g].det);
   */
 }

 if (*Ngauss<=0){
  printf("#ERROR:gaussfile \"%s\" does not contain any gaussians.\n",fname);
  exit(1);
 }

 return(chain);

} /* end of Read_Gaussian3D_File() */




int Number_Of_Gaussian3D_in_File(fname)
 char *fname;
{
 FILE *fp;
 char line[512];
 int  L,Ngauss;

 printf("#Number_Of_Gaussian3D_in_File('%s')\n",fname); fflush(stdout);
 Ngauss = 0;
 fp = fopen(fname,"r");

 if (fp==NULL) {printf("#ERROR:Can't open gaussfile \"%s\"\n",fname); exit(1);}
 while (feof(fp)==0){
  line[0] = '\0';
  fgets(line,511,fp);
  L = strlen(line);
  if ((L>0)&&(line[L-1]=='\n')) line[L-1] = '\0';

/* Count "REMARK GAUSS xxx M" line */
/* REMARK GAUSS   2 M 16.336505 24.889884 -3.055859 */
  if ((strncmp(line,"REMARK GAUSS ",13) == 0)&&(line[17]=='M')){
     Ngauss += 1; 
  }
 }
 fclose(fp);
 printf("#Ngauss %d\n",Ngauss);
 return(Ngauss);
} /* end of Number_Of_Gaussian3D_in_File() */








double Corr_Coeff_Bwn_Two_GAUSS3D_Arrays(NgaussA,gAarray,NgaussB,gBarray)
 int NgaussA;
 struct GAUSS3D *gAarray;
 int NgaussB;
 struct GAUSS3D *gBarray;
{
 double ovAB,ovAA,ovBB,CC;

 ovAB =  Overlap_Integral_Bwn_Two_GAUSS3D_Arrays(NgaussA,gAarray,NgaussB,gBarray);
 ovAA =  Overlap_Integral_Bwn_Two_GAUSS3D_Arrays(NgaussA,gAarray,NgaussA,gAarray);
 ovBB =  Overlap_Integral_Bwn_Two_GAUSS3D_Arrays(NgaussB,gBarray,NgaussB,gBarray);
 
 if ((ovAA>0.0)&&(ovBB>0.0)){
   CC = ovAB/sqrt(ovAA*ovBB);
 }
 else CC = 0.0;
 /* printf("ovAB %e ovAA %e ovBB %e CC %lf %e\n",ovAB,ovAA,ovBB,CC); */
 return(CC);
} /* end of Corr_Coeff_Bwn_Two_GAUSS3D_Arrays() */



                                                                                                                    
double Overlap_Integral_Bwn_Two_GAUSS3D_Arrays(NgaussA,gAarray,NgaussB,gBarray)
 int NgaussA;
 struct GAUSS3D *gAarray;
 int NgaussB;
 struct GAUSS3D *gBarray;
{
 int i,j;
 double ov, ov_all;
 ov_all = 0.0;
 for (i=0;i<NgaussA;++i){
  for (j=0;j<NgaussB;++j){
    ov = Overlap_Integral_Bwn_Two_GAUSS3Ds(&(gAarray[i]),&(gBarray[j]));
    ov_all += gAarray[i].Weight * gBarray[j].Weight * ov;
  }
 }
 return(ov_all);
} /* end of Overlap_Integral_Bwn_GAUSS3D_Arrays() */



double Overlap_Integral_Bwn_Two_GAUSS3Ds(gA,gB)
 struct GAUSS3D *gA,*gB;
{
 int i;
 double K[3][3],invK[3][3],DetK,ov;
 double nm[3];
 double qform;
 /*
   Overlap = 1/[(2pi)**3/2 * sqrt(|P+Q|)] *
             exp[-1/2 * tr_(n-m) inv(P+Q) (n-m)]
            = 1/(2pi)**3/2 /sqrt(|K|)*exp[-1/2 * (n-m) invK (n-m)]
    where P = gA->CovM, Q = gB->CovM, m = gA->M, n = gB->M.
     K   = P + Q
  invK   = inv(P + Q)
 */
 Add_Matrix3D(K,gA->CovM,gB->CovM);
 /* Cal_Inverse_Matrix3D_by_Cramer_Rule(invK,K,&DetK); */  
 Cal_Inverse_Matrix3D_by_Cramer_Rule_Symmetric(invK,K,&DetK); 
 for (i=0;i<3;++i){ nm[i] =  gB->M[i] - gA->M[i];}
 qform =  Quadratic_Form_3D(invK,nm);
 /* printf("nm %lf %lf %lf qform %lf\n",nm[0],nm[1],nm[2],qform);  */
 ov = exp(-0.5*qform) * one_over_2pi_32 / sqrt(DetK);
 return(ov);
} /* end of double Overlap_Integral_Bwn_Two_GAUSS3Ds() */


double Overlap_Integral_Bwn_Two_Isotropic_GAUSS3Ds(gA,gB)
 struct GAUSS3D *gA,*gB;
{
 int i;
 double qform,DetK,ov,Kdiag,invKdiag,nm[3];
 /*
   Overlap = 1/[(2pi)**3/2 * sqrt(|P+Q|)] *
             exp[-1/2 * tr_(n-m) inv(P+Q) (n-m)]
            = 1/(2pi)**3/2 /sqrt(|K|)*exp[-1/2 * (n-m) invK (n-m)]
    where P = gA->CovM, Q = gB->CovM, m = gA->M, n = gB->M.
     K   = P + Q
  invK   = inv(P + Q)
 */

 for (i=0;i<3;++i){ nm[i] =  gB->M[i] - gA->M[i];}
 Kdiag = gA->CovM[0][0] + gB->CovM[0][0];
 invKdiag = 1.0/Kdiag;
 DetK = Kdiag * Kdiag * Kdiag;
 qform = invKdiag * (nm[0]*nm[0] + nm[1]*nm[1] + nm[2]*nm[2]); 
 /*
 Add_Matrix3D(K,gA->CovM,gB->CovM);
 Cal_Inverse_Matrix3D_by_Cramer_Rule(invK,K,&DetK);
 qform =  Quadratic_Form_3D(invK,nm);
 */
 ov = exp(-0.5*qform) * one_over_2pi_32 / sqrt(DetK);
 
 return(ov);
} /* end of double Overlap_Integral_Bwn_Two_Isotropic_GAUSS3Ds() */





void Malloc_Voxel_From_Gaussians(vox,Ngauss,Ga,MaxSD)
 struct VOXEL *vox;
 int    Ngauss;
 struct GAUSS3D *Ga;
 double MaxSD;  /* SD for maximum density boundary. (about from 4.0 to 8.0) */
{
 double MeanG[3],VarG[3],SD,CovM[3][3],TotalWeight;
 int i,j,g;

 printf("#Malloc_Voxel_From_Gaussians(Ngauss %d MaxSD %lf gw %lf)\n",Ngauss,MaxSD,vox->grid_width);
 /*** (1) Calculate Mean and Var for Gaussian Array ***/
 for (i=0;i<3;++i) MeanG[i] = VarG[i] = 0.0;

 TotalWeight = 0.0;
 for (g=0;g<Ngauss;++g){ 
   TotalWeight += Ga[g].Weight;
 }

 for (g=0;g<Ngauss;++g){
   for (i=0;i<3;++i){ MeanG[i] += Ga[g].Weight * Ga[g].M[i]; }
 }
 for (i=0;i<3;++i){ MeanG[i] = MeanG[i]/TotalWeight;}



 for (i=0;i<3;++i){
   for (j=0;j<3;++j){CovM[i][j] = 0.0;}
 }

 for (g=0;g<Ngauss;++g){
  for (i=0;i<3;++i){
   for (j=0;j<3;++j){
    CovM[i][j] += Ga[g].Weight*(Ga[g].M[i]-MeanG[i])*(Ga[g].M[j]-MeanG[j]);
    CovM[i][j] += Ga[g].Weight*Ga[g].CovM[i][j];
   }
  }
 }

 for (i=0;i<3;++i){
   for (j=0;j<3;++j){
     CovM[i][j] = CovM[i][j]/TotalWeight;
   }
 }


 VarG[0] = CovM[0][0];
 VarG[1] = CovM[1][1];
 VarG[2] = CovM[2][2];
/*
  printf("CovM %lf %lf %lf\n",CovM[0][0],CovM[0][1],CovM[0][2]);
  printf("CovM %lf %lf %lf\n",CovM[1][0],CovM[1][1],CovM[1][2]);
  printf("CovM %lf %lf %lf\n",CovM[2][0],CovM[2][1],CovM[2][2]);
  */

 /*** (2) Setup Min and Max by MaxSD * sqrt(VarG) ***/

 for (i=0;i<3;++i){
   SD = sqrt(VarG[i]);
   vox->OrigPos[i] = MeanG[i] - MaxSD * SD;
   vox->N[i]   = (int)ceil((2.0*MaxSD * SD)/vox->grid_width);
   printf("#g %d SD %lf MIN %f MAX %f N %d\n",i,SD,vox->OrigPos[i],MeanG[i]+MaxSD*SD,vox->N[i]);
 }

 Malloc_Voxel(vox,vox->N[0],vox->N[1],vox->N[2]);


} /* end of Malloc_Voxel_From_Gaussians() */




void Set_Voxel_Value_By_Gaussians(vox,Ngauss,Garray)
 struct VOXEL *vox;
 int Ngauss;
 struct GAUSS3D *Garray;
{
 int g,x,y,z;
 float X[3];


 for (x=0;x< vox->N[0];++x){
   X[0] = vox->OrigPos[0] + vox->grid_width * x;
   for (y=0;y< vox->N[1];++y){
     X[1] = vox->OrigPos[1] + vox->grid_width * y;
     for (z=0;z< vox->N[2];++z){
       X[2] = vox->OrigPos[2] + vox->grid_width * z;
       vox->dat[x][y][z] = 0.0;
       for (g=0;g<Ngauss;++g)
         vox->dat[x][y][z] += (float)(Garray[g].Weight * Gaussian3D_Value(&(Garray[g]),X));
      }
    }
 }


} /* end of Set_Voxel_Value_By_Gaussians() */

void Set_Voxel_Value_By_Gaussian(vox,G)
 struct VOXEL *vox;
 struct GAUSS3D *G;
{
 int x,y,z;
 float X[3];


 for (x=0;x< vox->N[0];++x){
   X[0] = vox->OrigPos[0] + vox->grid_width * x;
   for (y=0;y< vox->N[1];++y){
     X[1] = vox->OrigPos[1] + vox->grid_width * y;
     for (z=0;z< vox->N[2];++z){
       X[2] = vox->OrigPos[2] + vox->grid_width * z;
       // set value
       vox->dat[x][y][z] = (float)(G->Weight * Gaussian3D_Value(G,X));
      }
    }
 }

}


void Transform_Gaussians_By_Rmat_Gorig_Gnew(Ngauss, garray, Rmat,Gorig, Gnew)
  int Ngauss;
  struct GAUSS3D *garray;
  double Rmat[3][3];
  double Gorig[3],Gnew[3]; 
{
  int g,i;
  double mgo[3],rmgo[3];
  double RCovM[3][3],det;
 
  for (g=0;g<Ngauss;++g){
    Sub_Vec3D(mgo, garray[g].M, Gorig);
    Multiply_Matrix3D_Vec3D(rmgo,Rmat,mgo);
    for (i=0;i<3;++i){ garray[g].M[i] = rmgo[i] + Gnew[i];}
    Transform_RAtR_Matrix3D(RCovM, garray[g].CovM,Rmat);
    Equal_Matrix3D(garray[g].CovM,RCovM);
    Cal_Inverse_Matrix3D_by_Cramer_Rule(garray[g].iCovM, garray[g].CovM,&det);
  }


} /* end of Transform_Gaussians_By_Rmat_Gorig_Gnew() */


int Check_NaN_Inf_Gaussian3D(G)
 struct GAUSS3D *G;
{
  int i,j;
 
  if ((isnan(G->Weight)!=0)||(isnan(G->Weight)!=0)){return(1);}
  for (i=0;i<3;++i){
    if ((isnan(G->M[i])!=0)||(isnan(G->M[i])!=0)){return(1);}
    if (isinf(G->M[i])!=0){return(1);}

    for (j=0;j<3;++j){
      if ((isnan(G->CovM[i][j])!=0)||(isnan(G->CovM[i][j])!=0)){return(1);}
      if (isinf(G->CovM[i][j])!=0){return(1);}
      if ((isnan(G->iCovM[i][j])!=0)||(isnan(G->iCovM[i][j])!=0)){return(1);}
      if (isinf(G->iCovM[i][j])!=0){return(1);}
    }
  } 

  return(0);
}


void show_Gaussian3D(G)
 struct GAUSS3D *G;
{
   printf("REMARK GAUSS W %.10lf\n",G->Weight);
   printf("REMARK GAUSS M %lf %lf %lf\n",G->M[0], G->M[1], G->M[2]);
   printf("REMARK GAUSS CovM  xx %15.10lf xy %15.10lf xz %15.10lf\n",
        G->CovM[0][0],G->CovM[0][1],G->CovM[0][2]);
   printf("REMARK GAUSS CovM  yy %15.10lf yz %15.10lf zz %15.10lf\n",
        G->CovM[1][1],G->CovM[1][2],G->CovM[2][2]);
}



void Initialize_Gauss_From_M_and_CovM(Ngauss,Garray,M,CovM)
 int Ngauss;
 struct GAUSS3D *Garray;
 double M[3];
 double CovM[3][3];
{
 int g,i,j;
 double loCovM[3][3]; 

 Setup_LowerTriangle_Matrix(loCovM,CovM);
 /* printf("#Initialize_Gauss_From_Randomly_Chosen_Voxel()\n"); */
 for (g=0;g<Ngauss;++g){
   Random_Point_from_GaussianFunction(Garray[g].M, M, CovM, loCovM);
   for (i=0;i<3;++i){
      for (j=0;j<3;++j){
         Garray[g].CovM[i][j] = CovM[i][j];
       }
    }

   Garray[g].Weight = 1.0/(double)Ngauss;
   Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));
   Garray[g].Cons = 1.0/(two_pi_32*sqrt(Garray[g].det));
/*
   printf("#g %d     M %lf %lf %lf\n",g,Garray[g].M[0],Garray[g].M[1],Garray[g].M[2]); 
   printf("#g %d CovM %lf %lf %lf\n",g,Garray[g].CovM[0][0],Garray[g].CovM[0][1],Garray[g].CovM[0][2]); 
   printf("#g %d CovM %lf %lf %lf\n",g,Garray[g].CovM[1][0],Garray[g].CovM[1][1],Garray[g].CovM[1][2]); 
   printf("#g %d CovM %lf %lf %lf\n",g,Garray[g].CovM[2][0],Garray[g].CovM[2][1],Garray[g].CovM[2][2]); 
   printf("#g %d iCovM %lf %lf %lf\n",g,Garray[g].iCovM[0][0],Garray[g].iCovM[0][1],Garray[g].iCovM[0][2]); 
   printf("#g %d iCovM %lf %lf %lf\n",g,Garray[g].iCovM[1][0],Garray[g].iCovM[1][1],Garray[g].iCovM[1][2]); 
   printf("#g %d iCovM %lf %lf %lf\n",g,Garray[g].iCovM[2][0],Garray[g].iCovM[2][1],Garray[g].iCovM[2][2]); 
 */ 
 }

} /* end of Initialize_Gauss_From_M_and_CovM() */





void Random_Point_from_GaussianFunction(Rpnt,M,CovM,loCovM)
 double Rpnt[3];    /* Random point  (to be calculated) */
 double M[3];       /* center of Gaussian function */
 double CovM[3][3]; /* covariance matrix of Gaussian function */
 double loCovM[3][3]; /* lower triangle matrix for CovM */
{
 int i;
 double U1,U2,X[3];

 for (i=0;i<3;++i){
   U1 = (double)(rand()+1)/(double)RAND_MAX;
   if (U1>1.0) U1 = 1.0;
   if (U1<0.0) {printf("#WARNING U1 is less than 0?\n"); fflush(stdout); U1 *= -1.0;}
   U2 = (double)rand()/(double)RAND_MAX;
   if (U2==0.0) {printf("#WARNING U2 is zero.\n"); fflush(stdout);}
   X[i] = sqrt(-2*log(U1))*cos(2*M_PI*U2);
 }


 Rpnt[0] = M[0] + loCovM[0][0]*X[0];
 Rpnt[1] = M[1] + loCovM[1][0]*X[0] + loCovM[1][1]*X[1];
 Rpnt[2] = M[2] + loCovM[2][0]*X[0] + loCovM[2][1]*X[1] + loCovM[2][2]*X[2];


} /* end of Random_Point_from_GaussianMatrix() */





void Setup_LowerTriangle_Matrix(loCovM,CovM)
  double  loCovM[3][3]; /* lower triangle matrix (to be calculated) */
  double  CovM[3][3];   /* covariance matrix */
{
 /**
    Cholesky decomposition of CovM(S)
 
    |L00 0   0  |   |L00 L10 L20|   |S00 S01 S02|
    |L10 L11 0  | x |0   L11 L21| = |S01 S11 S12|
    |L20 L21 L22|   |0   0   L22|   |S02 S12 S22|
 
        L00^2   = S00              --> L00 = sqrt(S00)
        L00*L10 = S01              --> L10 = S01/L00
        L10^2 + L11^2 = S11        --> L11 = sqrt(S11-L10*L10)
        L00*L20 = S02              --> L20 = S02/L00
        L10*L20 + L11*L21 = S12    --> L21 = (S12-L10*L20)/L11
        L20^2 + L21^2 + L22^2 =S22 --> K22 = sqrt(S22-L20^2-L21^2);

 **/
 double L00,L10,L11,L20,L21,L22;

 L00 = sqrt(CovM[0][0]);
 L10 = CovM[0][1]/L00;
 L11 = sqrt(CovM[1][1] - L10*L10);
 L20 = CovM[0][2]/L00;
 L21 = (CovM[1][2]-L10*L20)/L11;
 L22 = sqrt(CovM[2][2]-L20*L20-L21*L21);

 loCovM[0][0] = L00;
 loCovM[1][0] = L10;
 loCovM[1][1] = L11;
 loCovM[2][0] = L20;
 loCovM[2][1] = L21;
 loCovM[2][2] = L22;

} /* end of Setup_LowerTriangle_Matrix() */



void Delete_ZeroWeight_gdfs_of_GMM(Ngauss,Garray)
 int *Ngauss;
 struct GAUSS3D *Garray;
{
  struct GAUSS3D *GarrayOrig;
  int g,Ngauss_new;
  double tole;

  tole = 0.0000000001;
  /*
  printf("#Delete_ZeroWeight_gdfs_of_GMM()\n");
  */

  GarrayOrig = (struct GAUSS3D*)malloc(sizeof(struct GAUSS3D)*(*Ngauss));  
 
  for (g=0;g<(*Ngauss);++g){ 
    Copy_GAUSS3D(&(GarrayOrig[g]),&(Garray[g]));
  }

  Ngauss_new = 0;
  for (g=0;g<(*Ngauss);++g){ 
    if (fabs(GarrayOrig[g].Weight)>tole){
      Copy_GAUSS3D(&(Garray[Ngauss_new]),&(GarrayOrig[g]));
      Garray[Ngauss_new].num = Ngauss_new;
      Ngauss_new += 1;
    }
  }

  printf("#Ngauss %d --> %d\n",*Ngauss,Ngauss_new);
  *Ngauss = Ngauss_new;
  /*
  for (g=0;g<(*Ngauss);++g){ 
    printf("#Gauss %d W %lf\n",Garray[g].num,Garray[g].Weight);
  }
   */
} /* end of Delete_ZeroWeight_gdfs_of_GMM(G) */



void Delete_Identical_gdfs_of_GMM(Ngauss,Garray)
 int *Ngauss;
 struct GAUSS3D *Garray;
{
  int i,j,g,h;
  char identical;
  double tole;

  printf("#Delete_Identical_gdfs_of_GMM()\n");
  if (*Ngauss>=1){
    tole = 0.0001;
    /** (1) Set Weight=0.0 for similar gdfs **/
    for (g=0;g<(*Ngauss);++g){
      if (Garray[g].Weight>0.0){
        for (h=g+1;h<(*Ngauss);++h){
          if (Garray[h].Weight>0.0){
            identical = 1;
            for (i=0;i<3;++i){
              if (fabs(Garray[g].M[i] - Garray[h].M[i])>tole){identical = 0;}
              for (j=0;j<3;++j){
                if (fabs(Garray[g].CovM[i][j] - Garray[h].CovM[i][j])>tole){identical = 0;}
              }
            }
            /* printf("#g %d h %d identical %d\n",g,h,identical); */
            if (identical==1){
              Garray[g].Weight += Garray[h].Weight; /* merging weights */
              Garray[h].Weight = 0.0;
            }
          }
        }
      }
    }

   /** (2) Delete gdfs with Weight=0.0. **/
   Delete_ZeroWeight_gdfs_of_GMM(Ngauss,Garray);
 }

} /* end of Delete_Identical_gdfs_of_GMM() */





void Copy_GAUSS3D(gnew,gorig)
  struct GAUSS3D *gnew;
  struct GAUSS3D *gorig;
{
  int i,j;

  gnew->num = gorig->num;
  gnew->det = gorig->det;
  gnew->Cons = gorig->Cons;
  gnew->Weight = gorig->Weight;

  for (i=0;i<3;++i){
    gnew->M[i] = gorig->M[i]; 
    gnew->PCvar[i] = gorig->PCvar[i]; 
    for (j=0;j<3;++j){
      gnew->CovM[i][j]   = gorig->CovM[i][j]; 
      gnew->iCovM[i][j]  = gorig->iCovM[i][j]; 
      gnew->PCaxis[i][j] = gorig->PCaxis[i][j]; 
    }
  }

} /* end of Copy_GAUSS3D() */


void Copy_GAUSS3D_Array(Ngauss,Gnew,Gorig)
  int Ngauss;
  struct GAUSS3D *Gnew;
  struct GAUSS3D *Gorig;
{
  int g;

  for (g=0;g<Ngauss;++g){
    Copy_GAUSS3D(&(Gnew[g]),&(Gorig[g]));
  }

} /* end of Copy_GAUSS3D_Array() */


