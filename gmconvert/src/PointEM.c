/*

 <PointEM.c>

 for EM fit of Gaussian Mixture Model for Point Data (PDB atom data)

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

/*** Functions (GLOBAL) ***/
void Initialize_Gauss_From_Randomly_Chosen_Atoms();
void Initialize_Gauss_From_Mol_Distribution();
void EM_optimize_GaussMix_For_Atoms();
void Malloc_Voxel_From_Covariance_Matrix();
void Cal_Mean_and_CovarMat_from_Atoms();
void Cal_Memberships_for_Atoms();
void Write_PDB_File_With_GMM_and_Memberships();
void Read_PDB_File_With_GMM_and_Memberships();
void Estimate_M_and_CovM_from_Memberships();

/*****************/
/*** FUNCTIONS ***/
/*****************/

void Initialize_Gauss_From_Randomly_Chosen_Atoms(Ahead,Ngauss,Garray)
 struct ATOM *Ahead;
 int Ngauss;
 struct GAUSS3D *Garray;
{
 int g,i,j,Natom,Natom_rand,Natom_rand2;
 double G[3],Variance;
 struct ATOM *an;
 int occupied,hit,minNgauss_Natom,round,atom_number;


 /** Calculate Gravity center and Variance **/
 G[0] = G[1] = G[2] = 0.0; 
 Natom = 0;
 an = Ahead; 
 while (an->next != NULL){ 
   an = an->next;
   G[0] += an->Pos[0];  
   G[1] += an->Pos[1];  
   G[2] += an->Pos[2];
   ++Natom;
 }
 G[0] /= Natom; G[1] /= Natom; G[2] /= Natom;

 Variance = 0.0;
 an = Ahead; 
 while (an->next != NULL){ 
   an = an->next;
   Variance += 
     (   (an->Pos[0]-G[0])*(an->Pos[0]-G[0]) 
       + (an->Pos[1]-G[1])*(an->Pos[1]-G[1]) 
       + (an->Pos[2]-G[2])*(an->Pos[2]-G[2]) );  
   an->mark = 0;  
 }

 Variance /= Natom;

 if (Ngauss<Natom){ minNgauss_Natom = Ngauss;}
             else { minNgauss_Natom = Natom;}

 /*
 printf("#Initialize_Gauss_From_Randomly_Chosen_Atoms(Ngauss %d minNgauss_Natom %d)\n",Ngauss,minNgauss_Natom);
 */

 /* Assign Garray[g].M[] by random chosen atom  */
 for (g=0;g<minNgauss_Natom;++g){
   occupied = 0; 
   round  = 0;
   do{
      Natom_rand = rand() % Natom;
      /* printf("#g %d Natom_rand %d Natom %d\n",g,Natom_rand,Natom);  */
      an = Ahead;
      hit = 0;
      atom_number = 0;
      while ((an->next != NULL)&&(hit==0)){
        an = an->next;
        if (atom_number == Natom_rand){ 
           hit = 1; 
           if (an->mark == 1){occupied = 1;}
           else{
             occupied = 0;
             for (i=0;i<3;++i){
               Garray[g].M[i] = an->Pos[i];     
             }
             an->mark = 1;
           }
         }
        atom_number += 1;
      } 
     round += 1;
   } while ((occupied == 1) && (round < 10*minNgauss_Natom));
  
   if (occupied==1){
     exit(1);
   }
 
   /* printf("Natom_rand_chosen %d\n",Natom_rand); */
    
    for (i=0;i<3;++i){
      for (j=0;j<3;++j){ 
         if (i==j) Garray[g].CovM[i][j] = Variance/Natom;
         else Garray[g].CovM[i][j] = 0.0; 
      }
    }
   Garray[g].Weight = 1.0/(double)Ngauss;
   Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));
   Garray[g].Cons = 1/(pow(2*M_PI,1.5)*sqrt(Garray[g].det));

 } /* g */

 /* For the case for Ngauss is larger than Natom */
 for (g=minNgauss_Natom;g<Ngauss;++g){
   Natom_rand   = rand() % minNgauss_Natom;
   Natom_rand2  = rand() % minNgauss_Natom;
   for (i=0;i<3;++i){
     Garray[g].M[i] = (Garray[Natom_rand].M[i] + Garray[Natom_rand].M[i])/2.0;
   }
   for (i=0;i<3;++i){
     for (j=0;j<3;++j){ 
        if (i==j) Garray[g].CovM[i][j] = Variance/Natom;
        else Garray[g].CovM[i][j] = 0.0; 
     }
   }
   Garray[g].Weight = 1.0/(double)Ngauss;
   Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));
 } 


} /* end of Initialize_Gauss_From_Randomly_Chosen_Atoms() */





void Initialize_Gauss_From_Mol_Distribution(Ngauss,Garray,Mean,CovM,Min,Max)
 int Ngauss;
 struct GAUSS3D *Garray;
 double Mean[3],CovM[3][3],Min[3],Max[3];
{
 int g,i,j;
 double AveSig;

 AveSig = (CovM[0][0]+CovM[1][1]+CovM[2][2])/3.0/Ngauss;

 for (g=0;g<Ngauss;++g)
 {
  for (i=0;i<3;++i) 
   Garray[g].M[i] = Min[i] + (Max[i]-Min[i]) * (double)rand()/RAND_MAX;

 /*
  for (i=0;i<3;++i)
   for (j=0;j<3;++j) Garray[g].CovM[i][j] = CovM[i][j];
 */

 /*
  for (i=0;i<3;++i)
   for (j=0;j<3;++j) 
   { if (i==j) Garray[g].CovM[i][j] = CovM[i][j];
     else Garray[g].CovM[i][j] = 0.0; }
 */

  for (i=0;i<3;++i){
    for (j=0;j<3;++j){ 
      if (i==j) Garray[g].CovM[i][j] = AveSig;
      else Garray[g].CovM[i][j] = 0.0;
    }
  }

  Garray[g].Weight = 1.0/(double)Ngauss;
  Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));
  Garray[g].Cons = 1/(pow(2*M_PI,1.5)*sqrt(Garray[g].det));

 } /* g */

} /* end of Initialize_Gauss_From_Mol_Distribution() */








void EM_optimize_GaussMix_For_Atoms(Nrepeat,Ngauss,Garray,Ahead,logLike_fin,ConvThre)
 int Nrepeat;
 int Ngauss;
 struct GAUSS3D *Garray;
 struct ATOM *Ahead;
 double *logLike_fin; /* Final Log Likelihood */
 double ConvThre;  /* Threshold value for logLike convergence */
{
 int r,g,a,i,j;
 int Natom;
 double **H; /* H[0..Ngauss-1][0..Natom-1] (membership/posterior probability) */ 
 struct ATOM *an;
 double sumG,Hg_sum_a,P,logLike,logLikePre,DiffLogLike;
 /* FILE *fp; */

 printf("#EM_optimize_GaussMix_For_Atoms()\n");
 /*** [1] Count Natom and malloc matrix H ***/
 Natom = 0;
 an = Ahead;
 while (an->next != NULL){ 
   an = an->next;
   an->num = Natom;
   ++Natom;
 }

 H = (double **)malloc(sizeof(double*)*Ngauss);
 for (g=0;g<Ngauss;++g){ H[g] = (double *)malloc(sizeof(double)*Natom);}


  logLikePre = 0.0;
  an = Ahead;
  while (an->next != NULL){
    an = an->next;
    P = 0.0;
    for (g=0;g<Ngauss;++g){
      if (Garray[g].Weight>0.0){
        P += Garray[g].Weight * Gaussian3D_Value(&(Garray[g]),an->Pos);
      }
    }
    logLikePre += log(P); 
  }
 printf("#logLikeInit %lf %e\n",logLikePre,logLikePre);

 /*** [2] Iteration(r) for EM  ***/
 
 r = 0;
 do{
  /** (E-step) Calculate H[g][a] (posterior probablility )**/ 
  an = Ahead;
  while (an->next != NULL){ 
    an = an->next;
    sumG = 0.0;
    for (g=0;g<Ngauss;++g){
     if (Garray[g].Weight>0.0){ 
       H[g][an->num] = Garray[g].Weight * Gaussian3D_Value(&(Garray[g]),an->Pos);
       sumG += H[g][an->num];
      }
    } 
    for (g=0;g<Ngauss;++g){H[g][an->num] /= sumG;}
  } 

  /** (M-step) Update Weight[g],Mean[g],CovM[g] **/ 
  for (g=0;g<Ngauss;++g){ 

    /* update Garray[g].Weight */
    Hg_sum_a =  0.0;
    for (a=0;a<Natom;++a){Hg_sum_a += H[g][a];}
    Garray[g].Weight = Hg_sum_a/Natom;

    for (i=0;i<3;++i){Garray[g].M[i] = 0.0;}
    
    /* update Garray[g].M[] */
    an = Ahead;
    while (an->next != NULL){ 
      an = an->next; 
      for (i=0;i<3;++i){ Garray[g].M[i]  +=  H[g][an->num] * an->Pos[i];}
    }
    for (i=0;i<3;++i){ Garray[g].M[i]  /=  Hg_sum_a;}

    /* update Garray[g].CovM[] */
    for (i=0;i<3;++i){ 
      for (j=0;j<3;++j){ Garray[g].CovM[i][j] = 0.0;}
    }

    an = Ahead;
    while (an->next != NULL){ 
      an = an->next; 
      for (i=0;i<3;++i){ 
        for (j=0;j<3;++j){
          Garray[g].CovM[i][j] +=  
            H[g][an->num] * (an->Pos[i] - Garray[g].M[i]) * (an->Pos[j] - Garray[g].M[j]);
        }
      }
    }

    for (i=0;i<3;++i){
      for (j=0;j<3;++j){Garray[g].CovM[i][j] /= Hg_sum_a;}
    }

   Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));
   Garray[g].Cons = 1/(pow(2*M_PI,1.5)*sqrt(Garray[g].det));

  } /* g */

  /* Calculate Likelihood */
  logLike = 0.0;
  an = Ahead;
  while (an->next != NULL){
    an = an->next;
    P = 0.0;
    for (g=0;g<Ngauss;++g){
      if (Garray[g].Weight>0.0){
        P += Garray[g].Weight * Gaussian3D_Value(&(Garray[g]),an->Pos);
      }
    }
    logLike += log(P); 
  }
 
  DiffLogLike = logLike - logLikePre;
  printf("Step %d logLike %lf Diff %lf\n",r,logLike,DiffLogLike);
  logLikePre = logLike;
  *logLike_fin = logLike;
  ++r;
 
 } while ((r<Nrepeat)&&(DiffLogLike>ConvThre));

 printf("#Converged after %d steps (logLike_fin:%e)\n",r,*logLike_fin);

 for (g=0;g<Ngauss;++g){ free(H[g]); }
 free(H);

} /* end of EM_optimize_GaussMix_For_Atoms() */










void Malloc_Voxel_From_Covariance_Matrix(vox,Mean,CovM)
 struct VOXEL *vox;
 double Mean[3];
 double CovM[3][3];
{
 double SD,Mag;
 int i;

 Mag = 3;
 for (i=0;i<3;++i){
  SD = sqrt(CovM[i][i]);
  vox->OrigPos[i] = Mean[i] - Mag * SD;
  printf("MIN[%d] %f MAX %f\n",i,vox->OrigPos[i],Mean[i]+Mag*SD);
  vox->N[i]   = (int)ceil((2*Mag * SD)/vox->grid_width);
 }
 
 Malloc_Voxel(vox,vox->N[0],vox->N[1],vox->N[2]);

} /* end of Malloc_Voxel_From_Covariance_Matrix() */







void Cal_Mean_and_CovarMat_from_Atoms(Ahead,M,Cov,Min,Max)
 struct ATOM *Ahead;
 double  M[3];
 double Cov[3][3];
 double Min[3],Max[3];
{
 struct ATOM *an;
 int Natom;
 int i,j;
 double Dvec[3];

 /* Caliculate Mean,Min,Max */
 an = Ahead;
 for (i=0;i<3;++i) M[i] = 0.0;
 Natom = 0;
 while (an->next != NULL){
  an = an->next;
  for (i=0;i<3;++i) M[i] += an->Pos[i];
  
  if (Natom==0) {for (i=0;i<3;++i) Min[i] = Max[i] = an->Pos[i];}
  else{
   for (i=0;i<3;++i){ 
     if (an->Pos[i]<Min[i]) Min[i] = an->Pos[i];
     if (an->Pos[i]>Max[i]) Max[i] = an->Pos[i];
   }
 } 

   ++Natom;
 } /* an */

 for (i=0;i<3;++i){M[i] /= Natom;}

 /* Caliculate Covariance Matrix */

 an = Ahead;
 for (i=0;i<3;++i){ 
   for (j=0;j<3;++j){Cov[i][j] = 0.0;}
 }

 while (an->next != NULL){
   an = an->next;
   for (i=0;i<3;++i){Dvec[i] = an->Pos[i] - M[i];}
  
   for (i=0;i<3;++i){
     for (j=i;j<3;++j){ 
       Cov[i][j] += Dvec[i]*Dvec[j]/Natom;
       Cov[j][i] = Cov[i][j];
     }
   }
 } /* an */


} /* end of Cal_Mean_and_CovarMat_from_Atoms() */



void Cal_Memberships_for_Atoms(AtomMember,Ahead,Ngauss,Garray)
 struct MATRIX *AtomMember; 
 struct ATOM *Ahead;
 int    Ngauss;
 struct GAUSS3D *Garray;
{
 struct ATOM *an;
 int g,a;
 double sum;

 printf("#Cal_Memberships_for_Atoms(Nrow %d Ncol %d Ngauss %d)\n",AtomMember->Nrow,AtomMember->Ncol,Ngauss);

 /*
 an = Ahead;
  while (an->next != NULL){
    an = an->next;
    sumG = 0.0;
    for (g=0;g<Ngauss;++g){
     if (Garray[g].Weight>0.0){
       H[g][an->num] = Garray[g].Weight * Gaussian3D_Value(&(Garray[g]),an->Pos);
       sumG += H[g][an->num];
      }
    }
    for (g=0;g<Ngauss;++g){H[g][an->num] /= sumG;}
  } 
 */




 an = Ahead;
 a = 0;
 while (an->next != NULL){
   an = an->next;
   /** calculate memberships H[] **/
   sum = 0.0;
   if ((a>AtomMember->Nrow)||(an->num > AtomMember->Nrow)){
     printf("#ERROR:Natom %d %d is over AtomMember->Nrow %d.\n",a,an->num,AtomMember->Nrow);
     exit(1);
   }

   for (g=0;g<Ngauss;++g){   
     if (Garray[g].Weight>0.0){
       /* AtomMember->m[a][g] = Garray[g].Weight*Gaussian3D_Value(&(Garray[g]),an->Pos); */
       AtomMember->m[an->num][g] = Garray[g].Weight*Gaussian3D_Value(&(Garray[g]),an->Pos); 
       sum += AtomMember->m[an->num][g];
     }
   }
   for (g=0;g<Ngauss;++g){ AtomMember->m[an->num][g] /= sum;} 
   a += 1;
 }

 
} /* end Cal_Memberships_for_Atoms() */




void Write_PDB_File_With_GMM_and_Memberships(fname,Ahead,Ngauss, Garray,AtomMember,comment,command)
 char *fname;
 struct ATOM *Ahead;
 int Ngauss;
 struct GAUSS3D *Garray;
 struct MATRIX *AtomMember;
 char *comment;
 char *command;
{
 FILE *fp;
 struct ATOM *an;
 int a,g,max_g;
 char line[1024],subline[100];
 double maxH;

 printf("#Write_PDB_File_With_Memberships() -> \"%s\"\n",fname);
 if (fname[0]=='-') fp = stdout;
 else{
   fp= fopen(fname,"w");
  if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1);}
 }

 sprintf(line,"FILENAME \"%s\"",fname);
 Get_Part_Of_Line(subline,line,0,40);
 fprintf(fp,"HEADER    %-40s%s   0XXX      0XXX   1\n",subline,Get_Date_String_PDB());
 fprintf(fp,"REMARK    DATE %s\n",Get_Date_String());
 if (command[0] != '\0'){ fprintf(fp,"REMARK    COMMAND %s\n",command);}
 fprintf(fp,"REMARK    NOTE: THIS FILE CONTAINS ATOMIC COODINATES WITH MEMBERSHIP VALUES AND GMM MODEL\n");
 fprintf(fp,"REMARK #   %s\n",comment);
 
 /** Write GMM **/
 fprintf(fp,"REMARK NGAUSS %d\n",Ngauss);
 for (g=0;g<Ngauss;++g){
   fprintf(fp,"REMARK GAUSS%4d W %.10lf\n",g+1,Garray[g].Weight);
   fprintf(fp,"REMARK GAUSS%4d M %lf %lf %lf\n",g+1,Garray[g].M[0], Garray[g].M[1], Garray[g].M[2]);
   fprintf(fp,"REMARK GAUSS%4d CovM  xx %15.10lf xy %15.10lf xz %15.10lf\n",
        g+1,Garray[g].CovM[0][0],Garray[g].CovM[0][1],Garray[g].CovM[0][2]);
   fprintf(fp,"REMARK GAUSS%4d CovM  yy %15.10lf yz %15.10lf zz %15.10lf\n",
        g+1,Garray[g].CovM[1][1],Garray[g].CovM[1][2],Garray[g].CovM[2][2]);

 }
 fprintf(fp,"REMARK    Occupancy [55-60] : Radius of the atom.\n");
 fprintf(fp,"REMARK    tFactor   [61-65] : GDF number with the highest membership value.\n");

 fprintf(fp,"REMARK                                   MEMBERSHIP VALUE FOR GDF:");
 for (g=0;g<Ngauss;++g){
   fprintf(fp," [%3d]",g+1);
 }
 fprintf(fp,"\n"); 

/*** Write ATOM/HETATM ***/
 an = Ahead;
 a = 0;
 while (an->next != NULL){
   an = an->next;
   /** calculate max_g for maximum member ship H[] **/
   maxH = 0.0;
   max_g = 0;
   for (g=0;g<Ngauss;++g){   
     if ((g==0)||(AtomMember->m[a][g]>maxH)){ max_g = g+1; maxH = AtomMember->m[a][g];}
   }

   /** write ATOM/HETATM line **/
        if (an->AHtype=='A')  fprintf(fp,"ATOM  ");
   else if (an->AHtype=='H')  fprintf(fp,"HETATM");
   else fprintf(fp,"ATOM? ");

   fprintf(fp,"%5s %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f",
        an->Anum,an->Atom,an->Resi,an->Chain,an->Rnum,
        an->Pos[0],an->Pos[1],an->Pos[2],an->R,(float)max_g);
   for (g=0;g<Ngauss;++g){
     /* fprintf(fp," %5.3f",AtomMember->m[a][g]); */
     fprintf(fp," %5.3lf",AtomMember->m[a][g]); 
   }
   fprintf(fp,"\n");
   a += 1;
 }

 fprintf(fp,"TER\n");

 if (fp!=stdout) fclose(fp);
 
} /* end of Write_PDB_File_With_GMM_and_Memberships() */









void Read_PDB_File_With_GMM_and_Memberships(ifname,Ngauss,Garray,AtomMember)
  char  *ifname;
  int   Ngauss;
  struct GAUSS3D *Garray; /* already malloced (Garray[0..Ngauss-1]) */
  struct MATRIX  *AtomMember; 
{
 FILE *fp;
 int  g,h,gnum,gnum0,w,L;
 int Natom;
 char line[512],buff[512],word[10][100];
 int  Wsta[100],Wend[100],Nword;
 double Weight_all;


 fp = fopen(ifname,"r");
 if (fp==NULL) {printf("#ERROR:Can't open mpdbfile \"%s\"\n",ifname); exit(1);}

 printf("#Read_PDB_File_With_GMM_and_Memberships(\"%s\")\n",ifname);

/*
REMARK NGAUSS 4
REMARK GAUSS   1 W 0.3126277222
REMARK GAUSS   2 W 0.2191508529
:
REMARK GAUSS   4 W 0.2112233392

REMARK    Occupancy [55-60] : Radius of the atom.
REMARK    tFactor   [61-65] : GDF number with the highest membership value.
REMARK                                   MEMBERSHIP VALUE FOR GDF: [  1] [  2] [  3] [  4]
          1         2         3         4         5         6         7         8         9
01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM      1  N   ALA A   2      63.130 -67.331 -18.137  1.65  4.00 0.000 0.000 0.000 1.000
ATOM      2  CA  ALA A   2      63.914 -67.078 -16.894  1.87  4.00 0.000 0.000 0.000 1.000
ATOM      3  C   ALA A   2      63.317 -67.813 -15.708  1.76  4.00 0.000 0.000 0.000 1.000
ATOM      4  O   ALA A   2      62.179 -68.283 -15.768  1.40  4.00 0.000 0.000 0.000 1.000
:
ATOM    248  N   ARG A  36      36.608 -36.818 -11.210  1.65  2.00 0.000 0.560 0.424 0.016
ATOM    249  CA  ARG A  36      37.796 -35.940 -11.222  1.87  3.00 0.000 0.430 0.557 0.014
:
*/
 g = -1; gnum0 = -1; 
 Weight_all = 0.0;
 Natom = 0;
 while (feof(fp)==0){
   line[0] = '\0';
   fgets(line,511,fp);

   L = strlen(line);
   if ((L>0)&&(line[L-1]=='\n')) line[L-1] = '\0';

   /*** READ GMM INFORMATION ***/
   if (strncmp(line,"REMARK GAUSS",12)==0){
     Get_Part_Of_Line(buff,line,12,512);
     Split_to_Words(buff,' ',&Nword,Wsta,Wend,100);
     for (w=0;w<Nword;++w) Get_Part_Of_Line(word[w],buff,Wsta[w],Wend[w]);
     gnum = atoi(word[0]);
     if (gnum0!=gnum) { g += 1;}
     if (g > Ngauss){
       printf("#ERROR:Ngauss %d is over Ngauss_malloc %d.\n",g,Ngauss);
       exit(1);
     }
     gnum0 = gnum;
     /* printf("g %d gnum %d '%s'\n",g,gnum,line); */
         if (strncmp(word[1],"W",1)==0){
      Garray[g].Weight = atof(word[2]);
      Weight_all += Garray[g].Weight;
    }
   } /* REMARK GAUSS */

   /*** READ ATOM/HETATM INFORMATION ***/
   if ((strncmp(line,"ATOM",4)==0)||((PAR.HETtype=='T')&&(strncmp(line,"HETATM",6)==0))){
     for (h=0;h<Ngauss;++h){
       Get_Part_Of_Line(buff,line,66+6*h,66+6*h+6);
       AtomMember->m[Natom][h] = atof(buff);
     }
     /*
     printf("%s::",line);
     for (h=0;h<Ngauss;++h){ printf(" %5.3f",AtomMember->m[Natom][h]);}
     printf("\n");
     */
     Natom += 1;
     if (Natom>AtomMember->Nrow){
       printf("#ERROR:Natom %d is over premalloced AtomMember->Nrow %d\n",Natom,AtomMember->Nrow);
       exit(1);
     }
   }


 } /* while */
 
 for (g=0;g<Ngauss;++g){
   printf("#GDF[%d] Weight %f\n",g+1,Garray[g].Weight);
 }
 fclose(fp);

} /* end of Read_PDB_File_With_GMM_and_Memberships() */


void Estimate_M_and_CovM_from_Memberships(Ngauss,Garray,Ahead,AtomMember)
  int   Ngauss;
  struct GAUSS3D *Garray; /* already malloced (Garray[0..Ngauss-1]) */
  struct ATOM *Ahead;
  struct MATRIX  *AtomMember; 
{
  int i,j,a,g;
  struct ATOM *an;
  double sum_member_a; 
  printf("#Estimate_M_and_CovM_from_Memberships(Ngauss,Garray,Ahead,AtomMember)\n");
  /** Initialize Garray[g].M[] and CovM[][] **/ 
  for (g=0;g<Ngauss;++g){
    for (i=0;i<3;++i){
      Garray[g].M[i] = 0.0;
      for (j=0;j<3;++j){
        Garray[g].CovM[i][j] = 0.0;
      }
    }  
  }

  for (g=0;g<Ngauss;++g){
    /** cal sum_member **/
    an = Ahead;
    a = 0;
    sum_member_a = 0.0;
    while (an->next != NULL){ 
      an = an->next; 
      /* printf("#a %d g %d member %f\n",a,g,AtomMember->m[a][g]); */
      sum_member_a += AtomMember->m[a][g];
      ++a;
    }
    printf("#GAUSS[%d] sum_member %f\n",g+1,sum_member_a); 
    
   /* cal Garray[g].M[] */
    an = Ahead;
    a = 0;
    while (an->next != NULL){ 
      an = an->next; 
      for (i=0;i<3;++i){Garray[g].M[i]  +=  AtomMember->m[a][g] * an->Pos[i];}
      ++a;
    }
    for (i=0;i<3;++i){Garray[g].M[i]  /=  sum_member_a;}


    /* cal Garray[g].CovM[] */
    an = Ahead;
    a = 0;
    while (an->next != NULL){ 
      an = an->next; 
      for (i=0;i<3;++i){ 
        for (j=0;j<3;++j){
          Garray[g].CovM[i][j] +=  
            AtomMember->m[a][g] * (an->Pos[i]-Garray[g].M[i]) * (an->Pos[j]-Garray[g].M[j]);
        }
      }
      ++a;
    }

    for (i=0;i<3;++i){
     for (j=0;j<3;++j){Garray[g].CovM[i][j] /= sum_member_a;}
    }

   Cal_Inverse_Matrix3D_by_Cramer_Rule(Garray[g].iCovM, Garray[g].CovM, &(Garray[g].det));
   Garray[g].Cons = 1/(pow(2*M_PI,1.5)*sqrt(Garray[g].det));

  } /* g */
 

} /* end of Estimate_M_and_CovM_from_Memberships() */
