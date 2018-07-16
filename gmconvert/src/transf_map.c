/*
  <transf_map.c>

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "Voxel.h"
#include "MapCCP4.h"

static char LastModDate[] = "Feb 26, 2016";

#define MAX_FILE_LEN 1024

char MMCIF_DIR[MAX_FILE_LEN];
char PDB_GMM_DIR[MAX_FILE_LEN];
char EMDB_GMM_DIR[MAX_FILE_LEN];
char EMDB_DIR[MAX_FILE_LEN];
char TMPOUT_DIR[MAX_FILE_LEN];
char OMOKAGE_USERDATA_DIR[MAX_FILE_LEN];
char *ARG_STRING;

static void split_into_words();
static void get_part_of_line();
static void split_into_two_words();
static int read_environment_file();
static void get_Rmat_and_Tvec_from_strings();
static void get_Rmat_from_Raxis_str();
static void get_min_max_xyz_by_transformation();
static void setup_density_for_transformed_map();
static void transform_Rx_plus_tvec_3D();
static void inv_transform_Rx_plus_tvec_3D();

void split_into_words(str,splsym,Nword,Wsta,Wend,Nwordmax)
 char *str;          /* Input String */
 char splsym;        /* Symbol for split */
 int *Nword;         /* Number of words  */
 int Wsta[];         /* Start point of str for a word */
 int Wend[];         /* End point of str for a word */
 int Nwordmax;       /* Maxium number of Nword  */
{
 /* [Example]
     str = "//abc/d/ef//ghi//"
     splsym = '/'
      -->
     Nword 4
    (Wsta,Wend) = {(0,2), (4,4), (6,7), (10,12)}
 */
 int i,L;
 L = strlen(str);
 *Nword = 0; i = 0;
 while ((i<L)&&(*Nword < Nwordmax)){
  if (str[i]!=splsym){
     Wsta[*Nword] = i;
     while ((str[i]!=splsym)&&(i<=(L-1))) { ++i; }
     Wend[*Nword] = i-1;
     ++(*Nword);
   }
  ++i;
 }
} /* end of split_into_words() */


void get_part_of_line(part,line,s,e)
  char *part;
  char *line;
  int  s,e;
{
 int i,E,L;
 L = strlen(line)-1;
 /* printf("#Get_Part_Of_Line(line:'%s',s:%d e:%d L:%d)\n",line,s,e,L); fflush(stdout); */
 if (line[L] == '\n') L -= 1;
 if (s<0) s = 0;
 if (e>L) E = L; else E = e;
 for (i=s;i<=E;++i) part[i-s] = line[i];
 part[E-s+1] = '\0';

} /* end of get_part_of_line() */


void split_into_two_words(line,start_pos,splsym,first_word,second_word)
  char *line;
  int  start_pos;  /* start position of line[]. the default should be 0 */
  char splsym;     /* symbol for spliting. for example, '.', '\t',' ',.. */
  char *first_word;
  char *second_word;
{
 /*
    if line[] = '_audit_conform.dict_name'
     splsym = '.', start_pos = 1.
       ==> first_word[]  = 'audit_conform'
           second_word[] = 'dict_name'
 
     if (line[] = '_audit_conform.dict_name       mmcif_pdbx.dic'
        splsym = ' ', start_pos = 1.

     ==> first_word[]  =  'audit_conform.dict_name'
         second_word[] = 'mmcif_pdbx.dic'
 */

 int i,j,L;
  char findspl;

  L = strlen(line);
  i = start_pos;

  findspl = 0;
  while ((i<L)&&(findspl==0)){
    if (line[i]==splsym){
      findspl = 1;
      for (j=start_pos;j<i;++j){
        first_word[j-start_pos] = line[j];
      }
      first_word[i-start_pos] = '\0';
      while ((line[i]==splsym) && (i<L)){
        i += 1;
      }

      for (j=i;j<L;++j){
        second_word[j-i] = line[j];
      }
      second_word[L-i] = '\0';
    }
   ++i;
  }

  if (findspl==0){
    for (j=start_pos;j<L;++j){
      first_word[j-start_pos] = line[j];
    }
    first_word[L] = '\0';
    second_word[0] = '\0';
  }

} /* end of split_into_two_words() */







int read_environment_file(ienvfile)
  char *ienvfile;
{
/*
 ## FORMAT EXAMPLE ##
 # BASE_DIR  /usr/people/takawaba/work/STRMAT/Matras10
 # SCORE_DIR /usr/people/takawaba/work/STRMAT/Matras10/ROM-00Oct6
 # BSSP_DIR  /lab/DB/BSSP
 # PDB_DIR   /lab/PDB
 # EMDB_DIR      /kf1/PDBj/ftp/pdbj/pub/emdb/structures
 # OMOKAGE_USERDATA_DIR /home/pdbj/omokage/userdata
 # TMPOUT_DIR    /var/www/html/gmfit/TMPOUT

*/
  FILE *fp;
  char line[1000],key[1000],value[1000];
  int L;

  /* printf("#read_environment_file('%s')\n",ienvfile); */

  MMCIF_DIR[0]    = '\0'; 
  PDB_GMM_DIR[0]  = '\0'; 
  EMDB_GMM_DIR[0] = '\0'; 

  fp = fopen(ienvfile,"r");
  if (fp==NULL){
    /* printf("#WARNING: Can't open ienvfile '%s'.\n",ienvfile); */
    return(0);
  }

  while (feof(fp)==0){
    line[0] = '\0';
    fgets(line,1000,fp);
    L = strlen(line);
    /* printf("%s",line); */
    if (line[L-1]=='\n'){ 
      line[L-1] = '\0';
      L = strlen(line);
    }

    if ((L>5) && (line[0] != '#')){
       split_into_two_words(line,0,' ',key,value);  
       if (strcmp(key,"MMCIF_DIR")==0){ sprintf(MMCIF_DIR,"%s",value);}
       if (strcmp(key,"PDB_GMM_DIR")==0){ sprintf(PDB_GMM_DIR,"%s",value);}
       if (strcmp(key,"EMDB_GMM_DIR")==0){ sprintf(EMDB_GMM_DIR,"%s",value);}
       if (strcmp(key,"EMDB_DIR")==0){ sprintf(EMDB_DIR,"%s",value);}
       if (strcmp(key,"TMPOUT_DIR")==0){ sprintf(TMPOUT_DIR,"%s",value);}
       if (strcmp(key,"OMOKAGE_USERDATA_DIR")==0){ sprintf(OMOKAGE_USERDATA_DIR,"%s",value);}
    }
  }

  fclose(fp);
  return(1);
} /* end of read_environment_file() */






void get_Rmat_and_Tvec_from_strings(Rmat,Tvec,Rmat_str,Tvec_str)
  float Rmat[3][3],Tvec[3];
  char *Rmat_str, *Tvec_str;
{
  int Nword,Wsta[100],Wend[100],i,j;
  char value[100];

  split_into_words(Tvec_str,':',&Nword,Wsta,Wend,100);
  for (i=0;i<3;++i){
    get_part_of_line(value,Tvec_str,Wsta[i],Wend[i]);
    Tvec[i] = atof(value);
  }

  split_into_words(Rmat_str,':',&Nword,Wsta,Wend,100);

  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      get_part_of_line(value,Rmat_str,Wsta[3*i+j],Wend[3*i+j]);
      Rmat[i][j] = atof(value);
    }
  }

} /* end of get_Rmat_and_Tvec_from_strings() */


void get_Rmat_from_Raxis_str(Rmat,Raxis_str)
  float Rmat[3][3];
  char *Raxis_str;
{
  int Nword,Wsta[100],Wend[100],i,j;
  char value[100];
  float u[3],ulen; /* axis vector */
  float rangle,sin_t,cos_t; /* rotational angle (radian) */
  float U[3][3];
  float V[3][3]; /* U x U */
  split_into_words(Raxis_str,':',&Nword,Wsta,Wend,100);
  ulen = 0.0;
  for (i=0;i<3;++i){
    get_part_of_line(value,Raxis_str,Wsta[i],Wend[i]);
    u[i] = atof(value);
    ulen += u[i]*u[i];
  }

  for (i=0;i<3;++i){ u[i] /= sqrt(ulen);}

  get_part_of_line(value,Raxis_str,Wsta[3],Wend[3]);
  rangle = atof(value)/180.0*M_PI;

  /* printf("#axis %f %f %f rangle %f\n",u[0],u[1],u[2],rangle); */
  sin_t = sin(rangle);
  cos_t = cos(rangle);

  U[0][0] =   0.0; U[0][1] = -u[2]; U[0][2] =  u[1];
  U[1][0] =  u[2]; U[1][1] =   0.0; U[1][2] = -u[0];
  U[2][0] = -u[1]; U[2][1] =  u[0]; U[2][2] =  0.0;

  V[0][0] = u[0]*u[0]-1.0; V[0][1] = u[0]*u[1];     V[0][2] = u[0]*u[2];
  V[1][0] = u[0]*u[1];     V[1][1] = u[1]*u[1]-1.0; V[0][2] = u[1]*u[2];
  V[2][0] = u[0]*u[2];     V[2][1] = u[1]*u[2];     V[2][2] = u[2]*u[2]-1.0;

  for (i=0;i<3;++i){ 
    for (j=0;j<3;++j){ 
      Rmat[i][j]  =  sin_t * U[i][j] + (1.0-cos_t)*V[i][j];
      if (i==j) {Rmat[i][j] += 1.0;}
    }
  }
 /*
  printf("Rmat0 %+8.5f %+8.5f %+8.5f\n",Rmat[0][0],Rmat[0][1],Rmat[0][2]);
  printf("Rmat1 %+8.5f %+8.5f %+8.5f\n",Rmat[1][0],Rmat[1][1],Rmat[1][2]);
  printf("Rmat2 %+8.5f %+8.5f %+8.5f\n",Rmat[2][0],Rmat[2][1],Rmat[2][2]);
 */

} /* end of get_Rmat_and_Tvec_from_strings() */






void get_min_max_xyz_by_transformation(vox,R,T,min,max)
  struct VOXEL *vox;
  float R[3][3],T[3];
  float min[3],max[3];
{
  int i,j,k,n;
  float ver_ori[8][3];
  float ver_tra[8][3];
  char init[3];

  for (i=0;i<2;++i){
    for (j=0;j<2;++j){
      for (k=0;k<2;++k){
        n = 4*i + 2*j + k;
        ver_ori[n][0] = vox->OrigPos[0] + vox->grid_width * vox->N[0] * i;
        ver_ori[n][1] = vox->OrigPos[1] + vox->grid_width * vox->N[1] * j;
        ver_ori[n][2] = vox->OrigPos[2] + vox->grid_width * vox->N[2] * k;
      }
    }
  }

 init[0] = init[1] = init[2] = 1;
 
 for (n=0;n<8;++n){
   transform_Rx_plus_tvec_3D(ver_tra[n],ver_ori[n],R,T);

   for (i=0;i<3;++i){
     if ((init[i]==1) || (ver_tra[n][i] < min[i])){ min[i] = ver_tra[n][i];}
     if ((init[i]==1) || (ver_tra[n][i] > max[i])){ max[i] = ver_tra[n][i];}
     init[i] = 0;
   }
 } 


} /* end of get_min_max_xyz_by_transformation() */


void setup_density_for_transformed_map(Vnew,Vori,R,tvec,Algorithm)
  struct VOXEL *Vnew,*Vori;
  float R[3][3],tvec[3];
  char  Algorithm; /* 'N'aive,'A'verage 'L'inear interpolation */
{
  int i,j,k,io,jo,ko,i0,j0,k0,a,b,c;
  float pnew[3],pori[3],maxori[3],min_density,sumWeight,Weight;
  float delta[3]; /* dx,dy,dz */
  float wdelta[3][2]; /* w(dx),w(1-dx),w(dy),w(1-dy),w(dz),w(1-dz) */
  int Nnei,n;

/*
  printf("#setup_density_for_transformed_map(Vnew,Vori,R,tvec,Algorithm:%c)\n",Algorithm);
*/
 
 /* find min_density */
  min_density = Vori->dat[0][0][0];
  for (i=0;i<Vori->N[0];++i){
    for (j=0;j<Vori->N[1];++j){
      for (k=0;k<Vori->N[1];++k){
        if (Vori->dat[i][j][k]<min_density){
          min_density = Vori->dat[i][j][k];
        }
      }
    } 
  }

  for (i=0;i<3;++i){ maxori[i] = Vori->OrigPos[i] + Vori->grid_width * Vori->N[i]; }

 /* scan each voxel in Vnew */
  for (i=0;i<Vnew->N[0];++i){
    pnew[0] = Vnew->OrigPos[0] + Vnew->grid_width * i;
    for (j=0;j<Vnew->N[1];++j){
      pnew[1] = Vnew->OrigPos[1] + Vnew->grid_width * j;
      for (k=0;k<Vnew->N[2];++k){
        pnew[2] = Vnew->OrigPos[2] + Vnew->grid_width * k;
        Vnew->dat[i][j][k] = min_density; 
        inv_transform_Rx_plus_tvec_3D(pori,pnew,R,tvec);
        if (   (Vori->OrigPos[0]<=pori[0])&&(pori[0]<=maxori[0])
            && (Vori->OrigPos[1]<=pori[1])&&(pori[1]<=maxori[1])
            && (Vori->OrigPos[2]<=pori[2])&&(pori[2]<=maxori[2])){


          if (Algorithm == 'N'){  /* Naive (Shisha Gonyu ) */
            io = (int)((pori[0] - Vori->OrigPos[0])/Vori->grid_width);
            jo = (int)((pori[1] - Vori->OrigPos[1])/Vori->grid_width);
            ko = (int)((pori[2] - Vori->OrigPos[2])/Vori->grid_width);
            if ((0<=io)&&(io<Vori->N[0]) && (0<=jo)&&(jo<Vori->N[1]) && (0<=ko)&&(ko<Vori->N[2])){
              Vnew->dat[i][j][k] = Vori->dat[io][jo][ko];
            }
           }

          else if (Algorithm == 'A'){ /* Average */
            i0 = (int)floor((pori[0] - Vori->OrigPos[0])/Vori->grid_width);
            j0 = (int)floor((pori[1] - Vori->OrigPos[1])/Vori->grid_width);
            k0 = (int)floor((pori[2] - Vori->OrigPos[2])/Vori->grid_width);
            Vnew->dat[i][j][k] = 0.0;
            Nnei = 0;
            for (n=0;n<8;++n){
              if (n>=4)          {a = 1;} else {a = 0;}
              if ((n-4*a)>=2)    {b = 1;} else {b = 0;}
              if ((n-4*a-2*b)>=1){c = 1;} else {c = 0;}
              io = i0 + a;
              jo = j0 + b;
              ko = k0 + c;
              if ((0<=io)&&(io<Vori->N[0]) && (0<=jo)&&(jo<Vori->N[1]) && (0<=ko)&&(ko<Vori->N[2])){
                Vnew->dat[i][j][k] += Vori->dat[io][jo][ko];
                Nnei += 1;
              }
            }
            if (Nnei > 0) {Vnew->dat[i][j][k] /= Nnei;}
            /* printf("==> %f\n",Vnew->dat[i][j][k]);  */
          }
          else if (Algorithm == 'L'){  /* Linear Interpolation */
/*
 >> Interpolation of 2D image << 
 01------11 
  |       | 
  |----p  |
  dy   |  |
  |    |  |
 00-dx---10

  f(p) = (w(dx)w(dy)f00  + w(1-dx)w(dy)f10 + w(dx)w(1-dy)f01 + w(1-dx)w(1-dy)f11)/Z
  Z    =  w(dx)w(dy)  + w(1-dx)w(dy) + w(dx)w(1-dy) + w(1-dx)w(1-dy)
  w(x) = 1.0 -x

  dx,dy,dz ==> delta[3]
  w(dx),w(1-dx),w(dy),w(1-dy),w(dz),w(1-dz) ==> wdelta[3][2]
*/
            i0 = (int)floor((pori[0] - Vori->OrigPos[0])/Vori->grid_width);
            j0 = (int)floor((pori[1] - Vori->OrigPos[1])/Vori->grid_width);
            k0 = (int)floor((pori[2] - Vori->OrigPos[2])/Vori->grid_width);
            delta[0] = (pori[0] - (Vori->OrigPos[0] + Vori->grid_width * i0))/Vori->grid_width;
            delta[1] = (pori[1] - (Vori->OrigPos[1] + Vori->grid_width * j0))/Vori->grid_width;
            delta[2] = (pori[2] - (Vori->OrigPos[2] + Vori->grid_width * k0))/Vori->grid_width;
            for (a=0;a<3;++a){
              wdelta[a][0] = 1.0-delta[a]; 
              wdelta[a][1] = delta[a];
            }

            sumWeight = 0.0;
            Vnew->dat[i][j][k] = 0.0;
            for (a=0;a<=1;++a){
              for (b=0;b<=1;++b){
                for (c=0;c<=1;++c){
                  io = i0 + a;
                  jo = j0 + b;
                  ko = k0 + c;
                  if ((0<=io)&&(io<Vori->N[0]) && (0<=jo)&&(jo<Vori->N[1]) && (0<=ko)&&(ko<Vori->N[2])){
                    Weight = wdelta[0][a] * wdelta[1][b] *wdelta[2][c];
                    Vnew->dat[i][j][k] += Weight * Vori->dat[io][jo][ko];
                    sumWeight += Weight;
                  }
                }
              }
           } 
           if (sumWeight > 0.0) {Vnew->dat[i][j][k] /= sumWeight;}
            /* printf("==> %f\n",Vnew->dat[i][j][k]);  */
          }
        }

      }
    }
  }

} /* end of setup_density_for_transformed_map() */



void transform_Rx_plus_tvec_3D(xT,xO,R,tvec)    /* xT = R*xO + t */
 float xT[3],xO[3];
 float R[3][3],tvec[3];
{
 int i,j;

 for (i=0;i<3;++i){
   xT[i] = tvec[i];
   for (j=0;j<3;++j){xT[i] += R[i][j]*xO[j];}
 }

} /* end of transform_Rx_plus_tvec_3D() */


void inv_transform_Rx_plus_tvec_3D(xT,xO,R,tvec)    /* xT = R^T (x-t) */
 float xT[3],xO[3];
 float R[3][3],tvec[3];
{
 int i,j;

 for (i=0;i<3;++i){
   xT[i] = 0.0;
   for (j=0;j<3;++j){xT[i] += R[j][i]*(xO[j]-tvec[j]);}
 }

} /* end of inv_transform_Rx_plus_tvec_3D() */



/**************/
/**** MAIN ****/
/**************/

int main(argc,argv)
  int argc;
  char **argv;
{
  int i,k;
  char OutType[8],imapfile[MAX_FILE_LEN],omapfile[MAX_FILE_LEN],download_file[MAX_FILE_LEN];
  char id[MAX_FILE_LEN],R_str[MAX_FILE_LEN],T_str[MAX_FILE_LEN], Raxis_str[MAX_FILE_LEN],AlgoType, FileType[16];
  char word[128],key[128],value[MAX_FILE_LEN],OutByteOrder[32];
  int Nword,Wsta[100],Wend[100];
  int Lquery_string;
  float Rmat[3][3],Tvec[3];
  float MIN[3],MAX[3]; 
  struct VOXEL Vori,Vnew;
  double maxMbyte;

  sprintf(OutType,"plain");
  sprintf(FileType,"emdb");
  id[0] = '\0';
  sprintf(R_str,"1:0:0:0:1:0:0:0:1");
  sprintf(Raxis_str,"0:0:0:0");
  sprintf(T_str,"0:0:0");
  imapfile[0] = '\0';
  AlgoType = 'L';
  omapfile[0] = '\0';
  maxMbyte = -1.0;
  sprintf(OutByteOrder,"default");

  if ((argc<2)&&(getenv("QUERY_STRING")==NULL)){
    printf("Content-type: text/html;charset=utf-8\n\n");
    printf("<HTML><BODY><PRE>\n");

    printf("transf_map.cgi [key1]=[value1] [key2]=[value2] ... \n");
    printf(" for transform (rotate and translate) 3D density map.\n");
    printf(" mainly for the 'pairwise_gmfit' WEB service.\n");
    printf(" This program is a part of the 'gmconvert' program.\n");
    printf("  coded by T.Kawabata. LastModDate:%s.\n",LastModDate);
    printf("<keys>\n");
    printf(" out  : type of output. 'plain', 'html' or 'download'[%s]\n",OutType);
    printf(" imap : Input map file [%s]\n",imapfile);
    printf(" omap : Output map file. If not assigned output in stdout. [%s]\n",omapfile);
    printf(" type : type for map. 'emdb','gmfit_up', 'omokage_up' [%s]\n",FileType);
    printf(" id   : ID for the map (EMDB_ID, UPLOAD_ID). [%s]\n",id);
    printf(" R    : rotation matrix for target[%s]\n",R_str);
    printf(" Raxis: rotation axis and angle [%s]\n",Raxis_str);
    printf(" T    : translational vector for target[%s]\n",T_str);
    printf(" algo : algorithm for density estimation. 'N'aive, 'A'verage, 'L'inear_interpolation [%c]\n",AlgoType);
    printf(" maxMbyte : maximum memory (Mbyte) for input image. negative values means no limit. [%lf]\n",maxMbyte);
    printf(" obyteorder : Output Byte Order 'default', 'HtoN','NtoH'  [%s]\n",OutByteOrder);
    printf("<note> 'type=direct' is only for debugging.\n");
    printf("</PRE></BODY></HTML>\n");
    exit(1);
  }
 
 /** [1] Read options **/ 
  if (getenv("QUERY_STRING")!=NULL){
    Lquery_string = strlen(getenv("QUERY_STRING"));
    ARG_STRING = (char *)malloc(sizeof(char)*(Lquery_string+1));
    sprintf(ARG_STRING,"%s",getenv("QUERY_STRING"));
    split_into_words(getenv("QUERY_STRING"),'&',&Nword,Wsta,Wend,100);
  }
  else{
    Nword = argc-1;
    Lquery_string = 0;
    for (k=0;k<argc;++k){ Lquery_string += (strlen(argv[k])+1); }
    ARG_STRING = (char *)malloc(sizeof(char)*(Lquery_string+1));
    ARG_STRING[0] = '\0';
    for (k=0;k<argc;++k){ strcat(ARG_STRING,argv[k]); strcat(ARG_STRING," ");}
  }

  for (k=0;k<Nword;++k){
    if (getenv("QUERY_STRING")==NULL){
      split_into_two_words(argv[k+1],0,'=',key,value);
    }
    else{
      get_part_of_line(word,getenv("QUERY_STRING"),Wsta[k],Wend[k]);
      split_into_two_words(word,0,'=',key,value);
    }
         if (strcmp(key,"out")==0) {sprintf(OutType,"%s",value); }
    else if (strcmp(key,"id")==0)  {sprintf(id,"%s",value); }
    else if (strcmp(key,"imap")==0){sprintf(imapfile,"%s",value); }
    else if (strcmp(key,"R")==0) {sprintf(R_str,"%s",value); }
    else if (strcmp(key,"Raxis")==0) {sprintf(Raxis_str,"%s",value); }
    else if (strcmp(key,"T")==0) {sprintf(T_str,"%s",value); }
    else if (strcmp(key,"algo")==0) {AlgoType = value[0];}
    else if (strcmp(key,"omap")==0) {sprintf(omapfile,"%s",value); }
    else if (strcmp(key,"type")==0) {sprintf(FileType,"%s",value); }
    else if (strcmp(key,"maxMbyte")==0) {maxMbyte = atof(value); }
    else if (strcmp(key,"obyteorder")==0) {sprintf(OutByteOrder,"%s",value); }
    else { printf("#ERROR:Can't understand the option '%s'.\n",key); exit(1);}
  }

  if (strcmp("out","html")==0){
    printf("Content-Type: text/html;charset=utf-8\n\n");
    printf("<HTML><BODY><PRE>\n");
  }
  else if (strcmp(OutType,"download")==0){
    /* printf("Content-type: text/download\n"); */
    printf("Content-Type: application/force-download\n"); 
    download_file[0] = '\0';
    if (id[0] != '\0'){
      sprintf(download_file,"%s.map",id);
    }
    else{
      sprintf(download_file,"transf_map.map");
    }
    printf("Content-Disposition: attachment; filename=%s\n\n",download_file);
  }
  else{
    printf("Content-Type: text/html;charset=utf-8\n\n");
  }

  get_Rmat_and_Tvec_from_strings(Rmat, Tvec, R_str, T_str);


  if (strcmp(Raxis_str,"0:0:0:0")!=0){
    get_Rmat_from_Raxis_str(Rmat,Raxis_str);
  }

  if (strcmp(FileType,"direct")!=0){
     imapfile[0]='\0';
  }

  if (imapfile[0]=='\0'){
    read_environment_file("/var/www/html/gmfit/cgi-bin/CONF/webgmfit.conf");
    if ((strcmp(FileType,"emdb")==0) &&(strlen(id)==4)){
      sprintf(imapfile,"%s/EMD-%s/map/emd_%s.map.gz",EMDB_DIR,id,id);
    }
    else if ((strcmp(FileType,"gmfit_up")==0) &&(strlen(id)==8)){
      sprintf(imapfile,"%s/%s.map",TMPOUT_DIR,id);
    }
    else if ((strcmp(FileType,"omokage_up")==0) &&(strlen(id)==33)){
      sprintf(imapfile,"%s/%s.map",OMOKAGE_USERDATA_DIR,id);
    }
 }

  Read_Density_File(imapfile,&Vori,'F',maxMbyte);
  
/*
  printf("#Rmat0 %f %f %f T %f\n",Rmat[0][0],Rmat[0][1],Rmat[0][2],Tvec[0]);
  printf("#Rmat1 %f %f %f T %f\n",Rmat[1][0],Rmat[1][1],Rmat[1][2],Tvec[1]);
  printf("#Rmat2 %f %f %f T %f\n",Rmat[2][0],Rmat[2][1],Rmat[2][2],Tvec[2]);
  */

  get_min_max_xyz_by_transformation(&Vori,Rmat,Tvec,MIN,MAX);
  if (omapfile[0] != '\0'){
    printf("#Orig %f %f %f\n",Vori.OrigPos[0],Vori.OrigPos[1],Vori.OrigPos[2]); 
    printf("#X MIN %f MAX %f\n",MIN[0],MAX[0]); 
    printf("#Y MIN %f MAX %f\n",MIN[1],MAX[1]); 
    printf("#Z MIN %f MAX %f\n",MIN[2],MAX[2]); 
  }

  Vnew.grid_width = Vori.grid_width;

  for (i=0;i<3;++i) {
    Vnew.OrigPos[i] = MIN[i];
    Vnew.N[i] = (int)ceil((MAX[i]-MIN[i])/Vnew.grid_width);
  }


  if (omapfile[0] != '\0'){
    printf("#N %d %d %d -> %d %d %d\n",Vori.N[0], Vori.N[1], Vori.N[2], Vnew.N[0], Vnew.N[1], Vnew.N[2]); 
  }

  Malloc_Voxel(&Vnew,Vnew.N[0],Vnew.N[1],Vnew.N[2]);

  setup_density_for_transformed_map(&Vnew,&Vori,Rmat,Tvec,AlgoType); 
 
 
  /* printf("#Vnew %d %d %d\n",Vnew.N[0],Vnew.N[1],Vnew.N[2]); */

  if (omapfile[0] == '\0'){
    if (strcmp(OutByteOrder,"HtoN")==0){
      Write_MapCCP4_File_HtoN(&Vnew,"-");
    }
    if (strcmp(OutByteOrder,"NtoH")==0){
      Write_MapCCP4_File_NtoH(&Vnew,"-");
    }
    else{
      Write_MapCCP4_File(&Vnew,"-"); 
    }
  }
  else{
    if (strcmp(OutByteOrder,"HtoN")==0){
      Write_MapCCP4_File_HtoN(&Vnew,omapfile);
    }
    if (strcmp(OutByteOrder,"NtoH")==0){
      Write_MapCCP4_File_NtoH(&Vnew,omapfile);
    }
    else{
       Write_MapCCP4_File(&Vnew,omapfile); 
    }
  }
  return(1);
}
