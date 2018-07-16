/*

<MCubeIO.c>
 
 functions for input/output of surface model genrated by
 Marching-Cube algorithm

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

/** Connection Matrix **/
struct CMATRIX{
 int N;
 char **m; };


/*** FUNCTIONS (GLOBAL) ***/
int Output_MC_Faces_VRML();
int Output_MC_Faces_Object();
void Read_MC_Faces_Object();
int  Output_MC_Edges_VRML();
int  Output_MC_Faces_PDB();
int  Output_MC_Edges_PDB();


/*** FUNCTIONS (GLOBAL) ***/
static char *Get_Date_String_PDB();
static void Output_Comment();


int Output_MC_Faces_VRML(fname,fouttype,Vhead,Fhead,vox,RGBT,comment)
  char *fname;
  char  fouttype;  /* 'a'ttend 'w'rite */
  struct MC_VERTEX *Vhead;
  struct MC_FACE   *Fhead;
  struct VOXEL  *vox;
  float  RGBT[4];  /* Red, Green, Blue, Transparency */
  char  *comment;
{
  FILE *fp;
  struct MC_VERTEX *v;
  struct MC_FACE   *f;
  int Nver;



  Nver = Number_Of_MC_VERTEX(Vhead);
 
  if (fname[0]=='-') fp = stdout;
  else
  { if (fouttype=='a') fp = fopen(fname,"a");
                  else fp = fopen(fname,"w");
    if (fp==NULL){printf("#ERROR:Can't write to vrmlfile \"%s\"\n",fname); return(0);}
    printf("#Output_MC_Faces_VRML(Nver %d) --> \"%s\"\n",Nver,fname); 

  }

 
  /** Output Header ***/
  fprintf(fp,"#VRML V2.0 utf8\n");
  fprintf(fp,"#%s\n",fname);
  Output_Comment(fp,comment,"#[COMMENT]",100);

  fprintf(fp,"Shape{\n");
  fprintf(fp," appearance Appearance{\n");
  fprintf(fp," material Material{\n");
  fprintf(fp,"  diffuseColor %.1f %.1f %.1f\n",RGBT[0],RGBT[1],RGBT[2]);
  fprintf(fp,"  transparency %.1lf\n",RGBT[3]);
  fprintf(fp," }\n");
  fprintf(fp,"}\n");

 /*** Output IndexedFaceSet ***/
 /* Output Points (in stacked order (first in, last out) */
  fprintf(fp,"geometry IndexedFaceSet{\n");
  fprintf(fp," coord Coordinate{\n");
  fprintf(fp," point[\n");
 
  v = Vhead;
  while (v->next != NULL){
    v = v->next;
    fprintf(fp,"%.3f %.3f %.3f",
         vox->OrigPos[0] + vox->grid_width*v->pos[0],
         vox->OrigPos[1] + vox->grid_width*v->pos[1],
         vox->OrigPos[2] + vox->grid_width*v->pos[2]);
   if (v->next != NULL) fprintf(fp,","); 
   fprintf(fp,"\n"); 
  }
  fprintf(fp,"]}\n");
 
  /* Output MC_FACES  (Vertex nums are modified in the stacked order) */
  /*  [vertex_num] = [Nver] - [a->num] -1 */

  fprintf(fp,"  coordIndex[\n");
  f = Fhead;
  while (f->next != NULL){
    f = f->next;
    if ((f->a != NULL)&& (f->b != NULL) && (f->c != NULL)){
      fprintf(fp,"%d, %d, %d, -1", Nver - f->a->num -1, Nver - f->b->num -1, Nver - f->c->num -1); 
      if (f->next !=NULL) fprintf(fp,",");
      fprintf(fp,"\n");
    }
  }

  fprintf(fp,"]\n");
  /* fprintf(fp,"solid FALSE\n"); */
  fprintf(fp,"}\n");
  fprintf(fp,"}\n");
  if (fp != stdout) fclose(fp);
  return(1);
} /* end of Output_MC_Faces_VRML() */





int Output_MC_Faces_Object(fname,fouttype,Vhead,Fhead,vox,comment)
  char *fname;
  char  fouttype;  /* 'a'ttend 'w'rite */
  struct MC_VERTEX *Vhead;
  struct MC_FACE   *Fhead;
  struct VOXEL  *vox;
  char  *comment;
{
  FILE *fp;
  struct MC_VERTEX *v;
  struct MC_FACE   *f;
  int Nver;

  Nver = Number_Of_MC_VERTEX(Vhead);
 
  if (fname[0]=='-') fp = stdout;
  else
  { if (fouttype=='a') fp = fopen(fname,"a");
                  else fp = fopen(fname,"w");
    if (fp==NULL){printf("#ERROR:Can't write to object file \"%s\"\n",fname); return(0);}
    printf("#Output_MC_Faces_Object(Nver %d) --> \"%s\"\n",Nver,fname); 
  }

 
  /** Output Vertexes ***/
  v = Vhead;
  while (v->next != NULL){
    v = v->next;
    fprintf(fp,"v %.3f %.3f %.3f\n",
         vox->OrigPos[0] + vox->grid_width*v->pos[0],
         vox->OrigPos[1] + vox->grid_width*v->pos[1],
         vox->OrigPos[2] + vox->grid_width*v->pos[2]);
  }
 
  /** Output Faces ***/
  f = Fhead;
  while (f->next != NULL){
    f = f->next;
    fprintf(fp,"f %d %d %d\n", Nver - f->a->num, Nver - f->b->num, Nver - f->c->num);
  }
  return(1);
} /* end of Output_MC_Faces_Object() */





void Read_MC_Faces_Object(ifname,Vhead,Fhead)
  char *ifname;
  struct MC_VERTEX *Vhead;
  struct MC_FACE   *Fhead;
{
  FILE *fp;
  char line[100],word[100];
  int i,L,Nword,Wsta[32],Wend[32],a,b,c,Nvertex,nvertex;
  float pos[3]; 
  struct MC_VERTEX **ver_pointers; /* (malloc later) */
/* 
>> "*.obj" file format example <<
v 133.093 76.093 57.000
v 133.109 76.000 57.109
:
v 37.989 76.000 57.000
v 38.000 75.988 57.000
f 2 1 15
f 1 19 15
f 2 15 23
:
*/
  Vhead->next = NULL;
  Fhead->next = NULL;
  Vhead->num = -1;
  Fhead->num = -1;

  fp = fopen(ifname,"r");
  if (fp==NULL){
    printf("#ERROR:Can't open object file \"%s\"\n",ifname);
    exit(1);
  }

  /* just conunt Nvertex, and malloc ver_pointers[Nvertex] */
  Nvertex = 0;
  while (feof(fp)==0){
    line[0] = '\0';
    fgets(line,100,fp);
    if (line[0]=='v') { Nvertex += 1;} 
  }
  fclose(fp);

  printf("#Read_MC_Faces_Object('%s' Nvertex %d)\n",ifname,Nvertex);
  ver_pointers = (struct MC_VERTEX**)malloc(sizeof(struct MC_VERTEX*)*Nvertex);
  for (i=0;i<Nvertex;++i){ver_pointers[i] = NULL;}

  /* Read again */
  fp = fopen(ifname,"r");
  nvertex = 0;
  while (feof(fp)==0){
    line[0] = '\0';
    fgets(line,100,fp);
    L = strlen(line);
    if (line[L-1]=='\n'){ line[L-1] = '\0';}

    if ((line[0]!='#') && (L>5)){
      if (line[0]=='v'){
        Split_to_Words(line,' ',&Nword,Wsta,Wend,32);
        if (Nword>=4){
          Get_Part_Of_Line(word,line,Wsta[1],Wend[1]); pos[0] = atof(word);
          Get_Part_Of_Line(word,line,Wsta[2],Wend[2]); pos[1] = atof(word);
          Get_Part_Of_Line(word,line,Wsta[3],Wend[3]); pos[2] = atof(word);
          ver_pointers[nvertex] = Add_MC_Vertex_Stack(Vhead,pos);
          nvertex += 1;
        }
      }
      if (line[0]=='f'){
        Split_to_Words(line,' ',&Nword,Wsta,Wend,32);
        if (Nword>=4){
          Get_Part_Of_Line(word,line,Wsta[1],Wend[1]); a = atoi(word);
          Get_Part_Of_Line(word,line,Wsta[2],Wend[2]); b = atoi(word);
          Get_Part_Of_Line(word,line,Wsta[3],Wend[3]); c = atoi(word);
          /* printf("a %d b %d c %d\n",a,b,c); */
          if ((1<=a)&&(a<=Nvertex)&&(1<=b)&&(b<=Nvertex)&&(1<=c)&&(c<=Nvertex)){
            /* printf("a %d (%f %f %f)\n",a, ver_pointers[a-1]->pos[0], ver_pointers[a-1]->pos[1], ver_pointers[a-1]->pos[2]); */
            Add_MC_Face_Stack_Simple(Fhead, ver_pointers[a-1], ver_pointers[b-1], ver_pointers[c-1]);
          }
        }
      }
    }
  } 

  fclose(fp);
  free(ver_pointers);
} /* end of Read_MC_Faces_Object() */








int Output_MC_Edges_VRML(fname,fouttype,Vhead,Ehead,vox,RGBT,comment)
  char *fname;
  char  fouttype;  /* 'a'ttend 'w'rite */
  struct MC_VERTEX *Vhead;
  struct MC_EDGE   *Ehead;
  struct VOXEL  *vox;
  float  RGBT[4];  /* Red, Green, Blue, Transparency */
  char  *comment;
{
  FILE *fp;
  struct MC_VERTEX *v;
  struct MC_EDGE *e;
  int Nver;
  
  Nver = Number_Of_MC_VERTEX(Vhead);
 
  if (fname[0]=='-') fp = stdout;
  else{ 
    if (fouttype=='a') fp = fopen(fname,"a");
                  else fp = fopen(fname,"w");
    if (fp==NULL){printf("#ERROR:Can't write to vrmlfile \"%s\"\n",fname); return(0);}
    printf("#Output_MC_Edges_VRML(Nver %d) --> \"%s\"\n",Nver,fname); 
  }
 
  /** Output Header ***/
  fprintf(fp,"#VRML V2.0 utf8\n");
  fprintf(fp,"#%s\n",fname);
  Output_Comment(fp,comment,"#[COMMENT]",100);
  fprintf(fp,"#Nver %d\n\n",Nver);

  fprintf(fp,"Shape{\n");

 /*** Output IndexedLineSet ***/
 /* Output Points (in stacked order (first in, last out) */
  fprintf(fp,"geometry IndexedLineSet{\n");
  fprintf(fp," coord Coordinate{\n");
  fprintf(fp," point[\n");
 
  v = Vhead;
  while (v->next != NULL){
    v = v->next;
    fprintf(fp,"%.3f %.3f %.3f",
         vox->OrigPos[0] + vox->grid_width*v->pos[0],
         vox->OrigPos[1] + vox->grid_width*v->pos[1],
         vox->OrigPos[2] + vox->grid_width*v->pos[2]);
   if (v->next != NULL) fprintf(fp,",");
   fprintf(fp,"\n");
  }
  fprintf(fp,"]}\n");
 
  /* Output Edges  */
  fprintf(fp,"  coordIndex[\n");
  e = Ehead;
  while (e->next != NULL){
    e = e->next;
    fprintf(fp,"%d,%d,-1", Nver - e->a->num -1, Nver - e->b->num -1);
    if (e->next !=NULL) fprintf(fp,",");
    fprintf(fp,"\n");
  }

  fprintf(fp,"]\n");
  /* fprintf(fp,"solid FALSE\n"); */
  fprintf(fp,"colorPerVertex FALSE\n"); 
  fprintf(fp,"color Color{\n"); 
  fprintf(fp,"  color[%f %f %f]\n",RGBT[0],RGBT[1],RGBT[2]); 
  fprintf(fp,"}\n"); 
  fprintf(fp,"colorIndex[\n"); 
  e = Ehead;
  while (e->next != NULL){
    e = e->next;
    fprintf(fp,"0\n");
  }


  fprintf(fp,"]\n"); 
  
  fprintf(fp,"}\n");
  fprintf(fp,"}\n");
  if (fp != stdout) fclose(fp);
  return(1);

} /* end of Output_MC_Edges_VRML() */





















int Output_MC_Faces_PDB(fname,fouttype,Vhead,Fhead,vox,tFactor,comment,Chain)
 char *fname;
 char  fouttype;  /* 'a'ttend 'w'rite */
 struct MC_VERTEX *Vhead;
 struct MC_FACE   *Fhead;
 struct VOXEL  *vox;
 float  tFactor;
 char *comment;
 char Chain;
{
 FILE *fp;
 struct MC_VERTEX *v;
 struct MC_FACE   *f;
 struct CMATRIX Cmat;  /* Connection Matrix */
 int Nver,i,j,n,a,b,c,anum;
 int *Carray;
 char OutFace;
 
 /** [1] Open file **/
 if (fname[0]=='-') fp = stdout;
 else
  { if (fouttype=='a') fp = fopen(fname,"a");
                 else fp = fopen(fname,"w");
   if (fp==NULL){printf("#ERROR:Can't write to surf_pdb file\"%s\"\n",fname); return(0);}
   printf("#Output_MC_Faces_PDB() --> \"%s\"\n",fname); }

 /** [2] Output HEADER and REMARK **/ 
 fprintf(fp,"HEADER    Surface File by MarCube                 %s   %4s      %4s\n",
            Get_Date_String_PDB(),"0XXX","0XXX");
 Output_Comment(fp,comment,"REMARK  ",80);

 /** [3] Count Nvertex **/ 
 Nver = Number_Of_MC_VERTEX(Vhead);
 printf("#Nvertex %d\n",Nver);
 fprintf(fp,"REMARK  Nvertex %d\n",Nver);
 if (Nver >99999) 
 { printf("#WARNING(Output_MC_Faces_PDB)[%s]:Natom %d is too large for CONECT.\n",fname,Nver);
   fprintf(fp,"REMARK   #WARNING:Natom %d is too large for CONECT.\n",Nver);
   OutFace = 0;}
 else OutFace = 0;
 /* [4] Output Points (in the stacked order) as "HETATM" */

 v = Vhead; 
 while (v->next != NULL){
  v = v->next;
  anum = Nver - v->num-1;
  if (anum>99999) anum = anum % 100000;

  fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n",
       anum+1 ,"MAP","MAP",Chain,(anum/10)+1,
       vox->OrigPos[0] + vox->grid_width*v->pos[0],
       vox->OrigPos[1] + vox->grid_width*v->pos[1],
       vox->OrigPos[2] + vox->grid_width*v->pos[2],tFactor,tFactor);
 }

 if (OutFace == 1){  
 /* [5] Output Faces as "CONECT" */

 /** (1) Malloc Camt and Carray **/

 Malloc_CMATRIX(&Cmat,Nver);
 Carray = (int *)malloc(sizeof(int)*Nver);
 for (i=0;i<Nver;++i)
  { Carray[i] = 0;
    for (j=0;j<Nver;++j) Cmat.m[i][j] = 0; }

 /** (2) Make Camt and Carray **/
 f = Fhead;
 while (f->next != NULL){
  f = f->next;
  a = Nver - f->a->num-1; 
  b = Nver - f->a->num-1; 
  c = Nver - f->a->num-1; 
  Cmat.m[a][b] = Cmat.m[b][a] = 1;
  Cmat.m[b][c] = Cmat.m[c][b] = 1;
  ++Carray[a]; ++Carray[b]; ++Carray[c];
 }

 /** (3) Output CONECT **/
 fprintf(fp,"TER\n");
 for (i=0;i<Nver;++i){
  if (Carray[i]>0){ 
   fprintf(fp,"CONECT%5d",i+1);
   n = 0;
   for (j=0;j<Nver;++j) 
   { if (Cmat.m[i][j]==1) 
     { fprintf(fp,"%5d",j+1); ++n;
       if ((Carray[i]>n)&&((n%4)==0)) {fprintf(fp,"\nCONECT%5d",i+1);} 
      } 
     }    
    fprintf(fp,"\n"); 
   }
  }
 
 Free_CMATRIX(&Cmat);
 free(Carray);
 } 
fprintf(fp,"END\n"); 
 if (fp != stdout) fclose(fp);
 return(1);

} /* end of Output_MC_Faces_PDB() */




int Output_MC_Edges_PDB(fname,fouttype,Vhead,Ehead,vox,tFactor,comment,Chain,offset_atomnum)
 char *fname;
 char  fouttype;  /* 'a'ttend 'w'rite */
 struct MC_VERTEX *Vhead;
 struct MC_EDGE   *Ehead;
 struct VOXEL  *vox;
 float  tFactor;
 char *comment;
 char Chain;
 int  offset_atomnum;
{
 FILE *fp;
 struct MC_VERTEX *v;
 struct MC_EDGE   *e;
 int Nvertex,anum;
 char OutFace;
 int a,b,**ConnectAtomNums,*Nconnect;
 
 /** [1] Open file **/
 if (fname[0]=='-') fp = stdout;
 else
  { if (fouttype=='a') fp = fopen(fname,"a");
                 else fp = fopen(fname,"w");
   if (fp==NULL){printf("#ERROR:Can't write to surf_pdb file\"%s\"\n",fname); exit(1);}
   printf("#Output_MC_Edges_PDB() --> \"%s\"\n",fname); }

 /** [2] Output HEADER and REMARK **/ 
 fprintf(fp,"HEADER    Surface File by MarCube                 %s   %4s      %4s\n",
            Get_Date_String_PDB(),"0XXX","0XXX");
 fprintf(fp,"REMARK    COMMAND '%s'\n",PAR.COMMAND);
 fprintf(fp,"REMARK    grid_width %f tFactor %f offset_atomnum %d\n",vox->grid_width,tFactor,offset_atomnum);
 Output_Comment(fp,comment,"REMARK  ",80);

 /** [3] Count Nvertex **/ 
 Nvertex = Number_Of_MC_VERTEX(Vhead);
 printf("#Nvertex %d\n",Nvertex);
 fprintf(fp,"REMARK  Nvertex %d\n",Nvertex);
 if (Nvertex >99999){ 
   printf("#WARNING(Output_MC_Faces_PDB)[%s]:Natom %d is too large for CONECT.\n",fname,Nvertex);
   fprintf(fp,"REMARK   #WARNING:Natom %d is too large for CONECT.\n",Nvertex);
   OutFace = 0;
 }
 else OutFace = 1;

 printf("#OutFace %d\n",OutFace);
 /* [4] Output Points (in the stacked order) as "HETATM" */
 /* fprintf(fp,"MODEL        1\n"); */
 v = Vhead; 
 while (v->next != NULL){
  v = v->next;
  anum = Nvertex - v->num-1;
  if (anum>99999) anum = anum % 100000;

  fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n",
       anum+1 + offset_atomnum,"MAP","MAP",Chain,(anum/10)+1,
       vox->OrigPos[0] + vox->grid_width*v->pos[0],
       vox->OrigPos[1] + vox->grid_width*v->pos[1],
       vox->OrigPos[2] + vox->grid_width*v->pos[2],tFactor,tFactor);
 }

 /** (3) Output CONECT **/
 if (OutFace == 1){  
   ConnectAtomNums = (int **)malloc(sizeof(int*)*Nvertex);
   Nconnect = (int *)malloc(sizeof(int)*Nvertex);
   for (a=0;a<Nvertex;++a){
     Nconnect[a] = 0;
     ConnectAtomNums[a] = (int *)malloc(sizeof(int)*6);
     for (b=0;b<6;++b){
       ConnectAtomNums[a][b] = -1;
     } 
   }
   e = Ehead; 
   while (e->next != NULL){
     e = e->next;
     ConnectAtomNums[e->a->num][Nconnect[e->a->num]] = e->b->num; 
     Nconnect[e->a->num] += 1;
     ConnectAtomNums[e->b->num][Nconnect[e->b->num]] = e->a->num; 
     Nconnect[e->b->num] += 1;
   }
   for (a=0;a<Nvertex;++a){
     if (Nconnect[a]>0){
       fprintf(fp,"CONECT%5d",Nvertex - a + offset_atomnum);
       for (b=0;b<6;++b){
         if (ConnectAtomNums[a][b]>=0){
           fprintf(fp,"%5d",Nvertex - ConnectAtomNums[a][b] + offset_atomnum);
         }
       }
      fprintf(fp,"\n");
     }
   }
   /* fprintf(fp,"ENDMDL\n"); */
   fprintf(fp,"TER\n"); 
 }
 if (fp != stdout) fclose(fp);
 return(Nvertex);

} /* end of Output_MC_Edges_PDB() */









void Output_Comment(fp,comment,header,Nsymb_oneline)
 FILE *fp;
 char *comment,*header;
 int  Nsymb_oneline;  /* Number of symbol for one line */
{
 int Lcomment,i,j;
                                                                                
 Lcomment = strlen(comment);
 if (comment[Lcomment-1]=='\n') Lcomment -= 1;

 fprintf(fp,"%s ",header);
 i = 0; j = 0;
 while (i<Lcomment){
   if ((comment[i]=='\n')||(j==Nsymb_oneline))
    { fprintf(fp,"\n");
      if (i<Lcomment) fprintf(fp,"%s ",header);
      j = 0;}
   if (comment[i]!='\n') { fprintf(fp,"%c",comment[i]);  ++j;}
   ++i;
 }
 if (j>=0) fprintf(fp,"\n");
                                                                                
} /* end of Output_Comment() */


char *Get_Date_String_PDB()
{
 time_t      now_t;
 struct tm  *loc_t;
 static char Mon[][4] = {"JAN","FEB","MAR","APR","MAY","JUN",
                         "JUL","AUG","SEP","OCT","NOV","DEC"};
 static char string[64];
 char day[16],year[16];
 int  Nyear;
 now_t = time(NULL);
 loc_t = localtime(&now_t);
                                                                                
 if (loc_t->tm_mday <10) sprintf(day,"0%d",loc_t->tm_mday);
 else sprintf(day,"%2d",loc_t->tm_mday);
                                                                                
 if (loc_t->tm_year<100) Nyear = loc_t->tm_year;
   else                  Nyear = loc_t->tm_year - 100;
                                                                                
 if (Nyear < 10) sprintf(year,"0%d",Nyear);
            else sprintf(year,"%2d",Nyear);
 sprintf(string,"%s-%3s-%s",day, Mon[loc_t->tm_mon],year);
 return(string);
                                                                                
} /* end of Get_Date_String_PDB() */
                                                                                



