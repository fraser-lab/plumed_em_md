/*
 <Radius.c>
 
 Functions for assigning radius to atoms

*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "globalvar.h" 
#include "pdbstruct.h" 
#include "radius_vdw.h" 
#include "radius_cov.h" 
#include "BasicIO.h" 

/*** FUNCTIONS (GLOBAL) ***/
void Read_Radius_File();
void Set_Default_Radius();
void Assign_Radius();
void Assign_Uniform_Radius();
void Assign_Weight();
void Assign_Uniform_Weight();
void Assign_Radius_by_Residue();
void Assign_Weight_by_Residue();
float Atomic_Weight();
void  write_radius_weight_of_atoms();
void Set_Variance_of_Atoms();


/*** FUNCTIONS (LOCAL) ***/
static int Match();
static float Atomic_Weight_by_Residue();
static float Radius_by_Residue();


void Read_Radius_File(fname,Rhead)
 char *fname;
 struct RADIUS *Rhead;
{
 FILE *fp;
 char line[128],atom[5],resi[4],buff[32];
 int L;
 struct RADIUS *rn;

 fp = fopen(fname,"r");
 if (fp==NULL){printf("#ERROR:Can't open radius file \"%s\"\n",fname); exit(1);}

 rn = Rhead;
 rn->next = rn->prev = NULL;
 
 while (feof(fp)==0){ 
  line[0] = '\0';
  fgets(line,127,fp); 
  L = strlen(line);
  if ((line[0] != '#')&&(L>5)){
    Get_Part_Of_Line(atom,line,0,3);
    Get_Part_Of_Line(resi,line,5,7);
    Get_Part_Of_Line(buff,line,8,L-1);
    rn->next = (struct RADIUS *)malloc(sizeof(struct RADIUS));
    rn->next->prev = rn;
    rn->next->next = NULL;
    rn  = rn->next; 
    sprintf(rn->Atom,"%s",atom); 
    sprintf(rn->Resi,"%s",resi); 
    rn->R = atof(buff); 
  }  
 } 
 fclose(fp);

} /* end of Read_Radius_File() */







void Set_Default_Radius(Rhead,RadiusType)
 struct RADIUS *Rhead;
 char   RadiusType; /* 'V'dw:Rvdw, 'C':Rcovalent */
{
 char line[128],atom[5],resi[4],buff[32];
 int i,L,N;
 struct RADIUS *rn;

 rn = Rhead;
 rn->next = rn->prev = NULL;

   if (RadiusType=='C'){ N = N_RCOV_DATA_LINE;}
 else { N = N_RVDW_DATA_LINE;}

 for (i=0;i<N;++i){ 
   if (RadiusType=='C'){sprintf(line,"%s",RCOV_DATA_LINE[i]);}
                  else {sprintf(line,"%s",RVDW_DATA_LINE[i]);}
   L = strlen(line);
   if ((line[0] != '#')&&(L>5)){
     Get_Part_Of_Line(atom,line,0,3);
     Get_Part_Of_Line(resi,line,5,7);
     Get_Part_Of_Line(buff,line,8,L-1);
     rn->next = (struct RADIUS *)malloc(sizeof(struct RADIUS));
     rn->next->prev = rn;
     rn->next->next = NULL;
     rn  = rn->next; 
     sprintf(rn->Atom,"%s",atom); 
     sprintf(rn->Resi,"%s",resi); 
     rn->R = atof(buff); 
   }  
 } 

} /* end of Set_Default_Radius() */






void Assign_Radius(Ahead,Rhead)
 struct ATOM   *Ahead;
 struct RADIUS *Rhead;
{
 struct ATOM   *an;
 struct RADIUS *rn;
 char hit;

 an = Ahead;
 while (an->next != NULL){
   an = an->next;
   rn = Rhead;
   hit = 0;
   while ((hit==0)&&(rn->next != NULL)){
     rn = rn->next;
     if ((Match(rn->Resi,3,an->Resi,0)==1) &&
         (Match(rn->Atom,3,an->Atom,0)==1) ){
       hit = 1; 
       an->R = rn->R; 
     }
   }
   /* printf("#%s %s R %f\n",an->Resi,an->Atom,an->R); */
 } 

} /* end of Assign_Radius() */




void Assign_Uniform_Radius(Ahead,Runiform)
 struct ATOM   *Ahead;
 float  Runiform;
{
 struct ATOM   *an;

 printf("#Assign_Uniform_Radius(Runiform:%f)\n",Runiform);
 an = Ahead;
 while (an->next != NULL){
   an = an->next;
   an->R = Runiform;
 } 

} /* end of Assign_Uniform_Radius() */




int Match(Pat,plen,Str,offset)
 char *Pat;
 int  plen;
 char *Str;
 int offset;
{
  /*** Strict Pattern Matching ***/

  int i,j,ok,slen;

  slen = strlen(Str);
  if ((plen==0)||(slen==0)) return(0);

  if (slen>=plen){
   for (i= offset;i<(slen - plen+1);++i){ 
     ok = 1; j = 0;
     while ((ok==1) && (j<plen)){
       if ((Pat[j]!='x')&&(Pat[j] != Str[i+j])) { ok = 0; }
       ++j; 
     }
     if (ok==1) return(1);
   }
  }
 return(0);

} /* end of Match */


void Assign_Weight(Ahead)
 struct ATOM   *Ahead;
{
 struct ATOM   *an;
 char ele[4];

 an = Ahead;
 while (an->next != NULL){
   an = an->next;
   ele[0] = an->Atom[1];   
   ele[1] = '\0';
   an->Weight = Atomic_Weight(ele); 
   /* printf("#'%s' '%s' ele '%s' W %f\n",an->Resi,an->Atom,ele,an->Weight); */
 } 
} /* end of Assign_Weight() */



void Assign_Radius_by_Residue(Ahead)
 struct ATOM   *Ahead;
{
 struct ATOM   *an;

 an = Ahead;
 while (an->next != NULL){
   an = an->next;
   an->R = Radius_by_Residue(an->Resi); 
   /* printf("#'%s' '%s' '%s' R %f\n",an->Resi,an->Atom,an->Resi,an->R); */
 } 
} /* end of Assign_Weight_by_Residue() */




void Assign_Weight_by_Residue(Ahead)
 struct ATOM   *Ahead;
{
 struct ATOM   *an;

 an = Ahead;
 while (an->next != NULL){
   an = an->next;
   an->Weight = Atomic_Weight_by_Residue(an->Resi); 
   /* printf("#'%s' '%s' '%s' W %f\n",an->Resi,an->Atom,an->Resi,an->Weight); */
 } 
} /* end of Assign_Weight_by_Residue() */


void Assign_Uniform_Weight(Ahead)
 struct ATOM   *Ahead;
{
 struct ATOM   *an;

 printf("#Assign_Uniform_Weight()\n");
 an = Ahead;
 while (an->next != NULL){
   an = an->next;
   an->Weight = 1.0;
 } 
} /* end of Assign_Weight() */


float Atomic_Weight(ele)
  char *ele; /* element name */
{
 /*
  These values are taken from 
  Urabyoushi I. Genshiryou-hyou (1985)
  in "Rikagakujiten dai-4-ban".
 */
       if (strncmp(ele,"C",1)==0){ return(12.011);}
  else if (strncmp(ele,"N",1)==0){ return(14.00674);}
  else if (strncmp(ele,"O",1)==0){ return(15.9994);}
  else if (strncmp(ele,"S",1)==0){ return(32.066);}
  else if (strncmp(ele,"P",1)==0){ return(30.973762);}
  else { return(12.011);}
}


float Atomic_Weight_by_Residue(res)
  char *res; /* three-letter residue name */
{
 /*
 These weights were calculated by "src/res_rg_python/cal_mrg_pdb.py",
 using ideal 3D coordinates download from the LIGAND_EXPO web server. 
 Hydrogen atoms were NOT considered for the calculation.
 */
       if (strncmp(res,"ALA",3)==0){return(82.03854);}
  else if (strncmp(res,"ARG",3)==0){return(160.09176);}
  else if (strncmp(res,"ASN",3)==0){return(148.07768);}
  else if (strncmp(res,"ASP",3)==0){return(126.04834);}
  else if (strncmp(res,"CYS",3)==0){return(114.10454);}
  else if (strncmp(res,"GLN",3)==0){return(160.08868);}
  else if (strncmp(res,"GLU",3)==0){return(138.05934);}
  else if (strncmp(res,"GLY",3)==0){return(70.02754);}
  else if (strncmp(res,"HIS",3)==0){return(146.08502);}
  else if (strncmp(res,"ILE",3)==0){return(214.15954);}
  else if (strncmp(res,"LEU",3)==0){return(190.13754);}
  else if (strncmp(res,"LYS",3)==0){return(132.07828);}
  else if (strncmp(res,"MET",3)==0){return(138.12654);}
  else if (strncmp(res,"PHE",3)==0){return(154.10454);}
  else if (strncmp(res,"PRO",3)==0){return(106.06054);}
  else if (strncmp(res,"SER",3)==0){return(98.03794);}
  else if (strncmp(res,"THR",3)==0){return(146.08194);}
  else if (strncmp(res,"TRP",3)==0){return(192.13328);}
  else if (strncmp(res,"TYR",3)==0){return(170.10394);}
  else if (strncmp(res,"VAL",3)==0){return(178.12654);}
  else if (strncmp(res,"  A",3)==0){return(409.12186);}
  else if (strncmp(res,"  G",3)==0){return(425.12126);}
  else if (strncmp(res,"  C",3)==0){return(385.09678);}
  else if (strncmp(res,"  T",3)==0){return(379.11264);}
  else if (strncmp(res,"  U",3)==0){return(387.08944);}
  else                            {return(82.03854);}
} /* end of Atomic_Weight_by_Residue() */



float Radius_by_Residue(res)
  char *res; /* three-letter residue name */
{
 /*
 Raidius of gyration were calculated by "src/mrg_python/cal_mrg_pdb.py",
 using ideal 3D coordinates download from the LIGAND_EXPO web server. 
 Hydrogen atoms were not considered for the calculation.
 These radius are radius of gyration of GMM, assigning gdfs 
 with variance =Rvdw*Rvdw/5.0.

 For the sphere, Rg = 3/5 * Rvdw^2  => Rvdw^2 = 5/3 * Rg
             ==> Rvdw = sqrt(5/3) * Rg
 */
  double Rg,Rvdw;

/*
       if (strncmp(res,"ALA",3)==0){Rg = 2.187866;}
  else if (strncmp(res,"ARG",3)==0){Rg = 3.528325;}
  else if (strncmp(res,"ASN",3)==0){Rg = 2.840243;}
  else if (strncmp(res,"ASP",3)==0){Rg = 2.701685;}
  else if (strncmp(res,"CYS",3)==0){Rg = 2.622403;}
  else if (strncmp(res,"GLN",3)==0){Rg = 3.223137;}
  else if (strncmp(res,"GLU",3)==0){Rg = 3.046054;}
  else if (strncmp(res,"GLY",3)==0){Rg = 2.059648;}
  else if (strncmp(res,"HIS",3)==0){Rg = 2.841412;}
  else if (strncmp(res,"ILE",3)==0){Rg = 2.688200;}
  else if (strncmp(res,"LEU",3)==0){Rg = 2.856076;}
  else if (strncmp(res,"LYS",3)==0){Rg = 3.339634;}
  else if (strncmp(res,"MET",3)==0){Rg = 3.141510;}
  else if (strncmp(res,"PHE",3)==0){Rg = 3.154383;}
  else if (strncmp(res,"PRO",3)==0){Rg = 2.484154;}
  else if (strncmp(res,"SER",3)==0){Rg = 2.382992;}
  else if (strncmp(res,"THR",3)==0){Rg = 2.514618;}
  else if (strncmp(res,"TRP",3)==0){Rg = 3.413389;}
  else if (strncmp(res,"TYR",3)==0){Rg = 3.376502;}
  else if (strncmp(res,"VAL",3)==0){Rg = 2.479064;}
  else if (strncmp(res,"  A",3)==0){Rg = 4.388194;}
  else if (strncmp(res,"  G",3)==0){Rg = 4.497954;}
  else if (strncmp(res,"  C",3)==0){Rg = 4.007091;}
  else if (strncmp(res,"  T",3)==0){Rg = 3.757795;}
  else if (strncmp(res,"  U",3)==0){Rg = 4.016172;}
  else                             {Rg = 2.187866;}
 */

       if (strncmp(res,"ALA",3)==0){ Rg = 1.963913;}
  else if (strncmp(res,"ARG",3)==0){ Rg = 3.374007;}
  else if (strncmp(res,"ASN",3)==0){ Rg = 2.695111;}
  else if (strncmp(res,"ASP",3)==0){ Rg = 2.525241;}
  else if (strncmp(res,"CYS",3)==0){ Rg = 2.413249;}
  else if (strncmp(res,"GLN",3)==0){ Rg = 3.088783;}
  else if (strncmp(res,"GLU",3)==0){ Rg = 2.883527;}
  else if (strncmp(res,"GLY",3)==0){ Rg = 1.841949;}
  else if (strncmp(res,"HIS",3)==0){ Rg = 2.652737;}
  else if (strncmp(res,"ILE",3)==0){ Rg = 2.575828;}
  else if (strncmp(res,"LEU",3)==0){ Rg = 2.736953;}
  else if (strncmp(res,"LYS",3)==0){ Rg = 3.177825;}
  else if (strncmp(res,"MET",3)==0){ Rg = 2.959014;}
  else if (strncmp(res,"PHE",3)==0){ Rg = 2.979213;}
  else if (strncmp(res,"PRO",3)==0){ Rg = 2.266054;}
  else if (strncmp(res,"SER",3)==0){ Rg = 2.184637;}
  else if (strncmp(res,"THR",3)==0){ Rg = 2.366486;}
  else if (strncmp(res,"TRP",3)==0){ Rg = 3.248871;}
  else if (strncmp(res,"TYR",3)==0){ Rg = 3.217711;}
  else if (strncmp(res,"VAL",3)==0){ Rg = 2.351359;}
  else if (strncmp(res,"  A",3)==0){ Rg = 4.333750;}
  else if (strncmp(res,"  T",3)==0){ Rg = 3.700942;}
  else if (strncmp(res,"  G",3)==0){ Rg = 4.443546;}
  else if (strncmp(res,"  C",3)==0){ Rg = 3.954067;}
  else if (strncmp(res,"  U",3)==0){ Rg = 3.964129;}
  else                             { Rg = 1.963913;}

  Rvdw = sqrt(5.0/3.0) * Rg;
  return(Rvdw);

} /* end of Radius_by_Residue() */


void write_radius_weight_of_atoms(ofname,Ahead)
  char *ofname;
  struct ATOM *Ahead;
{
  FILE *fpo;
  struct ATOM *an;

  printf("#write_radius_weight_of_atoms()-->'%s'\n",ofname);
  fpo = fopen(ofname,"w");
  an = Ahead;
  fprintf(fpo,"#>> Radius and Weight Values for Atomic Model <<\n");
  fprintf(fpo,"#COMMAND %s\n",PAR.COMMAND);
  fprintf(fpo,"#DATE    %s\n",PAR.START_DATE);
  fprintf(fpo,"#[Anum] [Atom] [Resi] [Chain] [R] [W]\n");
  while (an->next != NULL){
    an = an->next;
    fprintf(fpo,"%s %s %s %c   %f   %f\n",an->Anum,an->Atom,an->Resi,an->Chain,an->R,an->Weight);
  }  
  fclose(fpo);
} /* end of write_radius_weight_of_atoms() */



void Set_Variance_of_Atoms(Ahead, VarType, Crr2var, Resolution)
  struct ATOM *Ahead;
  char   VarType;      /* 'A':Var = Crr2var * Rvdw*Rvdw, 'R':Var = [Reso/2]^2 */
  double Crr2var;      /* Constant bwn variance and Rvdw. Variance = Crr2var * Rvdw * Rvdw */
  double Resolution;   /* Resolution */
{
  struct ATOM *an;

  printf("#Set_Variance_of_Atoms(VarType:%c Crr2var:%lf Resolution:%lf)\n",VarType,Crr2var,Resolution);
  an = Ahead;

  while (an->next != NULL){ 
    an = an->next;
    if (VarType=='R'){
      an->variance = (Resolution*Resolution)/4.0;
    }
    else{
     an->variance = Crr2var * an->R * an->R;
   }
  }
 
} /* end of Set_Variance_of_Atoms() */
 


