/*

<PdbIO.c>

 functions for input/output of PDB files
 using list structure (ATOM and RESIDUE node). 

*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include "globalvar.h" 
#include "pdbstruct.h" 
#include "gauss.h" 
#include "BasicIO.h" 
#include "Matrix3D.h" 
#include "GaussIO.h" 

/*** Functions (GLOBAL) ***/
void Read_PDB_File();
void Make_Residue_List();
char *Get_Date_String_PDB();
int  Number_Of_Atom();
int  Number_Of_Residue();
int  Number_Of_Chain();
char Chain_Of_Atoms();
void Malloc_MATRIX();
void Free_MATRIX();
int  Free_ATOMs();
void Cal_Distance_MATRIX();
void Set_Constant_tFactor();
void Renumber_Atom_Number();
void Renumber_Residue_Number();
void Write_PDB_File();
void Write_PDB_File_Residue();
void Add_Atom_To_Residue_Head();
void Add_Atom_To_Residue_Tail();
void Transform_Atoms_By_Rmat_Gorig_Gnew();
float CorrCoeff_Bwn_Atoms_and_GMM();
char Judge_Atom_or_Residue_Model();


/*** Functions (LOCAL) ***/
static void Exclude_Doubling_altLoc_Atoms();
static int Find_Same_Res_Rnum_Atom_Before();
static void Change_N_CA_C_O_Hetatom_to_Atom();
static void Renumber_Atom_num();
static void Exclude_Hetero_Atoms();



void Read_PDB_File(ifname,Ahead,AtomSelect,ChainID)
 char *ifname;
 struct ATOM *Ahead;
 char  AtomSelect;  /* 'A'll atom except hydrogen, 'R'esidue-based ('CA' and 'P') */
 char  ChainID; /* Chain Identifier () */
{
 FILE *fp;
 char line[256],command[256],B[32],Aname[32],Rname[32],atomtype,chain,altloc;
 char atminfo[32],atminfo0[32],end;
 int accept,anum,L,model_num;
 struct ATOM *an;


 L = strlen(ifname);
 if ((L>3)&&(ifname[L-3]=='.')&&(ifname[L-2]=='g')&&(ifname[L-1]=='z')){
   sprintf(command,"zcat %s",ifname);
   fp = popen(command,"r");
 }
 else{
   fp = fopen(ifname,"r");
 }
 if (fp==NULL){printf("#ERROR:Can't open pdbfile \"%s\"\n",ifname); exit(1);} 

 model_num = 0;

 anum = 0;  
 Ahead->next = NULL;
 an = Ahead; end = 0;

 while ((feof(fp)==0)&&(end==0)){
   line[0] = '\0';
   fgets(line,255,fp);
   /* printf("%s",line); */
        if (strncmp(line,"ATOM",4)==0)   {atomtype = 'A';}
   else if (strncmp(line,"HETATM",6)==0) {atomtype = 'H';}
   else { atomtype = ' '; }

   if ((PAR.MODEL=='S')&&(strncmp(line,"ENDMDL",6)==0)){end = 1;}

   if (strncmp(line,"MODEL",5)==0){
/*
          1
012345678901234
MODEL        1
MODEL        2
MODEL       15
*/
      Get_Part_Of_Line(B,line,5,13);
      model_num = atoi(B);  
   }

   if (atomtype != ' '){ 
     accept = 1; 
     Get_Part_Of_Line(Aname,line,12,15);
     Get_Part_Of_Line(Rname,line,17,19);
     Get_Part_Of_Line(atminfo,line,12,25);
     atminfo[4] = ' '; 
    /* printf("\"%s\"\n",atminfo); */
     altloc = line[16]; 
    
     chain = line[21]; 

    /*** Decide Accept or Not ***/

 /*
    if ((PAR.HETtype == 'F') && (atomtype == 'H')) accept = 0;
  */  
    if (Aname[1]=='H') accept = 0;

      if (strncmp(Rname,"HOH",3)==0) accept = 0;
      if (strncmp(Rname,"DOD",3)==0) accept = 0;

      if (AtomSelect == 'R'){ 
        if ((strncmp(Aname," CA ",4)==0)||(strncmp(Aname," P  ",4)==0)){accept = 1;} 
        else {accept = 0;}
      }

      if ((ChainID!='-')&&(chain!=ChainID)){ accept = 0; }

      if (accept==1){ 
        an->next = (struct ATOM *)malloc(sizeof(struct ATOM));
        an->next->prev = an;
        an = an->next;
        an->next = NULL;
        Get_Part_Of_Line(B,line,30,37); an->Pos[0] = atof(B);
        Get_Part_Of_Line(B,line,38,45); an->Pos[1] = atof(B);
        Get_Part_Of_Line(B,line,46,53); an->Pos[2] = atof(B);
        Get_Part_Of_Line(an->Anum,line,6,10);
        Get_Part_Of_Line(an->Atom,line,12,15);
        Get_Part_Of_Line(an->Resi,line,17,19);
        Get_Part_Of_Line(an->Rnum,line,22,26);
        Get_Part_Of_Line(B,line,54,59); an->Occup  = atof(B);
        an->tFactor = 0.000;
        Get_Part_Of_Line(B,line,60,65); an->tFactor = atof(B);
        an->Chain  = chain; 
        an->altLoc = altloc; 
        an->num   = anum;
        an->model_num = model_num;
        an->res = NULL;
        an->rnext = an->rprev = NULL;
        an->AHtype = atomtype;
        ++anum; 
      }
     sprintf(atminfo0,"%s",atminfo); 
    }

 } /* while */

 fclose(fp);

 Exclude_Doubling_altLoc_Atoms(Ahead);
 Change_N_CA_C_O_Hetatom_to_Atom(Ahead);
 
 if (PAR.HETtype == 'F'){ Exclude_Hetero_Atoms(Ahead);}

} /* end of Read_PDB_File() */








void Write_PDB_File(fname,Ahead,mode,new_chainID,comment,command)
 char *fname;
 struct ATOM *Ahead;
 char mode; /* 'a'ppend, 'w'rite */
 char new_chainID;  /* if '-', use the original chainID */ 
 char *comment;
 char *command;
{
 FILE *fp;
 struct ATOM *an;
 int i,j,L,model_num_pre;
 char line[1024],subline[100],chainID;

 printf("#Write_PDB_File(newchain '%c') -> \"%s\"\n",new_chainID,fname);
 if (fname[0]=='-') fp = stdout;
 else{
    if (mode=='a') fp = fopen(fname,"a");
             else  fp = fopen(fname,"w");
  if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1);}
 }

 sprintf(line,"FILENAME \"%s\"",fname);
 Get_Part_Of_Line(subline,line,0,40);
 fprintf(fp,"HEADER    %-40s%s   0XXX      0XXX   1\n",subline,Get_Date_String_PDB());
 
 /*** Write Comments ***/
 if (comment[0] != '\0'){
  L = strlen(comment);
  i = j = 0;
  fprintf(fp,"REMARK #   ");
  while (i<L){ 
    if ((j>60)||(comment[i]=='\n')) { fprintf(fp,"\nREMARK #   "); j = 0;}
    if (comment[i]!='\n') fprintf(fp,"%c",comment[i]);
   ++i; ++j;
  } 
  fprintf(fp,"\n");
 }

 fprintf(fp,"REMARK    DATE %s\n",Get_Date_String());
 if (command[0] != '\0'){ fprintf(fp,"REMARK    COMMAND %s\n",command); }
 fprintf(fp,"REMARK    Occupancy [55-60] : Radius of the atom\n");

 /*** Write ATOM/HETATM ***/
 model_num_pre = 0;
 an = Ahead;
 while (an->next != NULL){
   an = an->next;
   if (an->model_num != model_num_pre){
     if (model_num_pre!=0){fprintf(fp,"ENDMDL\n");}
     fprintf(fp,"MODEL %8d\n",an->model_num);
   }
        if (an->AHtype=='A')  fprintf(fp,"ATOM  ");
   else if (an->AHtype=='H')  fprintf(fp,"HETATM");
   else fprintf(fp,"ATOM? ");
   if (new_chainID=='-'){chainID = an->Chain;}
                   else {chainID = new_chainID;}
  fprintf(fp,"%5s %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
        an->Anum,an->Atom,an->Resi,chainID,an->Rnum,
        an->Pos[0],an->Pos[1],an->Pos[2],an->R,an->tFactor);
  model_num_pre = an->model_num;
 }

 fprintf(fp,"TER\n");
 fprintf(fp,"ENDMDL\n");

 if (fp!=stdout) fclose(fp);

} /* end of Write_PDB_File() */
















void Write_PDB_File_Residue(fname,Rhead,mode,comment,command)
 char *fname;
 struct RESIDUE *Rhead;
 char mode; /* 'a'ppend, 'w'rite */
 char *comment;
 char *command;
{
 FILE *fp;
 struct RESIDUE *rn;
 struct ATOM *an;
 int i,j,L;
 char line[512],subline[100];

 printf("#Write_PDB_File_Residue -> \"%s\"\n",fname);
 if (fname[0]=='-') fp = stdout;
 else
 {
    if (mode=='a') fp = fopen(fname,"a");
             else  fp = fopen(fname,"w");
  if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1);}
 }

 sprintf(line,"FILENAME \"%s\"",fname);
 Get_Part_Of_Line(subline,line,0,40);
 fprintf(fp,"HEADER    %-40s%s   0XXX      0XXX   1\n",subline,Get_Date_String_PDB());

 if (comment[0] != '\0'){
  L = strlen(comment);
  i = j = 0;
  fprintf(fp,"REMARK #   ");
  while (i<L)
  { if ((j>60)||(comment[i]=='\n')) { fprintf(fp,"\nREMARK #   "); j = 0;}
    if (comment[i]!='\n') fprintf(fp,"%c",comment[i]);
   ++i; ++j;
  } 
  fprintf(fp,"\n");
 }

 if (command[0] != '\0') fprintf(fp,"REMARK COMMAND %s\n",command);
 rn = Rhead;
 while (rn->next != NULL){
  rn = rn->next;
  an = &(rn->Ahead);
  while (an->rnext != NULL){
   an = an->rnext;
        if (an->AHtype=='A')  fprintf(fp,"ATOM  ");
   else if (an->AHtype=='H')  fprintf(fp,"HETATM");
   else fprintf(fp,"ATOM? ");
   fprintf(fp,"%5s %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
        an->Anum,an->Atom,an->Resi,an->Chain,an->Rnum,
        an->Pos[0],an->Pos[1],an->Pos[2],an->R,an->tFactor);
  } /* an */

 } /* rn */

 fprintf(fp,"TER\n");
 if (fp!=stdout) fclose(fp);

} /* end of Write_PDB_File_Residue() */









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


int Number_Of_Atom(Ahead)
 struct ATOM *Ahead;
{ struct ATOM *an;
  int Natom;
  Natom = 0; an = Ahead;
  while (an->next != NULL){ 
    an = an->next; 
    ++Natom;
  }
  return(Natom);
} /* end of  Number_Of_Atom() */


int Number_Of_Residue(Rhead)
 struct RESIDUE *Rhead;
{ struct RESIDUE *rn;
  int Nres;

  Nres = 0; rn = Rhead;
  while (rn->next != NULL){
    rn = rn->next; ++Nres;
  }
  return(Nres);

} /* end of  Number_Of_Residue() */





int Number_Of_Chain(Ahead)
 struct ATOM *Ahead;
{
 struct ATOM *an;
 int Natom,Nchain;

 Natom = Nchain = 0;
 an = Ahead;
 while (an->next != NULL){
   an = an->next;
   if (Natom==0){ Nchain = 1;}

   if (Natom!=0){
     if (an->Chain != an->prev->Chain){
       ++Nchain;
     }
   }
   ++Natom;
 }

 return(Nchain);

} /* end of  Number_Of_Chain() */



char Chain_Of_Atoms(Ahead)
 struct ATOM *Ahead;
{
 struct ATOM *an;
 char chain;

 an = Ahead;
 if (an->next==NULL) return(' ');
 an = an->next;
 chain = an->Chain;
 while (an->next != NULL) {
  an = an->next;
  if (an->Chain != chain) return('X');
  }
 return(chain);

} /* end of  Number_Of_Chain() */






void Malloc_MATRIX(M,Nrow,Ncol)
 struct MATRIX *M;
 int Nrow,Ncol;
{ int i;
  double Mbyte;

  M->Nrow = Nrow;
  M->Ncol = Ncol;
/*
  Mbyte = (double)sizeof(float)*M->Nrow * M->Ncol/1000.0/1000.0;
  */ 
  Mbyte = (double)sizeof(double)*M->Nrow * M->Ncol/1000.0/1000.0;
  printf("#Malloc_MATRIX(Nrow %d Ncol %d) --> %.2f Mbyte\n",M->Nrow,M->Ncol,Mbyte);
  if (Mbyte > MAX_MEMORY_MEGA){
   printf("#ERROR:Malloc_MATRIX():%.2f Mbyte is larger than MAX (%d Mbyte).\n",
      Mbyte,MAX_MEMORY_MEGA);
   exit(1);  
  }  
 
 /* 
  M->m = (float **)malloc(sizeof(float *)*Nrow);
  for (i=0;i<Nrow;++i){
     M->m[i] = (float *)malloc(sizeof(float)*Ncol);
  }
 */
  M->m = (double **)malloc(sizeof(double *)*Nrow);
  for (i=0;i<Nrow;++i){
     M->m[i] = (double *)malloc(sizeof(double)*Ncol);
  }

} /* end of Malloc_MATRIX() */



void Free_MATRIX(M)
 struct MATRIX *M;
{ int i;
 
 for (i=0;i<M->Nrow;++i) free(M->m[i]);
 free(M->m);

} /* end of Free_MATRIX() */


int Free_ATOMs(Ahead)
 struct ATOM *Ahead;
{
 /* return "Number of free atom" */
 struct ATOM *an,*pn;
 int Nfree;

 printf("#Free_ATOMs()\n");
 if (Ahead->next == NULL) return(0);
                                                                                                        
 an = Ahead; pn = NULL; Nfree = 0;
 while (an->next != NULL){
   pn = an;
   an = an->next;
   if ((pn != Ahead)&&(pn != NULL)) { free(pn); pn = NULL; ++Nfree;}
  }
                                                                                                        
 if (an!=NULL) { free(an); an = NULL; ++Nfree;}
 Ahead->next = NULL;
 return(Nfree);

} /* end of Free_ATOMs() */





void Make_Residue_List(Ahead,Rhead)
 struct ATOM    *Ahead;
 struct RESIDUE *Rhead;
{
 struct ATOM *an;
 struct RESIDUE *rn;
 int  Nresidue;

 printf("#Make_Residue_List(Ahead,Rhead)\n");
 Rhead->next = NULL;

 an = Ahead;
 rn = Rhead;
 Nresidue = 0;
 while (an->next != NULL){
  an = an->next;
  
  if ((an->prev==NULL)||
      (strncmp(an->Rnum,an->prev->Rnum,5)!=0)||(an->Chain != an->prev->Chain)){
    rn->next = (struct RESIDUE *)malloc(sizeof(struct RESIDUE));
    rn->next->prev = rn;
    rn->next->next = NULL;
    rn = rn->next;
    rn->Chain = an->Chain;
    sprintf(rn->Resi,"%s",an->Resi);
    sprintf(rn->Rnum,"%s",an->Rnum);
    rn->Natom = 0; 
    ++Nresidue;  
    rn->num = Nresidue;
    rn->Ahead.rnext = NULL;  
}

   Add_Atom_To_Residue_Tail(rn,an); 
 }

 printf("#Nresidue %d\n",Nresidue);

} /* end of Make_Residue_List() */




void Add_Atom_To_Residue_Tail(rn,an)
 struct RESIDUE *rn;
 struct ATOM    *an;
{
 struct ATOM *bn;

 /*
 printf("#Add_Atom_To_Residue_Tail(rnum %d anum %d)\n",rn->num,an->num); fflush(stdout);
 */
 
 bn = &(rn->Ahead);
 while (bn->rnext != NULL) bn = bn->rnext;
 bn->rnext = an;
 an->rnext = NULL; 
 ++(rn->Natom);
 an->res = rn;
 an->rnum = rn->num; 
 sprintf(an->Rnum,"%s",rn->Rnum);
 sprintf(an->Resi,"%s",rn->Resi);

} /* end of Add_Atom_To_Residue_Tail() */




void Add_Atom_To_Residue_Head(rn,an)
 struct RESIDUE *rn;
 struct ATOM    *an;
{
 struct ATOM *bn;

 bn = rn->Ahead.rnext;
 rn->Ahead.rnext = an; 
 if (bn != NULL) { an->rnext = bn; } 
 ++(rn->Natom);
 an->res = rn;
 an->rnum = rn->num; 
 sprintf(an->Rnum,"%s",rn->Rnum);
 sprintf(an->Resi,"%s",rn->Resi);

} /* end of Add_Atom_To_Residue_Head() */







 

void Set_Constant_tFactor(Phead,val)
 struct ATOM *Phead;
 float val;
{  
 struct ATOM *an;
   
 an = Phead;
 while (an->next != NULL){
   an = an->next;
   an->tFactor = val;
 } 

} /* end of Set_Constant_tFactor() */



void Renumber_Atom_Number(Phead)
 struct ATOM *Phead;
{  
 struct ATOM *an;
 int Natom;
   
 an = Phead;
 Natom = 0;
 while (an->next != NULL){
   an = an->next;
   ++Natom;
   an->num = Natom;
   sprintf(an->Anum,"%5d",an->num);
   if (an->res != NULL) 
    { sprintf(an->Rnum,"%s",an->res->Rnum);
      an->rnum = an->res->num;}
   else 
   { sprintf(an->Rnum,"%4d ",an->num);
     an->rnum = an->num;}
 } 

} /* end of Renumber_Atom_Number() */



void Renumber_Residue_Number(Rhead)
 struct RESIDUE *Rhead;
{  
 struct RESIDUE *rn;
 int Nres;
   
 rn = Rhead;
 Nres = 0;
 while (rn->next != NULL){
   rn = rn->next;
   ++Nres;
   rn->num = Nres;
   sprintf(rn->Rnum,"%4d ",rn->num);
 } 

} /* end of Renumber_Residue_Number() */





void Exclude_Doubling_altLoc_Atoms(HeadAtom)
 struct ATOM *HeadAtom;
{
 struct ATOM     *an;
 char hit;
 int Nex;

 Nex = 0;
 an = HeadAtom;
 while (an->next != NULL){
   an = an->next;
   if (an->altLoc != ' '){
     hit = Find_Same_Res_Rnum_Atom_Before(HeadAtom,an);
     if (hit==1){ 
       ++Nex;
       if (an->prev != NULL)  an->prev->next   = an->next;
       if (an->next != NULL)  an->next->prev   = an->prev;
     }
   }

 } /* an */

  Renumber_Atom_num(HeadAtom);

} /* end of Exclude_Doubling_altLoc_Atoms() */




void Exclude_Hetero_Atoms(HeadAtom)
 struct ATOM *HeadAtom;
{
 struct ATOM     *an;
 int Nex;

 Nex = 0;
 an = HeadAtom;
 while (an->next != NULL){
   an = an->next;
   if (an->AHtype == 'H'){ 
     ++Nex;
     if (an->prev != NULL)  an->prev->next   = an->next;
     if (an->next != NULL)  an->next->prev   = an->prev; }

 } /* an */
    
  Renumber_Atom_num(HeadAtom);

} /* end of Exclude_Hetero_Atoms() */




void Renumber_Atom_num(HeadAtom)
 struct ATOM *HeadAtom;
{
 int Natom;
 struct ATOM     *an;

 Natom = 0;
 an = HeadAtom;
 while (an->next != NULL)
 { an = an->next; 
   an->num = Natom;
   ++Natom; }

} /* end of Renumber_Atom_num() */




int Find_Same_Res_Rnum_Atom_Before(HeadAtom,focus)
 struct ATOM  *HeadAtom;
 struct ATOM  *focus;
{
 struct ATOM *an;
 char hit;

 hit = 0;
 an = HeadAtom;
 while ((an->next != NULL)&&(hit==0)&&(an->next->num!=focus->num)){
  an = an->next;
  if ((an->altLoc != ' ')&&(an->num != focus->num)){
    if ( (focus->Chain == an->Chain) &&
         (strncmp(focus->Resi, an->Resi,3)==0) &&
         (strncmp(focus->Rnum, an->Rnum,5)==0) &&
         (strncmp(focus->Atom, an->Atom,4)==0)) { hit = 1; return(hit);}
   }

 } /* an */

  return(hit);
} /* end of Find_Same_Res_Rnum_Atom_Before() */







void Change_N_CA_C_O_Hetatom_to_Atom(HeadAtom)
 struct ATOM *HeadAtom;
{
 struct ATOM *atom,*satom;
 char Nok,Cok,CAok,Ook,stop;
 int Nhetatm;

 atom = HeadAtom;
 Nok = Cok = CAok = 0; Nhetatm = 0; satom = NULL;
 while (atom->next != NULL){
  if ((strcmp(atom->Resi,atom->next->Resi)!=0)||(strcmp(atom->Rnum,atom->next->Rnum)!=0)){

   if (((Nok*Cok*CAok*Ook)==1)&&(Nhetatm>0)){
     stop = 0;
     while ((satom!=NULL)&&(satom->num <= atom->num)){
      if (satom->AHtype != 'A') { satom->AHtype = 'A'; }
      satom = satom->next; }
     }

   satom = NULL;
   Nok = Cok = CAok = 0; Nhetatm = 0;
  }

  atom = atom->next;

  if (satom==NULL) satom = atom;
  if (atom->AHtype != 'A') ++Nhetatm;
  if (strncmp(atom->Atom," CA ",4)==0) CAok = 1;
  if (strncmp(atom->Atom," N  ",4)==0) Nok  = 1;
  if (strncmp(atom->Atom," C  ",4)==0) Cok  = 1;
  if (strncmp(atom->Atom," O  ",4)==0) Ook  = 1;
 }

} /* end of  Change_N_CA_C_O_Hetatom_to_Atom() */




void Transform_Atoms_By_Rmat_Gorig_Gnew(HeadAtom, Rmat,gorig, gref)
  struct ATOM *HeadAtom;
  double Rmat[3][3];
  double gorig[3], gref[3];
/*
 
  Xnew = R(Xorig-gorig) + gref 
 
 */
{
  struct ATOM *an;
  double Rx_gorig[3],x_gorig[3];
  int i;

  an = HeadAtom;
  while (an->next != NULL){
    an = an->next;
    for (i=0;i<3;++i){
      x_gorig[i] = an->Pos[i]-gorig[i];
    }
    Multiply_Matrix3D_Vec3D(Rx_gorig,Rmat,x_gorig);
    for (i=0;i<3;++i){
      an->Pos[i] = Rx_gorig[i] + gref[i];
    }
  }

} /* end of Transform_Atoms_By_Rmat_Gorig_Gnew() */




char Judge_Atom_or_Residue_Model(Ahead)
 struct ATOM *Ahead;
{ struct ATOM *an;
  int Natom,nCA,nP;
  float rCA_P;
  char atmrestype;
 
  Natom = 0;
  nCA = 0;
  nP  = 0;

  an = Ahead;
  while (an->next != NULL){ 
    an = an->next; 
    ++Natom; 
         if (strncmp(an->Atom," CA ",4)==0) {nCA += 1;}
    else if (strncmp(an->Atom," P  ",4)==0) {nP  += 1;}
  }

  rCA_P = (float)(nCA + nP)/(float)Natom;
  if (rCA_P>0.9) {atmrestype = 'R';}
              else {atmrestype = 'A';}
  printf("#Natom %d nCA %d nP %d rCA_P %f ==> Judged AtmResType '%c'\n",Natom,nCA,nP,rCA_P,atmrestype);
  return(atmrestype);

} /* end of Judge_Atom_or_Residue_Model() */





