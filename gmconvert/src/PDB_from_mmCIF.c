/*
  <PDB_from_mmCIF.c>

  functions for changing "struct BLOCK" to "struct PDB_DATA".
 
   LastModDate: Feb 4, 2015

  >> order for performing the functions <<
    make_ASSEMBLY_from_BLOCK(); => make_UNITMOL_from_BLOCK(); => make_OPERATION_from_BLOCK(); => make_ASMBLMOL_from_BLOCK();

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "io_mmCIF.h"
#include "PDB_from_mmCIF.h"


#define MAX_ASYM_ID 10000

/* FUNCTIONS (GLOBAL) */
int  make_ASSEMBLY_from_BLOCK();
int  make_OPERATION_from_BLOCK();
void make_UNITMOL_from_BLOCK();
void make_ASMBLMOL_from_BLOCK();
void write_ASSEMBLY_in_pdb();
void write_UNITMOLs_in_pdb();
struct ASSEMBLY* find_assembly();
struct UNITMOL* find_unitmol();
struct ASMBLMOL* find_asmblmol();
void setup_column_numbers_for_ATOM_HETATM_line();
void free_PDB_DATA();
int expand_oper_expression_string();
void split_into_words();
void get_part_of_line();

/* FUNCTIONS (LOCAL) */
static struct OPERATION* find_operation();
static void get_Rmat_and_Tvec_from_oper_expression();
static void free_ASSEMBLY();
static void transform_pos3D_by_Rmat_and_Tvec();
static void write_ATOM_HETATM_line();



int make_ASSEMBLY_from_BLOCK(P,blk)
  struct PDB_DATA *P;
  struct BLOCK *blk;
{
/*
  >> order for performing the functions <<
    make_ASSEMBLY_from_BLOCK(); => make_UNITMOL_from_BLOCK(); => make_OPERATION_from_BLOCK(); => make_ASMBLMOL_from_BLOCK();

>> examples for mmCIF file <<
data_1A5K
#
_pdbx_struct_assembly.id                   1
_pdbx_struct_assembly.details              author_and_software_defined_assembly
_pdbx_struct_assembly.method_details       PISA,PQS
_pdbx_struct_assembly.oligomeric_details   nonameric
_pdbx_struct_assembly.oligomeric_count     9
#
_pdbx_struct_assembly_gen.assembly_id       1
_pdbx_struct_assembly_gen.oper_expression   1,2,3
_pdbx_struct_assembly_gen.asym_id_list      A,D,B,E,C,F
#
loop_
_pdbx_struct_assembly_prop.biol_id
_pdbx_struct_assembly_prop.type
_pdbx_struct_assembly_prop.value
_pdbx_struct_assembly_prop.details
1 'ABSA (A^2)' 47440 ?
1 'SSA (A^2)'  55450 ?
1 MORE         -241  ?
#
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.name
_pdbx_struct_oper_list.symmetry_operation
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 'identity operation'         1_555 x,y,z 1.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 1.0000000000
0.0000000000 0.0000000000 0.0000000000 0.0000000000 1.0000000000 0.0000000000
2 'crystal symmetry operation' 5_555 z,x,y 0.0000000000 0.0000000000 1.0000000000 0.0000000000 1.0000000000 0.0000000000
0.0000000000 0.0000000000 0.0000000000 1.0000000000 0.0000000000 0.0000000000
3 'crystal symmetry operation' 9_555 y,z,x 0.0000000000 1.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
1.0000000000 0.0000000000 1.0000000000 0.0000000000 0.0000000000 0.0000000000
#

------------------------------------------------------------
data_2BUK
#
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
_pdbx_struct_assembly.oligomeric_details
_pdbx_struct_assembly.oligomeric_count
1   'complete icosahedral assembly'                60-MERIC 60
2   'icosahedral asymmetric unit'                  ?        ?
3   'icosahedral pentamer'                         ?        ?
4   'icosahedral 23 hexamer'                       ?        ?
PAU 'icosahedral asymmetric unit, std point frame' ?        ?
XAU 'crystal asymmetric unit, crystal frame'       ?        ?
#
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1   '(1-60)'           A,B,C,D,E
2   1                  A,B,C,D,E
3   '(1-5)'            A,B,C,D,E
4   '(1,2,6,10,23,24)' A,B,C,D,E
PAU P                  A,B,C,D,E
XAU '(X0)(1-60)'       A,B,C,D,E
#
#
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
P  'transform to point frame'   -0.46431533 -0.80595401 -0.36721848 51.39772  -0.85378970 0.29707252  0.42754073  43.44080
-0.23548765 0.51204107  -0.82605055 55.69698
X0 'transform to crystal frame' 1.00000000  0.00000000  0.00000000  0.00000   0.00000000  1.00000000  0.00000000  0.00000
0.00000000  0.00000000  1.00000000  0.00000
1  'point symmetry operation'   1.00000000  0.00000000  0.00000000  0.00000   0.00000000  1.00000000  0.00000000  -0.00000
0.00000000  0.00000000  1.00000000  0.00000
2  'point symmetry operation'   0.60022259  0.18907595  0.77716351  -6.37893  -0.71996118 0.55097588  0.42199701  33.78484
-0.34840886 -0.81281970 0.46683552  50.49749
3  'point symmetry operation'   -0.04663084 -0.41402987 0.90906811  35.42500  -0.97584571 -0.17556041 -0.13001419 78.30184
0.21342613  -0.89317289 -0.39584274 48.83300
:
59 'point symmetry operation'   -0.67599716 0.39670953  -0.62100676 152.89994 0.54001836  0.84009633  -0.05116959 -37.62950
0.50140604  -0.36994554 -0.78213316 45.39144
60 'point symmetry operation'   -0.47499997 0.59552910  -0.64785810 139.25557 -0.26287759 0.60656895  0.75031292  -15.27565
0.83980379  0.52670599  -0.13156898 -9.80131
#


*/
  struct TABLE *tb;
  struct ROWLINE *rn; 
  struct ASSEMBLY *an;

  int i,j,k;
  int NwordO,WstaO[MAX_ASYM_ID],WendO[MAX_ASYM_ID];
  int NwordA,WstaA[MAX_ASYM_ID],WendA[MAX_ASYM_ID];
  char *oper_string, *asym_string,wordO[1000],wordA[1000];
  int Loper_string;

  /* printf("#make_ASSEMBLY_from_BLOCK()\n"); */
  tb = TABLE_from_key(blk->table_hash_tab,"entry");
  if (tb==NULL){
    printf("#WARNING : Can't find the table '%s'\n","entry");
    return(0); }
  P->pdb_id = (char *)malloc(sizeof(char)*(strlen(data_from_rowline(tb,tb->head_rowline.next,"id"))+1));
  sprintf(P->pdb_id,"%s",data_from_rowline(tb,tb->head_rowline.next,"id"));

  /* (1) From the table "pdbx_struct_assembly" */
  tb = TABLE_from_key(blk->table_hash_tab,"pdbx_struct_assembly");
  if (tb==NULL){
    printf("#WARNING: Can't find the table '%s'\n","pdbx_struct_assembly");
    return(0);
  }
  P->Nassembly = tb->Nrow;
  /* printf("#Nassembly %d\n",P->Nassembly);  */
  P->head_assembly.next = NULL;

  an = &(P->head_assembly);
  rn = &(tb->head_rowline);
  i = 0; 
  while (rn->next != NULL){
    rn = rn->next;
    an->next = (struct ASSEMBLY*)malloc(sizeof(struct ASSEMBLY));
    an = an->next;
    an->next = NULL;
    an->id = (char *)malloc(sizeof(char)*((int)strlen(data_from_rowline(tb,rn,"id"))+1));
    sprintf(an->id,"%s",data_from_rowline(tb,rn,"id"));
  
    an->details = (char *)malloc(sizeof(char)*((int)strlen(data_from_rowline(tb,rn,"details"))+1));
    sprintf(an->details,"%s",data_from_rowline(tb,rn,"details"));
 
    an->method_details = (char *)malloc(sizeof(char)*((int)strlen(data_from_rowline(tb,rn,"method_details"))+1));
    sprintf(an->method_details,"%s",data_from_rowline(tb,rn,"method_details"));
 
    an->oligomeric_details = (char *)malloc(sizeof(char)*((int)strlen(data_from_rowline(tb,rn,"oligomreric_details"))+1));
    sprintf(an->oligomeric_details,"%s",data_from_rowline(tb,rn,"oligomeric_details"));
    an->Nasmblmol = 0;
    an->Natom     = 0;
    i += 1;
  }

/*
  an = &(P->head_assembly);
  while (an->next != NULL){
    an = an->next;
    printf("#'%s' '%s' '%s' '%s'\n",an->id,an->details,an->method_details,an->oligomeric_details);
  }
*/

  /* (2) From the table "pdbx_struct_assembly_gen" */
  tb = TABLE_from_key(blk->table_hash_tab,"pdbx_struct_assembly_gen");
  an = &(P->head_assembly);
  while (an->next != NULL){
    an = an->next;
    an->Nasmblmol = 0;
    /* (2-1) just count an->Nasmblmol */
    rn = &(tb->head_rowline);
    while (rn->next != NULL){
      rn = rn->next;
      if (strcmp(an->id,data_from_rowline(tb,rn,"assembly_id"))==0){

        oper_string = NULL;
        /* printf("oper_expression '%s'\n",data_from_rowline(tb,rn,"oper_expression")); */
        Loper_string = expand_oper_expression_string(oper_string,data_from_rowline(tb,rn,"oper_expression"),"length_only");
        oper_string  = (char *)malloc(sizeof(char)*(Loper_string+1));
        expand_oper_expression_string(oper_string,data_from_rowline(tb,rn,"oper_expression"),"");
        split_into_words(oper_string,',',&NwordO,WstaO,WendO,MAX_ASYM_ID);
        free(oper_string);

        asym_string = NULL;
        asym_string  = (char *)malloc(sizeof(char)*((int)strlen(data_from_rowline(tb,rn,"asym_id_list"))+1));
        sprintf(asym_string,"%s",data_from_rowline(tb,rn,"asym_id_list"));
        split_into_words(asym_string,',',&NwordA,WstaA,WendA,MAX_ASYM_ID);
        free(asym_string);

        an->Nasmblmol += (NwordO * NwordA);
      }
    }
    /* printf("#Nasmblmol %d\n",an->Nasmblmol); */
    /* (2-2) malloc an->OPER_EXPRESSION_LIST and an->ASYM_LIST */
    an->OPER_EXPRESSION_LIST = (char **)malloc(sizeof(char *)*an->Nasmblmol);  
    an->ASYM_ID_LIST = (char **)malloc(sizeof(char *)*an->Nasmblmol);  

    /* (2-3) set up an->OPER_EXPRESSION_LIST and an->ASYM_LIST */
    rn = &(tb->head_rowline);
    k = 0;
    while (rn->next != NULL){
      rn = rn->next;
      if (strcmp(an->id,data_from_rowline(tb,rn,"assembly_id"))==0){
        oper_string = NULL;
        Loper_string = expand_oper_expression_string(oper_string,data_from_rowline(tb,rn,"oper_expression"),"length_only");
        oper_string  = (char *)malloc(sizeof(char)*(Loper_string+1));
        expand_oper_expression_string(oper_string,data_from_rowline(tb,rn,"oper_expression"),"");
        split_into_words(oper_string,',',&NwordO,WstaO,WendO,MAX_ASYM_ID);

        asym_string = NULL;
        asym_string  = (char *)malloc(sizeof(char)*((int)strlen(data_from_rowline(tb,rn,"asym_id_list"))+1));
        sprintf(asym_string,"%s",data_from_rowline(tb,rn,"asym_id_list"));
        split_into_words(asym_string,',',&NwordA,WstaA,WendA,MAX_ASYM_ID);
        /* 
        printf("#asym_string '%s' NwordA %d  oper_string '%s' NwordO %d\n",asym_string,NwordA,oper_string,NwordO); 
        */
        for (i=0;i<NwordO;++i){
          get_part_of_line(wordO,oper_string,WstaO[i],WendO[i]);
          for (j=0;j<NwordA;++j){
            get_part_of_line(wordA,asym_string,WstaA[j],WendA[j]);
            an->OPER_EXPRESSION_LIST[k] = (char *)malloc(sizeof(char)*(strlen(wordO)+1));
            sprintf(an->OPER_EXPRESSION_LIST[k],"%s",wordO);
            an->ASYM_ID_LIST[k]         = (char *)malloc(sizeof(char)*(strlen(wordA)+1));
            sprintf(an->ASYM_ID_LIST[k],"%s",wordA);
            /* printf("#k %d asym_id '%s' oper '%s'\n",k,an->ASYM_ID_LIST[k],an->OPER_EXPRESSION_LIST[k]); */
            k += 1;
          }
        }
        free(asym_string);
        free(oper_string);
      }
    }
  }

  /* printf("#Nasmblmol %d\n",an->Nasmblmol); */
  /*
  for (k=0;k<an->Nasmblmol;++k){
    printf("k %d asym_id '%s' oper '%s'\n",k,an->ASYM_ID_LIST[k],an->OPER_EXPRESSION_LIST[k]);
  }
  */

  return(1);
} /* end of make_ASSEMBLY_from_BLOCK() */



int make_OPERATION_from_BLOCK(P,blk)
  struct PDB_DATA *P;
  struct BLOCK *blk;
{
/*
  >> order for performing the functions <<
    make_ASSEMBLY_from_BLOCK(); => make_UNITMOL_from_BLOCK(); => make_OPERATION_from_BLOCK(); => make_ASMBLMOL_from_BLOCK();
>> example taken from "2buk"  <<
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
P  'transform to point frame'   -0.46431533 -0.80595401 -0.36721848 51.39772  -0.85378970 0.29707252  0.42754073  43.44080
-0.23548765 0.51204107  -0.82605055 55.69698
X0 'transform to crystal frame' 1.00000000  0.00000000  0.00000000  0.00000   0.00000000  1.00000000  0.00000000  0.00000
0.00000000  0.00000000  1.00000000  0.00000
1  'point symmetry operation'   1.00000000  0.00000000  0.00000000  0.00000   0.00000000  1.00000000  0.00000000  -0.00000
0.00000000  0.00000000  1.00000000  0.00000
2  'point symmetry operation'   0.60022259  0.18907595  0.77716351  -6.37893  -0.71996118 0.55097588  0.42199701  33.78484
-0.34840886 -0.81281970 0.46683552  50.49749
*/
  struct TABLE *tb;
  struct ROWLINE *rn; 
  struct OPERATION *on; 
  int i;
 
  /* printf("#make_OPERATION_from_BLOCK()\n"); fflush(stdout); */

  tb = TABLE_from_key(blk->table_hash_tab,"pdbx_struct_oper_list");
  if (tb==NULL){
    printf("#WARNING: Can't find the table '%s'\n","pdbx_struct_oper_list");
    return(0);
  }
  P->Noperation = tb->Nrow;
  /* printf("#P->Noperation %d\n",P->Noperation); fflush(stdout); */
 
  rn = &(tb->head_rowline);
  P->head_operation.next = NULL;
  on = &(P->head_operation); 
  i = 0;
  while (rn->next != NULL){
    rn = rn->next;
    on->next = (struct OPERATION*)malloc(sizeof(struct OPERATION));
    on = on->next;
    on->next = NULL;
    on->id = (char *)malloc(sizeof(char)*(strlen(data_from_rowline(tb,rn,"id"))+1));
    sprintf(on->id,"%s",data_from_rowline(tb,rn,"id"));
 /*
    printf("type '%s' strlen(type) %d\n",data_from_rowline(tb,rn,"type"),strlen(data_from_rowline(tb,rn,"type"))); fflush(stdout);
 */
    on->type = (char *)malloc(sizeof(char)*(strlen(data_from_rowline(tb,rn,"type"))+1));
    sprintf(on->type,"%s",data_from_rowline(tb,rn,"type"));
    on->R[0][0] = atof(data_from_rowline(tb,rn,"matrix[1][1]"));
    on->R[0][1] = atof(data_from_rowline(tb,rn,"matrix[1][2]"));
    on->R[0][2] = atof(data_from_rowline(tb,rn,"matrix[1][3]"));
    on->R[1][0] = atof(data_from_rowline(tb,rn,"matrix[2][1]"));
    on->R[1][1] = atof(data_from_rowline(tb,rn,"matrix[2][2]"));
    on->R[1][2] = atof(data_from_rowline(tb,rn,"matrix[2][3]"));
    on->R[2][0] = atof(data_from_rowline(tb,rn,"matrix[3][1]"));
    on->R[2][1] = atof(data_from_rowline(tb,rn,"matrix[3][2]"));
    on->R[2][2] = atof(data_from_rowline(tb,rn,"matrix[3][3]"));

    on->T[0] = atof(data_from_rowline(tb,rn,"vector[1]"));
    on->T[1] = atof(data_from_rowline(tb,rn,"vector[2]"));
    on->T[2] = atof(data_from_rowline(tb,rn,"vector[3]"));

    i += 1;
  }

  return(1);

} /* end of make_OPERATION_from_BLOCK() */




void make_UNITMOL_from_BLOCK(P,blk,ExHOH)
  struct PDB_DATA *P;
  struct BLOCK    *blk;
  char  ExHOH;            /* if ExHOH = 'T', then reject  'HOH' and 'DOD'. */ 
{
/*
  >> order for performing the functions <<
    make_ASSEMBLY_from_BLOCK(); => make_UNITMOL_from_BLOCK(); => make_OPERATION_from_BLOCK(); => make_ASMBLMOL_from_BLOCK();
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x_esd
_atom_site.Cartn_y_esd
_atom_site.Cartn_z_esd
_atom_site.occupancy_esd
_atom_site.B_iso_or_equiv_esd
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1   N N   . MET A 1 1  ? 27.340 24.430 2.614  1.00 9.67  ? ? ? ? ? ? 1   MET A N   1
ATOM   2   C CA  . MET A 1 1  ? 26.266 25.413 2.842  1.00 10.38 ? ? ? ? ? ? 1   MET A CA  1
:
*/
  struct TABLE *tb;
  struct ROWLINE *rn,*rn_prev,*sn;
  struct UNITMOL *um;
  int cnum_label_asym_id, cnum_auth_asym_id,cnum_label_comp_id;
  char asym_id[16], asym_id_pre[16],auth_asym_id[16];
  
  /* printf("#make_UNITMOL_from_BLOCK(ExHOH '%c')\n",ExHOH); */
  P->Natom = 0;
  tb = TABLE_from_key(blk->table_hash_tab,"atom_site");   
  rn = &(tb->head_rowline);
  rn_prev = NULL;
  sn = NULL; 
  P->head_unitmol.next = NULL;
  um = &(P->head_unitmol);
  cnum_label_asym_id  = number_from_column(tb,"label_asym_id");
  cnum_auth_asym_id   = number_from_column(tb,"auth_asym_id");
  cnum_label_comp_id  = number_from_column(tb,"label_comp_id");

  asym_id_pre[0] = '\0';
  while (rn->next != NULL){
    rn_prev = rn;
    rn      = rn->next;
    sprintf(asym_id,"%s",rn->data[cnum_label_asym_id]);
    if ((ExHOH!='T') || ((strcmp(rn->data[cnum_label_comp_id],"HOH")!=0)&&(strcmp(rn->data[cnum_label_comp_id],"DOD")!=0))){
      if (strcmp(asym_id,asym_id_pre)!=0){
        um->next = (struct UNITMOL*)malloc(sizeof(struct UNITMOL));
        um = um->next;
        um->next = NULL;
        um->asym_id = (char *)malloc(sizeof(char)*(strlen(asym_id)+1));
        sprintf(um->asym_id,"%s",asym_id);
  
        sprintf(auth_asym_id,"%s",rn->data[cnum_auth_asym_id]);
        um->auth_asym_id = (char *)malloc(sizeof(char)*(strlen(auth_asym_id)+1));
        sprintf(um->auth_asym_id,"%s",auth_asym_id);

        um->Natom = 0;
        um->head_rowline_pointer = rn_prev;
      }
      um->tail_rowline_pointer = rn->next;
      um->Natom += 1;
      P->Natom += 1;
    } 

    sprintf(asym_id_pre,"%s",asym_id);
  }


} /* end of make_UNITMOL_from_BLOCK() */



void make_ASMBLMOL_from_BLOCK(P,blk,assembly_id)
  struct PDB_DATA *P;
  struct BLOCK *blk;
  char *assembly_id;
{
/*
  >> order for performing the functions <<
    make_ASSEMBLY_from_BLOCK(); => make_UNITMOL_from_BLOCK(); => make_OPERATION_from_BLOCK(); => make_ASMBLMOL_from_BLOCK();
*/

  struct ASSEMBLY *A;
  struct ASMBLMOL *am;
  struct UNITMOL   *um;
  struct TABLE     *tb;
  struct ROWLINE   *rn;
  int i,k;
  int cnum_Cartn_x, cnum_Cartn_y, cnum_Cartn_z;
  float Rmat[3][3],Tvec[3],pos[3],tra_pos[3];

  /* printf("#make_ASMBLMOL_from_ASSEMBLY(assembly_id '%s')\n",assembly_id); */

  P->head_asmblmol.next = NULL;
  am = &(P->head_asmblmol);
  tb = TABLE_from_key(blk->table_hash_tab,"atom_site");   
  cnum_Cartn_x  = number_from_column(tb,"Cartn_x");
  cnum_Cartn_y  = number_from_column(tb,"Cartn_y");
  cnum_Cartn_z  = number_from_column(tb,"Cartn_z");
  

  A = &(P->head_assembly);
  while (A->next != NULL){
    A = A->next;
    A->Natom = 0;
    if ((assembly_id[0]=='\0')||(strcmp(A->id,assembly_id)==0)){
      for (i=0;i<A->Nasmblmol;++i){
        get_Rmat_and_Tvec_from_oper_expression(Rmat,Tvec,P,A->OPER_EXPRESSION_LIST[i]);
        um = find_unitmol(P,A->ASYM_ID_LIST[i]);
        /* If ExHOH='T', then "water" um should be NULL!! */
        if (um!=NULL){
          am->next = (struct ASMBLMOL*)malloc(sizeof(struct ASMBLMOL));
          am = am->next;
          am->Natom  = um->Natom;
          A->Natom  += um->Natom;
          am->next = NULL;
          am->asym_id = (char *)malloc(sizeof(char)*(strlen(um->asym_id)+1));
          sprintf(am->asym_id,"%s",um->asym_id);
          am->oper_expression = (char *)malloc(sizeof(char)*(strlen(A->OPER_EXPRESSION_LIST[i])+1));
          sprintf(am->oper_expression,"%s",A->OPER_EXPRESSION_LIST[i]);
          am->Cartn_x = (float *)malloc(sizeof(float)*um->Natom);
          am->Cartn_y = (float *)malloc(sizeof(float)*um->Natom);
          am->Cartn_z = (float *)malloc(sizeof(float)*um->Natom);
          rn = um->head_rowline_pointer;
          k = 0;
          while ((rn->next != NULL)&&(rn->next != um->tail_rowline_pointer)){
            rn = rn->next;
            pos[0] = atof(rn->data[cnum_Cartn_x]);
            pos[1] = atof(rn->data[cnum_Cartn_y]);
            pos[2] = atof(rn->data[cnum_Cartn_z]);
            transform_pos3D_by_Rmat_and_Tvec(tra_pos,pos,Rmat,Tvec);
            am->Cartn_x[k] = tra_pos[0];
            am->Cartn_y[k] = tra_pos[1];
            am->Cartn_z[k] = tra_pos[2];
            k += 1;
          }
        }
      }
    }
  }

} /*end of make_ASMBLMOL_from_BLOCK */





void transform_pos3D_by_Rmat_and_Tvec(tra_pos,orig_pos,Rmat,Tvec)
/*
  tra_pos = Rmat * orig_pos + Tvec;
*/ 
  float tra_pos[3];
  float orig_pos[3];
  float Rmat[3][3];
  float Tvec[3];
{
  int i,j;
  for (i=0;i<3;++i){
    tra_pos[i] = Tvec[i];
    for (j=0;j<3;++j){
      tra_pos[i] += Rmat[i][j] * orig_pos[j];
    }
  }

} /* end of transform_pos3D_by_Rmat_and_Tvec() */



struct ASSEMBLY* find_assembly(P,assembly_id)
  struct PDB_DATA *P;
  char *assembly_id;
{
  struct ASSEMBLY *A;
  A = &(P->head_assembly);
  while (A->next != NULL){
    A = A->next;
    if (strcmp(A->id,assembly_id)==0){ return(A); }
  }
  return(NULL);
} /* end of find_assembly() */



struct UNITMOL* find_unitmol(P,asym_id)
  struct PDB_DATA *P;
  char *asym_id;
{
  struct UNITMOL *um;
  um = &(P->head_unitmol);
  while (um->next != NULL){
    um = um->next;
    if (strcmp(um->asym_id,asym_id)==0){ return(um); }
  }
  return(NULL);
} /* end of find_unitmol() */


struct ASMBLMOL* find_asmblmol(P,asym_id,oper_expression)
  struct PDB_DATA *P;
  char *asym_id;
  char *oper_expression;
{
  struct ASMBLMOL *am;
  
  am = &(P->head_asmblmol);
  while (am->next != NULL){
    am = am->next;
    if ((strcmp(am->asym_id,asym_id)==0)&&(strcmp(am->oper_expression,oper_expression)==0)) { 
      return(am);
    }
  }
  return(NULL);
} /* end of find_asmblmol() */


struct OPERATION* find_operation(P,id)
  struct PDB_DATA *P;
  char   *id;
{
  struct OPERATION *op;
  op = &(P->head_operation);
  while (op->next != NULL){
    op = op->next;
    if (strcmp(op->id,id)==0){ return(op); }
  }
  return(NULL);
} /* end of find_operation() */





void get_Rmat_and_Tvec_from_oper_expression(Rmat,Tvec,P,oper_expression)
  float Rmat[3][3];
  float Tvec[3];
  struct PDB_DATA *P;
  char   *oper_expression;
{
  struct OPERATION *op1,*op2;
  int i,j,k,hyphen_pos;
  char  oper_exp1[100],oper_exp2[100];
 
 /* (1) Single oper_expression */
  if (strstr(oper_expression,"-")==NULL){
    op1 = find_operation(P,oper_expression);
    for (i=0;i<3;++i){
      Tvec[i] = op1->T[i];
      for (j=0;j<3;++j){
         Rmat[i][j] = op1->R[i][j];
      }
    }
  }
 /* (2) Double oper_expression */
/*
 '[t2]-[t1]' = (R2,T2)-(R2,T1)

 x1  = R1*x0 + T1
 x2  = R2*x1 + T2
     = R2*(R1*x0 + T1) + T2
     = R2*R1*x0 + R2*T1 + T2
     = R*    x0 + T
 Threfore,
    R = R2*R1,   T = R2*T1 + T2
*/
  else{
    hyphen_pos = strstr(oper_expression,"-") - oper_expression;
    get_part_of_line(oper_exp2,oper_expression,0,hyphen_pos-1);
    get_part_of_line(oper_exp1,oper_expression,hyphen_pos+1,strlen(oper_expression)-1);
    op2 = find_operation(P,oper_exp2);
    op1 = find_operation(P,oper_exp1);
    /*  R = R2*R1 */
    for (i=0;i<3;++i){
      for (j=0;j<3;++j){
        Rmat[i][j] = 0.0;
        for (k=0;k<3;++k){
          Rmat[i][j] += op2->R[i][k] * op1->R[k][j];
        }
     }
   }
   /* T = R2*T1 + T2 */
    for (i=0;i<3;++i){
      Tvec[i] = op2->T[i];
      for (k=0;k<3;++k){
        Tvec[i] += op2->R[i][k] * op1->T[k];
      }
    }
  }

} /* end of get_Rmat_and_Tvec_from_oper_expression() */







void write_UNITMOLs_in_pdb(ofname,P,blk,residue_select,given_auth_asym_id,TransformType,Rmat,Tvec,My_tFactor)
  char *ofname;
  struct PDB_DATA *P;
  struct BLOCK *blk;
  char   residue_select;    /* if residue_select = 'T', then write only ' CA ' andn ' P '. */ 
  char   *given_auth_asym_id;
  char   TransformType;    /* if TransformType=='T', transform atomic XYZ by Rmat[][] and Tvec[] */
  float  Rmat[3][3],Tvec[3];
  float  My_tFactor;
{
  FILE *fpo;
  struct TABLE *tb;
  struct ROWLINE *rn;
  int group_PDB,id,auth_atom_id,label_alt_id,auth_comp_id,auth_asym_id,auth_seq_id,pdbx_PDB_ins_code;
  int Cartn_x,Cartn_y,Cartn_z,occupancy,B_iso_or_equiv,type_symbol,label_asym_id;
  struct UNITMOL  *um;
  int model_num,accept;
  float pos[3],tra_pos[3],tFactor;
 
  printf("#write_UNITMOLs_in_pdb(residue_select:%c auth_asym_id:%s TransformType:%c )\n",residue_select,given_auth_asym_id,TransformType); 
  if (ofname[0]=='-'){
    fpo = stdout;
  }
  else{
    fpo = fopen(ofname,"w");
    printf("#write_UNITMOLs_in_pdb()-->'%s'\n",ofname);
  } 

  tb = TABLE_from_key(blk->table_hash_tab,"atom_site");   
  if (tb==NULL){
     printf("#ERROR:Can't find table '%s'.\n","atom_site");
     exit(1); 
  }

  setup_column_numbers_for_ATOM_HETATM_line(
    tb, &group_PDB,&id, &auth_atom_id,&label_alt_id, &auth_comp_id,&auth_asym_id,&auth_seq_id,&pdbx_PDB_ins_code,
    &Cartn_x, &Cartn_y, &Cartn_z,&occupancy, &B_iso_or_equiv,&type_symbol,&label_asym_id);

  /** write REMARK lines **/
  fprintf(fpo,"REMARK     PDB_ID '%s'\n",P->pdb_id);
  fprintf(fpo,"REMARK     Natom  %d\n",P->Natom);
  fprintf(fpo,"REMARK     Nassembly  %d\n",P->Nassembly);
  um = &(P->head_unitmol);
  model_num = 0;
  while (um->next != NULL){
    um = um->next;
    model_num += 1;
    fprintf(fpo,"REMARK     MODEL_NUMBER %4d asym_id %2s auth_asym_id %2s Natom %6d\n", model_num,um->asym_id,um->auth_asym_id,um->Natom);
  }  

  /** write ATOM/HETATM lines **/
  um = &(P->head_unitmol);
  model_num = 0;
  while (um->next != NULL){
    um = um->next;
    rn = um->head_rowline_pointer; 

    model_num += 1;
    fprintf(fpo,"MODEL    %5d\n",model_num);

    while ((rn->next != NULL)&&(rn->next != um->tail_rowline_pointer)){
      rn = rn->next;
/*
      printf("#residue_select '%c'\n",residue_select);
      printf("#'%s' '%s' '%s' '%s' '%s' '%s'\n",rn->data[auth_atom_id],rn->data[type_symbol],rn->data[Cartn_x],rn->data[Cartn_y],rn->data[Cartn_z],rn->data[5]); 
*/
/*         
  >> residue select reject condition <<
    not [(atom_id='CA' and type='C') or (atom_id='P' and type='P')]
 =  not [(atom_id='CA' and type='C')] and not[(atom_id='P' and type='P')]
 =  (atom_id!='CA' or type!='C')] and (atom_id!='P' or type!='P')
*/ 
/* 
  >> residue select reject condition <<
    not [(atom_id='CA' and (type='C' or auth_atom_id='CA')) or (atom_id='P' and (type='P' or auth_atom_id='P')]
    not [(atom_id='CA' and (type='C' or auth_atom_id='CA'))] and not [(atom_id='P' and (type='P' or auth_atom_id='P')]
    (atom_id!'CA') or not (type='C' or auth_atom_id='CA')) and (atom_id!='P' or not(type='P' or auth_atom_id='P')
    (atom_id!'CA') or (type!='C' and auth_atom_id!='CA')) and (atom_id!='P' or (type!='P' and auth_atom_id!='P')
*/

      accept = 1;
      if (residue_select=='T'){   
        if ( ((strcmp(rn->data[auth_atom_id],"CA")!=0)|| (((strcmp(rn->data[type_symbol],"C")!=0)) && (strcmp(rn->data[auth_atom_id],"CA")!=0))) &&
             ((strcmp(rn->data[auth_atom_id],"P")!=0) || (((strcmp(rn->data[type_symbol],"P")!=0)) && (strcmp(rn->data[auth_atom_id],"P")!=0))) ){
          accept = 0; 
         }
      }

      if ((given_auth_asym_id[0] !='\0') && (given_auth_asym_id[0] != '-')){
        if (strcmp(rn->data[auth_asym_id],given_auth_asym_id)!=0){ accept = 0;}
      }

      if (accept == 1){
        pos[0] = atof(rn->data[Cartn_x]);
        pos[1] = atof(rn->data[Cartn_y]);
        pos[2] = atof(rn->data[Cartn_z]);
        tFactor = atof(rn->data[B_iso_or_equiv]);
        if (TransformType=='T'){
          transform_pos3D_by_Rmat_and_Tvec(tra_pos,pos,Rmat,Tvec);
          pos[0] = tra_pos[0];
          pos[1] = tra_pos[1];
          pos[2] = tra_pos[2];
        }

        if (My_tFactor>=0.0){ tFactor = My_tFactor; }
        write_ATOM_HETATM_line(fpo,
          rn->data[group_PDB], rn->data[id], rn->data[auth_atom_id], rn->data[label_alt_id], 
          rn->data[auth_comp_id], rn->data[auth_asym_id], rn->data[auth_seq_id], rn->data[pdbx_PDB_ins_code],
          pos[0],pos[1],pos[2],
          rn->data[occupancy], tFactor,rn->data[type_symbol]);
      }
    } 
    fprintf(fpo,"ENDMDL\n");
  }

  if (ofname[0]!='-'){
    fclose(fpo);
  }

} /* end of write_UNITMOLs_in_pdb() */








void write_ASSEMBLY_in_pdb(ofname,P,blk,assembly_id,residue_select,TransformType,Rmat,Tvec,My_tFactor)
  char *ofname;
  struct PDB_DATA *P;
  struct BLOCK *blk;
  char *assembly_id;
  char residue_select;    /* if residue_select = 'T', then write only ' CA ' and ' P '. */ 
  char  TransformType;    /* if TransformType=='T', transform atomic XYZ by Rmat[][] and Tvec[] */
  float Rmat[3][3],Tvec[3];
  float My_tFactor;
{
  FILE *fpo;
  struct TABLE *tb;
  struct ROWLINE *rn;
  int group_PDB,id,auth_atom_id,label_alt_id,auth_comp_id,auth_asym_id,auth_seq_id,pdbx_PDB_ins_code;
  int Cartn_x,Cartn_y,Cartn_z,occupancy,B_iso_or_equiv,type_symbol,label_asym_id;
  struct ASSEMBLY *A;
  struct UNITMOL  *um;
  struct ASMBLMOL *am;
  int i,k,model_num;
  float pos[3],tra_pos[3],tFactor;
  char accept;

  if (ofname[0]=='-'){
    fpo = stdout;
  }
  else{
    printf("#write_ASSEMBLY_in_pdb(residue_select %c TransformType %c)-->'%s'\n",residue_select,TransformType,ofname);
    fpo = fopen(ofname,"w");
  }

  tb = TABLE_from_key(blk->table_hash_tab,"atom_site");
  if (tb==NULL){
     printf("#ERROR:Can't find table '%s'.\n","atom_site");
     exit(1);
  }

  setup_column_numbers_for_ATOM_HETATM_line(
    tb, &group_PDB,&id, &auth_atom_id,&label_alt_id, &auth_comp_id,&auth_asym_id,&auth_seq_id,&pdbx_PDB_ins_code,
    &Cartn_x, &Cartn_y, &Cartn_z,&occupancy, &B_iso_or_equiv,&type_symbol,&label_asym_id);
/** write REMARK lines **/
  fprintf(fpo,"REMARK     PDB_ID '%s'\n",P->pdb_id);
  A = &(P->head_assembly);
  model_num = 0;
  while (A->next != NULL){
    A = A->next;
    if ((assembly_id[0]=='\0')||(strcmp(A->id,assembly_id)==0)){
      for (i=0;i<A->Nasmblmol;++i){
        um = find_unitmol(P,A->ASYM_ID_LIST[i]);
        am = find_asmblmol(P,A->ASYM_ID_LIST[i],A->OPER_EXPRESSION_LIST[i]);
        if ((um!=NULL)&&(am!=NULL)){
          model_num += 1;
          fprintf(fpo,"REMARK     MODEL %4d asym_id %2s oper %3s auth_asym_id %2s Natom %6d\n",
             model_num,am->asym_id,am->oper_expression,um->auth_asym_id,um->Natom);
        }
      }
    }
  }

 /** write ATOM/HETATM lines **/
  A = &(P->head_assembly);
  model_num = 0;
  while (A->next != NULL){
    A = A->next;
    if ((assembly_id[0]=='\0')||(strcmp(A->id,assembly_id)==0)){
      for (i=0;i<A->Nasmblmol;++i){
        um = find_unitmol(P,A->ASYM_ID_LIST[i]);
        am = find_asmblmol(P,A->ASYM_ID_LIST[i],A->OPER_EXPRESSION_LIST[i]);
        if ((um!=NULL)&&(am!=NULL)){
          rn = um->head_rowline_pointer;
          model_num += 1;
          fprintf(fpo,"MODEL    %5d\n",model_num);
          k = 0;
          while ((rn->next != NULL)&&(rn->next != um->tail_rowline_pointer)){
            rn = rn->next;
/*         
  >> residue select reject condition <<
    not [(atom_id='CA' and type='C') or (atom_id='P' and type='P')]
 =  not [(atom_id='CA' and type='C')] and not[(atom_id='P' and type='P')]
 =  (atom_id!='CA' or type!='C')] and (atom_id!='P' or type!='P')
*/ 
            accept = 1;
            if (residue_select=='T'){   
              if (((strcmp(rn->data[auth_atom_id],"CA")!=0)||(strcmp(rn->data[type_symbol],"C")!=0)) && 
                  ((strcmp(rn->data[auth_atom_id],"P")!=0) ||(strcmp(rn->data[type_symbol],"P")!=0))  ){
                accept = 0; 
               }
            }

            if (accept == 1){
              pos[0] = am->Cartn_x[k];
              pos[1] = am->Cartn_y[k];
              pos[2] = am->Cartn_z[k];

              tFactor = atof(rn->data[B_iso_or_equiv]);
              if (TransformType=='T'){
                transform_pos3D_by_Rmat_and_Tvec(tra_pos,pos,Rmat,Tvec);
                pos[0] = tra_pos[0];
                pos[1] = tra_pos[1];
                pos[2] = tra_pos[2];
              }
              if (My_tFactor>=0.0){ tFactor = My_tFactor; }
              write_ATOM_HETATM_line(fpo,
                rn->data[group_PDB], rn->data[id], rn->data[auth_atom_id], rn->data[label_alt_id],
                rn->data[auth_comp_id], rn->data[auth_asym_id], rn->data[auth_seq_id], rn->data[pdbx_PDB_ins_code],
                pos[0],pos[1],pos[2],
                rn->data[occupancy],tFactor,rn->data[type_symbol]);
           }
          k += 1;
        }
        fprintf(fpo,"ENDMDL\n");
        }
       }
     }
   }

  if (ofname[0]!='-'){
    fclose(fpo);
  }

} /* end of write_ASSEMBLY_in_pdb() */







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




int expand_oper_expression_string(output_string,input_string,length_only)
  char *output_string; /* expanded oper_expression_string considering '-' and ')(' */
  char *input_string;  /* input oper_expression_string */
  char *length_only;   /* if length_only[] = 'length_only', return only length */
{
/*
 >> example <<
 '1'                -> '1'
 'P'                -> 'P'
  '1,2,3'           -> '1,2,3'
 '(1-5)'            -> '1,2,3,4,5'
 '(1-60)'           -> '1,2,...,59,60'
  taken from PDBcode:'2buk'
 '(1,2,6,10,23,24)' -> '1,2,6,10,23,24'
 '(X0)(1-60)'       -> 'X0-1,X0-2,...,X0-59,X0-60'
  taken from PDBcode:'1m4x'
 '(1-60)(61-88)'       -> '1-61,1-62,...,60-87,60-88'
# >> caution <<
# '[t2]-[t1]' means that firstly, the transformation 't1' is done,
#  next, the transformation 't2' is done.
*/
  int i,j,k,f,L,hyphen_pos,blackets_pos;
  char field[2][1000],oper_string[2][1000],instr[1000];
  char buff0[1000],buff1[1000],buff2[1000],buff3[1000];
  int Nword0,Wsta0[1000],Wend0[1000];
  int Nword1,Wsta1[1000],Wend1[1000];
  int Nfield,Loutput_string,start,end; 

  /* printf("#expand_oper_expression_string('%s')\n",input_string); */
  /** (1) stript the first '(' and the last ')' **/

  L = strlen(input_string);
  sprintf(instr,"%s",input_string);
  if (instr[0]=='(') {
    for (i=1;i<L;++i){
      instr[i-1] = instr[i];
    }
    instr[L] = '\0';
    L = L -1;
  } 
  if (instr[L-1] == ')'){ 
    instr[L-1] = '\0'; L = L - 1;
  }
 
  /** (2) split into field[0] and field[1] by ')(' **/
  if (strstr(instr,")(")==NULL){
    Nfield = 1;
    sprintf(field[0],"%s",instr);
  }
  else{
    Nfield = 2;
    blackets_pos = strstr(instr,")(") - instr;
    get_part_of_line(field[0],instr,0,blackets_pos-1);
    get_part_of_line(field[1],instr,blackets_pos+2,strlen(instr)-1);
  }
 
  /* printf("#instr '%s' Nfield %d\n",instr,Nfield); */

  /** (3) field[f] --> oper_string[f] considering '[num1]-[num2]' patterns. **/
  for (f=0;f<Nfield;++f){
    j = 0;
    oper_string[f][0] = '\0';
    split_into_words(field[f],',',&Nword0,Wsta0,Wend0,1000); 

    for (i=0;i<Nword0;++i){
      get_part_of_line(buff0,field[f],Wsta0[i],Wend0[i]);
      if (strstr(buff0,"-")==NULL){
        if (oper_string[f][0]!='\0'){ strcat(oper_string[f],","); }
        strcat(oper_string[f],buff0);
      }
      else{
      /* special case for hyphen-included instr, such as '1-60'. */
        hyphen_pos = strstr(buff0,"-")-buff0;
        get_part_of_line(buff1,buff0,0,hyphen_pos-1);
        get_part_of_line(buff2,buff0,hyphen_pos+1,strlen(buff0)-1);
        start = atoi(buff1);
        end   = atoi(buff2);
        for (k=start;k<=end;++k){
          sprintf(buff3,"%d",k);
          if (oper_string[f][0]!='\0'){ strcat(oper_string[f],","); }
          strcat(oper_string[f],buff3);
        }
      }
    }
  }
  
  /** (4) if (Nfield=1), just return oper_string[0]. **/
  if (Nfield==1){
    if (strcmp(length_only,"length_only")!=0){
      sprintf(output_string,"%s",oper_string[0]);
    }
    return(strlen(oper_string[0]));
  }
  /** (5) if (Nfield=2), combine oper_string[0] and oper_string[1] to generate output_string[]. **/
  else if (Nfield==2){
    split_into_words(oper_string[0],',',&Nword0,Wsta0,Wend0,1000); 
    split_into_words(oper_string[1],',',&Nword1,Wsta1,Wend1,1000); 
    /* cal Loutput_string */ 
    Loutput_string = 0;
    for (i=0;i<Nword0;++i){
      get_part_of_line(buff0,oper_string[0],Wsta0[i],Wend0[i]);
      for (j=0;j<Nword1;++j){
        get_part_of_line(buff1,oper_string[1],Wsta1[j],Wend1[j]);
        Loutput_string += (Wend0[i]-Wsta0[i]+1) + (Wend1[j]-Wsta1[j]+1) + 2;
      }
    }
    if (strcmp(length_only,"length_only")==0){
      return(Loutput_string);    
    }
    /* make output_string[] */ 
    output_string[0] = '\0';
    for (i=0;i<Nword0;++i){
      get_part_of_line(buff0,oper_string[0],Wsta0[i],Wend0[i]);
      for (j=0;j<Nword1;++j){
        get_part_of_line(buff1,oper_string[1],Wsta1[j],Wend1[j]);
        if (output_string[0]!='\0'){ strcat(output_string,","); }
        sprintf(buff3,"%s-%s",buff0,buff1);
        strcat(output_string,buff3);
        /* printf("output_string '%s' buff3 '%s'\n",output_string,buff3); */
      }
    }
   return(Loutput_string);    
  }

  
 return(0);

} /* end of expand_oper_expression_string() */




void setup_column_numbers_for_ATOM_HETATM_line(
  tb,
  group_PDB, id, auth_atom_id, label_alt_id, auth_comp_id, auth_asym_id, auth_seq_id, pdbx_PDB_ins_code,
  Cartn_x, Cartn_y, Cartn_z,occupancy, B_iso_or_equiv, type_symbol,label_asym_id
)
  struct TABLE *tb;
/* these are column number to be calculated */
  int *group_PDB,*id,*auth_atom_id,*label_alt_id,*auth_comp_id,*auth_asym_id,*auth_seq_id,*pdbx_PDB_ins_code;
  int *Cartn_x, *Cartn_y, *Cartn_z,*occupancy,*B_iso_or_equiv,*type_symbol,*label_asym_id;
{
/*
>> Items of TABLE "atom_site" <<
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x_esd
_atom_site.Cartn_y_esd
_atom_site.Cartn_z_esd
_atom_site.occupancy_esd
_atom_site.B_iso_or_equiv_esd
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1   N N   . MET A 1 1  ? 27.340 24.430 2.614  1.00 9.67  ? ? ? ? ? ? 1   MET A N   1
ATOM   2   C CA  . MET A 1 1  ? 26.266 25.413 2.842  1.00 10.38 ? ? ? ? ? ? 1   MET A CA  1
:
>>PDB to PDBx/mmCIF Data Item Correspondences <<
 [Field Name] 	[mmCIF Data Item ]	 
  Section   	  _atom_site.group_PDB   	 
  Serial_No   	  _atom_site.id   	 
  Atom_Name   	  _atom_site.auth_atom_id   	 
  Alt_Loc   	  _atom_site.label_alt_id   	 
  Residue_Name    _atom_site.auth_comp_id   	 
  Strand_ID   	  _atom_site.auth_asym_id   	 
  Residue_No   	  _atom_site.auth_seq_id   	 
  Ins_Code   	  _atom_site.pdbx_PDB_ins_code   	 
  X   	          _atom_site.Cartn_x   	 
  Y   	          _atom_site.Cartn_y   	 
  Z   	          _atom_site.Cartn_z   	 
  Occupancy   	  _atom_site.occupancy   	 
  T_Factor   	  _atom_site.B_iso_or_equiv
  Symbol   	  _atom_site.type_symbol   	 
*/

  *group_PDB     = number_from_column(tb,"group_PDB");
  *id            = number_from_column(tb,"id");
  *auth_atom_id  = number_from_column(tb,"auth_atom_id");
  *label_alt_id  = number_from_column(tb,"label_alt_id");
  *auth_comp_id  = number_from_column(tb,"auth_comp_id");
  *auth_asym_id  = number_from_column(tb,"auth_asym_id");
  *auth_seq_id   = number_from_column(tb,"auth_seq_id");
  *pdbx_PDB_ins_code   = number_from_column(tb,"pdbx_PDB_ins_code");
  *Cartn_x = number_from_column(tb,"Cartn_x"); 
  *Cartn_y = number_from_column(tb,"Cartn_y"); 
  *Cartn_z = number_from_column(tb,"Cartn_z"); 
  *occupancy = number_from_column(tb,"occupancy");
  *B_iso_or_equiv = number_from_column(tb,"B_iso_or_equiv");
  *type_symbol    = number_from_column(tb,"type_symbol");
  *label_asym_id  = number_from_column(tb,"label_asym_id");

} /* end of setup_column_numbers_for_ATOM_HETATM_line() */




void write_ATOM_HETATM_line(
  fpo, 
  group_PDB, id, auth_atom_id, label_alt_id, auth_comp_id, auth_asym_id, auth_seq_id, pdbx_PDB_ins_code,
  Cartn_x, Cartn_y, Cartn_z, occupancy, B_iso_or_equiv,type_symbol
)
  FILE *fpo;
  char *group_PDB, *id, *auth_atom_id, *label_alt_id, *auth_comp_id, *auth_asym_id,*auth_seq_id, *pdbx_PDB_ins_code;
  float Cartn_x, Cartn_y, Cartn_z;
  char *occupancy;
  float B_iso_or_equiv;
  char *type_symbol;
{
/*
>> OLD PDB FORMAT <<

COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.

>> example << 
         1         2         3         4         5         6         7
1234567890123456789012345678901234567890123456789012345678901234567890123456789
ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67           N
ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 10.38           C
ATOM     34  HG3 GLN A   2     -25.797  27.111   8.654  1.00  0.00           H
ATOM     35 HE21 GLN A   2     -28.648  26.923   9.929  1.00  0.00           H
ATOM     36 HE22 GLN A   2     -28.281  26.708  11.607  1.00  0.00           H
TER     603      GLY A  76
HETATM  604  O   HOH A  77      45.747  30.081  19.708  1.00 12.43           O
HETATM  605  O   HOH A  78      19.168  31.868  17.050  1.00 12.65           O
HETATM 2765  HN1 GDP A 180     140.538   7.338  14.020  1.00  1.32           H
HETATM 2766 HN21 GDP A 180     140.460   7.807  16.187  1.00  1.44           H
HETATM 2767 HN22 GDP A 180     139.262   7.113  17.256  1.00  1.31           H
ATOM     18  N  AASN A  15      52.924 -36.416  13.354  0.25 50.46           N
ATOM     19  CA AASN A  15      52.188 -35.814  14.462  0.25 50.58           C
ATOM     18  N  BASN A  15      52.924 -36.416  13.354  0.25 50.46           N
ATOM     19  CA BASN A  15      52.188 -35.814  14.462  0.25 50.58           C
*/
  fprintf(fpo,"%-6s",group_PDB);
  fprintf(fpo,"%5d",atoi(id)%100000);

  if (strlen(auth_atom_id)<4){
    fprintf(fpo,"  %-3s",auth_atom_id);
  }
  else{ fprintf(fpo," %-4s",auth_atom_id); }

  if ((label_alt_id[0] != '.')&&(label_alt_id[0] != '?')){
    fprintf(fpo,"%s",label_alt_id);
  }
  else{ fprintf(fpo," "); }

  fprintf(fpo,"%3s",auth_comp_id);

  fprintf(fpo," %c", auth_asym_id[0]);
  fprintf(fpo,"%4s", auth_seq_id);

  if ((pdbx_PDB_ins_code[0] != '.')&&(pdbx_PDB_ins_code[0] != '?')){
    fprintf(fpo,"%s", pdbx_PDB_ins_code);
  }
  else{ fprintf(fpo," "); }

  fprintf(fpo,"   ");

/*
  fprintf(fpo,"%8s", Cartn_x);
  fprintf(fpo,"%8s", Cartn_y);
  fprintf(fpo,"%8s", Cartn_z);
*/
  fprintf(fpo,"%8.3f", Cartn_x);
  fprintf(fpo,"%8.3f", Cartn_y);
  fprintf(fpo,"%8.3f", Cartn_z);

  fprintf(fpo,"%6s", occupancy);
  /* fprintf(fpo,"%6s", B_iso_or_equiv); */
  fprintf(fpo,"%6.2f", B_iso_or_equiv);
  fprintf(fpo,"           ");
  fprintf(fpo,"%6s", type_symbol);
  fprintf(fpo,"\n");

} /* end of write_ATOM_HETATM_line() */





void free_ASSEMBLY(A)
  struct ASSEMBLY *A;
{
  int i;

  free(A->id); 
  free(A->details); 
  free(A->method_details); 
  free(A->oligomeric_details); 

  for (i=0;i<A->Nasmblmol;++i){
    free(A->OPER_EXPRESSION_LIST[i]); 
    free(A->ASYM_ID_LIST[i]); 
  }
  free(A->OPER_EXPRESSION_LIST);
  free(A->ASYM_ID_LIST);


} /* end of free_ASSEMBLY() */


void free_PDB_DATA(P)
  struct PDB_DATA *P;
{
  struct ASSEMBLY  *A_curr,*A_next;
  struct UNITMOL   *u_curr,*u_next;
  struct ASMBLMOL  *a_curr,*a_next;
  struct OPERATION *o_curr,*o_next;

  /* printf("#free_PDB_DATA()\n"); */
  a_curr = &(P->head_asmblmol);
  a_next = a_curr->next;
  while (a_next != NULL){
    a_curr = a_next;
    a_next = a_curr->next;
    free(a_curr->asym_id);
    free(a_curr->oper_expression);
    free(a_curr->Cartn_x);
    free(a_curr->Cartn_y);
    free(a_curr->Cartn_z);
    free(a_curr);
  }

 
  A_curr = &(P->head_assembly);
  A_next = A_curr->next;
  while (A_next != NULL){
    A_curr = A_next;
    A_next = A_curr->next;
    free_ASSEMBLY(A_curr);
  }


  o_curr = &(P->head_operation);
  o_next = o_curr->next;
  while (o_next != NULL){
    o_curr = o_next;
    o_next = o_curr->next;
    free(o_curr->id);
    free(o_curr->type);
    free(o_curr);
  }

  u_curr = &(P->head_unitmol);
  u_next = u_curr->next;
  while (u_next != NULL){
    u_curr = u_next;
    u_next = u_curr->next;
    free(u_curr->asym_id);
    free(u_curr->auth_asym_id);
    free(u_curr);
  }

} /* end of free_PDB_DATA() */
