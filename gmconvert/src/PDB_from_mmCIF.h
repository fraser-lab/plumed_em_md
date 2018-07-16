/*
  <PDB_from_mmCIF.h>

*/


struct UNITMOL{
  char *asym_id; 
  char *auth_asym_id; 
  int    Natom;   /* in other words, number of ROWLINE from head to tail_rowline_pointer */
  struct ROWLINE *head_rowline_pointer; /* pointer to rowline of the TABLE 'atom_site' */
  struct ROWLINE *tail_rowline_pointer; /* pointer to rowline of the TABLE 'atom_site' */
  struct UNITMOL *next;
};



struct ASMBLMOL{
  char *asym_id; 
  char *oper_expression; 
  int    Natom;   /* equal to the Natom of the corresponding UNITMOL. */
  float *Cartn_x; /* transformed x-coordinate.[Natom] (malloc later) */ 
  float *Cartn_y; /* transformed x-coordinate.[Natom] (malloc later) */ 
  float *Cartn_z; /* transformed x-coordinate.[Natom] (malloc later) */ 
  struct ASMBLMOL *next; 
};


struct OPERATION{
  char *id;       /* malloc later (corresponds to oper_expression)*/
  char *type;     /* malloc later */
  float R[3][3];
  float T[3];
  struct OPERATION *next;
};


struct ASSEMBLY{
  char *id;                     /* (malloc later) */
  char *details;                /* (malloc later) */
  char *method_details;         /* (malloc later) */
  char *oligomeric_details;     /* (malloc later) */
  int  Nasmblmol;               /* number of ASMBLMOL */
  int  Natom;                   /* total number of atoms of the all the ASMBLMOL */ 
  char **OPER_EXPRESSION_LIST;  /* array of oper_expression [Nasmblmol] */
  char **ASYM_ID_LIST;          /* array of asym_id         [Nasmblmol] */
  struct ASSEMBLY *next;
};


struct PDB_DATA{
  char *pdb_id;
  int Natom;    /* number of atoms in all the UNITMOLs */
  int Nassembly;
  struct ASSEMBLY head_assembly;
  int Noperation;
  struct OPERATION head_operation;
  struct UNITMOL  head_unitmol;
  struct ASMBLMOL head_asmblmol;
};

extern int make_ASSEMBLY_from_BLOCK();
extern int make_OPERATION_from_BLOCK();
extern void make_UNITMOL_from_BLOCK();
extern void make_ASMBLMOL_from_BLOCK();
extern void write_ASSEMBLY_in_pdb();
extern void write_UNITMOLs_in_pdb();
extern struct ASSEMBLY* find_assembly();
extern struct UNITMOL* find_unitmol();
extern struct ASMBLMOL* find_asmblmol();
extern void setup_column_numbers_for_ATOM_HETATM_line();
extern void free_PDB_DATA();
extern int  expand_oper_expression_string();
extern void split_into_words();
extern void get_part_of_line();

