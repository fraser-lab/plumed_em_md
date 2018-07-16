/*

<pdbstructs.h>

Definition of basic structures

*/

struct RADIUS{
 char   Atom[5];
 char   Resi[4];
 float  R;
 struct RADIUS *next,*prev;
};




/** ATOM : for both protein atoms and probe spheres */
struct ATOM{
 float Pos[3];        /* X,Y,Z coordinates */
 float R;             /* Radius (van der Waals)   */ 
 float Weight;        /* Atomic Weight     */ 
 /* float RR;             Radius^2         */ 
 int   num;           /* Atom Number (0,1,2,....) */
 int   rnum;          /* Residue Number */
 char  Anum[7];       /* Atom Number String */
 char  Atom[5];       /* Name of Atom */
 char  Resi[5];       /* Name of Residue */
 char  Rnum[6];       /* Residue Number String */
 char  Chain;         /* Chain ID */
 char  altLoc;        /* Alternate location indicator */
 float Occup;         /* Occupancy */ 
 float tFactor;       /* Temparature Factor */ 
 char  AHtype;        /* 'A'tom, 'H'etatm */
 int   model_num;     /* MODEL number     */
 struct ATOM *next,*prev;    /* Pointer for atom list (bi-direction) */
 struct RESIDUE *res;        /* Pointer to corresponding residues */
 struct ATOM *rnext,*rprev;  /* Pointer for the list from RESIDUE (bi-direction)*/
 char   region;
 char   mark;           /* Mark for various purpose */
 float  variance;     /* variance of Atom GDF */
};



/** RESIDUE : for both protein residues and probe clusters */
struct RESIDUE{
 int    num;           /* Residue Number */
 char   Resi[5];       /* Name of Residue */
 char   Rnum[6];       /* Residue Number String */
 char   Chain;         /* Chain ID */
 int    Natom;         /* Number of Atoms */ 
 float  Vprobe;        /* Volume of Contacted Probes */
 struct RESIDUE *next,*prev; /* Pointer for RESIDUE list (bi-direction) */
 struct ATOM    Ahead; /* Head of ATOM list, belonged to the residue (using ATOM->rnext, rprev) */
 float  value;         /* Value for various purpose such as Sorting */  
};




/** MATRIX : simple 2D matrix for general purpose **/
struct MATRIX{
 int   Nrow;  /* Number of rows            */
 int   Ncol;  /* Number of columns         */
 double **m;  /* 2D valiable: M[0..Nrow-1][1..Ncol-1]   */
};

