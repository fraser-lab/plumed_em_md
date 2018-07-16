/*

 <PdbIO.h>

*/

extern void Read_PDB_File();
extern void Make_Residue_List();
extern void Write_PDB_File();
extern void Write_PDB_File_Residue();
extern void Find_Filename_Core();
extern void Get_Part_Of_Line();
extern char *Get_Date_String();
extern char *Get_Date_String_PDB();
extern int Number_Of_Atom();
extern int Number_Of_Residue();
extern int Number_Of_Chain();
extern char Chain_Of_Atoms();
extern void Malloc_MATRIX();
extern void Free_MATRIX();
extern int  Free_ATOMs();
extern void Cal_Distance_MATRIX();
extern void Set_Region_Of_Atoms();
extern void Set_Constant_tFactor();
extern void Renumber_Atom_Number();
extern void Renumber_Residue_Number();
extern void Add_Atom_To_Residue_Head();
extern void Add_Atom_To_Residue_Tail();
extern void Transform_Atoms_By_Rmat_Gorig_Gnew();
extern float CorrCoeff_Bwn_Atoms_and_GMM();
extern char Judge_Atom_or_Residue_Model();


