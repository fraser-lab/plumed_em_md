/*

 <PointEM.h>

*/

/*** Functions (GLOBAL) ***/
extern void Initialize_Gauss_From_Randomly_Chosen_Atoms();
extern void Initialize_Gauss_From_Mol_Distribution();
extern void EM_optimize_GaussMix_For_Atoms();
extern void Malloc_Voxel_From_Covariance_Matrix();
extern void Cal_Mean_and_CovarMat_from_Atoms();
extern void Cal_Memberships_for_Atoms();
extern void Write_PDB_File_With_GMM_and_Memberships();
extern void Read_PDB_File_With_GMM_and_Memberships();
extern void Estimate_M_and_CovM_from_Memberships();
