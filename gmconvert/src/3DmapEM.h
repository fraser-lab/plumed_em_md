/*

 <3DmapEM.h>

*/

/*** Functions (GLOBAL) ***/
extern void Set_Less_Than_Threshold_Voxels_to_Zero();
extern void Cal_Mean_and_CovarMat_from_Voxel();
extern int  EM_optimize_GaussMix_for_3D_Map_SmallMemory();
extern int  EM_optimize_GaussMix_for_3D_Map();
extern double Corr_Coeff_Bwn_Two_Voxels();
extern double Corr_Coeff_Bwn_Voxel_and_GMM();
extern double Log_Likelihood_Of_GaussMix_For_3Dmap();
extern void Set_GaussMix_Density_To_3Dmap();
