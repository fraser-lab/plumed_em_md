/*

 <gaussIO.h>

*/

/*** Functions (GLOBAL) ***/
extern double Gaussian3D_Value();
extern double Gaussian3D_Value_Table();
extern int Write_Gaussian3D_File();
extern int Write_Gaussian3D_File_with_PCaxis();
extern int Write_Gaussian3D_Center_Points();
extern char Read_Gaussian3D_File();
extern int  Number_Of_Gaussian3D_in_File();
extern double Corr_Coeff_Bwn_Two_GAUSS3D_Arrays();
extern double Overlap_Integral_Bwn_Two_GAUSS3Ds();
extern double Overlap_Integral_Bwn_Two_Isotropic_GAUSS3Ds();
extern double Overlap_Integral_Bwn_Two_GAUSS3D_Arrays();
extern void Malloc_Voxel_From_Gaussians();
extern void Set_Voxel_Value_By_Gaussians();
extern void Set_Voxel_Value_By_Gaussian();
extern void Transform_Gaussians_By_Rmat_Gorig_Gnew();
extern int  Check_NaN_Inf_Gaussian3D();
extern void show_Gaussian3D();
extern void Initialize_Gauss_From_M_and_CovM();
extern void Delete_ZeroWeight_gdfs_of_GMM();
extern void Delete_Identical_gdfs_of_GMM();
extern void Copy_GAUSS3D();
extern void Copy_GAUSS3D_Array();
