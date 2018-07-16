/*

 <Voxel.h>

*/

struct VOXEL{
 int   N[3];           /* Number of XYZ column 0:X, 1:Y, 2:Z */
 float grid_width;     /* Grid Length  */
 float OrigPos[3];     /* Real XYZ coordinates on the grid (0,0,0) */
 float min,max,ave,sd; /* 'min'imum, 'max'imum, 'ave'rage, 'sd':standard deviation of voxel values */
 char mal;             /* 1:already malloced 0: not malloced */
 float ***dat;         /* [0..X-1][0..Y-1][0..Z-1] */
};
                                                                                                    

extern int  Malloc_Voxel();
extern void Free_Voxel();
extern void Copy_Voxel();
extern void Initialize_Voxel();
extern void Read_Voxel_File();
extern void Write_Voxel_File();
extern void Write_Voxel_File_In_PDB();
extern void Voxel_Statistics();
extern float Find_Threshold_Density_For_NumVoxel();
extern void Output_Voxel_Histgram();
extern void Make_Reduced_Size_Voxel();
