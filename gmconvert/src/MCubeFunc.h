/*
 <MCubeFunc.h>
*/

extern void Marching_Cube_Tetrahedral();
extern void Free_MC_Vertex_Stack();
extern void Free_MC_Face_Stack();
extern struct MC_VERTEX *Add_MC_Vertex_Stack();
extern struct MC_FACE *Add_MC_Face_Stack();
extern struct MC_FACE *Add_MC_Face_Stack_Simple();
extern void Malloc_MIDMAP();


/** these functions are external mainly for 'MCubeFuncHex.c' **/
extern void Sub_Vec_F();
extern void Sub_Vec_FI();
extern void Cross_Prod_F();
extern float Dot_Prod_F();
extern void Make_Mid_Point();
extern struct MC_VERTEX *Add_Pnt_to_Ver_Map();
extern struct MC_VERTEX *Find_Ver_Map();
extern void Malloc_CMATRIX();
extern void Free_CMATRIX();
extern void Malloc_MIDMAP_YZ();
extern void Free_MIDMAP_YZ();
extern void Reset_MIDMAP_YZ();
extern void Free_MC_Edge_Stack();
extern struct MC_EDGE *Add_MC_Edge_Stack();
extern int  check_MC_Vertex_neighbors();
extern int Check_Overlap_MC_Edge();
extern int common_xyz_in_two_vertexes();
extern int Number_Of_MC_VERTEX();
extern int Number_Of_MC_FACE();
