/*

 <Matrix3D.h>
 
*/

/** FUNCTIONS (GLOBAL) **/
extern double Dot_Prod3D();
extern void Cross_Prod3D();
extern void Add_Vec3D();
extern void Sub_Vec3D();
extern void Equal_Vec3D();
extern void Multiply_Scalar_Vec3D();
extern double Norm_Vec3D();
extern double Length_Vec3D();
extern double Normalize_Vec3D();
extern double Length_Normalize_Vec3D();
extern void Print_Vec3D();
extern int Cal_Inverse_Matrix3D_by_Cramer_Rule();
extern int Cal_Inverse_Matrix3D_by_Cramer_Rule_Symmetric();
extern void Set_Unit_Matrix3D();
extern void Equal_Matrix3D();
extern void Add_Matrix3D();
extern void Sub_Matrix3D();
extern void Multiply_Scalar_Matrix3D();
extern void Multiply_Scalar_Matrix3D_Self();
extern void Multiply_Scalar_Matrix3D_Symmetric();
extern void Trace_Matrix3D();
extern void Multiply_Matrix3D();
extern void Multiply_Matrix3D_Vec3D();
extern void Multiply_Transpose_Matrix3D();
extern void Transform_RAtR_Matrix3D();
extern void Transform_RAtQ_Matrix3D();
extern void Transform_Rx_plus_tvec_3D();
extern void Transform_Rx_minus_x0_plus_tvec_3D();
extern void Transform_Rmat_around_tvec_3D();
extern void Multiply_Matrix_By_Vector3D();
extern void Print_Matrix3D();
extern void Transpose_Matrix3D();
extern void Make_Matrix3D_By_Prod_Vector3D();
extern double Determinant_Matrix3D();
extern double Quadratic_Form_3D();
extern void Cal_EigenVectors_For_Matrix3D();
