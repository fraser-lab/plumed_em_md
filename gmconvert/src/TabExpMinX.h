/*
 <TabExpMinX.h>

*/

struct TAB_EXP_MIN_X{
 double *F;          /* F[i] = F(Xi) = exp(-Xi) [Ndivide] (malloc later) */
 double *X;          /* Input Value of exp(-x)  [Ndivide] (malloc later) */
 double *Slope;      /* Slope between i+1 and i := (F[i+1]-F[i])/Wid [Ndivide] (malloc later) */
 double *Cons;       /* Cons[i] := - Slp[i]*x[i] + F[i]  [Ndivide] (malloc later) */
 double  MaxX;        /* Maximum Input X value for exponential  */
 double  Width;      /* Width of Input value for exponential */
 int     Ndivide;     /* number of divition (= MaxInp_iex/WidInp_iex) */

 /*
 >> Basic Principle of interploation  <<
 i = floor(X/Width)
 F(X) - F[i] = [(F[i+1]-F[i])/Width] (X - Xi);
 F(X) = [(F[i+1]-F[i])/Width] (X - Xi) + F[i];
 F(X) = Slope[i](X - Xi) + F[i];
 F(X) = Slope[i]* X - Slope[i]*Xi+F[i];

 Therefore,  F(X) = Slope[i]* X + Cons[i];

 Slope[i] and Cons[i] were also calculated and stored.
  Slope[i] = (F[i+1]-F[i])/Width;
  Cons[i]  = -Slope[i]*Xi + F[i];
 */
};



extern void Make_Exp_minus_X_Table();
extern double Exp_minus_X_Table();
