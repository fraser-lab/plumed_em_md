/*

 <TabExpMinX.c>
  Functions for making table for fast calculation of exp(-x).
 
 >> Discretization  <<
 range        :[0..MaxX]
 width of bin : Width
 Ndivide      : MaxX/Width
 lattce_range : [0,1,2.... Ndivide]

 if (x>=tab.MaxX) return(tab->F[tab->Ndivide]);
 if (x<=0.0)      return(tab->F[0]);
 
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TabExpMinX.h" 


/*** FUNCTIONS ***/
void Make_Exp_minus_X_Table();
double Exp_minus_X_Table();

void Make_Exp_minus_X_Table(tab,MaxX,Width)
 struct TAB_EXP_MIN_X *tab;
 double MaxX,Width;
{
 int i;

 tab->MaxX    = MaxX;
 tab->Width   = Width;
 tab->Ndivide = (int)ceil(tab->MaxX/tab->Width);
 printf("#Make_Exp_minus_X_Table(MaxX %lf Width %lf Ndivide %d)\n",tab->MaxX, tab->Width, tab->Ndivide);
 
 /** (1) Malloc Variables **/
 tab->F = (double *)malloc(sizeof(double)*(tab->Ndivide+1));
 tab->X     = (double *)malloc(sizeof(double)*(tab->Ndivide+1));
 tab->Slope  = (double *)malloc(sizeof(double)*(tab->Ndivide+1));
 tab->Cons = (double *)malloc(sizeof(double)*(tab->Ndivide+1));
                                                                                                        
 /** (2) Calulate X[i] and F[i]  **/
 for (i=0;i<=tab->Ndivide;++i){
  tab->X[i] = i * tab->Width;
  tab->F[i] = exp(-tab->X[i]);
 /*
  printf("i %d X %lf exp %lf\n",i,tab->X    [i],tab->F[i]);
 */
 }
                                                                                                        
 /** (3) Calulate Slope[i] and Cons[i]  **/
 for (i=0;i<tab->Ndivide;++i){
  tab->Slope[i] = (tab->F[i+1] - tab->F[i])/tab->Width;
  tab->Cons[i]  = -tab->Slope[i]*tab->X[i] + tab->F[i];
 }
                                                                                                        
 tab->Slope[tab->Ndivide] = tab->Cons[tab->Ndivide] =  0.0;
} /* end of Make_Exp_minus_X_Table() */



double Exp_minus_X_Table(x,tab)
 double x;
 struct TAB_EXP_MIN_X *tab;
{
 int ind;
 
 if (x>=tab->MaxX) return(tab->F[tab->Ndivide]);
 if (x<=0.0)       return(tab->F[0]);

 ind = x/tab->Width; 
 
 if (ind<tab->Ndivide){
  return(x*tab->Slope[ind] + tab->Cons[ind] );
 }
 else {
  return(tab->F[tab->Ndivide]); 
 }
                                                                                                        
} /* end of Exp_minus_X_Table() */

