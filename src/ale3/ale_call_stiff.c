/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale_keku' which calculates the stiffness
matrix at one integration point for a 2d ale element

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale3.h"

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates usual stiffness matrix in total lagrangian formulation 

<pre>                                                              mn 06/02
This routine calculates usual stiffness matrix in total lagrangian
formulation.

</pre>
\param **s   double    (o)  element stiffness matrix 
\param **bs  double    (i)  derivative operator
\param **d   double    (i)  constitutive matrix
\param fac   double    (i)  integration factor
\param nd    int       (i)  total number degrees of freedom of element
\param neps  int       (i)  actual number of strain components

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; caled by: ale2_static_ke(), ale3_static_ke()

*----------------------------------------------------------------------*/
void ale_keku(double  **s, 
              double  **bs, 
              double  **d, 
              double    fac, 
              int       nd,
              int       neps)
{
int            i, j, k, l, m;
double         dum;
double         db[24];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("ale_keku");
#endif
/*----------------------------------------------------------------------*/
   for (j=0; j<nd; j++)
   {
     for (k=0; k<neps; k++)
     {
      db[k] = 0.0 ;                                                              
       for (l=0; l<neps; l++)
       {
       db[k] = db[k] + d[k][l]*bs[l][j]*fac ;
       }
     }
     for (i=0; i<nd; i++)
     {
       dum = 0.0 ;                                                                
       for (m=0; m<neps; m++)
       {
        dum = dum + bs[m][i]*db[m] ;
       }
        s[i][j] = s[i][j] + dum ;
     }
   }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale_keku */
#endif
/*! @} (documentation module close)*/
