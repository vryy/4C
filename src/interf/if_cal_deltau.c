/*!-----------------------------------------------------------------------
\file
\brief contains the routine 'if_jumpu' which calculates [un],[ut],
       DELTA[un] andDELTA[ut]
*-----------------------------------------------------------------------*/
#ifdef D_INTERF
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h"

/*! 
\addtogroup INTERF
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  contains the dispacement-jump calculation for the interface element

<pre>                                                              mn 05/03 
This routine computes the dispacement-jump of the interface element

</pre>
\param *mat          STVENANT   (i)   blabal

\warning There is nothing special to this routine
\return void                                               
\sa calling:   ---; 
    called by: if_mat();

*----------------------------------------------------------------------*/
void if_jumpu(ELEMENT  *ele, 
              DOUBLE  **bop,
              DOUBLE   *disjump,
              DOUBLE   *Deltadisjump) 
{
INT i,j,nodei,iel;
DOUBLE d[16];
DOUBLE Deltad[16];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("if_deltau");
#endif
/* get nodal displ. and nodal increm. displ (=Delta d + sum delta d) ---*/

iel     = ele->numnp;

for (nodei=0; nodei<iel; nodei++)
{
  d[2*nodei]        = ele->node[nodei]->sol.a.da[0][0];
  d[2*nodei+1]      = ele->node[nodei]->sol.a.da[0][1];
  Deltad[2*nodei]   = ele->node[nodei]->sol_increment.a.da[0][0];
  Deltad[2*nodei+1] = ele->node[nodei]->sol_increment.a.da[0][1];
}
/*-- calculate tangential and normal displacement jumps at interface ---*/
for (i=0; i<2; i++)
{
  Deltadisjump[i] = 0.0;
  disjump[i]      = 0.0;
  for (j=0; j<2*iel; j++)
  {
    disjump[i]      += bop[i][j] * d[j];
    Deltadisjump[i] += bop[i][j] * Deltad[j];
  }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of if_jumpu */

/*----------------------------------------------------------------------*/
#endif /*D_INTERF*/
/*! @} (documentation module close)*/
