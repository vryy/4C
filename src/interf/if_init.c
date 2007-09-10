/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ifinit' which initializes the element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_INTERF
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h"

/*!
\addtogroup INTERF
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  initialization routine for the interface element

<pre>                                                              ah 05/03
This routine inits the integration points and allocates stresses.

</pre>
\param *actpart      PARTITION   (I)   my partition

\warning There is nothing special to this routine.
\return void

*----------------------------------------------------------------------*/
void ifinit(PARTITION *actpart)
{
INT             i,k,size_i,size_j;
ELEMENT     *actele;

#ifdef DEBUG
dstrc_enter("ifinit");
#endif

for (i=0; i<actpart->pdis[0].numele; i++)
{
  actele = actpart->pdis[0].element[i];
  if (actele->eltyp != el_interf) continue;
  /*---------------------------------------- init integration points ---*/
  /*ifintg(&data);*/
  /*-------------------------------- allocate the space for stresses ---*/
 am4def("stress_GP",&(actele->e.interf->stress_GP),1,3,4,0,"D3");
 am4def("stress_ND",&(actele->e.interf->stress_ND),1,3,8,0,"D3");
 size_i = 1;
 actele->e.interf->elewa = (IF_ELE_WA*)CCACALLOC(size_i,sizeof(IF_ELE_WA));
 size_j = actele->e.interf->nGP;
 actele->e.interf->elewa[0].ipwa = (IF_IP_WA*)CCACALLOC(size_j,sizeof(IF_IP_WA));
 for (k=0; k<size_j; k++)
 {
    actele->e.interf->elewa[0].ipwa[k].yip = -1;
    actele->e.interf->elewa[0].ipwa[k].jump_ut_pl = 0.0;
    actele->e.interf->elewa[0].ipwa[k].Q[0][0] = 0.0;
    actele->e.interf->elewa[0].ipwa[k].Q[0][1] = 0.0;
    actele->e.interf->elewa[0].ipwa[k].Q[1][0] = 0.0;
    actele->e.interf->elewa[0].ipwa[k].Q[1][1] = 0.0;
    actele->e.interf->elewa[0].ipwa[k].dn = 0.0;
    actele->e.interf->elewa[0].ipwa[k].dt = 0.0;
    actele->e.interf->elewa[0].ipwa[k].Tn = 0.0;
    actele->e.interf->elewa[0].ipwa[k].Tt = 0.0;
    actele->e.interf->elewa[0].ipwa[k].kappa_n = 0.0;
    actele->e.interf->elewa[0].ipwa[k].kappa_t = 0.0;
 }

}
return;
} /* end of ifinit */
/*----------------------------------------------------------------------*/
#endif /*D_INTERF*/
/*! @} (documentation module close)*/
#endif
