/*!----------------------------------------------------------------------
\file
\brief service routines for fluid2 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID2_PRO
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2_PRO
#include "../headers/standardtypes.h"
#include "fluid2pro_prototypes.h"
/*!----------------------------------------------------------------------
\brief positions of physical values in node arrays

<pre>                                                        chfoe 11/04

This structure contains the positions of the various fluid solutions 
within the nodal array of sol_increment.a.da[ipos][dim].

extern variable defined in fluid_service.c
</pre>

------------------------------------------------------------------------*/
extern struct _FLUID_POSITION ipos;
/*!---------------------------------------------------------------------
\brief set all arrays for element calculation

<pre>                                                         basol 11/02

get the element velocities and the pressure at different times

</pre>
\param   *elevel      ELEMENT	   (i)    actual element for velocity
\param   *elepre      ELEMENT	   (i)    actual element for pressure
\param   **xyze       DOUBLE	   (o)    coordinates of the element
\param   **eveln      DOUBLE	   (o)    ele velocities at time n
\param   *epren       DOUBLE	   (o)    ele pressures at time n
\return void

------------------------------------------------------------------------*/
void f2pro_calset(
	        ELEMENT         *elevel,
		ELEMENT         *elepre,
		DOUBLE         **xyze,
                DOUBLE         **eveln,
	        DOUBLE          *epren
	      )
{
INT i, veln;         /* simply some counters                            */
NODE *actnode;       /* actual node for element                         */

#ifdef DEBUG
dstrc_enter("f2pro_calset");
#endif

veln = ipos.veln;
/*--------------------------------------------- set element coordinates */
for(i=0;i<elevel->numnp;i++)
{
   xyze[0][i]=elevel->node[i]->x[0];
   xyze[1][i]=elevel->node[i]->x[1];
}
/*---------------------------------------------------------------------*
 | position of the different solutions:                                |
 | node->sol_incement: solution history used for calculations          |
 |       sol_increment[ipos][i]: solution at some time level           |
 |       ipos.velnp .. solution at time n+1                            |
 |       ipos.veln  .. solution at time n                              |
 |       ipos.velnm .. solution at time n-1                            |
 *---------------------------------------------------------------------*/

/* -> computation of time forces -------------------
   -> velocities and pressure at (n) are needed ----*/

      for(i=0;i<elevel->numnp;i++) /* loop nodes of element for velocity  */
      {
         actnode=elevel->node[i];
/*------------------------------------- set element velocities at (n) */
         eveln[0][i]=actnode->sol_increment.a.da[veln][0];
	 eveln[1][i]=actnode->sol_increment.a.da[veln][1];
/*------------------------------------------------- set pressures (n) */
      } /* end of loop over nodes of element for velocity */

/*----------------------------------------------------------------------*
 | REMARK::in projection method that we implement here the pressure     |
 | element and velocity element have different number of nodes and      |
 | different shape functions. Therefore below loop runs up to 4 nodes.  |
 | Where it runs up to 9 in the upper loop.                             |
/-----------------------------------------------------------------------*/
      for(i=0;i<elepre->numnp;i++) /* loop nodes of element for pressure  */
      {
         actnode=elepre->node[i];
	 epren[i]   =actnode->sol_increment.a.da[veln][0];
      } /* end of loop over nodes of element  for pressure */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f2pro_calset */

#endif
/*! @} (documentation module close)*/
