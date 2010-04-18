/*-----------------------------------------------------------------------*/
/*!
\file
\brief contains the routine 'f3f_out_gid_sol_str' which
  writes data for a 3D fluid element for gid post processing


<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

#ifndef CCADISCRET
/*!
\addtogroup Fluid3_fast
*//*! @{ (documentation module open)*/


#ifdef D_FLUID3_F

#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"
#include "f3f_prototypes.h"

/*!
  \addtogroup FLUID3_FAST
  *//*! @{ (documentation module open)*/




/*-----------------------------------------------------------------------*/
/*!
  \brief output of data for a 3D fluid element

  This routine writes the viscouse stresses, extrapolated from GP to nodal
  points to the flavia.res file. The stresses are averaged between the
  different elements, so be careful with the interpretation! Test future
  versions of "Gid" for integration point value visualization!

  \param out      *FILE  (i)   file to be written on (flavia.res)
  \param actfield *FIELD (i)   actual field
  \param init      INT   (i)   allocate/free memory

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3f_out_gid_sol_str(
    FILE               *out,           /* File pointer to flavia.res */
    FIELD              *actfield,      /* active field */
    INT                 disnum,
    INT                 init           /* allocate/free memory */
    )
{

  INT                 i,j,k;         /* counter */
  INT                 numnp;         /* number of nodal points in actfield */
  INT                numele;
  INT                 count;
  DOUBLE           invcount;
  DOUBLE           **stress;
  DOUBLE             str[6];
  ELEMENT           *actele;
  NODE            *sactnode;
  NODE            *mactnode;


#ifdef DEBUG
  dstrc_enter("f3f_out_gid_sol_str");
#endif


  numnp  = actfield->dis[disnum].numnp;
  /* write stresses */
  for (i=0; i<numnp; i++)
  {
    mactnode = &(actfield->dis[disnum].node[i]);
    count = 0;


#ifdef SUBDIV
    if (actfield->subdivide > 0)
      sactnode = mactnode->slave_node;
    else
#endif
      sactnode = mactnode;


    for (j=0; j<6; j++) str[j] = 0.0;
    numele  = sactnode->numele;


    for (j=0; j<numele; j++)
    {
      actele = sactnode->element[j];
      if (actele->eltyp != el_fluid3_fast || actele->numnp !=8) continue;
      count++;
      stress=actele->e.f3->stress_ND.a.da;
      for(k=0; k<8; k++)
        if (actele->node[k] == sactnode) break;
#ifndef PARALLEL
      /* Stress output does not work in parallel execution! */
      str[0] += stress[k][0];
      str[1] += stress[k][1];
      str[2] += stress[k][2];
      str[3] += stress[k][3];
      str[4] += stress[k][4];
      str[5] += stress[k][5];
#endif
    }

    invcount = 1.0/count;
    for (j=0; j<6; j++) str[j] = str[j] * invcount;

    fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
        mactnode->Id+1,
        str[0],
        str[1],
        str[2],
        str[3],
        str[4],
        str[5]
        );
  }


#ifdef DEBUG
    dstrc_exit();
#endif


    return;
} /* end of f3f_out_gid_sol_str */


#endif /* ifdef D_FLUID3_F */

/*! @} (documentation module close)*/


#endif
