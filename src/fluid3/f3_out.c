
/*!----------------------------------------------------------------------
  \file
  \brief contains the routine 'f3_out_gid_sol_str' which
  writes data for a 3D fluid element for gid post processing

  <pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0771 - 685-6121
  </pre>

 *----------------------------------------------------------------------*/
#ifdef D_FLUID3

#include "../headers/standardtypes.h"
#include "fluid3.h"
#include "fluid3_prototypes.h"

/*!
  \addtogroup FLUID3
  *//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
  \brief output of data for a 3D fluid element

  <pre>                                                             mn 03/04
  This routine writes the viscouse stresses, extrapolated from GP to nodal
  points to the flavia.res file. The stresses are averaged between the
  different elements, so be careful with the interpretation! Test future
  versions of "Gid" for integration point value visualization!
  </pre>

  \param  out       FILE*   (i)   file to be written on (flavia.res)
  \param  actfield FIELD*   (i)   actual field
  \param  init       INT*   (i)   allocate/free memory

  \warning Check if it is possible to give "Gid" integration point values.
  If it is possible change this code!

  \return void
  \sa calling: ---; called by: out_gid_sol()

 *----------------------------------------------------------------------*/

void f3_out_gid_sol_str(
    FILE       *out,     /* File pointer to flavia.res */
    FIELD *actfield,     /* active field */
    INT        type,     /* viscous (1) or press (0) stresses */
    INT         init     /* allocate/free memory */
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
  NODE             *actnode;


#ifdef DEBUG
  dstrc_enter("f3_out_gid_sol_str");
#endif


  numnp  = actfield->dis[0].numnp;
  /* write stresses */
  for (i=0; i<numnp; i++)
  {
    actnode = &(actfield->dis[0].node[i]);
    count = 0;
    for (j=0; j<6; j++) str[j] = 0.0;
    numele  = actnode->numele;

    for (j=0; j<numele; j++)
    {
      actele = actnode->element[j];
      if (actele->eltyp != el_fluid3 || actele->numnp !=8) continue;
      count++;
      stress=actele->e.f3->stress_ND.a.da;
      for(k=0; k<8; k++)
        if (actele->node[k] == actnode) break;
      str[0] += stress[k][type*6+0];
      str[1] += stress[k][type*6+1];
      str[2] += stress[k][type*6+2];
      str[3] += stress[k][type*6+3];
      str[4] += stress[k][type*6+4];
      str[5] += stress[k][type*6+5];
    }

    invcount = 1.0/count;
    for (j=0; j<6; j++) str[j] = str[j] * invcount;

    fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
        actnode->Id+1,
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
} /* end of c1_out_gid_sol_str */
/*----------------------------------------------------------------------*/
#endif /*D_FLUID3*/
/*! @} (documentation module close)*/
