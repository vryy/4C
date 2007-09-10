/*======================================================================*/
/*!
\file
\brief Print THERM2 element data to Ccarat output file

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 10/06
*/

#ifndef CCADISCRET
/*----------------------------------------------------------------------*/
#ifdef D_THERM2

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm2.h"

/*----------------------------------------------------------------------*/
/*!
\brief General problem data

\author bborn
\date 03/06
*/
extern GENPROB genprob;

/*----------------------------------------------------------------------*/
/*!
\brief Fields

vector of numfld FIELDs, defined in global_control.c

\author bborn
\date 10/06
*/
extern FIELD *field;


/*======================================================================*/
/*!
\brief Write the heatflux in element to Ccarat output file

\param   *actele        ELEMENT     (i)   pointer to current element
\param   *out           FILE        (i/o) pointer to output file

\return void

\author bborn
\date 08/06
*/
void th2_out_hflux(ELEMENT *actele,
                   FILE *out)
{
  const INT place = 0;
  INT i;

  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("th2_out_hflux");
#endif

  /*--------------------------------------------------------------------*/
  /* print title */
  fprintf(out, "_________________________________________________"
               "_______________________________\n");
  fprintf(out, "Element glob_Id %d loc_Id %d"
               "                THERM2\n", actele->Id, actele->Id_loc);
  fprintf(out, "\n");

  /*--------------------------------------------------------------------*/
  /* check heat flux type */
  switch (actele->e.th2->hfluxtype)
  {
    /* none */
    case th2_hflux_none:
      break;
    /* at Gauss points */
    case th2_hflux_gpxy:
    case th2_hflux_gp12:
    case th2_hflux_gpxy12:
      fprintf(out, "Gaussian   "
                   "  heatflux-x   heatflux-y   heatflux-z"
                   "   abs.hflux        angle   heatflux-z\n");
      /* loop all Gauss points */
      for (i=0; i<actele->e.th2->gptotal; i++)
      {
        fprintf(out, "Gauss %2d   "
                     "%12.3E %12.3E %12.3E"
                     "%12.3E %12.3E %12.3E \n",
                i,
                actele->e.th2->hflux_gp.a.d3[place][0][i],
                actele->e.th2->hflux_gp.a.d3[place][1][i],
                actele->e.th2->hflux_gp.a.d3[place][2][i],
                actele->e.th2->hflux_gp.a.d3[place][3][i],
                actele->e.th2->hflux_gp.a.d3[place][4][i],
                actele->e.th2->hflux_gp.a.d3[place][5][i]);
      }
      break;
    /* at element nodes */
    case th2_hflux_ndxy:
    case th2_hflux_nd12:
    case th2_hflux_ndxy12:
      fprintf(out, "Nodal      "
                   "  heatflux-x   heatflux-y   heatflux-z"
                   "   abs.hflux        angle   heatflux-z\n");
      /* loop all element nodes */
      for (i=0; i<actele->numnp; i++)
      {
        fprintf(out, "Node  %2d   "
                     "%12.3E %12.3E %12.3E"
                     "%12.3E %12.3E %12.3E \n",
                i,
                actele->e.th2->hflux_nd.a.d3[place][0][i],
                actele->e.th2->hflux_nd.a.d3[place][1][i],
                actele->e.th2->hflux_nd.a.d3[place][2][i],
                actele->e.th2->hflux_nd.a.d3[place][3][i],
                actele->e.th2->hflux_nd.a.d3[place][4][i],
                actele->e.th2->hflux_nd.a.d3[place][5][i]);
      }
      break;
    default:
      dserror("Heat flux type is not available for printing!");
  }  /* end switch */

  /*====================================================================*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_out_hflux */

#endif /* end of #ifdef D_THERM2 */
/*! @} (documentation module close)*/
#endif
