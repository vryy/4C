/*======================================================================*/
/*!
\file
\brief Print THERM3 element data to Ccarat output file

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 10/06
*/

/*----------------------------------------------------------------------*/
#ifdef D_THERM3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm3.h"

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
void th3_out_hflux(ELEMENT *actele,
                   FILE *out)
{
  const INT place = 0;
  INT gptotal = 0;  /* total number of Gauss points in element domain */
  INT nelenod;  /* number of element nodes */
  INT igp;  /* Gauss point counter */
  INT inod;  /* node counter */
  
  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("th3_out_hflux");
#endif

  /*--------------------------------------------------------------------*/
  /* print title */
  fprintf(out, "_________________________________________________"
               "_______________________________\n");
  fprintf(out, "Element glob_Id %d loc_Id %d"
               "                THERM3\n", actele->Id, actele->Id_loc);
  fprintf(out, "\n");

  /*--------------------------------------------------------------------*/
  /* get total number of Gauss points */
  if ( (actele->distyp == hex8) 
       || (actele->distyp == hex20)
       || (actele->distyp == hex27) )
  {
    gptotal = actele->e.th3->gpnum[0]
      * actele->e.th3->gpnum[1]
      * actele->e.th3->gpnum[2];
  }
  else if ( (actele->distyp == tet4) 
            || (actele->distyp == tet10) )
  {
    gptotal = actele->e.th3->gpnum[0];
  }

  /*--------------------------------------------------------------------*/
  /* number of element nodes */
  nelenod = actele->numnp;

  /*--------------------------------------------------------------------*/
  /* check heat flux type */
  switch (actele->e.th3->hfluxtype)
  {
    /*------------------------------------------------------------------*/
    /* none */
    case th3_hflux_none:
      break;
    /*------------------------------------------------------------------*/
    /* at Gauss points */
    case th3_hflux_gpxyz:
      fprintf(out, "Gaussian   "
                   "  heatflux-x   heatflux-y   heatflux-z\n");
      /* loop all Gauss points */
      for (igp=0; igp<gptotal; igp++)
      {
        fprintf(out, "      %2d"
                     "   %12.3E %12.3E %12.3E\n",
                igp,
                actele->e.th3->hflux_gp_xyz.a.d3[place][igp][0],
                actele->e.th3->hflux_gp_xyz.a.d3[place][igp][1],
                actele->e.th3->hflux_gp_xyz.a.d3[place][igp][2]);
      }
      break;
    case th3_hflux_gp123:
      fprintf(out, "Gaussian   "
                   "  abs.hflux"
                   "  x-angle      y-angle      z-angle\n");
      /* loop all Gauss points */
      for (igp=0; igp<gptotal; igp++)
      {
        fprintf(out, "      %2d   "
                     "%12.3E %12.3E %12.3E "
                     "%12.3E \n",
                igp,
                actele->e.th3->hflux_gp_123.a.d3[place][igp][3],  /* mod */
                actele->e.th3->hflux_gp_123.a.d3[place][igp][0],
                actele->e.th3->hflux_gp_123.a.d3[place][igp][1],
                actele->e.th3->hflux_gp_123.a.d3[place][igp][2]);
      }
      break;
    case th3_hflux_gpxyz123:
      fprintf(out, "Gaussian   "
                   "  heatflux-x   heatflux-y   heatflux-z  "
                   "  abs.hflux   "
                   "   x-angle      y-angle      z-angle  \n");
      /* loop all Gauss points */
      for (igp=0; igp<gptotal; igp++)
      {
        fprintf(out, "      %2d   "
                     "%12.3E %12.3E %12.3E "
                     "%12.3E %12.3E %12.3E "
                     "%12.3E \n",
                igp,
                actele->e.th3->hflux_gp_xyz.a.d3[place][igp][0],
                actele->e.th3->hflux_gp_xyz.a.d3[place][igp][1],
                actele->e.th3->hflux_gp_xyz.a.d3[place][igp][2],
                actele->e.th3->hflux_gp_123.a.d3[place][igp][3],  /* mod */
                actele->e.th3->hflux_gp_123.a.d3[place][igp][0],
                actele->e.th3->hflux_gp_123.a.d3[place][igp][1],
                actele->e.th3->hflux_gp_123.a.d3[place][igp][2]);
      }
      break;
      break;
    case th3_hflux_gprst:
      fprintf(out, "Gaussian   "
                   "  heatflux-r   heatflux-s   heatflux-t\n");
      /* loop all Gauss points */
      for (igp=0; igp<gptotal; igp++)
      {
        fprintf(out, "      %2d   "
                     "%12.3E %12.3E %12.3E",
                igp,
                actele->e.th3->hflux_gp_rst.a.d3[place][igp][0],
                actele->e.th3->hflux_gp_rst.a.d3[place][igp][1],
                actele->e.th3->hflux_gp_rst.a.d3[place][igp][2]);
      }
      break;
    /*------------------------------------------------------------------*/
    /* at element nodes */
    case th3_hflux_ndxyz:
      fprintf(out, "Node       "
                   "  heatflux-x   heatflux-y   heatflux-z\n");
      /* loop all Gauss points */
      for (inod=0; inod<nelenod; inod++)
      {
        fprintf(out, "      %2d   "
                     "%12.3E %12.3E %12.3E\n",
                inod,
                actele->e.th3->hflux_nd_xyz.a.d3[place][inod][0],
                actele->e.th3->hflux_nd_xyz.a.d3[place][inod][1],
                actele->e.th3->hflux_nd_xyz.a.d3[place][inod][2]);
      }
      break;
    case th3_hflux_nd123:
      fprintf(out, "Node       "
                   "  abs.hflux"
                   "  x-angle      y-angle      z-angle\n");
      /* loop all Gauss points */
      for (inod=0; inod<nelenod; inod++)
      {
        fprintf(out, "Node  %2d   "
                     "%12.3E %12.3E %12.3E "
                     "%12.3E \n",
                inod,
                actele->e.th3->hflux_nd_123.a.d3[place][inod][3],  /* mod */
                actele->e.th3->hflux_nd_123.a.d3[place][inod][0],
                actele->e.th3->hflux_nd_123.a.d3[place][inod][1],
                actele->e.th3->hflux_nd_123.a.d3[place][inod][2]);
      }
      break;
    case th3_hflux_ndxyz123:
      fprintf(out, "Node       "
                   "  heatflux-x   heatflux-y   heatflux-z  "
                   "  abs.hflux   "
                   "   x-angle      y-angle      z-angle  \n");
      /* loop all Gauss points */
      for (inod=0; inod<nelenod; inod++)
      {
        fprintf(out, "      %2d   "
                     "%12.3E %12.3E %12.3E"
                     "%12.3E %12.3E %12.3E "
                     "%12.3E \n",
                inod,
                actele->e.th3->hflux_nd_xyz.a.d3[place][inod][0],
                actele->e.th3->hflux_nd_xyz.a.d3[place][inod][1],
                actele->e.th3->hflux_nd_xyz.a.d3[place][inod][2],
                actele->e.th3->hflux_nd_123.a.d3[place][inod][3],  /* mod */
                actele->e.th3->hflux_nd_123.a.d3[place][inod][0],
                actele->e.th3->hflux_nd_123.a.d3[place][inod][1],
                actele->e.th3->hflux_nd_123.a.d3[place][inod][2]);
      }
      break;
    case th3_hflux_ndrst:
      dserror("Parametric heatflux in (r,s,t) co-ordinates "
              "is not implemented!\n");
      fprintf(out, "Node       "
                   "  heatflux-r   heatflux-s   heatflux-t"
                   "  x-angle      y-angle      z-angle   "
                   "  abs.hflux\n");
      /* loop all Gauss points */
      for (inod=0; inod<nelenod; inod++)
      {
        fprintf(out, "      %2d   "
                     "%12.3E %12.3E %12.3E",
                inod,
                actele->e.th3->hflux_nd_rst.a.d3[place][inod][0],
                actele->e.th3->hflux_nd_rst.a.d3[place][inod][1],
                actele->e.th3->hflux_nd_rst.a.d3[place][inod][2]);
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
}  /* end of th3_out_hflux */

#endif /* end of #ifdef D_THERM3 */
/*! @} (documentation module close)*/
