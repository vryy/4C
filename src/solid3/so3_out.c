/*======================================================================*/
/*!
\file
\brief Output of SOLID3 to classic output *.out file 

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 01/07
*/

/*----------------------------------------------------------------------*/
#ifdef D_SOLID3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Write the stress in element to Ccarat output file

\param   *actele        ELEMENT     (i)   pointer to current element
\param    place         INT         (i)   first index in ARRAY4D stress 
\param   *out           FILE        (i/o) pointer to output file

\return void

\author bborn
\date 01/07
*/
void so3_out_stress(ELEMENT *actele,
                    INT place,
                    FILE *out)
{
  SOLID3 *actso3 = actele->e.so3;  /* pointer to SOLID3 element contents */
  INT ngauss;  /* total number of Gauss points */
  INT igp, inod;  /* loop indices */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_out_stress");
#endif

  /*--------------------------------------------------------------------*/
  /* print header */
  fprintf(out, "_________________________________________________________"
               "_______________________\n");
  fprintf(out, "Element glob_Id %d loc_Id %d                SOLID3\n",
          actele->Id, actele->Id_loc);
  fprintf(out,"\n");
  /*--------------------------------------------------------------------*/
  /* total number of Gauss points ind element domain */
  ngauss = actso3->gptot;
  /*--------------------------------------------------------------------*/
  /* select acc. to output stress type */
  switch(actso3->stresstype)
  {
    /* no stress output */
    case so3_stress_none:
      break;
    /* output stress XYZ-oriented at Gauss points */
    case so3_stress_gpxyz:
      fprintf(out, "INT.point   x-coord.     y-coord.     z-coord."
                   "     stress-xx    stress-yy    stress-zz"
                   "stress-xy    stress-yz    stress-zx\n");
      for (igp=0; igp<ngauss; igp++)
      {
        fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                igp,
                0.0,  /* x-coord */
                0.0,  /* y-coord */
                0.0,  /* z-ccord */
                actso3->stress_gpxyz.a.d3[place][igp][0],  /* stress-xx */
                actso3->stress_gpxyz.a.d3[place][igp][1],  /* stress-yy */
                actso3->stress_gpxyz.a.d3[place][igp][2],  /* stress-zz */
                actso3->stress_gpxyz.a.d3[place][igp][3],  /* stress-xy */
                actso3->stress_gpxyz.a.d3[place][igp][4],  /* stress-yz */
                actso3->stress_gpxyz.a.d3[place][igp][5]);  /* stress-zx */
      }
      break;
    /* output stress rst-oriented at Gauss point */
    case so3_stress_gprst:
      fprintf(out, "r,s,t    ---> local system on element level \n");
      fprintf(out, "rr,ss,tt ---> normal-stresses               \n");
      fprintf(out, "rs,st,tr ---> shear -stresses               \n\n");
      fprintf(out, "INT.point   x-coord.     y-coord.     z-coord."
                   "     stress-rr    stress-ss    stress-tt"
                   "    stress-rs    stress-st    stress-tr\n");
      for (igp=0; igp<ngauss; igp++)
      {
        fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                igp,
                0.0,  /* x-coord */
                0.0,  /* y-coord */
                0.0,  /* z-coord */
                actso3->stress_gprst.a.d3[place][igp][0],  /* stress-rr */
                actso3->stress_gprst.a.d3[place][igp][1],  /* stress-ss */
                actso3->stress_gprst.a.d3[place][igp][2],  /* stress-tt */
                actso3->stress_gprst.a.d3[place][igp][3],  /* stress-rs */
                actso3->stress_gprst.a.d3[place][igp][4],  /* stress-st */
                actso3->stress_gprst.a.d3[place][igp][5]); /* stress-tr */
      }
      break;
    /* output stress principals at Gauss point */
    case so3_stress_gp123:
      fprintf(out, "11,22,33 ---> principal-stresses\n");
      fprintf(out, "r1,s1,t1 ---> angles to the first  principal direction\n");
      fprintf(out, "r2,s2,t2 ---> angles to the second principal direction\n");
      fprintf(out, "r3,s3,t3 ---> angles to the third  principal direction\n\n");
      fprintf(out, "INT.point   stress-11    stress-22    stress-33"
                   "  ang-r1  ang-s1   ang-t1  "
                   "  ang-r2  ang-s2   ang-t2  "
                   "  ang-r3  ang-s3   ang-t3\n");
      for (igp=0; igp<ngauss; igp++)
      {
        fprintf(out,"  %-6d %12.3E %12.3E %12.3E %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f \n",
                igp,
                actso3->stress_gp123.a.d3[place][igp][0],  /* stress-11 */
                actso3->stress_gp123.a.d3[place][igp][1],  /* stress-22 */
                actso3->stress_gp123.a.d3[place][igp][2],  /* stress-33 */
                actso3->stress_gp123.a.d3[place][igp][3],  /* ang-r1 */
                actso3->stress_gp123.a.d3[place][igp][4],  /* ang-s1 */
                actso3->stress_gp123.a.d3[place][igp][5],  /* ang-t1 */
                actso3->stress_gp123.a.d3[place][igp][6],  /* ang-r2 */
                actso3->stress_gp123.a.d3[place][igp][7],  /* ang-s2 */
                actso3->stress_gp123.a.d3[place][igp][8],  /* ang-t2 */
                actso3->stress_gp123.a.d3[place][igp][9],  /* ang-r3 */
                actso3->stress_gp123.a.d3[place][igp][10], /* ang-s3 */
                actso3->stress_gp123.a.d3[place][igp][11]);/* ang-t3 */
      }
      break;
    /* output stress XYZ-oriented at element nodes */
    case so3_stress_ndxyz:
      fprintf(out, "elenode "
                   "    stress-XX    stress-YY    stress-ZZ"
                   "    stress-XY    stress-YZ    stress-ZX\n");
      for (inod=0; inod<actele->numnp; inod++)
      {
        fprintf(out, "  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                inod,
                actso3->stress_ndxyz.a.d3[place][inod][0],  /* stress-XX */
                actso3->stress_ndxyz.a.d3[place][inod][1],  /* stress-YY */
                actso3->stress_ndxyz.a.d3[place][inod][2],  /* stress-ZZ */
                actso3->stress_ndxyz.a.d3[place][inod][3],  /* stress-XY */
                actso3->stress_ndxyz.a.d3[place][inod][4],  /* stress-YZ */
                actso3->stress_ndxyz.a.d3[place][inod][5]); /* stress-ZX */
      }
      break;
    /* output stress rst-oriented at element nodes */
    case so3_stress_ndrst:
      fprintf(out, "elenode "
                   "    stress-rr    stress-ss    stress-tt"
                   "    stress-rs    stress-st    stress-tr\n");
      for (inod=0; inod<actele->numnp; inod++)
      {
        fprintf(out, "  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                inod,
                actso3->stress_ndrst.a.d3[place][inod][0],  /* stress-rr */
                actso3->stress_ndrst.a.d3[place][inod][1],  /* stress-ss */
                actso3->stress_ndrst.a.d3[place][inod][2],  /* stress-tt */
                actso3->stress_ndrst.a.d3[place][inod][3],  /* stress-rs */
                actso3->stress_ndrst.a.d3[place][inod][4],  /* stress-st */
                actso3->stress_ndrst.a.d3[place][inod][5]); /* stress-tr */
      }
      break;
    /* output stress principals at element nodes */
    case so3_stress_nd123:
      fprintf(out, "elenode"
                   "     stress-11    stress-22    stress-33"
                   "  ang-r1  ang-s1  ang-t1  "
                   "  ang-r2  ang-s2  ang-t2  "
                   "  ang-r3  ang-s3  ang-t3\n");
      for (inod=0; inod<actele->numnp; inod++)
      {
        fprintf(out, "  %-6d %12.3E %12.3E %12.3E %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f \n",
                inod,
                actso3->stress_nd123.a.d3[place][inod][0],  /* stress-11 */
                actso3->stress_nd123.a.d3[place][inod][1],  /* stress-22 */
                actso3->stress_nd123.a.d3[place][inod][2],  /* stress-33 */
                actso3->stress_nd123.a.d3[place][inod][3],  /* ang-r1 */
                actso3->stress_nd123.a.d3[place][inod][4],  /* ang-s1 */
                actso3->stress_nd123.a.d3[place][inod][5],  /* ang-t1 */
                actso3->stress_nd123.a.d3[place][inod][6],  /* ang-r2 */
                actso3->stress_nd123.a.d3[place][inod][7],  /* ang-s2 */
                actso3->stress_nd123.a.d3[place][inod][8],  /* ang-t2 */
                actso3->stress_nd123.a.d3[place][inod][9],  /* ang-r3 */
                actso3->stress_nd123.a.d3[place][inod][10],  /* ang-s3 */
                actso3->stress_nd123.a.d3[place][inod][11]);  /* ang-t3 */
      }
      break;
    /* catch it if you must */
    default:
      fprintf(out, "no stresses available\n");
  }  /* end switch */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of so3_out_stress */

#endif /* D_SOLID3 */
/*! @} (documentation module close)*/
