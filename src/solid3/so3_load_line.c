/*======================================================================*/
/*!
\file
\brief Spatial integration of loads (ie body forces/traction) applied
       to edges (lines)

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/
#ifdef D_SOLID3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*======================================================================*/
/*!
\brief Spatial integration of 
       traction on element edges (lines) [force/length]

\param  ele       ELEMENT*      (i)  current element
\param  data      SO3_DATA*     (i)  Gauss point data
\param  ex        DOUBLE[]      (i)  material element node coord.
\param  exs       DOUBLE[]      (i)  spatial element node coord.s
\param  ngline    INT           (i)  number of geometry lines
\param  gline     GLINE*[]      (i)   element geometry lines (is changed)
\param  eload     DOUBLE[][]    (io) element load vector (is changed)
\author bborn
\date 04/07
*/
void so3_load_line_int(ELEMENT* ele,
                       SO3_DATA* data,
                       DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                       DOUBLE exs[MAXNOD_SOLID3][NDIM_SOLID3],
                       INT ngline,
                       GLINE* gline[MAXEDG_SOLID3],
                       DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3])
{
  /* element */
  const DIS_TYP distyp = ele->distyp;  /* local copy of discretisation type */

  /* line */
  INT igline;  /* line index */
  INT idim;  /* dimension index */
  DOUBLE gpcidim;  /* dummy Gauss point coordinate */

  /* integration */
  INT igp[NDIM_SOLID3];  /* Gauss point indices */
  INT gpnum[NDIM_SOLID3];  /* Gauss point numbers */
  INT gpintc[NDIM_SOLID3];  /* Gauss point integration cases */
  DOUBLE fac;  /* integration factors */
  DOUBLE gpc[NDIM_SOLID3];  /* r,s,t-coord current GP */
  DOUBLE shape[MAXNOD_SOLID3];  /* shape functions */
  DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3];  /* derivatives of shape fct. */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load_line_int");
#endif

  /*--------------------------------------------------------------------*/
  /* distinguish hex and tet */
  switch (distyp)
  {
    /*==================================================================*/
    /* hexahedron elements */
    case hex8: case hex20: case hex27:
      /*----------------------------------------------------------------*/
      /* loop all element lines */
      for (igline=0; igline<ngline; igline++)
      {
        /*--------------------------------------------------------------*/
        /* check if current line is subjected to a load */
        if (gline[igline]->neum != NULL)
        {
          /*------------------------------------------------------------*/
          /* get number of Gauss points for current line
           * these could change depending on the set Gauss point numbers
           * in r-, s- and t-direction */
          for (idim=0; idim<NDIM_SOLID3; idim++)
          {
            /* set number of Gauss points for side integration
             * the anchor ancsidh is normal on the edge direction */
            if ((INT) data->ancedgh[igline][idim] == 0)
            {
              gpnum[0] = ele->e.so3->gpnum[idim];
              gpintc[0] = ele->e.so3->gpintc[idim];
            }
          }  /* end for */
          /*------------------------------------------------------------*/
          /* integration loops 
           * For each side one of the following loops is not repeated,
           * but as all sides are generally integrated here. */
          for (igp[0]=0; igp[0]<gpnum[0]; igp[0]++)
          {
            /*----------------------------------------------------------*/
            /* obtain current Gauss coordinates and weights */
            for (idim=0; idim<NDIM_SOLID3; idim++)
            {
              /* set anchor point in current side
               * This is the intersection of edge with its normal 
               * parameter coordinate plane */
              gpcidim = data->ancedgh[igline][idim];
              /* add coordinate components */
              gpcidim = gpcidim
                + data->rededgh[igline][idim] * data->ghlc[gpintc[0]][igp[0]];
              /* final set of idim-component */
              gpc[idim] = gpcidim;
            }  /* end for */
            /*----------------------------------------------------------*/
            /* set weight */
            fac = data->ghlw[gpintc[0]][igp[0]];
            /*----------------------------------------------------------*/
            /* Shape functions at Gauss point */
            so3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1, 
                            shape, deriv);
            /*----------------------------------------------------------*/
            /* add surface load contribution ==> eload modified */
            so3_load_line_valh(ele, data, igline, gline[igline], 
                               ex, exs, gpc, shape, deriv, fac, 
                               eload);
          }  /* end for */
          /*------------------------------------------------------------*/
          /* the line load of this element has been done,
           * -> switch off the line load condition */
          gline[igline]->neum = NULL;
        }  /* end if */
      }  /* end for */
      break;  /* end of cases */
    /*==================================================================*/
    /* tetrahedron elements */
    case tet4: case tet10:
      /*----------------------------------------------------------------*/
      /* get number of Gauss points for all sides */
      gpnum[0] = ele->e.so3->gpnum[2];
      gpintc[0] = ele->e.so3->gpintc[2];
      /*----------------------------------------------------------------*/
      /* loop all element lines */
      for (igline=0; igline<ngline; igline++)
      {
        /*--------------------------------------------------------------*/
        /* check if current line is subjected to a load */
        if (gline[igline]->neum != NULL)
        {
          /*------------------------------------------------------------*/
          /* integration loops 
           * For each side one of the following loops is not repeated,
           * but as all sides are generally integrated here. */
          for (igp[0]=0; igp[0]<gpnum[0]; igp[0]++)
          {
            /*----------------------------------------------------------*/
            /* obtain current Gauss coordinates and weights */
            for (idim=0; idim<NDIM_SOLID3; idim++)
            {
              /* set anchor point in current side
               * This is the intersection of edge with its normal 
               * parameter coordinate plane */
              gpcidim = data->ancedgt[igline][idim];
              /* add coordinate components */
              gpcidim = gpcidim
                + data->rededgt[igline][idim] * data->gtlc[gpintc[0]][igp[0]];
              /* final set of idim-component */
              gpc[idim] = gpcidim;
            }  /* end for */
            /*----------------------------------------------------------*/
            /* shape functions at Gauss point */
            so3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1, 
                            shape, deriv);
            /*----------------------------------------------------------*/
            /* add surface load contribution ==> eload modified */
            so3_load_line_valt(ele, data, igline, gline[igline], 
                               ex, exs, gpc, shape, deriv, fac, 
                               eload);
          }  /* end for */
          /*------------------------------------------------------------*/
          /* the line load of this element has been done,
           * -> switch off the line load condition */
          gline[igline]->neum = NULL;
        }  /* end if */
      }  /* end for */
      break;  /* end of cases */
    default:
      dserror("ele->distyp unknown!");
  }  /* end switch */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_load_line_int(...) */


/*======================================================================*/
/*!
\brief Determine load due to tractions on element edges (lines)

\param   *ele      ELEMENT      (i)    actual element
\param    nelenod  INT          (i)    number of element nodes
\param   *igline   GLINE        (i)    current geometry line
\param    shape[]  DOUBLE       (i)    shape function at Gauss point
\param    fac      DOUBLE       (i)    integration factor
\param    eload[][]DOUBLE       (io)   element load vector contribution
\return void

\author bborn
\date 01/07
*/
void so3_load_line_valh(ELEMENT* ele,
                        SO3_DATA* data,
                        INT igline,
                        GLINE* gline,
                        DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE exs[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE gpc[NDIM_SOLID3],
                        DOUBLE shape[MAXNOD_SOLID3],
                        DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE fac,
                        DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3])
{
  const INT nelenod = ele->numnp;  /* number of element nodes */
  const DIS_TYP distyp = ele->distyp;  /* local copy of discretisation type */
  INT onoff[NUMDOF_SOLID3];
  DOUBLE traction[NUMDOF_SOLID3];  /* traction on edge [force/length] */
  INT idof, jdof, inode;  /* loopers */
  DOUBLE metr;  /* metric */
  DOUBLE gamt[NDIM_SOLID3];  /* base vector physically oriented */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load_line_val");
#endif

  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(gline->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform prescribed surface load */
    case neum_dead:
      /* metric at Gauss point */
      so3_metr_line(ele, nelenod, ex, deriv, data->rededgh[igline], 
                    NULL, &metr);
      for (idof=0; idof<NUMDOF_SOLID3; idof++)
      {
        onoff[idof] = gline->neum->neum_onoff.a.iv[idof];
        traction[idof] = gline->neum->neum_val.a.dv[idof];
      }
      /* add load vector component to element load vector */
      for (inode=0; inode<nelenod; inode++)
      {
        for (idof=0; idof<NUMDOF_SOLID3; idof++)
        {
          /* if load is switched on : apply */
          if (onoff[idof] == 1)
          {
            eload[inode][idof] += shape[inode] * traction[idof] * fac * metr;
          }
        }
      }
      break;
    /*------------------------------------------------------------------*/
    /* uniform tangential line load */
#if 0
    /* quick hack to apply tangential 'line load' on a cylindrical 
     * cantilever beam (length 120) subjected to a torsional torque
     * at its tip. (bborn/mgit 04/07) */
    case neum_orthopressure:
      idof = 0;  /* only 1st index relevant */
      onoff[idof] = gline->neum->neum_onoff.a.iv[idof];
      if (onoff[idof] == 1)
      {
        /* traction value */
        traction[idof] = gline->neum->neum_val.a.dv[idof];
        /* metric at Gauss point */
        so3_metr_line(ele, nelenod, ex, deriv, data->rededgh[igline], 
                      gamt, &metr);
        /* make unit base vector */
        so3_tns3_unitvct(gamt);
        /* check orientation of line */
        INT node0 = data->nodedghl[igline][0];
        INT node1 = data->nodedghl[igline][1];
        DOUBLE cc[NDIM_SOLID3] = {120.0, 0.0, 0.0};
        DOUBLE c0[NDIM_SOLID3];
        DOUBLE c1[NDIM_SOLID3];
        DOUBLE c2[NDIM_SOLID3];
        INT idim;
        for (idim=0; idim<NDIM_SOLID3; idim++)
        {
          c0[idim] = ex[node0][idim] - cc[idim];
          c1[idim] = ex[node1][idim] - cc[idim];
        }
        so3_tns3_unrm(c0, c1, 0, c2);
        DOUBLE dir;
        if (c2[0] > 0.0)
        {
          dir = +1.0;
        }
        else
        {
          dir = -1.0;
        }
        /* */
        for (inode=0; inode<nelenod; inode++)
        {
          for (jdof=0; jdof<NUMDOF_SOLID3; jdof++)
          {
            eload[inode][jdof] += shape[inode] * traction[idof] 
              * dir * gamt[jdof] * fac * metr;
          }
        }
      }
      break;
#endif
    /*------------------------------------------------------------------*/
    default:
      dserror("load case unknown");
      break;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_load_line_valh */

/*======================================================================*/
/*!
\brief Determine load due to tractions on element edges (lines)

\param   *ele      ELEMENT      (i)    actual element
\param    nelenod  INT          (i)    number of element nodes
\param   *igline   GLINE        (i)    current geometry line
\param    shape[]  DOUBLE       (i)    shape function at Gauss point
\param    fac      DOUBLE       (i)    integration factor
\param    eload[][]DOUBLE       (io)   element load vector contribution
\return void

\author bborn
\date 01/07
*/
void so3_load_line_valt(ELEMENT* ele,
                        SO3_DATA* data,
                        INT igline,
                        GLINE* gline,
                        DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE exs[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE gpc[NDIM_SOLID3],
                        DOUBLE shape[MAXNOD_SOLID3],
                        DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE fac,
                        DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3])
{
  const INT nelenod = ele->numnp;  /* number of element nodes */
  INT onoff[NUMDOF_SOLID3];
  DOUBLE traction[NUMDOF_SOLID3];  /* traction on edge [force/length] */
  INT idof, inode;  /* loopers */
  DOUBLE metr;  /* metric */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load_line_valt");
#endif

  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(gline->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform prescribed surface load */
    case neum_dead:
      /* metric at Gauss point */
      so3_metr_line(ele, nelenod, ex, deriv, data->rededgt[igline], 
                    NULL, &metr);
      for (idof=0; idof<NUMDOF_SOLID3; idof++)
      {
        onoff[idof] = gline->neum->neum_onoff.a.iv[idof];
        traction[idof] = gline->neum->neum_val.a.dv[idof];
      }
      /* add load vector component to element load vector */
      for (inode=0; inode<nelenod; inode++)
      {
        for (idof=0; idof<NUMDOF_SOLID3; idof++)
        {
          /* if load is switched on : apply */
          if (onoff[idof] == 1)
          {
            eload[inode][idof] += shape[inode] * traction[idof] * fac * metr;
          }
        }
      }
      break;
    /*------------------------------------------------------------------*/
    default:
      dserror("load case unknown");
      break;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_load_line_valt */

/*======================================================================*/
#endif  /*end of #ifdef D_SOLID3 */
