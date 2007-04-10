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
\author bborn
\date 04/07
*/
void so3_load_line_int(ELEMENT* ele,  /*!< current element */
                       SO3_DATA* data,  /*!< Gauss point data */
                       DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],  /*!< deform. element node coord. */
                       INT ngline,  /*!< number of geometry lines */
                       GLINE* gline[MAXEDG_SOLID3],  /*!< element geometry lines (is changed) */
                       DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3])  /*!< element load vector (is changed) */
{
  /* element */
  const DIS_TYP distyp = ele->distyp;  /* local copy of discretisation type */
  const INT nelenod = ele->numnp;  /* number of element nodes */

  /* line */
  INT igline;  /* line index */
  INT idim;  /* dimension index */
  DOUBLE linredv[NDIM_SOLID3];  /* dimension reduction vector */
  DOUBLE metr;  /* metric */
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
          /* get dimension reduction matrix */
          for (idim=0; idim<NDIM_SOLID3; idim++)
          {
            linredv[idim] = data->rededgh[igline][idim];
          }
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
                + linredv[idim] * data->ghlc[gpintc[0]][igp[0]];
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
            /* Jacobi matrix and determinant */
            so3_metr_line(ele, nelenod, ex, deriv, linredv, &metr);
            /*----------------------------------------------------------*/
            /* integration factor */
            fac = fac * metr;
            /*----------------------------------------------------------*/
            /* add surface load contribution ==> eload modified */
            so3_load_line_val(ele, nelenod, gline[igline], shape, fac, 
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
          /* get dimension reduction matrix */
          for (idim=0; idim<NDIM_SOLID3; idim++)
          {
            linredv[idim] = data->rededgt[igline][idim];
          }
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
                + linredv[idim] * data->gtlc[gpintc[0]][igp[0]];
              /* final set of idim-component */
              gpc[idim] = gpcidim;
            }  /* end for */
            /*----------------------------------------------------------*/
            /* set weight */
            fac = data->gtlw[gpintc[0]][igp[0]];
            /*----------------------------------------------------------*/
            /* shape functions at Gauss point */
            so3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1, 
                            shape, deriv);
            /*----------------------------------------------------------*/
            /* Jacobi matrix and determinant */
            so3_metr_line(ele, nelenod, ex, deriv, linredv, &metr);
            /*----------------------------------------------------------*/
            /* integration factor */
            fac = fac * metr;
            /*----------------------------------------------------------*/
            /* add surface load contribution ==> eload modified */
            so3_load_line_val(ele, nelenod, gline[igline], shape, fac, 
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
void so3_load_line_val(ELEMENT* ele,
                       INT nelenod,
                       GLINE* gline,
                       DOUBLE shape[MAXNOD_SOLID3],
                       DOUBLE fac,
                       DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3])
{
  INT onoff[NUMDOF_SOLID3];
  DOUBLE traction[NUMDOF_SOLID3];  /* traction on edge [force/length] */
  INT idof, inode;  /* loopers */

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
            eload[inode][idof] += shape[inode] * traction[idof] * fac;
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
} /* end of so3_load_line_val */

/*======================================================================*/
#endif  /*end of #ifdef D_SOLID3 */
