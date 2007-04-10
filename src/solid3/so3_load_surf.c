/*======================================================================*/
/*!
\file
\brief Integration of tractions on Neumann sturfaces

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


/*----------------------------------------------------------------------*/
/*!
\brief Spatial integration of
       surface tractions on element von Neumann boundaries
\author bborn
\date 04/07
*/
void so3_load_surf_int(ELEMENT* ele,  /*!< current element */
                        SO3_DATA* data,  /*!< Gauss point data */
                        DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],  /*!< deform. element node coord. */
                        INT ngsurf,  /*!< number of element surfaces */
                        GSURF* gsurf[MAXSID_SOLID3],  /*!< geometry surf. */
                        DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3])  /*!< element load vector (is changed) */

{
  /* element */
  const DIS_TYP distyp = ele->distyp;  /* local copy of discretisation type */
  const INT nelenod = ele->numnp;  /* number of element nodes */

  /* surface */
  INT igsurf;  /* element surface index */
  INT idim;  /* dimension index */
  INT idimsid;  /* dimension of surface index */
  DOUBLE sidredm[DIMSID_SOLID3][NDIM_SOLID3];  /* dimens. reduct. matrix */
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
  dstrc_enter("so3_load_surf_int");
#endif

  /*--------------------------------------------------------------------*/
  /* distinguish hex and tet */
  switch (distyp)
  {
    /*==================================================================*/
    /* hexahedron elements */
    case hex8: case hex20: case hex27:
      /*----------------------------------------------------------------*/
      /* loop all element surfaces */
      for (igsurf=0; igsurf<ngsurf; igsurf++)
      {
        /*--------------------------------------------------------------*/
        /* check if current side is subjected to a load */
        if (gsurf[igsurf]->neum != NULL)
        {
          /*------------------------------------------------------------*/
          /* get dimension reduction matrix */
          /* sidredm = &&(data->redsidh[igsurf][0][0]); */
          for (idimsid=0; idimsid<DIMSID_SOLID3; idimsid++)
          {
            for (idim=0; idim<NDIM_SOLID3; idim++)
            {
              sidredm[idimsid][idim] 
                = data->redsidh[igsurf][idimsid][idim];
            }
          }
          /*------------------------------------------------------------*/
          /* get number of Gauss points in each direction */
          idimsid = 0;  /* initialise dimension index */
          for (idim=0; idim<NDIM_SOLID3; idim++)
          {
            /* set number of Gauss points for side integration */
            if ((INT) data->ancsidh[igsurf][idim] == 0)
            {
              gpnum[idimsid] = ele->e.so3->gpnum[idim];
              gpintc[idimsid] = ele->e.so3->gpintc[idim];
              idimsid++;
            }
          }  /* end for */
          /*------------------------------------------------------------*/
          /* integration loops */
          for (igp[0]=0; igp[0]<gpnum[0]; igp[0]++)
          {
            for (igp[1]=0; igp[1]<gpnum[1]; igp[1]++)
            {
              /*--------------------------------------------------------*/
              /* obtain current Gauss coordinates */
              for (idim=0; idim<NDIM_SOLID3; idim++)
              {
                /* set anchor point in current side
                 * This is the intersection of side plane with its
                 * normal parameter coordinate axis */
                gpcidim = data->ancsidh[igsurf][idim];
                for (idimsid=0; idimsid<DIMSID_SOLID3; idimsid++)
                {
                  /* add coordinate components */
                  gpcidim = gpcidim 
                    + sidredm[idimsid][idim] 
                      * data->ghlc[gpintc[idimsid]][igp[idimsid]];
                }  /* end for */
                /* final set of idim-component */
                gpc[idim] = gpcidim;
              }  /* end for */
              /*--------------------------------------------------------*/
              /* Gauss weight */
              fac = 1.0;  /* initialise integration factor */
              for (idimsid=0; idimsid<DIMSID_SOLID3; idimsid++)
              {
                /* multiply weight */
                fac = fac * data->ghlw[gpintc[idimsid]][igp[idimsid]];
              }
              /*--------------------------------------------------------*/
              /* Shape functions at Gauss point */
              so3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1, 
                              shape, deriv);
              /*--------------------------------------------------------*/
              /* Jacobi matrix and determinant */
              so3_metr_surf(ele, nelenod, ex, deriv, sidredm, &metr);
              /*--------------------------------------------------------*/
              /* integration factor */
              fac = fac * metr;
              /*--------------------------------------------------------*/
              /* add surface load contribution ==> eload modified */
              so3_load_surf_val(ele, nelenod, gsurf[igsurf], shape, fac, 
                            eload);
            }  /* end for */
          }  /* end for */
          /*------------------------------------------------------------*/
          /* the surface load of this element has been done,
           * -> switch off the surface load condition */
          gsurf[igsurf]->neum = NULL;
        }  /* end if */
      } /* end for */
      break;  /* end cases */
    /*==================================================================*/
    /* tetrahedron elements */
    case tet4: case tet10:
      /*----------------------------------------------------------------*/
      /* get number of Gauss points for all sides */
      gpnum[0] = ele->e.so3->gpnum[1];
      gpintc[0] = ele->e.so3->gpintc[1];
      /*----------------------------------------------------------------*/
      /* loop all element surfaces */
      for (igsurf=0; igsurf<ngsurf; igsurf++)
      {
        /*--------------------------------------------------------------*/
        /* check if current side is subjected to a load */
        if (gsurf[igsurf]->neum != NULL)
        {
          /*------------------------------------------------------------*/
          /* set dimension reduction matrix */
          for (idimsid=0; idimsid<DIMSID_SOLID3; idimsid++)
          {
            for (idim=0; idim<NDIM_SOLID3; idim++)
            {
              sidredm[idimsid][idim] 
                = data->redsidt[igsurf][idimsid][idim];
            }
          }
          /*------------------------------------------------------------*/
          /* integration loop */
          for (igp[0]=0; igp[0]<gpnum[0]; igp[0]++)
          {
            /*----------------------------------------------------------*/
            /* create Gauss point (r,s,t) coordinate in side */
            for (idim=0; idim<NDIM_SOLID3; idim++)
            {
              /* set anchor point in current side */
              gpcidim = data->ancsidt[igsurf][idim];
              /* add position sideways */
              for (idimsid=0; idimsid<DIMSID_SOLID3; idimsid++)
              {
                gpcidim = gpcidim
                  + sidredm[idim][idimsid]
                    * data->gtsc[gpintc[0]][igp[0]][idimsid];
              }  /* end for */
              /* final set of idim-component */
              gpc[idim] = gpcidim;
            }  /* end for */
            /*----------------------------------------------------------*/
            /* multiply weight */
            fac = data->ghlw[gpintc[0]][igp[0]];
            /*----------------------------------------------------------*/
            /* Shape functions at Gauss point */
            so3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1, 
                            shape, deriv);
            /*----------------------------------------------------------*/
            /* Jacobi matrix and determinant */
            so3_metr_surf(ele, nelenod, ex, deriv, sidredm, &metr);
            /*----------------------------------------------------------*/
            /* integration factor */
            fac = fac * metr;
            /*----------------------------------------------------------*/
            /* add surface load contribution ==> eload modified */
            so3_load_surf_val(ele, nelenod, gsurf[igsurf], shape, fac, 
                          eload);
          }  /* end for */
        }  /* end for */
        /*--------------------------------------------------------------*/
        /* the surface load of this element has been done,
         * -> switch off the surface load condition */
        gsurf[igsurf]->neum = NULL;
      }  /* end if */
      break;  /* end cases */
    /*------------------------------------------------------------------*/
    default:
      dserror("ele->distyp unknown!");
  }  /* end switch */

  /*--------------------------------------------------------------------*/
  /* finish */
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_load_surf_itgr */

/*======================================================================*/
/*!
\brief Determine load due to tractions on element faces (surfaces)

\param   *ele      ELEMENT      (i)    actual element
\param    nelenod  INT          (i)    number of element nodes
\param   *gsurf    GSURF        (i)    current geometry surface
\param    shape[]  DOUBLE       (i)    shape function at Gauss point
\param    fac      DOUBLE       (i)    integration factor at Gauss point
\param    eload[][]DOUBLE       (io)   element load vector contribution
\return void

\author bborn
\date 01/07
*/
void so3_load_surf_val(ELEMENT* ele,
                       INT nelenod,
                       GSURF* gsurf,
                       DOUBLE shape[MAXNOD_SOLID3],
                       DOUBLE fac,
                       DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3])
{
  INT onoff[NUMDOF_SOLID3];
  DOUBLE traction[NUMDOF_SOLID3];  /* elem. bound. traction [force/area] */
  INT idof, inode;  /* loopers */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load_surf_val");
#endif

  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(gsurf->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform traction applied in undeformed configuration */
    case neum_dead:
      for (idof=0; idof<NUMDOF_SOLID3; idof++)
      {
        onoff[idof] = gsurf->neum->neum_onoff.a.iv[idof];
        traction[idof] = gsurf->neum->neum_val.a.dv[idof];
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
} /* end of so3_load_surf_val(...) */

/*======================================================================*/
#endif  /*end of #ifdef D_SOLID3 */
