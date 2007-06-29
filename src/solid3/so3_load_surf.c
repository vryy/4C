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
\brief General problem data
\author bborn
\date 05/07
*/
extern GENPROB genprob;

/*----------------------------------------------------------------------*/
/*!
\brief Spatial integration of
       surface tractions on element von Neumann boundaries

\param  ele       ELEMENT*      (i)  current element
\param  data      SO3_DATA*     (i)  Gauss point data
\param  ex        DOUBLE[]      (i)  material element node coord.
\param  exs       DOUBLE[]      (i)  spatial element node coord.s
\param  timen     DOUBLE        (i)  curr. load factor/curr. time
\param  ngsurf    INT           (i)  number of element surfaces
\param  gsurf     GSURF*[]      (i)  geometry surf.
\param  eload     DOUBLE[][]    (io) element load vector (is changed)
\return void

\author bborn
\date 04/07
*/
void so3_load_surf_int(ELEMENT* ele,
                       SO3_DATA* data,
                       DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                       DOUBLE exs[MAXNOD_SOLID3][NDIM_SOLID3],
                       const DOUBLE timen,
                       INT ngsurf,
                       GSURF* gsurf[MAXSID_SOLID3],
                       DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3])
{
  /* element */
  const DIS_TYP distyp = ele->distyp;  /* local copy of discretisation type */

  INT curve;  /* curve index */
  DOUBLE cfac;  /* curve factor */

  /* surface */
  INT igsurf;  /* element surface index */
  INT idim;  /* dimension index */
  INT idimsid;  /* dimension of surface index */
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
          /* curve factor */
          curve = gsurf[igsurf]->neum->curve - 1;
          if (curve < 0)
          {
            cfac = 1.0;
          }
          else
          {
            if (genprob.timetyp == time_static)
            {
              /* not implemented, could be something like: */
              /* dyn_facfromcurve(curve, timen, &(cfac)); */
              cfac = 1.0;  /* remove this */
            }
            else if (genprob.timetyp == time_dynamic)
            {
              dyn_facfromcurve(curve, timen, &(cfac));
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
          }
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
                  gpcidim += data->redsidh[igsurf][idimsid][idim] 
                    * data->ghlc[gpintc[idimsid]][igp[idimsid]];
                }
                /* final set of idim-component */
                gpc[idim] = gpcidim;
              }
              /*--------------------------------------------------------*/
              /* Gauss weight */
              fac = cfac;  /* initialise integration factor */
              for (idimsid=0; idimsid<DIMSID_SOLID3; idimsid++)
              {
                /* multiply weight */
                fac *= data->ghlw[gpintc[idimsid]][igp[idimsid]];
              }
              /*--------------------------------------------------------*/
              /* Shape functions at Gauss point */
              so3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1, 
                              shape, deriv);
              /*--------------------------------------------------------*/
              /* add surface load contribution ==> eload modified */
              so3_load_surf_valh(ele, data, igsurf, gsurf[igsurf], 
                                ex, exs, gpc, shape, deriv,
                                fac, eload);
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
          /* curve factor */
          curve = gsurf[igsurf]->neum->curve - 1;
          if (curve < 0)
          {
            cfac = 1.0;
          }
          else
          {
            if (genprob.timetyp == time_static)
            {
              /* not implemented, could be something like: */
              /* dyn_facfromcurve(curve, timen, &(cfac)); */
              cfac = 1.0;  /* remove this */
            }
            else if (genprob.timetyp == time_dynamic)
            {
              dyn_facfromcurve(curve, timen, &(cfac));
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
                  + data->redsidt[igsurf][idimsid][idim]
                    * data->gtsc[gpintc[0]][igp[0]][idimsid];
              }  /* end for */
              /* final set of idim-component */
              gpc[idim] = gpcidim;
            }  /* end for */
            /*----------------------------------------------------------*/
            /* multiply weight */
            fac = cfac * data->gtsw[gpintc[0]][igp[0]];
            /*----------------------------------------------------------*/
            /* Shape functions at Gauss point */
            so3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1, 
                            shape, deriv);
            /*----------------------------------------------------------*/
            /* add surface load contribution ==> eload modified */
            so3_load_surf_valt(ele, data, igsurf, gsurf[igsurf], 
                              ex, exs, gpc, shape, deriv, 
                              fac, eload);
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
} /* end of so3_load_surf_int */


/*======================================================================*/
/*!
\brief Determine load due to tractions on element faces (surfaces)
       (Hexahedra)

\param    ele      ELEMENT*     (i)    actual element
\param    data     SO3_DATA*    (i)    Gauss point data
\param    igsurf   INT          (i)    curr. geom. surf. index
\param    gsurf    GSURF*       (i)    current geometry surface
\param    ex       DOUBLE[][]   (i)    material (undeformed) node coord.s
\param    exs      DOUBLE[][]   (i)    spatial (deformed) node coord.s
\param    gpc      DOUBLE[]     (i)    r,s,t-coord current GP
\param    shape    DOUBLE[]     (i)    shape function at Gauss point
\param    deriv    DOUBLE[][]   (i)    derivatives of shape fct. at Gauss point
\param    fac      DOUBLE       (i)    integration factor at Gauss point
\param    eload    DOUBLE[][]   (io)   element load vector contribution
\return void

\author bborn
\date 04/07
*/
void so3_load_surf_valh(ELEMENT* ele,
                        SO3_DATA* data,
                        INT igsurf,
                        GSURF* gsurf,
                        DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE exs[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE gpc[NDIM_SOLID3],
                        DOUBLE shape[MAXNOD_SOLID3],
                        DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE fac,
                        DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3])
{
  const DOUBLE pressureminus = -1.0;
  const INT nelenod = ele->numnp;  /* number of element nodes */
  DOUBLE metr;  /* metric */
  INT onoff[NUMDOF_SOLID3];
  DOUBLE traction[NUMDOF_SOLID3];  /* elem. bound. traction [force/area] */
  INT idof, jdof, inode;  /* loopers */
  DOUBLE gamt[DIMSID_SOLID3][NDIM_SOLID3];  /* base vectors physically oriented */
  DOUBLE unrm[NDIM_SOLID3];  /* unit normal */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load_surf_valh");
#endif

  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(gsurf->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform traction applied in undeformed configuration */
    case neum_dead:
      /* metric at Gauss point */
      so3_metr_surf(ele, nelenod, ex, deriv, data->redsidh[igsurf], 
                    NULL, &metr);
      /* copy switch and values */
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
            eload[inode][idof] += shape[inode] * traction[idof] * fac * metr;
          }
        }
      }
      break;
    /*------------------------------------------------------------------*/
    /* uniform orthogonal pressure on deformed configuration */
    case neum_orthopressure:
      idof = 0;  /* only 1st index relevant */
      onoff[idof] = gsurf->neum->neum_onoff.a.iv[idof];
      if (onoff[idof] == 1)
      {
        /* pressure value */
        traction[idof] = pressureminus * gsurf->neum->neum_val.a.dv[idof];
        /* physically oriented parametric base vectors */
        so3_metr_surf(ele, nelenod, exs, deriv, 
                      data->redsidh[igsurf], gamt, &metr);
        /* unit normal in deformed space */
        so3_tns3_unrm(gamt[0], gamt[1], data->nrmsidh[igsurf], unrm);
        /* apply load to element load vector */
        for (inode=0; inode<nelenod; inode++)
        {
          for (jdof=0; jdof<NUMDOF_SOLID3; jdof++)
          {
            eload[inode][jdof] += shape[inode] * traction[idof] 
              * unrm[jdof] * fac * metr;
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
/*!
\brief Determine load due to tractions on element faces (surfaces)
       (Tetrahedra)

\param    ele      ELEMENT*     (i)    actual element
\param    data     SO3_DATA*    (i)    Gauss point data
\param    igsurf   INT          (i)    curr. geom. surf. index
\param    gsurf    GSURF*       (i)    current geometry surface
\param    ex       DOUBLE[][]   (i)    material (undeformed) node coord.s
\param    exs      DOUBLE[][]   (i)    spatial (deformed) node coord.s
\param    gpc      DOUBLE[]     (i)    r,s,t-coord current GP
\param    shape    DOUBLE[]     (i)    shape function at Gauss point
\param    deriv    DOUBLE[][]   (i)    derivatives of shape fct. at Gauss point
\param    fac      DOUBLE       (i)    integration factor at Gauss point
\param    eload    DOUBLE[][]   (io)   element load vector contribution

\author bborn
\date 04/07
*/
void so3_load_surf_valt(ELEMENT* ele,
                        SO3_DATA* data,
                        INT igsurf,
                        GSURF* gsurf,
                        DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE exs[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE gpc[NDIM_SOLID3],
                        DOUBLE shape[MAXNOD_SOLID3],
                        DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE fac,
                        DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3])
{
  const DOUBLE pressureminus = -1.0;
  const INT nelenod = ele->numnp;  /* number of element nodes */
  DOUBLE metr;  /* metric */
  INT onoff[NUMDOF_SOLID3];
  DOUBLE traction[NUMDOF_SOLID3];  /* elem. bound. traction [force/area] */
  INT idof, jdof, inode;  /* loopers */
  DOUBLE gamt[DIMSID_SOLID3][NDIM_SOLID3];  /* base vectors physically oriented */
  DOUBLE unrm[NDIM_SOLID3];  /* unit normal */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load_surf_valt");
#endif

  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(gsurf->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform traction applied in undeformed configuration */
    case neum_dead:
      /* metric at Gauss point */
      so3_metr_surf(ele, nelenod, ex, deriv, data->redsidt[igsurf], 
                    NULL, &metr);
      /* copy switch and values */
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
            eload[inode][idof] += shape[inode] * traction[idof] * fac * metr;
          }
        }
      }
      break;
    /*------------------------------------------------------------------*/
    /* uniform orthogonal pressure on deformed configuration */
    case neum_orthopressure:
      idof = 0;  /* only 1st index relevant */
      onoff[idof] = gsurf->neum->neum_onoff.a.iv[idof];
      if (onoff[idof] == 1)
      {
        /* pressure value */
        traction[idof] = pressureminus * gsurf->neum->neum_val.a.dv[idof];
        /* physically oriented parametric base vectors */
        so3_metr_surf(ele, nelenod, exs, deriv, 
                      data->redsidt[igsurf], gamt, &metr);
        /* unit normal in deformed space */
        so3_tns3_unrm(gamt[0], gamt[1], data->nrmsidt[igsurf], unrm);
        /* apply load to element load vector */
        for (inode=0; inode<nelenod; inode++)
        {
          for (jdof=0; jdof<NUMDOF_SOLID3; jdof++)
          {
            eload[inode][jdof] += shape[inode] * traction[idof] 
              * unrm[jdof] * fac * metr;
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
