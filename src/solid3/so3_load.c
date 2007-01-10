/*======================================================================*/
/*!
\file
\brief Spatial integration of loads (ie body forces/traction) applied
       to element domain (volume), sides (faces) and edges (lines)

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

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

global variable GENPROB genprob is defined in global_control.c

\author bborn
\date 01/07
*/
extern GENPROB genprob;


/*======================================================================*/
/*!
\brief Spatial integration of 
       (i)    body force in element domain (volume) [force/volume]
       (ii)   traction on element sides (faces) [force/area]
       (iii)  traction on element edges (lines) [force/length]

The traction is integrated over the element domain (surface).
The traction is integrated along the elements boundaries (edges/lines).
The integration results in the external element load vector.

The parameter space is defined by the triple (r,s,t)
Hexahedra triunit cube  { (r,s,t) | -1<=r<=1, -1<=s<=1, -1<=t<=1 }
Tetrahedra  { (r,s,t) | -1<=r<=1, -1<=s<=1-r, -1<=t<=1-r-s }

\param   *ele           ELEMENT     (i)  pointer to current element
\param   *data          SO3_DATA    (i)  common element data
\param 
\param    imyrank       INT         (i)  ??????? parallel stuff
\param   *loadvec       DOUBLE      (o)  global element load vector fext
\return void

\author bborn
\date 09/06
*/
void so3_load(ELEMENT *ele,  /* actual element */
              SO3_DATA *data,
              INT imyrank,
              DOUBLE *loadvec) /* global element load vector fext */
{

  /* general variables */
  INT idim, inode, idof;  /* some counters */
  INT nelenod;  /* number of element nodes */
  INT neledof;  /* element DOF */
  DIS_TYP distyp;  /* local copy of discretisation type */
  DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3];  /* material coord. of element */

  /* volume load */
  GVOL *gvol;  /* local pointer to geometry volume of element */
  INT foundgvolneum;  /* flag for identifying loaded volume */
  DOUBLE shape[MAXNOD_SOLID3];  /* shape functions */
  DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3];
  DOUBLE xjm[NDIM_SOLID3][NDIM_SOLID3];  /* Jacobian matrix */
  DOUBLE det;  /* Jacobi determinant */
  DOUBLE xji[NDIM_SOLID3][NDIM_SOLID3];  /* inverse Jacobian matrix */

  /* surface load */
  INT foundgsurfneum;  /* flag for identifying loaded surfaces */
  INT ngsurf;  /* number of geometry surfaces of volumetric element */
  INT igsurf; /* surface index */
  INT idimsid; /* surface dimension index */
  GSURF *gsurf[MAXSID_SOLID3];
  DOUBLE sidredm[DIMSID_SOLID3][NDIM_SOLID3];  /* dimens. reduct. matrix */
  DOUBLE metr;
  DOUBLE gpcidim;  /* dummy Gauss point coordinate */

  /* line load */
  INT foundglineneum;  /* flag for identifying loaded line */
  INT ngline;  /* number of geometry line of volumetric element */
  INT igline;  /* line index */
  GLINE *gline[MAXEDG_SOLID3];  /* local pointers to lines */
  DOUBLE linredv[NDIM_SOLID3];  /* dimension reduction vector */

  /* integration */
  INT jgp, igp[NDIM_SOLID3];  /* Gauss point indices */
  INT ngp, gpnum[NDIM_SOLID3];  /* Gauss point numbers */
  INT gpintc[NDIM_SOLID3];  /* Gauss point integration cases */
  DOUBLE fac;  /* integration factors */
  DOUBLE gpc[NDIM_SOLID3];  /* r,s,t-coord current GP */

  /* result */
  DOUBLE eload[NUMDOF_SOLID3][MAXNOD_SOLID3];  /* element load */
  

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load");
#endif

  /*--------------------------------------------------------------------*/
  /* element properties */
  nelenod = ele->numnp;  /* element nodes */
  neledof = nelenod * NUMDOF_SOLID3;  /* total number of element DOFs */
  distyp = ele->distyp;
  gvol = ele->g.gvol;

  /*--------------------------------------------------------------------*/
  /* initialize vectors */
  memset(eload, 0, sizeof(eload));  /* set eload to zero */

  /*====================================================================*/
  /* check if external load is applied */
  /*--------------------------------------------------------------------*/
  /* check for presence of heat sources in domain (volume) */
  if (gvol->neum == NULL)
  {
    foundgvolneum = 0;
  }
  else
  {
    switch (ele->g.gvol->neum->neum_type)
    {
      case neum_dead:
        foundgvolneum = 1;
        break;
      default:
        dserror("load type not implemented");
        foundgvolneum = 0;
        break;
    }
  }
  /*--------------------------------------------------------------------*/
  /* number of geom. surfaces */
  ngsurf = gvol->ngsurf;
  /* initialise flag for applied heat fluxes */
  foundgsurfneum = 0;
  /* check if heat fluxes are applied */
  for (igsurf=0; igsurf<ngsurf; igsurf++)
  {
    gsurf[igsurf] = gvol->gsurf[igsurf];
    if (gsurf[igsurf]->neum != NULL)
    {
      switch (gsurf[igsurf]->neum->neum_type)
      {
        case neum_dead:
          foundgsurfneum = 1;
          break;
        default:
          dserror("Neumann BC type is not available on side!");
          break;
      }
    }
  }
  /*--------------------------------------------------------------------*/
  /* number of geom. lines */
  ngline = gvol->ngline;
  /* initialise flag for applied heat fluxes */
  foundglineneum = 0;
  /* check if heat fluxes are applied */
  for (igline=0; igline<ngline; igline++)
  {
    gline[igline] = gvol->gline[igline];
    if (gline[igline]->neum != NULL)
    {
      switch (gline[igline]->neum->neum_type)
      {
        case neum_dead:
          foundglineneum = 1;
          break;
        default:
          dserror("Neumann BC type is not available on edge!");
          break;
      }
    }
  }

  
  /*--------------------------------------------------------------------*/
  /* material coordinates of element nodes */
  if ( (foundgvolneum > 0) 
       || (foundgsurfneum > 0) 
       || (foundglineneum > 0) )
  {
    for (inod=0; inod<nelenod; inod++)
    {
      actnode = ele->node[inod];
      for (jdim=0; jdim<NDIM_SOLID3; jdim++)
      {
        ex[inod][jdim] = actnode->x[jdim];
      }
    }
  }


  /*====================================================================*/
  /* domain load ==> volume heat load, heat source */
  /*------------------------------------------------------------------- */
  /* integrate volume load */
  if ( (foundgvolneum > 0) && (imyrank == ele->proc) )
  {
    /*------------------------------------------------------------------*/
    /* total number of Gauss points */
    ngp = gpshade->gptot;
    /*------------------------------------------------------------------*/
    /* integration loop */
    for (jgp=0; jgp<ngp; jgp++)
    {
      /* initialise intgration factor */
      fac = gpshade->gpwg[jgp];  /* Gauss weight */
      /* compute Jacobian matrix, its determinant
       * inverse Jacobian is not calculated */
      so3_metr_jaco(ele, nelenod, ex, gpshade->gpderiv[igp], 1, 
                    xjm, &det, xji);
      /* integration (quadrature) factor */
      fac = fac * det;
      /* volume-load  ==> eload modified */
      so3_load_vol(ele, nelenod, shape, fac, eload);
    }  /* end of for */
    /*------------------------------------------------------------------*/
    /* the volume load of this element has been done,
     * -> switch off the volume load condition
     * will be switched on again at beginning of next time step */
    gvol->neum = NULL;
  }


  /*====================================================================*/
  /* side loads ==> surface heat load, heat fluxes */
  /*--------------------------------------------------------------------*/
  /* loop all element sides (surfaces) */
  if (foundgsurfneum > 0)
  {
    /*------------------------------------------------------------------*/
    /* distinguish hex and tet */
    switch (distyp)
    {
      /*================================================================*/
      /* hexahedron elements */
      case hex8: case hex20: case hex27:
        /*--------------------------------------------------------------*/
        /* loop all element surfaces */
        for (igsurf=0; igsurf<ngsurf; igsurf++)
        {
          /*------------------------------------------------------------*/
          /* check if current side is subjected to a load */
          if (gsurf[igsurf]->neum != NULL)
          {
            /*----------------------------------------------------------*/
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
            /*----------------------------------------------------------*/
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
            /*----------------------------------------------------------*/
            /* integration loops */
            for (igp[0]=0; igp[0]<gpnum[0]; igp[0]++)
            {
              for (igp[1]=0; igp[1]<gpnum[1]; igp[1]++)
              {
                /*------------------------------------------------------*/
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
                /*------------------------------------------------------*/
                /* Gauss weight */
                fac = 1.0;  /* initialise integration factor */
                for (idimsid=0; idimsid<DIMSID_SOLID3; idimsid++)
                {
                  /* multiply weight */
                  fac = fac * data->ghlw[gpintc[idimsid]][igp[idimsid]];
                }
                /*------------------------------------------------------*/
                /* Shape functions at Gauss point */
                so3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1, 
                                shape, deriv);
                /*------------------------------------------------------*/
                /* Jacobi matrix and determinant */
                so3_metr_surf(ele, nelenod, ex, deriv, sidredm, &metr);
                /*------------------------------------------------------*/
                /* integration factor */
                fac = fac * metr;
                /*------------------------------------------------------*/
                /* add surface load contribution ==> eload modified */
                so3_load_surf(ele, nelenod, gsurf[igsurf], shape, fac, 
                              eload);
              }  /* end for */
            }  /* end for */
            /*----------------------------------------------------------*/
            /* the surface load of this element has been done,
             * -> switch off the surface load condition */
            gsurf[igsurf]->neum = NULL;
          }  /* end if */
        } /* end for */
        break;  /* end cases */
      /*================================================================*/
      /* tetrahedron elements */
      case tet4: case tet10:
        /*--------------------------------------------------------------*/
        /* get number of Gauss points for all sides */
        gpnum[0] = ele->e.so3->gpnum[1];
        gpintc[0] = ele->e.so3->gpintc[1];
        /*--------------------------------------------------------------*/
        /* loop all element surfaces */
        for (igsurf=0; igsurf<ngsurf; igsurf++)
        {
          /*------------------------------------------------------------*/
          /* check if current side is subjected to a load */
          if (gsurf[igsurf]->neum != NULL)
          {
            /*----------------------------------------------------------*/
            /* set dimension reduction matrix */
            for (idimsid=0; idimsid<DIMSID_SOLID3; idimsid++)
            {
              for (idim=0; idim<NDIM_SOLID3; idim++)
              {
                sidredm[idimsid][idim] 
                  = data->redsidt[igsurf][idimsid][idim];
              }
            }
            /*----------------------------------------------------------*/
            /* integration loop */
            for (igp[0]=0; igp[0]<gpnum[0]; igp[0]++)
            {
              /*--------------------------------------------------------*/
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
              /*--------------------------------------------------------*/
              /* multiply weight */
              fac = data->ghlw[gpintc[0]][igp[0]];
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
              so3_load_surf(ele, nelenod, gsurf[igsurf], shape, fac, 
                            eload);
            }  /* end for */
          }  /* end for */
          /*------------------------------------------------------------*/
          /* the surface load of this element has been done,
           * -> switch off the surface load condition */
          gsurf[igsurf]->neum = NULL;
        }  /* end if */
        break;  /* end cases */
      /*----------------------------------------------------------------*/
      default:
        dserror("ele->distyp unknown!");
    }  /* end switch */
  }  /* end if */


  /*====================================================================*/
  /* edge loads ==> edge heat load, heat fluxes */
  /*--------------------------------------------------------------------*/
  /* loop all element sides (surfaces) */
  if (foundglineneum > 0)
  {
    /*------------------------------------------------------------------*/
    /* distinguish hex and tet */
    switch (distyp)
    {
      /*================================================================*/
      /* hexahedron elements */
      case hex8: case hex20: case hex27:
        /*--------------------------------------------------------------*/
        /* loop all element lines */
        for (igline=0; igline<ngline; igline++)
        {
          /*------------------------------------------------------------*/
          /* check if current line is subjected to a load */
          if (gline[igline]->neum != NULL)
          {
            /*----------------------------------------------------------*/
            /* get dimension reduction matrix */
            for (idim=0; idim<NDIM_SOLID3; idim++)
            {
              linredv[idim] = data->rededgh[igline][idim];
            }
            /*----------------------------------------------------------*/
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
            /*----------------------------------------------------------*/
            /* integration loops 
             * For each side one of the following loops is not repeated,
             * but as all sides are generally integrated here. */
            for (igp[0]=0; igp[0]<gpnum[0]; igp[0]++)
            {
              /*--------------------------------------------------------*/
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
              /*--------------------------------------------------------*/
              /* set weight */
              fac = data->ghlw[gpintc[0]][igp[0]];
              /*--------------------------------------------------------*/
              /* Shape functions at Gauss point */
              so3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1, 
                              shape, deriv);
              /*--------------------------------------------------------*/
              /* Jacobi matrix and determinant */
              so3_metr_line(ele, nelenod, ex, deriv, linredv, &metr);
              /*--------------------------------------------------------*/
              /* integration factor */
              fac = fac * metr;
              /*--------------------------------------------------------*/
              /* add surface load contribution ==> eload modified */
              so3_load_line(ele, nelenod, gline[igline], shape, fac, 
                            eload);
            }  /* end for */
            /*----------------------------------------------------------*/
            /* the line load of this element has been done,
             * -> switch off the line load condition */
            gline[igline]->neum = NULL;
          }  /* end if */
        }  /* end for */
        break;  /* end of cases */
      /*================================================================*/
      /* tetrahedron elements */
      case tet4: case tet10:
        /*--------------------------------------------------------------*/
        /* get number of Gauss points for all sides */
        gpnum[0] = ele->e.so3->gpnum[2];
        gpintc[0] = ele->e.so3->gpintc[2];
        /*--------------------------------------------------------------*/
        /* loop all element lines */
        for (igline=0; igline<ngline; igline++)
        {
          /*------------------------------------------------------------*/
          /* check if current line is subjected to a load */
          if (gline[igline]->neum != NULL)
          {
            /*----------------------------------------------------------*/
            /* get dimension reduction matrix */
            for (idim=0; idim<NDIM_SOLID3; idim++)
            {
              linredv[idim] = data->rededgt[igline][idim];
            }
            /*----------------------------------------------------------*/
            /* integration loops 
             * For each side one of the following loops is not repeated,
             * but as all sides are generally integrated here. */
            for (igp[0]=0; igp[0]<gpnum[0]; igp[0]++)
            {
              /*--------------------------------------------------------*/
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
              /*--------------------------------------------------------*/
              /* set weight */
              fac = data->gtlw[gpintc[0]][igp[0]];
              /*--------------------------------------------------------*/
              /* shape functions at Gauss point */
              so3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1, 
                              shape, deriv);
              /*--------------------------------------------------------*/
              /* Jacobi matrix and determinant */
              so3_metr_line(ele, nelenod, ex, deriv, linredv, &metr);
              /*--------------------------------------------------------*/
              /* integration factor */
              fac = fac * metr;
              /*--------------------------------------------------------*/
              /* add surface load contribution ==> eload modified */
              so3_load_line(ele, nelenod, gline[igline], shape, fac, 
                            eload);
            }  /* end for */
            /*----------------------------------------------------------*/
            /* the line load of this element has been done,
             * -> switch off the line load condition */
            gline[igline]->neum = NULL;
          }  /* end if */
        }  /* end for */
        break;  /* end of cases */
      default:
        dserror("ele->distyp unknown!");
    }  /* end switch */
  }  /* end if */

  /*====================================================================*/
  /* add eload to global load vector */
  if ( (foundgvolneum > 0) 
       || (foundgsurfneum > 0) 
       || (foundglineneum > 0) )
  {
    for (inode=0; inode<nelenod; inode++)
    {
      for (idof=0; idof<NUMDOF_SOLID3; idof++)
      {
        loadvec[inode*NUMDOF_SOLID3+idof] = eload[idof][inode];
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* finish */
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_load_heat */



/*======================================================================*/
/*!
\brief Determine load due to body force on element domain (volume)

\param    ele      ELEMENT*     (i)    actual element
\param    nelenod  INT          (i)    number of element nodes
\param    shape    DOUBLE[]     (i)    shape function at Gauss point
\param    fac      DOUBLE       (i)    integration factor
\param    eload    DOUBLE[][]   (io)   element load vector contribution
\return void

\author bborn
\date 01/07
*/
void so3_load_vol(ELEMENT *ele,
                  INT nelenod,
                  DOUBLE shape[MAXNOD_SOLID3],
                  DOUBLE fac,
                  DOUBLE eload[NUMDOF_SOLID3][MAXNOD_SOLID3])
{
  DOUBLE source[NUMDOF_SOLID3];  /* body force [force/vol] */
  INT idof, inode;  /* loopers (i = loaddirection x or y)(j=node) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load_vol");
#endif
  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(ele->g.gvol->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform (density-proportional) dead load */
    case neum_dead:
      for (idof=0; idof<NUMDOF_SOLID3; idof++)
      {
        source[idof] = ele->g.gvol->neum->neum_val.a.dv[idof];
      }
      /* add load vector component to element load vector */
      for (inode=0; inode<nelenod; inode++)
      {
        for (idof=0; idof<NUMDOF_SOLID3; idof++)
        {
          eload[idof][inode] += shape[inode] * source[idof] * fac;
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
} /* end of so3_load_vol(...) */


/*======================================================================*/
/*!
\brief Determine load due to heat fluxes on element sides (surfaces)

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
void so3_load_surf(ELEMENT *ele,
                   INT nelenod,
                   GSURF *gsurf,
                   DOUBLE shape[MAXNOD_SOLID3],
                   DOUBLE fac,
                   DOUBLE eload[NUMDOF_SOLID3][MAXNOD_SOLID3])
{
  INT onoff[NUMDOF_SOLID3];
  DOUBLE traction[NUMDOF_SOLID3];  /* elem. bound. traction [force/area] */
  INT idof, inode;  /* loopers */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load_surf");
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
            eload[idof][inode] += shape[inode] * traction[idof] * fac;
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
} /* end of so3_load_vol(...) */


/*======================================================================*/
/*!
\brief Determine load due to heat fluxes on element edges (lines)

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
void so3_load_line(ELEMENT *ele,
                   INT nelenod,
                   GLINE *gline,
                   DOUBLE shape[MAXNOD_SOLID3],
                   DOUBLE fac,
                   DOUBLE eload[NUMDOF_SOLID3][MAXNOD_SOLID3])
{
  INT onoff[NUMDOF_SOLID3];
  DOUBLE traction[NUMDOF_SOLID3];  /* traction on edge [force/length] */
  INT idof, inode;  /* loopers */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load_line");
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
            eload[idof][inode] += shape[inode] * traction[idof] * fac;
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
} /* end of so3_load_line */


/*======================================================================*/
#endif  /*end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
