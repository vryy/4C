/*======================================================================*/
/*!
\file
\brief Spatial integration of loads (ie heat sources/fluxes) applied
       to element domain (volume), sides and edges

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/
#ifndef CCADISCRET
#ifdef D_THERM3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm3.h"

/*!
\addtogroup THERM3
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data
       global variable GENPROB genprob is defined in global_control.c
\author bborn
\date 03/06
*/
extern GENPROB genprob;

/*----------------------------------------------------------------------*/
/*!
\brief Alldyn dynamic control
       pointer to allocate dynamic variables if needed
       dedfined in global_control.c
\auther bborn
\date 05/07
*/
extern ALLDYNA* alldyn;

/*----------------------------------------------------------------------*/
/*!
\brief Spatial distributions
\author bborn
\date 10/07
*/
extern INT numfunct;
extern FUNCT* funct;


/*======================================================================*/
/*!
\brief Spatial integration of
       (i)    heat source in element domain (volume)
       (ii)   heat flux on element sides (surfaces)
       (iii)  heat flux on element edges (lines)

The heat source is integration over the element domain (surface).
The heat flux is integratted along the elements boundaries (edges/lines).
The integration results in the external element heat load vector.

The parameter space is defined by the triple (r,s,t)
Hexahedra biunit cube  { (r,s,t) | -1<=r<=1, -1<=s<=1, -1<=t<=1 }
Tetrahedra  { (r,s,t) | -1<=r<=1, -1<=s<=1-r, -1<=t<=1-r-s }

\param   *ele           ELEMENT     (i)  pointer to current element
\param   *data          TH3_DATA    (i)  common element data
\param    imyrank       INT         (i)  ??????? parallel stuff
\param   *loadvec       DOUBLE      (o)  global element load vector fext
\return void

\author bborn
\date 09/06
*/
void th3_load_heat(CONTAINER* container,
                   ELEMENT *ele,
                   TH3_DATA *data,
                   INT imyrank,
                   ARRAY* intforce_global)
{
  DOUBLE* loadvec = intforce_global->a.dv;  /* internal load vector */
  const INT itdyn = genprob.numtf;  /* index of therm dynamics data */
  const THERM_DYNAMIC* tdyn = alldyn[itdyn].tdyn;  /* therm. dyn. data */
  const DOUBLE timen = tdyn->time;  /* current time */
  INT curve;  /* curve index */
  DOUBLE cfac;  /* curve factor */

  /* general variables */
  INT idim, inode, idof;  /* some counters */
  const INT nelenod = ele->numnp;  /* number of element nodes */
  /* const INT neledof = nelenod * NUMDOF_THERM3;  /\* element DOF *\/ */
  const DIS_TYP distyp = ele->distyp;  /* elem. discretisation type */

  /* volume load */
  GVOL* gvol = ele->g.gvol;  /* local pointer to geometry volume of element */
  INT foundgvolneum;  /* flag for identifying loaded volume */
  DOUBLE shape[MAXNOD_THERM3];  /* shape functions */
  DOUBLE deriv[MAXNOD_THERM3][NDIM_THERM3];
  DOUBLE xjm[NDIM_THERM3][NDIM_THERM3];  /* Jacobian matrix */
  DOUBLE det;  /* Jacobi determinant */
  DOUBLE xji[NDIM_THERM3][NDIM_THERM3];  /* inverse Jacobian matrix */

  /* surface load */
  INT foundgsurfneum;  /* flag for identifying loaded surfaces */
  INT ngsurf;  /* number of geometry surfaces of volumetric element */
  INT igsurf;  /* surface index */
  INT idimsid;  /* surface dimension index */
  GSURF* gsurf[MAXSID_THERM3];  /* geometry surfaces of element */
  DOUBLE sidredm[DIMSID_THERM3][NDIM_THERM3];  /* dimens. reduct. matrix */
  DOUBLE metr;  /* metric */
  DOUBLE gpcidim;  /* dummy Gauss point coordinate */

  /* line load */
  INT foundglineneum;  /* flag for identifying loaded line */
  INT ngline;  /* number of geometry line of volumetric element */
  INT igline;  /* line index */
  GLINE* gline[MAXEDG_THERM3];  /* local pointers to lines */
  DOUBLE linredv[NDIM_THERM3];  /* dimension reduction vector */

  /* integration */
  INT igp[NDIM_THERM3];  /* Gauss point indices */
  INT gpnum[NDIM_THERM3];  /* Gauss point numbers */
  INT gpintc[NDIM_THERM3];  /* Gauss point integration cases */
  DOUBLE fac;  /* integration factors */
  DOUBLE gpc[NDIM_THERM3];  /* r,s,t-coord current GP */

  /* result */
  DOUBLE eload[NUMDOF_THERM3][MAXNOD_THERM3];  /* element load */


  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_load_heat");
#endif

  /*--------------------------------------------------------------------*/
  /* initialize vectors */
  {
    INT idof;
    for (idof=0; idof<NUMDOF_THERM3; idof++)
    {
      INT jnod;
      for (jnod=0; jnod<MAXNOD_THERM3; jnod++)
      {
        eload[idof][jnod] = 0.0;  /* set eload to zero */
      }
    }
  }

  /*====================================================================*/
  /* check if external heat is applied (load) */
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
        case neum_heatconvection:
          foundgsurfneum = 1;
          break;
        case neum_heatconv_fct1:
          foundgsurfneum = 1;
          break;
        case neum_heatconv_fct2:
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

  /*====================================================================*/
  /* domain load ==> volume heat load, heat source */
  /*------------------------------------------------------------------- */
  /* integrate volume load */
  if ( (foundgvolneum > 0) && (imyrank == ele->proc) )
  {
    /*------------------------------------------------------------------*/
    /* number of Gauss points */
    if ( (distyp == hex8) || (distyp == hex20) || (distyp == hex27) )
    {
      for (idim=0; idim<NDIM_THERM3; idim++)
      {
        gpnum[idim] = ele->e.th3->gpnum[idim];
        gpintc[idim] = ele->e.th3->gpintc[idim];
      }
    }
    else if ( (distyp == tet4) || (distyp == tet10) )  /* tet */
    {
      gpnum[0] = 1;
      gpintc[0] = 1;
      gpnum[1] = 1;
      gpintc[1] = 1;
      gpnum[2] = ele->e.th3->gpnum[0];
      gpintc[2] = ele->e.th3->gpintc[0];
    }
    else
    {
      dserror("ele->distyp unknown!");
    } /* end of if */
    /*------------------------------------------------------------------*/
    /* load curve factor at current time */
    curve = gvol->neum->curve - 1;  /* load curve index */
    if (curve < 0)
    {
      cfac = 1.0;  /* default */
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
    /*------------------------------------------------------------------*/
    /* integration loops */
    for (igp[0]=0; igp[0]<gpnum[0]; igp[0]++)
    {
      for (igp[1]=0; igp[1]<gpnum[1]; igp[1]++)
      {
        for (igp[2]=0; igp[2]<gpnum[2]; igp[2]++)
        {
          /*------------------------------------------------------------*/
          /* initialise intgration factor */
          fac = 1.0;
          /*------------------------------------------------------------*/
          /* obtain current Gauss coordinates and weights */
          switch (distyp)
          {
            /* hexahedra */
            case hex8: case hex20: case hex27:
              for (idim=0; idim<NDIM_THERM3; idim++)
              {
                /* r,s,t-coordinate */
                gpc[idim] = data->ghlc[gpintc[idim]][igp[idim]];
                /* weight */
                fac *= data->ghlw[gpintc[idim]][igp[idim]];
              }
              break;
            /* tetrahedra */
            case tet4: case tet10:
              gpc[0] = data->gtdc[gpintc[2]][igp[2]][0];  /* r-coordinate */
              gpc[1] = data->gtdc[gpintc[2]][igp[2]][1];  /* s-coordinate */
              gpc[2] = data->gtdc[gpintc[2]][igp[2]][2];  /* t-coordinate */
              fac = data->gtdw[gpintc[2]][igp[2]];  /* weight */
              break;
            default:
              dserror("ele->distyp unknown!");
          }  /* end of switch */
          /*------------------------------------------------------------*/
          /* shape functions */
          th3_shape_deriv(ele->distyp, gpc[0], gpc[1], gpc[2], 1,
                          shape, deriv);
          /*------------------------------------------------------------*/
          /* compute Jacobian matrix, its determinant
           * inverse Jacobian is not calculated */
          th3_metr_jaco(ele, nelenod, deriv, 0, xjm, &det, xji);
          /*------------------------------------------------------------*/
          /* integration (quadrature) factor */
          fac *= det;
          /*------------------------------------------------------------*/
          /* volume-load  ==> eload modified */
          th3_load_vol(ele, nelenod, shape, fac, cfac, eload);
        }  /* end of for */
      }  /* end of for */
    }  /* end of for */
    /*------------------------------------------------------------------*/
    /* the volume load of this element has been done,
     * -> switch off the volume load condition
     * will be switched on again at beginning of next time step */
    gvol->neum = NULL;
  } /* end of if ((foundsurface > 0) && (imyrank==ele->proc)) */


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
            /*----------------------------------------------------------*/
            /* get dimension reduction matrix */
            /* sidredm = &&(data->redsidh[igsurf][0][0]); */
            for (idimsid=0; idimsid<DIMSID_THERM3; idimsid++)
            {
              for (idim=0; idim<NDIM_THERM3; idim++)
              {
                sidredm[idimsid][idim]
                  = data->redsidh[igsurf][idimsid][idim];
              }
            }
            /*----------------------------------------------------------*/
            /* get number of Gauss points in each direction */
            idimsid = 0;  /* initialiase dimension index */
            for (idim=0; idim<NDIM_THERM3; idim++)
            {
              /* set number of Gauss points for side integration */
              if ((INT) data->ancsidh[igsurf][idim] == 0)
              {
                gpnum[idimsid] = ele->e.th3->gpnum[idim];
                gpintc[idimsid] = ele->e.th3->gpintc[idim];
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
                for (idim=0; idim<NDIM_THERM3; idim++)
                {
                  /* set anchor point in current side
                   * This is the intersection of side plane with its
                   * normal parameter coordinate axis */
                  gpcidim = data->ancsidh[igsurf][idim];
                  for (idimsid=0; idimsid<DIMSID_THERM3; idimsid++)
                  {
                    /* add coordinate components */
                    gpcidim += sidredm[idimsid][idim]
                        * data->ghlc[gpintc[idimsid]][igp[idimsid]];
                  }  /* end for */
                  /* final set of idim-component */
                  gpc[idim] = gpcidim;
                }  /* end for */
                /*------------------------------------------------------*/
                /* Gauss weight */
                fac = 1.0;  /* initialise integration factor */
                for (idimsid=0; idimsid<DIMSID_THERM3; idimsid++)
                {
                  /* multiply weight */
                  fac *= data->ghlw[gpintc[idimsid]][igp[idimsid]];
                }
                /*------------------------------------------------------*/
                /* Shape functions at Gauss point */
                th3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1,
                                shape, deriv);
                /*------------------------------------------------------*/
                /* Jacobi matrix and determinant */
                th3_metr_surf(ele, nelenod, deriv, sidredm, &metr);
                /*------------------------------------------------------*/
                /* integration factor */
                fac *= metr;
                /*------------------------------------------------------*/
                /* add surface load contribution ==> eload modified */
                th3_load_surf(container, ele, nelenod, gsurf[igsurf],
                              gpc, shape, fac, cfac, eload);
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
        gpnum[0] = ele->e.th3->gpnum[1];
        gpintc[0] = ele->e.th3->gpintc[1];
        /*--------------------------------------------------------------*/
        /* loop all element surfaces */
        for (igsurf=0; igsurf<ngsurf; igsurf++)
        {
          /*------------------------------------------------------------*/
          /* check if current side is subjected to a load */
          if (gsurf[igsurf]->neum != NULL)
          {
            /*----------------------------------------------------------*/
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
            /*----------------------------------------------------------*/
            /* set dimension reduction matrix */
            for (idimsid=0; idimsid<DIMSID_THERM3; idimsid++)
            {
              for (idim=0; idim<NDIM_THERM3; idim++)
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
              for (idim=0; idim<NDIM_THERM3; idim++)
              {
                /* set anchor point in current side */
                gpcidim = data->ancsidt[igsurf][idim];
                /* add position sideways */
                for (idimsid=0; idimsid<DIMSID_THERM3; idimsid++)
                {
                  gpcidim += sidredm[idim][idimsid]
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
              th3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1,
                              shape, deriv);
              /*--------------------------------------------------------*/
              /* Jacobi matrix and determinant */
              th3_metr_surf(ele, nelenod, deriv, sidredm, &metr);
              /*--------------------------------------------------------*/
              /* integration factor */
              fac *= metr;
              /*--------------------------------------------------------*/
              /* add surface load contribution ==> eload modified */
              th3_load_surf(container, ele, nelenod, gsurf[igsurf], 
                            gpc, shape, fac, cfac, eload);
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
            /* curve factor */
            curve = gline[igline]->neum->curve - 1;
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
            /*----------------------------------------------------------*/
            /* get dimension reduction matrix */
            for (idim=0; idim<NDIM_THERM3; idim++)
            {
              linredv[idim] = data->rededgh[igline][idim];
            }
            /*----------------------------------------------------------*/
            /* get number of Gauss points for current line
             * these could change depending on the set Gauss point numbers
             * in r-, s- and t-direction */
            for (idim=0; idim<NDIM_THERM3; idim++)
            {
              /* set number of Gauss points for side integration
               * the anchor ancsidh is normal on the edge direction */
              if ((INT) data->ancedgh[igline][idim] == 0)
              {
                gpnum[0] = ele->e.th3->gpnum[idim];
                gpintc[0] = ele->e.th3->gpintc[idim];
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
              for (idim=0; idim<NDIM_THERM3; idim++)
              {
                /* set anchor point in current side
                 * This is the intersection of edge with its normal
                 * parameter coordinate plane */
                gpcidim = data->ancedgh[igline][idim];
                /* add coordinate components */
                gpcidim += linredv[idim] * data->ghlc[gpintc[0]][igp[0]];
                /* final set of idim-component */
                gpc[idim] = gpcidim;
              }  /* end for */
              /*--------------------------------------------------------*/
              /* set weight */
              fac = data->ghlw[gpintc[0]][igp[0]];
              /*--------------------------------------------------------*/
              /* Shape functions at Gauss point */
              th3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1,
                              shape, deriv);
              /*--------------------------------------------------------*/
              /* Jacobi matrix and determinant */
              th3_metr_line(ele, nelenod, deriv, linredv, &metr);
              /*--------------------------------------------------------*/
              /* integration factor */
              fac *= metr;
              /*--------------------------------------------------------*/
              /* add surface load contribution ==> eload modified */
              th3_load_line(ele, nelenod, gline[igline], shape, fac, cfac,
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
        gpnum[0] = ele->e.th3->gpnum[2];
        gpintc[0] = ele->e.th3->gpintc[2];
        /*--------------------------------------------------------------*/
        /* loop all element lines */
        for (igline=0; igline<ngline; igline++)
        {
          /*------------------------------------------------------------*/
          /* check if current line is subjected to a load */
          if (gline[igline]->neum != NULL)
          {
            /*----------------------------------------------------------*/
            /* curve factor */
            curve = gline[igline]->neum->curve - 1;
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
            /*----------------------------------------------------------*/
            /* get dimension reduction matrix */
            for (idim=0; idim<NDIM_THERM3; idim++)
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
              for (idim=0; idim<NDIM_THERM3; idim++)
              {
                /* set anchor point in current side
                 * This is the intersection of edge with its normal
                 * parameter coordinate plane */
                gpcidim = data->ancedgt[igline][idim];
                /* add coordinate components */
                gpcidim += linredv[idim] * data->gtlc[gpintc[0]][igp[0]];
                /* final set of idim-component */
                gpc[idim] = gpcidim;
              }  /* end for */
              /*--------------------------------------------------------*/
              /* set weight */
              fac = data->gtlw[gpintc[0]][igp[0]];
              /*--------------------------------------------------------*/
              /* shape functions at Gauss point */
              th3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1,
                              shape, deriv);
              /*--------------------------------------------------------*/
              /* Jacobi matrix and determinant */
              th3_metr_line(ele, nelenod, deriv, linredv, &metr);
              /*--------------------------------------------------------*/
              /* integration factor */
              fac *= metr;
              /*--------------------------------------------------------*/
              /* add surface load contribution ==> eload modified */
              th3_load_line(ele, nelenod, gline[igline], shape, fac, cfac,
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
      for (idof=0; idof<NUMDOF_THERM3; idof++)
      {
        loadvec[inode*NUMDOF_THERM3+idof] = eload[idof][inode];
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* finish */
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th3_load_heat */



/*======================================================================*/
/*!
\brief Determine load due to heat source on element domain (volume)

\param   *ele      ELEMENT      (i)    actual element
\param    nelenod  INT          (i)    number of element nodes
\param   *shape    DOUBLE       (i)    shape function at Gauss point
\param    fac      DOUBLE       (i)    integration factor
\param  **eload    DOUBLE       (io)   element load vector contribution
\return void

\author bborn
\date 09/06
*/
void th3_load_vol(ELEMENT *ele,
                  INT nelenod,
                  DOUBLE shape[MAXNOD_THERM3],
                  DOUBLE fac,
                  DOUBLE cfac,
                  DOUBLE eload[NUMDOF_THERM3][MAXNOD_THERM3])
{
  const INT heatminus = -1.0;  /* heat minus occuring in heat conduction */
  DOUBLE heatsource[NUMDOF_THERM3];
  INT idof, inode;  /* loopers (i = loaddirection x or y)(j=node) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_load_vol");
#endif
  /*--------------------------------------------------------------------*/
  /* apply heat minus */
  fac *= heatminus;
  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(ele->g.gvol->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform prescribed surface load */
    case neum_dead:
      for (idof=0; idof<NUMDOF_THERM3; idof++)
      {
        heatsource[idof] = ele->g.gvol->neum->neum_val.a.dv[idof];
      }
      /* add load vector component to element load vector */
      for (inode=0; inode<nelenod; inode++)
      {
        for (idof=0; idof<NUMDOF_THERM3; idof++)
        {
          eload[idof][inode] += shape[inode] * heatsource[idof] * fac * cfac;
        }
      }
      break;
    /*------------------------------------------------------------------*/
    case neum_live:
      dserror("load case unknown");
      break;
    /*------------------------------------------------------------------*/
    case neum_consthydro_z:
      dserror("load case unknown");
      break;
    /*------------------------------------------------------------------*/
    case neum_increhydro_z:
      dserror("load case unknown");
      break;
    /*------------------------------------------------------------------*/
    default:
      dserror("load case unknown");
      break;
  }  /* end of switch(ele->g.gsurf->neum->neum_type) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th3_load_vol(...) */


/*======================================================================*/
/*!
\brief Determine load due to heat fluxes on element sides (surfaces)

\param   container  CONTAINER*  (i)    container
\param   ele        ELEMENT*    (i)    actual element
\param   nelenod    INT         (i)    number of element nodes
\param   gsurf      GSURF*      (i)    current geometry surface
\param   gpc        DOUBLE[]    (i)    r,s,t-coord current GP
\param   shape      DOUBLE[]    (i)    shape function at Gauss point
\param   fac        DOUBLE      (i)    integration factor at Gauss point
\param   cfac       DOUBLE      (i)    load curve scale
\param   eload      DOUBLE[][]  (io)   element load vector contribution
\return void

\author bborn
\date 09/06
*/
void th3_load_surf(CONTAINER* container,
                   ELEMENT* ele,
                   INT nelenod,
                   GSURF* gsurf,
                   DOUBLE gpc[NDIM_THERM3],
                   DOUBLE shape[MAXNOD_THERM3],
                   DOUBLE fac,
                   DOUBLE cfac,
                   DOUBLE eload[NUMDOF_THERM3][MAXNOD_THERM3])
{
  const INT heatminus = -1.0;
  INT onoff[NUMDOF_THERM3];
  DOUBLE heatflux[NUMDOF_THERM3];
  INT idof, inode;  /* loopers */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_load_surf");
#endif
  /*--------------------------------------------------------------------*/
  /* apply heat minus */
  fac *= heatminus;
  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(gsurf->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform prescribed surface load */
    case neum_dead:
      for (idof=0; idof<NUMDOF_THERM3; idof++)  /* mark NUMDOF_THERM3==1 */
      {
        onoff[idof] = gsurf->neum->neum_onoff.a.iv[idof];
        heatflux[idof] = gsurf->neum->neum_val.a.dv[idof];
      }
      /* add load vector component to element load vector */
      for (inode=0; inode<nelenod; inode++)
      {
        for (idof=0; idof<NUMDOF_THERM3; idof++)
        {
          /* if load is switched on : apply */
          if (onoff[idof] == 1)
          {
            eload[idof][inode] += shape[inode] * heatflux[idof] * fac * cfac;
          }
        }
      }
      break;
    /*------------------------------------------------------------------*/
    /* heat convection
     * Strictly speaking this is a non-linear heat load type,
     * due to its dependence on the current heat state (follower load-like).
     * Here, we will use the _last converged_ temperature to avoid
     * the effort of a linearisation. */
    case neum_heatconvection:
    case neum_heatconv_fct1:
    case neum_heatconv_fct2:
    {
#if 1
      const INT jdof = 0;
      DOUBLE hflx;
      th3_load_heatconv(container, ele, nelenod, gsurf,
                        gpc, shape, fac, cfac, &hflx);
      /* add load vector component to element load vector */
      for (inode=0; inode<nelenod; inode++)
      {
        eload[jdof][inode] += shape[inode] * hflx * fac;
      }
#else
      const INT jdof = 0;
      INT have_val0 = gsurf->neum->neum_onoff.a.iv[0];
      INT have_val1 = gsurf->neum->neum_onoff.a.iv[1];
      INT have_val2 = gsurf->neum->neum_onoff.a.iv[2];
      INT have_val3 = gsurf->neum->neum_onoff.a.iv[3];
      INT have_val4 = gsurf->neum->neum_onoff.a.iv[4];
      if ( (have_val0)
           && (have_val1) && (have_val2)
           && (have_val3) && (have_val4) )
      {
        /* environmental temperature T_infty */
        DOUBLE envtem = gsurf->neum->neum_val.a.dv[0];
        /* 1st base temperature */
        DOUBLE t1 = gsurf->neum->neum_val.a.dv[1];
        /* 1st convection coefficient */
        DOUBLE c1 = gsurf->neum->neum_val.a.dv[2];
        /* 2nd base temperature */
        DOUBLE t2 = gsurf->neum->neum_val.a.dv[3];
        /* 2nd convection coefficient */
        DOUBLE c2 = gsurf->neum->neum_val.a.dv[4];
        /* temperature of last converged state */
        DOUBLE tem;
        th3_temper_sh(container, ele, shape, &tem);
        /* base temperature */
        DOUBLE bastem = cfac*envtem;
        /* check denominator */
        if (fabs(t2-t1) < EPS12)
        {
          dserror("Division by zero");
        }
        /* convection coefficient */
        DOUBLE convcoeff = ((c2-c1)*bastem + c1*t2 - c2*t1)/(t2-t1);
        /* convective boundary flux */
#if 0
        DOUBLE hflx = convcoeff * bastem;
#else
        DOUBLE hflx = convcoeff * (bastem - tem);
#endif
#if 0
        if (envtem < 0.0)
        {
          DOUBLE hflx = 0.0;
        }
#endif
#if 0
        if ( ele->Id_loc == 713)
          printf("Ele %d: T %g; T_infty %g; Coeff %g; Flux %g\n", ele->Id_loc, tem, bastem, convcoeff, hflx);
        if ( ele->Id_loc == 709)
          printf("Ele %d: T %g; T_infty %g; Coeff %g; Flux %g\n", ele->Id_loc, tem, bastem,  convcoeff, hflx);
#endif
        /* add load vector component to element load vector */
        for (inode=0; inode<nelenod; inode++)
        {
          eload[jdof][inode] += shape[inode] * hflx * fac;
        }
      }
      else if ( (have_val0) && (have_val1) )
      {
        /* environmental temperature T_infty */
        DOUBLE envtem = gsurf->neum->neum_val.a.dv[0];
        /* convection coefficient */
        DOUBLE convcoeff = gsurf->neum->neum_val.a.dv[1];
        /* temperature of last converged state */
        DOUBLE tem;
        th3_temper_sh(container, ele, shape, &tem);
        /* convective boundary flux */
        DOUBLE hflx = convcoeff * (envtem - tem);
        /* add load vector component to element load vector */
        for (inode=0; inode<nelenod; inode++)
        {
          eload[jdof][inode] += shape[inode] * hflx * fac;
        }
      }
      else
      {
#if 0
        dserror("Not enough parameters for heat convection");
#endif
      }
#endif
      break;
    }
    /*------------------------------------------------------------------*/
    case pres_domain_load:
      dserror("load case unknown");
      break;
    /*------------------------------------------------------------------*/
    case neum_consthydro_z:
      dserror("load case unknown");
      break;
    /*------------------------------------------------------------------*/
    case neum_increhydro_z:
      dserror("load case unknown");
      break;
    /*------------------------------------------------------------------*/
    default:
      dserror("load case unknown");
      break;
  }  /* end of switch(ele->g.gsurf->neum->neum_type) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th3_load_vol(...) */


/*======================================================================*/
/*!
\brief Determine load due to heat fluxes on element edges (lines)

\param   *ele      ELEMENT      (i)    actual element
\param    nelenod  INT          (i)    number of element nodes
\param   *igline   GLINE        (i)    current geometry line
\param   *shape    DOUBLE       (i)    shape function at Gauss point
\param    fac      DOUBLE       (i)    integration factor
\param  **eload    DOUBLE       (io)   element load vector contribution
\return void

\author bborn
\date 09/06
*/
void th3_load_line(ELEMENT *ele,
                   INT nelenod,
                   GLINE *gline,
                   DOUBLE shape[MAXNOD_THERM3],
                   DOUBLE fac,
                   DOUBLE cfac,
                   DOUBLE eload[NUMDOF_THERM3][MAXNOD_THERM3])
{
  const INT heatminus = -1.0;
  INT onoff[NUMDOF_THERM3];
  DOUBLE heatflux[NUMDOF_THERM3];
  INT idof, inode;  /* loopers */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_load_line");
#endif
  /*--------------------------------------------------------------------*/
  /* apply heat minus */
  fac *= heatminus;
  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(gline->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform prescribed surface load */
    case neum_dead:
      for (idof=0; idof<NUMDOF_THERM3; idof++)
      {
        onoff[idof] = gline->neum->neum_onoff.a.iv[idof];
        heatflux[idof] = gline->neum->neum_val.a.dv[idof];
      }
      /* add load vector component to element load vector */
      for (inode=0; inode<nelenod; inode++)
      {
        for (idof=0; idof<NUMDOF_THERM3; idof++)
        {
          /* if load is switched on : apply */
          if (onoff[idof] == 1)
          {
            eload[idof][inode] += shape[inode] * heatflux[idof] * fac * cfac;
          }
        }
      }
      break;
    /*------------------------------------------------------------------*/
    case pres_domain_load:
      dserror("load case unknown");
      break;
    /*------------------------------------------------------------------*/
    case neum_consthydro_z:
      dserror("load case unknown");
      break;
    /*------------------------------------------------------------------*/
    case neum_increhydro_z:
      dserror("load case unknown");
      break;
    /*------------------------------------------------------------------*/
    default:
      dserror("load case unknown");
      break;
  }  /* end of switch(ele->g.gsurf->neum->neum_type) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th3_load_line */


/*======================================================================*/
#endif  /*end of #ifdef D_THERM3 */
/*! @} (documentation module close)*/
#endif
