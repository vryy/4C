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
\param    imyrank       INT         (i)  ???????
\param   *loadvec       DOUBLE      (o)  global element load vector fext
\return void

\author bborn
\date 09/06
*/
void th3_load_heat(ELEMENT *ele,  /* actual element */
                   TH3_DATA *data,
                   INT imyrank,
                   DOUBLE *loadvec) /* global element load vector fext */
{

  /* general variables */
  const INT heatminus = -1.0;  /* minus sign occuring in heat conduction */
  INT idim, inode, idof;  /* some counters */
  INT nelenod;  /* number of element nodes */
  INT neledof;  /* element DOF */
  DIS_TYP distyp;  /* local copy of discretisation type */

  /* volume load */
  GVOL *gvol;  /* local pointer to geometry volume of element */
  INT foundgvolneum;  /* flag for identifying loaded volume */
  DOUBLE shape[MAXNOD_THERM3];  /* shape functions */
  DOUBLE deriv[MAXNOD_THERM3][NDIM_THERM3];
  DOUBLE xjm[NDIM_THERM3][NDIM_THERM3];  /* Jacobian matrix */
  DOUBLE det;  /* Jacobi determinant */
  DOUBLE xji[NDIM_THERM3][NDIM_THERM3];  /* inverse Jacobian matrix */

  /* surface load */
  INT foundgsurfneum;  /* flag for identifying loaded surfaces */
  INT ngsurf;  /* number of geometry surfaces of volumetric element */
  INT igsurf; /* surface index */
  INT idimsid; /* surface dimension index */
  GSURF *gsurf[MAXSID_THERM3];
  DOUBLE sidredm[DIMSID_THERM3][NDIM_THERM3];  /* dimens. reduct. matrix */
  DOUBLE metr;

  /* line load */
  INT foundglineneum;  /* flag for identifying loaded line */
  INT ngline;  /* number of geometry line of volumetric element */
  INT igline;  /* line index */
  GLINE *gline[MAXEDG_THERM3];  /* local pointers to lines */
  DOUBLE linredv[NDIM_THERM3];  /* dimension reduction vector */

  /* integration */
  INT igp[NDIM_THERM3];  /* Gauss point indices */
  INT gpnum[NDIM_THERM3];  /* Gauss point numbers */
  INT gpintc[NDIM_THERM3];  /* Gauss point integration cases */
  DOUBLE fac;  /* integration factors */
  DOUBLE gpc[NDIM_THERM3];  /* r,s,t-coord current GP */
  INT gpsubnum[NDIM_THERM3];
  INT gpsubintc[NDIM_THERM3];

  /* result */
  DOUBLE eload[NUMDOF_THERM3][MAXNOD_THERM3];  /* element load */
  

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_load_heat");
#endif

  /*--------------------------------------------------------------------*/
  /* element properties */
  nelenod = ele->numnp;  /* element nodes */
  neledof = nelenod * NUMDOF_THERM3;  /* total number of element DOFs */
  distyp = ele->distyp;
  gvol = ele->g.gvol;

  /*--------------------------------------------------------------------*/
  /* initialize vectors */
  memset(eload, 0, sizeof(eload));
/*   for (idim=0; idim<NUMDOF_THERM3; idim++)  /\* element load *\/ */
/*   { */
/*     for (inode=0; inode<nelenod; inode++) */
/*     { */
/*       eload[idim][inode] = 0.0; */
/*     } */
/*   } */

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
      case pres_domain_load:
        foundgvolneum = 1;
        break;
      default:
        dserror("load case not implemented");
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
        case pres_domain_load:
          foundgsurfneum = 1;
          break;
        default:
          dserror("Neumann BC type is not available!");
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
        case pres_domain_load:
          foundglineneum = 1;
          break;
        default:
          dserror("Neumann BC type is not available!");
          break;
      }
    }
  }
  /*--------------------------------------------------------------------*/
  /* Gauss integraton data */
  if ( (foundgvolneum > 0) 
       || (foundgsurfneum > 0) 
       || (foundglineneum > 0) )
  {
    switch (distyp)
    {
      /* hexahedra elements */
      case hex8: case hex20: case hex27:
        for (idim=0; idim<NDIM_THERM3; idim++)
        {
          gpnum[idim] = ele->e.th3->gpnum[idim];
          gpintc[idim] = ele->e.th3->gpintc[idim];
        }
        break;
      /* tetrahedra elements */
      case tet4: case tet10:
        gpnum[0] = 1;
        gpnum[1] = 1;
        /* gpnum[2] has to be set according to domain/side/edge */
        gpintc[0] = 1;
        gpintc[1] = 1;
        /* gpintc[2] has to be set according to domain/side/edge */
        break;
      default:
        dserror("ele->distyp unknown!");
        break;
    }  /* end of switch */
  }  /* end of if */

  /*====================================================================*/
  /* domain load ==> volume heat load, heat source */
  /*------------------------------------------------------------------- */
  /* integrate volume load */
  if ( (foundgvolneum > 0) && (imyrank == ele->proc) )
  {
    /*------------------------------------------------------------------*/
    /* Gauss integraton data for tetrahedron domain */
    if ( (distyp == tet4) || (distyp == tet10) )
    {
        gpnum[2] = ele->e.th3->gpnum[0];
        gpintc[2] = ele->e.th3->gpintc[0];
    }  /* end of if */
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
                fac = fac * data->ghlw[gpintc[idim]][igp[idim]];
              }
              break;
            /* tetrahedra */
            case tet4: case tet10:
              gpc[0] = data->gtdcr[gpintc[2]][igp[2]];  /* r-coordinate */
              gpc[1] = data->gtdcs[gpintc[2]][igp[2]];  /* s-coordinate */
              gpc[2] = data->gtdct[gpintc[2]][igp[2]];  /* t-coordinate */
              fac = data->gtdw[gpintc[2]][igp[2]];  /* weight */
              break;
            default:
              dserror("ele->distyp unknown!");
          }  /* end of switch */
          /*------------------------------------------------------------*/
          /* shape functions */
          th3_shape_deriv(ele->distyp, gpc[0], gpc[1], gpc[2], 0, 
                          shape, deriv);
          /*------------------------------------------------------------*/
          /* compute Jacobian matrix, its determinant
           * inverse Jacobian is not calculated */
          th3_metr_jaco(ele, nelenod, deriv, 0, xjm, &det, xji);
          /*------------------------------------------------------------*/
          /* integration (quadrature) factor */
          fac = heatminus * fac * det;
          /*------------------------------------------------------------*/
          /* volume-load  ==> eload modified */
          th3_load_vol(ele, nelenod, shape, fac, eload);
        }  /* end of for */
      }  /* end of for */
    }  /* end of for */
    /*------------------------------------------------------------------*/
    /* WHY, WHY, WHY ??????????????? */
    /* the volume load of this element has been done,
     * -> switch off the volume load condition */
    gvol->neum = NULL;
  } /* end of if ((foundsurface > 0) && (imyrank==ele->proc)) */


  /*====================================================================*/
  /* side loads ==> surface heat load, heat fluxes */
  /*--------------------------------------------------------------------*/
  /* loop all element sides (surfaces) */
  if (foundgsurfneum > 0)
  {
    /*------------------------------------------------------------------*/
    /* Gauss integraton data for tetrahedron sides */
    if ( (distyp == tet4) || (distyp == tet10) )
    {
        gpnum[2] = ele->e.th3->gpnum[1];
        gpintc[2] = ele->e.th3->gpintc[1];
    }  /* end of if */
    /*------------------------------------------------------------------*/
    /* loop all element surfaces */
    for (igsurf=0; igsurf<ngsurf; igsurf++)
    {
      /*----------------------------------------------------------------*/
      /* get number of Gauss points for current surface
       * these could change depending on the set Gauss point numbers
       * in r-, s- and t-direction */
      switch (distyp)
      {
        /* hexahedron elements */
        case hex8: case hex20: case hex27:
          for (idim=0; idim<NDIM_THERM3; idim++)
          {
            /* set number of Gauss points for side integration */
            if (data->dirsidh[igsurf][idim] == 1)
            {
              gpsubnum[idim] = gpnum[idim];
              gpsubintc[idim] = gpintc[idim];
            }
            else
            {
              gpsubnum[idim] = 1;
              gpsubintc[idim] = 1;
            }
          }
          break;
        /* tetrahedron elements */
        case tet4: case tet10:
          dserror("surface loads for tet.s are not available");
          break;
        /* catch dubious discretisation types */
        default:
          break;
      }
      /*----------------------------------------------------------------*/
      /* integration loops 
       * For each side one of the following loops is not repeated,
       * but as all sides are generally integrated here. */
      for (igp[0]=0; igp[0]<gpsubnum[0]; igp[0]++)
      {
        for (igp[1]=0; igp[1]<gpsubnum[1]; igp[1]++)
        {
          for (igp[2]=0; igp[2]<gpsubnum[2]; igp[2]++)
          {
            /*----------------------------------------------------------*/
            /* set anchor point in current side
             * This is the intersection of side plane with its normal 
             * parameter coordinate axis */
            for (idim=0; idim<NDIM_THERM3; idim++)
            {
              gpc[idim] = data->ancsidh[igsurf][idim];
            }
            /* initialise integration factor */
            fac = 1.0;
            /* initialise dimension reduction matrix to zero */
            memset(sidredm, 0, sizeof(sidredm));
            /* obtain current Gauss coordinates and weights 
             * and dimension reduction matrix */
            switch (ele->distyp)
            {
              /* hexahedra */
              case hex8: case hex20: case hex27:
                /* initialise dimension reduction counter */
                idimsid = 0;
                /* loop potential basis vector directions */
                for (idim=0; idim<NDIM_THERM3; idim++)
                {
                  if (data->dirsidh[igsurf][idim] == 1)
                  {
                    /* add coordinate components */
                    gpc[idim] = gpc[idim] 
                              + data->ghlc[gpsubintc[idim]][igp[idim]];
                    /* add weight */
                    fac = fac * data->ghlw[gpsubintc[idim]][igp[idim]];
                    /* dimension reduction matrix */
                    sidredm[idimsid][idim] = 1.0;
                    idimsid = idimsid + 1;
                  }
                }
                break;
              /* tetrahedra */
              case tet4: case tet10:
                dserror("tets are not available");
                break;
              default:
                dserror("ele->distyp unknown!");
            }  /* end of switch (ele->distyp) */
            /*----------------------------------------------------------*/
            /* Shape functions at Gauss point */
            th3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1, 
                            shape, deriv);
            /*----------------------------------------------------------*/
            /* Jacobi matrix and determinant */
            th3_metr_surf(ele, nelenod, deriv, sidredm, &metr);
            /*----------------------------------------------------------*/
            /* integration factor */
            fac = fac * metr;
            /*----------------------------------------------------------*/
            /* add surface load contribution ==> eload modified */
            th3_load_surf(ele, nelenod, gsurf[igsurf], shape, fac, 
                          eload);
          }  /* end of for */
        }  /* end of for */
      }  /* end of for */
      /*----------------------------------------------------------------*/
      /* WHY, WHY AND WHY ??? */
      /* the surface load of this element has been done,
       * -> switch off the surface load condition */
      gsurf[igsurf]->neum = NULL;
    }
  }
  
  /*====================================================================*/
  /* edge loads ==> edge heat load, heat fluxes */
  
  /*--------------------------------------------------------------------*/
  /* loop all element sides (surfaces) */
  if (foundglineneum > 0)
  {
    /*------------------------------------------------------------------*/
    /* Gauss integraton data for tetrahedron edges */
    if ( (distyp == tet4) || (distyp == tet10) )
    {
        gpnum[2] = ele->e.th3->gpnum[2];
        gpintc[2] = ele->e.th3->gpintc[2];
    }  /* end of if */
    /*------------------------------------------------------------------*/
    /* loop all element lines */
    for (igline=0; igline<ngline; igline++)
    {
      /*----------------------------------------------------------------*/
      /* get number of Gauss points for current line
       * these could change depending on the set Gauss point numbers
       * in r-, s- and t-direction */
      switch (distyp)
      {
        /* hexahedron elements */
        case hex8: case hex20: case hex27:
          for (idim=0; idim<NDIM_THERM3; idim++)
          {
            /* set number of Gauss points for side integration */
            if (data->dirsidh[igsurf][idim] == 1)
            {
              gpsubnum[idim] = gpnum[idim];
              gpsubintc[idim] = gpintc[idim];
            }
            else
            {
              gpsubnum[idim] = 1;
              gpsubintc[idim] = 1;
            }
          }
          break;
        /* tetrahedron elements */
        case tet4: case tet10:
          dserror("tets are not available");
          break;
        /* catch dubious discretisation types */
        default:
          break;
      }
      /*----------------------------------------------------------------*/
      /* integration loops 
       * For each side one of the following loops is not repeated,
       * but as all sides are generally integrated here. */
      for (igp[0]=0; igp[0]<gpsubnum[0]; igp[0]++)
      {
        for (igp[1]=0; igp[1]<gpsubnum[1]; igp[1]++)
        {
          for (igp[2]=0; igp[2]<gpsubnum[2]; igp[2]++)
          {
            /*----------------------------------------------------------*/
            /* set anchor point in current side
             * This is the intersection of side plane with its normal 
             * parameter coordinate axis */
            for (idim=0; idim<NDIM_THERM3; idim++)
            {
              gpc[idim] = data->ancedgh[igline][idim];
            }
            /* initialise integration factor */
            fac = 1.0;
            /* initialise dimension reduction vector to zero */
            memset(linredv, 0, sizeof(linredv));
            /* obtain current Gauss coordinates and weights 
             * and dimension reduction matrix */
            switch (distyp)
            {
              /* hexahedra */
              case hex8: case hex20: case hex27:
                /* loop potential basis vector directions */
                for (idim=0; idim<NDIM_THERM3; idim++)
                {
                  if ( (data->dirsidh[igline][idim] == 1)
                       || (data->dirsidh[igline][idim] == -1) )
                  {
                    /* add coordinate components */
                    gpc[idim] = gpc[idim] 
                              + data->ghlc[gpsubintc[idim]][igp[idim]];
                    /* add weight */
                    fac = fac * data->ghlw[gpsubintc[idim]][igp[idim]];
                    /* dimension reduction vector */
                    linredv[idim] = 1.0;
                  }
                }
                break;
              /* tetrahedra */
              case tet4: case tet10:
                dserror("tets are not available");
                break;
              default:
                dserror("ele->distyp unknown!");
            }  /* end of switch (ele->distyp) */
            /*----------------------------------------------------------*/
            /* Shape functions at Gauss point */
            th3_shape_deriv(distyp, gpc[0], gpc[1], gpc[2], 1, 
                            shape, deriv);
            /*----------------------------------------------------------*/
            /* Jacobi matrix and determinant */
            th3_metr_line(ele, nelenod, deriv, linredv, &metr);
            /*----------------------------------------------------------*/
            /* integration factor */
            fac = fac * metr;
            /*----------------------------------------------------------*/
            /* add surface load contribution ==> eload modified */
            th3_load_line(ele, nelenod, gline[igline], shape, fac, eload);
          }  /* end of for */
        }  /* end of for */
      }  /* end of for */
      /*----------------------------------------------------------------*/
      /* WHY, WHY AND WHY ??? */
      /* the line load of this element has been done,
       * -> switch off the line load condition */
      gline[igline]->neum = NULL;
    }
  }

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
        loadvec[inode*NUMDOF_THERM3+idof] += eload[idof][inode];
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

\param	 *ele      ELEMENT      (i)    actual element
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
                  DOUBLE eload[NUMDOF_THERM3][MAXNOD_THERM3])
{
  DOUBLE heatsource[NUMDOF_THERM3];
  INT idof, inode;  /* loopers (i = loaddirection x or y)(j=node) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_load_vol");
#endif
  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(ele->g.gvol->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform prescribed surface load */
    case pres_domain_load:
      for (idof=0; idof<NUMDOF_THERM3; idof++)
      {
        heatsource[idof] = ele->g.gvol->neum->neum_val.a.dv[idof];
      }
      /* add load vector component to element load vector */
      for (inode=0; inode<nelenod; inode++)
      {
        for (idof=0; idof<NUMDOF_THERM3; idof++)
        {
          eload[idof][inode] += shape[inode] * heatsource[idof] * fac;
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

\param	 *ele      ELEMENT      (i)    actual element
\param    nelenod  INT          (i)    number of element nodes
\param   *gsurf    GSURF        (i)    current geometry surface
\param   *shape    DOUBLE       (i)    shape function at Gauss point
\param    fac      DOUBLE       (i)    integration factor at Gauss point
\param  **eload    DOUBLE       (io)   element load vector contribution
\return void

\author bborn
\date 09/06
*/
void th3_load_surf(ELEMENT *ele,
                   INT nelenod,
                   GSURF *gsurf,
                   DOUBLE shape[MAXNOD_THERM3],
                   DOUBLE fac,
                   DOUBLE eload[NUMDOF_THERM3][MAXNOD_THERM3])
{
  INT onoff[NUMDOF_THERM3];
  DOUBLE heatflux[NUMDOF_THERM3];
  INT idof, inode;  /* loopers */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_load_surf");
#endif
  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(gsurf->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform prescribed surface load */
    case pres_domain_load:
      for (idof=0; idof<NUMDOF_THERM3; idof++)
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
            eload[idof][inode] += shape[inode] * heatflux[idof] * fac;
          }
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
\brief Determine load due to heat fluxes on element edges (lines)

\param	 *ele      ELEMENT      (i)    actual element
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
                   DOUBLE eload[NUMDOF_THERM3][MAXNOD_THERM3])
{
  INT onoff[NUMDOF_THERM3];
  DOUBLE heatflux[NUMDOF_THERM3];
  INT idof, inode;  /* loopers */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_load_line");
#endif
  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(gline->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform prescribed surface load */
    case pres_domain_load:
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
            eload[idof][inode] += shape[inode] * heatflux[idof] * fac;
          }
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
} /* end of th3_load_line */


/*======================================================================*/
#endif  /*end of #ifdef D_THERM3 */
/*! @} (documentation module close)*/
