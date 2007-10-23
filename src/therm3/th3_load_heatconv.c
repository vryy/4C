/*======================================================================*/
/*!
\file
\brief Boundary condition : Heat convection

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
\brief Heat convection BC

Strictly speaking this is a non-linear heat load type,
due to its dependence on the current heat state (follower load-like).
Here, we will use the _last converged_ temperature to avoid
the effort of a linearisation.

\param container CONTAINER* (i)   container data
\param ele       ELEMENT*   (i)   the element
\param nelenod   INT        (i)   number of element nodes
\param gsurf     GSURF*     (i)   current geometry surface
\param gpc       DOUBLE[]   (i)   parametric Gauss point co-ordinates
\param shape     DOUBLE[]   (i)   shape functions at Gauss point
\param fac       DOUBLE     (i)   integration weight -> factor
\param cfac      DOUBLE     (i)   load curve factor
\param hflx      DOUBLE*    (o)   heat flux to be applied
\return void

\author bborn
\date 10/07
*/
void th3_load_heatconv(CONTAINER* container,
                       ELEMENT* ele,
                       INT nelenod,
                       GSURF* gsurf,
                       DOUBLE gpc[NDIM_THERM3],
                       DOUBLE shape[MAXNOD_THERM3],
                       DOUBLE fac,
                       DOUBLE cfac,
                       DOUBLE* hflx)
{
  DOUBLE sfac = 1.0;  /* spatial scale due to FUNCT */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_load_heatconv");
#endif

  /*--------------------------------------------------------------------*/
  /* scale spatially by FUNCTION */
  if (gsurf->neum->neum_type == neum_heatconvection)
  {
    /* do nothing hflux has final value */
    sfac = 1.0;
  }
  else if (gsurf->neum->neum_type == neum_heatconv_fct1)
  {
    /* physical position of Gauss point */
    DOUBLE gpx[NDIM_THERM3];
    th3_metr_refpos(ele, nelenod, gpc, gpx);
#if 0
    /* debug */
    printf("Ele %d : RST %g %g %g : XYZ %g %g %g\n", ele->Id,
           gpc[0], gpc[1], gpc[2],
           gpx[0], gpx[1], gpx[2]);
#endif
    /* get spatial scale */
    th3_load_heatconv_fct(0, gpx, &sfac);
#if 0
    /* debug */
    printf("Ele %d : RST %g %g %g : XYZ %g %g %g : Scale %g\n", 
           ele->Id,
           gpc[0], gpc[1], gpc[2],
           gpx[0], gpx[1], gpx[2],
           sfac);
#endif
  }
  else if (gsurf->neum->neum_type == neum_heatconv_fct2)
  {
    /* physical position of Gauss point */
    DOUBLE gpx[NDIM_THERM3];
    th3_metr_refpos(ele, nelenod, gpc, gpx);
    /* get spatial scale */
    th3_load_heatconv_fct(1, gpx, &sfac);
  }
  else
  {
    dserror("Wrong load type at this place");
  }

  /*--------------------------------------------------------------------*/
  /* determine heat flux value at Gauss point */
  {
    INT have_val0 = gsurf->neum->neum_onoff.a.iv[0];
    INT have_val1 = gsurf->neum->neum_onoff.a.iv[1];
    INT have_val2 = gsurf->neum->neum_onoff.a.iv[2];
    INT have_val3 = gsurf->neum->neum_onoff.a.iv[3];
    INT have_val4 = gsurf->neum->neum_onoff.a.iv[4];
    /* type 1 has 5 parameters */
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
      DOUBLE bastem = cfac * envtem;
      /* apply spatial distribution */
      if (bastem > t1)
      {
        bastem *= sfac;
      }
      /* check denominator */
      if (fabs(t2-t1) < EPS12)
      {
        dserror("Division by zero");
      }
      /* convection coefficient linearly interpolated */
      DOUBLE convcoeff = ((c2-c1)*bastem + c1*t2 - c2*t1)/(t2-t1);
      /* convective boundary flux */
      *hflx = convcoeff * (bastem - tem);
#if 0
      /* debug : print out */
      if ( ele->Id_loc == 713)
      {
        printf("Ele %d: T %g; T_infty %g; Coeff %g; Flux %g\n", 
               ele->Id_loc, tem, bastem, convcoeff, hflx);
      }
      if ( ele->Id_loc == 709)
      {
          printf("Ele %d: T %g; T_infty %g; Coeff %g; Flux %g\n", 
                 ele->Id_loc, tem, bastem,  convcoeff, hflx);
      }
#endif
    }
    /* type 2 has 2 parameters */
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
      *hflx = convcoeff * (envtem - tem);
    }
    else
    {
      dserror("Not enough parameters for heat convection");
    }
  }



  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of  th3_load_heatconv(...) */


/*======================================================================*/
/*!
\brief Spatial scale

\param ifct      INT        (i)   FUNCT index
\param gpx       DOUBLE[]   (i)   physical Gauss point co-ordinates
\param sfac      DOUBLE*    (o)   spatial scale factor
\return void

\author bborn
\date 10/07
*/
void th3_load_heatconv_fct(INT ifct,
                           DOUBLE gpx[NDIM_THERM3],
                           DOUBLE* sfac)
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_load_heatconv_fct");
#endif

  /*--------------------------------------------------------------------*/
  /* spatial scale */
  if (numfunct > ifct)
  {
    if (funct[ifct].functtyp == funct_none)
    {
      dserror("FUNCT%d was not defined in input file", ifct+1);
    }
    else if (funct[ifct].functtyp == funct_line_lin)
    {
      /* linear line function */
      FUNCT_LINE_LIN* fll = funct[ifct].typ.funct_line_lin;
      /* check if line is oriented parallel to an XYZ-axis */
      INT idir = -1;  /* index of relevant direction */
      INT ndir = 0;  /* counter of non-zero directions */
      INT idim;
      for (idim=0; idim<NDIM_THERM3; idim++)
      {
        DOUBLE dist = fabs(fll->x2[idim] - fll->x1[idim]);
        if (dist > EPS12)
        {
          ndir += 1;
          idir = idim;
        }
      }
      /* determine scale */
      if (ndir == 1)
      { 
        /* Linear interpolation along line/axis:
         *                y2 - y1
         *   y  =  y1  +  ------- (x - x1)
         *                x2 - x1
         *
         *        y1 (x2 - x1) + (y2 - y1) (x - x1)
         *      = ---------------------------------
         *                   x2 - x1
         */
        *sfac = ( fll->val1*(fll->x2[idir] - fll->x1[idir]) 
                  + (fll->val2 - fll->val1)*(gpx[idir] - fll->x1[idir]) ) 
                / (fll->x2[idir]-fll->x1[idir]);
      }
      else
      {
        dserror("Line must be oriented parallel to X-, Y- or Z-axis");
      }
    }
    else
    {
      dserror("Function type is wrong.");
    }
  }
  else
  {
    dserror("No function FUNCT%d defined.", ifct+1);
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
#endif  /*end of #ifdef D_THERM3 */
/*! @} (documentation module close)*/
