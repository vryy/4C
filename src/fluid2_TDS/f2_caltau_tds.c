/*!----------------------------------------------------------------------
\file
\brief stabilisation parameter for time dependent subscales

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/gammi/
            +49-(0)89-289-15235

</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2

#ifdef D_FLUID2_TDS
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2_TDS_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*!---------------------------------------------------------------------
\brief routine to evaluate the stabilisation paramter for time dependent
subscales

<pre>                                                         gammi 11/06

This routine evaluates the stabilisation parameters tau_M and tau_C for
the current and the last timestep according to a combination of the
Codina tau for time dependent subscales and the Franca version
implemented in the USEFEM part.
The tau_C is computed out of tau_M exactly like in the Codina reference.
tau_M is the one of Codina with one of the switches from the Franca tau_M.



see also: Barrenechea, G.R. and Valentin, F.: An unusual stabilized finite
          element method for a generalized Stokes problem. Numerische
          Mathematik, Vol. 92, pp. 652-677, 2002.
          http://www.lncc.br/~valentin/publication.htm
and:      Franca, L.P. and Valentin, F.: On an Improved Unusual Stabilized
          Finite Element Method for the Advective-Reactive-Diffusive
          Equation. Computer Methods in Applied Mechanics and Enginnering,
          Vol. 190, pp. 1785-1800, 2000.
          http://www.lncc.br/~valentin/publication.htm
and:      Codina, R: Time dependent subscales in the stabilized finite
          element approximation of incompressible flow problems.
	  (preprint)

</pre>

\param    *ele	        ELEMENT	   (i)    actual element
\param   **xyze	        DOUBLE 	   (i)    elemental coordinates
\param   *funct         DOUBLE 	   (-)    shape functions
\param  **deriv         DOUBLE 	   (-)    deriv. of shape funcs
\param    visc    	DOUBLE 	   (i)    viscosity


\return void

------------------------------------------------------------------------*/
void f2_get_time_dependent_sub_tau(ELEMENT *ele,
				   DOUBLE  **xyze,
				   DOUBLE   *funct,
				   DOUBLE  **deriv,
				   DOUBLE  **evelng,
				   DOUBLE  **eveln,
 				   DOUBLE    visc
    )
{

INT       i;              /* a counter                                  */
INT       iel;            /* number of element nodal points             */
DOUBLE    norm_p;
DOUBLE    mk=0,hk=0;
DOUBLE    re, xi2;

DOUBLE    velint[2],velint_old[2];

/* nach Moeglichkeit zu entfernen:*/
DOUBLE    xyz[2][MAXNOD];
DOUBLE    area;

DIS_TYP        typ;

FLUID_DYNAMIC   *fdyn;

#ifdef DEBUG
dstrc_enter("f2_get_time_dependent_sub_tau");
#endif

/*--------------------------------------------------- initialisation ---*/
fdyn    = alldyn[genprob.numff].fdyn;

typ = ele->distyp;
iel = ele->numnp;


/*------------ get shape functions and derivatives at element center ---*/
switch (ele->distyp)
{
  case quad4: case quad8: case quad9:
    f2_rec(funct,deriv,NULL,0.0,0.0,typ,2);
    break;
  case tri3: case tri6:
    f2_tri(funct,deriv,NULL,0.0,0.0,typ,2);
    break;
  default:
    dserror("typ unknown!");
} /* end switch(typ) */


/*--------------------------------- get velocities at element center ---*/
f2_veci(velint    ,funct,evelng,iel);
f2_veci(velint_old,funct,eveln ,iel);


/* tau for time dependent subscales --- a combination
 * of Franca and Codina tau for time dependent subscales */

/* be careful with this definition --- it's just a try......  */

/*--- get proper constant mk ---*/
switch(typ)
{
    case tri3:
	mk = 0.333333333333333333333;
	break;
    case tri6:
	mk = 0.083333333333333333333;
	break;
    case quad4:
	mk = 0.333333333333333333333;
	break;
    case quad8:
	mk = 0.083333333333333333333;
	break;
    case quad9:
	mk = 0.083333333333333333333;
	break;
    default: dserror("element type not implemented!");
}
   
/* square root of area for element length calculation */
/*--------------------- rewrite array of elemental coordinates ---*/
for(i=0; i<iel; i++)
{
    xyz[0][i] = xyze[0][i];
    xyz[1][i] = xyze[1][i];
}

/*-------------------------------- get area and element length ---*/
area = area_lin_2d(ele,xyz);

switch(typ)
{
    case tri3:
    case tri6:
	hk = sqrt(area);
	/* this is a rough estimate .... */
	break;
    case quad4:
    case quad8:
    case quad9:
	hk = sqrt(area);
	break;
    default: dserror("element type not implemented!");
}

/*---------------------------------------------------- get p-norm ---*/
norm_p = sqrt(DSQR(velint[0]) + DSQR(velint[1]));

re = mk * norm_p * hk / (2.0 * visc);  /* advective : viscous forces */

xi2 = DMAX(re,1.0);

fdyn->tau[0] = DSQR(hk) / (2 * visc/mk + (4.0 * visc/mk) * xi2);

fdyn->tau[2]=  DSQR(hk) /(fdyn->tau[0]*2./mk);

/*---------------------------------------------------- get p-norm ---*/
norm_p = sqrt(DSQR(velint_old[0]) + DSQR(velint_old[1]));

re = mk * norm_p * hk / (2.0 * visc);  /* advective : viscous forces */

xi2 = DMAX(re,1.0);

fdyn->tau_old[0] = DSQR(hk) / (2 * visc/mk + (4.0 * visc/mk) * xi2);

fdyn->tau_old[2]=  DSQR(hk) /(fdyn->tau_old[0]*2./mk);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
}
#endif /*D_FLUID2_TDS*/
#endif /*D_FLUID2*/
/*! @} (documentation module close)*/
