/*!----------------------------------------------------------------------
\file
\brief Calculate stabilisation parameter

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0711 - 685-6575
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID2TU
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
#include "fluid2_tu.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
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
\brief routine to calculate stabilisation parameter per element at center

<pre>                                                         he  02/03

</pre>
\param  *ele     ELEMENT	       (i)   actual element
\param  *elev    ELEMENT	       (i)   actual element for vel.
\param  *funct   DOUBLE 	       (-)   shape functions
\param **deriv   DOUBLE 	       (-)   deriv. of shape funcs
\param **deriv2  DOUBLE 	       (-)   2nd deriv. of sh. funcs
\param **evel    DOUBLE 	       (i)   element velocities
\param  *eddyg   DOUBLE 	       (i)   eddy-viscosity
\param  *velint  DOUBLE 	       (-)   vel. at integr. point
\param  *velint_dc  DOUBLE 	       (-)   vel. at integr. point for D.C.
\param  *kapomen    DOUBLE 	       (i)   kapomen at nodes
\param **xjm        DOUBLE 	       (-)   jacobian  matrix
\param **xzye       DOUBLE           (-)   nodal coordinates
\param **derxy      DOUBLE 	       (-)   global deriv.
\param  *kapomederxy   DOUBLE 	 (-)   kapome global deriv.
\param **cutp         DOUBLE 	       (-)   array
\return void

------------------------------------------------------------------------*/
void f2_calelesize_tu_1(
	           ELEMENT         *ele,
		     ELEMENT         *elev,
	           DOUBLE          *funct,
	           DOUBLE         **deriv,
	           DOUBLE         **deriv2,
                 DOUBLE         **evel,
                 DOUBLE          *eddyg,
                 DOUBLE          *velint,
                 DOUBLE          *velint_dc,
                 DOUBLE          *kapomen,
	           DOUBLE         **xjm,
	           DOUBLE         **xyze,
	           DOUBLE         **derxy,
                 DOUBLE          *kapomederxy,
                 DOUBLE         **cutp
                 )
{
INT     actmat;         /* number of actual material		            */
INT     iel;            /* number of nodes of actual element            */
INT     icode=2;        /* flag for eveluation of shape functions       */
INT     ihoel=0;        /* flag for higher order elements               */
DOUBLE  e1,e2;          /* natural coordinates of inegration point      */
DOUBLE  visc;           /* fluid viscosity                              */
DOUBLE  eddyint;
DOUBLE  facr,facs,det;
DOUBLE  gcoor[2];       /* global coordinates                           */
DIS_TYP typ;
FLUID_DYNAMIC  *fdyn;
FLUID_DATA    *data;

#ifdef DEBUG
dstrc_enter("f2_calelesize_tu_1");
#endif

/*---------------------------------------------------------- initialise */
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;
actmat = ele->mat-1;
visc   = mat[actmat].m.fluid->viscosity;
iel    = ele->numnp;
typ    = ele->distyp;

/*------ get values of integration parameters, shape functions and their
         derivatives ---------------------------------------------------*/
   switch(typ)
   {
   case quad4: case quad8: case quad9:    /* --> quad - element */
      icode= 3;
      ihoel= 1;
      e1   = data->qxg[0][0];
      facr = data->qwgt[0][0];
      e2   = data->qxg[0][0];
      facs = data->qwgt[0][0];
      f2_rec(funct,deriv,deriv2,e1,e2,typ,icode);
   break;
   case tri3: case tri6:       /* --> tri - element */
     if (iel>3)
     {
      icode   = 3;
      ihoel   = 1;
     }
      e1   = data->txgr[0][0];
      facr = data->twgt[0][0];
      e2   = data->txgs[0][0];
      facs = ONE;
      f2_tri(funct,deriv,deriv2,e1,e2,typ,icode);
   break;
   default:
      dserror("typ unknown!\n");
   } /*end switch(typ) */

/*-------------------------------------------- compute Jacobian matrix */
   f2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);

/*------------------------------------------- compute global derivates */
   f2_gder(derxy,deriv,xjm,det,iel);

/*--------------------------------- get eddy-visc. at center of element */
   f2_eddyi(&eddyint,funct,eddyg,iel);

/*----------------------------------- get velocity at center of element */
   f2_veci(velint,funct,evel,iel);

/*------------------- calculate stabilisation parameter for DISC. CAPT. */
   f2_kapomeder(kapomederxy,derxy,kapomen,iel);
   f2_vel_dc_1(velint,velint_dc,kapomederxy);

/*---------------- get streamlenght for elementlenght for DISC. CAPT.  */
   f2_gcoor(xyze,funct,iel,gcoor);
   f2_calstrlen_tu(velint_dc,ele,gcoor,cutp,typ);

/*----------------------------------- calculate stabilisation parameter */
   f2_calstabpar_tu_1(ele,elev,eddyint,velint,velint_dc,visc);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calelesize */

#endif
