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

<pre>                                                         he  12/02

</pre>
\param  *ele     ELEMENT	       (i)   actual element
\param  *elev    ELEMENT	       (i)   actual element for vel.
\param  *data    FLUID_DATA	       (i)
\param  *funct   DOUBLE 	       (-)   shape functions
\param **deriv   DOUBLE 	       (-)   deriv. of shape funcs
\param **deriv2  DOUBLE 	       (-)   2nd deriv. of shape funcs
\param **evel    DOUBLE 	       (i)   element velocities
\param  *eddyg    DOUBLE 	       (i)   eddy-viscosity
\param  *velint   DOUBLE 	       (-)   vel. at integr. point
\param  *velint_dc DOUBLE 	       (-)   vel. for D.C. at integr. point
\param  *kapepsn   DOUBLE 	       (i)   kapepsn at nodes
\param **xjm        DOUBLE 	       (-)   jacobian  matrix
\param **xzye       DOUBLE           (-)   nodal coordinates
\param **derxy      DOUBLE 	       (-)   global deriv.
\param  *kapepsderxy  DOUBLE 	       (-)   kapeps global deriv.
\param **cutp         DOUBLE 	       (-)   array
\return void

------------------------------------------------------------------------*/
void f2_calelesize_tu(
	           ELEMENT         *ele,
		     ELEMENT         *elev,
	           DOUBLE          *funct,
	           DOUBLE         **deriv,
	           DOUBLE         **deriv2,
                 DOUBLE         **evel,
                 DOUBLE          *eddyg,
                 DOUBLE          *velint,
                 DOUBLE          *velint_dc,
                 DOUBLE          *kapepsn,
	           DOUBLE         **xjm,
	           DOUBLE         **xyze,
	           DOUBLE         **derxy,
                 DOUBLE          *kapepsderxy,
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
FLUID_DYNAMIC *fdyn;
FLUID_DATA    *data;

#ifdef DEBUG
dstrc_enter("f2_calelesize_tu");
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

/*-------------------------------- get eddy-visc. at center of element */
   f2_eddyi(&eddyint,funct,eddyg,iel);

/*---------------------------------- get velocity at center of element */
   f2_veci(velint,funct,evel,iel);

/*------------------- calculate stabilisation parameter for DISC. CAPT. */
   f2_kapepsder(kapepsderxy,derxy,kapepsn,iel);
   f2_vel_dc(velint,velint_dc,kapepsderxy);

/*---------------- get streamlenght for elementlenght for DISC. CAPT.  */
   f2_gcoor(xyze,funct,iel,gcoor);
   f2_calstrlen_tu(velint_dc,ele,gcoor,cutp,typ);

/*----------------------------------- calculate stabilisation parameter */
   f2_calstabpar_tu(ele,elev,eddyint,velint,velint_dc,visc);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calelesize */

/*!---------------------------------------------------------------------
\brief routine to calculate streamlength

<pre>                                                         he 03/03

in this routine the the streamlength, used for calculation of
stabilisation parameter is calculated. for higher order element this
is only an approximation, since the boundaries are assumed to be
straight.

</pre>
\param  *velint_dc    DOUBLE   (i)    velocities at integr. point for D.C.
\param  *ele 	    ELEMENT  (i)    actual element
\param  *gcoor        DOUBLE   (i)    global coord. of INT. point
\param **cutp         DOUBLE   (-)    cutting points
\param   typ	    DIS_TYP     (i)    flag for element type
\return void

------------------------------------------------------------------------*/
void f2_calstrlen_tu(
		         DOUBLE   *velint_dc,
		         ELEMENT  *ele,
                     DOUBLE   *gcoor,
		         DOUBLE  **cutp,
		         DIS_TYP       typ
                     )
{
INT     nodcut=-1;
INT     nodmax;
INT     inod;
DOUBLE dl,dx,dy,dxh,dyh;
DOUBLE dsub,dval;

#ifdef DEBUG
dstrc_enter("f2_calstrlen_tu");
#endif

dval = FABS(velint_dc[0])+FABS(velint_dc[1]);
if (dval == ZERO)  /* no flow at this point - take some arbitr. measure for streamlength */
{
   dx = ele->node[2]->x[0] - ele->node[0]->x[0];
   dy = ele->node[2]->x[1] - ele->node[0]->x[1];
   goto calc2;
} /* enidf (dval == ZERO) */

/*----------------------------------------------------------------------*
   streamlength is calculated via cutting points of velocity vector
   with straight boundaries
*/
switch(typ)
{
case quad4: case quad8: case quad9: /* max number of nodes for quad: 4 --> C-numbering nodmax = 4-1 */
   nodmax = 3;
break;
case tri3: case tri6:  /* max number of nodes for tri: 3 --> C-numbering nodmax = 3-1 */
   nodmax = 2;
break;
default:
   dserror("typ unknown!\n");
} /* end switch(typ) */
 /*------------------------------------------------- get cutting points */
for (inod=0;inod<nodmax;inod++)
{
   dxh = ele->node[inod+1]->x[0] - ele->node[inod]->x[0];
   dyh = ele->node[inod+1]->x[1] - ele->node[inod]->x[1];
   dsub = dxh*velint_dc[1]-dyh*velint_dc[0];
   if (dsub==ZERO)  /* check for parallel vectors */
      continue;
   dl = ((ele->node[inod]->x[1]-gcoor[1])*velint_dc[0] -	\
	 (ele->node[inod]->x[0]-gcoor[0])*velint_dc[1])/dsub;
   if (dl>=ZERO && dl<=ONE)
   {
      nodcut++;
      cutp[0][nodcut]=ele->node[inod]->x[0]+dl*dxh;
      cutp[1][nodcut]=ele->node[inod]->x[1]+dl*dyh;
      if (nodcut==1)
	 goto calc1;
   } /* endif (dl>=ZERO && dl<=ONE) */
} /* end loop over inod */
/*------------------------------------------------- the last boundary */
dxh = ele->node[0]->x[0]-ele->node[nodmax]->x[0];
dyh = ele->node[0]->x[1]-ele->node[nodmax]->x[1];

dsub = dxh*velint_dc[1] - dyh*velint_dc[0];
if (dsub==ZERO)
   dserror("Couldn't find two cutting points!\n");
dl = ((ele->node[nodmax]->x[1]-gcoor[1])*velint_dc[0] -	\
      (ele->node[nodmax]->x[0]-gcoor[0])*velint_dc[1])/dsub;
if (dl>=ZERO && dl <= ONE)
{
   nodcut++;
   cutp[0][nodcut]=ele->node[nodmax]->x[0]+dl*dxh;
   cutp[1][nodcut]=ele->node[nodmax]->x[1]+dl*dyh;
   if(nodcut==1)
      goto calc1;
} /* endif  (dl>=ZERO && dl <= ONE) */

dserror("Couldn't find two cutting points!\n");

calc1:
dx = cutp[0][1]-cutp[0][0];
dy = cutp[1][1]-cutp[1][0];

calc2:
ele->e.f2_tu->strom_dc = sqrt(dx*dx+dy*dy);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calstrlen_tu */

#endif
