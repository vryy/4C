/*!----------------------------------------------------------------------
\file
\brief Calculate stabilisation parameter for second fluid implementation

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6572
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID3 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"
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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!---------------------------------------------------------------------
\brief routine to calculate element size and stabilisation parameter

<pre>                                                        chfoe 10/04
This routine evaluates the stabilisation parameter for USFEM stabilised 
fluid elements in 3D. Different ways of calculating the tau are 
implemented. 
Here everything is done once at the element center (center of its local
coordinates) and used for the entire element.
</pre>
\param  *ele     ELEMENT	       (i)   actual element
\param **xzye    DOUBLE                (i)   nodal coordinates
\param  *funct   DOUBLE 	       (-)   shape functions
\param **deriv   DOUBLE 	       (-)   deriv. of shape funcs
\param **deriv2  DOUBLE 	       (-)   2nd deriv. of sh. funcs
\param **xjm     DOUBLE 	       (-)   jacobian matrix
\param  *visc     DOUBLE 	       (-)   viscosity
\return void             

------------------------------------------------------------------------*/
void f3_caltau(			     
	       ELEMENT         *ele, 
	       DOUBLE         **xyze,
	       DOUBLE          *funct,  
	       DOUBLE         **deriv,  
	       DOUBLE         **derxy,              
	       DOUBLE         **xjm,  
	       DOUBLE         **evelng,
               DOUBLE         **wa1,
               DOUBLE           visc  
              )
{

INT            i, iel;

DOUBLE         hk;
DOUBLE         velint[3];
DOUBLE         det, vol;
DOUBLE         e1, e2, e3, facr, facs, fact;
DOUBLE         mk, norm_p, pe, re, xi1, xi2, timefac;
DOUBLE         val, strle;
DOUBLE         velino[3]; /* normed velocity at element centre */

FLUID_DYNAMIC *fdyn;
DIS_TYP        typ;
FLUID_DATA    *data;

#ifdef DEBUG 
dstrc_enter("f3_caltau");
#endif		

/*------------------------------------------------------- initialise ---*/
fdyn    = alldyn[genprob.numff].fdyn;
data    = fdyn->data;
timefac = fdyn->thsl;

/*----------------------------------- set parameter for this element ---*/
iel    = ele->numnp;
typ    = ele->distyp;

/*---------------------- shape functions and derivs at element center --*/
switch(typ)
{
case hex8: case hex20: case hex27:   /* --> hex - element */
   e1   = data->qxg[0][0];
   facr = data->qwgt[0][0];
   e2   = data->qxg[0][0];
   facs = data->qwgt[0][0];
   e3   = data->qxg[0][0];
   fact = data->qwgt[0][0];
   f3_hex(funct,deriv,NULL,e1,e2,e3,typ,2);
break;
case tet4: case tet10:   /* --> tet - element */
   e1   = data->txgr[0][0];
   facr = data->twgt[0][0];
   e2   = data->txgs[0][0];
   facs = ONE;
   e3   = data->txgs[0][0];
   fact = ONE;
   f3_tet(funct,deriv,NULL,e1,e2,e3,typ,2);
break;
default:
   dserror("type unknown!\n");
} /*end switch(typ) */

/*------------------------------- get element type constant for tau ---*/
switch(typ)
   {
   case tet4:
   case hex8:
      mk = 0.333333333333333333333;
      mk = 0.333333333333333333333;
   break;
   case hex20:
   case hex27:
   case tet10:
      mk = 0.083333333333333333333;
   break;
   default: dserror("type unknown!\n");
   }
/*--------------------------------- get velocities at element center ---*/
f3_veci(velint,funct,evelng,iel);

/*----------------------------- get Jacobian matrix and determinante ---*/
f3_jaco(xyze,deriv,xjm,&det,ele,iel);
vol=facr*facs*fact*det;

/*----------------------------------------------- get element length ---*/
hk = pow((SIX*vol/PI),(ONE/THREE))/sqrt(THREE);

/*------------------------------------------------- get streamlength ---*/
f3_gder(derxy,deriv,xjm,wa1,det,iel);
val = ZERO;
/* get p-norm */
norm_p=sqrt( velint[0]*velint[0]
           + velint[1]*velint[1]
           + velint[2]*velint[2]);
if(norm_p>=EPS6)
{
   velino[0] = velint[0]/norm_p;
   velino[1] = velint[1]/norm_p;
   velino[2] = velint[2]/norm_p;
}
else
{
   velino[0] = ONE;
   velino[1] = ZERO;
   velino[2] = ZERO;
}
for (i=0;i<iel;i++) /* loop element nodes */
{
   val += FABS(velino[0]*derxy[0][i] \
              +velino[1]*derxy[1][i] \
              +velino[2]*derxy[2][i]);
} /* end of loop over elements */
strle=TWO/val;

/*--------------------------------------------------------- get tau_Mu ---*/
pe = 4.0 * timefac * visc / (mk * DSQR(strle)); /* viscous : reactive forces */
re = mk * norm_p * strle / (2.0 * visc);       /* advective : viscous forces */

xi1 = DMAX(pe,1.0);
/*xi1 = pe; */
xi2 = DMAX(re,1.0);

fdyn->tau[0] = DSQR(strle) / (DSQR(strle)*xi1+(4.0*timefac*visc/mk)*xi2);

/*--------------------------------------------------------- get tau_Mp ---*/
pe = 4.0 * timefac * visc / (mk * DSQR(hk)); /* viscous : reactive forces */
re = mk * norm_p * hk / (2.0 * visc);       /* advective : viscous forces */

xi1 = DMAX(pe,1.0);
/*xi1 = pe; */
xi2 = DMAX(re,1.0);

fdyn->tau[1] = DSQR(hk) / (DSQR(hk) * xi1 + (4.0*timefac * visc/mk) * xi2);

/*---------------------------------------------------------- get tau_C ---*/
fdyn->tau[2] = norm_p * hk * 0.5 * xi2;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_caltau */

#endif
