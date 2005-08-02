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
#ifdef D_FLUID3_F
#include "../headers/standardtypes.h"
#include "../fluid3_fast/f3f_prototypes.h"
#include "../fluid3/fluid3.h"
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
void f3fcaltau(			     
	       ELEMENT         *ele[LOOPL], 
	       DOUBLE          *elecord,
	       DOUBLE          *funct,  
	       DOUBLE          *deriv,  
	       DOUBLE          *derxy, 
               DOUBLE          *velint,
	       DOUBLE          *xjm,  
	       DOUBLE          *evelng,
               DOUBLE          *tau,
               DOUBLE          *wa1,
               INT              sizevec[6]
              )
{

INT            i,l;
INT            actmat;             /* number of actual material */

DOUBLE         visc;               /* fluid viscosity */
DOUBLE         hk[LOOPL];
DOUBLE         det[LOOPL];
DOUBLE         vol[LOOPL];
DOUBLE         val[LOOPL];
DOUBLE         re[LOOPL];
DOUBLE         pe[LOOPL];
DOUBLE         norm_p[LOOPL];
DOUBLE         strle[LOOPL];
DOUBLE         xi1[LOOPL];
DOUBLE         xi2[LOOPL];
DOUBLE         e1, e2, e3;
DOUBLE         facr=0.0, facs=0.0, fact=0.0;
DOUBLE         mk=0.0, timefac;
DOUBLE         velino[3]; /* normed velocity at element centre */

FLUID_DYNAMIC *fdyn;
DIS_TYP        typ;
FLUID_DATA    *data;

/*fortran variables */
INT typint;
INT icode;

#ifdef DEBUG 
dstrc_enter("f3fcaltau");
#endif		

/*------------------------------------------------------- initialise ---*/
fdyn    = alldyn[genprob.numff].fdyn;
data    = fdyn->data;
timefac = fdyn->thsl;
actmat  = ele[0]->mat-1;
visc    = mat[actmat].m.fluid->viscosity;

icode = 2;

for(l=0;l<sizevec[4];l++)
{
  vol[l]   = ZERO;
  strle[l] = ZERO;
}/*loop*/

/*------------------------------------ set general element parameter ---*/
typ    = ele[0]->distyp;

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

   if(typ==hex8) typint=8;
   else if(typ==hex20) typint=20;
   else if(typ==hex27) typint=27;

   f3fhex(funct, deriv, NULL, &e1, &e2, &e3, &typint,
          &icode, sizevec);
break;
case tet4: case tet10:   /* --> tet - element */
   e1   = data->txgr[0][0];
   facr = data->twgt[0][0];
   e2   = data->txgs[0][0];
   facs = ONE;
   e3   = data->txgs[0][0];
   fact = ONE;

   if(typ==tet4) typint=4;
   else if(typ==tet10) typint=10;

   f3ftet(funct, deriv, NULL, &e1, &e2, &e3, &typint,
       &icode, sizevec);
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

/*----------------------------- get Jacobian matrix and determinante ---*/
f3fjaco(funct, deriv, xjm, det, elecord, sizevec);

/*------------------------------------ get volume and element length ---*/
for(l=0;l<sizevec[4];l++)
{
   vol[l]  = facr*facs*fact*det[l];
   hk[l]   = pow((SIX*vol[l]/PI),(ONE/THREE))/sqrt(THREE);
} /*loop*/

/*--------------------------------- get velocities at element center ---*/
f3fveli(velint, funct, evelng, sizevec);

/*------------------------------------------------- get streamlength ---*/
f3fgder(derxy, deriv, xjm, wa1, det, sizevec);

for(l=0;l<sizevec[4];l++) /* loop elements of this set */
{
   val[l] = ZERO;
   /* get p-norm */
   norm_p[l] = sqrt( velint[l]*velint[l]
                   + velint[LOOPL+l]*velint[LOOPL+l]
                   + velint[2*LOOPL+l]*velint[2*LOOPL+l]);
   if( norm_p[l] >= EPS6 )
   {
      velino[0] = velint[l]/norm_p[l];
      velino[1] = velint[LOOPL+l]/norm_p[l];
      velino[2] = velint[2*LOOPL+l]/norm_p[l];
   }
   else
   {
      velino[0] = ONE;
      velino[1] = ZERO;
      velino[2] = ZERO;
   }
   for (i=0;i<sizevec[1];i++) /* loop element nodes */
   {
      val[l] += FABS(velino[0]*derxy[                   i*LOOPL+l]
                    +velino[1]*derxy[  sizevec[0]*LOOPL+LOOPL*i+l]
                    +velino[2]*derxy[2*sizevec[0]*LOOPL+LOOPL*i+l]);
   } /* end of loop over nodes */

   strle[l] = TWO/val[l];

   /*------------------------------------------------------ get tau_Mu ---*/
   /* viscous : reactive forces */
   pe[l] = 4.0 * timefac * visc / (mk * DSQR(strle[l]));
   /* advective : viscous forces */
   re[l] = mk * norm_p[l] * strle[l] / (2.0 * visc);       

   xi1[l] = DMAX(pe[l],1.0);
   xi2[l] = DMAX(re[l],1.0);

   /* stab parameter for convective stabilisation */
   tau[l] = DSQR(strle[l]) / 
           (DSQR(strle[l])*xi1[l]+(4.0*timefac*visc/mk)*xi2[l]);

   /*------------------------------------------------------ get tau_Mp ---*/
   /* viscous : reactive forces */
   pe[l] = 4.0 * timefac * visc / (mk * DSQR(hk[l]));
   /* advective : viscous forces */
   re[l] = mk * norm_p[l] * hk[l] / (2.0 * visc);

   xi1[l] = DMAX(pe[l],1.0);
   xi2[l] = DMAX(re[l],1.0);

   /* stab parameter for non-convective stabilisation terms */
   tau[sizevec[3]+l] = DSQR(hk[l]) 
                  / (DSQR(hk[l])*xi1[l]+(4.0*timefac*visc/mk)*xi2[l]);

   /*------------------------------------------------------- get tau_C ---*/
   tau[2*sizevec[3]+l] = norm_p[l] * hk[l] * 0.5 * xi2[l];

} /* loop */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3fcaltau */

#endif
