/*!----------------------------------------------------------------------
\file
\brief Calculate stabilisation parameter

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
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

/*!---------------------------------------------------------------------
\brief routine to calculate stabilisation parameter per element at center

<pre>                                                         he  02/03

</pre>
\param  *ele     ELEMENT	       (i)   actual element
\param  *elev    ELEMENT	       (i)   actual element for vel.
\param  *dynvar  FLUID_DYN_CALC      (i/o)
\param  *data    FLUID_DATA	       (i)
\param  *funct   double 	       (-)   shape functions
\param **deriv   double 	       (-)   deriv. of shape funcs
\param **deriv2  double 	       (-)   2nd deriv. of sh. funcs
\param **evel    double 	       (i)   element velocities
\param  *eddyg   double 	       (i)   eddy-viscosity
\param  *velint  double 	       (-)   vel. at integr. point
\param  *velint_dc  double 	       (-)   vel. at integr. point for D.C.
\param  *kapomen    double 	       (i)   kapomen at nodes 
\param **xjm        double 	       (-)   jacobian  matrix
\param **xzye       double           (-)   nodal coordinates
\param **derxy      double 	       (-)   global deriv. 
\param  *kapomederxy   double 	 (-)   kapome global deriv. 
\param **cutp         double 	       (-)   array 
\return void             

------------------------------------------------------------------------*/
void f2_calelesize_tu_1(			     
	           ELEMENT         *ele,    
		     ELEMENT         *elev,    
                 FLUID_DYN_CALC  *dynvar,
		     FLUID_DATA      *data, 
	           double          *funct,  
	           double         **deriv,  
	           double         **deriv2,  		 
                 double         **evel, 
                 double          *eddyg, 
                 double          *velint, 
                 double          *velint_dc, 
                 double          *kapomen, 
	           double         **xjm,     
	           double         **xyze,     
	           double         **derxy,   
                 double          *kapomederxy,  
                 double         **cutp
                 )
{
int     ntyp;           /* element type (TRI or QUAD)  		      */
int     actmat;         /* number of actual material		            */
int     iel;            /* number of nodes of actual element            */
int     icode=2;        /* flag for eveluation of shape functions       */     
int     ihoel=0;        /* flag for higher order elements               */
double  e1,e2;          /* natural coordinates of inegration point      */
double  visc;           /* fluid viscosity                              */
double  eddyint;
double  facr,facs,det;
double  gcoor[2];       /* global coordinates                           */
DIS_TYP typ;

#ifdef DEBUG 
dstrc_enter("f2_calelesize_tu_1");
#endif		

/*---------------------------------------------------------- initialise */
actmat = ele->mat-1;
visc   = mat[actmat].m.fluid->viscosity;
ntyp   = ele->e.f2_tu->ntyp;
iel    = ele->numnp;
typ    = ele->distyp;

/*------ get values of integration parameters, shape functions and their
         derivatives ---------------------------------------------------*/
   switch(ntyp)
   {
   case 1:    /* --> quad - element */
      icode= 3;
      ihoel= 1;
      e1   = data->qxg[0][0];
      facr = data->qwgt[0][0];
      e2   = data->qxg[0][0];
      facs = data->qwgt[0][0];
      f2_rec(funct,deriv,deriv2,e1,e2,typ,icode);
   break;
   case 2:       /* --> tri - element */              
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
      dserror("ntyp unknown!\n");      
   } /*end switch(ntyp) */
     
/*-------------------------------------------- compute Jacobian matrix */
   f2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);

/*------------------------------------------- compute global derivates */
   f2_gder(derxy,deriv,xjm,det,iel);

/*--------------------------------- get eddy-visc. at center of element */               
   f2_eddyi(&eddyint,funct,eddyg,iel);

/*----------------------------------- get velocity at center of element */               
   f2_veli(velint,funct,evel,iel);    

/*------------------- calculate stabilisation parameter for DISC. CAPT. */
   f2_kapomeder(kapomederxy,derxy,kapomen,iel);
   f2_vel_dc_1(dynvar,velint,velint_dc,kapomederxy);

/*---------------- get streamlenght for elementlenght for DISC. CAPT.  */
   f2_gcoor(xyze,funct,iel,gcoor);
   f2_calstrlen_tu(velint_dc,ele,gcoor,cutp,ntyp);       

/*----------------------------------- calculate stabilisation parameter */               
   f2_calstabpar_tu_1(ele,elev,dynvar,eddyint,velint,velint_dc,visc); 

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calelesize */

#endif
