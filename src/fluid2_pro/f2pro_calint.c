/*!----------------------------------------------------------------------
\file
\brief integration loop for one fluid2 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2_PRO 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2_PRO 
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2pro_prototypes.h"
#include "fluid2pro.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 | integration loop for one fluid element                               |
 |                                                                      |
 |                                                           genk 03/02 |
 *----------------------------------------------------------------------*/
/*!---------------------------------------------------------------------
\brief integration loop for one fluid2pro element

<pre>                                                         genk 04/02
                                                             basol 11/02
In this routine the element stiffness matrix, iteration-RHS and
time-RHS for one fluid2pro element is calculated
      
</pre>
\param  *data      FLUID_DATA	   (i)	  integration data
\param  *elev	   ELEMENT	   (i)    actual element for velocity
\param  *elep	   ELEMENT	   (i)    actual element for pressure
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     DOUBLE	   (o)    element stiffness matrix
\param **emass     DOUBLE	   (o)    element mass matrix
\param **gradopr    DOUBLE	   (o)    gradient operator
\param  *etforce   DOUBLE	   (o)    element time force vector part1
\param  *eiforce   DOUBLE	   (o)    element iter force vector
\param  *funct     DOUBLE	   (-)    natural shape functions for velocity
\param  **xyze     DOUBLE	   (i)    coordinates of the element
\param  *functpr   DOUBLE	   (-)    natural shape functions for pressure
\param **deriv     DOUBLE	   (-)	  deriv. of velocity nat. shape funcs
\param **derivpr   DOUBLE	   (-)	  deriv. of pressure nat. shape funcs
\param **xjm	   DOUBLE	   (-)    jacobian matrix
\param **derxy     DOUBLE	   (-)	  global derivatives for velocity shape fnc.
\param **derxypr   DOUBLE	   (-)	  global derivatives for pressure shape fnc.
\param **eveln     DOUBLE	   (i)    ele vel. at time n
\param  *epren     DOUBLE	   (-)    ele pres. at time n
\param  *velint    DOUBLE	   (-)    vel at integration point
\param  *covint    DOUBLE	   (-)    conv. vel. at integr. point
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\param  *pderxy    DOUBLE	   (-)    global pres. derivatives
\param **wa1       DOUBLE	   (-)    working array
\param  *dirich    DOUBLE	   (-)    dirichlet vector
\param **deriv2    DOUBLE	   (-)    second der. of the nat. shape funcs. 
\param  *dirich_onoff DOUBLE	   (-)    element dirich_onoff values
\return void                                                   

------------------------------------------------------------------------*/
void f2pro_calint(
               FLUID_DATA      *data,     
	       ELEMENT         *elev,
	       ELEMENT         *elep,     
	       FLUID_DYN_CALC  *dynvar, 
               DOUBLE         **estif,   
	       DOUBLE         **emass,
	       DOUBLE         **gradopr,
	       DOUBLE         *etforce,
	       DOUBLE         *eiforce,
	       DOUBLE         **xyze, 
	       DOUBLE         *funct,
	       DOUBLE         *functpr,   
	       DOUBLE         **deriv,
	       DOUBLE         **derivpr,   
	       DOUBLE         **xjm,     
	       DOUBLE         **derxy,
	       DOUBLE         **derxypr,   
	       DOUBLE         **eveln,   
	       DOUBLE         *epren,   
	       DOUBLE         *velint,  
	       DOUBLE         *covint,  
	       DOUBLE         **vderxy,  
	       DOUBLE         *pderxy,
	       DOUBLE         **wa1,  
	       DOUBLE         *dirich,
	       DOUBLE         **deriv2,  
	       INT             *dirich_onoff
               )
{ 
INT       i,j;          /* simply a counter                               */
INT       iel;          /* number of nodes                                */
INT       ielp;         /* number of nodes for pressure element           */
INT       ntyp;         /* element type: 1 - quad; 2 - tri                */
INT       intc;         /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c              */
INT       nir,nis;      /* number of integration nodesin r,s direction    */
INT       actmat;       /* material number of the element                 */
INT       icode=2;      /* flag for eveluation of shape functions         */     
INT       lr, ls;       /* counter for integration                        */
DOUBLE    dens;         /* density                                        */
DOUBLE    visc;         /* viscosity                                      */
DOUBLE    fac;
DOUBLE    facr, facs;   /* integration weights                            */
DOUBLE    det;          /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;        /* natural coordinates of integr. point           */
DOUBLE    preint;       /* pressure at integration point                  */
DIS_TYP   typ1;	        /* element type 1 (quad9)                         */
DIS_TYP   typ2;	        /* element type 2 (quad4)                         */
DISMODE   mode;         /* element discretization mode q2q1 or none       */

#ifdef DEBUG 
dstrc_enter("f2pro_calint");
#endif

/*----------------------------------------------------- initialisation */
iel =elev->numnp;
ielp=elep->numnp;                             
actmat=elev->mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity;
ntyp = elev->e.f2pro->ntyp; 
mode = elev->e.f2pro->dm;

if (mode==dm_q2q1)
{
   typ1= quad9;
   typ2= quad4;
}
else
{
   dserror ("unknown discretization mode");
}/*end of if(mode==dm_q2q1)*/

/*--------------------------------- an additional plausibility check */
dsassert(typ1==elev->distyp,"wrong element typ!\n");
dsassert(typ2==elep->distyp,"wrong element typ!\n");

/* initialise integration */
nir = elev->e.f2pro->nGP[0];
nis = elev->e.f2pro->nGP[1];

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{    
   for (ls=0;ls<nis;ls++)
   {
      e1   = data->qxg[lr][nir-1];
      facr = data->qwgt[lr][nir-1];
      e2   = data->qxg[ls][nis-1];
      facs = data->qwgt[ls][nis-1];
      f2_rec(funct,deriv,deriv2,e1,e2,typ1,icode);
      /*-------------------------- we need the pressure shape fnct. too */
      f2_rec(functpr,derivpr,deriv2,e1,e2,typ2,icode); 
      /*--------------------------------------- compute Jacobian matrix */
      f2_jaco(xyze,funct,deriv,xjm,&det,iel,elev);
      fac = facr*facs*det;
      /*---------- compute global derivates for velocity shape function */
      f2_gder(derxy,deriv,xjm,det,iel);
      /*----------------------- get the velocities at integration point */
      f2_veli(velint,funct,eveln,iel);
      /*----------- get velocity (n,i) derivatives at integration point */
      f2_vder(vderxy,derxy,eveln,iel);
      /*----------------------------------------------------------------*/
      if (dynvar->pro_calmat==1)
      {
      /*---------------------------- compute standard Galerkin matrices */
      /*-------------- standard Galerkin mass matrix is stored in emass */
      /*--------------------- element stiffness matrix is stored in kvv */
         if (dynvar->pro_kvv==1) 
	 f2pro_calkvv(estif,derxy,fac,visc,dynvar->dta,iel);
         /*-- calculate the Balancing Diffusivity Tensor and add to the */
         /*-------------------------------------------------- estif term*/
         f2pro_calbdt(elev,dynvar,estif,velint,derxy,vderxy,funct,fac,visc,iel);
         /*-------------------------------- get the gradient operator C */
         if (dynvar->pro_gra==1)
	 f2pro_gradopr(gradopr,derxy,functpr,fac,ielp,iel); 
         /*-------------------------- calculate the element mass matrix */
         if (dynvar->pro_mvv==1)
	 f2_calmvv(emass,funct,fac,iel);         
         /*---- end of Galerkin matrices calculation for left hand side * 
	  |                           (M+dt*K)u~(n+1)                   | 
         /*-------------------------------------------------------------*/
      }/*end of if(cal_mat==1)*/
      /*----------------------------------------------------------------*/
      if (dynvar->pro_calrhs==1) 
      {
         /*---------------------- get pressure (n) at integration point */
         f2_prei(&preint,functpr,epren,ielp);
         /*-------------------- get velocities (n) at integration point */
	 f2_veli(velint,funct,eveln,iel);
         /*---------------------------- get global velocity derivatives */
         f2_vder(vderxy,derxy,eveln,iel);
         /*--------- get convective velocities (n) at integration point */
         f2_covi(vderxy,velint,covint); 
         /*------------ calculate galerkin part of "Time-RHS" (vel-dofs)*/
         f2pro_calgaltfv(dynvar,etforce,eiforce,velint,covint,vderxy,
	                 funct,derxy,preint,visc,fac,dynvar->dta,iel);
      }/*end of if(dynvar->pro_calrhs==1)*/
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG     
dstrc_exit(); 
#endif

return; 
} /* end of f2pro_calint */

#endif
/*! @} (documentation module close)*/
