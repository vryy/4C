/*!----------------------------------------------------------------------
\file
\brief calculate submesh element size for fluid3

<pre>
Maintainer:  name
             e-mail
             homepage
             telephone number


</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3_prototypes.h"
#include "fluid3ml_prototypes.h"
#include "../fluid3/fluid3.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!---------------------------------------------------------------------
\brief routine to calculate submesh element size for fluid3

<pre>                                                       gravem 07/03

In this routine, the characteristic submesh element length is calculated
and the routine for the calculation of the stabilization parameter or 
the subgrid viscosity (depending on the respective flag) is called. 

</pre>
\param  *ele       ELEMENT	     (i)   actual element
\param  *data      FLUID_DATA	     (i)
\param  *dynvar    FLUID_DYN_CALC    (i/o)
\param  *mlvar     FLUID_DYN_ML      (i)
\param  *funct     DOUBLE 	     (-)   l-s shape functions
\param **deriv     DOUBLE 	     (-)   deriv. of l-s shape funcs
\param **deriv2    DOUBLE 	     (-)   2nd deriv. of l-s sh. funcs
\param  *smfunct   DOUBLE 	     (-)   submesh shape functions
\param **smderiv   DOUBLE 	     (-)   deriv. of submesh shape fun.
\param **smderiv2  DOUBLE 	     (-)   2nd deriv. of sm sh. fun.
\param **derxy     DOUBLE 	     (-)   global l-s sh. fun. deriv.
\param **xjm       DOUBLE 	     (-)   jacobian matrix
\param **evel      DOUBLE 	     (i)   l-s element velocities
\param  *velint    DOUBLE 	     (-)   l-s vel. at integr. point
\param **vderxy    DOUBLE 	     (-)   l-s vel. deriv. at int. point
\param **smxyze    DOUBLE 	     (i)   submesh element coordinates
\param **smxyzep   DOUBLE 	     (i)   sm ele. coord. on parent dom.
\param **wa1       DOUBLE 	     (-)   working array
\return void             

------------------------------------------------------------------------*/
void f3_smelesize(ELEMENT         *ele,    
		  FLUID_DATA	  *data, 
		  FLUID_DYN_CALC  *dynvar,
		  FLUID_DYN_ML    *mlvar,
	          DOUBLE	  *funct,  
	          DOUBLE	 **deriv,  
	          DOUBLE	 **deriv2,		
	          DOUBLE	  *smfunct,  
	          DOUBLE	 **smderiv,  
	          DOUBLE	 **smderiv2,		
	          DOUBLE	 **derxy,  
		  DOUBLE	 **xjm,    
		  DOUBLE	 **evel,		 
		  DOUBLE	  *velint, 
	          DOUBLE	 **vderxy,  
	          DOUBLE	 **smxyze,  
	          DOUBLE	 **smxyzep,  
		  DOUBLE	 **wa1)
{

INT     i,ilen,inod;    /* simply some counters	        		*/
INT     ntyp,nsmtyp;    /* l-s and sm element type (TRI or QUAD)        */
INT     actmat;         /* number of actual material		        */
INT     iel,smiel;      /* l-s and sm number of nodes of actual element */
DOUBLE  visc;           /* fluid viscosity                              */
DOUBLE  vol;            /* element volume                               */
DOUBLE  det;            /* determinant of jacobian                      */
DOUBLE  val;            /* temporary calculation value                  */
DOUBLE  velno;          /* velocity norm                                */
DOUBLE  strle;          /* streamlength                                 */
DOUBLE  e1,e2,e3;       /* natural coordinates of inegration point      */
DOUBLE  fac,facr;       /* factors                                      */
DOUBLE  facs,fact;      /* factors                                      */
DOUBLE  velino[3];      /* normed velocity vector at integration point  */
DOUBLE  gcoor[3],coor[3];/* global coordinates                           */
DIS_TYP typ,smtyp;

#ifdef DEBUG 
dstrc_enter("f3_smelesize");
#endif		

/*---------------------------------------------------------- initialize */
ntyp   = ele->e.f3->ntyp;
nsmtyp = mlvar->submesh.ntyp;
typ    = ele->distyp;
smtyp  = mlvar->submesh.typ;
iel    = ele->numnp;
smiel  = mlvar->submesh.numen;
vol   = ZERO;
strle = ZERO;
 
actmat=ele->mat-1;
visc = mat[actmat].m.fluid->viscosity;

/*------ get values of integration parameters, shape functions and their
         derivatives for submesh element -------------------------------*/
switch(nsmtyp)
{
case 1:    /* --> hex - element */
   e1	= data->qxg[0][0];
   facr = data->qwgt[0][0];
   e2	= data->qxg[0][0];
   facs = data->qwgt[0][0];
   e3	= data->qxg[0][0];
   fact = data->qwgt[0][0]; 
   f3_hex(smfunct,smderiv,smderiv2,e1,e2,e3,smtyp,2);
break;
case 2:       /* --> tet - element */		   
   e1	= data->txgr[0][0];
   facr = data->twgt[0][0];
   e2	= data->txgs[0][0];
   facs = ONE;
   e3	= data->txgs[0][0]; 
   fact = ONE;
   f3_tet(smfunct,smderiv,smderiv2,e1,e2,e3,smtyp,2);	   
break;
default:
   dserror("nsmtyp unknown!\n");      
} /*end switch(nsmtyp) */

if (mlvar->smesize<4)     /* compute volume */
{
/* -------------------------------------------- compute jacobian matrix */      
  f3_jaco3(smxyze,smfunct,smderiv,xjm,&det,smiel,ele);
  fac=facr*facs*fact*det;
  vol += fac;
}/* endif (mlvar->smesize<4) */
  
/* compute diagonal based diameter */
if (mlvar->smesize==4) dserror("no diagonal-based diameter in 3D yet!\n");

if (mlvar->smesize==5)    /* compute streamlength based on l-s velocity */
{
  f3_gcoor2(smfunct,smxyzep,smiel,coor);
  switch(ntyp)
  {
  case 1:    /* --> quad - element */
    f3_hex(funct,deriv,deriv2,coor[0],coor[1],coor[2],typ,2);
  break;
  case 2:	/* --> tri - element */ 	     
    f3_tet(funct,deriv,deriv2,coor[0],coor[1],coor[2],typ,2);
  break;
  default:
    dserror("ntyp unknown!\n");      
  } /*end switch(ntyp) */
  f3_veli(velint,funct,evel,iel);
  f3_jaco(funct,deriv,xjm,&det,ele,iel);
  f3_gder(derxy,deriv,xjm,wa1,det,iel); 
  val = ZERO;
  velno=sqrt( velint[0]*velint[0] \
  	    + velint[1]*velint[1] \
    	    + velint[2]*velint[2]);
  if(velno>=EPS6)
  {
     velino[0] = velint[0]/velno;
     velino[1] = velint[1]/velno;
     velino[2] = velint[2]/velno;
  }
  else
  {
     velino[0] = ONE;
     velino[1] = ZERO;
     velino[2] = ZERO;
  }	    
  for (inod=0;inod<iel;inod++) /* loop element nodes */
  {
     val += FABS(velino[0]*derxy[0][inod] \
    		+velino[1]*derxy[1][inod] \
    		+velino[2]*derxy[2][inod]);
  } /* end of loop over element nodes */
  strle=TWO/val;      
} /* endif (mlvar->smesize==5) */

/*----------------------------------- set characteristic element length */
if      (mlvar->smesize==1) ele->e.f3->smcml = pow(vol,(ONE/THREE));
else if (mlvar->smesize==2) ele->e.f3->smcml = pow((SIX*vol/PI),(ONE/THREE));
else if (mlvar->smesize==3) ele->e.f3->smcml = pow((SIX*vol/PI),(ONE/THREE))/sqrt(THREE);
else if (mlvar->smesize==5) ele->e.f3->smcml = strle;  

if (mlvar->smesize<5)     /* compute l-s velocity */
{
  f3_gcoor2(smfunct,smxyzep,smiel,coor);
  switch(ntyp)
  {
  case 1:    /* --> hex - element */
    f3_hex(funct,deriv,deriv2,coor[0],coor[1],coor[2],typ,2);
  break;
  case 2:	/* --> tet - element */ 	     
    f3_tet(funct,deriv,deriv2,coor[0],coor[1],coor[2],typ,2);
  break;
  default:
    dserror("ntyp unknown!\n");      
  } /*end switch(ntyp) */
  f3_veli(velint,funct,evel,iel);
}  

/*----------------------------------- calculate stabilization parameter */               
if (mlvar->smstabi>0) f3_smstabpar(ele,dynvar,mlvar,velint,visc,smiel,ntyp); 

/*--------------------------------------------------- subgrid viscosity */               
if (mlvar->smsgvi==1 || mlvar->smsgvi==2)
{ 
  f3_jaco(funct,deriv,xjm,&det,ele,iel);
/*-------------------------------------------- compute global derivates */
  f3_gder(derxy,deriv,xjm,wa1,det,iel);
/*----------------------- get velocity derivatives at integration point */
  f3_vder(vderxy,derxy,evel,iel);
/*----------------------------------------- calculate subgrid viscosity */               
  f3_smsgvisc(ele,mlvar,velint,vderxy,visc,smiel,ntyp);
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_smelesize */

#endif
