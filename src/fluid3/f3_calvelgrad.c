/*!----------------------------------------------------------------------
\file
\brief calculation of components of the velocity gradient

<pre>
Maintainer:  name
             e-mail
             homepage
             telephone number


</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"

/*!--------------------------------------------------------------------- 
\brief calculation of vel. derivative extrapolated to nodes for fluid3

<pre>                                                       gravem 10/02

In this routine, the velocity derivative dux/dy (or various other 
components of the velocity gradient) at the integration points is 
calculated, extrapolated to the nodes and stored at the nodes. 
The extrapolation is performed using a virtual element. The integration
points of the actual element represent the nodes of the virtual element.
The values are averaged depending on the number of elements the actual 
node belongs to.
			     
</pre>

\param  *data       FLUID_DATA  (i)
\param  *ele        ELMENT      (i/o)  actual element
\param **evel       DOUBLE      (-)    element velocities
\param  *funct      DOUBLE      (-)    shape functions
\param  *vefunct    DOUBLE      (-)    virtual elem. shape functions
\param **deriv      DOUBLE      (-)    natural deriv. of shape funct.
\param **derxy      DOUBLE      (-)    global deriv. of sape funct.       
\param **vderxy     DOUBLE      (-)    global vel. deriv
\param **xjm        DOUBLE      (-)    jacobian matrix
\param **xyze       DOUBLE      (-)    element coordinates
\param **vderint    DOUBLE      (-)    velocity derivatives at GAUSS point 
\param **wa1        DOUBLE      (-)    working array 
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_calvelgrad(FLUID_DATA	  *data, 
       	           ELEMENT	  *ele,
		   DOUBLE	 **evel, 
		   DOUBLE	  *funct,
		   DOUBLE	  *vefunct,
		   DOUBLE	 **deriv,
		   DOUBLE	 **derxy,
		   DOUBLE	 **vderxy,
		   DOUBLE	 **xjm,
		   DOUBLE	 **xyze,
		   DOUBLE	 **vderint,
	           DOUBLE        **wa1)
{
INT     i,j,node,lr,ls,lt;   /* some counters                          */
INT     iel,nir,nis,nit;     /* number of nodes/integr. points         */
INT     intc,icode;          /* flags                                  */
INT     ntyp;                /* flag for element type                  */
INT     iv;                  /* counter for GAUSS points               */
DOUBLE  det,val;             /* element values                         */
DOUBLE  e1,e2,e3,r,s,t;         
DIS_TYP typ,vetyp;	     /* element displacement type              */
NODE   *actnode;             /* actual node                            */

#ifdef DEBUG 
dstrc_enter("f3_calvelgrad");
#endif

/*----------------------------------------------------- initialisation */
iel  = ele->numnp;
ntyp = ele->e.f3->ntyp; 
typ  = ele->distyp;

/*------------------------------------------------ get integraton data  */
switch (ntyp)
{
case 1:  /* --> quad - element */
  icode = 2;
  nir = ele->e.f3->nGP[0];
  nis = ele->e.f3->nGP[1];
  nit = ele->e.f3->nGP[2];
break;
case 2: /* --> tri - element */  
  icode = 2;
  nir  = ele->e.f3->nGP[0];
  nis  = 1;
  nit  = 1; 
  intc = ele->e.f3->nGP[1];  
break;
default:
  dserror("ntyp unknown!");
} /* end switch(ntyp) */

/*----------------------------- set element velocities and coordinates */
for (j=0;j<iel;j++)
{
  evel[0][j] = ele->node[j]->sol_increment.a.da[3][0];
  evel[1][j] = ele->node[j]->sol_increment.a.da[3][1];
  evel[2][j] = ele->node[j]->sol_increment.a.da[3][2];
  xyze[0][j] = ele->node[j]->x[0];
  xyze[1][j] = ele->node[j]->x[1];
  xyze[2][j] = ele->node[j]->x[2];
}/*end for (j=0;j<iel;j++) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
iv=0;
for (lr=0;lr<nir;lr++)
{    
for (ls=0;ls<nis;ls++)
{
for (lt=0;lt<nit;lt++)
{
/*--------------- get values of  shape functions and their derivatives */
  switch(ntyp)  
  {
  case 1:   /* --> quad - element */
    e1 = data->qxg[lr][nir-1];           
    e2 = data->qxg[ls][nis-1];
    e3 = data->qxg[lt][nit-1];
    f3_hex(funct,deriv,NULL,e1,e2,e3,typ,icode);
  break;
  case 2:   /* --> tri - element */              
    e1 = data->txgr[lr][intc];
    e2 = data->txgs[lr][intc];
    e3 = data->txgt[lr][intc]; 
    f3_tet(funct,deriv,NULL,e1,e2,e3,typ,icode); 
  break;
  default:
    dserror("ntyp unknown!");
  } /* end switch(ntyp) */
  
/*-------------------------------------------- compute Jacobian matrix */
  f3_jaco2(xyze,deriv,xjm,&det,iel);
/*------------------------------------------- compute global derivates */
  f3_gder(derxy,deriv,xjm,wa1,det,iel);
/*--------------------- get velocity  derivatives at integration point */
  f3_vder(vderxy,derxy,evel,iel);
/*---------------
             | Ux,x  ->Ux,y   Ux,z |
             |                     |
    vderxy = | Uy,x    Uy,y   Uy,z |  
             |                     |
             | Uz,x    Uz,y   Uz,z |    */
/*  calculate velocity derivative Ux,y  */  
    vderint[1][iv] = vderxy[0][1];
    iv++;			     
} /* end loop over nit */
} /* end loop over nis */
} /* end loop over nir */

/*------ extrapolate velocity derivatives to the nodes using a virtual  */
/*                      element made up of integrations points as nodes */
for (node=0;node<iel;node++)
{
  actnode = ele->node[node];
/*------------------- determine coordinates of nodes in virtual element */
  r = f3_rstve(node,0,iel,iv);
  s = f3_rstve(node,1,iel,iv);
  t = f3_rstve(node,2,iel,iv);
    
/*------ calculate value at node using virtual element with integration */
/*                                points representing the virtual nodes */
  switch (iv)
  {
  case 8: 
    vetyp=hex8; 
  break;
  case 27: 
    vetyp=hex27; 
  break;
  default:
    dserror("number of integration points unknown");
  }/* end of switch(iv) */    
  if      (ntyp==1) f3_hex(vefunct,NULL,NULL,r,s,t,vetyp,1);
  else if (ntyp==2) dserror("Tetraeder not implemented yet!\n");
  else dserror("ntyp unknown!\n");
    
  val=ZERO;
  for (j=0;j<iv;j++)
  {
    val+=vefunct[j]*vderint[1][j];
  }  
    	 	 
/*------ gather derivatives at node as mean values of elements at node */
  actnode->deruxy+=(val/actnode->numele);    
} /* end for (node=0;node<iel;node++) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_calvelgrad */
#endif
