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
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"

static DOUBLE Q12=(ONE/TWO);

/*!--------------------------------------------------------------------- 
\brief calculation of vel. derivative extrapolated to nodes for fluid2

<pre>                                                       gravem 10/02

In this routine, the velocity derivative dux/dy (or various other 
components of the velocity gradient) at the integration points is 
calculated, extrapolated to the nodes and stored at the nodes. The 
values are averaged depending on the number of elements the actual 
node belongs to.
			     
</pre>

\param  *data       FLUID_DATA  (i)
\param  *ele        ELMENT      (i/o)  actual element
\param **evel       DOUBLE      (-)    element velocities
\param  *funct      DOUBLE      (-)    shape functions
\param **deriv      DOUBLE      (-)    natural deriv. of shape funct.
\param **derxy      DOUBLE      (-)    global deriv. of sape funct.       
\param **vderxy     DOUBLE      (-)    global vel. deriv
\param **xjm        DOUBLE      (-)    jacobian matrix
\param **xyze       DOUBLE      (-)    element coordinates
\param **vderint    DOUBLE      (-)    velocity derivatives at GAUSS point 
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_calvelgrad(FLUID_DATA	  *data, 
       	           ELEMENT	  *ele,
		   DOUBLE	 **evel, 
		   DOUBLE	  *funct,
		   DOUBLE	 **deriv,
		   DOUBLE	 **derxy,
		   DOUBLE	 **vderxy,
		   DOUBLE	 **xjm,
		   DOUBLE	 **xyze,
		   DOUBLE	 **vderint)
{
INT     i,j,node,lr,ls;      /* some counters                          */
INT     iel,nir,nis;         /* number of nodes/integr. points         */
INT     intc,icode;          /* flags                                  */
INT     ntyp;                /* flag for element type                  */
INT     iv;                  /* counter for GAUSS points               */
DOUBLE  det,val;             /* element values                         */
DOUBLE  e1,e2,r,s;         
DIS_TYP typ;	             /* element displacement type              */
DOUBLE  fpar[MAXGAUSS];      /* working array                          */
NODE   *actnode;             /* actual node                            */

#ifdef DEBUG 
dstrc_enter("f2_calvelgrad");
#endif

/*----------------------------------------------------- initialisation */
iel  = ele->numnp;
ntyp = ele->e.f2->ntyp; 
typ  = ele->distyp;

/*---------------------------------------------- get integration data  */
switch (ntyp)
{
case 1:  /* --> quad - element */
  icode = 2;
  nir = ele->e.f2->nGP[0];
  nis = ele->e.f2->nGP[1];
break;
case 2: /* --> tri - element */  
  icode = 2;
  nir  = ele->e.f2->nGP[0];
  nis  = 1;
  intc = ele->e.f2->nGP[1];  
break;
default:
  dserror("ntyp unknown!");
} /* end switch(ntyp) */

/*----------------------------- set element velocities and coordinates */
for (j=0;j<iel;j++)
{
  evel[0][j] = ele->node[j]->sol_increment.a.da[3][0];
  evel[1][j] = ele->node[j]->sol_increment.a.da[3][1];
  xyze[0][j] = ele->node[j]->x[0];
  xyze[1][j] = ele->node[j]->x[1];
}/*end for (j=0;j<iel;j++) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
iv=0;
for (lr=0;lr<nir;lr++)
{    
  for (ls=0;ls<nis;ls++)
  {
/*--------------- get values of  shape functions and their derivatives */
    switch(ntyp)  
    {
    case 1:   /* --> quad - element */
      e1 = data->qxg[lr][nir-1];           
      e2 = data->qxg[ls][nis-1];
      f2_rec(funct,deriv,NULL,e1,e2,typ,icode);
    break;
    case 2:   /* --> tri - element */              
      e1 = data->txgr[lr][intc];
      e2 = data->txgs[lr][intc];
      f2_tri(funct,deriv,NULL,e1,e2,typ,icode);
    break;
    default:
      dserror("ntyp unknown!");
    } /* end switch(ntyp) */
  
/*-------------------------------------------- compute Jacobian matrix */
    f2_jaco2(xyze,funct,deriv,xjm,&det,iel,ele);
/*------------------------------------------- compute global derivates */
    f2_gder(derxy,deriv,xjm,det,iel);
/*--------------------- get velocity  derivatives at integration point */
    f2_vder(vderxy,derxy,evel,iel);
/*---------------
             | Ux,x  ->Ux,y |
    vderxy = |              |
             | Uy,x    Uy,y |  */
/*  Alternative 1: calculate velocity derivative Ux,y   
    vderint[1][iv] = vderxy[0][1];  */ 	      
/*  Alternative 2: calculate 2xvorticity e.g. for plane mixing layer  */
    vderint[1][iv] = vderxy[1][0]-vderxy[0][1]; 
    iv++;			     
  } /* end loop over nis */
} /* end loop over nir */

/*----------------------- extrapolate velocity derivative to the nodes */
if      (ntyp==1) f2_recex(NULL,fpar,ZERO,ZERO,&(vderint[1][0]),iv,0);      
else if (ntyp==2) f2_triex(NULL,fpar,ZERO,ZERO,&(vderint[1][0]),iv,0); 
else dserror("ntyp unknown!\n");
        
for (node=0;node<iel;node++)
{
  actnode = ele->node[node];
  r = f2_rsn(node,0,iel);
  s = f2_rsn(node,1,iel);
  if      (ntyp==1) f2_recex(&val,fpar,r,s,NULL,iv,2);
  else if (ntyp==2) f2_triex(&val,fpar,r,s,NULL,iv,2); 
  else dserror("ntyp unknown!\n");	 	 
/*------ gather derivatives at node as mean values of elements at node */
  actnode->deruxy+=(val/actnode->numele);    
} /* end for (node=0;node<iel;node++) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calvelgrad */

/*!---------------------------------------------------------------------
\brief routine to calculate the vorticity

<pre>                                                         genk 12/02
   calculation of vorticity (for visualisation):
      vort = 1/2 (Ux,y + Uy,x)
   when extrapolating the vorticity from the gauss points to the nodes
   the results are averaged with the number of elements the actual node
   belongs to.
   An alternative would be to average with the element areas;
   smaller elements = better values!     
</pre>

\param  *data	            FLUID_DATA       (i)
\param  *dynvar	            FLUID_DYN_CALC   (i)
\param  *ele	            ELEMENT	     (i)   actual element
\param   init	            INT              (i)   init flag
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_calvort(FLUID_DATA     *data, 
                FLUID_DYN_CALC *dynvar, 
	        ELEMENT        *ele,                
       	        INT             init)
{
static ARRAY     evel_a;     /* element velocities at (n)                 */
static DOUBLE  **evel;
static ARRAY     funct_a;    /* shape functions                           */
static DOUBLE   *funct;
static ARRAY     deriv_a;    /* first natural derivatives                 */
static DOUBLE  **deriv;
static ARRAY     xjm_a;      /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     vderxy_a;   /* vel - derivatives                         */
static DOUBLE  **vderxy;
static ARRAY     derxy_a;    /* coordinate - derivatives                  */
static DOUBLE  **derxy;
static ARRAY     vort_a;     /* vorticity                                 */
static DOUBLE   *vort;      
static ARRAY     xyze_a;     /* element coordinates                       */
static DOUBLE  **xyze;   

INT j,icol,kk;
INT iv=0;		     /* index for gauss points                    */ 	
INT ivmax;		     /* index for maximum number of gauss points  */
INT icode=2;                 /* flag for eveluation of shape functions    */ 
INT iel=0;		     /* number of nodes                           */
INT lr,ls;		     /* index to loop over integration points     */
INT ntyp;		     /* element type: 1 - quad; 2 - tri           */
INT node;		     /* index to loop over nodal points           */
INT nir,nis;                 /* number of gauss points in r&s directions  */
INT intc;
INT ncols;                   /* number of columns                         */
INT numele;                  /* number of elements to the actual element  */
DOUBLE det;                  /* determinant                               */
DOUBLE e1,e2;                /* coordinates of the current gauss points   */
DOUBLE r;                    /* local coord. of gauss points in r-dir.    */
DOUBLE s;                    /* local coord. of gauss points in s-dir     */ 
DOUBLE f;                    /* vorticity value at the nodes              */
DOUBLE fpar[MAXGAUSS];
NODE *actnode;               /* actual node                               */
DIS_TYP typ;                 /* element type                              */

#ifdef DEBUG 
dstrc_enter("f2_calvort");
#endif

if (init==1) /* allocate working arrays and set pointers */
{
   evel    = amdef("evel"   ,&evel_a  ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
   funct   = amdef("funct"  ,&funct_a  ,MAXNOD_F2,1,"DV");
   deriv   = amdef("deriv"  ,&deriv_a  ,2,MAXNOD_F2,"DA");
   xjm     = amdef("xjm"    ,&xjm_a    ,2,2        ,"DA");
   vderxy  = amdef("vderxy" ,&vderxy_a ,2,2        ,"DA");
   derxy   = amdef("derxy"  ,&derxy_a  ,2,MAXNOD_F2,"DA");
   vort    = amdef("vort"   ,&vort_a   ,MAXGAUSS,1 ,"DV");   
   xyze    = amdef("xyze"   ,&xyze_a   ,2,MAXNOD_F2,"DA");
   goto end;
} /* endif (init==1) */

/*------------------------------------------------- calculate vorticity */
/*	                vort = 1/2*(Uy,x - Ux,y)                        */ 
/*----------------------------------------------------------------------*/

/*-------------------------------------------- loop over each time step */
ncols = dynvar->ncols;
for (icol=0;icol<ncols;icol++)  
{

/*------------------------------------------------------ initialization */
iel = ele->numnp;
ntyp= ele->e.f2->ntyp;
typ = ele->distyp;

/*------------------------------------------------- get integraton data */
switch (ntyp)
{
case 1: /* -----> quad element */
   icode = 2;
	/* initialise integration */
   nir = ele->e.f2->nGP[0];
   nis = ele->e.f2->nGP[1];
break;
case 2: /* triangular element */	
      icode = 2;
      nir  = ele->e.f2->nGP[0];
      nis  = 1;
      intc = ele->e.f2->nGP[1];  
break;	
default:
   dserror("ntyp unknown!");
} /* end switch(ntyp) */


/*---------------------- set element velocities, and coordinates -------*/
for (j=0;j<iel;j++)
{
   evel[0][j] = ele->node[j]->sol.a.da[icol][0];
   evel[1][j] = ele->node[j]->sol.a.da[icol][1];
   xyze[0][j] = ele->node[j]->x[0];
   xyze[1][j] = ele->node[j]->x[1];
}/*end for (j=0;j<iel;j++) */

/*---------------------------------- start loop over integration points */
iv = 0;
for (lr=0;lr<nir;lr++)
{
   for (ls=0;ls<nis;ls++)
      {
/*---------------- get values of  shape functions and their derivatives */
      switch(ntyp)
      {
      case 1:  /* quad element*/
         e1 = data->qxg[lr][nir-1];
         e2 = data->qxg[ls][nis-1];
         f2_rec(funct,deriv,NULL,e1,e2,typ,icode);
      break;
      case 2:   /* --> tri - element */              
	 e1   = data->txgr[lr][intc];
	 e2   = data->txgs[lr][intc];
	 f2_tri(funct,deriv,NULL,e1,e2,typ,icode);
      break;
      default:
         dserror("ntyp unknown!");
      }/*end of switch(ntyp) */		
      
/*--------------------------------------------- compute Jacobian matrix */
      f2_jaco2(xyze,funct,deriv,xjm,&det,iel,ele);
/*--------------------------------- compute global derivative ----------*/
      f2_gder(derxy,deriv,xjm,det,iel);
/*-------- compute velocity derivatives at integration points ----------*/
      f2_vder(vderxy,derxy,evel,iel);
/*----------------------- calculate vorticity --------------------------*/
      vort[iv] = Q12*(vderxy[1][0]-vderxy[0][1]); 
      iv++;
   }/* end of for (ls=0;ls<nis;ls++) */	
}/* end of for (lr=0;lr<nir;lr++) */

ivmax = iv;

/*------------------------------start loop over nodal points------------*/
for (node=0;node<iel;node++)
{
   actnode = ele->node[node];                          /*--actual node--*/ 
   numele  = actnode->numele; /*--number of elements to the actual node-*/
   if (actnode->sol.sdim < 4)
   {
       amredef(&(actnode->sol),actnode->sol.fdim,4,"DA"); 
       for (kk=0;kk<ncols;kk++)
       actnode->sol.a.da[kk][3] = ZERO;
   }
/*------ --------get local coordinates of nodes-------------------------*/
   r = f2_rsn(node,0,iel);
   s = f2_rsn(node,1,iel);
/*---------------------------------- extrapolate vorticity to the nodes */   
   if      (ntyp==1) f2_recex(&f,fpar,r,s,vort,ivmax,1); 
   else if (ntyp==2) f2_triex(&f,fpar,r,s,vort,ivmax,1); 
   else dserror("ntyp unknown!\n");
          
/*------------------- store the vorticity value in the solution history-*/   
   actnode->sol.a.da[icol][3] += f/numele;                    

}/* end of loop for (node=0;node<iel;node++) */
}/* end of loop over columns of solution history */
  

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calvort */

#endif
