/*!----------------------------------------------------------------------
\file
\brief integration loop for multilevel fluid3 elements

<pre>
Maintainer: Volker Gravemeier
            vgravem@stanford.edu
            
            
</pre>

------------------------------------------------------------------------*/
#ifdef FLUID3_ML 
#include "../headers/standardtypes.h"
#include "../fluid_full/fluid_prototypes.h"
#include "../fluid3/fluid3_prototypes.h"
#include "fluid3ml_prototypes.h"
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

static FLUID_DYNAMIC *fdyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!---------------------------------------------------------------------
\brief integration loop for submesh element for fluid3

<pre>                                                       gravem 07/03

In this routine, the element stiffness matrix, mass matrix, VMM-RHS and
Time-RHS for one submesh element is calculated.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *mlvar       FLUID_DYN_ML  (i)
\param  *submesh     FLUID_ML_SMESH(i)   
\param **smestif      DOUBLE       (o)	sm element stiffness matrix
\param **smemass      DOUBLE       (o)	sm element mass matrix
\param  *smevfor      DOUBLE       (o)	sm element VMM force vector
\param  *smetfor      DOUBLE       (o)	sm element time force vector
\param  *smxyze       DOUBLE       (i)	submesh element coordinates
\param  *smxyzep      DOUBLE       (i)	sm ele. coord. on parent dom.
\param  *funct        DOUBLE       (-)	natural shape functions
\param **deriv        DOUBLE       (-)	deriv. of nat. shape funcs
\param **deriv2       DOUBLE       (-)	2nd deriv. of nat. shape f.
\param **xjm	      DOUBLE       (-)	jacobian matrix
\param **derxy        DOUBLE       (-)	global derivatives
\param **derxy2       DOUBLE       (-)	2nd global derivatives
\param  *smfunct      DOUBLE       (-)	sm natural shape functions
\param **smderiv      DOUBLE       (-)  sm deriv. of nat. shape funcs
\param **smderiv2     DOUBLE       (-)  sm 2nd deriv. of nat. shape f.
\param **smxjm	      DOUBLE       (-)  sm jacobian matrix
\param **smderxy      DOUBLE       (-)  sm global derivatives
\param **smderxy2     DOUBLE       (-)  sm 2nd global derivatives
\param **eveln        DOUBLE       (i)  ele vel. at time step n
\param **evel         DOUBLE       (i)  ele vel. at time step n+1
\param  *epren        DOUBLE       (i)  ele pres. at time step n
\param  *epre         DOUBLE       (i)  ele pres. at time step n+1
\param **evbub        DOUBLE       (i)  sm ele vel. bubble functions
\param **epbub        DOUBLE       (i)  sm ele pre. bubble functions
\param **efbub        DOUBLE       (i)  sm ele rhs bubble functions
\param **evbubn       DOUBLE       (i)  sm ele vel. bubble fun. at n
\param **epbubn       DOUBLE       (i)  sm ele pre. bubble fun. at n
\param **efbubn       DOUBLE       (i)  sm ele rhs bubble fun. at n
\param  *vbubint      DOUBLE       (-)  vel. bubble fun. at int. p.
\param  *vbubderxy    DOUBLE       (-)  vel. bub. fun. der. at int. p.
\param  *vbubderxy2   DOUBLE       (-)  2nd vel. bub. fun. der. at i.p.
\param  *pbubint      DOUBLE       (-)  pre. bubble fun. at int. p.
\param  *pbubderxy    DOUBLE       (-)  pre. bub. fun. der. at int. p.
\param  *pbubderxy2   DOUBLE       (-)  2nd pre. bub. fun. der. at i.p.
\param  *vbubintn     DOUBLE       (-)  vel. bubble fun. at int. p. at n
\param  *vbubderxyn   DOUBLE       (-)  vel. bub. fun. der. at int. p. at n
\param  *vbubderxy2n  DOUBLE       (-)  2nd vel. bub. fun. der. at i.p. at n
\param  *pbubintn     DOUBLE       (-)  pre. bubble fun. at int. p. at n
\param  *pbubderxyn   DOUBLE       (-)  pre. bub. fun. der. at int. p. at n
\param  *pbubderxy2n  DOUBLE       (-)  2nd pre. bub. fun. der. at i.p. at n
\param  *velint       DOUBLE       (-)  vel at integration point
\param  *velintn      DOUBLE       (-)  vel at integration point at n
\param  *velintnt     DOUBLE       (-)  'temporal' vel at int. p. at n
\param  *velintnc     DOUBLE       (-)  'convective' vel at int. p. at n
\param **vderxy       DOUBLE       (-)  global vel. derivatives
\param **vderxyn      DOUBLE       (-)  global vel. derivatives at n
\param **vderxync     DOUBLE       (-)  global 'convective' vel. der. at n
\param **vderxynv     DOUBLE       (-)  global 'viscous' vel. der. at n
\param **vderxy2      DOUBLE       (-)  2nd global vel. deriv.
\param **vderxy2n     DOUBLE       (-)  2nd global vel. derivatives at n
\param **vderxy2nv    DOUBLE       (-)  2nd global 'viscous' vel. der. at n
\param  *pderxyn      DOUBLE       (-)  global pres. derivatives at n
\param  *smvelint     DOUBLE       (-)  sm vel at integration point
\param **smvderxy     DOUBLE       (-)  sm global vel. derivatives
\param  *smpreint     DOUBLE       (-)  sm pre at integration point
\param **smpderxy     DOUBLE       (-)  sm global pre. derivatives
\param  *smvelintn    DOUBLE       (-)  sm vel at integration point at n
\param **smvderxyn    DOUBLE       (-)  sm global vel. derivatives at n
\param **smvderxy2n   DOUBLE       (-)  2nd sm global vel. derivatives at n
\param  *smpreintn    DOUBLE       (-)  sm pre at integration point at n
\param **smpderxyn    DOUBLE       (-)  sm global pre. derivatives at n
\param **smpderxy2n   DOUBLE       (-)  2nd sm global pre. derivatives at n
\param  *smfint       DOUBLE       (-)  sm rhs at integration point
\param **smfderxy     DOUBLE       (-)  sm global rhs. derivatives
\param  *smfintn      DOUBLE       (-)  sm rhs at integration point at n
\param **smfderxyn    DOUBLE       (-)  sm global rhs. derivatives at n
\param **smfderxy2n   DOUBLE       (-)  2nd sm global rhs. derivatives at n
\param **wa1	      DOUBLE       (-)  working array
\param **wa2	      DOUBLE       (-)  working array
\return void                                                   

------------------------------------------------------------------------*/
void f3_smint(FLUID_DATA      *data,     
	      ELEMENT	      *ele,
	      FLUID_DYN_ML    *mlvar, 
	      FLUID_ML_SMESH  *submesh, 
              DOUBLE	     **smestif,   
	      DOUBLE	     **smemass,   
	      DOUBLE	     **smevfor, 
	      DOUBLE	     **smetfor, 
	      DOUBLE	     **smxyze, 
	      DOUBLE	     **smxyzep, 
	      DOUBLE	      *funct,	
	      DOUBLE	     **deriv,	
	      DOUBLE	     **deriv2,  
	      DOUBLE	     **xjm,	
	      DOUBLE	     **derxy,	
	      DOUBLE	     **derxy2,  
	      DOUBLE	      *smfunct,   
	      DOUBLE	     **smderiv,   
	      DOUBLE	     **smderiv2,  
	      DOUBLE	     **smxjm,	  
	      DOUBLE	     **smderxy,   
	      DOUBLE	     **smderxy2,  
	      DOUBLE	     **eveln,	
	      DOUBLE	     **evel,  
	      DOUBLE	      *epren,	
	      DOUBLE	      *epre,
	      DOUBLE	     **evbub,	
	      DOUBLE	     **epbub,	
	      DOUBLE	     **efbub,	
	      DOUBLE	     **evbubn,   
	      DOUBLE	     **epbubn,   
	      DOUBLE	     **efbubn,   
              DOUBLE	      *vbubint,    
              DOUBLE	     **vbubderxy,  
              DOUBLE	     **vbubderxy2, 
              DOUBLE	     **pbubint,    
              DOUBLE	    ***pbubderxy,  
              DOUBLE	    ***pbubderxy2, 
              DOUBLE	      *vbubintn,   
              DOUBLE	     **vbubderxyn, 
              DOUBLE	     **vbubderxy2n,
              DOUBLE	     **pbubintn,   
              DOUBLE	    ***pbubderxyn, 
              DOUBLE	    ***pbubderxy2n,
	      DOUBLE	      *velint,  
              DOUBLE	      *velintn,   
              DOUBLE	      *velintnt,  
              DOUBLE	      *velintnc,  
	      DOUBLE	     **vderxy,  
              DOUBLE	     **vderxyn,   
              DOUBLE	     **vderxync,  
              DOUBLE	     **vderxynv,  
	      DOUBLE	     **vderxy2, 
              DOUBLE	     **vderxy2n,  
              DOUBLE	     **vderxy2nv,  
              DOUBLE	      *pderxyn,   
              DOUBLE	      *smvelint,  
              DOUBLE	     **smvderxy,  
              DOUBLE	      *smpreint,  
              DOUBLE	     **smpderxy,  
              DOUBLE	      *smvelintn, 
              DOUBLE	     **smvderxyn, 
              DOUBLE	     **smvderxy2n, 
              DOUBLE	      *smpreintn,  
              DOUBLE	     **smpderxyn,  
              DOUBLE	     **smpderxy2n, 
              DOUBLE	      *smfint,	 
              DOUBLE	     **smfderxy,  
              DOUBLE	      *smfintn,    
              DOUBLE	     **smfderxyn,  
              DOUBLE	     **smfderxy2n,
	      DOUBLE	     **wa1,	
	      DOUBLE	     **wa2)
{ 
INT       i,j;        /* simply some counters                           */
INT       iel,smiel;  /* large-scale and submesh number of nodes        */
INT       nsmtyp;     /* l-s and submesh element type: 1 - hex; 2 - tet */
INT       intc;       /* "integration case" for tet for further infos
                          see f3_inpele.c and f3_intg.c                 */
INT       nir,nis,nit;/* number of integration nodes in r,s direction   */
INT       actmat;     /* material number of the element                 */
INT       ihoel=0;    /* flag for higher order large-scale elements     */
INT       ihoelsm=0;  /* flag for higher order submesh elements         */
INT       icode=2;    /* flag for eveluation of l-s shape functions     */     
INT       icodesm=2;  /* flag for eveluation of sm shape functions      */     
INT       lr,ls,lt;   /* counter for integration                        */
DOUBLE    dens;       /* density                                        */
DOUBLE    visc;       /* viscosity                                      */
DOUBLE    fac;
DOUBLE    facr,facs,fact;/* integration weights                         */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2,e3;   /* natural coordinates of integr. point           */
DOUBLE    coor[3];    /* coordinates                                    */
DIS_TYP   typ,smtyp;  /* large-scale and submesh element type           */

#ifdef DEBUG 
dstrc_enter("f3_smint");
#endif

/*----------------------------------------------------- initialization */
fdyn = alldyn[genprob.numff].fdyn;

iel    = ele->numnp;
smiel  = submesh->numen;
actmat = ele->mat-1;
dens   = mat[actmat].m.fluid->density;
visc   = mat[actmat].m.fluid->viscosity;
typ    = ele->distyp;
nsmtyp = submesh->ntyp; 
smtyp  = submesh->typ;

/*--- get integration data and check if sm-elements are "higher order" */
switch (nsmtyp)
{
case 1:  /* --> hex - element */
   icodesm = 3;
   ihoelsm = 1;
   /* initialize integration */
   nir = submesh->ngpr;
   nis = submesh->ngps;
   nit = submesh->ngpt;
break;
case 2: /* --> tet - element */  
   if (smiel>4)
   {
      icodesm = 3;
      ihoelsm = 1;
   }
   /* initialize integration */
   nir  = submesh->ngpr;
   nis  = 1;
   nit  = 1;
   intc = submesh->ngps;  
break;
default:
   dserror("nsmtyp unknown!");
} /* end switch(nsmtyp) */

/*---------------------------- check if ls-elements are "higher order" */
switch (typ)
{
case hex8: case hex20: case hex27:  /* --> hex - element */
   icode   = 3;
   ihoel   = 1;
break;
case tet4: case tet10: /* --> tet - element */  
   if (iel>4)
   {
     icode   = 3;
     ihoel   = 1;
   }
break;
default:
   dserror("typ unknown!");
} /* end switch(typ) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{    
   for (ls=0;ls<nis;ls++)
   {
   for (lt=0;lt<nit;lt++)
   {
/*-------- get values of submesh shape functions and their derivatives */
      switch(nsmtyp)  
      {
      case 1:   /* --> hex - element */
	 e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
         e3   = data->qxg[lt][nit-1];
         fact = data->qwgt[lt][nit-1];
         f3_hex(smfunct,smderiv,smderiv2,e1,e2,e3,smtyp,icodesm);
      break;
      case 2:   /* --> tet - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
         e3   = data->txgt[lr][intc]; 
         fact = ONE;
         f3_tet(smfunct,smderiv,smderiv2,e1,e2,e3,smtyp,icodesm); 
      break;
      default:
         dserror("typ unknown!");
      } /* end switch(nsmtyp) */
/*-------------------------------- compute Jacobian matrix for submesh */
      f3_mljaco3(smxyze,smfunct,smderiv,smxjm,&det,smiel,ele);
      fac = facr*facs*fact*det;
/*------------------------------- compute global derivates for submesh */
      f3_gder(smderxy,smderiv,smxjm,wa1,det,smiel);

/*- get coordinates of current int. p. on parent domain of l-s element */
      f3_mlgcoor2(smfunct,smxyzep,smiel,coor);

/*---- get values of large-scale shape functions and their derivatives */
      switch(typ)
      {
      case hex8: case hex20: case hex27:    /* --> hex - element */
        f3_hex(funct,deriv,deriv2,coor[0],coor[1],coor[2],typ,icode);
      break;
      case tet4: case tet10:	/* --> tet - element */ 	     
        f3_tet(funct,deriv,deriv2,coor[0],coor[1],coor[2],typ,icode); 
      break;
      default:
        dserror("typ unknown!\n");      
      } /*end switch(typ) */
/*------------ compute Jacobian matrix for large-scale shape functions */
      f3_mljaco(funct,deriv,xjm,&det,ele,iel);
/*----------- compute global derivates for large-scale shape functions */
      f3_gder(derxy,deriv,xjm,wa1,det,iel);
/*------- compute 2nd global derivates for large-scale shape functions */
      if (ihoel!=0) f3_mlgder2(ele,xjm,wa2,derxy,derxy2,deriv2,iel);

/*------------ get large-scale velocities (n+1,i) at integration point */
      f3_veci(velint,funct,evel,iel);
/*-- get large-scale velocity (n+1,i) derivatives at integration point */
      f3_vder(vderxy,derxy,evel,iel);
      
      if (mlvar->convel==0)
      { 
/*-----------------------------------------------------------------------
     get values of bubble functions and their derivatives 
----------------------------------------------------------------------- */
/*------------------- get velocity bubble functions at integraton point */
        fluid_bubint (vbubint,smfunct,evbub,smiel,mlvar->nvbub);            
/*-------- get velocity bubble function derivatives at integraton point */
        fluid_bubder (vbubderxy,smderxy,evbub,smiel,mlvar->nvbub,3);            
/*----------------- get 'pressure' bubble functions at integraton point */
        fluid_pbubint (pbubint,smfunct,epbub,smiel,iel,3);            
/*------ get 'pressure' bubble function derivatives at integraton point */
        fluid_pbubder (pbubderxy,smderxy,epbub,smiel,iel,3,3);            
/*--------------------------- get small-scale 'rhs' at integraton point */
        fluid_bubint (smfint,smfunct,efbub,smiel,3);            
/*--------------- get small-scale 'rhs' derivatives at integraton point */
        fluid_bubder (smfderxy,smderxy,efbub,smiel,3,3);            
	 
/*---------------------- get small-scale velocities at integraton point */
        f3_veci (smvelint,vbubint,evel,iel);	      
/*------------ get small-scale velocity derivatives at integraton point */
        f3_vder (smvderxy,vbubderxy,evel,iel);	      
/*--------------------- get small-scale 'pressures' at integraton point */
        f3_smprei (smpreint,pbubint,epre,iel);	      
/*---------- get small-scale 'pressure' derivatives at integraton point */
        f3_smpder (smpderxy,pbubderxy,epre,iel);	      

/*--- calculate velocities and velocity derivatives at integraton point */
        for (i=0; i<3; i++)
	{
	  velint[i] += smvelint[i] + smpreint[i] + smfint[i];
	  for (j=0; j<3; j++)
	  {
	    vderxy[i][j] += smvderxy[i][j] + smpderxy[i][j] + smfderxy[j][i];
	  }
	}
      }	    

/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrices for submesh             |
 *----------------------------------------------------------------------*/
/*-------------------------------------------------  compute matrix SMK */      
      f3_calsmk(mlvar,smestif,velint,vderxy,smfunct,smderxy,fac,
                visc,smiel);
/*-------------------------------------------------- compute matrix SMM */
      if (fdyn->nis==0 && mlvar->quastabub==0) 
        f3_calsmm(smemass,smfunct,fac,smiel);
      
/*----------------------------------------------------------------------*
 |         compute Stabilization matrices for submesh                   |
 *----------------------------------------------------------------------*/
      if (mlvar->smstabi>0)
      { 
/*----------------------- compute second global derivatives for submesh */ 
        if (ihoelsm!=0) 
	  f3_mlcogder2(smxyze,smxjm,wa2,smderxy,smderxy2,smderiv2,smiel);
/*---------------------------------------- stabilization for matrix SMK */
        f3_calstabsmk(mlvar,smestif,velint,vderxy,smfunct,smderxy,
                      smderxy2,fac,visc,smiel,ihoelsm);
/*---------------------------------------- stabilization for matrix SMM */
        if (fdyn->nis==0) 
	  f3_calstabsmm(mlvar,smemass,velint,vderxy,smfunct,smderxy,
	                smderxy2,fac,visc,smiel,ihoelsm); 
      } /* endif (mlvar->smstabi>0) */

/*----------------------------------------------------------------------*
 |         compute "VMM" Force Vectors                                  |
 *----------------------------------------------------------------------*/ 
/*----------------------- standard Galerkin part for "VMM" force vector */
      f3_calsmfv(mlvar,smevfor,velint,vderxy,smfunct,funct,derxy,
                 derxy2,fac,visc,smiel,iel,ihoel); 
/*--------------------------- stabilization part for "VMM" force vector */
      if (mlvar->smstabi>0)
        f3_calstabsmfv(mlvar,smevfor,velint,vderxy,smfunct,smderxy,
	               smderxy2,funct,derxy,derxy2,fac,visc,smiel,iel,
		       ihoelsm,ihoel);

 /*----------------------------------------------------------------------*
 |         compute "Time" Force Vector                                  |
 *----------------------------------------------------------------------*/
/*------- cancel "Time" force vector calculation for certain parameters */
      if (mlvar->transterm>1 && fdyn->thsr==ZERO && mlvar->quastabub!=0) 
        continue;
	
     if (fdyn->nis==0)
      {
/*------- get large-scale pressure derivatives (n) at integration point */
        if (fdyn->iprerhs>0) f3_pder(pderxyn,derxy,epren,iel);
/*------------------ get large-scale velocities (n) at integraton point */
        f3_veci(velintn,funct,eveln,iel);
/*------- get large-scale velocity derivatives (n) at integration point */
        f3_vder(vderxyn,derxy,eveln,iel);
/*--- get large-scale 2nd velocity derivatives (n) at integration point */
        if (ihoel!=0) f3_vder2(vderxy2n,derxy2,eveln,iel);
      
        if (mlvar->quastabub==0)
        { 
/*--------------- get velocity bubble functions (n) at integraton point */
          fluid_bubint (vbubintn,smfunct,evbubn,smiel,mlvar->nvbub);            
/*---- get velocity bubble function derivatives (n) at integraton point */
          fluid_bubder (vbubderxyn,smderxy,evbubn,smiel,mlvar->nvbub,3);            
/* get 2nd velocity bubble function derivatives (n) at integraton point */
          if (ihoelsm!=0) 
	  fluid_bubder (vbubderxy2n,smderxy2,evbubn,smiel,mlvar->nvbub,6);
	              
/*------------- get 'pressure' bubble functions (n) at integraton point */
          fluid_pbubint (pbubintn,smfunct,epbubn,smiel,iel,3);            
/*-- get 'pressure' bubble function derivatives (n) at integraton point */
          fluid_pbubder (pbubderxyn,smderxy,epbubn,smiel,iel,3,3);            
/*--- get 2nd 'pressure' bubble function deriv. (n) at integraton point */
          if (ihoelsm!=0) 
	  fluid_pbubder (pbubderxy2n,smderxy2,epbubn,smiel,iel,3,6);
	  
/*----------------------- get small-scale 'rhs' (n) at integraton point */
          fluid_bubint (smfintn,smfunct,efbubn,smiel,3);            
/*----------- get small-scale 'rhs' derivatives (n) at integraton point */
          fluid_bubder (smfderxyn,smderxy,efbubn,smiel,3,3);            
/*------- get small-scale 2nd 'rhs' derivatives (n) at integraton point */
          if (ihoelsm!=0) 
	  fluid_bubder (smfderxy2n,smderxy2,efbubn,smiel,3,6);
	   
/*------------------ get small-scale velocities (n) at integraton point */
          f3_veci (smvelintn,vbubintn,eveln,iel);	      
/*-------- get small-scale velocity derivatives (n) at integraton point */
          f3_vder (smvderxyn,vbubderxyn,eveln,iel);	      
/*---- get small-scale 2nd velocity derivatives (n) at integraton point */
          if (ihoelsm!=0) f3_vder2(smvderxy2n,vbubderxy2n,eveln,iel);
	  	      
/*----------------- get small-scale 'pressures' (n) at integraton point */
          f3_smprei (smpreintn,pbubintn,epren,iel);	      
/*------ get small-scale 'pressure' derivatives (n) at integraton point */
          f3_smpder (smpderxyn,pbubderxyn,epren,iel);	      
/*-- get small-scale 2nd 'pressure' derivatives (n) at integraton point */
          if (ihoelsm!=0) f3_smpder2(smpderxy2n,pbubderxy2n,epren,iel);
	}
	  
/*----- calculate various velocities for "Time" force vector evaluation */
        for (i=0; i<3; i++)
	{
/*------------------------------------------------ quasi-static bubbles */
          if (mlvar->quastabub!=0)
	  {
/*----------------- velocities for temporal part of "Time" force vector */
	    if (mlvar->transterm>1) velintnt[i] = ZERO;
	    else                    velintnt[i] = velintn[i];
/*--------------- velocities for convective part of "Time" force vector */
	    if (mlvar->convel==0) velintnc[i] = velintn[i] + smvelint[i]\
	                                      + smpreint[i] + smfint[i];
	    else                  velintnc[i] = velintn[i];
	    for (j=0; j<3; j++)
	    {
/*----- velocity derivatives for convective part of "Time" force vector */
	      if (mlvar->convel==0) vderxync[i][j] = vderxyn[i][j]\
	          + smvderxy[i][j] + smpderxy[i][j] + smfderxy[j][i];
	      else                  vderxync[i][j] = vderxyn[i][j];
	    }
/*---- 2nd velocity derivatives for viscous part of "Time" force vector */
	    if (ihoel!=0 && ihoelsm!=0)
	    {
	      for (j=0; j<6; j++)
	      {
	        vderxy2nv[i][j] = vderxy2n[i][j];
	      }
	    }
	  }
/*--------------------------------------------- no quasi-static bubbles */
	  else
	  {
/*----------------- velocities for temporal part of "Time" force vector */
	    if (mlvar->transterm>1) velintnt[i] = smvelintn[i]\
	                                        + smpreintn[i] + smfintn[i];
	    else                    velintnt[i] = velintn[i] + smvelintn[i]\
	                                        + smpreintn[i] + smfintn[i];
/*--------------- velocities for convective part of "Time" force vector */
	    if (mlvar->convel==0) velintnc[i] = velintn[i] + smvelintn[i]\
	                                      + smpreintn[i] + smfintn[i];
	    else                  velintnc[i] = velintn[i];
/*-------------------------- general velocities for "Time" force vector */
	    velintn[i] += smvelintn[i] + smpreintn[i] + smfintn[i];
	    for (j=0; j<3; j++)
	    {
/*----- velocity derivatives for convective part of "Time" force vector */
	      if (mlvar->convel==0) vderxync[i][j] = vderxyn[i][j]\
	        + smvderxyn[i][j] + smpderxyn[i][j] + smfderxyn[j][i];
	      else                  vderxync[i][j] = vderxyn[i][j];
/*-------- velocity derivatives for viscous part of "Time" force vector */
	      vderxynv[i][j] = smvderxyn[i][j] + smpderxyn[i][j] + smfderxyn[j][i];
/*-- general velocity deriv. for convective part of "Time" force vector */
	      vderxyn[i][j] += smvderxyn[i][j] + smpderxyn[i][j] + smfderxyn[j][i];
	    }
	    if (ihoel!=0 && ihoelsm!=0)
	    {
	      for (j=0; j<6; j++)
	      {
/* general 2nd velocity der. for convective part of "Time" force vector */
	        vderxy2n[i][j] += smvderxy2n[i][j] + smpderxy2n[i][j]\
	                        + smfderxy2n[j][i];
/*---- 2nd velocity derivatives for viscous part of "Time" force vector */
	        vderxy2nv[i][j] = vderxy2n[i][j];
	      }
	    }
	  }
	}
	
/*---------------------- standard Galerkin part for "Time" force vector */
        f3_calsmft(mlvar,smetfor,velintn,velintnt,velintnc,vderxyn, 
		   vderxync,vderxynv,vderxy2nv,pderxyn,smfunct,smderxy,fac,
		   visc,smiel,iel,ihoel);
/*-------------------------- stabilization part for "Time" force vector */
        if (mlvar->smstabi>0)
          f3_calstabsmft(mlvar,smetfor,velintn,velintnt,velintnc,
	                 vderxyn,vderxync,vderxynv,vderxy2n,pderxyn,smfunct,
			 smderxy,smderxy2,fac,visc,smiel,iel,ihoelsm,ihoel);
      }  
   } /* end of loop over integration points lt*/
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_smint */

/*!---------------------------------------------------------------------
\brief integration loop for bubble funct. on submesh element for fluid3

<pre>                                                       gravem 07/03

In this routine, the bubble function part of the large-scale element 
stiffness matrix, mass matrix, VMM-RHS and Time-RHS for one submesh 
element is calculated.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *mlvar       FLUID_DYN_ML  (i)
\param  *submesh     FLUID_ML_SMESH(i)   
\param **estif        DOUBLE	   (o)  element stiffness matrix
\param **emass        DOUBLE	   (o)  element mass matrix
\param  *eiforce      DOUBLE	   (o)  element iteration force vector
\param  *smxyze	      DOUBLE	   (i)  submesh element coordinates
\param  *smxyzep      DOUBLE	   (i)  sm ele. coord. on parenyt dom.
\param  *funct        DOUBLE       (-)	natural shape functions
\param **deriv        DOUBLE       (-)	deriv. of nat. shape funcs
\param **deriv2       DOUBLE       (-)	2nd deriv. of nat. shape f.
\param **xjm	      DOUBLE       (-)	jacobian matrix
\param **derxy        DOUBLE       (-)	global derivatives
\param  *smfunct      DOUBLE       (-)	sm natural shape functions
\param **smderiv      DOUBLE       (-)  sm deriv. of nat. shape funcs
\param **smderiv2     DOUBLE       (-)  sm 2nd deriv. of nat. shape f.
\param **smxjm	      DOUBLE       (-)  sm jacobian matrix
\param **smderxy      DOUBLE       (-)  sm global derivatives
\param **evel         DOUBLE       (i)  ele vel. at time step n+1
\param  *epre         DOUBLE       (i)  ele pres. at time step n+1
\param **evbub        DOUBLE       (i)  sm ele vel. bubble functions
\param **epbub        DOUBLE       (i)  sm ele pre. bubble functions
\param **efbub        DOUBLE       (i)  sm ele rhs bubble functions
\param  *vbubint      DOUBLE       (-)  vel. bubble fun. at int. p.
\param  *vbubderxy    DOUBLE       (-)  vel. bub. fun. der. at int. p.
\param  *pbubint      DOUBLE       (-)  pre. bubble fun. at int. p.
\param  *pbubderxy    DOUBLE       (-)  pre. bub. fun. der. at int. p.
\param  *covint       DOUBLE       (-)  convective vel at int. point
\param  *velint       DOUBLE       (-)  vel at integration point
\param **vderxy       DOUBLE       (-)  global vel. derivatives
\param  *smvelint     DOUBLE       (-)  sm vel at integration point
\param **smvderxy     DOUBLE       (-)  sm global vel. derivatives
\param  *smpreint     DOUBLE       (-)  sm pre at integration point
\param **smpderxy     DOUBLE       (-)  sm global pre. derivatives
\param  *smfint       DOUBLE       (-)  sm rhs at integration point
\param **smfderxy     DOUBLE       (-)  sm global rhs. derivatives
\param **wa1	      DOUBLE       (-)  working array
\param **wa2	      DOUBLE       (-)  working array
\return void                                                   

------------------------------------------------------------------------*/
void f3_bubint(FLUID_DATA      *data,     
	       ELEMENT         *ele,
	       FLUID_DYN_ML    *mlvar, 
	       FLUID_ML_SMESH  *submesh, 
               DOUBLE	      **estif,	
	       DOUBLE	      **emass,	
	       DOUBLE	       *eiforce, 
	       DOUBLE         **smxyze, 
	       DOUBLE         **smxyzep, 
	       DOUBLE          *funct,   
	       DOUBLE         **deriv,   
	       DOUBLE         **deriv2,   
	       DOUBLE         **xjm,     
	       DOUBLE         **derxy,   
	       DOUBLE          *smfunct,   
	       DOUBLE         **smderiv,   
	       DOUBLE         **smderiv2,   
	       DOUBLE         **smxjm,     
	       DOUBLE         **smderxy,   
	       DOUBLE         **evel,  
	       DOUBLE          *epre,
	       DOUBLE         **evbub,   
	       DOUBLE         **epbub,   
	       DOUBLE         **efbub,   
               DOUBLE          *vbubint,    
               DOUBLE         **vbubderxy,  
               DOUBLE         **pbubint,    
               DOUBLE        ***pbubderxy,  
	       DOUBLE	       *covint,  
	       DOUBLE          *velint,  
	       DOUBLE         **vderxy,  
               DOUBLE          *smvelint,  
               DOUBLE         **smvderxy,  
               DOUBLE          *smpreint,  
               DOUBLE         **smpderxy,  
               DOUBLE          *smfint,    
               DOUBLE         **smfderxy,  
	       DOUBLE         **wa1,     
	       DOUBLE         **wa2)
{ 
INT       i,j;        /* simply some counters                           */
INT       iel,smiel;  /* number of nodes                                */
INT       nsmtyp;     /* l-s and submesh element type: 1 - hex; 2 - tet */
INT       intc;       /* "integration case" for tet for further infos
                          see f3_inpele.c and f3_intg.c                 */
INT       nir,nis,nit;/* number of integration nodesin r,s direction    */
INT       actmat;     /* material number of the element                 */
INT       ihoel=0;    /* flag for higher order large-scale elements     */
INT       ihoelsm=0;  /* flag for higher order submesh elements         */
INT       icode=2;    /* flag for eveluation of l-s shape functions     */     
INT       icodesm=2;  /* flag for eveluation of sm shape functions      */     
INT       lr,ls,lt;   /* counter for integration                        */
DOUBLE    dens;       /* density                                        */
DOUBLE    visc;       /* viscosity                                      */
DOUBLE    fac;
DOUBLE    facr,facs,fact;/* integration weights                         */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2,e3;   /* natural coordinates of integr. point           */
DOUBLE    coor[3];    /* coordinates                                    */
DIS_TYP   typ,smtyp;  /* large-scale and submesh element type           */

#ifdef DEBUG 
dstrc_enter("f3_bubint");
#endif

/*----------------------------------------------------- initialization */
fdyn = alldyn[genprob.numff].fdyn;

iel    = ele->numnp;
smiel  = submesh->numen;
actmat = ele->mat-1;
dens   = mat[actmat].m.fluid->density;
visc   = mat[actmat].m.fluid->viscosity;
typ    = ele->distyp;
nsmtyp = submesh->ntyp; 
smtyp  = submesh->typ;

/*---- get integraton data and check if sm-elements are "higher order" */
switch (nsmtyp)
{
case 1:  /* --> hex - element */
   icodesm = 3;
   ihoelsm = 1;
   /* initialixe integration */
   nir = submesh->ngpr;
   nis = submesh->ngps;
   nit = submesh->ngpt;
break;
case 2: /* --> tet - element */  
   if (smiel>4)
   {
      icodesm = 3;
      ihoelsm = 1;
   }
   /* initialize integration */
   nir  = submesh->ngpr;
   nis  = 1;
   nit  = 1;
   intc = submesh->ngps;  
break;
default:
   dserror("nsmtyp unknown!");
} /* end switch(nsmtyp) */

/*---------------------------- check if ls-elements are "higher order" */
switch (typ)
{
case hex8: case hex20: case hex27:  /* --> hex - element */
   icode   = 3;
   ihoel   = 1;
break;
case tet4: case tet10: /* --> tet - element */  
   if (iel>4)
   {
      icode   = 3;
      ihoel   = 1;
   }
break;
default:
   dserror("typ unknown!");
} /* end switch(typ) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{    
   for (ls=0;ls<nis;ls++)
   {
   for (lt=0;lt<nit;lt++)
   {
/*-------- get values of submesh shape functions and their derivatives */
      switch(nsmtyp)  
      {
      case 1:   /* --> hex - element */
	 e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
         e3   = data->qxg[lt][nit-1];
         fact = data->qwgt[lt][nit-1];
         f3_hex(smfunct,smderiv,smderiv2,e1,e2,e3,smtyp,icodesm);
      break;
      case 2:   /* --> tet - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
         e3   = data->txgt[lr][intc]; 
         fact = ONE;
	 f3_tet(smfunct,smderiv,smderiv2,e1,e2,e3,smtyp,icodesm);
      break;
      default:
         dserror("nsmtyp unknown!");
      } /* end switch(nsmtyp) */
/*-------------------------------- compute Jacobian matrix for submesh */
      f3_mljaco3(smxyze,smfunct,smderiv,smxjm,&det,smiel,ele);
      fac = facr*facs*fact*det;
/*------------------------------- compute global derivates for submesh */
      f3_gder(smderxy,smderiv,smxjm,wa1,det,smiel);

/*- get coordinates of current int. p. on parent domain of l-s element */
      f3_mlgcoor2(smfunct,smxyzep,smiel,coor);

/*---- get values of large-scale shape functions and their derivatives */
      switch(typ)
      {
      case hex8: case hex20: case hex27:    /* --> hex - element */
        f3_hex(funct,deriv,deriv2,coor[0],coor[1],coor[2],typ,icode);
      break;
      case tet4: case tet10:	/* --> tri - element */ 	     
        f3_tet(funct,deriv,deriv2,coor[0],coor[1],coor[2],typ,icode);   
      break;
      default:
        dserror("typ unknown!\n");      
      } /*end switch(typ) */
/*------------ compute Jacobian matrix for large-scale shape functions */
      f3_mljaco(funct,deriv,xjm,&det,ele,iel);
/*----------- compute global derivates for large-scale shape functions */
      f3_gder(derxy,deriv,xjm,wa1,det,iel);

/*------------ get large-scale velocities (n+1,i) at integration point */
      f3_veci(velint,funct,evel,iel);
/*-- get large-scale velocity (n+1,i) derivatives at integration point */
      f3_vder(vderxy,derxy,evel,iel);
      
/*-----------------------------------------------------------------------
     get values of bubble functions and their derivatives 
----------------------------------------------------------------------- */
/*------------------ get velocity bubble functions at integration point */
      fluid_bubint (vbubint,smfunct,evbub,smiel,mlvar->nvbub);            
/*------- get velocity bubble function derivatives at integration point */
      fluid_bubder (vbubderxy,smderxy,evbub,smiel,mlvar->nvbub,3);            
/*---------------- get 'pressure' bubble functions at integration point */
      fluid_pbubint (pbubint,smfunct,epbub,smiel,iel,3);            
/*----- get 'pressure' bubble function derivatives at integration point */
      fluid_pbubder (pbubderxy,smderxy,epbub,smiel,iel,3,3);            
/*-------------------------- get small-scale 'rhs' at integration point */
      fluid_bubint (smfint,smfunct,efbub,smiel,3);            
/*-------------- get small-scale 'rhs' derivatives at integration point */
      fluid_bubder (smfderxy,smderxy,efbub,smiel,3,3);            
	 
      if (mlvar->convel==0)
      { 
/*--------------------- get small-scale velocities at integration point */
        f3_veci (smvelint,vbubint,evel,iel);	      
/*----------- get small-scale velocity derivatives at integration point */
        f3_vder (smvderxy,vbubderxy,evel,iel);	      
/*-------------------- get small-scale 'pressures' at integration point */
        f3_smprei (smpreint,pbubint,epre,iel);	      
/*--------- get small-scale 'pressure' derivatives at integration point */
        f3_smpder (smpderxy,pbubderxy,epre,iel);	      

/*-- calculate velocities and velocity derivatives at integration point */
        for (i=0; i<3; i++)
	{
	  velint[i] += smvelint[i] + smpreint[i] + smfint[i];
	  for (j=0; j<3; j++)
	  {
	    vderxy[i][j] += smvderxy[i][j] + smpderxy[i][j] + smfderxy[j][i];
	  }
	}
      }	    

/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrix                           |
 | NOTE:                                                                |
 |  Standard Galerkin matrices are all stored in one matrix "estif"     |
 |  Standard Galerkin mass matrix is stored in "emass"                  |
 *----------------------------------------------------------------------*/
      if (mlvar->convel==0)
      {
/*------------------------ compute standard Galerkin part of matrix Kvv */      
        f3_lscalkvv(estif,velint,vderxy,funct,derxy,fac,visc,iel);
      }	 

/*----------------------------------------------------------------------*
 |         compute "Bubble" matrix                                      |
 | NOTE:                                                                |
 |  Bubble matrices are all stored in one matrix "estif"                |
 |  Bubble mass matrix is stored in "emass"                             |
 *----------------------------------------------------------------------*/
/*-------------------------- compute bubble function part of matrix Kvv */      
      f3_calbkvv(estif,velint,vderxy,funct,derxy,vbubint,vbubderxy,
                 fac,visc,iel);
/*-------------------------- compute bubble function part of matrix Kvp */      
      f3_calbkvp(estif,velint,vderxy,funct,derxy,pbubint,pbubderxy,
                 fac,visc,iel);
/*-------------------------- compute bubble function part of matrix Kpv */      
      f3_calbkpv(estif,funct,vbubderxy,fac,iel);
/*-------------------------- compute bubble function part of matrix Kpp */      
      f3_calbkpp(estif,funct,pbubderxy,fac,iel);
	
      if ((fdyn->nis==0 && mlvar->transterm==0) || 
          (fdyn->nis==0 && mlvar->transterm==2))
      {	   
/*-------------------------- compute bubble function part of matrix Mvv */      
	f3_calbmvv(emass,funct,vbubint,fac,iel);
/*-------------------------- compute bubble function part of matrix Mvp */      
	f3_calbmvp(emass,funct,pbubint,fac,iel);
      }    	
/*----------------------------------------------------------------------*
 |         compute "VMM" Force Vectors                                  |
 *----------------------------------------------------------------------*/ 
      if (fdyn->nis==0)
      {
/*------------------------ compute "VMM" force vector for velocity dofs */      
        f3_calbfv(mlvar,eiforce,velint,vderxy,funct,derxy,smfint,
                  smfderxy,fac,visc,iel); 
/*------------------------ compute "VMM" force vector for pressure dofs */      
        f3_calbfp(eiforce,funct,smfderxy,fac,iel); 
      }
      
/*----------------------------------------------------------------------*
 |        compute "Iteration" Force Vectors                             |
 *----------------------------------------------------------------------*/ 
      if (fdyn->nii!=0)
      {
        if (mlvar->convel!=0)
        {
/*--------------------- get small-scale velocities at integration point */
          f3_veci (smvelint,vbubint,evel,iel);	      
/*----------- get small-scale velocity derivatives at integration point */
          f3_vder (smvderxy,vbubderxy,evel,iel);	      
/*-------------------- get small-scale 'pressures' at integration point */
          f3_smprei (smpreint,pbubint,epre,iel);	      
/*--------- get small-scale 'pressure' derivatives at integration point */
          f3_smpder (smpderxy,pbubderxy,epre,iel);	      

/*-- calculate velocities and velocity derivatives at integration point */
          for (i=0; i<3; i++)
	  {
	    velint[i] = smvelint[i] + smpreint[i] + smfint[i];
	    for (j=0; j<3; j++)
	    {
	      vderxy[i][j] = smvderxy[i][j] + smpderxy[i][j] + smfderxy[j][i];
	    }
	  }
	}
/*-------------- get convective velocities (n+1,i) at integration point */
        f3_covi(vderxy,velint,covint);
/*- calculate "Iteration" force vector for velocity dofs (no pre. dofs) */
        f3_lscalgalifv(eiforce,covint,velint,vderxy,funct,fac,iel);
      } /* endif (fdyn->nii!=0) */
   } /* end of loop over integration points lt*/
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_bubint */

/*!---------------------------------------------------------------------
\brief integration loop for one large-scale element for fluid3

<pre>                                                       gravem 07/03

In this routine the large-scale part of the large-scale element 
stiffness matrix, iteration-RHS, time-RHS and external-RHS for one 
large-scale element is calculated
      
</pre>
\param  *data      FLUID_DATA	   (i)	  integration data
\param  *ele	   ELEMENT	   (i)    actual element
\param  *mlvar     FLUID_DYN_ML    (i)
\param  *hasext    int             (i)    element flag
\param **estif     DOUBLE	   (o)    element stiffness matrix
\param **emass     DOUBLE	   (o)    element mass matrix
\param  *eiforce   DOUBLE	   (o)    element iter force vector
\param  *etforce   DOUBLE	   (o)    element time force vector
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **deriv2    DOUBLE	   (-)    2nd deriv. of nat. shape f.
\param **xjm	   DOUBLE	   (-)    jacobian matrix
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **derxy2    DOUBLE	   (-)    2nd global derivatives
\param **evel      DOUBLE	   (i)    ele vel. at time n+1
\param **eveln     DOUBLE	   (i)    ele vel. at time n
\param  *epren     DOUBLE	   (-)    ele pres. at time n
\param  *edeadn    DOUBLE	   (-)    ele dead load at n 
\param  *edead     DOUBLE	   (-)    ele dead load at n+1
\param  *velint    DOUBLE	   (-)    vel at integration point
\param  *velintn   DOUBLE	   (-)    vel at integration point at n
\param  *covint    DOUBLE	   (-)    conv. vel. at integr. point
\param  *covintn   DOUBLE	   (-)    conv. vel. at integr. point at n
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\param **vderxyn   DOUBLE	   (-)    global vel. derivatives at n
\param  *pderxy    DOUBLE	   (-)    global pres. derivatives
\param **wa1	   DOUBLE	   (-)    working array
\param **wa2	   DOUBLE	   (-)    working array
\return void                                                   

------------------------------------------------------------------------*/
void f3_lsint(FLUID_DATA      *data,     
	      ELEMENT	      *ele,
	      FLUID_DYN_ML    *mlvar, 
              INT             *hasext,
              DOUBLE	     **estif,	
	      DOUBLE	     **emass,	
	      DOUBLE	      *eiforce, 
	      DOUBLE	      *etforce, 
	      DOUBLE	      *funct,	
	      DOUBLE	     **deriv,	
	      DOUBLE	     **deriv2,  
	      DOUBLE	     **xjm,	
	      DOUBLE	     **derxy,	
	      DOUBLE	     **derxy2,  
	      DOUBLE	     **evel,  
	      DOUBLE	     **eveln,  
	      DOUBLE          *epren,   
	      DOUBLE          *edeadn,   
	      DOUBLE          *edead,   
	      DOUBLE	      *velint,  
              DOUBLE          *velintn,   
	      DOUBLE	      *covint,  
	      DOUBLE	      *covintn,  
	      DOUBLE	     **vderxy,  
              DOUBLE         **vderxyn,  
	      DOUBLE	     **wa1,	
	      DOUBLE	     **wa2)
{ 
INT      iel;	      /* number of nodes 			        */
INT      intc;        /* "integration case" for tet for further infos
         		 see f3_inpele.c and f3_intg.c                  */
INT      nir,nis,nit; /* number of integration nodes in r,s,t direction */
INT      actmat;      /* material number of the element                 */
INT      ihoel=0;     /* flag for higher order elements                 */
INT      icode=2;     /* flag for eveluation of shape functions         */   
INT      lr, ls, lt;  /* counter for integration                        */
DOUBLE   dens;        /* density                                        */
DOUBLE   visc;        /* viscosity                                      */
DOUBLE   fac;
DOUBLE   facr,facs,fact; /* integration weights                         */
DOUBLE   det;	      /* determinant of jacobian matrix                 */
DOUBLE   e1,e2,e3;    /* natural coordinates of integr. point           */
DOUBLE   preintn;     /* pressure at integration point at times step n  */
DIS_TYP  typ;         /* element type                                   */
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f3_lsint");
#endif

/*----------------------------------------------------- initialization */
fdyn = alldyn[genprob.numff].fdyn;

iel=ele->numnp;
actmat=ele->mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity;
typ  = ele->distyp;
gls    = ele->e.f3->stabi.gls;

if (ele->e.f3->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");
   
/*------- get integraton data and check if elements are "higher order" */
switch (typ)
{
case hex8: case hex20: case hex27:  /* --> hex - element */
   icode   = 3;
   ihoel   = 1;
   /* initialize integration */
   nir = ele->e.f3->nGP[0];
   nis = ele->e.f3->nGP[1];
   nit = ele->e.f3->nGP[2];
   break;
case tet4: case tet10: /* --> tet - element */  
   if (iel>4)
   {
      icode   = 3;
      ihoel   = 1;
   }
   /* initialize integration */
   nir  = ele->e.f3->nGP[0];
   nis  = 1;
   nit  = 1; 
   intc = ele->e.f3->nGP[1];  
   break;
default:
   dserror("typ unknown!");
} /* end switch (typ) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{
for (ls=0;ls<nis;ls++)
{
for (lt=0;lt<nit;lt++)
{
/*---------------- get values of  shape functions and their derivatives */
   switch(typ)  
   {
   case hex8: case hex20: case hex27:   /* --> hex - element */
      e1   = data->qxg[lr][nir-1];
      facr = data->qwgt[lr][nir-1];
      e2   = data->qxg[ls][nis-1];
      facs = data->qwgt[ls][nis-1];
      e3   = data->qxg[lt][nit-1];
      fact = data->qwgt[lt][nit-1];
      f3_hex(funct,deriv,deriv2,e1,e2,e3,typ,icode);
      break;
   case tet4: case tet10:   /* --> tet - element */		  	
      e1   = data->txgr[lr][intc];
      facr = data->twgt[lr][intc];
      e2   = data->txgs[lr][intc];
      facs = ONE;
      e3   = data->txgt[lr][intc]; 
      fact = ONE;
      f3_tet(funct,deriv,deriv2,e1,e2,e3,typ,icode); 
      break;
   default:
      dserror("typ unknown!");
   } /* end switch (typ) */
/*-------------------------------------------- compute Jacobian matrix */  
   f3_mljaco(funct,deriv,xjm,&det,ele,iel);
   fac = facr*facs*fact*det;
/*------------------------------------------- compute global derivates */
   f3_gder(derxy,deriv,xjm,wa1,det,iel);

/*------------------------ get velocities (n+1,i) at integration point */
   f3_veci(velint,funct,evel,iel);
/*-------------- get velocity derivatives (n+1,i) at integration point */
   f3_vder(vderxy,derxy,evel,iel);
 
/*---- compute stab. par. or subgrid viscosity during integration loop */
      if (gls->iduring!=0 || 
          (fdyn->sgvisc>0    && gls->iduring!=0))
        f3_mlcalelesize2(ele,velint,vderxy,derxy,visc,iel,typ);
/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrices                         |
 | NOTE:                                                                |
 |  Standard Galerkin matrices are all stored in one matrix "estif"     |
 |  Standard Galerkin mass matrix is stored in "emass"                  |
 *----------------------------------------------------------------------*/
/*------------------------ compute standard Galerkin part of matrix Kvv */      
   if (mlvar->convel!=0) 
      f3_lscalkvv(estif,velint,vderxy,funct,derxy,fac,visc,iel);
/*-------------- compute standard Galerkin part of matrices Kvp and Kpv */      
   f3_lscalkvp(estif,funct,derxy,fac,iel);
/*------------------------ compute standard Galerkin part of matrix Mvv */      
   if (fdyn->nis==0) f3_lscalmvv(emass,funct,fac,iel);
   
/*----------------------------------------------------------------------*
 |         compute Stabilization matrices                               |
 |  (in most of the cases only the bulk viscosity term is calculated)   |
 | NOTE:                                                                |
 |  Stabilization matrices are all stored in one matrix "estif"         |
 |  Stabilization mass matrices are all stored in one matrix "emass"    |
 *----------------------------------------------------------------------*/
/*------------------------------------ compute second global derivative */ 
   if (ihoel!=0) f3_mlgder2(ele,xjm,wa2,derxy,derxy2,deriv2,iel);
   
/*---------------------------------------- stabilization for matrix Kvv */
   f3_lscalstabkvv(ele,estif,velint,vderxy,
                          funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilization for matrix Kvp */
   f3_lscalstabkvp(ele,estif,velint,vderxy,
                       funct,derxy,derxy2,fac,visc,iel,ihoel); 
/*---------------------------------------- stabilization for matrix Mvv */
   if (fdyn->nis==0) 
      f3_lscalstabmvv(ele,emass,velint,vderxy,
	                     funct,derxy,derxy2,fac,visc,iel,ihoel);
   if (gls->ipres!=0)	        
   {
/*---------------------------------------- stabilization for matrix Kpv */ 
      f3_lscalstabkpv(estif,velint,vderxy,
	                  funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilization for matrix Mpv */
      if (fdyn->nis==0)
         f3_lscalstabmpv(emass,funct,derxy,fac,iel);
   } /* endif (ele->e.f3->ipres!=0) */
   /*---------------------------------------- stabilization for matrix Kpp */ 
   if (gls->ipres!=0)
      f3_lscalstabkpp(estif,derxy,fac,iel);  

/*----------------------------------------------------------------------*
 |         compute "Iteration" Force Vectors                            |
 *----------------------------------------------------------------------*/ 
   if (fdyn->nii!=0 && mlvar->convel!=0)
   {
/*-------------- get convective velocities (n+1,i) at integration point */
     f3_covi(vderxy,velint,covint);
/*- calculate "Iteration" force vector for velocity dofs (no pre. dofs) */
     f3_lscalgalifv(eiforce,covint,velint,vderxy,funct,fac,iel);
   } /* endif (fdyn->nii!=0) */

/*----------------------------------------------------------------------*
 |         compute "Time" Force Vector                                  |
 *----------------------------------------------------------------------*/
    if (fdyn->nis==0)
    {
 /*------------------------------ get pressure (n) at integration point */
      if (fdyn->iprerhs>0) 
         preintn=f3_scali(funct,epren,iel);
/*----------------------------- get velocities (n) at integration point */
      f3_veci(velintn,funct,eveln,iel);
/*------------------- get velocity derivatives (n) at integration point */
      f3_vder(vderxyn,derxy,eveln,iel);
/*------------------ get convective velocities (n) at integration point */
      f3_covi(vderxyn,velintn,covintn);        	    
/*--------------------- calculate "Time" force vector for velocity dofs */
      f3_lscalgaltfv(etforce,velintn,velintn,covintn,funct,derxy,
	             vderxyn,preintn,visc,fac,iel);		       
/*--------------------- calculate "Time" force vector for pressure dofs */
      if (fdyn->thsr!=ZERO) 
	f3_lscalgaltfp(&(etforce[3*iel]),funct,vderxyn,fac,iel);
      }  

/*----------------------------------------------------------------------*
 |         compute "External" Force Vector                              |
 *----------------------------------------------------------------------*/
    if (*hasext!=0)
/*----------------- calculate "External" force vector for velocity dofs */
      f3_lscalgalexfv(etforce,funct,edeadn,edead,fac,iel);
}
}
} /* end of loop over integration points */
 
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_lsint */

/*!---------------------------------------------------------------------
\brief integration loop for dynamic procedure for sm ele. for fluid3

<pre>                                                       gravem 07/03

In this routine, the integral of the diffusive part of the lhs as well
as the normalized rhs are determined.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *submesh     FLUID_ML_SMESH(i)   
\param **smiediff     DOUBLE       (o)	sm ele. int. of diffusive lhs
\param  *smierhs      DOUBLE       (o)	sm ele. int. of rhs
\param  *smxyze       DOUBLE       (i)	submesh element coordinates
\param  *smfunct      DOUBLE       (-)	sm natural shape functions
\param **smderiv      DOUBLE       (-)  sm deriv. of nat. shape funcs
\param **smderiv2     DOUBLE       (-)  sm 2nd deriv. of nat. shape f.
\param **smxjm	      DOUBLE       (-)  sm jacobian matrix
\param **smderxy      DOUBLE       (-)  sm global derivatives
\return void                                                   

------------------------------------------------------------------------*/
void f3_smint2(FLUID_DATA      *data,     
	       ELEMENT	       *ele,	
	       FLUID_ML_SMESH  *submesh, 
               DOUBLE	      **smediff,   
	       DOUBLE	       *smerhs,   
	       DOUBLE	      **smxyze, 
	       DOUBLE	       *smfunct,   
	       DOUBLE	      **smderiv,   
	       DOUBLE	      **smderiv2,  
	       DOUBLE	      **smxjm,	  
	       DOUBLE	      **smderxy,
	       DOUBLE	      **wa1)
{ 
INT       smiel;      /* submesh number of nodes                        */
INT       nsmtyp;     /* submesh element type: 1 - hex; 2 - tri        */
INT       intc;       /* "integration case" for tet for further infos
                          see f3_inpele.c and f3_intg.c                 */
INT       nir,nis,nit;/* number of integration nodesin r,s direction    */
INT       lr,ls,lt;   /* counter for integration                        */
DOUBLE    fac;
DOUBLE    facr,facs,fact;/* integration weights                         */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2,e3;   /* natural coordinates of integr. point           */
DIS_TYP   typ,smtyp;  /* large-scale and submesh element type           */

#ifdef DEBUG 
dstrc_enter("f3_smint2");
#endif

/*------------------------------------------------------ initialization */
smiel  = submesh->numen;
nsmtyp = submesh->ntyp; 
smtyp  = submesh->typ;

/*--------------------------- get integration data for submesh elements */
switch (nsmtyp)
{
case 1:  /* --> hex - element */
   nir = submesh->ngpr;
   nis = submesh->ngps;
   nit = submesh->ngpt;
break;
case 2: /* --> tet - element */  
   nir  = submesh->ngpr;
   nis  = 1;
   nit  = 1;
   intc = submesh->ngps;  
break;
default:
   dserror("nsmtyp unknown!");
} /* end switch(nsmtyp) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{    
   for (ls=0;ls<nis;ls++)
   {
   for (lt=0;lt<nit;lt++)
   {
/*--------- get values of submesh shape functions and their derivatives */
      switch(nsmtyp)  
      {
      case 1:   /* --> hex - element */
	 e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
         e3   = data->qxg[lt][nit-1];
         fact = data->qwgt[lt][nit-1];
         f3_hex(smfunct,smderiv,smderiv2,e1,e2,e3,smtyp,2);
      break;
      case 2:   /* --> tri - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
         e3   = data->txgt[lr][intc]; 
         fact = ONE;
	 f3_tet(smfunct,smderiv,smderiv2,e1,e2,e3,smtyp,2);
      break;
      default:
         dserror("nsmtyp unknown!");
      } /* end switch(nsmtyp) */
/*--------------------------------- compute Jacobian matrix for submesh */
      f3_mljaco3(smxyze,smfunct,smderiv,smxjm,&det,smiel,ele);
      fac = facr*facs*fact*det;
/*-------------------------------- compute global derivates for submesh */
      f3_gder(smderxy,smderiv,smxjm,wa1,det,smiel);
      
/*-------------- compute diffusive part of submesh stiffness matrix SMK */      
      f3_calsmkd(smediff,smderxy,fac,smiel);
/*---------------------------- compute integral of normalized RHS force */
      fluid_calnofo(smerhs,smfunct,fac,smiel); 
  }
  }
}
        
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_smint2 */

/*!---------------------------------------------------------------------
\brief integration loop for sub-submesh element for fluid2

<pre>                                                       gravem 07/03

In this routine, the element stiffness matrix and normalized RHS for 
one sub-submesh element is calculated.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *mlvar       FLUID_DYN_ML  (i)
\param  *ssmesh      FLUID_ML_SMESH(i)   
\param **ssestif      DOUBLE       (o)	ssm element stiffness matrix
\param  *ssenfor      DOUBLE       (o)	ssm element norm. force vector
\param  *ssxyze       DOUBLE       (i)	sub-submesh element coordinates
\param  *ssxyzep      DOUBLE       (i)	ssm ele. coord. on parent dom.
\param  *funct        DOUBLE       (-)	natural shape functions
\param **deriv        DOUBLE       (-)	deriv. of nat. shape funcs
\param **deriv2       DOUBLE       (-)	2nd deriv. of nat. shape f.
\param **xjm	      DOUBLE       (-)	jacobian matrix
\param **derxy        DOUBLE       (-)	global derivatives
\param  *ssfunct      DOUBLE       (-)	sm natural shape functions
\param **ssderiv      DOUBLE       (-)  sm deriv. of nat. shape funcs
\param **ssderiv2     DOUBLE       (-)  sm 2nd deriv. of nat. shape f.
\param **ssxjm	      DOUBLE       (-)  sm jacobian matrix
\param **ssderxy      DOUBLE       (-)  sm global derivatives
\param **evel         DOUBLE       (i)  ele vel. at time step n+1
\param  *velint       DOUBLE       (-)  vel at integration point
\param **vderxy       DOUBLE       (-)  global vel. derivatives
\param **wa1          DOUBLE       (-)  working array
\return void                                                   

------------------------------------------------------------------------*/
void f3_ssint(FLUID_DATA      *data,     
	      ELEMENT	      *ele,
	      FLUID_DYN_ML    *mlvar, 
	      FLUID_ML_SMESH  *ssmesh, 
              DOUBLE	     **ssestif,   
	      DOUBLE	      *ssenfor,   
	      DOUBLE	     **ssxyze, 
	      DOUBLE	     **ssxyzep, 
	      DOUBLE	      *funct,	
	      DOUBLE	     **deriv,	
	      DOUBLE	     **deriv2,  
	      DOUBLE	     **xjm,	
	      DOUBLE	     **derxy,	
	      DOUBLE	      *ssfunct,   
	      DOUBLE	     **ssderiv,   
	      DOUBLE	     **ssderiv2,  
	      DOUBLE	     **ssxjm,	  
	      DOUBLE	     **ssderxy,   
	      DOUBLE	     **evel,  
	      DOUBLE	      *velint,  
	      DOUBLE	     **vderxy,
	      DOUBLE	     **wa1)
{ 
INT       iel,ssiel;  /* large-scale and sub-submesh number of nodes    */
INT       nsstyp;     /* l-s and ssm element type: 1 - hex; 2 - tri    */
INT       intc;       /* "integration case" for tet for further infos
                          see f3_inpele.c and f3_intg.c                 */
INT       nir,nis,nit;/* number of integration nodesin r,s direction    */
INT       actmat;     /* material number of the element                 */
INT       lr,ls,lt;   /* counter for integration                        */
DOUBLE    visc;       /* viscosity                                      */
DOUBLE    fac;
DOUBLE    facr,facs,fact;/* integration weights                         */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2,e3;   /* natural coordinates of integr. point           */
DOUBLE    coor[3];    /* coordinates                                    */
DIS_TYP   typ,sstyp;  /* large-scale and sub-submesh element type       */

#ifdef DEBUG 
dstrc_enter("f3_ssint");
#endif

/*------------------------------------------------------ initialization */
iel    = ele->numnp;
ssiel  = ssmesh->numen;
actmat = ele->mat-1;
visc   = mat[actmat].m.fluid->viscosity;
typ    = ele->distyp;
nsstyp = ssmesh->ntyp; 
sstyp  = ssmesh->typ;

/*---------------------- get integration data for sub-submesh elements */
switch (nsstyp)
{
case 1:  /* --> hex - element */
   nir = ssmesh->ngpr;
   nis = ssmesh->ngps;
   nit = ssmesh->ngpt;
break;
case 2: /* --> tet - element */  
   nir  = ssmesh->ngpr;
   nis  = 1;
   nit  = 1;
   intc = ssmesh->ngps;  
break;
default:
   dserror("nsstyp unknown!");
} /* end switch(nsstyp) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{    
   for (ls=0;ls<nis;ls++)
   {
   for (lt=0;lt<nit;lt++)
   {
/*---- get values of sub-submesh shape functions and their derivatives */
      switch(nsstyp)  
      {
      case 1:   /* --> hex - element */
	 e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
	 e3   = data->qxg[ls][nis-1];
	 fact = data->qwgt[ls][nis-1];
         f3_hex(ssfunct,ssderiv,ssderiv2,e1,e2,e3,sstyp,2);
      break;
      case 2:   /* --> tet - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
	 f3_tet(ssfunct,ssderiv,ssderiv2,e1,e2,e3,sstyp,2);
      break;
      default:
         dserror("nsstyp unknown!");
      } /* end switch(nsstyp) */
/*---------------------------- compute Jacobian matrix for sub-submesh */
      f3_mljaco3(ssxyze,ssfunct,ssderiv,ssxjm,&det,ssiel,ele);
      fac = facr*facs*fact*det;
/*--------------------------- compute global derivates for sub-submesh */
      f3_gder(ssderxy,ssderiv,ssxjm,wa1,det,ssiel);

/*- get coordinates of current int. p. on parent domain of l-s element */
      f3_mlgcoor2(ssfunct,ssxyzep,ssiel,coor);
      
/*---- get values of large-scale shape functions and their derivatives */
      switch(typ)
      {
      case hex8: case hex20: case hex27:    /* --> hex - element */
        f3_hex(funct,deriv,deriv2,coor[0],coor[1],coor[2],typ,2);
      break;
      case tet4: case tet10:	/* --> tet - element */ 	     
        f3_tet(funct,deriv,deriv2,coor[0],coor[1],coor[2],typ,2);   
      break;
      default:
        dserror("typ unknown!\n");      
      } /*end switch(typ) */
/*------------ compute Jacobian matrix for large-scale shape functions */
      f3_mljaco(funct,deriv,xjm,&det,ele,iel);
/*----------- compute global derivates for large-scale shape functions */
      f3_gder(derxy,deriv,xjm,wa1,det,iel);

/*------------ get large-scale velocities (n+1,i) at integration point */
      f3_veci(velint,funct,evel,iel);
/*-- get large-scale velocity (n+1,i) derivatives at integration point */
      f3_vder(vderxy,derxy,evel,iel);
      
/*----------------------------------------------------------------------*
 |       compute matrix and rhs for sub-submesh                         |
 *----------------------------------------------------------------------*/
/*------------------------------------------------  compute matrix SSK */      
      f3_calsmk(mlvar,ssestif,velint,vderxy,ssfunct,ssderxy,fac,
                visc,ssiel);
/*--------------------------------------------- compute normalized rhs */
      fluid_calnofo(ssenfor,ssfunct,fac,ssiel); 
   } /* end of loop over integration points lt*/
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_ssint */

/*!---------------------------------------------------------------------
\brief integration loop for norm. bubble on sub-submesh ele. for fluid3

<pre>                                                       gravem 07/03

In this routine, the integral of the normalized bubble function for 
one sub-submesh element is calculated.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *ssmesh      FLUID_ML_SMESH(i)   
\param  *ssinbu       DOUBLE	   (o)  sub-submesh bubble integral
\param  *ebub         DOUBLE       (i)	ssm element bubble function
\param **ssxyze       DOUBLE       (i)	sub-submesh element coordinates
\param  *ssfunct      DOUBLE       (-)	sm natural shape functions
\param **ssderiv      DOUBLE       (-)  sm deriv. of nat. shape funcs
\param **ssderiv2     DOUBLE       (-)  sm 2nd deriv. of nat. shape f.
\param **ssxjm	      DOUBLE       (-)  sm jacobian matrix
\return void                                                   

------------------------------------------------------------------------*/
void f3_inbu(FLUID_DATA      *data,     
	     ELEMENT	     *ele,   
	     FLUID_ML_SMESH  *ssmesh, 
             DOUBLE	     *ssinbu,   
             DOUBLE	     *ebub,   
	     DOUBLE	    **ssxyze, 
	     DOUBLE	     *ssfunct,   
	     DOUBLE	    **ssderiv,   
	     DOUBLE	    **ssderiv2,  
	     DOUBLE	    **ssxjm)
{ 
INT       i;          /* simply a counter                               */
INT       ssiel;      /* sub-submesh number of nodes                    */
INT       nsstyp;     /* ssm element type: 1 - hex; 2 - tet             */
INT       intc;       /* "integration case" for tet for further infos
                          see f3_inpele.c and f3_intg.c                 */
INT       nir,nis,nit;/* number of integration nodesin r,s direction    */
INT       lr,ls,lt;   /* counter for integration                        */
DOUBLE    fac;
DOUBLE    facr,facs,fact;/* integration weights                         */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2,e3;   /* natural coordinates of integr. point           */
DOUBLE    bubnint;    /* normalized bubble function at int. point       */
DIS_TYP   sstyp;      /* sub-submesh element type                       */

#ifdef DEBUG 
dstrc_enter("f3_inbu");
#endif

/*------------------------------------------------------ initialization */
ssiel  = ssmesh->numen;
nsstyp = ssmesh->ntyp; 
sstyp  = ssmesh->typ;

/*----------------------- get integration data for sub-submesh elements */
switch (nsstyp)
{
case 1:  /* --> hex - element */
   nir = ssmesh->ngpr;
   nis = ssmesh->ngps;
   nit = ssmesh->ngpt;
break;
case 2: /* --> tet - element */  
   nir  = ssmesh->ngpr;
   nis  = 1;
   nit  = 1;
   intc = ssmesh->ngps;  
break;
default:
   dserror("nsstyp unknown!");
} /* end switch(nsstyp) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{    
   for (ls=0;ls<nis;ls++)
   {
   for (lt=0;lt<nit;lt++)
   {
/*---- get values of sub-submesh shape functions and their derivatives */
      switch(nsstyp)  
      {
      case 1:   /* --> hex - element */
	 e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
         e3   = data->qxg[lt][nit-1];
         fact = data->qwgt[lt][nit-1];
         f3_hex(ssfunct,ssderiv,ssderiv2,e1,e2,e3,sstyp,2);
      break;
      case 2:   /* --> tet - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
	 f3_tet(ssfunct,ssderiv,ssderiv2,e1,e2,e3,sstyp,2);
      break;
      default:
         dserror("nsstyp unknown!");
      } /* end switch(nsstyp) */
/*---------------------------- compute Jacobian matrix for sub-submesh */
      f3_mljaco3(ssxyze,ssfunct,ssderiv,ssxjm,&det,ssiel,ele);
      fac = facr*facs*fact*det;

/*------------ compute normalized bubble function at integration point */
      bubnint=ZERO; 
      for (i=0;i<ssiel;i++) /* loop all nodes of the ssm element */
      {
        bubnint += ssfunct[i]*ebub[i];
      } /* end loop all nodes of the ssm element */
      
/*----- compute sub-submesh integral of the normalized bubble function */
      *ssinbu += bubnint*fac; 
   } /* end of loop over integration points lt*/
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_inbu */

#endif
