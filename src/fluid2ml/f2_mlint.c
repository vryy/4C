/*!----------------------------------------------------------------------
\file
\brief integration routines for multi-level fluid2

<pre>
Maintainer:  name
             e-mail
             homepage
             telephone number


</pre>

------------------------------------------------------------------------*/
#ifdef FLUID2_ML
#include "../headers/standardtypes.h"
#include "../fluid_full/fluid_prototypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2ml_prototypes.h"
#include "../fluid2/fluid2.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!---------------------------------------------------------------------
\brief integration loop for submesh element for fluid2

<pre>                                                       gravem 07/03

In this routine, the element stiffness matrix, mass matrix, VMM-RHS and
Time-RHS for one submesh element is calculated.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *dynvar      FLUID_DYN_CALC(i)
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
void f2_smint(FLUID_DATA      *data,     
	      ELEMENT	      *ele,	
	      FLUID_DYN_CALC  *dynvar, 
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
INT       ntyp,nsmtyp;/* l-s and submesh element type: 1-quad; 2-tri    */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis;    /* number of integration points in r,s direction  */
INT       actmat;     /* material number of the element                 */
INT       ihoel=0;    /* flag for higher order large-scale elements     */
INT       ihoelsm=0;  /* flag for higher order submesh elements         */
INT       icode=2;    /* flag for eveluation of l-s shape functions     */     
INT       icodesm=2;  /* flag for eveluation of sm shape functions      */     
INT       lr,ls;      /* counter for integration                        */
DOUBLE    dens;       /* density                                        */
DOUBLE    visc;       /* viscosity                                      */
DOUBLE    fac;
DOUBLE    facr, facs; /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DOUBLE    coor[2];    /* coordinates                                    */
DIS_TYP   typ,smtyp;  /* large-scale and submesh element type           */

#ifdef DEBUG 
dstrc_enter("f2_smint");
#endif

/*----------------------------------------------------- initialization */
iel    = ele->numnp;
smiel  = submesh->numen;
actmat = ele->mat-1;
dens   = mat[actmat].m.fluid->density;
visc   = mat[actmat].m.fluid->viscosity;
ntyp   = ele->e.f2->ntyp; 
typ    = ele->distyp;
nsmtyp = submesh->ntyp; 
smtyp  = submesh->typ;

/*--- get integration data and check if sm-elements are "higher order" */
switch (nsmtyp)
{
case 1:  /* --> quad - element */
   if (smiel>4)
   {
     icodesm = 3;
     ihoelsm = 1;
   }  
   /* initialize integration */
   nir = submesh->ngpr;
   nis = submesh->ngps;
break;
case 2: /* --> tri - element */  
   if (smiel>3)
   {
      icodesm = 3;
      ihoelsm = 1;
   }
   /* initialize integration */
   nir  = submesh->ngpr;
   nis  = 1;
   intc = submesh->ngps;  
break;
default:
   dserror("nsmtyp unknown!");
} /* end switch(nsmtyp) */

/*---------------------------- check if ls-elements are "higher order" */
switch (ntyp)
{
case 1:  /* --> quad - element */
   icode   = 3;
   ihoel   = 1;
break;
case 2: /* --> tri - element */  
   if (iel>3)
   {
     icode   = 3;
     ihoel   = 1;
   }
break;
default:
   dserror("ntyp unknown!");
} /* end switch(ntyp) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{    
   for (ls=0;ls<nis;ls++)
   {
/*-------- get values of submesh shape functions and their derivatives */
      switch(nsmtyp)  
      {
      case 1:   /* --> quad - element */
	 e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
         f2_rec(smfunct,smderiv,smderiv2,e1,e2,smtyp,icodesm);
      break;
      case 2:   /* --> tri - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
	 f2_tri(smfunct,smderiv,smderiv2,e1,e2,smtyp,icodesm);
      break;
      default:
         dserror("nsmtyp unknown!");
      } /* end switch(nsmtyp) */
/*-------------------------------- compute Jacobian matrix for submesh */
      f2_mljaco3(smxyze,smfunct,smderiv,smxjm,&det,smiel,ele);
      fac = facr*facs*det;
/*------------------------------- compute global derivates for submesh */
      f2_gder(smderxy,smderiv,smxjm,det,smiel);

/*- get coordinates of current int. p. on parent domain of l-s element */
      f2_mlgcoor2(smfunct,smxyzep,smiel,coor);

/*---- get values of large-scale shape functions and their derivatives */
      switch(ntyp)
      {
      case 1:    /* --> quad - element */
        f2_rec(funct,deriv,deriv2,coor[0],coor[1],typ,icode);
      break;
      case 2:	/* --> tri - element */ 	     
        f2_tri(funct,deriv,deriv2,coor[0],coor[1],typ,icode);   
      break;
      default:
        dserror("ntyp unknown!\n");      
      } /*end switch(ntyp) */
/*------------ compute Jacobian matrix for large-scale shape functions */
      f2_mljaco(funct,deriv,xjm,&det,ele,iel);
/*----------- compute global derivates for large-scale shape functions */
      f2_gder(derxy,deriv,xjm,det,iel);
/*------- compute 2nd global derivates for large-scale shape functions */
      if (ihoel!=0) f2_mlgder2(ele,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);

/*------------ get large-scale velocities (n+1,i) at integration point */
      f2_veli(velint,funct,evel,iel);
/*-- get large-scale velocity (n+1,i) derivatives at integration point */
      f2_vder(vderxy,derxy,evel,iel);
      
      if (mlvar->convel==0)
      { 
/*-----------------------------------------------------------------------
     get values of bubble functions and their derivatives 
----------------------------------------------------------------------- */
/*------------------- get velocity bubble functions at integraton point */
        fluid_bubint (vbubint,smfunct,evbub,smiel,mlvar->nvbub);            
/*-------- get velocity bubble function derivatives at integraton point */
        fluid_bubder (vbubderxy,smderxy,evbub,smiel,mlvar->nvbub,2);            
/*----------------- get 'pressure' bubble functions at integraton point */
        fluid_pbubint (pbubint,smfunct,epbub,smiel,iel,2);            
/*------ get 'pressure' bubble function derivatives at integraton point */
        fluid_pbubder (pbubderxy,smderxy,epbub,smiel,iel,2,2);            
/*--------------------------- get small-scale 'rhs' at integraton point */
        fluid_bubint (smfint,smfunct,efbub,smiel,2);            
/*--------------- get small-scale 'rhs' derivatives at integraton point */
        fluid_bubder (smfderxy,smderxy,efbub,smiel,2,2);            
	 
/*---------------------- get small-scale velocities at integraton point */
        f2_veli (smvelint,vbubint,evel,iel);	      
/*------------ get small-scale velocity derivatives at integraton point */
        f2_vder (smvderxy,vbubderxy,evel,iel);	      
/*--------------------- get small-scale 'pressures' at integraton point */
        f2_smprei (smpreint,pbubint,epre,iel);	      
/*---------- get small-scale 'pressure' derivatives at integraton point */
        f2_smpder (smpderxy,pbubderxy,epre,iel);	      

/*--- calculate velocities and velocity derivatives at integraton point */
        for (i=0; i<2; i++)
	{
	  velint[i] += smvelint[i] + smpreint[i] + smfint[i];
	  for (j=0; j<2; j++)
	  {
	    vderxy[i][j] += smvderxy[i][j] + smpderxy[i][j] + smfderxy[j][i];
	  }
	}
      }	    

/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrices for submesh             |
 *----------------------------------------------------------------------*/
/*-------------------------------------------------  compute matrix SMK */      
      f2_calsmk(dynvar,mlvar,smestif,velint,vderxy,smfunct,smderxy,fac,
                visc,smiel);
/*-------------------------------------------------- compute matrix SMM */
      if (dynvar->nis==0 && mlvar->quastabub==0) 
        f2_calsmm(smemass,smfunct,fac,smiel);
      
/*----------------------------------------------------------------------*
 |         compute Stabilization matrices for submesh                   |
 *----------------------------------------------------------------------*/
      if (mlvar->smstabi>0)
      { 
/*----------------------- compute second global derivatives for submesh */ 
        if (ihoelsm!=0) 
	  f2_mlcogder2(smxyze,smxjm,wa1,wa2,smderxy,smderxy2,smderiv2,smiel);
/*---------------------------------------- stabilization for matrix SMK */
        f2_calstabsmk(dynvar,mlvar,smestif,velint,vderxy,smfunct,smderxy,
                      smderxy2,fac,visc,smiel,ihoelsm);
/*---------------------------------------- stabilization for matrix SMM */
        if (dynvar->nis==0 && mlvar->quastabub==0) 
	  f2_calstabsmm(dynvar,mlvar,smemass,velint,vderxy,smfunct,smderxy,
	                smderxy2,fac,visc,smiel,ihoelsm); 
      } /* endif (mlvar->smstabi>0) */

/*----------------------------------------------------------------------*
 |         compute "VMM" Force Vectors                                  |
 *----------------------------------------------------------------------*/ 
/*----------------------- standard Galerkin part for "VMM" force vector */
      f2_calsmfv(dynvar,mlvar,smevfor,velint,vderxy,smfunct,funct,derxy,
                 derxy2,fac,visc,smiel,iel,ihoel); 
/*--------------------------- stabilization part for "VMM" force vector */
      if (mlvar->smstabi>0)
        f2_calstabsmfv(dynvar,mlvar,smevfor,velint,vderxy,smfunct,smderxy,
	               smderxy2,funct,derxy,derxy2,fac,visc,smiel,iel,
		       ihoelsm,ihoel);

/*----------------------------------------------------------------------*
 |         compute "Time" Force Vector                                  |
 *----------------------------------------------------------------------*/
/*------- cancel "Time" force vector calculation for certain parameters */
      if (mlvar->transterm>1 && dynvar->thsr==ZERO && mlvar->quastabub!=0) 
        continue;
	
      if (dynvar->nis==0)
      {
/*------- get large-scale pressure derivatives (n) at integration point */
        if (dynvar->iprerhs>0) f2_pder(pderxyn,derxy,epren,iel);
/*------------------ get large-scale velocities (n) at integraton point */
        f2_veli(velintn,funct,eveln,iel);
/*------- get large-scale velocity derivatives (n) at integration point */
        f2_vder(vderxyn,derxy,eveln,iel);
/*--- get large-scale 2nd velocity derivatives (n) at integration point */
        if (ihoel!=0) f2_vder2(vderxy2n,derxy2,eveln,iel);
      
        if (mlvar->quastabub==0)
        { 
/*--------------- get velocity bubble functions (n) at integraton point */
          fluid_bubint (vbubintn,smfunct,evbubn,smiel,mlvar->nvbub);            
/*---- get velocity bubble function derivatives (n) at integraton point */
          fluid_bubder (vbubderxyn,smderxy,evbubn,smiel,mlvar->nvbub,2);            
/* get 2nd velocity bubble function derivatives (n) at integraton point */
          if (ihoelsm!=0) 
	  fluid_bubder (vbubderxy2n,smderxy2,evbubn,smiel,mlvar->nvbub,3);
	              
/*------------- get 'pressure' bubble functions (n) at integraton point */
          fluid_pbubint (pbubintn,smfunct,epbubn,smiel,iel,2);            
/*-- get 'pressure' bubble function derivatives (n) at integraton point */
          fluid_pbubder (pbubderxyn,smderxy,epbubn,smiel,iel,2,2);            
/*--- get 2nd 'pressure' bubble function deriv. (n) at integraton point */
          if (ihoelsm!=0) 
	  fluid_pbubder (pbubderxy2n,smderxy2,epbubn,smiel,iel,2,3);
	  
/*----------------------- get small-scale 'rhs' (n) at integraton point */
          fluid_bubint (smfintn,smfunct,efbubn,smiel,2);            
/*----------- get small-scale 'rhs' derivatives (n) at integraton point */
          fluid_bubder (smfderxyn,smderxy,efbubn,smiel,2,2);            
/*------- get small-scale 2nd 'rhs' derivatives (n) at integraton point */
          if (ihoelsm!=0) 
	  fluid_bubder (smfderxy2n,smderxy2,efbubn,smiel,2,3);
	   
/*------------------ get small-scale velocities (n) at integraton point */
          f2_veli (smvelintn,vbubintn,eveln,iel);	      
/*-------- get small-scale velocity derivatives (n) at integraton point */
          f2_vder (smvderxyn,vbubderxyn,eveln,iel);	      
/*---- get small-scale 2nd velocity derivatives (n) at integraton point */
          if (ihoelsm!=0) f2_vder2(smvderxy2n,vbubderxy2n,eveln,iel);
	  	      
/*----------------- get small-scale 'pressures' (n) at integraton point */
          f2_smprei (smpreintn,pbubintn,epren,iel);	      
/*------ get small-scale 'pressure' derivatives (n) at integraton point */
          f2_smpder (smpderxyn,pbubderxyn,epren,iel);	      
/*-- get small-scale 2nd 'pressure' derivatives (n) at integraton point */
          if (ihoelsm!=0) f2_smpder2(smpderxy2n,pbubderxy2n,epren,iel);
	}    
	  
/*----- calculate various velocities for "Time" force vector evaluation */
        for (i=0; i<2; i++)
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
	    for (j=0; j<2; j++)
	    {
/*----- velocity derivatives for convective part of "Time" force vector */
	      if (mlvar->convel==0) vderxync[i][j] = vderxyn[i][j]\
	          + smvderxy[i][j] + smpderxy[i][j] + smfderxy[j][i];
	      else                  vderxync[i][j] = vderxyn[i][j];
	    }
/*---- 2nd velocity derivatives for viscous part of "Time" force vector */
	    if (ihoel!=0 && ihoelsm!=0)
	    {
	      for (j=0; j<3; j++)
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
	    for (j=0; j<2; j++)
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
	      for (j=0; j<3; j++)
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
        f2_calsmft(dynvar,mlvar,smetfor,velintn,velintnt,velintnc,vderxyn, 
		   vderxync,vderxynv,vderxy2nv,pderxyn,smfunct,smderxy,fac,
		   visc,smiel,iel,ihoelsm,ihoel);
/*-------------------------- stabilization part for "Time" force vector */
        if (mlvar->smstabi>0)
          f2_calstabsmft(dynvar,mlvar,smetfor,velintn,velintnt,velintnc,
	                 vderxyn,vderxync,vderxynv,vderxy2n,pderxyn,smfunct,
			 smderxy,smderxy2,fac,visc,smiel,iel,ihoelsm,ihoel);
      }  
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_smint */

/*!---------------------------------------------------------------------
\brief integration loop for bubble funct. on submesh element for fluid2

<pre>                                                       gravem 07/03

In this routine, the bubble function part of the large-scale element 
stiffness matrix, mass matrix, VMM-RHS and Time-RHS for one submesh 
element is calculated.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *dynvar      FLUID_DYN_CALC(i)
\param  *mlvar       FLUID_DYN_ML  (i)
\param  *submesh     FLUID_ML_SMESH(i)   
\param **estif        DOUBLE	   (o)  element stiffness matrix
\param **emass        DOUBLE	   (o)  element mass matrix
\param  *eiforce      DOUBLE	   (o)  element iteration force vector
\param  *smxyze	      DOUBLE	   (i)  submesh element coordinates
\param  *smxyzep      DOUBLE	   (i)  sm ele. coord. on parent dom.
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
void f2_bubint(FLUID_DATA      *data,     
	       ELEMENT         *ele,     
	       FLUID_DYN_CALC  *dynvar, 
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
INT       ntyp,nsmtyp;/* l-s and submesh element type: 1-quad; 2-tri    */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis;    /* number of integration nodesin r,s direction    */
INT       actmat;     /* material number of the element                 */
INT       ihoel=0;    /* flag for higher order large-scale elements     */
INT       ihoelsm=0;  /* flag for higher order submesh elements         */
INT       icode=2;    /* flag for eveluation of l-s shape functions     */     
INT       icodesm=2;  /* flag for eveluation of sm shape functions      */     
INT       lr,ls;      /* counter for integration                        */
DOUBLE    dens;       /* density                                        */
DOUBLE    visc;       /* viscosity                                      */
DOUBLE    fac;
DOUBLE    facr,facs;  /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DOUBLE    coor[2];    /* coordinates                                    */
DIS_TYP   typ,smtyp;  /* large-scale and submesh element type           */

#ifdef DEBUG 
dstrc_enter("f2_bubint");
#endif

/*----------------------------------------------------- initialization */
iel    = ele->numnp;
smiel  = submesh->numen;
actmat = ele->mat-1;
dens   = mat[actmat].m.fluid->density;
visc   = mat[actmat].m.fluid->viscosity;
ntyp   = ele->e.f2->ntyp; 
typ    = ele->distyp;
nsmtyp = submesh->ntyp; 
smtyp  = submesh->typ;

/*---- get integraton data and check if sm-elements are "higher order" */
switch (nsmtyp)
{
case 1:  /* --> quad - element */
   if (smiel>4)
   {
     icodesm = 3;
     ihoelsm = 1;
   }  
   /* initialize integration */
   nir = submesh->ngpr;
   nis = submesh->ngps;
break;
case 2: /* --> tri - element */  
   if (smiel>3)
   {
      icodesm = 3;
      ihoelsm = 1;
   }
   /* initialize integration */
   nir  = submesh->ngpr;
   nis  = 1;
   intc = submesh->ngps;  
break;
default:
   dserror("nsmtyp unknown!");
} /* end switch(nsmtyp) */

/*---------------------------- check if ls-elements are "higher order" */
switch (ntyp)
{
case 1:  /* --> quad - element */
   icode   = 3;
   ihoel   = 1;
break;
case 2: /* --> tri - element */  
   if (iel>3)
   {
      icode   = 3;
      ihoel   = 1;
   }
break;
default:
   dserror("ntyp unknown!");
} /* end switch(ntyp) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{    
   for (ls=0;ls<nis;ls++)
   {
/*-------- get values of submesh shape functions and their derivatives */
      switch(nsmtyp)  
      {
      case 1:   /* --> quad - element */
	 e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
         f2_rec(smfunct,smderiv,smderiv2,e1,e2,smtyp,icodesm);
      break;
      case 2:   /* --> tri - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
	 f2_tri(smfunct,smderiv,smderiv2,e1,e2,smtyp,icodesm);
      break;
      default:
         dserror("nsmtyp unknown!");
      } /* end switch(nsmtyp) */
/*-------------------------------- compute Jacobian matrix for submesh */
      f2_mljaco3(smxyze,smfunct,smderiv,smxjm,&det,smiel,ele);
      fac = facr*facs*det;
/*------------------------------- compute global derivates for submesh */
      f2_gder(smderxy,smderiv,smxjm,det,smiel);

/*- get coordinates of current int. p. on parent domain of l-s element */
      f2_mlgcoor2(smfunct,smxyzep,smiel,coor);

/*---- get values of large-scale shape functions and their derivatives */
      switch(ntyp)
      {
      case 1:    /* --> quad - element */
        f2_rec(funct,deriv,deriv2,coor[0],coor[1],typ,icode);
      break;
      case 2:	/* --> tri - element */ 	     
        f2_tri(funct,deriv,deriv2,coor[0],coor[1],typ,icode);   
      break;
      default:
        dserror("ntyp unknown!\n");      
      } /*end switch(ntyp) */
/*------------ compute Jacobian matrix for large-scale shape functions */
      f2_mljaco(funct,deriv,xjm,&det,ele,iel);
/*----------- compute global derivates for large-scale shape functions */
      f2_gder(derxy,deriv,xjm,det,iel);

/*------------ get large-scale velocities (n+1,i) at integration point */
      f2_veli(velint,funct,evel,iel);
/*-- get large-scale velocity (n+1,i) derivatives at integration point */
      f2_vder(vderxy,derxy,evel,iel);
      
/*-----------------------------------------------------------------------
     get values of bubble functions and their derivatives 
----------------------------------------------------------------------- */
/*------------------- get velocity bubble functions at integraton point */
      fluid_bubint (vbubint,smfunct,evbub,smiel,mlvar->nvbub);            
/*-------- get velocity bubble function derivatives at integraton point */
      fluid_bubder (vbubderxy,smderxy,evbub,smiel,mlvar->nvbub,2);            
/*----------------- get 'pressure' bubble functions at integraton point */
      fluid_pbubint (pbubint,smfunct,epbub,smiel,iel,2);            
/*------ get 'pressure' bubble function derivatives at integraton point */
      fluid_pbubder (pbubderxy,smderxy,epbub,smiel,iel,2,2);            
/*--------------------------- get small-scale 'rhs' at integraton point */
      fluid_bubint (smfint,smfunct,efbub,smiel,2);            
/*--------------- get small-scale 'rhs' derivatives at integraton point */
      fluid_bubder (smfderxy,smderxy,efbub,smiel,2,2);            
	 
      if (mlvar->convel==0)
      { 
/*---------------------- get small-scale velocities at integraton point */
        f2_veli (smvelint,vbubint,evel,iel);	      
/*------------ get small-scale velocity derivatives at integraton point */
        f2_vder (smvderxy,vbubderxy,evel,iel);	      
/*--------------------- get small-scale 'pressures' at integraton point */
        f2_smprei (smpreint,pbubint,epre,iel);	      
/*---------- get small-scale 'pressure' derivatives at integraton point */
        f2_smpder (smpderxy,pbubderxy,epre,iel);	      

/*--- calculate velocities and velocity derivatives at integraton point */
        for (i=0; i<2; i++)
	{
	  velint[i] += smvelint[i] + smpreint[i] + smfint[i];
	  for (j=0; j<2; j++)
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
      if (dynvar->nik>0 && mlvar->convel==0)
      {
/*------------------------ compute standard Galerkin part of matrix Kvv */      
        f2_lscalkvv(dynvar,estif,velint,vderxy,funct,derxy,fac,visc,iel);
      }	 

/*----------------------------------------------------------------------*
 |         compute "Bubble" matrix                                      |
 | NOTE:                                                                |
 |  Bubble matrices are all stored in one matrix "estif"                |
 |  Bubble mass matrix is stored in "emass"                             |
 *----------------------------------------------------------------------*/
/*-------------------------- compute bubble function part of matrix Kvv */      
      f2_calbkvv(dynvar,estif,velint,vderxy,funct,derxy,vbubint,vbubderxy,
                 fac,visc,iel);
/*-------------------------- compute bubble function part of matrix Kvp */      
     f2_calbkvp(dynvar,estif,velint,vderxy,funct,derxy,pbubint,pbubderxy,
                 fac,visc,iel);
/*-------------------------- compute bubble function part of matrix Kpv */      
      f2_calbkpv(estif,funct,vbubderxy,fac,iel);
/*-------------------------- compute bubble function part of matrix Kpp */      
      f2_calbkpp(estif,funct,pbubderxy,fac,iel);
	
      if (dynvar->nis==0 && mlvar->transterm==0 || 
          dynvar->nis==0 && mlvar->transterm==2)
      {	   
/*-------------------------- compute bubble function part of matrix Mvv */      
	f2_calbmvv(emass,funct,vbubint,fac,iel);
/*-------------------------- compute bubble function part of matrix Mvp */      
	f2_calbmvp(emass,funct,pbubint,fac,iel);
      }    	
/*----------------------------------------------------------------------*
 |         compute "VMM" Force Vectors                                  |
 *----------------------------------------------------------------------*/ 
      if (dynvar->nis==0)
      {
/*------------------------ compute "VMM" force vector for velocity dofs */      
        f2_calbfv(dynvar,mlvar,eiforce,velint,vderxy,funct,derxy,smfint,
                  smfderxy,fac,visc,iel); 
/*------------------------ compute "VMM" force vector for pressure dofs */      
        f2_calbfp(dynvar,eiforce,funct,smfderxy,fac,iel); 
      }
      
/*----------------------------------------------------------------------*
 |        compute "Iteration" Force Vectors                             |
 *----------------------------------------------------------------------*/ 
      if (dynvar->nii!=0)
      {
        if (mlvar->convel!=0)
        {
/*--------------------- get small-scale velocities at integration point */
          f2_veli (smvelint,vbubint,evel,iel);	      
/*----------- get small-scale velocity derivatives at integration point */
          f2_vder (smvderxy,vbubderxy,evel,iel);	      
/*-------------------- get small-scale 'pressures' at integration point */
          f2_smprei (smpreint,pbubint,epre,iel);	      
/*--------- get small-scale 'pressure' derivatives at integration point */
          f2_smpder (smpderxy,pbubderxy,epre,iel);	      

/*-- calculate velocities and velocity derivatives at integration point */
          for (i=0; i<2; i++)
	  {
	    velint[i] = smvelint[i] + smpreint[i] + smfint[i];
	    for (j=0; j<2; j++)
	    {
	      vderxy[i][j] = smvderxy[i][j] + smpderxy[i][j] + smfderxy[j][i];
	    }
	  }
	}
/*-------------- get convective velocities (n+1,i) at integration point */
        f2_covi(vderxy,velint,covint);
/*- calculate "Iteration" force vector for velocity dofs (no pre. dofs) */
        f2_lscalgalifv(dynvar,eiforce,covint,velint,vderxy,funct,fac,iel);
      } /* endif (dynvar->nii!=0) */
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_bubint */

/*!---------------------------------------------------------------------
\brief integration loop for one large-scale element for fluid2

<pre>                                                       gravem 07/03

In this routine the large-scale part of the large-scale element 
stiffness matrix, iteration-RHS, time-RHS and external-RHS for one 
large-scale element is calculated
      
</pre>
\param  *data      FLUID_DATA	   (i)	  integration data
\param  *ele	   ELEMENT	   (i)    actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param  *hasext    INT             (i)    element flag
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
void f2_lsint(FLUID_DATA      *data,     
	      ELEMENT	      *ele,	
	      FLUID_DYN_CALC  *dynvar,
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
INT       iel;        /* number of nodes                                */
INT       ntyp;       /* element type: 1 - quad; 2 - tri                */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis;    /* number of integration nodesin r,s direction    */
INT       actmat;     /* material number of the element                 */
INT       ihoel=0;    /* flag for higher order elements                 */
INT       icode=2;    /* flag for eveluation of shape functions         */     
INT       lr, ls;     /* counter for integration                        */
DOUBLE    dens;       /* density                                        */
DOUBLE    visc;       /* viscosity                                      */
DOUBLE    fac;
DOUBLE    facr, facs; /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DOUBLE    preintn;    /* pressure at integration point at times step n  */
DIS_TYP   typ;	      /* element type                                   */

#ifdef DEBUG 
dstrc_enter("f2_lsint");
#endif

/*----------------------------------------------------- initialization */
iel=ele->numnp;
actmat=ele->mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity;
ntyp = ele->e.f2->ntyp; 
typ  = ele->distyp;

/*------ get integration data and check if elements are "higher order" */
switch (ntyp)
{
case 1:  /* --> quad - element */
   icode   = 3;
   ihoel   = 1;
   /* initialize integration */
   nir = ele->e.f2->nGP[0];
   nis = ele->e.f2->nGP[1];
break;
case 2: /* --> tri - element */  
   if (iel>3)
   {
      icode   = 3;
      ihoel   = 1;
   }
   /* initialize integration */
   nir  = ele->e.f2->nGP[0];
   nis  = 1;
   intc = ele->e.f2->nGP[1];  
break;
default:
   dserror("ntyp unknown!");
} /* end switch(ntyp) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{    
   for (ls=0;ls<nis;ls++)
   {
/*----------------- get values of shape functions and their derivatives */
      switch(ntyp)  
      {
      case 1:   /* --> quad - element */
	 e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
         f2_rec(funct,deriv,deriv2,e1,e2,typ,icode);
      break;
      case 2:   /* --> tri - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
	 f2_tri(funct,deriv,deriv2,e1,e2,typ,icode);
      break;
      default:
         dserror("ntyp unknown!");
      } /* end switch(ntyp) */
/*-------------------------------------------- compute Jacobian matrix */
      f2_mljaco(funct,deriv,xjm,&det,ele,iel);
      fac = facr*facs*det;
/*------------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);

/*------------------------ get velocities (n+1,i) at integration point */
      f2_veli(velint,funct,evel,iel);
/*-------------- get velocity derivatives (n+1,i) at integration point */
      f2_vder(vderxy,derxy,evel,iel);
      
 /*--- compute stab. par. or subgrid viscosity during integration loop */
      if (ele->e.f2->istabi>0 && ele->e.f2->iduring!=0 || 
          dynvar->sgvisc>0    && ele->e.f2->iduring!=0)
        f2_mlcalelesize2(ele,dynvar,funct,velint,vderxy,wa1,visc,iel,ntyp);
	
/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrices                         |
 | NOTE:                                                                |
 |  Standard Galerkin matrices are all stored in one matrix "estif"     |
 |  Standard Galerkin mass matrix is stored in "emass"                  |
 *----------------------------------------------------------------------*/
      if(dynvar->nik>0)
      {
/*------------------------ compute standard Galerkin part of matrix Kvv */      
         if (mlvar->convel!=0)
	   f2_lscalkvv(dynvar,estif,velint,vderxy,funct,derxy,fac,visc,iel);
/*-------------- compute standard Galerkin part of matrices Kvp and Kpv */      
	 f2_lscalkvp(estif,funct,derxy,fac,iel);
/*------------------------ compute standard Galerkin part of matrix Mvv */      
	 if (dynvar->nis==0) f2_lscalmvv(emass,funct,fac,iel);
      } /* endif (dynvar->nik>0) */
      
/*----------------------------------------------------------------------*
 |         compute Stabilization matrices                               |
 |  (in most of the cases only the bulk viscosity term is calculated)   |
 | NOTE:                                                                |
 |  Stabilization matrices are all stored in one matrix "estif"         |
 |  Stabilization mass matrices are all stored in one matrix "emass"    |
 *----------------------------------------------------------------------*/
      if (ele->e.f2->istabi>0)
      { 
/*----------------------------------- compute second global derivatives */ 
         if (ihoel!=0) f2_mlgder2(ele,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
   
         if (dynvar->nie==0)
         {
/*---------------------------------------- stabilization for matrix Kvv */
            f2_lscalstabkvv(ele,dynvar,estif,velint,vderxy,
                          funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilization for matrix Kvp */
            f2_lscalstabkvp(ele,dynvar,estif,velint,vderxy,
                          funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilization for matrix Mvv */
            if (dynvar->nis==0) 
               f2_lscalstabmvv(ele,dynvar,emass,velint,vderxy, 
	                     funct,derxy,derxy2,fac,visc,iel,ihoel); 
            if (ele->e.f2->ipres!=0)	        
            {
/*---------------------------------------- stabilization for matrix Kpv */
               f2_lscalstabkpv(dynvar,estif,velint,vderxy,
	                     funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilization for matrix Mpv */
	       if (dynvar->nis==0)
		  f2_lscalstabmpv(dynvar,emass,funct,derxy,fac,iel);
            } /* endif (ele->e.f2->ipres!=0) */
         } /* endif (dynvar->nie==0) */
/*---------------------------------------- stabilization for matrix Kpp */
         if (ele->e.f2->ipres!=0)
	    f2_lscalstabkpp(dynvar,estif,derxy,fac,iel);  
      } /* endif (ele->e.f2->istabi>0) */

/*----------------------------------------------------------------------*
 |         compute "Iteration" Force Vectors                            |
 *----------------------------------------------------------------------*/ 
      if (dynvar->nii!=0 && mlvar->convel!=0)
      {
/*-------------- get convective velocities (n+1,i) at integration point */
        f2_covi(vderxy,velint,covint);
/*- calculate "Iteration" force vector for velocity dofs (no pre. dofs) */
        f2_lscalgalifv(dynvar,eiforce,covint,velint,vderxy,funct,fac,iel);
      } /* endif (dynvar->nii!=0) */

/*----------------------------------------------------------------------*
 |         compute "Time" Force Vector                                  |
 *----------------------------------------------------------------------*/
      if (dynvar->nis==0)
      {
 /*------------------------------ get pressure (n) at integration point */
        if (dynvar->iprerhs>0) f2_prei(&preintn,funct,epren,iel);
/*----------------------------- get velocities (n) at integration point */
        f2_veli(velintn,funct,eveln,iel);
/*------------------- get velocity derivatives (n) at integration point */
        f2_vder(vderxyn,derxy,eveln,iel);
/*------------------ get convective velocities (n) at integration point */
        f2_covi(vderxyn,velintn,covintn);        	    
/*--------------------- calculate "Time" force vector for velocity dofs */
	f2_lscalgaltfv(dynvar,etforce,velintn,velintn,covintn,funct,derxy,
	             vderxyn,preintn,visc,fac,iel);		       
/*--------------------- calculate "Time" force vector for pressure dofs */
        if (dynvar->thsr!=ZERO) 
	  f2_lscalgaltfp(dynvar,&(etforce[2*iel]),funct,vderxyn,fac,iel);
      }  

/*----------------------------------------------------------------------*
 |         compute "External" Force Vector                              |
 *----------------------------------------------------------------------*/
      if (*hasext!=0)
/*----------------- calculate "External" force vector for velocity dofs */
        f2_lscalgalexfv(dynvar,etforce,funct,edeadn,edead,fac,iel);
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_lsint */

/*!---------------------------------------------------------------------
\brief integration loop for dynamic procedure for sm ele. for fluid2

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
void f2_smint2(FLUID_DATA      *data,     
	       ELEMENT	       *ele,	
	       FLUID_ML_SMESH  *submesh, 
               DOUBLE	      **smiediff,   
	       DOUBLE	       *smierhs,   
	       DOUBLE	      **smxyze, 
	       DOUBLE	       *smfunct,   
	       DOUBLE	      **smderiv,   
	       DOUBLE	      **smderiv2,  
	       DOUBLE	      **smxjm,	  
	       DOUBLE	      **smderxy)
{ 
INT       smiel;      /* submesh number of nodes                        */
INT       nsmtyp;     /* submesh element type: 1 - quad; 2 - tri        */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis;    /* number of integration nodes in r,s direction   */
INT       lr, ls;     /* counter for integration                        */
DOUBLE    fac;
DOUBLE    facr, facs; /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DIS_TYP   typ,smtyp;  /* large-scale and submesh element type           */

#ifdef DEBUG 
dstrc_enter("f2_smint2");
#endif

/*------------------------------------------------------ initialization */
smiel  = submesh->numen;
nsmtyp = submesh->ntyp; 
smtyp  = submesh->typ;

/*--------------------------- get integration data for submesh elements */
switch (nsmtyp)
{
case 1:  /* --> quad - element */
   nir = submesh->ngpr;
   nis = submesh->ngps;
break;
case 2: /* --> tri - element */  
   nir  = submesh->ngpr;
   nis  = 1;
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
/*--------- get values of submesh shape functions and their derivatives */
      switch(nsmtyp)  
      {
      case 1:   /* --> quad - element */
	 e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
         f2_rec(smfunct,smderiv,smderiv2,e1,e2,smtyp,2);
      break;
      case 2:   /* --> tri - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
	 f2_tri(smfunct,smderiv,smderiv2,e1,e2,smtyp,2);
      break;
      default:
         dserror("nsmtyp unknown!");
      } /* end switch(nsmtyp) */
/*--------------------------------- compute Jacobian matrix for submesh */
      f2_mljaco3(smxyze,smfunct,smderiv,smxjm,&det,smiel,ele);
      fac = facr*facs*det;
/*-------------------------------- compute global derivates for submesh */
      f2_gder(smderxy,smderiv,smxjm,det,smiel);
      
/*-------------- compute diffusive part of submesh stiffness matrix SMK */      
      f2_calsmkd(smiediff,smderxy,fac,smiel);
/*---------------------------- compute integral of normalized RHS force */
      fluid_calnofo(smierhs,smfunct,fac,smiel); 
  }
}
        
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_smint2 */

/*!---------------------------------------------------------------------
\brief integration loop for sub-submesh element for fluid2

<pre>                                                       gravem 07/03

In this routine, the element stiffness matrix and normalized RHS for 
one sub-submesh element is calculated.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *dynvar      FLUID_DYN_CALC(i)
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
\return void                                                   

------------------------------------------------------------------------*/
void f2_ssint(FLUID_DATA      *data,     
	      ELEMENT	      *ele,	
	      FLUID_DYN_CALC  *dynvar, 
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
	      DOUBLE	     **vderxy)
{ 
INT       iel,ssiel;  /* large-scale and sub-submesh number of nodes    */
INT       ntyp,nsstyp;/* l-s and ssm element type: 1 - quad; 2 - tri    */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis;    /* number of integration nodes in r,s direction   */
INT       actmat;     /* material number of the element                 */
INT       lr, ls;     /* counter for integration                        */
DOUBLE    visc;       /* viscosity                                      */
DOUBLE    fac;
DOUBLE    facr, facs; /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DOUBLE    coor[2];    /* coordinates                                    */
DIS_TYP   typ,sstyp;  /* large-scale and sub-submesh element type       */

#ifdef DEBUG 
dstrc_enter("f2_ssint");
#endif

/*------------------------------------------------------ initialization */
iel    = ele->numnp;
ssiel  = ssmesh->numen;
actmat = ele->mat-1;
visc   = mat[actmat].m.fluid->viscosity;
ntyp   = ele->e.f2->ntyp; 
typ    = ele->distyp;
nsstyp = ssmesh->ntyp; 
sstyp  = ssmesh->typ;

/*---------------------- get integration data for sub-submesh elements */
switch (nsstyp)
{
case 1:  /* --> quad - element */
   nir = ssmesh->ngpr;
   nis = ssmesh->ngps;
break;
case 2: /* --> tri - element */  
   nir  = ssmesh->ngpr;
   nis  = 1;
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
/*---- get values of sub-submesh shape functions and their derivatives */
      switch(nsstyp)  
      {
      case 1:   /* --> quad - element */
	 e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
         f2_rec(ssfunct,ssderiv,ssderiv2,e1,e2,sstyp,2);
      break;
      case 2:   /* --> tri - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
	 f2_tri(ssfunct,ssderiv,ssderiv2,e1,e2,sstyp,2);
      break;
      default:
         dserror("nsstyp unknown!");
      } /* end switch(nsstyp) */
/*---------------------------- compute Jacobian matrix for sub-submesh */
      f2_mljaco3(ssxyze,ssfunct,ssderiv,ssxjm,&det,ssiel,ele);
      fac = facr*facs*det;
/*--------------------------- compute global derivates for sub-submesh */
      f2_gder(ssderxy,ssderiv,ssxjm,det,ssiel);

/*- get coordinates of current int. p. on parent domain of l-s element */
      f2_mlgcoor2(ssfunct,ssxyzep,ssiel,coor);
      
/*---- get values of large-scale shape functions and their derivatives */
      switch(ntyp)
      {
      case 1:    /* --> quad - element */
        f2_rec(funct,deriv,deriv2,coor[0],coor[1],typ,2);
      break;
      case 2:	/* --> tri - element */ 	     
        f2_tri(funct,deriv,deriv2,coor[0],coor[1],typ,2);   
      break;
      default:
        dserror("ntyp unknown!\n");      
      } /*end switch(ntyp) */
/*------------ compute Jacobian matrix for large-scale shape functions */
      f2_mljaco(funct,deriv,xjm,&det,ele,iel);
/*----------- compute global derivates for large-scale shape functions */
      f2_gder(derxy,deriv,xjm,det,iel);

/*------------ get large-scale velocities (n+1,i) at integration point */
      f2_veli(velint,funct,evel,iel);
/*-- get large-scale velocity (n+1,i) derivatives at integration point */
      f2_vder(vderxy,derxy,evel,iel);
      
/*----------------------------------------------------------------------*
 |       compute matrix and rhs for sub-submesh                         |
 *----------------------------------------------------------------------*/
/*------------------------------------------------  compute matrix SSK */      
      f2_calsmk(dynvar,mlvar,ssestif,velint,vderxy,ssfunct,ssderxy,fac,
                visc,ssiel);
/*--------------------------------------------- compute normalized rhs */
      fluid_calnofo(ssenfor,ssfunct,fac,ssiel); 
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_ssint */

/*!---------------------------------------------------------------------
\brief integration loop for norm. bubble on sub-submesh ele. for fluid2

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
void f2_inbu(FLUID_DATA      *data,     
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
INT       nsstyp;     /* ssm element type: 1 - quad; 2 - tri            */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis;    /* number of integration nodesin r,s direction    */
INT       lr, ls;     /* counter for integration                        */
DOUBLE    fac;
DOUBLE    facr, facs; /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DOUBLE    bubnint;    /* normalized bubble function at int. point       */
DIS_TYP   sstyp;      /* sub-submesh element type                       */

#ifdef DEBUG 
dstrc_enter("f2_inbu");
#endif

/*------------------------------------------------------ initialization */
ssiel  = ssmesh->numen;
nsstyp = ssmesh->ntyp; 
sstyp  = ssmesh->typ;

/*----------------------- get integration data for sub-submesh elements */
switch (nsstyp)
{
case 1:  /* --> quad - element */
   nir = ssmesh->ngpr;
   nis = ssmesh->ngps;
break;
case 2: /* --> tri - element */  
   nir  = ssmesh->ngpr;
   nis  = 1;
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
/*---- get values of sub-submesh shape functions and their derivatives */
      switch(nsstyp)  
      {
      case 1:   /* --> quad - element */
	 e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
         f2_rec(ssfunct,ssderiv,ssderiv2,e1,e2,sstyp,2);
      break;
      case 2:   /* --> tri - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
	 f2_tri(ssfunct,ssderiv,ssderiv2,e1,e2,sstyp,2);
      break;
      default:
         dserror("nsstyp unknown!");
      } /* end switch(nsstyp) */
/*---------------------------- compute Jacobian matrix for sub-submesh */
      f2_mljaco3(ssxyze,ssfunct,ssderiv,ssxjm,&det,ssiel,ele);
      fac = facr*facs*det;

/*------------ compute normalized bubble function at integration point */
      bubnint=ZERO; 
      for (i=0;i<ssiel;i++) /* loop all nodes of the ssm element */
      {
        bubnint += ssfunct[i]*ebub[i];
      } /* end loop all nodes of the ssm element */
      
/*----- compute sub-submesh integral of the normalized bubble function */
     *ssinbu += bubnint*fac; 
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_inbu */

#endif
