/*!----------------------------------------------------------------------
\file
\brief time RHS for fluid2_xfem element

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_XFEM
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid2/fluid2.h"
#include "xfem_prototypes.h"
/*! 
\addtogroup XFEM 
*//*! @{ (documentation module open)*/



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



/*!--------------------------------------------------------------------- 
\brief galerkin part of time forces for vel dofs

<pre>                                                            irhan 05/04

In this routine the galerkin part of the time forces for vel dofs
is calculated:

EULER/ALE:
      /
   + |  v * u     d_omega
    /  
    
EULER:
                      /
   (-) (1-THETA)*dt  |  v * u * grad(u)     d_omega
                    /

EULER/ALE:
                      /
   (-) (1-THETA)*dt  |  2*nue * eps(v) : eps(u)  d_omega
                    /  
		  
                    /
   + (1-THETA)*dt  |  div(v) * p  d_omega
                  /		  		      

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE: vel2int = U(n)	     
      EULER: covint = U(n) * grad(U(n))
      ALE:   covint = C(n) * grad(U(n))
   
</pre>
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param   *vel2int     DOUBLE	      (i)    vel. at integr. point
\param   *covint      DOUBLE	      (i)    conv. vel. at integr. p.
\param	 *funct       DOUBLE	      (i)    nat. shape functions      
\param	**derxy       DOUBLE	      (i)    global derivatives
\param	**vderxy      DOUBLE	      (i/    global vel. deriv.
\param	  preint      DOUBLE	      (i)    pres. at integr. point
\param	  visc	      DOUBLE	      (i)    fluid viscosity
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\param    index	      INT  	      (i)    index for local assembly
\param    DENS	      DOUBLE  	      (i)    fluid density
\return void                                                                       

------------------------------------------------------------------------*/
void xfem_f2_calgaltfv(
  DOUBLE          *eforce,    
  DOUBLE          *vel2int,    
  DOUBLE          *covint,   
  DOUBLE          *funct,    
  DOUBLE         **derxy,    
  DOUBLE         **vderxy,   
  DOUBLE           preint,   
  DOUBLE           visc,     
  DOUBLE           fac,      
  INT              iel,
  INT             *index,
  DOUBLE           DENS
  )  
{
  INT        j,irow,isd,inode;  
  DOUBLE     c;
  DOUBLE     aux;
  DOUBLE     facsr;
  DOUBLE     facpr;
  DOUBLE     fact[2];
 
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calgaltfv");
#endif
/*---------------------------------------------------------------------*/

  /* initialize */
  fdyn = alldyn[genprob.numff].fdyn;
  
  /* set some factors */
  facsr = fac * fdyn->thsr;
  facpr = fac * fdyn->thpr;
  c     = facsr * visc;
/*----------------------------------------------------------------------*
   Calculate intertia forces of time force vector:
      /
   + |  v * u     d_omega
    /  
 *----------------------------------------------------------------------*/
  fact[0] = vel2int[0]*fac;
  fact[1] = vel2int[1]*fac;
  for (inode=0; inode<TWO*iel; inode++)
  {
    irow = index[inode];
    for (isd=0; isd<2; isd++)
    {
      eforce[irow] += funct[inode]*fact[isd]*DENS;
      irow++;
    }
  }
/*----------------------------------------------------------------------*
   Calculate convective forces of time force vector:
EULER:
                    /
   - (1-THETA)*dt  |  v * u * grad(u)     d_omega
                  / 
 *----------------------------------------------------------------------*/
  for (inode=0; inode<TWO*iel; inode++)
  {
    irow = index[inode];
    for (isd=0; isd<2; isd++)
    {
      eforce[irow] -= funct[inode]*covint[isd]*facsr*DENS;
      irow++;
    }
  }
/*----------------------------------------------------------------------*
   Calculate viscous forces of time force vector:
                    /
   - (1-THETA)*dt  |  2*nue * eps(v) : eps(u)  d_omega
                  /  
 *----------------------------------------------------------------------*/
  for (inode=0; inode<TWO*iel; inode++)
  {
    irow = index[inode];
    for (isd=0; isd<2; isd++)
    {
      for (j=0;j<2;j++)
     {
       eforce[irow] -= (derxy[j][inode]*(vderxy[isd][j]+vderxy[j][isd])*c);
     } 
      irow++;
    } 
  }
/*----------------------------------------------------------------------*
   Calculate pressure forces of time force vector:
                    /
   + (1-THETA)*dt  |  div(v) * p  d_omega
                  /
 *----------------------------------------------------------------------*/ 
  if (fdyn->iprerhs>0)
  {
   aux = preint * facpr;
   for (inode=0; inode<TWO*iel; inode++)
   {
     irow = index[inode];
     for (isd=0; isd<2; isd++)
     {
       eforce[irow] += derxy[isd][inode]*aux;
       irow++;
     }
   }
  }
  
/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calgaltfv */



/*!--------------------------------------------------------------------- 
\brief stabilisation part of time forces for vel dofs

<pre>                                                            irhan 05/04

In this routine the stabilisation part of the time forces for vel dofs
is calculated:

EULER:		 		                   	  		      
      /
   + |  tau_mu * u * grad(v) * u  d_omega
    /

EULER/ALE:
        /
   -/+ |  tau_mp * 2*nue * div( eps(v) ) * u  d_omega
      /
      
EULER:
                     /
   (-) (1-THETA)*dt |  tau_mu * u * grad(v) * u * grad(u) d_omega
                   /
		   
EULER:
                     /
   +/- (1-THETA)*dt |  tau_mp * 2*nue * div( eps(v) ) * u * grad(u)  d_omega
                   /
   
EULER:		   
                   /
   + (1-THETA)*dt |  tau_mu * 2*nue * u *grad(v) * div( eps(u) )  d_omega
                 /

EULER/ALE:		 
                     /
   -/+ (1-THETA)*dt |  tau_mp * 4*nue^2 * div( eps(v) ) * div ( eps(u) )  d_omega
                   /		 
		 	
       /
  (-) |  tau_c * div(v) * div(u)  d_omega
     /
   
EULER:
                    /
  (-) (1-THETA)*dt | tau_mu * u * grad(v) * grad(p)   d_omega
                  /

EULER/ALE:
                    /
  -/+ (1-THETA)*dt | tau_mp * 2*nue * div( eps(v) ) * grad(p)    d_omega
                  /		  
				       
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE: for EULER
      in ONESTEP methods: velint = vel2int = U(n);  
      in TWOSTEP methods: velint = U(n+gamma)      
      in TWOSTEP methods: vel2int= U(n)
      in ONESTEP methods: covint = U(n) * grad(U(n))
      in TWOSTEP methods: covint = U(n+gamma) * grad(U(n+gamma))
       
NOTE: for ALE
      only ONESTEP method implemented!!!
        velint = C(n)
	vel2int= U(n)
	covint = C(n) * grad(U(n))
</pre>
\param   *ele        ELEMENT	      (i)    actual element
\param   *eforce     DOUBLE	      (i/o)  element force vector
\param   *velint     DOUBLE	      (i)    vel. at integr. point
\param   *vel2int    DOUBLE	      (i)    vel. at integr. point
\param	 *covint     DOUBLE	      (i)    conv. vel. at integr. p.
\param	**derxy      DOUBLE	      (i)    global derivatives
\param	**derxy2     DOUBLE	      (i)    2nd global derivatives
\param	**vderxy     DOUBLE	      (i)    global vel. deriv.
\param	**vderxy2    DOUBLE	      (i)    2nd global vel. deriv.
\param	 *pderxy     DOUBLE	      (i)    global pressure deriv.
\param	  fac	     DOUBLE	      (i)    weighting factor
\param	  visc	     DOUBLE	      (i)    fluid viscosity
\param	  ihoel      INT	      (i)    flag for higer ord. ele
\param	  iel	     INT	      (i)    num. of nodes in ele
\param    index	     INT  	      (i)    index for local assembly
\param    DENS	     DOUBLE  	      (i)    fluid density
\return void                                                                       

------------------------------------------------------------------------*/
void xfem_f2_calstabtfv(
  ELEMENT         *ele,      
  DOUBLE          *eforce,  
  DOUBLE          *velint,  
  DOUBLE          *vel2int, 
  DOUBLE          *covint,  
  DOUBLE         **derxy,   
  DOUBLE         **derxy2,  
  DOUBLE         **vderxy,  
  DOUBLE         **vderxy2, 
  DOUBLE          *pderxy,  
  DOUBLE           fac,     
  DOUBLE           visc,    
  INT              ihoel,   
  INT              iel,
  INT             *index,
  DOUBLE           DENS
  )
{
  INT        irow,isd,inode;
  DOUBLE     c,cc;
  DOUBLE     aux;
  DOUBLE     taumu,taump,tauc;
  DOUBLE     facsr;
  DOUBLE     facpr;
  DOUBLE     fvts,fvtsr,fvvtsr;
  DOUBLE     fact[2];
  DOUBLE     sign;

  STAB_PAR_GLS  *gls;	/* pointer to GLS stabilisation parameters	*/  
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabtfv");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */
  fdyn = alldyn[genprob.numff].fdyn;
  gls = ele->e.f2->stabi.gls;
  
  /* set some factors */
  taumu = fdyn->tau[0];
  taump = fdyn->tau[1];
  tauc  = fdyn->tau[2];
  
  facsr = fac * fdyn->thsr;
  facpr = fac * fdyn->thpr;
  c     = facsr * visc;
/*----------------------------------------------------------------------*
   Calculate inertia/convective stab-forces of time force vector:
EULER:		 		                   	  		      
      /
   + |  tau_mu * u * grad(v) * u  d_omega
    /
 *----------------------------------------------------------------------*/
  if (gls->iadvec!=0)
  {
    fact[0] = vel2int[0]*fac*taumu;
    fact[1] = vel2int[1]*fac*taumu;
    for (inode=0; inode<TWO*iel; inode++)
    {
      irow = index[inode];
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1];
      for (isd=0; isd<2; isd++)
      {
        eforce[irow] += aux*fact[isd]*DENS*DENS;
        irow++;
      }
    }
  }
/*----------------------------------------------------------------------*
   Calculate inertia/viscous stab-forces of time force vector:
        /
   -/+ |  tau_mp * 2*nue * div( eps(v) ) * u  d_omega
      /
 *----------------------------------------------------------------------*/
  if (gls->ivisc!=0 && ihoel!=0)
  {
    switch (gls->ivisc) /* choose stabilisation type --> sign */
    {
        case 1: /* GLS- */
          sign = ONE;
          break;
        case 2: /* GLS+ */
          sign = -ONE;
          break;
        default:
          dserror("viscous stabilisation parameter unknown: IVISC");
    }
    
    fvts = fac*visc*taump*sign;
    for (inode=0; inode<TWO*iel; inode++)
    {
      irow = index[inode];
      eforce[irow  ] -= ((TWO*derxy2[0][inode] + derxy2[1][inode])*vel2int[0] +
                              derxy2[2][inode]*vel2int[1])*fvts*DENS;
      eforce[irow+1] -= ((TWO*derxy2[1][inode] + derxy2[0][inode])*vel2int[1] +
                              derxy2[2][inode]*vel2int[0])*fvts*DENS;
    }
  }
/*----------------------------------------------------------------------*
   Calculate convective/convective stab-forces of time force vector:
EULER
                   /
   - (1-THETA)*dt |  tau_mu * u * grad(v) * u * grad(u) d_omega
                 /
 *----------------------------------------------------------------------*/ 
  if (gls->iadvec!=0)
  {
    fact[0] = taumu*covint[0]*facsr;
    fact[1] = taumu*covint[1]*facsr;
    for (inode=0; inode<TWO*iel; inode++)
    {
      irow = index[inode];      
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1];
      for (isd=0; isd<2; isd++)
      {
        eforce[irow] -= aux*fact[isd]*DENS*DENS;
        irow++;
      }
    }
  }
/*----------------------------------------------------------------------*
   Calculate convective/viscous stab-forces of time force vector:
EULER:
                     /
   +/- (1-THETA)*dt |  tau_mp * 2*nue * div( eps(v) ) * u * grad(u)  d_omega
                   /

ALE:
                     /
   +/- (1-THETA)*dt |  tau_mp * 2*nue * div( eps(v) ) * c * grad(u)  d_omega
                   /
 *----------------------------------------------------------------------*/
  if (gls->ivisc!=0 && ihoel!=0)
  {
    fvtsr = fvts * fdyn->thsr;
    for (inode=0; inode<TWO*iel; inode++)
    {
      irow = index[inode];            
      aux = derxy2[0][inode] + derxy2[1][inode];
      eforce[irow  ] += ((derxy2[0][inode] + aux)*covint[0] +
                          derxy2[2][inode]*covint[1])*fvtsr*DENS;
      eforce[irow+1] += ((derxy2[1][inode] + aux)*covint[1] +
                          derxy2[2][inode]*covint[0])*fvtsr*DENS;
    }
  }
/*----------------------------------------------------------------------*
   Calculate viscous/convective stab-forces of time force vector:
EULER:
                   /
   + (1-THETA)*dt |  tau_mu * 2*nue * u *grad(v) * div( eps(u) )  d_omega
                 /
 *----------------------------------------------------------------------*/ 
  if (gls->iadvec!=0 && ihoel!=0)
  {
    cc = c*taumu;
    fact[0] = (TWO*vderxy2[0][0] + vderxy2[0][1] + vderxy2[1][2])*cc;
    fact[1] = (TWO*vderxy2[1][1] + vderxy2[1][0] + vderxy2[0][2])*cc;
    for (inode=0; inode<TWO*iel; inode++)
    {
      irow = index[inode];                  
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1];
      eforce[irow  ] += aux*fact[0]*DENS;
      eforce[irow+1] += aux*fact[1]*DENS;
    }
  }
/*----------------------------------------------------------------------*
   Calculate viscous/viscous stab-forces of time force vector:
                     /
   -/+ (1-THETA)*dt |  tau_mp * 4*nue^2 * div( eps(v) ) * div ( eps(u) )  d_omega
                   /
 *----------------------------------------------------------------------*/ 
  if (gls->ivisc!=0 && ihoel!=0)
  {
    fvvtsr = fvtsr*visc;
    fact[0] = (TWO*vderxy2[0][0] + vderxy2[0][1] + vderxy2[1][2])*fvvtsr;
    fact[1] = (TWO*vderxy2[1][1] + vderxy2[1][0] + vderxy2[0][2])*fvvtsr;
    for (inode=0; inode<TWO*iel; inode++)
    {
      irow = index[inode];                        
      aux = derxy2[0][inode] + derxy2[1][inode];
      eforce[irow  ] -= (derxy2[0][inode]*fact[0] +
                         derxy2[2][inode]*fact[1] + aux*fact[0]);
      eforce[irow+1] -= (derxy2[2][inode]*fact[0] +
                         derxy2[1][inode]*fact[1] + aux*fact[1]);       
    }
  }
/*----------------------------------------------------------------------*
   Calculate continuity stab-forces of time force vector:
     /
  - |  tau_c * div(v) * div(u)  d_omega
   /
 *----------------------------------------------------------------------*/  
  if (gls->icont!=0) 
  {
    aux = tauc*facsr*(vderxy[0][0] + vderxy[1][1]);
    for (inode=0; inode<TWO*iel; inode++)
    {
      irow = index[inode];
      for (isd=0;isd<2;isd++)
      {
        eforce[irow] -= derxy[isd][inode]*aux;
        irow++;
      }
    }
  }
/*----------------------------------------------------------------------*
   Calculate pressure/convective stab-forces of time force vector:
EULER:
                    /
  (-) (1-THETA)*dt | tau_mu * u * grad(v) * grad(p)   d_omega
                  /
 *----------------------------------------------------------------------*/  
  if (gls->iadvec!=0 && fdyn->iprerhs>0)
  {
    fact[0] = taumu*pderxy[0]*facpr;
    fact[1] = taumu*pderxy[1]*facpr;
    for (inode=0; inode<TWO*iel; inode++)
    {
      irow = index[inode];      
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1];
      for (isd=0;isd<2;isd++)
      {
        eforce[irow] -= aux*fact[isd]*DENS;
        irow++;
      }
    }
  }
/*----------------------------------------------------------------------*
   Calculate pressure/viscous stab-forces of time force vector:
                    /
  -/+ (1-THETA)*dt | tau_mp * 2*nue * div( eps(v) ) * grad(p)    d_omega
                  /
 *----------------------------------------------------------------------*/  
  if (gls->ivisc!=0 && ihoel!=0 && fdyn->iprerhs>0)
  {
    cc = facpr*visc*taump*sign;
    for (inode=0; inode<TWO*iel; inode++)
    {
      irow = index[inode];
      aux = derxy2[0][inode] + derxy2[1][inode];
      eforce[irow  ] += ((derxy2[0][inode] + aux)*pderxy[0] +
                          derxy2[2][inode]*pderxy[1])*cc;
      eforce[irow+1] += ((derxy2[1][inode] + aux)*pderxy[1] +
                          derxy2[2][inode]*pderxy[0])*cc;
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calstabtfv */



/*!--------------------------------------------------------------------- 
\brief stabilisation part of time forces for pre dofs

<pre>                                                         genk 04/02
                                             modified for ALE genk 10/02

In this routine the stabilisation part of the time forces for vel dofs
is calculated:

EULER/ALE:		 		                   	  		      
          /
   (-)   |  tau_mp * grad(q) * u  d_omega
        /
      
EULER:
                     /
   + (1-THETA)*dt   |  tau_mp * grad(q) * u * grad(u)  d_omega
                   /
		   		   
EULER/ALE:
                     /
   + (1-THETA)*dt   |  tau_mp * grad(q) * grad(p)  d_omega
                   /		         
				       
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:						       
    there's only one full element force vector         
    for pre-dofs the pointer eforce points to the entry
    eforce[2*iel]				       

NOTE: for EULER
      velint = U(n)
      covint = U(n) * grad(U(n))

NOTE: for ALE:
      velint = C(n) (ale-convective velocity)
      covint = C(n) * grad(U(n))      
      
</pre>
\param   *eforce,     DOUBLE	       (i/o)  element force vector
\param  **derxy,      DOUBLE	       (i)    global derivatives
\param  **vderxy2,    DOUBLE	       (i)    2nd global vel. deriv.
\param   *velint,     DOUBLE	       (i)    vel. at integr. point
\param	 *covint,     DOUBLE	       (i)    conv. vel. at integr. p.
\param	 *pderxy,     DOUBLE	       (i)    global pressure deriv.
\param	  visc,       DOUBLE	       (i)    fluid viscosity
\param	  fac,        DOUBLE	       (i)    weighting factor
\param	  ihoel,      INT	       (i)    flag for higer ord. ele
\param	  iel	      INT	       (i)    num. of nodes in ele
\param    DENS	      DOUBLE  	       (i)    fluid density
\return void                                                                       

------------------------------------------------------------------------*/
void xfem_f2_calstabtfp(
  DOUBLE          *eforce,    
  DOUBLE         **derxy,   
  DOUBLE         **vderxy2, 
  DOUBLE          *velint,  
  DOUBLE          *covint,  
  DOUBLE          *pderxy,  
  DOUBLE           visc,    
  DOUBLE           fac,     
  INT              ihoel,   
  INT              iel,
  DOUBLE           DENS
  ) 
{
INT      inode;
DOUBLE   fact[2];
DOUBLE   facsr,facpr;
DOUBLE   taump;

#ifdef DEBUG 
dstrc_enter("xfem_f2_calstabtfp");
#endif

/*--------------------------------------------------- set some factors */
fdyn  = alldyn[genprob.numff].fdyn;

facsr = fac * fdyn->thsr;
facpr = fac * fdyn->thpr;
taump = fdyn->tau[1];

/*----------------------------------------------------------------------*
   Calculate inertia/pressure stab forces of time force vector:
        /
   -   |  tau_mp * grad(q) * u  d_omega
      /
 *----------------------------------------------------------------------*/
fact[0] = velint[0]*taump*fac;
fact[1] = velint[1]*taump*fac;
for (inode=0;inode<iel;inode++)
{
   eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1])*DENS;
}
 
/*----------------------------------------------------------------------*
   Calculate convective/pressure stab forces of time force vector:
EULER:
                     /
   + (1-THETA)*dt   |  tau_mp * grad(q) * u * grad(u)  d_omega
                   /
 *----------------------------------------------------------------------*/
fact[0] = covint[0]*taump*facsr;
fact[1] = covint[1]*taump*facsr;
for (inode=0;inode<iel;inode++)
{
   eforce[inode] += (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1])*DENS;
}
/*----------------------------------------------------------------------*
   Calculate vsicous/pressure stab forces of time force vector:
                     /
   - (1-THETA)*dt   |  tau_mp * 2*nue * grad(q) * div( eps(u) )  d_omega
                   /
 *----------------------------------------------------------------------*/
if (ihoel!=0)
{
   fact[0] = (TWO*vderxy2[0][0] + vderxy2[0][1] + vderxy2[1][2]) \
             *taump*visc*facsr;   
   fact[1] = (TWO*vderxy2[1][1] + vderxy2[1][0] + vderxy2[0][2]) \
             *taump*visc*facsr;
   for (inode=0;inode<iel;inode++)
   {
      eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1]);
   }
} 
 
/*----------------------------------------------------------------------*
   Calculate pressure/pressure stab forces of time force vector:
                     /
   + (1-THETA)*dt   |  tau_mp * grad(q) * grad(p)  d_omega
                   /
 *----------------------------------------------------------------------*/
if (fdyn->iprerhs!=0)
{
   fact[0] = taump*pderxy[0]*facpr;
   fact[1] = taump*pderxy[1]*facpr;
   for (inode=0;inode<iel;inode++)
   {
      eforce[inode] += (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1]);
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of xfem_f2_calstabtfp */
/*! @} (documentation module close)*/
#endif
