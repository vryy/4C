/*!----------------------------------------------------------------------
\file
\brief time RHS for fluid3 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"
/*!--------------------------------------------------------------------- 
\brief galerkin part of time forces for vel dofs

<pre>                                                         genk 05/02

In this routine the galerkin part of the time forces for vel dofs
is calculated:

      /
   + |  v * u     d_omega
    /  
    
                      /
   (-) (1-THETA)*dt  |  v * u * grad(u)     d_omega
                    /
		  
                      /
   (-) (1-THETA)*dt  |  2*nue * eps(v) : eps(u)  d_omega
                    /  
		  
                    /
   + (1-THETA)*dt  |  div(v) * p  d_omega
                  /		  		      

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:					     
   in ONESTEP methods: velint  = vel2int = U(n)
   in TWOSTEP methods: velint  = U(n+gamma)  
   in TWOSTEP methods: vel2int = U(n)	     
      
</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param   *vel2int     DOUBLE	      (i)    vel. at integr. point
\param   *covint      DOUBLE	      (i)    conv. vel. at integr. p.
\param	 *funct       DOUBLE	      (i)    nat. shape functions      
\param	**derxy       DOUBLE	      (i)    global derivatives
\param	**vderxy      DOUBLE	      (i/    global vel. deriv.
\param	  preint      DOUBLE	      (i)    pres. at integr. point
\param	  visc	      DOUBLE	      (i)    fluid viscosity
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calgaltfv(
                  FLUID_DYN_CALC  *dynvar, 
                  DOUBLE          *eforce,    
		  DOUBLE          *velint,
		  DOUBLE          *vel2int,
		  DOUBLE          *covint,
		  DOUBLE          *funct,
		  DOUBLE         **derxy,
		  DOUBLE         **vderxy,
		  DOUBLE           preint,
		  DOUBLE           visc,
		  DOUBLE           fac,
		  INT              iel
              )  
{
INT    j,irow,isd,inode;  
DOUBLE c;
DOUBLE aux;
DOUBLE facsr;
DOUBLE facpr;
DOUBLE fact[3];
 
 #ifdef DEBUG 
dstrc_enter("f3_calgaltfv");
#endif

/*--------------------------------------------------- set some factors */
facsr = fac * dynvar->thsr;
facpr = fac * dynvar->thpr;
c     = facsr * visc;

/*----------------------------------------------------------------------*
   Calculate intertia forces of time force vector:
      /
   + |  v * u     d_omega
    /  
 *----------------------------------------------------------------------*/ 

fact[0] = vel2int[0]*fac;
fact[1] = vel2int[1]*fac;
fact[2] = vel2int[2]*fac;
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<3;isd++)
   {
      irow++;
      eforce[irow] += funct[inode]*fact[isd];
   } /* end of loop over isd */
} /* end of loop over inode */

/*----------------------------------------------------------------------*
   Calculate convective forces of time force vector:
                    /
   - (1-THETA)*dt  |  v * u * grad(u)     d_omega
                  / 
 *----------------------------------------------------------------------*/
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<3;isd++)
   {
      irow++;
      eforce[irow] -= funct[inode]*covint[isd]*facsr;
   } /* end of loop over isd */
} /* end of loop over inode */

/*----------------------------------------------------------------------*
   Calculate viscous forces of time force vector:
                    /
   - (1-THETA)*dt  |  2*nue * eps(v) : eps(u)  d_omega
                  /  
 *----------------------------------------------------------------------*/ 
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<3;isd++)
   {
      irow++;
      for (j=0;j<3;j++)
      {
         eforce[irow] -= (derxy[j][inode]*(vderxy[isd][j]+vderxy[j][isd])*c);  	                 
      } /* end of loop over j */
   } /* end of loop over isd */
} /* end of loop over inode */

/*----------------------------------------------------------------------*
   Calculate pressure forces of time force vector:
                    /
   + (1-THETA)*dt  |  div(v) * p  d_omega
                  /
 *----------------------------------------------------------------------*/ 
if (dynvar->iprerhs>0)
{
   aux = preint * facpr;
   irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      for (isd=0;isd<3;isd++)
      {
         irow++;
         eforce[irow] += derxy[isd][inode]*aux;
      } /* end of loop over isd */
   } /* end of loop over inode */
} /*endif (dynvar->iprerhs>0) */
/*----------------------------------------------------------------------*
   Calculate external forces of time force vector (due to self-weight):
                    /
   + (1-THETA)*dt  |  v * b  d_omega
                  / 
 *----------------------------------------------------------------------*/
/* NOT IMPLEMENTED YET!!!!! */

/*----------------------------------------------------------------------*
   Calculate external forces of time force vector
                         (due to surface tension):
                    /
   + (1-THETA)*dt  |  v * h  d_gamma
                  / 
 *----------------------------------------------------------------------*/
/* NOT IMPLEMENTED YET!!!!!!!!!!!!!!!
   maybe its better to do the integration over the element edge in 
   a seperate loop (different gauss-points!!!!)
*/
 
/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calgaltfv */

/*!--------------------------------------------------------------------- 
\brief galerkin part of time forces for pre dofs

<pre>                                                         genk 05/02

In this routine the galerkin part of the time forces for pre dofs
is calculated:

                    /
   + (1-THETA)*dt  |  q * div(u)  d_omega
                  /

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:						       
    there's only one full element force vector         
    for pre-dofs the pointer eforce points to the entry
    eforce[3*iel]				        	    
      
</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param   *funct       DOUBLE	      (i)    nat. shape functions
\param  **vderxy      DOUBLE	      (i)    global vel. deriv.
\param    fac         DOUBLE	      (i)    weighting factor	     
\param	  iel	      INT	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calgaltfp(
                   FLUID_DYN_CALC  *dynvar, 
                   DOUBLE          *eforce,    
		   DOUBLE          *funct,
		   DOUBLE         **vderxy,
		   DOUBLE           fac,
		   INT              iel
                  ) 
{
INT      inode;
DOUBLE   aux;
DOUBLE   facsr;

#ifdef DEBUG 
dstrc_enter("f3_calgaltfp");
#endif

/*--------------------------------------------------- set some factors */
facsr = fac * dynvar->thsr;

/*----------------------------------------------------------------------*
   Calculate continuity forces of time force vector:
                    /
   + (1-THETA)*dt  |  q * div(u)  d_omega
                  /
 *----------------------------------------------------------------------*/
aux = facsr * (vderxy[0][0] + vderxy[1][1] + vderxy[2][2]);
for(inode=0;inode<iel;inode++)
{
   eforce[inode] += funct[inode]*aux;
} /* end of loop over inode */

/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calgaltfp */

/*!--------------------------------------------------------------------- 
\brief stabilisation part of time forces for vel dofs

<pre>                                                         genk 05/02

In this routine the stabilisation part of the time forces for vel dofs
is calculated:
		 		                   	  		      
      /
   + |  tau_mu * u * grad(v) * u  d_omega
    /
    
        /
   -/+ |  tau_mp * 2*nue * div( eps(v) ) * u  d_omega
      /
      
                     /
   (-) (1-THETA)*dt |  tau_mu * u * grad(v) * u * grad(u) d_omega
                   /          
    
                     /
   +/- (1-THETA)*dt |  tau_mp * 2*nue * div( eps(v) ) * u * grad(u)  d_omega
                   /

                   /
   + (1-THETA)*dt |  tau_mu * 2*nue * u *grad(v) * div( eps(u) )  d_omega
                 /

                     /
   -/+ (1-THETA)*dt |  tau_mp * 4*nue^2 * div( eps(v) ) * div ( eps(u) )  d_omega
                   /		 
		 	
       /
  (-) |  tau_c * div(v) * div(u)  d_omega
     /
   
                    /
  (-) (1-THETA)*dt | tau_mu * u * grad(v) * grad(p)   d_omega
                  /
		
                    /
  -/+ (1-THETA)*dt | tau_mp * 2*nue * div( eps(v) ) * grad(p)    d_omega
                  /
		  
                  /
  + (1-THETA)*dt | tau_mu * u * grad(v) * b    d_omega
                /
				       
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
      
</pre>
\param   *dynvar     FLUID_DYN_CALC   (i)
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
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabtfv(
                   FLUID_DYN_CALC  *dynvar, 
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
		   INT              iel
                  ) 
{
INT    irow,isd,inode;
DOUBLE c,cc;
DOUBLE aux;
DOUBLE taumu,taump,tauc; 
DOUBLE facsr;
DOUBLE facpr;
DOUBLE fvts,fvtsr,fvvtsr;
DOUBLE fact[3];
DOUBLE sign;
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f3_calstabtfv");
#endif	

/*---------------------------------------------------------- initialise */
gls    = ele->e.f3->stabi.gls;

if (ele->e.f3->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");
  
/*--------------------------------------------------- set some factors */
taumu = dynvar->tau[0];
taump = dynvar->tau[1];
tauc  = dynvar->tau[2];

facsr = fac * dynvar->thsr;
facpr = fac * dynvar->thpr;
c     = facsr * visc;

/*----------------------------------------------------------------------*
   Calculate inertia/convective stab-forces of time force vector:
      /
   + |  tau_mu * u * grad(v) * u  d_omega
    /
 *----------------------------------------------------------------------*/
if (gls->iadvec!=0)
{
   fact[0] = vel2int[0]*fac*taumu;
   fact[1] = vel2int[1]*fac*taumu;
   fact[2] = vel2int[2]*fac*taumu;
   irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1] \
          + derxy[2][inode]*velint[2];
      for (isd=0;isd<3;isd++)
      {
         irow++;
	 eforce[irow] += aux*fact[isd];
      } /* end of loop over isd */
   } /* end of loop over inode */
} /* endif (ele->e.f3->iadvec!=0) */

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
   } /* end switch (ele->e.f3->ivisc) */

   fvts = fac*visc*taump*sign;
   irow = 0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy2[0][inode] + derxy2[1][inode] + derxy2[2][inode];
      eforce[irow]   -= ( derxy2[0][inode]*vel2int[0]  \
                        + derxy2[3][inode]*vel2int[1]  \
		        + derxy2[4][inode]*vel2int[2] + aux*vel2int[0])*fvts;
      eforce[irow+1] -= ( derxy2[3][inode]*vel2int[0]  \
                        + derxy2[1][inode]*vel2int[1]  \
		        + derxy2[5][inode]*vel2int[2] + aux*vel2int[1])*fvts;
      eforce[irow+2] -= ( derxy2[4][inode]*vel2int[0]  \
                        + derxy2[5][inode]*vel2int[1]  \
		        + derxy2[2][inode]*vel2int[2] + aux*vel2int[2])*fvts;
      irow += 3;			 
   } /* end of loop over inode */   
} /* endif (ele->e.f3->ivisc!=0 && ihoel!=0) */
 
/*----------------------------------------------------------------------*
   Calculate convective/convective stab-forces of time force vector:
                   /
   - (1-THETA)*dt |  tau_mu * u * grad(v) * u * grad(u) d_omega
                 /
 *----------------------------------------------------------------------*/ 
if (gls->iadvec!=0)
{
   fact[0] = taumu*covint[0]*facsr;
   fact[1] = taumu*covint[1]*facsr;
   fact[2] = taumu*covint[2]*facsr;
   irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1] \
          + derxy[2][inode]*velint[2];
      for (isd=0;isd<3;isd++)
      {
         irow++;
	 eforce[irow] -= aux*fact[isd];
      } /* end of loop over isd */
   } /* end of loop over inode */
} /* endif (ele->e.f3->iadvec!=0) */

/*----------------------------------------------------------------------*
   Calculate convective/viscous stab-forces of time force vector:
                     /
   +/- (1-THETA)*dt |  tau_mp * 2*nue * div( eps(v) ) * u * grad(u)  d_omega
                   /
 *----------------------------------------------------------------------*/
if (gls->ivisc!=0 && ihoel!=0)
{
   fvtsr = fvts * dynvar->thsr;
   irow = 0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy2[0][inode] + derxy2[1][inode] + derxy2[2][inode];
      eforce[irow]   += (derxy2[0][inode]*covint[0] \
                       + derxy2[3][inode]*covint[1] \
		       + derxy2[4][inode]*covint[2] + aux*covint[0])*fvtsr;
      eforce[irow+1] += (derxy2[3][inode]*covint[0] \
                       + derxy2[1][inode]*covint[1] \
		       + derxy2[5][inode]*covint[2] + aux*covint[1])*fvtsr;
      eforce[irow+2] += (derxy2[4][inode]*covint[0] \
                       + derxy2[5][inode]*covint[1] \
		       + derxy2[2][inode]*covint[2] + aux*covint[2])*fvtsr;
      irow +=3;
   } /* end of loop over inode */
} /* endif (ele->e.f3->ivisc!=0 && ihoel!=0) */
 
/*----------------------------------------------------------------------*
   Calculate viscous/convective stab-forces of time force vector:
                   /
   + (1-THETA)*dt |  tau_mu * 2*nue * u *grad(v) * div( eps(u) )  d_omega
                 /
 *----------------------------------------------------------------------*/ 
if (gls->iadvec!=0 && ihoel!=0)
{
   cc = c*taumu;
   fact[0] = (TWO*vderxy2[0][0] + vderxy2[0][1] + vderxy2[0][2] \
                + vderxy2[1][3] + vderxy2[2][4])*cc;
   fact[1] = (TWO*vderxy2[1][1] + vderxy2[1][0] + vderxy2[1][2] \
                + vderxy2[0][3] + vderxy2[2][5])*cc;
   fact[2] = (TWO*vderxy2[2][2] + vderxy2[2][0] + vderxy2[2][1] \
                + vderxy2[0][4] + vderxy2[1][5])*cc;
   irow=0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1]  \
          + derxy[2][inode]*velint[2];
      eforce[irow]   += aux*fact[0];
      eforce[irow+1] += aux*fact[1];
      eforce[irow+2] += aux*fact[2];
      irow += 3;
   } /* end of loop over inode */
} /* endif (ele->e.f3->iadvec!=0 && ihoel!=0) */

/*----------------------------------------------------------------------*
   Calculate viscous/viscous stab-forces of time force vector:
                     /
   -/+ (1-THETA)*dt |  tau_mp * 4*nue^2 * div( eps(v) ) * div ( eps(u) )  d_omega
                   /
 *----------------------------------------------------------------------*/ 
if (gls->ivisc!=0 && ihoel!=0)
{
   fvvtsr = fvtsr*visc;
   fact[0] = (TWO*vderxy2[0][0] + vderxy2[1][3] + vderxy2[2][4] \
                + vderxy2[0][1] + vderxy2[0][2])*fvvtsr;
   fact[1] = (TWO*vderxy2[1][1] + vderxy2[0][3] + vderxy2[2][5] \
                + vderxy2[1][0] + vderxy2[1][2])*fvvtsr;
   fact[2] = (TWO*vderxy2[2][2] + vderxy2[0][4] + vderxy2[1][5] \
                + vderxy2[2][0] + vderxy2[2][1])*fvvtsr;
   irow=0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy2[0][inode] + derxy2[1][inode] + derxy2[2][inode];
      eforce[irow]   -= (derxy2[0][inode]*fact[0] \
                       + derxy2[3][inode]*fact[1] \
		       + derxy2[4][inode]*fact[2] + aux*fact[0]);
      eforce[irow+1] -= (derxy2[3][inode]*fact[0] \
                       + derxy2[1][inode]*fact[1] \
		       + derxy2[5][inode]*fact[2] + aux*fact[1]);       
      eforce[irow+2] -= (derxy2[4][inode]*fact[0] \
                       + derxy2[5][inode]*fact[1] \
		       + derxy2[2][inode]*fact[2] + aux*fact[2]);
   irow += 3;
   } /* end of loop over inode */
} /* endif (ele->e.f3->ivisc!=0 && ihoel!=0) */
   
/*----------------------------------------------------------------------*
   Calculate continuity stab-forces of time force vector:
     /
  - |  tau_c * div(v) * div(u)  d_omega
   /
 *----------------------------------------------------------------------*/  
if (gls->icont!=0)
{
   aux = tauc*facsr*(vderxy[0][0] + vderxy[1][1] + vderxy[2][2]);
   irow = -1;
   for (inode=0;inode<iel;inode++)
   {
      for (isd=0;isd<3;isd++)
      {
         irow++;
	 eforce[irow] -= derxy[isd][inode]*aux;
      } /* end of loop over isd */
   } /* end of loop over inode */
} /* endif (ele->e.f3->ivisc!=0) */
 
/*----------------------------------------------------------------------*
   Calculate pressure/convective stab-forces of time force vector:
                  /
  - (1-THETA)*dt | tau_mu * u * grad(v) * grad(p)   d_omega
                /
 *----------------------------------------------------------------------*/  
if (gls->iadvec!=0 && dynvar->iprerhs>0)
{
   fact[0] = taumu*pderxy[0]*facpr;
   fact[1] = taumu*pderxy[1]*facpr;
   fact[2] = taumu*pderxy[2]*facpr;
   irow = -1;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1] \
          + derxy[2][inode]*velint[2];
      for (isd=0;isd<3;isd++)
      {
         irow++;
	 eforce[irow] -= aux*fact[isd];	 
      } /* end of loop over isd */
   } /* end of loop over inode */
} /* endif (ele->e.f3->iadvec!=0 && dynvar->iprerhs>0) */

/*----------------------------------------------------------------------*
   Calculate pressure/viscous stab-forces of time force vector:
                    /
  -/+ (1-THETA)*dt | tau_mp * 2*nue * div( eps(v) ) * grad(p)    d_omega
                  /
 *----------------------------------------------------------------------*/  
if (gls->ivisc!=0 && ihoel!=0 && dynvar->iprerhs>0)
{
   cc = facpr*visc*taump*sign;
   irow=0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy2[0][inode] + derxy2[1][inode] + derxy2[2][inode];
      eforce[irow]   += (derxy2[0][inode]*pderxy[0]  
                       + derxy2[3][inode]*pderxy[1] 
		       + derxy2[4][inode]*pderxy[2] + aux*pderxy[0])*cc;
      eforce[irow+1] += (derxy2[3][inode]*pderxy[0]  
                       + derxy2[1][inode]*pderxy[1] 
		       + derxy2[5][inode]*pderxy[2] + aux*pderxy[1])*cc;
      eforce[irow+2] += (derxy2[4][inode]*pderxy[0] 
                       + derxy2[5][inode]*pderxy[1] 
		       + derxy2[2][inode]*pderxy[2] + aux*pderxy[2])*cc;
      irow += 3;		
   } /* end of loop over inode */
} /* endif (ele->e.f3->ivisc!=0 && ihoel!=0 && dynvar->iprerhs>0) */
 
/*----------------------------------------------------------------------*
   Calculate external load/convective stab-forces of time force vector:
                  /
  + (1-THETA)*dt | tau_mu * u * grad(v) * b    d_omega
                /
 *----------------------------------------------------------------------*/
/* NOT IMPLEMENTED YET!!!!*/

/*----------------------------------------------------------------------*
   Calculate external load/viscous stab-forces of time force vector:
                    /
  -/+ (1-THETA)*dt | tau_mp * 2*nue * div( eps(v) ) * b    d_omega
                  /
 *----------------------------------------------------------------------*/
/* NOT IMPLEMENTED YET!!!!*/ 

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabtfv */

/*!--------------------------------------------------------------------- 
\brief stabilisation part of time forces for pre dofs

<pre>                                                         genk 05/02

In this routine the stabilisation part of the time forces for vel dofs
is calculated:
		 		                   	  		      
          /
   (-)   |  tau_mp * grad(q) * u  d_omega
        /
      
                     /
   + (1-THETA)*dt   |  tau_mp * grad(q) * u * grad(u)  d_omega
                   /
		   
                     /
   + (1-THETA)*dt   |  tau_mp * grad(q) * grad(p)  d_omega
                   /		         
				       
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:						       
    there's only one full element force vector         
    for pre-dofs the pointer eforce points to the entry
    eforce[3*iel]				       
      
</pre>
\param   *dynvar,     FLUID_DYN_CALC   (i)
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
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabtfp(
                   FLUID_DYN_CALC  *dynvar, 
                   DOUBLE          *eforce,    
     		   DOUBLE         **derxy,
		   DOUBLE         **vderxy2,
		   DOUBLE          *velint,
		   DOUBLE          *covint,
		   DOUBLE          *pderxy,
		   DOUBLE           visc,
		   DOUBLE           fac,
		   INT              ihoel,
		   INT              iel
                  ) 
{
INT      inode;
DOUBLE   fact[3];
DOUBLE   facsr,facpr;
DOUBLE   taump;

#ifdef DEBUG 
dstrc_enter("f3_calstabtfp");
#endif

/*--------------------------------------------------- set some factors */
facsr = fac * dynvar->thsr;
facpr = fac * dynvar->thpr;
taump = dynvar->tau[1];

/*----------------------------------------------------------------------*
   Calculate inertia/pressure stab forces of time force vector:
        /
   -   |  tau_mp * grad(q) * u  d_omega
      /
 *----------------------------------------------------------------------*/
fact[0] = velint[0]*taump*fac;
fact[1] = velint[1]*taump*fac;
fact[2] = velint[2]*taump*fac;
for (inode=0;inode<iel;inode++)
{
   eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1] \
                   + derxy[2][inode]*fact[2]);
} /* end of loop over inode */

/*----------------------------------------------------------------------*
   Calculate convective/pressure stab forces of time force vector:
                     /
   + (1-THETA)*dt   |  tau_mp * grad(q) * u * grad(u)  d_omega
                   /
 *----------------------------------------------------------------------*/
fact[0] = covint[0]*taump*facsr;
fact[1] = covint[1]*taump*facsr;
fact[2] = covint[2]*taump*facsr;
for (inode=0;inode<iel;inode++)
{
   eforce[inode] += (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1] 
                   + derxy[2][inode]*fact[2]);
} /* end of loop over inode */


/*----------------------------------------------------------------------*
   Calculate vsicous/pressure stab forces of time force vector:
                     /
   - (1-THETA)*dt   |  tau_mp * 2*nue * grad(q) * div( eps(u) )  d_omega
                   /
 *----------------------------------------------------------------------*/
if (ihoel!=0)
{
   fact[0] = (TWO*vderxy2[0][0] + vderxy2[1][3] + vderxy2[2][4] \
                + vderxy2[0][1] + vderxy2[0][2])*taump*visc*facsr;   
   fact[1] = (TWO*vderxy2[1][1] + vderxy2[0][3] + vderxy2[2][5] \
                + vderxy2[1][0] + vderxy2[1][2])*taump*visc*facsr;
   fact[2] = (TWO*vderxy2[2][2] + vderxy2[0][4] + vderxy2[1][5] \
                + vderxy2[2][0] + vderxy2[2][1])*taump*visc*facsr;
   for (inode=0;inode<iel;inode++)
   {
      eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1] \
                      + derxy[2][inode]*fact[2]);
   } /* end of loop over inode */
} /* endif (ihoel!=0) */

/*----------------------------------------------------------------------*
   Calculate pressure/pressure stab forces of time force vector:
                     /
   + (1-THETA)*dt   |  tau_mp * grad(q) * grad(p)  d_omega
                   /
 *----------------------------------------------------------------------*/
if (dynvar->iprerhs!=0)
{
   fact[0] = taump*pderxy[0]*facpr;
   fact[1] = taump*pderxy[1]*facpr;
   fact[2] = taump*pderxy[2]*facpr;
   for (inode=0;inode<iel;inode++)
   {
      eforce[inode] += (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1] \
                      + derxy[2][inode]*fact[2]);
   } /* end of loop over inode */
} /* endif (dynvar->iprerhs!=0) */

/*----------------------------------------------------------------------*
   Calculate external load/pressure stab forces of time force vector:
                     /
   - (1-THETA)*dt   |  tau_mp * grad(q) * b  d_omega
                   /
 *----------------------------------------------------------------------*/
/* NOT IMPLEMENTED YET!!!! */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabtfp */



#endif
 
