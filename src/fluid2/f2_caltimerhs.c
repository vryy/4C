#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#define ONE (1.0)
#define TWO (2.0)
/*----------------------------------------------------------------------*
 | routine to calculate galerkin part of time forces for vel dofs       |
 | NOTE:							        |
 |    in ONESTEP methods: velint  = vel2int = U(n)		        |
 |    in TWOSTEP methods: velint  = U(n+gamma)			        |
 |    in TWOSTEP methods: vel2int = U(n)			        |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calgaltfv(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,   /* element force vector */
		  double          *velint,
		  double          *vel2int,
		  double          *covint,
		  double          *funct,
		  double         **derxy,
		  double         **vderxy,
		  double           preint,
		  double           visc,
		  double           fac,
		  int              iel
              )  
{
int    j,irow,isd,inode;  
double c;
double aux;
double facsr;
double facpr;
double fact[2];
 
 
#ifdef DEBUG 
dstrc_enter("f2_calgaltfv");
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
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      eforce[irow] += funct[inode]*fact[isd];
   }
}

/*----------------------------------------------------------------------*
   Calculate convective forces of time force vector:
                    /
   - (1-THETA)*dt  |  v * u * grad(u)     d_omega
                  / 
 *----------------------------------------------------------------------*/
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      eforce[irow] -= funct[inode]*covint[isd]*facsr;
   }
}

/*----------------------------------------------------------------------*
   Calculate viscous forces of time force vector:
                    /
   - (1-THETA)*dt  |  2*nue * eps(v) : eps(u)  d_omega
                  /  
 *----------------------------------------------------------------------*/ 
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      for (j=0;j<2;j++)
      {
         eforce[irow] -= (derxy[j][inode]*(vderxy[isd][j]+vderxy[j][isd])*c);  	                 
      }
   }
}

/*----------------------------------------------------------------------*
   Calculate pressure forces of time force vector:
                    /
   + (1-THETA)*dt  |  div(v) * p  d_omega
                  /
 *----------------------------------------------------------------------*/ 
aux = preint * facpr;
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      eforce[irow] += derxy[isd][inode]*aux;
   }
}

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
 
 
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calgaltfv */

/*----------------------------------------------------------------------*
 | routine to calculate galerkin part of time forces for pre dofs       |
 | NOTE:                                                                |
 |     there's only one full element force vector                       |  
 |     for pre-dofs the pointer eforce points to the entry              |
 |     eforce[2*iel]                                                    |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calgaltfp(
                   FLUID_DYN_CALC  *dynvar, 
                   double          *eforce,   /* element force vector */
		   double          *funct,
		   double         **vderxy,
		   double           fac,
		   int              iel
                  ) 
{
int      inode;
double   aux;
double   facsr;

#ifdef DEBUG 
dstrc_enter("f2_calgaltfp");
#endif

/*--------------------------------------------------- set some factors */
facsr = fac * dynvar->thsr;

/*----------------------------------------------------------------------*
   Calculate continuity forces of time force vector:
                    /
   + (1-THETA)*dt  |  q * div(u)  d_omega
                  /
 *----------------------------------------------------------------------*/
aux = facsr * (vderxy[0][0] + vderxy[1][1]);
for(inode=0;inode<iel;inode++)
{
   eforce[inode] += funct[inode]*aux;
}

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calgaltfp */

/*----------------------------------------------------------------------*
 | routine to calculate stabilisation part of time forces for vel dofs  |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calstabtfv(
                   FLUID_DYN_CALC  *dynvar, 
                   ELEMENT         *ele,
	           double          *eforce,   /* element force vector */
	 	   double          *velint,
		   double          *vel2int,
		   double          *covint,
		   double         **derxy,
		   double         **derxy2,
		   double         **vderxy,
		   double         **vderxy2,
		   double          *pderxy,
		   double           fac,
		   double           visc,
		   int              ihoel,
		   int              iel
                  ) 
{
int    j,irow,isd,inode;
double c,cc;
double aux;
double taumu,taump,tauc;
double facsr;
double facpr;
double fvts,fvtsr,fvvtsr;
double fact[2];
double sign;
 
 
#ifdef DEBUG 
dstrc_enter("f2_calstabtfv");
#endif
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
if (ele->e.f2->iadvec!=0)
{
   fact[0] = vel2int[0]*fac*taumu;
   fact[1] = vel2int[1]*fac*taumu;
   irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1];
      for (isd=0;isd<2;isd++)
      {
         irow++;
	 eforce[irow] += aux*fact[isd];
      }
   }
}

/*----------------------------------------------------------------------*
   Calculate inertia/viscous stab-forces of time force vector:
        /
   -/+ |  tau_mp * 2*nue * div( eps(v) ) * u  d_omega
      /
 *----------------------------------------------------------------------*/
if (ele->e.f2->ivisc!=0 && ihoel!=0)
{
   switch (ele->e.f2->ivisc) /* choose stabilisation type --> sign */
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
   irow = 0;
   for (inode=0;inode<iel;inode++)
   {
      eforce[irow]   -= ((TWO*derxy2[0][inode] + derxy2[1][inode])*vel2int[0] \
                            + derxy2[2][inode]*vel2int[1])*fvts;
      eforce[irow+1] -= ((TWO*derxy2[1][inode] + derxy2[0][inode])*vel2int[1] \
                            + derxy2[2][inode]*vel2int[0])*fvts;
      irow += 2;			 
   }   
} 

/*----------------------------------------------------------------------*
   Calculate convective/convective stab-forces of time force vector:
                   /
   - (1-THETA)*dt |  tau_mu * u * grad(v) * u * grad(u) d_omega
                 /
 *----------------------------------------------------------------------*/ 
if (ele->e.f2->iadvec!=0)
{
   fact[0] = taumu*covint[0]*facsr;
   fact[1] = taumu*covint[1]*facsr;
   irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1];
      for (isd=0;isd<2;isd++)
      {
         irow++;
	 eforce[irow] -= aux*fact[isd];
      }
   }
} 

/*----------------------------------------------------------------------*
   Calculate convective/viscous stab-forces of time force vector:
                     /
   +/- (1-THETA)*dt |  tau_mp * 2*nue * div( eps(v) ) * u * grad(u)  d_omega
                   /
 *----------------------------------------------------------------------*/
if (ele->e.f2->ivisc!=0 && ihoel!=0)
{
   fvtsr = fvts * dynvar->thsr;
   irow = 0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy2[0][inode] + derxy[1][inode];
      eforce[irow]   += ((derxy2[0][inode] + aux)*covint[0] 
                        + derxy2[2][inode]*covint[1])*fvtsr;
      eforce[irow+1] += ((derxy2[1][inode] + aux)*covint[1] 
                        + derxy2[2][inode]*covint[0])*fvtsr;
      irow +=2;
   }   
}

/*----------------------------------------------------------------------*
   Calculate viscous/convective stab-forces of time force vector:
                   /
   + (1-THETA)*dt |  tau_mu * 2*nue * u *grad(v) * div( eps(u) )  d_omega
                 /
 *----------------------------------------------------------------------*/ 
if (ele->e.f2->iadvec!=0 && ihoel!=0)
{
   cc = c*taumu;
   fact[0] = (TWO*vderxy2[0][0] + vderxy2[0][1] + vderxy2[1][2])*cc;
   fact[1] = (TWO*vderxy2[1][1] + vderxy2[1][0] + vderxy2[0][2])*cc;
   irow=0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1];
      eforce[irow]   += aux*fact[0];
      eforce[irow+1] += aux*fact[1];
   }
   irow += 2;
}
 
  
/*----------------------------------------------------------------------*
   Calculate viscous/viscous stab-forces of time force vector:
                     /
   -/+ (1-THETA)*dt |  tau_mp * 4*nue^2 * div( eps(v) ) * div ( eps(u) )  d_omega
                   /
 *----------------------------------------------------------------------*/ 
if (ele->e.f2->ivisc!=0 && ihoel!=0)
{
   fvvtsr = fvtsr*visc;
   fact[0] = (TWO*vderxy2[0][0] + vderxy2[0][1] + vderxy2[1][2])*fvvtsr;
   fact[1] = (TWO*vderxy2[1][1] + vderxy2[1][0] + vderxy2[0][2])*fvvtsr;
   irow=0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy2[0][inode] + derxy2[1][inode];
      eforce[irow]   -= (derxy2[0][inode]*fact[0] \
                       + derxy2[2][inode]*fact[1] + aux*fact[0]);
      eforce[irow+1] -= (derxy2[2][inode]*fact[0] \
                       + derxy2[1][inode]*fact[1] + aux*fact[1]);       
   irow += 2;
   }
} 
  
/*----------------------------------------------------------------------*
   Calculate continuity stab-forces of time force vector:
     /
  - |  tau_c * div(v) * div(u)  d_omega
   /
 *----------------------------------------------------------------------*/  
if (ele->e.f2->ivisc!=0)
{
   aux = tauc*facsr*(vderxy[0][0] + vderxy[1][1]);
   irow = -1;
   for (inode=0;inode<iel;inode++)
   {
      for (isd=0;isd<2;isd++)
      {
         irow++;
	 eforce[irow] -= derxy[isd][inode]*aux;
      }
   }
}

/*----------------------------------------------------------------------*
   Calculate pressure/convective stab-forces of time force vector:
                  /
  - (1-THETA)*dt | tau_mu * u * grad(v) * grad(p)   d_omega
                /
 *----------------------------------------------------------------------*/  
if (ele->e.f2->iadvec!=0)
{
   fact[0] = taumu*pderxy[0]*facpr;
   fact[1] = taumu*pderxy[1]*facpr;
   irow = -1;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1];
      for (isd=0;isd<2;isd++)
      {
         irow++;
	 eforce[irow] -= aux*fact[isd];	 
      }
   }
} 
 
	 
/*----------------------------------------------------------------------*
   Calculate pressure/viscous stab-forces of time force vector:
                    /
  -/+ (1-THETA)*dt | tau_mp * 2*nue * div( eps(v) ) * grad(p)    d_omega
                  /
 *----------------------------------------------------------------------*/  
if (ele->e.f2->ivisc!=0 && ihoel!=0)
{
   cc = facpr*visc*taump*sign;
   irow=0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy2[0][inode] + derxy2[1][inode];
      eforce[irow]   += ((derxy2[0][inode] + aux)*pderxy[0] \
                        + derxy2[2][inode]*pderxy[1])*cc;
      eforce[irow+1] += ((derxy2[1][inode] + aux)*pderxy[1] \
                        + derxy2[2][inode]*pderxy[0])*cc;
      irow += 2;		
   }
}

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

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calstabtfv */


/*----------------------------------------------------------------------*
 | routine to calculate stabilisation part of time forces for pre dofs  |
 | NOTE:                                                                |
 |     there's only one full element force vector                       |  
 |     for pre-dofs the pointer eforce points to the entry              |
 |     eforce[2*iel]                                                    |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calstabtfp(
                   FLUID_DYN_CALC  *dynvar, 
                   double          *eforce,   /* element force vector */ 
     		   double         **derxy,
		   double         **vderxy2,
		   double          *velint,
		   double          *covint,
		   double          *pderxy,
		   double           visc,
		   double           fac,
		   int              ihoel,
		   int              iel
                  ) 
{
int      inode;
double   aux;
double   fact[2];
double   facsr,facpr;
double   taump;

#ifdef DEBUG 
dstrc_enter("f2_calstabtfp");
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
for (inode=0;inode<iel;inode++)
{
   eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1]);
}

/*----------------------------------------------------------------------*
   Calculate convective/pressure stab forces of time force vector:
                     /
   + (1-THETA)*dt   |  tau_mp * grad(q) * u * grad(u)  d_omega
                   /
 *----------------------------------------------------------------------*/
fact[0] = covint[0]*taump*facsr;
fact[1] = covint[1]*taump*facsr;
for (inode=0;inode<iel;inode++)
{
   eforce[inode] += (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1]);
}

/*----------------------------------------------------------------------*
   Calculate vsicous/pressure stab forces of time force vector:
                     /
   - (1-THETA)*dt   |  tau_mp * 2*nue * grad(q) * div( eps(u) )  d_omega
                   /
 *----------------------------------------------------------------------*/
if (ihoel!=0)
{
   fact[0] = (TWO*vderxy2[0][0] + vderxy2[0][1] + vderxy2[1][2] \
             *taump*visc*facsr);   
   fact[0] = (TWO*vderxy2[1][1] + vderxy2[1][0] + vderxy2[0][2] \
             *taump*visc*facsr);
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
if (dynvar->iprerhs!=0)
{
   fact[0] = taump*pderxy[0]*facpr;
   fact[1] = taump*pderxy[1]*facpr;
   for (inode=0;inode<iel;inode++)
   {
      eforce[inode] += (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1]);
   }
}

/*----------------------------------------------------------------------*
   Calculate external load/pressure stab forces of time force vector:
                     /
   - (1-THETA)*dt   |  tau_mp * grad(q) * b  d_omega
                   /
 *----------------------------------------------------------------------*/
/* NOT IMPLEMENTED YET!!!! */


#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calstabtfp */



#endif
