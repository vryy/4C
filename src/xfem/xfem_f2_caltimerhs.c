#ifdef D_XFEM
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid2/fluid2.h"
#include "xfem_prototypes.h"



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





/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
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
  /* calculate intertia forces of time force vector */
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
  /* calculate convective forces of time force vector */
  for (inode=0; inode<TWO*iel; inode++)
  {
    irow = index[inode];
    for (isd=0; isd<2; isd++)
    {
      eforce[irow] -= funct[inode]*covint[isd]*facsr*DENS;
      irow++;
    }
  }
  /* calculate viscous forces of time force vector */
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
  /* calculate pressure forces of time force vector */
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



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
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
  /* calculate inertia/convective stab-forces of time force vector */
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
  /* calculate inertia/viscous stab-forces of time force vector */
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
  /* calculate convective/convective stab-forces of time force vector */
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
  /* calculate convective/viscous stab-forces of time force vector */
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
  /* calculate viscous/convective stab-forces of time force vector */
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
  /* calculate viscous/viscous stab-forces of time force vector */
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
  /* calculate continuity stab-forces of time force vector */
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
  /* calculate pressure/convective stab-forces of time force vector */
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
  /* calculate pressure/viscous stab-forces of time force vector */
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
dstrc_enter("f2_calstabtfp");
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
ALE:
                     /
   + (1-THETA)*dt   |  tau_mp * grad(q) * c * grad(u)  d_omega
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
#endif
