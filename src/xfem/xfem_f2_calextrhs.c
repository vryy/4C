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
extern ALLDYNA            *alldyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

static FLUID_DYNAMIC      *fdyn;





/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calgalexfv(
  DOUBLE             *eforce,     
  DOUBLE             *funct,       
  DOUBLE             *edeadn,
  DOUBLE             *edeadng,
  DOUBLE              fac,      
  INT                 iel,
  INT                *index,
  DOUBLE              DENS
  ) 
{
  DOUBLE  facsl,facsr;
  INT     inode,irow,isd;

#ifdef DEBUG 
  dstrc_enter("xfem_f2_calgalexfv");
#endif  
/*----------------------------------------------------------------------*/

  /* initialize */
  fdyn = alldyn[genprob.numff].fdyn;
  
  /* set some factors */
  facsl = fac*fdyn->thsl;
  facsr = fac*fdyn->thsr;

  /* calculate galerkin part of external forces => */
  for (inode=0; inode<TWO*iel; inode++)
  {
    irow = index[inode];
    for (isd=0;isd<2;isd++)
    {
      eforce[irow] += funct[inode]*(edeadn[isd]*facsr+edeadng[isd]*facsl)*DENS;
      irow++;      
    }
  }
 
/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calgalexfv */ 



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calstabexfv(
  ELEMENT            *ele,  
  DOUBLE             *eforce,     
  DOUBLE           **derxy,
  DOUBLE           **derxy2,      
  DOUBLE            *edead,
  DOUBLE            *velint,  
  DOUBLE             fac,      
  DOUBLE             visc,
  INT                iel,
  INT                ihoel,
  INT                flag,
  INT               *index,
  DOUBLE             DENS
  ) 
{
  INT         irow,inode,isd;
  DOUBLE      sign,aux;
  DOUBLE      taumu,taump;
  DOUBLE      c,fvts; 
  DOUBLE      fact[2];

  STAB_PAR_GLS  *gls;	/* pointer to GLS stabilisation parameters	*/  
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabexfv");
#endif    
/*---------------------------------------------------------------------*/

  /* initialize */
  fdyn = alldyn[genprob.numff].fdyn;
  gls = ele->e.f2->stabi.gls;

  /* set some factors */
  taumu = fdyn->tau[0];
  taump = fdyn->tau[1];
  switch (flag)
  {
      case 0: /* evaluation at n */
        c = fdyn->thsr;
        break;
      case 1: /* evaluation at n+1 */
        c = fdyn->thsl;
        break;
      default:
        dserror("value of flag not valid!!!\n");
  }

  /* calculate external/convective stab-forces of time force vector => */
  if (gls->iadvec!=0)
  {
    fact[0] = edead[0]*fac*taumu*c*DENS*DENS;
    fact[1] = edead[1]*fac*taumu*c*DENS*DENS;
    for (inode=0;inode<TWO*iel;inode++)
    {
      irow = index[inode];
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1];
      for (isd=0;isd<2;isd++)
      {
        eforce[irow] += aux*fact[isd];
        irow++;
      }
    } 
  }

  /* calculate external/viscous stab-forces of time force vector => */
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
    
    fvts = fac*visc*taump*sign*c*DENS;
    for (inode=0;inode<TWO*iel;inode++)
    {
      irow = index[inode];
      eforce[irow]   -= ((TWO*derxy2[0][inode] + derxy2[1][inode])*edead[0]
                         + derxy2[2][inode]*edead[1])*fvts;
      eforce[irow+1] -= ((TWO*derxy2[1][inode] + derxy2[0][inode])*edead[1]
                         + derxy2[2][inode]*edead[0])*fvts;
    }
  }
  
/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calstabexfv */ 



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calstabexfp(
  DOUBLE          *eforce,     
  DOUBLE         **derxy,       
  DOUBLE          *edead,  
  DOUBLE           fac,      
  INT              iel,
  INT              flag,
  DOUBLE           DENS
  ) 
{
  INT        inode;
  DOUBLE     c;
  DOUBLE     taump;
  DOUBLE     fact[2];
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabexfp");
#endif    
/*---------------------------------------------------------------------*/

  /* initialize */
  fdyn = alldyn[genprob.numff].fdyn;
  
  /* set some factors */
  taump = fdyn->tau[1];
  switch (flag)
  {
      case 0: /* evaluation at n */
        c = fdyn->thpr;
        break;
      case 1: /* evaluation at n+1 */
        c = fdyn->thpl;
        break;
      default:
        dserror("value of flag not valid!!!\n");
  }

  /* calculate external/pressure stab forces of time force vector */
  fact[0] = edead[0]*taump*fac*c*DENS;
  fact[1] = edead[1]*taump*fac*c*DENS;
  for (inode=0;inode<iel;inode++)
  {
    eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1]);
  }
  
/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_stabexfp */ 
#endif
