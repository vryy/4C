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
void xfem_f2_calstabkvv(			      
  ELEMENT         *ele,    
  DOUBLE         **estif,  
  DOUBLE          *velint,
  DOUBLE          *vel2int, 
  DOUBLE          *gridvint,
  DOUBLE         **vderxy, 
  DOUBLE          *funct,  
  DOUBLE         **derxy,  
  DOUBLE         **derxy2, 
  DOUBLE           fac,    
  DOUBLE           visc,   
  INT              iel,    
  INT              ihoel,
  INT             *index 
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |  
 *----------------------------------------------------------------------*/
  INT        irow,icol,irn,icn;
  DOUBLE     taumu;
  DOUBLE     taump;
  DOUBLE     tauc;
  DOUBLE     c,cc;
  DOUBLE     aux,auxr,auxc;
  DOUBLE     sign;

  STAB_PAR_GLS  *gls;	/* pointer to GLS stabilisation parameters	*/  
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabkvv");
#endif
/*----------------------------------------------------------------------*/

  /* initialise */
  fdyn = alldyn[genprob.numff].fdyn;
  gls  = ele->e.f2->stabi.gls;

  /* set stabilisation parameter */
  taumu = fdyn->tau[0];
  taump = fdyn->tau[1];
  tauc  = fdyn->tau[2];

  /* calculate continuity stabilisation part */
  if (gls->icont!=0)
  {
    c = fac*tauc;
    for(icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];
      for(irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];
        estif[irow  ][icol  ] += derxy[0][icn]*derxy[0][irn]*c;
        estif[irow+1][icol  ] += derxy[0][icn]*derxy[1][irn]*c;
        estif[irow  ][icol+1] += derxy[1][icn]*derxy[0][irn]*c;
        estif[irow+1][icol+1] += derxy[1][icn]*derxy[1][irn]*c;
      }
    }
  }
  
  c = fac*taumu;
  cc = c;
  /* calculate advection stabilisation part */
  if (gls->iadvec!=0)
  {
    /* evaluate for Newton- and fixed-point-like-iteration */
    if (fdyn->nic!=0)
    {
      for (icn=0; icn<TWO*iel; icn++)
      {
        auxc = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn])*cc;
        icol = index[icn];
        for (irn=0; irn<TWO*iel; irn++)
        {
          irow = index[irn];
          aux = (vel2int[0]*derxy[0][irn] + vel2int[1]*derxy[1][irn])*auxc;
          estif[irow  ][icol  ] += aux;
          estif[irow+1][icol+1] += aux;
        }
      }
    }
    
    /* evaluate for Newton iteration */
    if (fdyn->nir!=0) 
    {
      for (icn=0; icn<TWO*iel; icn++)
      {
        auxc = funct[icn]*cc;
        icol = index[icn];
        for (irn=0; irn<TWO*iel; irn++)
        {
          irow = index[irn];          
          aux = (vel2int[0]*derxy[0][irn] + vel2int[1]*derxy[1][irn])*auxc;
          estif[irow  ][icol  ] += aux*vderxy[0][0];
          estif[irow+1][icol  ] += aux*vderxy[1][0];
          estif[irow  ][icol+1] += aux*vderxy[0][1];
          estif[irow+1][icol+1] += aux*vderxy[1][1];
        }
      }
    }

    /* calculate advection stabilisation part for higher order elements */
    if (ihoel!=0)
    {
      cc = c*visc;
      for (icn=0; icn<TWO*iel; icn++)
      {
        icol = index[icn];        
        auxc = derxy2[0][icn] + derxy2[1][icn];
        for (irn=0; irn<TWO*iel; irn++)
        {
          irow = index[irn];                    
          aux = (vel2int[0]*derxy[0][irn] + vel2int[1]*derxy[1][irn])*cc;
          estif[irow  ][icol  ] -= aux*(derxy2[0][icn] + auxc);
          estif[irow+1][icol  ] -= aux* derxy2[2][icn];
          estif[irow+1][icol+1] -= aux*(derxy2[1][icn] + auxc);
          estif[irow  ][icol+1] -= aux* derxy2[2][icn];
        }
      }
    }
  }
  
  /* calculate viscous stabilisation part */
  if (ihoel!=0 && gls->ivisc!=0)
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

    /* calculate viscous stabilisation part for higher order elements */
    cc = fac * taump * visc*visc * sign;
    
    for (icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];              
      auxc = derxy2[0][icn] + derxy2[1][icn];
      for (irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];        
        auxr = derxy2[0][irn] + derxy2[1][irn];
        aux  = auxc*auxr;
        estif[irow  ][icol  ] +=  (derxy2[0][icn]*(derxy2[0][irn]+auxr) +
                                   derxy2[2][icn]* derxy2[2][irn] +
                                   derxy2[0][irn]* auxc + aux)*cc;
        estif[irow+1][icol  ] +=  (derxy2[0][icn]* derxy2[2][irn] +
                                   derxy2[2][icn]*(derxy2[1][irn]+auxr) +
                                   derxy2[2][irn]* auxc)*cc;  
        estif[irow+1][icol+1] +=  (derxy2[2][icn]* derxy2[2][irn] +
                                   derxy2[1][icn]*(derxy2[1][irn]+auxr) +
                                   derxy2[1][irn]* auxc + aux)*cc;
        estif[irow  ][icol+1] +=  (derxy2[2][icn]*(derxy2[0][irn]+auxr) +
                                   derxy2[1][icn]* derxy2[2][irn] +
                                   derxy2[2][irn]* auxc)*cc;
      }
    }

    /* calculate viscous stabilisation part Nc(u) for higher order elements */
    cc = fac * taump * visc * sign;
    
    if (fdyn->nic!=0) /* evaluate for Newton- and fixed-point-like-iteration */
    {
      for (icn=0; icn<TWO*iel; icn++)
      {
        icol = index[icn];              
        aux = velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn];
        for (irn=0; irn<TWO*iel; irn++)
        {
          irow = index[irn];                  
          auxr = derxy2[0][irn] + derxy2[1][irn];
          estif[irow  ][icol  ] -= (derxy2[0][irn]+auxr)*aux*cc;
          estif[irow+1][icol  ] -=  derxy2[2][irn]*aux*cc;
          estif[irow  ][icol+1] -=  derxy2[2][irn]*aux*cc;
          estif[irow+1][icol+1] -= (derxy2[1][irn]+auxr)*aux*cc;
        }
      }
    }

    /* calculate viscous stabilisation part Nr(u) for higher order elements */
    if (fdyn->nir!=0)  /* evaluate for Newton iteraton */
    {
      for (icn=0; icn<TWO*iel; icn++)
      {
        icol = index[icn];              
        aux = funct[icn]*cc;
        for (irn=0; irn<TWO*iel; irn++)
        {
          irow = index[irn];                  
          auxr = derxy2[0][irn]*derxy2[1][irn];
          estif[irow  ][icol  ] -=  (derxy2[0][irn]*vderxy[0][0] +
                                     derxy2[2][irn]*vderxy[1][0] +
                                     auxr*vderxy[0][0])*aux;
          estif[irow+1][icol  ] -=  (derxy2[2][irn]*vderxy[0][0] +
                                     derxy2[1][irn]*vderxy[1][0] +
                                     auxr*vderxy[1][0])*aux;
          estif[irow  ][icol+1] -=  (derxy2[0][irn]*vderxy[0][1] +
                                     derxy2[2][irn]*vderxy[1][1] +
                                     auxr*vderxy[0][1])*aux;
          estif[irow+1][icol+1] -=  (derxy2[2][irn]*vderxy[0][1] +
                                     derxy2[1][irn]*vderxy[1][1] +
                                     auxr*vderxy[1][1])*aux;
        }
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of f2_calstabkvv */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calstabkvp(
  ELEMENT         *ele,    
  DOUBLE         **estif, 
  DOUBLE          *velint,
  DOUBLE          *funct, 
  DOUBLE         **derxy, 
  DOUBLE         **derxy2,
  DOUBLE           fac,   
  DOUBLE           visc,  
  INT              iel,   
  INT              ihoel,
  INT             *index   
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |
 |   posc - since there's only one full element stiffness matrix the    |
 |          column number has to be changed!                            |   
/*----------------------------------------------------------------------*/
  INT        irow,icol,irn,posc;
  DOUBLE     taumu;
  DOUBLE     taump;
  DOUBLE     c;
  DOUBLE     aux;
  DOUBLE     sign;

  STAB_PAR_GLS  *gls;	/* pointer to GLS stabilisation parameters	*/  
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabkvp");
#endif
/*----------------------------------------------------------------------*/

  /* initialise */
  fdyn = alldyn[genprob.numff].fdyn;
  gls  = ele->e.f2->stabi.gls;

  /* set stabilisation parameter */
  taumu = fdyn->tau[0];
  taump = fdyn->tau[1];

  c = fac * taumu;
  /* calculate advection stabilisation part */
  if (gls->iadvec!=0)
  {
    for (icol=0; icol<iel; icol++)
    {
      posc=2*iel+icol;
      for (irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];
        aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*c;
        estif[irow  ][posc] += derxy[0][icol]*aux;
        estif[irow+1][posc] += derxy[1][icol]*aux;
      }
    }
  }
  /* calculate viscous stabilisation part */
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
    c = fac * taump * visc * sign;
    
    for (icol=0; icol<iel; icol++)
    {
      posc=2*iel+icol;
      for (irn=0; irn<TWO*iel; irn++)
      {
        irow=index[irn];
        aux = derxy2[0][irn] + derxy2[1][irn];
        estif[irow  ][posc] -=  (derxy2[0][irn]*derxy[0][icol] +
                                 derxy2[2][irn]*derxy[1][icol] +
                                 aux*derxy[0][icol])*c;
        estif[irow+1][posc] -=  (derxy2[2][irn]*derxy[0][icol] +
                                 derxy2[1][irn]*derxy[1][icol] +
                                 aux*derxy[1][icol])*c;
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calstabkvp */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calstabmvv(
  ELEMENT         *ele,     
  DOUBLE         **emass,  
  DOUBLE          *velint, 
  DOUBLE          *funct,  
  DOUBLE         **derxy,  
  DOUBLE         **derxy2, 
  DOUBLE           fac,    
  DOUBLE           visc,   
  INT              iel,    
  INT              ihoel,
  INT             *index,
  DOUBLE           dens
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |   
/*----------------------------------------------------------------------*/
  INT        irow,icol,irn,icn;
  DOUBLE     taumu;
  DOUBLE     taump;
  DOUBLE     c,cc;
  DOUBLE     aux,auxc;
  DOUBLE     sign;

  STAB_PAR_GLS  *gls;	/* pointer to GLS stabilisation parameters	*/  
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabmvv");
#endif
/*----------------------------------------------------------------------*/

  /* initialise */
  fdyn = alldyn[genprob.numff].fdyn;
  gls  = ele->e.f2->stabi.gls;

  /* set stabilisation parameter */
  taumu = fdyn->tau[0];
  taump = fdyn->tau[1];
  
  c = fac * taumu * dens;
  cc = c;

  /* calculate advection stabilisation part */
  if (gls->iadvec!=0)
  {
    for (icn=0; icn<TWO*iel; icn++)
    {
      auxc = funct[icn]*cc;
      icol = index[icn];
      for (irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];
        aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*auxc;
        emass[irow  ][icol  ] += aux;
        emass[irow+1][icol+1] += aux;
      }
    }
  }

  /* calculate viscous stabilisation part */
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
    c = fac * taump * visc * sign;
    
    for (icn=0; icn<TWO*iel; icn++)
    {      
      aux = funct[icn]*c;
      icol = index[icn];
      for(irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];
        emass[irow  ][icol  ] -= (TWO*derxy2[0][irn] + derxy2[1][irn])*aux;
        emass[irow+1][icol  ] -=      derxy2[2][irn]*aux;
        emass[irow+1][icol+1] -= (TWO*derxy2[1][irn] + derxy2[0][irn])*aux;
        emass[irow  ][icol+1] -=      derxy2[2][irn]*aux;
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
 return;
} /* end of xfem_f2_calstabmvv */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calstabkpv(
  ELEMENT         *ele,
  DOUBLE         **estif,   
  DOUBLE          *velint,
  DOUBLE          *gridvint, 
  DOUBLE         **vderxy, 
  DOUBLE          *funct,  
  DOUBLE         **derxy,  
  DOUBLE         **derxy2, 
  DOUBLE           fac,    
  DOUBLE           visc,   
  INT              iel,    
  INT              ihoel,
  INT             *index   
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |
 |   posr - since there's only one full element stiffness matrix the    |
 |          row number has to be changed!                               |
/*----------------------------------------------------------------------*/
  INT        irow,icol,icn,posr;
  DOUBLE     c;
  DOUBLE     aux;
  DOUBLE     taump;

#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabkpv");
#endif
/*----------------------------------------------------------------------*/

  /* initialise */
  fdyn = alldyn[genprob.numff].fdyn;
  
  /* set stabilisation parameter */
  taump = fdyn->tau[1];
  
  c = fac * taump;
  /* calculate stabilisation part Nc(u) */
  /* evaluate for Newton- and fixed-point-like-iteration */
  if (fdyn->nic!=0)
  {
    for (icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];
      aux = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn])*c;
      for (irow=0; irow<iel; irow++)
      {
        posr = irow + 2*iel;
        estif[posr][icol  ] -= derxy[0][irow]*aux;
        estif[posr][icol+1] -= derxy[1][irow]*aux;
      }
    }
  }
  /* calculate stabilisation part Nr(u) */
  if (fdyn->nir!=0) /* evaluate for Newton iteration */
  {
    for (icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];
      aux = funct[icn]*c;
      for (irow=0; irow<iel; irow++)
      {
        posr = irow + 2*iel;
        estif[posr][icol  ] -= aux*(derxy[0][irow]*vderxy[0][0] +
                                    derxy[1][irow]*vderxy[1][0]);
        estif[posr][icol+1] -= aux*(derxy[0][irow]*vderxy[0][1] +
                                    derxy[1][irow]*vderxy[1][1]);
      }
    }
  }

  /* calculate stabilisation part for higher order elements */
  if (ihoel!=0)
  {
    c = c * visc;
    
    for (icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];
      aux = derxy2[0][icn] + derxy2[1][icn];
      for (irow=0; irow<iel; irow++)
      {
        posr = irow + 2*iel;
        estif[posr][icol  ] += ((derxy2[0][icn]+aux)*derxy[0][irow] +
                                derxy2[2][icn]     *derxy[1][irow])*c;
        estif[posr][icol+1] +=  (derxy2[2][icn]     *derxy[0][irow] +
                                 (derxy2[1][icn]+aux)*derxy[1][irow])*c;
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calstabkpv */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calstabmpv(
  DOUBLE         **emass,   
  DOUBLE          *funct,  
  DOUBLE         **derxy,  
  DOUBLE           fac,    
  INT              iel,
  INT             *index,
  DOUBLE           dens
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |
 |   posr - since there's only one full element stiffness matrix the    |
 |          row number has to be changed!                               |
/*----------------------------------------------------------------------*/
  INT        irow,icol,icn,posr;
  DOUBLE     c;
  DOUBLE     taump;
  DOUBLE     auxc;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabmpv");
#endif
/*----------------------------------------------------------------------*/

  /* initialise */
  fdyn = alldyn[genprob.numff].fdyn;
  
  /* set stabilisation parameter */
  taump = fdyn->tau[1];
  
  c = fac * taump * dens;
  /* calculate stabilisation part for matrix Mpv */
  for (icn=0; icn<TWO*iel; icn++)
  {
    icol = index[icn];
    auxc = funct[icn]*c;
    for (irow=0; irow<iel; irow++)
    {
      posr = irow + 2*iel;
      emass[posr][icol  ] -= derxy[0][irow]*auxc;
      emass[posr][icol+1] -= derxy[1][irow]*auxc;
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calstabmpv */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calstabkpp(
  DOUBLE         **estif,   
  DOUBLE         **derxy,  
  DOUBLE           fac,    
  INT              iel
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   posr - since there's only one full element stiffness matrix the    |
 |          row number has to be changed!                               |
 |   posc - since there's only one full element stiffness matrix the    |
 |          column number has to be changed!                            |
/*----------------------------------------------------------------------*/
  INT        irow,icol,posc,posr;
  DOUBLE     c;
  DOUBLE     taump;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabkpp");
#endif
/*----------------------------------------------------------------------*/

  /* initialise */
  fdyn = alldyn[genprob.numff].fdyn;

  /* set stabilisation parameter */
  taump = fdyn->tau[1];
  
  c = fac * taump;
  /* calculate stabilisation part for matrix Kpp */
  for (icol=0; icol<iel; icol++)
  {
    posc = icol + 2*iel;
    for (irow=0; irow<iel; irow++)
    {
      posr = irow + 2*iel;
      estif[posr][posc] -= (derxy[0][irow]*derxy[0][icol] +
                            derxy[1][irow]*derxy[1][icol])*c;
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calstabkpp */
#endif
