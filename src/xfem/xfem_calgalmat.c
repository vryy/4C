#ifdef D_XFEM
#include "../headers/standardtypes.h"
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
void xfem_f2_calkvv(
  ELEMENT         *ele,
  DOUBLE         **estif,   
  DOUBLE          *velint,
  DOUBLE          *gridvint,
  DOUBLE         **vderxy, 
  DOUBLE          *funct,  
  DOUBLE         **derxy,  
  DOUBLE           fac,    
  DOUBLE           visc,   
  INT              iel,
  INT             *index,
  DOUBLE           DENS
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |  
/*----------------------------------------------------------------------*/
  INT        irow,icol,irn,icn;
  DOUBLE     c,aux;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calkvv");
#endif		
/*----------------------------------------------------------------------*/

  /* initialize */
  fdyn   = alldyn[genprob.numff].fdyn;
  
  /*
   * => part I
   */

  /* set c */
  c=fac*visc;
  
  for (icn=0; icn<TWO*iel; icn++) 
  {
    icol = index[icn];
    for (irn=0; irn<TWO*iel; irn++)
    {
      irow = index[irn];      
      estif[irow  ][icol  ] += c*(TWO*derxy[0][irn]*derxy[0][icn] +
                                      derxy[1][irn]*derxy[1][icn]); 	     
      estif[irow+1][icol  ] += c*(    derxy[0][irn]*derxy[1][icn]);
      estif[irow+1][icol+1] += c*(TWO*derxy[1][irn]*derxy[1][icn] +
                                      derxy[0][irn]*derxy[0][icn]);
      estif[irow  ][icol+1] += c*(    derxy[1][irn]*derxy[0][icn]);
    }
  }

  /*
   * => part II
   */

  /* evaluate for Newton- and fixed-point-like-iteration */
  if(fdyn->nic!=0)
  {
    for (icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];  
      for (irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];
        aux = (velint[0]*derxy[0][icn] +
               velint[1]*derxy[1][icn])*funct[irn]*fac*DENS;
        estif[irow  ][icol  ] += aux;
        estif[irow+1][icol+1] += aux;
      }
    }
  }
  
  /*
   * => part III
   */

  /* evaluate for Newton iteration */
  if (fdyn->nir!=0) 
  {
    for (icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];  
      for (irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];
        aux = funct[irn]*funct[icn]*fac*DENS;
        estif[irow  ][icol  ] += aux*vderxy[0][0];
        estif[irow+1][icol  ] += aux*vderxy[1][0];
        estif[irow+1][icol+1] += aux*vderxy[1][1];
        estif[irow  ][icol+1] += aux*vderxy[0][1];
      }
    }
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calkvv */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calkvp(
  DOUBLE         **estif,   
  DOUBLE          *funct,  
  DOUBLE         **derxy,  
  DOUBLE           fac,    
  INT              iel,
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
  INT        irow,icol,irn,ird;  
  INT        posc;
  DOUBLE     aux;

#ifdef DEBUG 
  dstrc_enter("xfem_f2_calkvp");
#endif		
/*----------------------------------------------------------------------*/

  for (icol=0; icol<iel; icol++)
  {
    posc = icol + 2*iel;
    for (irn=0; irn<TWO*iel; irn++)
    {
      irow = index[irn];
      for(ird=0; ird<2; ird++)
      {      
	aux = -funct[icol]*derxy[ird][irn]*fac;
	estif[irow][posc] += aux;
	estif[posc][irow] += aux;
        irow++;
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calkvp */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calmvv(
  DOUBLE         **emass,  
  DOUBLE          *funct, 
  DOUBLE           fac,   
  INT              iel,
  INT             *index,
  DOUBLE           DENS
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |  
/*----------------------------------------------------------------------*/
  INT        irow,icol,irn,icn;  
  DOUBLE     aux;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calmvv");
#endif		
/*----------------------------------------------------------------------*/

  /* calculate full Galerkin part of matrix Mvv */
  for (icn=0; icn<TWO*iel; icn++) 
  {
    icol = index[icn];
    for (irn=0; irn<TWO*iel; irn++)
    {
      irow = index[irn];      
      aux = funct[icn]*funct[irn]*fac*DENS;
      emass[irow  ][icol  ] += aux;
      emass[irow+1][icol+1] += aux;   
    }
  }

/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of xfem_f2_calmvv */
#endif
