/*!----------------------------------------------------------------------
\file
\brief evaluate galerkin part of submesh matrices for fluid3

<pre>
Maintainer:  name
             e-mail
             homepage
             telephone number


</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3ml_prototypes.h"
/*!---------------------------------------------------------------------                                         
\brief evaluate Galerkin part of submesh stiffness matrix SMK for fluid3

<pre>                                                       gravem 07/03

In this routine, the Galerkin part of the submesh stiffness matrix SMK 
is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smestif   DOUBLE	   (i/o)  submesh ele stiffness matrix
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun. 
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calsmk(FLUID_DYN_CALC  *dynvar,
	       FLUID_DYN_ML    *mlvar, 
	       DOUBLE         **smestif,   
	       DOUBLE          *velint, 
	       DOUBLE         **vderxy, 
	       DOUBLE          *smfunct,  
	       DOUBLE         **smderxy,  
	       DOUBLE           fac,    
	       DOUBLE           visc,   
	       INT              smiel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
/*----------------------------------------------------------------------*/
INT     irow,icol;
DOUBLE  con,aux,beta,divv;

#ifdef DEBUG 
dstrc_enter("f3_calsmk");
#endif		

/*-------------------------------------------------- subgrid viscosity */
if (mlvar->smsgvi>0) con=fac*(visc+mlvar->smsgvisc);
else con=fac*visc;

/*----------------------------------------------------------------------*
   Calculate viscous part of matrix SMK (including subgrid viscosity):
    /
   |  (nue+nueT) * grad(w) : grad(bub)   d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<smiel;icol++) 
{
  for (irow=0;irow<smiel;irow++)
  {
    smestif[irow][icol] += con*(smderxy[0][irow]*smderxy[0][icol]\
                               +smderxy[1][irow]*smderxy[1][icol]\
			       +smderxy[2][irow]*smderxy[2][icol]);		   
  } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*
    Calculate convective part of matrix SMK:
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  w * u_old[ls_u_old] * grad(bub)     d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<smiel;icol++)
{
  aux = (velint[0]*smderxy[0][icol]+velint[1]*smderxy[1][icol]\
        +velint[2]*smderxy[2][icol])*fac;
  for (irow=0;irow<smiel;irow++)
  {
    smestif[irow][icol] += aux*smfunct[irow];
  } /* end loop over irow */
}/* end loop over icol */ 
  
if (dynvar->conte!=0)  
{
/*----------------------------------------------------------------------*
    /
   | beta * w * bub * div(u_old[ls_u_old])   d_omega
  /
 *----------------------------------------------------------------------*/
  if (dynvar->conte==1) beta = ONE;
  else beta = ONE/TWO;
/*---------------------------------------------- divergence of velocity */
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
  
  for (icol=0;icol<smiel;icol++)
  {
    aux = beta*smfunct[icol]*divv*fac;
    for (irow=0;irow<smiel;irow++)
    {
      smestif[irow][icol] += aux*smfunct[irow];
    } /* end loop over irow */
  } /* end loop over icol */
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calsmk */

/*!---------------------------------------------------------------------                                         
\brief evaluate Galerkin part of submesh mass matrix SMM for fluid3

<pre>                                                       gravem 07/03

In this routine, the Galerkin part of the submesh mass matrix SMM is 
calculated.

</pre>
\param **smemass   DOUBLE	   (i/o)  submesh element mass matrix
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calsmm(DOUBLE         **smemass,   
	       DOUBLE          *smfunct,  
	       DOUBLE           fac,    
	       INT              smiel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
/*----------------------------------------------------------------------*/
INT     irow,icol;

#ifdef DEBUG 
dstrc_enter("f3_calsmm");
#endif		

/*----------------------------------------------------------------------*
   Calculate matrix SMM:
    /
   |  v * bub  d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<smiel;icol++) 
{
  for (irow=0;irow<smiel;irow++)
  {
    smemass[irow][icol] += smfunct[irow]*smfunct[icol]*fac;		   
  } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calsmm */

/*!---------------------------------------------------------------------                                         
\brief evaluate diffusive part of submesh stiffn. matrix SMK for fluid3

<pre>                                                       gravem 07/03

In this routine, the diffusive Galerkin part of the submesh stiffness 
matrix SMK is calculated.

</pre>
\param **smiediff  DOUBLE	   (i/o)  sm ele stiffn. matrix (diff.)
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun. 
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calsmkd(DOUBLE         **smiediff,   
        	DOUBLE         **smderxy,  
	        DOUBLE           fac,    
	        INT              smiel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
/*----------------------------------------------------------------------*/
INT     irow,icol;

#ifdef DEBUG 
dstrc_enter("f3_calsmkd");
#endif		

/*----------------------------------------------------------------------*
   Calculate diffusive (viscous) part of matrix SMK (inc. subgrid visc.):
    /
   |  grad(w) : grad(bub)   d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<smiel;icol++) 
{
  for (irow=0;irow<smiel;irow++)
  {
    smiediff[irow][icol] += fac*(smderxy[0][irow]*smderxy[0][icol]\
                                +smderxy[1][irow]*smderxy[1][icol]\
			        +smderxy[2][irow]*smderxy[2][icol]);		   
  } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calsmkd */

#endif
