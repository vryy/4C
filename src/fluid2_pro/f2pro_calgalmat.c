/*!----------------------------------------------------------------------
\file
\brief evaluate galerkin part of stiffness matrix

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2_PRO 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2_PRO 
#include "../headers/standardtypes.h"
#include "fluid2pro_prototypes.h"
/*!----------------------------------------------------------------------                
\brief evaluate galerkin part of Kvv

<pre>                                                         basol 11/02

In this routine the galerkin part of matrix Kvv is calculated:
    /
   |  2 * nue * eps(v) : eps(u)   d_omega
  /
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
NOTE: there's only one elestif
      --> Kvv is stored in estif[0..(2*iel-1)][0..(2*iel-1)]
</pre>
\param **estif     DOUBLE	   (i/o)  ele stiffness matrix
\param **derxy     DOUBLE	   (i)    global coord. deriv.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   dt        DOUBLE	   (i)    incremental time step	
\param   iel	   INT  	   (i)	  number of nodes of act. ele
\return void                                                                       
------------------------------------------------------------------------*/
void f2pro_calkvv(
                
		DOUBLE         **estif,   
		DOUBLE         **derxy,  
		DOUBLE           fac,    
		DOUBLE           visc,
		DOUBLE           dt,   
		INT              iel     
              ) 
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |  
/*----------------------------------------------------------------------*/
INT     irow, icol,irn,icn;
DOUBLE  c,aux;

#ifdef DEBUG 
dstrc_enter("f2pro_calkvv");
#endif		

c=fac*visc;

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix K:
    /
   |  2 * nue * eps(v) : eps(u)   d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
for (icn=0;icn<iel;icn++) 
{
   irow=0;
   for (irn=0;irn<iel;irn++)
   {
      estif[irow][icol]     += c*(TWO*derxy[0][irn]*derxy[0][icn] \
                                    + derxy[1][irn]*derxy[1][icn]); 	     
      estif[irow+1][icol]   += c*(    derxy[0][irn]*derxy[1][icn]);
      estif[irow+1][icol+1] += c*(TWO*derxy[1][irn]*derxy[1][icn] \
                                    + derxy[0][irn]*derxy[0][icn]);
      estif[irow][icol+1]   += c*(    derxy[1][irn]*derxy[0][icn]);
      irow += 2;
   } /* end loop over irn */
   icol += 2;
} /* end loop over icn */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2pro_calkvv */

/*!----------------------------------------------------------------------                 
\brief evaluate the lumped mass matrix
<pre>                                                         basol 11/02

In this routine lumped mass matrix is calculated

</pre>
\param **lmass     DOUBLE	   (o)    lumped ele mass matrix
\param **emass     DOUBLE	   (i)    ele mass matrix
\param   iel	   INT  	   (i)	  number of nodes
\return void                                                                       

------------------------------------------------------------------------*/
void f2pro_lmass(        
		 DOUBLE         **lmass,   
		 DOUBLE         **emass,  
		 INT                iel     
              )
{
INT     i,j;
DOUBLE sum;

#ifdef DEBUG 
dstrc_enter("f2pro_lmass");
#endif		

for (i=0;i<2*iel;i++)
{
   sum=ZERO;
   for(j=0;j<2*iel;j++)
   {
      sum +=emass[i][j];
      lmass[i][j]=ZERO; 
   }
   lmass[i][i]=sum;
}   

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2pro_lmass */

/*!---------------------------------------------------------------------                 
\brief evaluate the gradient operator                         

<pre>                                                   basol 11/02      
In this routine the gradient operator is calculated

psi,x = derivative of shape function wrt x for velocities
psi,y = derivative of shape function wrt y for velocities
phi   = shape function for pressure 
i,j   = index
        
           /
 Cxx = (-)|  (psi,x)(i)*phi(j) d_omega
         /
           /
 Cyy = (-)|  (psi,y)(i)*phi(j) d_omega
         /
    
 C = [Cxx, 
      Cyy]   
 C --->18x4
 CT --->4x18
  
</pre>
\param **gradopr    DOUBLE	   (o)    grad. operator C-->(18x4)
\param **derxy     DOUBLE	   (i)    global coord. deriv.
\param  *functpr    DOUBLE	   (i)    pressure shape functions	      
\param   fac 	   DOUBLE	   (i)    weighting factor	     
\param   ielp	   INT  	   (i)	  number of nodes for pressure ele.
\param   iel	   INT  	   (i)	  number of nodes for velocity ele.
\return void                                                           
------------------------------------------------------------------------*/
void f2pro_gradopr(
                
		DOUBLE          **gradopr,   
		DOUBLE          **derxy,  
		DOUBLE           *functpr,    
		DOUBLE            fac,
		INT               ielp,    
		INT               iel   
              )
{
INT     i,j;
INT     irow,icol;

#ifdef DEBUG 
dstrc_enter("f2pro_ gradopr");
#endif		

irow=0;
for (i=0;i<iel;i++)
{
   icol=0;
   for (j=0;j<ielp;j++)
   {
      gradopr[irow][icol] -= derxy[0][i]*functpr[j]*fac;
      gradopr[irow+1][icol] -= derxy[1][i]*functpr[j]*fac;
      icol ++;
   }
   irow +=2;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2pro_gradopr */

#endif
/*! @} (documentation module close)*/
