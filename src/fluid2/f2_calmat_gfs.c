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
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
/*!---------------------------------------------------------------------
\brief evaluate galerkin part of Kgv and Kgg

<pre>                                                         genk 01/03

      
</pre>
\param **estif     DOUBLE     (i/o)  ele stiffness matrix
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calmat_gfs(DOUBLE       *funct,
                   DOUBLE       *vn,
                   DOUBLE      **estif,
                   DOUBLE        fac,
                   DOUBLE        length,
                   INT           iel,
                   INT           ngnode,
                   INT          *iedgnod
                 )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |  
 *----------------------------------------------------------------------*/
INT     irow,icolv,icolg1,icolg2,irn,icn;  
INT     nd;                 
DOUBLE  aux;

#ifdef DEBUG 
dstrc_enter("f2_calmat_gfs");
#endif		 

nd = NUMDOF_FLUID2*iel;

/*----------------------------------------------------------------------*
       /
    + |  psi * u_x*n_x    d_gamma_FS
     /
       /
    + |  psi * u_y*n_y    d_gamma_FS
     /

 *----------------------------------------------------------------------*/
for(icn=0;icn<ngnode;icn++)
{
   icolv = NUM_F2_VELDOF*iedgnod[icn];
   for(irn=0;irn<ngnode;irn++)
   {
      irow = iedgnod[irn]+nd;
      estif[irow][icolv]   += funct[irn]*funct[icn]*vn[0]*fac;
      estif[irow][icolv+1] += funct[irn]*funct[icn]*vn[1]*fac;
   } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*
       /
    - |  psi * u_G_x*n_x    d_gamma_FS
     /
       /
    - |  psi * u_G_y*n_y    d_gamma_FS
     /

 *----------------------------------------------------------------------*/
for(icn=0;icn<ngnode;icn++)
{
   icolg1 = nd+iedgnod[icn];
   icolg2 = nd+iel+iedgnod[icn];
   for(irn=0;irn<ngnode;irn++)
   {
      irow = iedgnod[irn]+nd;
      estif[irow][icolg1] -= funct[irn]*funct[icn]*vn[0]*fac;
      estif[irow][icolg2] -= funct[irn]*funct[icn]*vn[1]*fac;
   } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*
       /
    + |  zeta * c * uG_x    d_gamma_FS
     /
 *----------------------------------------------------------------------*/
for(icn=0;icn<ngnode;icn++)
{
   icolg1 = nd+iedgnod[icn];
   for(irn=0;irn<ngnode;irn++)
   {
      irow = iedgnod[irn]+nd+iel;
      estif[irow][icolg1] += funct[irn]*funct[icn]/DSQR(length)*fac;
   } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calkgedge */

#endif
/*! @} (documentation module close)*/
