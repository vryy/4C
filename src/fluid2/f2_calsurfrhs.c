/*!----------------------------------------------------------------------
\file
\brief rhs part due to surface tension effects

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
/*!---------------------------------------------------------------------
\brief rhs part due to surface tension effects

<pre>                                                         genk 02/03

In this routine the galerkin part of the forces due to surface tension
is calulated:

                 /
     + (facs)*  |  v * h^   d_gamma
               /  
 	  		      

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:
   h^ is either h at (n) or h at (n+1)
   facs is either (1-THETA)*dt at (n) or THETA*dt at (n+1)
		     
</pre>
\param  *eforce    double   (o)    element force vector
\param  *funct     double   (i)    natural shape functions
\param  *vn 	   double   (i)    normal vector on surface at gp
\param   sigmaint  double   (i)    surface tension at gauss point
\param   facs 	   double   (i)    time factor
\param   fac 	   double   (i)    weighting factor
\param   ngnode    int      (i)    number of nodes on actual edge
\param   iedgnod   int	    (i)    actual edge nodes
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calsurftenfv( 
                     double   *eforce, 
		     double   *funct, 
		     double   *vn, 
		     double    sigmaint,
		     double    facs, 
		     double    fac,
                     int       ngnode,
		     int      *iedgnod     
		    )
{
int irow, inode,isd;
double c;

#ifdef DEBUG 
dstrc_enter("f2_calsurftenfv");
#endif
/*--------------------------------------------------- set some factors */
c = facs*fac*sigmaint;

/*----------------------------------------------------------------------*
   Calculate galerkin part of external forces:
   
           /                      /
   + facs |  v * h   d_gamma  =  | v * n * sigma d_gamma 
         /                      /
   
 *----------------------------------------------------------------------*/
for (inode=0;inode<ngnode;inode++)
{
   irow = NUM_F2_VELDOF*iedgnod[inode];
   for (isd=0;isd<2;isd++)
   {
      eforce[irow] += funct[inode]*vn[isd]*c;
      irow++;
   } /* end of loop over isd */
} /* end of loop over inode */ 

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calsurftenfv */
#endif
/*! @} (documentation module close)*/
