/*!----------------------------------------------------------------------
\file
\brief routine to evaluate iteration rhs on element base resulting from
mass and acceleration in generalised alpha time integration

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID3 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"
#include "../fluid_full/fluid_prototypes.h"

/*!--------------------------------------------------------------------- 
\brief calculates mass-acceleration-part of the rhs vector for 
       generalised-alpha time integration

<pre>                                                        chfoe 11/03

This routine performs the multiplication emass*acc(n) and writes its
negative result on the elemental iteration force vector. This gives a 
part of the rhs required in the generalised alpha time integration scheme
for fluids.
			     
</pre>
\param	 *ele		ELEMENT		(i)	actual element
\param	**emass		DOUBLE		(i)	element mass matrix
\param  **eaccn		DOUBLE		(i)	element nodal accelerations
\param 	 *eiforce	DOUBLE	        (o)	element iteration force vec

\warning This routine has not yet been tested!!!

\return void                                               
\sa calling: 
             called by: f3_calele()
                     
------------------------------------------------------------------------*/
void f3_massrhs(ELEMENT *ele, DOUBLE **emass, DOUBLE **eaccn, DOUBLE *eiforce)
{
INT 	i,j;		/* counters					*/
INT 	iel;		/* element node number				*/
INT	dim=3;		/* dimension, number of velocity-dofs		*/
INT	dofpern=dim+1;	/* degrees of freedom per node			*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("f3_massrhs");
#endif

/*------------------------------------------------- set initial data ---*/
iel=ele->numnp;

/*------------------------------------------- perform multiplication ---*/
for (i=0; i<iel; i++)
{
   for (j=0; j<iel; j++)
   {
      eiforce[i*dofpern]   += (emass[i*dim][j*dim]   * eaccn[0][j]
                             + emass[i*dim][j*dim+1] * eaccn[1][j]
                             + emass[i*dim][j*dim+2] * eaccn[2][j]);
      eiforce[i*dofpern+1] += (emass[i*dim+1][j*dim]   * eaccn[0][j]
                             + emass[i*dim+1][j*dim+1] * eaccn[1][j]
                             + emass[i*dim+1][j*dim+2] * eaccn[2][j]);
      eiforce[i*dofpern+2] += (emass[i*dim+2][j*dim]   * eaccn[0][j]
                             + emass[i*dim+2][j*dim+1] * eaccn[1][j]
                             + emass[i*dim+2][j*dim+2] * eaccn[2][j]);
      eiforce[i*dofpern+3] += (emass[iel*dim+i][j*dim]   * eaccn[0][j]
                             + emass[iel*dim+i][j*dim+1] * eaccn[1][j]
                             + emass[iel*dim+i][j*dim+2] * eaccn[2][j]);
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_massrhs */



#endif
/*! @} (documentation module close)*/
