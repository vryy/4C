/*!----------------------------------------------------------------------
\file
\brief iteration RHS for fluid2 element

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0771 - 685-6575
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID2TU
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2_tu.h"
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

static FLUID_DYNAMIC   *fdyn;
/*!---------------------------------------------------------------------
\brief galerkin part of iteration forces for kapome dofs

<pre>                                                        he  03/02

In this routine the galerkin part of the iteration forces for kapome dofs
is calculated:

                    /
          THETA*dt | factor * (kapome_old)^2 * psi  d_omega
                  /


</pre>
\param  *eforce      DOUBLE	    (i/o)   element force vector
\param   kapomeint   DOUBLE	     (i)    kapome at integr. point
\param  *funct       DOUBLE	     (i)    nat. shape funcs
\param   fac 	   DOUBLE	     (i)    weighting factor
\param   factor2 	   DOUBLE	     (i)    factor
\param   iel	   INT           (i)	num. of nodes of act. ele
\return void

------------------------------------------------------------------------*/
void f2_calgalifkapome(
                  DOUBLE          *eforce,
                  DOUBLE           kapomeint,
                  DOUBLE          *funct,
		      DOUBLE           fac,
                  DOUBLE           factor2,
                  INT              iel
                 )
{
INT    inode,isd;
INT    irow;
DOUBLE facsl;

#ifdef DEBUG
dstrc_enter("f2_calgalifkapome");
#endif

/*----------------------------------------------- set some factors */
fdyn  = alldyn[genprob.numff].fdyn;
facsl = fac * fdyn->thsl;

/*------------------------------------------------------------------*
   Calculate forces of iteration force vector:
             /
   THETA*dt |  factor * (kapome_old)^2 * psi  d_omega
           /
 *------------------------------------------------------------------*/
irow = -1;
for (inode=0;inode<iel;inode++)
{
      irow++;
      eforce[irow] += factor2 * pow(kapomeint,2) * funct[inode] * facsl;
} /* end loop over inode */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calgalifv */

/*!---------------------------------------------------------------------
\brief stabilisation part of iteration forces for kapome dofs

<pre>                                                         he  03/02

In this routine the stabilisation part of the iteration forces for kapome


           /
 THETA*dt | tau_tu * factor * (kapome_old)^2 * grad(psi) * u  d_omega + D.C.
         /



</pre>
\param   *ele      ELEMENT	   (i)    actual element
\param   *eforce   DOUBLE	   (i/o)  element force vector
\param    kapomeint DOUBLE	   (i)    kapome at integr. point
\param   *velint   DOUBLE	   (i)    vel at integr. point
\param   *velint_dc  DOUBLE	   (i)    vel at integr. point for D.C.
\param   *funct    DOUBLE	   (i)    nat. shape funcs
\param  **derxy    DOUBLE	   (i)    global derivative
\param    fac 	 DOUBLE	   (i)    weighting factor
\param    factor2  DOUBLE	   (i)    factor
\param    iel	   INT	   (i)    num. of nodes of act. ele
\return void

------------------------------------------------------------------------*/
void f2_calstabifkapome(
                        ELEMENT         *ele,
	                DOUBLE          *eforce,
		        DOUBLE           kapomeint,
                        DOUBLE          *velint,
                        DOUBLE          *velint_dc,
                        DOUBLE          *funct,
		        DOUBLE         **derxy,
		        DOUBLE           fac,
                        DOUBLE           factor2,
		        INT              iel
                      )
{
INT    inode,isd;
INT    irow;
DOUBLE facsl,facsl_dc;
DOUBLE aux,aux_dc;
DOUBLE taumu,taumu_dc;

#ifdef DEBUG
dstrc_enter("f2_calstabifkapome");
#endif

/*--------------------------------------------------- set some factors */
fdyn     = alldyn[genprob.numff].fdyn;
taumu    = fdyn->tau_tu;
taumu_dc = fdyn->tau_tu_dc;
facsl    = fac * fdyn->thsl * taumu;
facsl_dc = fac * fdyn->thsl * taumu_dc;
/*----------------------------------------------------------------------
   Calculate  stabilastion of iteration force vector:

           /
 THETA*dt |tau_tu * factor * (kapome_old)^2 * grad(psi) * u  d_omega
         /

*----------------------------------------------------------------------*/
 irow=0;
   for (inode=0;inode<iel;inode++)
   {
      aux    = derxy[0][inode]*velint[0]    + derxy[1][inode]*velint[1];
      aux_dc = derxy[0][inode]*velint_dc[0] + derxy[1][inode]*velint_dc[1];

      eforce[irow] += factor2*pow(kapomeint,2)*aux   *facsl;
      eforce[irow] += factor2*pow(kapomeint,2)*aux_dc*facsl_dc;

      irow ++;
   } /* end loop over inode */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calgalifkapeps */

#endif
