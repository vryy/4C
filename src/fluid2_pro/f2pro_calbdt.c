/*!----------------------------------------------------------------------
\file
\brief balancing diffusitivity tensor

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

/*!---------------------------------------------------------------------
\brief calculate the balancing diffusivity tensor

<pre>                                                       basol 03/03

    /
   |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega
  /
  
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param **derxy     DOUBLE	   (i)     global derivatives
\param  *func      DOUBLE	   (i)     shape functions
\param   fac	   DOUBLE	   (i)	   weighting factor	   
\param   visc	   DOUBLE	   (i)	   fluid viscosity
\param   iel	   INT		   (i)	   num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2pro_calbdt(			      
                    ELEMENT         *ele, 
		    DOUBLE         **estif,  
		    DOUBLE          *velint, 
		    DOUBLE         **derxy,  
		    DOUBLE         **vderxy,
                    DOUBLE          *funct,
                    DOUBLE           fac,    
		    DOUBLE          visc,   
		    INT               iel    
                   )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
/*----------------------------------------------------------------------*/
INT    irow,icol,irn,icn;
DOUBLE taumu;
DOUBLE c;
DOUBLE aux,auxc;
FLUID_DYNAMIC *fdyn;    /* pointer to fluid dynamic variables           */

#ifdef DEBUG 
dstrc_enter("f2pro_calbtd");
#endif

/*---------------------------------------- set stabilisation parameter */
/*---stabilization parameter is taken to be dt/2 for the BTD method----*/
fdyn = alldyn[genprob.numff].fdyn;
taumu = fdyn->dta/TWO;
c = fac*taumu;

/*----------------------------------------------------------------------*
   Calculate advection stabilisation part Nc(u):
    /
   |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega
  /
*----------------------------------------------------------------------*/
icol=0;
for (icn=0;icn<iel;icn++)
{
   auxc = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn])*c;
   irow=0;
   for (irn=0;irn<iel;irn++)
   {
      aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*auxc;
      estif[irow][icol]     += aux;
      estif[irow+1][icol+1] += aux;
      irow += 2;
   } /* end of loop over irn */
   icol += 2;
} /* end of loop over icn */

/*---------------------------------------------------------------------*/   
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2pro_calbdt */

#endif
/*! @} (documentation module close)*/
