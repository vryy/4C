/*!----------------------------------------------------------------------
\file
\brief routine to evaluate mass right hand side

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
#include "../fluid_full/fluid_prototypes.h"
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
\brief calculates mass-rhs of history values

<pre>                                                        chfoe 09/03
                                    free surface included    chfoe 02/05

This routine performs the multiplication emass*(hist+deadload) and writes 
its result on the elemental iteration force vector. This replaces the 
time right hand side of WAW.

</pre>
\param	 *ele		ELEMENT		(i)	actual element
\param	**emass		DOUBLE		(i)	element mass matrix
\param  **hist		DOUBLE		(i)	elemental nodal history data
\param   *edeadng       DOUBLE          (i)     elemental body force at t=n+1
\param 	 *eiforce	DOUBLE	        (o)	element iteration force vec
\param   *hasext        INT             (i/o)   flag, if there's external load

\warning This routine sets *hasext to true in order to ensure assembly!

\return void
\sa calling:
             called by: f2_calele()

------------------------------------------------------------------------*/
void f2_massrhs(ELEMENT *ele, 
                DOUBLE **emass, 
                DOUBLE **hist,  
                DOUBLE  *edeadng, 
                DOUBLE  *eiforce,
                INT     *hasext)
{
INT 	i,j;		/* counters					*/
INT 	iel;		/* element node number				*/
INT	dim=2;		/* dimension, number of velocity-dofs		*/
INT	dofpern=dim+1;	/* degrees of freedom per node			*/
INT     row[MAXNOD_F2]; /* index field for free surface elements        */
DOUBLE  timefac;        /* factor from time integration                 */
FLUID_DYNAMIC   *fdyn;  /* pointer to fluid dynamic variables           */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("f2_massrhs");
#endif

/*------------------------------------------------- set initial data ---*/
iel     = ele->numnp;
fdyn    = alldyn[genprob.numff].fdyn;
timefac = fdyn->thsl;

/*----------------------------------------------------------------------*
 |       compute "external" Force Vector (b)                            |
 |  dead load may vary over time, but stays constant over               |
 |  the whole domain --> no nodal dependence                            |
 *----------------------------------------------------------------------*/
/*--------------------- add body force load to right hand side terms ---*/
if(*hasext)
   for(i=0; i<iel; i++)
   {
      hist[0][i] += timefac * edeadng[0];
      hist[1][i] += timefac * edeadng[1];
   }
else *hasext = 1; /* this assures rhs assembly! */

switch(ele->e.f2->fs_on)
{
case 0: case 1: case 3:
   /*---------------------------------------- perform multiplication ---*/
   for (i=0; i<iel; i++)
   {
      for (j=0; j<iel; j++)
      {
         eiforce[i*dofpern]   += (emass[i*dim][j*dim]   * hist[0][j]
                                + emass[i*dim][j*dim+1] * hist[1][j]);
         eiforce[i*dofpern+1] += (emass[i*dim+1][j*dim]   * hist[0][j]
                                + emass[i*dim+1][j*dim+1] * hist[1][j]);
         eiforce[i*dofpern+2] += (emass[iel*dim+i][j*dim]   * hist[0][j]
                                + emass[iel*dim+i][j*dim+1] * hist[1][j]);
      }
   }
break;
case 2: /* element has free surface contribution */
   row[0] = 0;
   for (i=0; i<iel; i++) row[i+1] = row[i]+ele->node[i]->numdf;
   /*---------------------------------------- perform multiplication ---*/
   for (i=0; i<iel; i++)
   {
      for (j=0; j<iel; j++)
      {
         eiforce[row[i]]   += (emass[i*dim][j*dim]   * hist[0][j]
                             + emass[i*dim][j*dim+1] * hist[1][j]);
         eiforce[row[i]+1] += (emass[i*dim+1][j*dim]   * hist[0][j]
                             + emass[i*dim+1][j*dim+1] * hist[1][j]);
         eiforce[row[i]+2] += (emass[iel*dim+i][j*dim]   * hist[0][j]
                             + emass[iel*dim+i][j*dim+1] * hist[1][j]);
      }
   }
break;
default:
   dserror("parameter fs_on out of range!\n");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_massrhs */



#endif
/*! @} (documentation module close)*/
