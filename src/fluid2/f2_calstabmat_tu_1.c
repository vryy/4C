/*!----------------------------------------------------------------------
\file
\brief stabilisation part of element stiffness matrix for fluid2

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
#include "fluid2.h"
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

static FLUID_DYNAMIC *fdyn;
/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Kkapome

<pre>                                                         he  02/03

In this routine the stabilisation part of matrix Kvv is calculated:

    /
   |  tau_tu * u * grad(kapome) * grad(psi) * u   d_omega + D.C.
  /

    /
   |  -tau_tu * div[(nue+nue_t*sig)*grad(kapome)] * grad(psi) * u   d_omega + D.C.
  /

    /
   |  tau_tu * factor * kapome_old * kapome * grad(psi) * u    d_omega + D.C.
  /



NOTE: there's only one elestif
      --> Kkapome is stored in estif[0..(iel-1)][0..(iel-1)]

</pre>
\param  *ele	     ELEMENT	   (i)	   actual element
\param  *elev	     ELEMENT	   (i)	   actual element for vel.
\param **estif          DOUBLE	   (i/o)   ele stiffness matrix
\param   kapomeint      DOUBLE	   (i)     kapome at integr. point
\param  *velint         DOUBLE	   (i)     vel. at integr. point
\param  *velint_dc      DOUBLE	   (i)     vel. at integr. point for D.C.
\param   eddyint        DOUBLE	   (i)     eddy at integr. point
\param  *funct          DOUBLE	   (i)     nat. shape functions
\param **derxy          DOUBLE	   (i)     global derivatives
\param **derxy2         DOUBLE	   (i)     2nd global derivatives
\param   fac	      DOUBLE	   (i)     weighting factor
\param   visc	      DOUBLE	   (i)     fluid viscosity
\param   factor	      DOUBLE	   (i)     factor
\param   sig	      DOUBLE	   (i)     factor
\param   iel	      INT		   (i)     num. of nodes in ele
\return void

------------------------------------------------------------------------*/
void f2_calstabkkapome(
                ELEMENT         *ele,
		    ELEMENT         *elev,
		    DOUBLE         **estif,
		    DOUBLE           kapomeint,
		    DOUBLE          *velint,
		    DOUBLE          *velint_dc,
                DOUBLE           eddyint,
                DOUBLE          *funct,
		    DOUBLE         **derxy,
		    DOUBLE         **derxy2,
		    DOUBLE           fac,
		    DOUBLE           visc,
		    DOUBLE           factor,
		    DOUBLE           sig,
                INT              iel
                   )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |
/-----------------------------------------------------------------------*/
INT    irow,icol,irn,icn,ird;
DOUBLE taumu,taumu_dc;
DOUBLE c,c_dc;
DOUBLE auxc,aux;

#ifdef DEBUG
dstrc_enter("f2_calstabkkapome");
#endif

/*---------------------------------------- set stabilisation parameter */
fdyn      = alldyn[genprob.numff].fdyn;

taumu     = fdyn->tau_tu;
taumu_dc  = fdyn->tau_tu_dc;
c    = fac*taumu;
c_dc = fac*taumu_dc;
/*----------------------------------------------------------------------*
    /
   | tau_tu * u * grad(kapome) * grad(psi) * u  d_omega
  /
 *----------------------------------------------------------------------*/
 icol=0;
   for(icn=0;icn<iel;icn++)
   {
     irow=0;
     auxc = velint[0]*derxy[0][icn]+velint[1]*derxy[1][icn];

      for(irn=0;irn<iel;irn++)
      {
         estif[irow][icol] += auxc*c   *(velint[0]   *derxy[0][irn]+velint[1]   *derxy[1][irn]);
         estif[irow][icol] += auxc*c_dc*(velint_dc[0]*derxy[0][irn]+velint_dc[1]*derxy[1][irn]);
         irow += 1;
      } /* end of loop over irn */
      icol += 1;
   } /* end of loop over icn */

/*----------------------------------------------------------------------*
    /
   |  -tau_tu * div[(nue+nue_t*sig)*grad(kapome)] * grad(psi) * u   d_omega =
  /

    /
   |  -tau_tu *  (nue+nue_t*sig) *  div grad(kapome) * grad(psi) * u   d_omega
  /

 *----------------------------------------------------------------------*/
   icol=0;
      for (icn=0;icn<iel;icn++)
      {
 	 irow = 0;
       auxc = (visc+eddyint*sig)*(derxy2[0][icn]+derxy2[1][icn]);

	 for (irn=0;irn<iel;irn++)
	 {
          estif[irow][icol] -= auxc*c   *(derxy[0][irn]*velint[0]   +derxy[1][irn]*velint[1]);
          estif[irow][icol] -= auxc*c_dc*(derxy[0][irn]*velint_dc[0]+derxy[1][irn]*velint_dc[1]);
	    irow += 1;
	 } /* end of loop over irn */
	 icol += 1;
      } /* end of loop over icn */

/*----------------------------------------------------------------------*
    /
   |  tau_tu * factor * kapome_old * kapome * grad(psi) * u   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
	 irow = 0;
       auxc = factor*kapomeint*funct[icn];

	 for (irn=0;irn<iel;irn++)
	 {
          estif[irow][icol]     += auxc*c   *(velint[0]   *derxy[0][irn] + velint[1]   *derxy[1][irn]);
          estif[irow][icol]     += auxc*c_dc*(velint_dc[0]*derxy[0][irn] + velint_dc[1]*derxy[1][irn]);
          irow += 1;
	 } /* end of loop over irn */
	 icol += 1;
      } /* end of loop over icn */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calstabkvv */


/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Mkapome

<pre>                                                         he  02/03

In this routine the stabilisation part of matrix Mkapome is calculated:

    /
   |   tau_tu * grad(psi) * u * kapome  d_omega + D.C.
  /


NOTE: there's only one elestif
      --> Mkapeps is stored in emass[0..(iel-1)][0..(iel-1)]

</pre>
\param  *ele	 ELEMENT	     (i)	   actual element
\param **emass     DOUBLE	   (i/o)   ele mass matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param  *velint_dc DOUBLE	   (i)     vel. at integr. point for D.C.
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   iel	   INT  	   (i)	   num. of nodes in ele
\return void

------------------------------------------------------------------------*/
void f2_calstabmkapome(
                    ELEMENT         *ele,
		        DOUBLE         **emass,
    		        DOUBLE          *velint,
    		        DOUBLE          *velint_dc,
                    DOUBLE          *funct,
		        DOUBLE         **derxy,
		        DOUBLE           fac,
		        INT              iel
                     )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |
/-----------------------------------------------------------------------*/
INT    irow,icol,irn,icn;
DOUBLE taumu,taumu_dc;
DOUBLE c,c_dc;
DOUBLE auxc;

#ifdef DEBUG
dstrc_enter("f2_calstabmkapome");
#endif

/*---------------------------------------- set stabilisation parameter */
fdyn       = alldyn[genprob.numff].fdyn;

taumu      = fdyn->tau_tu;
taumu_dc   = fdyn->tau_tu_dc;

c    = fac * taumu;
c_dc = fac * taumu_dc;

/*----------------------------------------------------------------------*
   Calculate convection stabilisation part:
    /
   |   tau_tu * grad(psi) * u * kapome  d_omega
  /
 *----------------------------------------------------------------------*/
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      irow = 0;
      auxc = funct[icn];

      for (irn=0;irn<iel;irn++)
      {
	 emass[irow][icol] += auxc*c   *(velint[0]   *derxy[0][irn]+velint[1]   *derxy[1][irn]);
 	 emass[irow][icol] += auxc*c_dc*(velint_dc[0]*derxy[0][irn]+velint_dc[1]*derxy[1][irn]);
	 irow += 1;
      } /* end loop over irn */
      icol += 1;
   } /* end loop over icn */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calstabmvv */

#endif
