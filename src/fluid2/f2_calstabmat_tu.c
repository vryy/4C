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
\brief evaluate stabilisaton part of Kkapeps

<pre>                                                         he  12/02

In this routine the stabilisation part of matrix Kvv is calculated:

    /
   |  tau_tu * u * grad(kapeps) * grad(psi) * u   d_omega    +  D. C.
  /

    /
   |  -tau_tu * div[(nue+nue_t/sig)*grad(kapeps)] * grad(psi) * u   d_omega  +  D. C.
  /

    /
   |  tau_tu * factor * kapeps_old * kapeps * grad(psi) * u    d_omega   +  D. C.
  /

LOW-REYNOLD's MODEL only for kappa:

    /
   |  tau_tu *2.0*visc*[ 2*grad(k_old)/(4*k_old)  * grad (k) - grad(k_old)*grad(k_old)/(4*k_old^2) * k ] *  grad(psi) * ud_omega  +  D. C.
  /


NOTE: there's only one elestif
      --> Kkapeps is stored in estif[0..(iel-1)][0..(iel-1)]

</pre>
\param  *ele	     ELEMENT	   (i)	   actual element
\param  *elev	     ELEMENT	   (i)	   actual element for vel.
\param **estif          DOUBLE	   (i/o)   ele stiffness matrix
\param   kapepsint      DOUBLE	   (i)     kapeps at integr. point
\param  *velint         DOUBLE	   (i)     vel. at integr. point
\param  *velint_dc      DOUBLE	   (i)     vel. at integr. point for DISC. CAPT.
\param   eddyint        DOUBLE	   (i)     eddy at integr. point
\param  *kapepsderxy    DOUBLE	   (i)     kapeps deriv. at integr. point
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
void f2_calstabkkapeps(
                ELEMENT         *ele,
		    ELEMENT         *elev,
		    DOUBLE         **estif,
		    DOUBLE           kapepsint,
		    DOUBLE          *velint,
		    DOUBLE          *velint_dc,
                DOUBLE           eddyint,
                DOUBLE          *kapepsderxy,
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
dstrc_enter("f2_calstabkkapeps");
#endif

/*---------------------------------------- set stabilisation parameter */
fdyn      = alldyn[genprob.numff].fdyn;

taumu     = fdyn->tau_tu;
taumu_dc  = fdyn->tau_tu_dc;
c    = fac*taumu;
c_dc = fac*taumu_dc;
/*----------------------------------------------------------------------*
    /
   | tau_tu * u * grad(kapeps) * grad(psi) * u  d_omega
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
   |  -tau_tu * div[(nue+nue_t/sig)*grad(kapeps)] * grad(psi) * u   d_omega =
  /

    /
   |  -tau_tu *  (nue+nue_t/sig) *  div grad(kapeps) * grad(psi) * u   d_omega
  /

 *----------------------------------------------------------------------*/
   icol=0;
      for (icn=0;icn<iel;icn++)
      {
 	 irow = 0;
       auxc = (visc+eddyint/sig)*(derxy2[0][icn]+derxy2[1][icn]);

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
   |  tau_tu * factor * kapeps_old * kapeps * grad(psi) * u   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
	 irow = 0;
       auxc = factor*kapepsint*funct[icn];

	 for (irn=0;irn<iel;irn++)
	 {
          estif[irow][icol]     += auxc*c   *(velint[0]   *derxy[0][irn] + velint[1]   *derxy[1][irn]);
          estif[irow][icol]     += auxc*c_dc*(velint_dc[0]*derxy[0][irn] + velint_dc[1]*derxy[1][irn]);
          irow += 1;
	 } /* end of loop over irn */
	 icol += 1;
      } /* end of loop over icn */


if(fdyn->kapeps_flag==0)
{
/*----------------------------------------------------------------------*
LOW-REYNOLD's MODEL:

    /
   |  tau_tu * 2*grad(k_old)/(4*k_old) * 2.0 * visc * grad k* grad(psi) * u  d_omega
  /
*----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
	 irow = 0;

       auxc  =  2/(4*kapepsint)*
               (kapepsderxy[0]*derxy[0][icn]+kapepsderxy[1]*derxy[1][icn]);
       auxc -= (pow(kapepsderxy[0],2)+pow(kapepsderxy[1],2))/
               (4*kapepsint*kapepsint) * funct[icn];

	 for (irn=0;irn<iel;irn++)
	 {
          estif[irow][icol]   += 2.0*visc*c   *auxc*(velint[0]   *derxy[0][irn] + velint[1]   *derxy[1][irn]);
          estif[irow][icol]   += 2.0*visc*c_dc*auxc*(velint_dc[0]*derxy[0][irn] + velint_dc[1]*derxy[1][irn]);
          irow += 1;
	 } /* end of loop over irn */
	 icol += 1;
      } /* end of loop over icn */

} /* endif */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calstabkvv */


/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Mkapeps

<pre>                                                         he  12/02

In this routine the stabilisation part of matrix Mkapeps is calculated:

    /
   |   tau_tu * grad(psi) * u * kapeps  d_omega + D. C.
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
void f2_calstabmkapeps(
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
dstrc_enter("f2_calstabmkapeps");
#endif

/*---------------------------------------- set stabilisation parameter */
fdyn      = alldyn[genprob.numff].fdyn;

taumu     = fdyn->tau_tu;
taumu_dc  = fdyn->tau_tu_dc;

c    = fac * taumu;
c_dc = fac * taumu_dc;

/*----------------------------------------------------------------------*
   Calculate convection stabilisation part:
    /
   |   tau_tu * grad(psi) * u * kapeps  d_omega
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
