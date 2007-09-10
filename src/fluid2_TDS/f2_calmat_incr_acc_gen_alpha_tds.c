/*!----------------------------------------------------------------------
\file
\brief evaluate 2D fluid coefficient matrix for the generalised alpha
       time integration.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/gammi/
            +49-(0)89-289-15235
</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
/*!
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2.h"
#ifdef D_FLUID2_TDS
#include "fluid2_TDS_prototypes.h"

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
\brief evaluate fluid coefficient matrix

<pre>                                                        gammi 11/06

In this routine the Gauss point contributions to the Galerkin part of
the elemental coefficient matrix of a stabilised fluid2 element are
calculated.

The procedure is based on the linearisation discribed in dis. Whiting.

integration schemes:

Gen-Alpha:


for further comments see comment lines within code.


</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *velint     DOUBLE        (i)   last intermediate vel at INT point
\param  *funct      DOUBLE        (i)   nat. shape funcs
\param **derxy      DOUBLE        (i)   global coord. deriv.
\param **derxy2     DOUBLE        (i)   2nd global coord. deriv.
\param   fac 	    DOUBLE        (i)   weighting factor
\param   visc       DOUBLE        (i)   fluid viscosity
\param   iel	    INT           (i)   number of nodes of act. ele
                                        derivative


\return void

------------------------------------------------------------------------*/
void f2_calgalmat_gen_alpha_tds(
                DOUBLE **estif,
		DOUBLE  *velint,
		DOUBLE  *funct,
		DOUBLE **derxy,
		DOUBLE **derxy2,
		DOUBLE   fac,
		DOUBLE   visc,
		int      iel
              )
{
 int       i;        /* default counter                             */
 int       ri,ci;    /* counters for row and column index           */


 double    alpha_M;  /* generalised alpha parameter                 */
                     /* accelerations -> intermediate accelerations */
 double    alpha_F;  /* generalised alpha parameter                 */
                     /* velocities -> intermediate velocities       */
 double    theta;    /* generalised alpha parameter                 */
                     /* accelerations -> velocities                 */
 double    dt;       /* timestepsize                                */

 double    aftdt;    /* alpha_F * theta * dt                        */

 FLUID_DYNAMIC   *fdyn;

 /* viscous term partially integrated */
 DOUBLE  viscous[2][2][2*MAXNOD];

 double    aux;

#ifdef DEBUG
 dstrc_enter("f2_calgalmat_gen_alpha_tds");
#endif


 /*----------------------- get time integration control information */
 fdyn    = alldyn[genprob.numff].fdyn;

 /*----------------------------- set constants for time integration */
 dt      = fdyn->dt;
 theta   = fdyn->theta;
 alpha_M = fdyn->alpha_m;
 alpha_F = fdyn->alpha_f;

 aftdt   = alpha_F*theta*dt;



 for (i=0; i<iel; i++) /* loop over nodes of element */
 {
     /*--- viscous term (after integr. by parts) -------------------------*/
     /*   /				\
	  1 |  2 N_x,x    N_x,y + N_y,x  |    with N_x .. x-line of N
	  - |                            |         N_y .. y-line of N
	  2 |  N_y,x + N_x,y    2 N_y,y  |
	  \                             /                                 */
     viscous[0][0][2*i]   = derxy[0][i];
     viscous[0][0][2*i+1] = 0.0;                /* 1st index:             */
     viscous[0][1][2*i]   = 0.5 * derxy[1][i];  /*   line of epsilon      */
     viscous[0][1][2*i+1] = 0.5 * derxy[0][i];  /* 2nd index:             */
     viscous[1][0][2*i]   = 0.5 * derxy[1][i];  /*   column of epsilon    */
     viscous[1][0][2*i+1] = 0.5 * derxy[0][i];  /* 3rd index:             */
     viscous[1][1][2*i+1] = derxy[1][i];        /*   elemental vel dof    */
     viscous[1][1][2*i]   = 0.0;
 }

 for (ri=0; ri<iel; ri++)       /* row index    */
 {
  for (ci=0; ci<iel; ci++)   /* column index */
  {
   /*--------------------------------------------------------------*/
   /*                       GALERKIN MVV                           */
   /*--------------------------------------------------------------*/

   /* 'mass matrix' alpha_M * (Dacc,v) */
   aux = alpha_M * fac ;

   estif[ri*3  ][ci*3  ] += funct[ri] * funct[ci] * aux;
   estif[ri*3+1][ci*3+1] += funct[ri] * funct[ci] * aux;

   /*------------------------------------------------------------------*/
   /*                         GALERKIN KVV                             */
   /*------------------------------------------------------------------*/

   /*  alpha_F * theta * dt * (2 * nu * epsilon(Dacc), epsilon(v)) */
   aux = 2 * visc * aftdt * fac;

   estif[ri*3  ][ci*3  ] += (viscous[0][0][ri*2  ]*viscous[0][0][ci*2  ]
			    +viscous[0][1][ri*2  ]*viscous[1][0][ci*2  ]
			    +viscous[1][0][ri*2  ]*viscous[0][1][ci*2  ]
			    +viscous[1][1][ri*2  ]*viscous[1][1][ci*2  ]
                            ) * aux;
   estif[ri*3  ][ci*3+1] += (viscous[0][0][ri*2  ]*viscous[0][0][ci*2+1]
			    +viscous[0][1][ri*2  ]*viscous[1][0][ci*2+1]
			    +viscous[1][0][ri*2  ]*viscous[0][1][ci*2+1]
			    +viscous[1][1][ri*2  ]*viscous[1][1][ci*2+1]
                            ) * aux;
   estif[ri*3+1][ci*3  ] += (viscous[0][0][ri*2+1]*viscous[0][0][ci*2  ]
			    +viscous[0][1][ri*2+1]*viscous[1][0][ci*2  ]
			    +viscous[1][0][ri*2+1]*viscous[0][1][ci*2  ]
			    +viscous[1][1][ri*2+1]*viscous[1][1][ci*2  ]
                            ) * aux;
   estif[ri*3+1][ci*3+1] += (viscous[0][0][ri*2+1]*viscous[0][0][ci*2+1]
			    +viscous[0][1][ri*2+1]*viscous[1][0][ci*2+1]
			    +viscous[1][0][ri*2+1]*viscous[0][1][ci*2+1]
			    +viscous[1][1][ri*2+1]*viscous[1][1][ci*2+1]
                            ) * aux;

   /*-------------------------------------------------------------------*/
   /*                         GALERKIN KPV                              */
   /*-------------------------------------------------------------------*/

   /*  alpha_F * theta * dt * (div Dacc, q) */
   aux = aftdt * fac;


   estif[ri*3+2][ci*3  ] += funct[ri] * derxy[0][ci] * aux;
   estif[ri*3+2][ci*3+1] += funct[ri] * derxy[1][ci] * aux;


   /*-------------------------------------------------------------------*/
   /*                         GALERKIN KVP                              */
   /*-------------------------------------------------------------------*/

   /* - ( Dp , div v ) */
   aux = fac;

   estif[ri*3  ][ci*3+2] -= derxy[0][ri] * funct[ci] * aux;
   estif[ri*3+1][ci*3+2] -= derxy[1][ri] * funct[ci] * aux;
  } /* end loop over column index ci */
 } /* end loop over row index ri */


 /*----------------------------------------------------------------------*/
#ifdef DEBUG
 dstrc_exit();
#endif

 return;
}


/*!---------------------------------------------------------------------
\brief evaluate fluid stabilisation matrix

<pre>                                                        gammi 11/06

In this routine the Gauss point contributions to the stabilisation part
of the elemental coefficient matrix of a stabilised fluid2 element are
calculated.

The procedure is based on the linearisation discribed in dis. Whiting.

integration schemes:

Gen-Alpha:


for further comments see comment lines within code.


</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *velint     DOUBLE        (i)   last intermediate vel at INT point
\param  *funct      DOUBLE        (i)   nat. shape funcs
\param **derxy      DOUBLE        (i)   global coord. deriv.
\param **derxy2     DOUBLE        (i)   2nd global coord. deriv.
\param   fac 	    DOUBLE        (i)   weighting factor
\param   visc       DOUBLE        (i)   fluid viscosity
\param   iel	    INT           (i)   number of nodes of act. ele
                                        derivative


\return void

------------------------------------------------------------------------*/
void f2_calstabmat_gen_alpha_tds(
                DOUBLE **estif,
		DOUBLE  *velint,
		DOUBLE  *funct,
		DOUBLE **derxy,
		DOUBLE **derxy2,
		DOUBLE   fac,
		DOUBLE   visc,
		int      iel
              )
{
 int       i;        /* default counter                             */
 int       ri,ci;    /* counters for row and column index           */


 double    alpha_M;  /* generalised alpha parameter                 */
                     /* accelerations -> intermediate accelerations */
 double    alpha_F;  /* generalised alpha parameter                 */
                     /* velocities -> intermediate velocities       */
 double    theta;    /* generalised alpha parameter                 */
                     /* accelerations -> velocities                 */
 double    dt;       /* timestepsize                                */

 double    aftdt;    /* alpha_F * theta * dt                        */


 double    tau_C;    /* stabilisation parameter --- continuity      */
 double    tau_M;    /* stabilisation parameter --- momentum        */

 DOUBLE  viscs2 [2][2*MAXNOD]; /* viscous term incluiding 2nd derivatives */
 DOUBLE  div[2*MAXNOD];             /* divergence of u or v              */

 FLUID_DYNAMIC   *fdyn;

 double    aux;




#ifdef DEBUG
 dstrc_enter("f2_calstabmat_gen_alpha_tds");
#endif


 /*----------------------- get time integration control information */
 fdyn    = alldyn[genprob.numff].fdyn;

 /*----------------------------- set constants for time integration */
 dt      = fdyn->dt;
 theta   = fdyn->theta;
 alpha_M = fdyn->alpha_m;
 alpha_F = fdyn->alpha_f;

 aftdt   = alpha_F*theta*dt;

 tau_M   = fdyn->tau[0];
 tau_C   = fdyn->tau[2];

 for (i=0; i<iel; i++) /* loop over nodes of element */
 {
   /*--- viscous term  - grad * epsilon(u): ----------------------------*/
   /*   /                              \
      1 |  2 N_x,xx + N_x,yy + N_y,xy  |    with N_x .. x-line of N
    - - |                              |         N_y .. y-line of N
      2 |  N_y,xx + N_x,yx + 2 N_y,yy  |
        \                             /                                 */
   viscs2[0][2*i]   = - 0.5 * ( 2.0 * derxy2[0][i] + derxy2[1][i] );
   viscs2[0][2*i+1] = - 0.5 * ( derxy2[2][i] );
   viscs2[1][2*i]   = - 0.5 * ( derxy2[2][i] );
   viscs2[1][2*i+1] = - 0.5 * ( derxy2[0][i] + 2.0 * derxy2[1][i] );

   /*--- divergence u term ---------------------------------------------*/
   div[2*i]   = derxy[0][i];
   div[2*i+1] = derxy[1][i];
 }

 for (ri=0; ri<iel; ri++)       /* row index    */
 {
  for (ci=0; ci<iel; ci++)   /* column index */
  {
   /*---------------------------------------------------------------*/
   /*                     STABILISATION KPP                         */
   /*---------------------------------------------------------------*/


      /* term  : (u_old * grad u, grad q) */



      /* factor: tau_M /(alpha_M * tau_M + dt * alpha_F * theta)
       *          * alpha_M                                         */
      /* term  : (Dacc, grad q)                                     */
      aux = tau_M/(alpha_M * tau_M + aftdt) *alpha_M * aftdt * fac;

      estif[ri*3+2][ci*3  ] += derxy[0][ri] * funct[ci] * aux;
      estif[ri*3+2][ci*3+1] += derxy[1][ri] * funct[ci] * aux;

      /* (div epsilon(Dacc), grad q) */
      /* viscs already contains - sign!!                                   */
      aux = 2 * visc * tau_M/(alpha_M * tau_M + aftdt) * aftdt * aftdt *fac;

      estif[ri*3+2][ci*3]   += (derxy[0][ri] * viscs2[0][2*ci]
                               +derxy[1][ri] * viscs2[1][2*ci]) * aux;
      estif[ri*3+2][ci*3+1] += (derxy[0][ri] * viscs2[0][2*ci+1]
                               +derxy[1][ri] * viscs2[1][2*ci+1]) * aux;


      /* factor: tau_M /(alpha_M * tau_M + dt * alpha_F * theta)
         term  : (grad Dp, grad q)                                   */
      aux = tau_M/(alpha_M * tau_M + aftdt) * aftdt * fac;

      estif[ri*3+2][ci*3+2] += (derxy[0][ri] * derxy[0][ci]
			       +derxy[1][ri] * derxy[1][ci]) * aux;


      /* factor:                                                    */
      /* term  : (div Dacc, div v)                                  */
      aux = theta*dt*(tau_C/(alpha_M*tau_C+aftdt))*aftdt*fac;

      estif[ri*3  ][ci*3  ] += div[ri*2  ] * div[ci*2  ] * aux;
      estif[ri*3  ][ci*3+1] += div[ri*2  ] * div[ci*2+1] * aux;
      estif[ri*3+1][ci*3  ] += div[ri*2+1] * div[ci*2  ] * aux;
      estif[ri*3+1][ci*3+1] += div[ri*2+1] * div[ci*2+1] * aux;

      /* factor:                                                    */
      /* term  : (grad Dp, div eps(v))                                  */
      /* viscs already contains - sign!!                              */
      aux = tau_M/(alpha_M * tau_M + aftdt) * 2 * visc *fac * aftdt;

      estif[ri*3  ][ci*3+2] -= (viscs2[0][2*ri  ] * derxy[0][ci]
                               +viscs2[1][2*ri  ] * derxy[1][ci]) * aux;
      estif[ri*3+1][ci*3+2] -= (viscs2[0][2*ri+1] * derxy[0][ci]
                               +viscs2[1][2*ri+1] * derxy[1][ci]) * aux;



      /* factor:                                                     */
      /* term  : (Dacc, div eps(v))                                  */
      /* viscs already contains - sign!!                             */
      aux = funct[ci] * alpha_M * 2 * visc * aftdt *
	  tau_M/(alpha_M * tau_M + aftdt) * fac;

      estif[ri*3  ][ci*3  ] -= viscs2[0][2*ri  ] * aux;
      estif[ri*3  ][ci*3+1] -= viscs2[1][2*ri  ] * aux;
      estif[ri*3+1][ci*3  ] -= viscs2[0][2*ri+1] * aux;
      estif[ri*3+1][ci*3+1] -= viscs2[1][2*ri+1] * aux;


      /* factor:                                                     */
      /* term  : (div epsilon(Dacc), div epsilon(v))                 */
      /* viscs already contains - sign!!                             */
      aux = tau_M/(alpha_M * tau_M + aftdt) *
	  4 * visc * visc * aftdt * aftdt * fac;

      estif[ri*3  ][ci*3  ] -= (viscs2[0][2*ri  ] * viscs2[0][2*ci  ]
                               +viscs2[1][2*ri  ] * viscs2[1][2*ci  ]) * aux;
      estif[ri*3+1][ci*3  ] -= (viscs2[0][2*ri+1] * viscs2[0][2*ci  ]
                               +viscs2[1][2*ri+1] * viscs2[1][2*ci  ]) * aux;
      estif[ri*3  ][ci*3+1] -= (viscs2[0][2*ri  ] * viscs2[0][2*ci+1]
                               +viscs2[1][2*ri  ] * viscs2[1][2*ci+1]) * aux;
      estif[ri*3+1][ci*3+1] -= (viscs2[0][2*ri+1] * viscs2[0][2*ci+1]
                               +viscs2[1][2*ri+1] * viscs2[1][2*ci+1]) * aux;


      /* factor:                                                    */
      /* term  : (grad Dp, v)                                       */
      aux = alpha_M * tau_M/(alpha_M * tau_M + aftdt) * fac;

      estif[ri*3  ][ci*3+2] -= funct[ri] * derxy[0][ci] * aux;
      estif[ri*3+1][ci*3+2] -= funct[ri] * derxy[1][ci] * aux;



      /* factor:                                                    */
      /* term (Dacc,v) */
      aux = alpha_M * alpha_M * tau_M/(alpha_M * tau_M + aftdt) * fac;

      estif[ri*3  ][ci*3  ] -= funct[ri] * funct[ci] * aux;
      estif[ri*3+1][ci*3+1] -= funct[ri] * funct[ci] * aux;


      /* factor:                                                     */
      /* term  : (div eps(Dacc), v)                                  */
      /* viscs already contains - sign!!                             */
      aux = 2 * visc * aftdt * alpha_M * tau_M/(alpha_M * tau_M + aftdt) * fac;

      estif[ri*3  ][ci*3  ]   -= funct[ri] * viscs2[0][2*ci  ] * aux;
      estif[ri*3  ][ci*3+1]   -= funct[ri] * viscs2[1][2*ci  ] * aux;
      estif[ri*3+1][ci*3  ]   -= funct[ri] * viscs2[0][2*ci+1] * aux;
      estif[ri*3+1][ci*3+1]   -= funct[ri] * viscs2[1][2*ci+1] * aux;
  } /* end loop over column index ci */
 } /* end loop over row index ri */


 /*-----------------------------------------------------------------*/
#ifdef DEBUG
 dstrc_exit();
#endif

 return;
}





#endif /*D_FLUID2_TDS*/
#endif /*D_FLUID2*/
/*! @} (documentation module close)*/
#endif
