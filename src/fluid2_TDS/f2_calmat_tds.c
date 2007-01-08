/*!----------------------------------------------------------------------
\file
\brief evaluate 2D fluid coefficient matrix for one step theta

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/gammi/
            +49-(0)89-289-15235
</pre>

------------------------------------------------------------------------*/
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

static FLUID_DYNAMIC *fdyn;

/*!---------------------------------------------------------------------
\brief evaluate fluid coefficient matrix

<pre>                                                        gammi 11/06

In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilised fluid2 element are calculated. The procedure is
based on the Rothe method of first integrating in time. Hence the
resulting terms include coefficients containing time integration variables
such as theta or delta t which are represented by 'timefac'.

The stabilisation is based on the time evolution of the residuum.

integration schemes:

One-step-Theta:
rhsint = u_old + Theta dt f + (1-Theta) acc_old +
         +u_subscale_old + (1-Theta) acc_subscale_old 

NOTE: Galerkin and stabilisation matrices are calculated within one
      routine.


Notational remarks:

                   /              \
                  | u_x,x   u_x,y |
vderxy = grad u = |               |
                  | u_y,x   u_y,y |
                  \               /

           /                         \
          | u_x,xx   u_x,yy   u_x,xy |
vderxy2 = |                          |
          | u_y,xx   u_y,yy   u_y,xy |
          \                          /

for further comments see comment lines within code.


</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\param  *velint     DOUBLE        (i)   vel at INT point
\param  *histvec    DOUBLE        (i)   rhs at INT point
\param  *gridvint   DOUBLE        (i)   gridvel at INT point
\param **vderxy     DOUBLE        (i)   global vel derivatives
\param  *vderxy2    DOUBLE        (i)   2nd global vel derivatives
\param  *funct      DOUBLE        (i)   nat. shape funcs
\param **derxy      DOUBLE        (i)   global coord. deriv.
\param **derxy2     DOUBLE        (i)   2nd global coord. deriv.
\param  *edeadng    DOUBLE        (i)   dead load at time n+1
\param   fac 	    DOUBLE        (i)   weighting factor
\param   visc       DOUBLE        (i)   fluid viscosity
\param   iel	    INT           (i)   number of nodes of act. ele
\param  *hasext     INT           (i)   flag, if element has volume load
\param   isale      INT           (i)   flag, if ALE or EULER
\param   is_relax   INT           (i)
\param	 sub_pres   DOUBLE        (i)   old subscale pressure at GP
\param	 divu_old   DOUBLE        (i)   old divergence at GP
\param	 sub_vel[2] DOUBLE        (i)   old subscale velocities at GP
\param   sub_vel_trial_wo_facMtau[2]
                    DOUBLE        (i)   estimdated subscale velocities at GP
                                        (still have to be multiplied by
                                         facMtau)
\param	 old_vel[2] DOUBLE        (i)   old velocities at GP
\param	 res_old[2] DOUBLE        (i)   part of old residual without time
                                        derivative


\return void

------------------------------------------------------------------------*/
void f2_calmat_tds(
                DOUBLE **estif,
		DOUBLE  *eforce,
		DOUBLE  *velint,
		DOUBLE   histvec[2],
		DOUBLE   gridvint[2],
		DOUBLE   press,
		DOUBLE **vderxy,
                DOUBLE **vderxy2,
                DOUBLE   gradp[2],
		DOUBLE  *funct,
		DOUBLE **derxy,
		DOUBLE **derxy2,
                DOUBLE  *edeadng,
		DOUBLE   fac,
		DOUBLE   visc,
		INT      iel,
                INT     *hasext,
                INT      isale,
                INT      is_relax,
		DOUBLE   sub_pres,
		DOUBLE   divu_old,
		DOUBLE   sub_vel[2],
		DOUBLE   sub_vel_trial_wo_facMtau[2],
		DOUBLE   old_vel[2],
		DOUBLE   old_acc[2],
		DOUBLE   res_old[2]
              )
{
INT     i, j, ri, ci;
DOUBLE  timefac;                  /* One-step-Theta: timefac = theta*dt */
DOUBLE  dt;                                           /* time step size */
DOUBLE  aux;
DOUBLE  auxmat[2][2];
DOUBLE  tau_M, tau_C;                        /* stabilisation parameter */
DOUBLE  tau_Mp;                              /* stabilisation parameter */
DOUBLE  facC,facCtau;
DOUBLE  facM,facMtau;
DOUBLE  old_facC,old_facCtau;
DOUBLE  old_facM,old_facMtau;
DOUBLE  viscs2[2][2*MAXNOD]; /* viscous term incluiding 2nd derivatives */
DOUBLE  viscous[2][2][2*MAXNOD];   /* viscous term partially integrated */
DOUBLE  conv_c[MAXNOD];    /* linearisation of convect, convective part */
DOUBLE  conv_g[MAXNOD];          /* linearisation of convect, grid part */
DOUBLE  conv_r[2][2*MAXNOD]; /* linearisation of convect, reactive part */
DOUBLE  div[2*MAXNOD];                          /* divergence of u or v */
DOUBLE  ugradv[MAXNOD][2*MAXNOD];        /* linearisation of u * grad v */
DOUBLE  conv_old[2]; /* convective term evalaluated with old velocities */
DOUBLE  visc_old[2]; /* viscous term evaluated with old velocities      */
DOUBLE  rhsint[2];   /* total right hand side terms at int.-point       */

DOUBLE  time2nue, timetauM, timetauMp, ttimetauM, ttimetauMp, timefacfac;
DOUBLE  theta;

INT     cross_stress   =1;
INT     reynolds_stress=1;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("f2_calmat_tds");
#endif
/*========================== initialisation ============================*/
fdyn = alldyn[genprob.numff].fdyn;

timefac = fdyn->thsl;
theta   = fdyn->theta;
dt      = fdyn->dt;

tau_M  = fdyn->tau[0]*fac;
tau_Mp = fdyn->tau[0]*fac;
tau_C  = fdyn->tau[2]*fac;

facC       = 1./(fdyn->tau[2]+theta*dt);
facCtau    = fdyn->tau[2]*facC;


facM       = 1./(fdyn->tau[0]+theta*dt);
facMtau    = fdyn->tau[0]*facM;

old_facC   = 1./(fdyn->tau_old[2]+theta*dt);
old_facCtau= fdyn->tau_old[2]*facC;


old_facM   = 1./(fdyn->tau_old[0]+theta*dt);
old_facMtau= fdyn->tau_old[0]*facM;


/* integration factors and koefficients of single terms */
time2nue  = timefac * 2.0 * visc;
timetauM   = timefac * tau_M;
timetauMp  = timefac * tau_Mp;

ttimetauM  = timefac * timetauM;
ttimetauMp = timefac * timetauMp;
timefacfac = timefac * fac;


if(isale)
{
    dserror("ALE for time dependent subscales not implemented yet!");
}

/*------------------------- evaluate rhs vector at integration point ---*/
if (!is_relax)
{
if (*hasext)
{
   rhsint[0] = timefac * edeadng[0] + histvec[0];
   rhsint[1] = timefac * edeadng[1] + histvec[1];
}
else
{
   rhsint[0] = histvec[0];
   rhsint[1] = histvec[1];
}
}
else
{
  rhsint[0]=0;
  rhsint[1]=0;
}

/*----------------- get numerical representation of single operators ---*/

/* Convective term  u_old * grad u_old: */
conv_old[0] = vderxy[0][0] * velint[0] + vderxy[0][1] * velint[1];
conv_old[1] = vderxy[1][0] * velint[0] + vderxy[1][1] * velint[1];

/* Viscous term  div epsilon(u_old) */
visc_old[0] = 0.5 * (2.0*vderxy2[0][0] + vderxy2[0][1] + vderxy2[1][2]);
visc_old[1] = 0.5 * (2.0*vderxy2[1][1] + vderxy2[1][0] + vderxy2[0][2]);

for (i=0; i<iel; i++) /* loop over nodes of element */
{
   /* Reactive term  u:  funct */
   /* linearise convective term */

   /*--- convective part u_old * grad (funct) --------------------------*/
   /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
   conv_c[i] = derxy[0][i] * velint[0] + derxy[1][i] * velint[1] ;

   /*--- convective grid part u_G * grad (funct) -----------------------*/
   /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
   if(isale)
   {
     conv_g[i] = - derxy[0][i] * gridvint[0] - derxy[1][i] * gridvint[1];
   }
   else
   {
     conv_g[i] = 0;
   }

   /*--- reactive part funct * grad (u_old) ----------------------------*/
   /* /                          \
      |  u_old_x,x   u_old_x,y   |
      |                          | * N   with  N .. form function matrix
      |  u_old_y,x   u_old_y,y   |
      \                         /                                       */
   conv_r[0][2*i]   = vderxy[0][0]*funct[i];
   conv_r[0][2*i+1] = vderxy[0][1]*funct[i];
   conv_r[1][2*i]   = vderxy[1][0]*funct[i];
   conv_r[1][2*i+1] = vderxy[1][1]*funct[i];

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

   /*--- viscous term (after integr. by parts) -------------------------*/
   /*   /                            \
      1 |  2 N_x,x    N_x,y + N_y,x  |    with N_x .. x-line of N
      - |                            |         N_y .. y-line of N
      2 |  N_y,x + N_x,y    2 N_y,y  |
        \                           /                                   */
   viscous[0][0][2*i]   = derxy[0][i];
   viscous[0][0][2*i+1] = 0.0;                /* 1st index:             */
   viscous[0][1][2*i]   = 0.5 * derxy[1][i];  /*   line of epsilon      */
   viscous[0][1][2*i+1] = 0.5 * derxy[0][i];  /* 2nd index:             */
   viscous[1][0][2*i]   = 0.5 * derxy[1][i];  /*   column of epsilon    */
   viscous[1][0][2*i+1] = 0.5 * derxy[0][i];  /* 3rd index:             */
   viscous[1][1][2*i+1] = derxy[1][i];        /*   elemental vel dof    */
   viscous[1][1][2*i]   = 0.0;

   /* pressure gradient term derxy, funct without or with integration   *
    * by parts, respectively                                            */

   /*--- divergence u term ---------------------------------------------*/
   div[2*i]   = derxy[0][i];
   div[2*i+1] = derxy[1][i];

   /*--- ugradv-Term ---------------------------------------------------*/
   /*
     /                                                          \
     |  N1*N1,x  N1*N1,y  N2*N1,x  N2*N1,y  N3*N1,x ...       . |
     |                                                          |
     |  N1*N2,x  N1*N2,y  N2*N2,x  N2*N2,y  N3*N2,x ...       . |
     |                                                          |
     |  N1*N3,x  N1*N3,y  N2*N3,x  N2*N3,y  N3*N3,x ...       . |
     |                                           .              |
     |  . . .                                        .          |
     |                                                  Ni*Ni,y |
     \                                                          /       */
   /* remark: vgradu = ugradv^T */
   for (j=0; j<iel; j++)
   {
      ugradv[i][2*j]   = derxy[0][i] * funct[j];
      ugradv[i][2*j+1] = derxy[1][i] * funct[j];
   }

}

/*--------------------------------- now build single stiffness terms ---*/
for (ri=0; ri<iel; ri++)      /* row index */
{
   for (ci=0; ci<iel; ci++)   /* column index */
   {
      /************** integrate element coefficient matrix **************/
/*===================== GALERKIN part of the matrix ====================*/

      /* a concentration of the following terms: */
      /* 'mass matrix' (u,v) */
      /* N_c (u_old * grad u, v) */
      /* N_r (u * grad u_old, v) */
      aux = funct[ri] * ( funct[ci]*fac + timefacfac*conv_c[ci] );
      estif[ri*3][ci*3]     += funct[ri] * conv_r[0][2*ci] * timefacfac + aux;
      estif[ri*3][ci*3+1]   += funct[ri] * conv_r[0][2*ci+1] * timefacfac;
      estif[ri*3+1][ci*3]   += funct[ri] * conv_r[1][2*ci] * timefacfac;
      estif[ri*3+1][ci*3+1] += funct[ri] * conv_r[1][2*ci+1] * timefacfac + aux;
      /* ALE: N_c (-u_G * grad u, v) */
      if(isale)
      {
         aux = timefacfac * funct[ri] * conv_g[ci];
         estif[ri*3][ci*3]     += aux;
         estif[ri*3+1][ci*3+1] += aux;
      }

      /* K (2 * nu * epsilon(u), epsilon(v)) */
      auxmat[0][0] = viscous[0][0][ri*2]   * viscous[0][0][ci*2]
                   + viscous[0][1][ri*2]   * viscous[1][0][ci*2]
                   + viscous[1][0][ri*2]   * viscous[0][1][ci*2]
                   + viscous[1][1][ri*2]   * viscous[1][1][ci*2];
      auxmat[0][1] = viscous[0][0][ri*2]   * viscous[0][0][ci*2+1]
                   + viscous[0][1][ri*2]   * viscous[1][0][ci*2+1]
                   + viscous[1][0][ri*2]   * viscous[0][1][ci*2+1]
                   + viscous[1][1][ri*2]   * viscous[1][1][ci*2+1];
      auxmat[1][0] = viscous[0][0][ri*2+1] * viscous[0][0][ci*2]
                   + viscous[0][1][ri*2+1] * viscous[1][0][ci*2]
                   + viscous[1][0][ri*2+1] * viscous[0][1][ci*2]
                   + viscous[1][1][ri*2+1] * viscous[1][1][ci*2];
      auxmat[1][1] = viscous[0][0][ri*2+1] * viscous[0][0][ci*2+1]
                   + viscous[0][1][ri*2+1] * viscous[1][0][ci*2+1]
                   + viscous[1][0][ri*2+1] * viscous[0][1][ci*2+1]
                   + viscous[1][1][ri*2+1] * viscous[1][1][ci*2+1];
      aux = time2nue * fac;
      estif[ri*3][ci*3]     += auxmat[0][0]*aux;
      estif[ri*3][ci*3+1]   += auxmat[0][1]*aux;
      estif[ri*3+1][ci*3]   += auxmat[1][0]*aux;
      estif[ri*3+1][ci*3+1] += auxmat[1][1]*aux;
      /* G (- div v, p) */
      estif[ri*3][ci*3+2]   -= timefacfac * derxy[0][ri] * funct[ci];
      estif[ri*3+1][ci*3+2] -= timefacfac * derxy[1][ri] * funct[ci];
      /* G^T ( div u, q) */
      estif[ri*3+2][ci*3]   += timefacfac * funct[ri] * derxy[0][ci];
      estif[ri*3+2][ci*3+1] += timefacfac * funct[ri] * derxy[1][ci];

/*=================== Stabilisation part of the matrix =================*/

#if 0
      /* ALE: -tau_M*timefac*timefac*(-u_G * grad u, u_old * grad v) */
      if(isale)
      {
         aux = ttimetauM * conv_c[ri] * conv_g[ci];
         estif[ri*3][ci*3]     += aux;
         estif[ri*3+1][ci*3+1] += aux;
      }

      /*--- ALE only: CONVECTIVE GRID stabilisation ---*/
      if(isale)
      {
         /* a concentration of the following terms: */
         /* -tau_M*timefac*(u, -u_G * grad v) */
         /* -tau_M*timefac*timefac*(u_old * grad u, -u_G * grad v) */
         /* -tau_M*timefac*timefac*(-u_G * grad u, -u_G * grad v) */
         aux = conv_g[ri] *
              (ttimetauM*(conv_c[ci]+conv_g[ci]) + timetauM*funct[ci]);
         estif[ri*3][ci*3]     += aux;
         estif[ri*3+1][ci*3+1] += aux;
         /* a concentration of the following two terms: */
         /* -tau_M*timefac*timefac*(u * grad u_old, -u_G * grad v) */
         /* tau_M*timefac*timefac*2*nu*(div epsilon(u), -u_G * grad v) */
         aux = timetauM * time2nue;
         estif[ri*3][ci*3]     += conv_g[ri] * ( conv_r[0][2*ci]*ttimetauM
                                                +viscs2[0][2*ci] * aux );
         estif[ri*3][ci*3+1]   += conv_g[ri] * ( conv_r[0][2*ci+1]*ttimetauM
                                                +viscs2[0][2*ci+1] * aux );
         estif[ri*3+1][ci*3]   += conv_g[ri] * ( conv_r[1][2*ci]*ttimetauM
                                                +viscs2[1][2*ci] * aux );
         estif[ri*3+1][ci*3+1] += conv_g[ri] * ( conv_r[1][2*ci+1]*ttimetauM
                                                +viscs2[1][2*ci+1] * aux );
         /* -tau_M*timefac*timefac*(grad p, -u_G * grad v) */
         estif[ri*3][ci*3+2]   += conv_g[ri] * derxy[0][ci] * ttimetauM;
         estif[ri*3+1][ci*3+2] += conv_g[ri] * derxy[1][ci] * ttimetauM;
      }

      /*ALE -tau_M*timefac*timefac*(-u_G * grad u, grad q) */
      if(isale)
      {
	  estif[ri*3+2][ci*3]   += derxy[0][ri] * conv_g[ci] * ttimetauMp;
	  estif[ri*3+2][ci*3+1] += derxy[1][ri] * conv_g[ci] * ttimetauMp;
      }

      /*ALE: tau_M*timefac*timefac*2*nu*(-u_G * grad u, div epsilon(v)) */
      if(isale)
      {
         aux = timetauMp * time2nue * conv_g[ci];
         estif[ri*3][ci*3]     += viscs2[0][2*ri] * aux;
         estif[ri*3][ci*3+1]   += viscs2[1][2*ri] * aux;
         estif[ri*3+1][ci*3]   += viscs2[0][2*ri+1] * aux;
         estif[ri*3+1][ci*3+1] += viscs2[1][2*ri+1] * aux;
      }

#endif

      /*--- TIME DEPENDENT part of stabilisation --- HIGHER ORDER TERMS */
      /*--- DIFFUSION part of stabilisation ---*/

      /* facMtau*timefac*timefac*2*nu*(grad p, div epsilon(v)) */
      aux = facMtau * timefac * timefac * 2 * visc *fac;
      estif[ri*3][ci*3+2]   -= (viscs2[0][2*ri] * derxy[0][ci]
                               +viscs2[1][2*ri] * derxy[1][ci]) * aux;
      estif[ri*3+1][ci*3+2] -= (viscs2[0][2*ri+1] * derxy[0][ci]
                               +viscs2[1][2*ri+1] * derxy[1][ci]) * aux;

      /* facMtau*timefac*2*nu*(u, div epsilon(v)) */
      aux = funct[ci] * facMtau * timefac * 2 * visc *fac;
      estif[ri*3][ci*3]     -= viscs2[0][2*ri] * aux;
      estif[ri*3][ci*3+1]   -= viscs2[1][2*ri] * aux;
      estif[ri*3+1][ci*3]   -= viscs2[0][2*ri+1] * aux;
      estif[ri*3+1][ci*3+1] -= viscs2[1][2*ri+1] * aux;


      /* -facMtau*timefac*timefac*4*nu^2(div epsilon(u), div epsilon(v)) */
      aux = facMtau * timefac * timefac * 4 * visc * visc * fac;
      estif[ri*3][ci*3]     -= (viscs2[0][2*ri]   * viscs2[0][2*ci]
                               +viscs2[1][2*ri]   * viscs2[1][2*ci]) * aux;
      estif[ri*3+1][ci*3]   -= (viscs2[0][2*ri+1] * viscs2[0][2*ci]
                               +viscs2[1][2*ri+1] * viscs2[1][2*ci]) * aux;
      estif[ri*3][ci*3+1]   -= (viscs2[0][2*ri]   * viscs2[0][2*ci+1]
                               +viscs2[1][2*ri]   * viscs2[1][2*ci+1]) * aux;
      estif[ri*3+1][ci*3+1] -= (viscs2[0][2*ri+1] * viscs2[0][2*ci+1]
                               +viscs2[1][2*ri+1] * viscs2[1][2*ci+1]) * aux;


      /* facMtau*timefac*timefac*2*nu*(u_old * grad u, div epsilon(v)) */
      aux = conv_c[ci] * facMtau * timefac * timefac * 2 * visc * fac;
      estif[ri*3][ci*3]     -= viscs2[0][2*ri] * aux;
      estif[ri*3][ci*3+1]   -= viscs2[1][2*ri] * aux;
      estif[ri*3+1][ci*3]   -= viscs2[0][2*ri+1] * aux;
      estif[ri*3+1][ci*3+1] -= viscs2[1][2*ri+1] * aux;

      
      /* facMtau*timefac*timefac*2*nu*(u * grad u_old, div epsilon(v)) */
      aux = facMtau * timefac * timefac * 2 * visc * fac;
      estif[ri*3][ci*3]     -= (viscs2[0][2*ri]   * conv_r[0][2*ci]
                               +viscs2[1][2*ri]   * conv_r[1][2*ci]) * aux;
      estif[ri*3+1][ci*3]   -= (viscs2[0][2*ri+1] * conv_r[0][2*ci]
                               +viscs2[1][2*ri+1] * conv_r[1][2*ci]) * aux;
      estif[ri*3][ci*3+1]   -= (viscs2[0][2*ri]   * conv_r[0][2*ci+1]
                               +viscs2[1][2*ri]   * conv_r[1][2*ci+1]) * aux;
      estif[ri*3+1][ci*3+1] -= (viscs2[0][2*ri+1] * conv_r[0][2*ci+1]
                               +viscs2[1][2*ri+1] * conv_r[1][2*ci+1]) * aux;


      /*- CROSS STRESS part of stabilisation --- TIME DEPENDENT FORMULATION -*/
      if(cross_stress==1)
      {
	  /* facMtau * timefac * timefac *(u, (((u_old * grad) u_old ) *grad ) v) */
	  aux = facMtau * timefac * timefac * fac;      
      
	  estif[ri*3  ][ci*3  ] += funct[ci] * ( conv_old[0] * derxy[0][ri]+ conv_old[1] * derxy[1][ri]) * aux;
	  estif[ri*3+1][ci*3+1] += funct[ci] * ( conv_old[0] * derxy[0][ri]+ conv_old[1] * derxy[1][ri]) * aux;

	  /* facMtau * timefac * timefac *(u_old, (((u * grad) u_old ) *grad ) v) */
	  aux = facMtau * timefac * timefac * fac;      
	  
	  estif[ri*3  ][ci*3  ] += velint[0] * funct[ci] * (vderxy[0][0] * derxy[0][ri] + vderxy[1][0] * derxy[1][ri]) * aux;
	  estif[ri*3  ][ci*3+1] += velint[0] * funct[ci] * (vderxy[0][1] * derxy[0][ri] + vderxy[1][1] * derxy[1][ri]) * aux;
	  estif[ri*3+1][ci*3  ] += velint[1] * funct[ci] * (vderxy[0][0] * derxy[0][ri] + vderxy[1][0] * derxy[1][ri]) * aux;
	  estif[ri*3+1][ci*3+1] += velint[1] * funct[ci] * (vderxy[0][1] * derxy[0][ri] + vderxy[1][1] * derxy[1][ri]) * aux;

	  /* facMtau * timefac * timefac *(u_old, (((u_old * grad) u ) *grad ) v) */
	  aux = facMtau * timefac * timefac * fac;      
	  
	  estif[ri*3  ][ci*3  ] += velint[0] * conv_c[ci] * derxy[0][ri] * aux;
	  estif[ri*3  ][ci*3+1] += velint[0] * conv_c[ci] * derxy[1][ri] * aux;
	  estif[ri*3+1][ci*3  ] += velint[1] * conv_c[ci] * derxy[0][ri] * aux;
	  estif[ri*3+1][ci*3+1] += velint[1] * conv_c[ci] * derxy[1][ri] * aux;
	  
	  /* facMtau*timefac*(u, u_old * grad v) */
	  aux = conv_c[ri] * (funct[ci]) * facMtau * timefac * fac;
	  estif[ri*3][ci*3]     += aux;
	  estif[ri*3+1][ci*3+1] += aux;
	  
	  /* facMtau*timefac*(u_old, u * grad v) */
	  aux = facMtau * timefac * fac;
	  
	  estif[ri*3][ci*3]     += velint[0] * ugradv[ri][2*ci]  * aux;
	  estif[ri*3][ci*3+1]   += velint[0] * ugradv[ri][2*ci+1]* aux;
	  estif[ri*3+1][ci*3]   += velint[1] * ugradv[ri][2*ci]  * aux;
	  estif[ri*3+1][ci*3+1] += velint[1] * ugradv[ri][2*ci+1]* aux;

	  /* facMtau*timefac*timefac*(u_old, grad p * grad v) */
	  aux = facMtau * timefac * timefac * fac;
	  
	  estif[ri*3  ][ci*3+2]   += velint[0] * (derxy[0][ci] * derxy[0][ri] + derxy[1][ci] *derxy[1][ri])*aux;
	  estif[ri*3+1][ci*3+2]   += velint[1] * (derxy[0][ci] * derxy[0][ri] + derxy[1][ci] *derxy[1][ri])*aux;
	  
	  /* facMtau*timefac*timefac*(u, grad p_old * grad v) */
	  aux = facMtau * timefac * timefac * fac ;
	  
	  estif[ri*3  ][ci*3  ]     += funct[ci] * (gradp[0] * derxy[0][ri] + gradp[1] * derxy[1][ri]) * aux;
	  estif[ri*3+1][ci*3+1]     += funct[ci] * (gradp[0] * derxy[0][ri] + gradp[1] * derxy[1][ri]) * aux;

	  /* -facMtau*timefac*timefac*((u, (div epsilon (u_old)) * grad )v) */
	  aux = facMtau * timefac * timefac * 2 * visc * fac;

	  estif[ri*3  ][ci*3  ]   -= funct[ci] * (visc_old[0] * derxy[0][ri] + visc_old[1] * derxy[1][ri]) * aux;
	  estif[ri*3  ][ci*3+1]   -= funct[ci] * (visc_old[0] * derxy[0][ri] + visc_old[1] * derxy[1][ri]) * aux;
	  estif[ri*3+1][ci*3  ]   -= funct[ci] * (visc_old[0] * derxy[0][ri] + visc_old[1] * derxy[1][ri]) * aux;
	  estif[ri*3+1][ci*3+1]   -= funct[ci] * (visc_old[0] * derxy[0][ri] + visc_old[1] * derxy[1][ri]) * aux;
	  

	  /* -facMtau*timefac*timefac*((u_old, (div epsilon (u)) * grad )v) */
	  aux = facMtau * timefac * timefac * 2 * visc * fac;

	  estif[ri*3  ][ci*3  ]   += velint[0] * (viscs2[0][2*ci  ] * derxy[0][ri] + viscs2[1][2*ci  ] * derxy[1][ri]) * aux;
	  estif[ri*3  ][ci*3+1]   += velint[0] * (viscs2[0][2*ci+1] * derxy[0][ri] + viscs2[1][2*ci+1] * derxy[1][ri]) * aux;
	  estif[ri*3+1][ci*3  ]   += velint[1] * (viscs2[0][2*ci  ] * derxy[0][ri] + viscs2[1][2*ci  ] * derxy[1][ri]) * aux;
	  estif[ri*3+1][ci*3+1]   += velint[1] * (viscs2[0][2*ci+1] * derxy[0][ri] + viscs2[1][2*ci+1] * derxy[1][ri]) * aux;
	  
	  
	  /* -facMtau*timefac*(u, MRHS * grad v) */
	  aux = facMtau*timefac * fac;
	  
	  estif[ri*3  ][ci*3  ]     -= funct[ci] * ((sub_vel[0]+old_vel[0]+timefac*edeadng[0]) * derxy[0][ri]
						    +
						    (sub_vel[1]+old_vel[1]+timefac*edeadng[1]) * derxy[1][ri]) * aux;
	  estif[ri*3+1][ci*3+1]     -= funct[ci] * ((sub_vel[0]+old_vel[0]+timefac*edeadng[0]) * derxy[0][ri]
						    +
						    (sub_vel[1]+old_vel[1]+timefac*edeadng[1]) * derxy[1][ri]) * aux;

	  aux = (1.-theta) * dt * facMtau * timefac * fac;
	  
	  estif[ri*3  ][ci*3  ]     += funct[ci] * ((sub_vel[0]/fdyn->tau_old[0]+res_old[0]) * derxy[0][ri]
						    +
						    (sub_vel[1]/fdyn->tau_old[0]+res_old[1]) * derxy[1][ri]) * aux;
	  estif[ri*3+1][ci*3+1]     += funct[ci] * ((sub_vel[0]/fdyn->tau_old[0]+res_old[0]) * derxy[0][ri]
						    +
						    (sub_vel[1]/fdyn->tau_old[0]+res_old[1]) * derxy[1][ri]) * aux;
      }
      
      /*------- AGLS part of stabilisation --- TIME DEPENDENT FORMULATION ---*/

      /* facMtau*timefac*timefac*(grad p, u_old * grad v) */
      estif[ri*3][ci*3+2]   += conv_c[ri] * derxy[0][ci]
	                       * facMtau * timefac * timefac * fac;
      estif[ri*3+1][ci*3+2] += conv_c[ri] * derxy[1][ci]
	                       * facMtau * timefac * timefac * fac;      


      /* facMtau*timefac*timefac*(grad p_old, u * grad v) */
      aux = facMtau * timefac * timefac * fac;

      estif[ri*3][ci*3]     += gradp[0] * ugradv[ri][2*ci  ] * aux;
      estif[ri*3][ci*3+1]   += gradp[0] * ugradv[ri][2*ci+1] * aux;
      estif[ri*3+1][ci*3]   += gradp[1] * ugradv[ri][2*ci  ] * aux;
      estif[ri*3+1][ci*3+1] += gradp[1] * ugradv[ri][2*ci+1] * aux;


      /* facMtau*timefac*(u, u_old * grad v) */
      aux = conv_c[ri] * (funct[ci]) * facMtau * timefac * fac;
      estif[ri*3][ci*3]     += aux;
      estif[ri*3+1][ci*3+1] += aux;


      /* facMtau*timefac*(u_old, u * grad v) */
      aux = facMtau * timefac * fac;
      
      estif[ri*3][ci*3]     += velint[0] * ugradv[ri][2*ci]  * aux;
      estif[ri*3][ci*3+1]   += velint[0] * ugradv[ri][2*ci+1]* aux;
      estif[ri*3+1][ci*3]   += velint[1] * ugradv[ri][2*ci]  * aux;
      estif[ri*3+1][ci*3+1] += velint[1] * ugradv[ri][2*ci+1]* aux;

      
      /* facMtau*timefac*timefac*(u_old * grad u, u_old * grad v) */
      aux = conv_c[ri] * (conv_c[ci]) * facMtau * timefac * timefac * fac;
      
      estif[ri*3][ci*3]     += aux;
      estif[ri*3+1][ci*3+1] += aux;

      
      /* facMtau*timefac*timefac*(u * grad u_old, u_old * grad v) */
      aux = facMtau * timefac * timefac * fac;
      
      estif[ri*3][ci*3]     += conv_c[ri] * ( conv_r[0][2*ci]  )* aux;
      estif[ri*3][ci*3+1]   += conv_c[ri] * ( conv_r[0][2*ci+1])* aux;
      estif[ri*3+1][ci*3]   += conv_c[ri] * ( conv_r[1][2*ci]  )* aux;
      estif[ri*3+1][ci*3+1] += conv_c[ri] * ( conv_r[1][2*ci+1])* aux;


      /* facMtau*timefac*timefac*(u_old * grad u_old, u * grad v) */
      aux = facMtau * timefac * timefac * fac;
      
      estif[ri*3][ci*3]     += conv_old[0] * ugradv[ri][2*ci  ]* aux;
      estif[ri*3][ci*3+1]   += conv_old[0] * ugradv[ri][2*ci+1]* aux;
      estif[ri*3+1][ci*3]   += conv_old[1] * ugradv[ri][2*ci  ]* aux;
      estif[ri*3+1][ci*3+1] += conv_old[1] * ugradv[ri][2*ci+1]* aux;

      
      /* -facMtau*timefac*timefac*2*nu*(div epsilon(u), u_old * grad v) */
      aux = facMtau * timefac * timefac * 2.0 * visc * fac;
      
      estif[ri*3][ci*3]     += conv_c[ri] * ( viscs2[0][2*ci]  *aux );
      estif[ri*3][ci*3+1]   += conv_c[ri] * ( viscs2[0][2*ci+1]*aux );
      estif[ri*3+1][ci*3]   += conv_c[ri] * ( viscs2[1][2*ci]  *aux );
      estif[ri*3+1][ci*3+1] += conv_c[ri] * ( viscs2[1][2*ci+1]*aux );

      
      /* -facMtau*timefac*timefac*2*nu*(div epsilon(u_old), u * grad v) */
      aux = - facMtau * timefac * time2nue * fac;
      
      estif[ri*3][ci*3]     += ( visc_old[0] * aux ) * ugradv[ri][2*ci  ];
      estif[ri*3][ci*3+1]   += ( visc_old[0] * aux ) * ugradv[ri][2*ci+1];
      estif[ri*3+1][ci*3]   += ( visc_old[1] * aux ) * ugradv[ri][2*ci  ];
      estif[ri*3+1][ci*3+1] += ( visc_old[1] * aux ) * ugradv[ri][2*ci+1];


      /* -facMtau*timefac*(MRHS, u * grad v) */
      aux = facMtau  * timefac * fac;
      
      estif[ri*3][ci*3]     -= (sub_vel[0]+old_vel[0]+timefac*edeadng[0]) * ugradv[ri][2*ci  ] * aux;
      estif[ri*3][ci*3+1]   -= (sub_vel[0]+old_vel[0]+timefac*edeadng[0]) * ugradv[ri][2*ci+1] * aux;
      estif[ri*3+1][ci*3]   -= (sub_vel[1]+old_vel[1]+timefac*edeadng[1]) * ugradv[ri][2*ci  ] * aux;
      estif[ri*3+1][ci*3+1] -= (sub_vel[1]+old_vel[1]+timefac*edeadng[1]) * ugradv[ri][2*ci+1] * aux;

      aux = (1.-theta) * dt * facMtau * timefac * fac;
      
      estif[ri*3][ci*3]     += (sub_vel[0]/fdyn->tau_old[0]+res_old[0]) * ugradv[ri][2*ci  ] * aux;
      estif[ri*3][ci*3+1]   += (sub_vel[0]/fdyn->tau_old[0]+res_old[0]) * ugradv[ri][2*ci+1] * aux;
      estif[ri*3+1][ci*3]   += (sub_vel[1]/fdyn->tau_old[0]+res_old[1]) * ugradv[ri][2*ci  ] * aux;
      estif[ri*3+1][ci*3+1] += (sub_vel[1]/fdyn->tau_old[0]+res_old[1]) * ugradv[ri][2*ci+1] * aux;


      /*--- CONTINUITY STABILISATION  --- TIME DEPENDENT FORMULATION ---*/

      /* 1/(1+dt/tau_M)*timefac*timefac*(u_old * grad u, grad q) */
      aux = facMtau * timefac * timefac * fac;
      
      estif[ri*3+2][ci*3]   += derxy[0][ri] * conv_c[ci]* aux;
      estif[ri*3+2][ci*3+1] += derxy[1][ri] * conv_c[ci]* aux;

      
      /* 1/(1+dt/tau_M)*timefac*timefac*(u * grad u_old, grad q) */
      aux = facMtau * timefac * timefac * fac;
      
      estif[ri*3+2][ci*3]   += (derxy[0][ri] * conv_r[0][2*ci]
                               +derxy[1][ri] * conv_r[1][2*ci])
	                       *aux;
      estif[ri*3+2][ci*3+1] += (derxy[0][ri] * conv_r[0][2*ci+1]
                               +derxy[1][ri] * conv_r[1][2*ci+1])
	                       * aux;

      
      /* -1/(1+theta*dt/tau_M)*timefac*timefac*2*nu*(div epsilon(u), grad q) */
      /* viscs already contains - sign!!                                   */
      aux = timefac * facMtau * fac * timefac * 2.0 * visc;
      estif[ri*3+2][ci*3]   += (derxy[0][ri] * viscs2[0][2*ci]
                               +derxy[1][ri] * viscs2[1][2*ci]) * aux;
      estif[ri*3+2][ci*3+1] += (derxy[0][ri] * viscs2[0][2*ci+1]
                               +derxy[1][ri] * viscs2[1][2*ci+1]) * aux;


      /* 1/(1+theta*dt/tau_M)*timefac*(u, grad q) */
      aux = timefac * facMtau * fac;

      estif[ri*3+2][ci*3]   += derxy[0][ri] * funct[ci] * aux;
      estif[ri*3+2][ci*3+1] += derxy[1][ri] * funct[ci] * aux;

      
      /* 1/(1+theta*dt/tau_M)*timefac*timefac* fac * (grad p, grad q) */
      aux = timefac * timefac * facMtau * fac;
      
      estif[ri*3+2][ci*3+2] += (derxy[0][ri] * derxy[0][ci]
				+derxy[1][ri] * derxy[1][ci]) * aux;


      /* TIME DEPENDENT STABILISATION --- SUBSCALE PRESSURE STABILISATION ---*/
      /* facCtau*timefac*timefac*(div u, div v) */

      aux = timefac * timefac * facCtau * fac;

      estif[ri*3  ][ci*3  ] += derxy[0][ri]*(derxy[0][ci])* aux;
      estif[ri*3  ][ci*3+1] += derxy[0][ri]*(derxy[1][ci])* aux;
      estif[ri*3+1][ci*3  ] += derxy[1][ri]*(derxy[0][ci])* aux;
      estif[ri*3+1][ci*3+1] += derxy[1][ri]*(derxy[1][ci])* aux;


      /* TIME DEPENDENT STABILISATION --- TIME DERIVATIVE */

      /* -(u * grad u_old , v) * timefac * facMtau */
      aux = timefac * facMtau * fac ;
      estif[ri*3  ][ci*3  ] -=  funct[ri] * conv_r[0][2*ci  ] * aux;
      estif[ri*3  ][ci*3+1] -=  funct[ri] * conv_r[0][2*ci+1] * aux;
      estif[ri*3+1][ci*3  ] -=  funct[ri] * conv_r[1][2*ci  ] * aux;
      estif[ri*3+1][ci*3+1] -=  funct[ri] * conv_r[1][2*ci+1] * aux;


      /* -(u_old *grad u , v) * timefac * facMtau*/
      aux = timefac * facMtau * fac ;
      estif[ri*3  ][ci*3  ] -=  funct[ri] * conv_c[ci] * aux;
      estif[ri*3+1][ci*3+1] -=  funct[ri] * conv_c[ci] * aux;


      /* -(grad p , v) * timefac * facMtau */
      aux = timefac * facMtau * fac ;
      estif[ri*3  ][ci*3+2] -=  funct[ri] * derxy[0][ci] * aux;
      estif[ri*3+1][ci*3+2] -=  funct[ri] * derxy[1][ci] * aux;


      /* -(u , v) * facMtau  */
      aux = facMtau * fac ;
      estif[ri*3  ][ci*3  ] -=  funct[ri] * funct[ci] * aux ;
      estif[ri*3+1][ci*3+1] -=  funct[ri] * funct[ci] * aux ;


      /* facMtau * timefac * 2.0 * visc *(div epsilon(u), v) */
      aux = facMtau * timefac * 2.0 * visc * fac ;
      estif[ri*3  ][ci*3  ] -= funct[ri] * ( viscs2[0][2*ci  ]*aux );
      estif[ri*3  ][ci*3+1] -= funct[ri] * ( viscs2[0][2*ci+1]*aux );
      estif[ri*3+1][ci*3  ] -= funct[ri] * ( viscs2[1][2*ci  ]*aux );
      estif[ri*3+1][ci*3+1] -= funct[ri] * ( viscs2[1][2*ci+1]*aux );

      
   }  /* end column loop (ci) */


   /**************** integrate element force vector *********************/
   /*================== Galerkin part of the RHS =======================*/
   /*--- 'Original' RHS, concentrated ---*/
   /* this expression contains (1-theta)*dt*(d/dt vel)_n (massrhs), the */
   /* old velocity and theta*dt times the body forces at the new        */
   /* timestep.                                                         */
   /* (rhsint, v) */
   eforce[ri*3]   += funct[ri] * ( rhsint[0]*fac );
   eforce[ri*3+1] += funct[ri] * ( rhsint[1]*fac );

   /* the following expression contains (d/dt sub_vel)_n !!!            */
   /* it's nothing but the right hand side of the evolution equation of */
   /* the subscales.                                                    */

   /* This is the part of the MASS RIGHT HAND SIDE  associated with the */
   /* velocity small scales.                                            */
   aux = fac * (1.-theta) * dt;
   
   eforce[ri*3]   += (-1./fdyn->tau_old[0]*sub_vel[0]-old_acc[0]-res_old[0]) * funct[ri] * aux;
   eforce[ri*3+1] += (-1./fdyn->tau_old[0]*sub_vel[1]-old_acc[1]-res_old[1]) * funct[ri] * aux;


   /*--- from Nonlinearity of Galerkin stiffness ---*/
   /* timefac*(u_old * grad u_old, v) */
   eforce[ri*3]   += funct[ri] * ( conv_old[0]*timefacfac);
   eforce[ri*3+1] += funct[ri] * ( conv_old[1]*timefacfac);

   
#if 0   
   /* -tau_M*timefac*(rhsint, -u_G * grad v) */
   if(isale)
   {
      eforce[ri*3]   += rhsint[0] * conv_g[ri] * timetauM;
      eforce[ri*3+1] += rhsint[1] * conv_g[ri] * timetauM;
   }
   /* ALE: -tau_M*timefac*timefac*(u_old * grad u_old, u_old * grad v) */
   if(isale)
   {
      eforce[ri*3]   += conv_old[0] * conv_g[ri] * ttimetauM;
      eforce[ri*3+1] += conv_old[1] * conv_g[ri] * ttimetauM;
   }
#endif

/*----------------------------------------------------------------------*/
/* ADDITIONAL TIME DEPENDENT STABILISATION PART                         */
/*----------------------------------------------------------------------*/

   /* fac*(sub_vel,v)                                      */
   aux = fac;
   eforce[ri*3]   += (sub_vel[0]) * funct[ri] * aux;
   eforce[ri*3+1] += (sub_vel[1]) * funct[ri] * aux;


   /* TIME DEPENDENT STABILISATION --- SUBSCALE PRESSURE STABILISATION */

   /* facCtau * timefac * (CRHS, div v) */
   aux = timefac * fac * facCtau;
   
   eforce[ri*3  ] +=  div[2*ri  ] * sub_pres * aux;
   eforce[ri*3+1] +=  div[2*ri+1] * sub_pres * aux;


   aux = timefac * fac * facCtau * (1.-theta) * dt;

   eforce[ri*3  ] -=  div[2*ri  ]*(sub_pres/fdyn->tau_old[2]+divu_old ) * aux;
   eforce[ri*3+1] -=  div[2*ri+1]*(sub_pres/fdyn->tau_old[2]+divu_old ) * aux;
   
   /* TIME DEPENDENT STABILISATION --- CONTINUITY STABILISATION */

   /* 1/(1+theta*dt/tau_M)*timefac*(MRHS, grad q) */
   aux = facMtau * timefac * fac;
   eforce[ri*3+2] +=
       ((  sub_vel[0]+old_vel[0]+timefac*edeadng[0])*derxy[0][ri]
 	 +(sub_vel[1]+old_vel[1]+timefac*edeadng[1])*derxy[1][ri])*aux;

   aux=(1.-theta)*dt*facMtau*timefac*fac;
   eforce[ri*3+2] -=(
        (sub_vel[0]/fdyn->tau_old[0]+res_old[0])*derxy[0][ri]*aux
       +(sub_vel[1]/fdyn->tau_old[0]+res_old[1])*derxy[1][ri]*aux
       );

   /* from linearisation of convective term                    */
   /* 1/(1+dt/tau_M)*timefac*timefac*(u_old * grad u_old, grad q) */
   aux = facMtau * timefac * timefac * fac;
   
   eforce[ri*3+2] += (conv_old[0] * derxy[0][ri]
                     +conv_old[1] * derxy[1][ri])*aux;


   /* TIME DEPENDENT STABILISATION --- HIGHER ORDER PART */
   
   /* facMtau*timefac*2*visc*(MRHS,div epsilon(v)) */
   aux = facMtau * timefac * 2 * visc * fac;
   eforce[ri*3  ] -= ((sub_vel[0]+old_vel[0]+timefac*edeadng[0]) * viscs2[0][2*ri]
       	             +(sub_vel[1]+old_vel[1]+timefac*edeadng[1]) * viscs2[1][2*ri]  ) * aux;
   eforce[ri*3+1] -= ((sub_vel[0]+old_vel[0]+timefac*edeadng[0]) * viscs2[0][2*ri+1]
		     +(sub_vel[1]+old_vel[1]+timefac*edeadng[1]) * viscs2[1][2*ri+1]) * aux;

   aux = (1.-theta) *dt * facMtau * timefac * 2 * visc * fac;

   eforce[ri*3  ] += ((sub_vel[0]/fdyn->tau_old[0]+res_old[0]) * viscs2[0][2*ri]
       	             +(sub_vel[1]/fdyn->tau_old[0]+res_old[1]) * viscs2[1][2*ri]) * aux;
   eforce[ri*3+1] += ((sub_vel[0]/fdyn->tau_old[0]+res_old[0]) * viscs2[0][2*ri+1]
		     +(sub_vel[1]/fdyn->tau_old[0]+res_old[1]) * viscs2[1][2*ri+1]) * aux;

   
   /* from linearisation of convective term */
   /* facMtau*timefac*timefac*2*nu*(u_old * grad u_old, div epsilon(v)) */
   aux = facMtau * timefac * timefac * 2 * visc * fac;
   eforce[ri*3  ] -= (conv_old[0] * viscs2[0][2*ri]
                     +conv_old[1] * viscs2[1][2*ri]) * aux;
   eforce[ri*3+1] -= (conv_old[0] * viscs2[0][2*ri+1]
                     +conv_old[1] * viscs2[1][2*ri+1]) * aux;

   /* TIME DEPENDENT STABILISATION --- TIME DERIVATIVE OF SUBSCALES */

   /* -facMtau*(MRHS,v)                                */
   aux = facMtau * fac;
   eforce[ri*3]   -= (sub_vel[0]+old_vel[0]+timefac*edeadng[0]) * funct[ri] * aux;
   eforce[ri*3+1] -= (sub_vel[1]+old_vel[1]+timefac*edeadng[1]) * funct[ri] * aux;

   aux = (1.-theta) * dt * facMtau * fac ;
   eforce[ri*3]   += (sub_vel[0]/fdyn->tau_old[0]+res_old[0]) * funct[ri] * aux;
   eforce[ri*3+1] += (sub_vel[1]/fdyn->tau_old[0]+res_old[1]) * funct[ri] * aux;

   /* from linearisation of convective term */
   /* -timefac * facMtau * fac*(u_old * grad u_old,v)                   */
   aux = timefac * facMtau * fac ;
   eforce[ri*3]   -= conv_old[0] * funct[ri] * aux;
   eforce[ri*3+1] -= conv_old[1] * funct[ri] * aux;


   /* TIME DEPENDENT STABILISATION --- LINEARISATION OF AGLS PART */

   /* facMtau*timefac*(u_old, u_old * grad v) */
   eforce[ri*3]   += conv_c[ri] * velint[0] * facMtau * timefac * fac;
   eforce[ri*3+1] += conv_c[ri] * velint[1] * facMtau * timefac * fac;

   
   /* facMtau*2*timefac*timefac*(u_old * grad u_old, u_old * grad v) */
   aux = 2.0 * facMtau * timefac * timefac * fac;
   eforce[ri*3]   += conv_old[0] * conv_c[ri] * aux;
   eforce[ri*3+1] += conv_old[1] * conv_c[ri] * aux;

   
   /* -facMtau*timefac*timefac*2*nu*(div epsilon(u_old), u_old * grad v) */
   aux = timefac * facMtau * fac * timefac * 2.0 * visc;
   eforce[ri*3]   -= conv_c[ri] * visc_old[0]*aux ;
   eforce[ri*3+1] -= conv_c[ri] * visc_old[1]*aux ;

   
   /* facMtau*timefac*timefac*(grad p_old, u_old * grad v) */
   aux = facMtau * timefac * timefac * fac;
   
   eforce[ri*3]   += conv_c[ri] * gradp[0] * aux;
   eforce[ri*3+1] += conv_c[ri] * gradp[1] * aux;
   
   /* TIME DEPENDENT FORMULATION ---                               */
   /*      --- linearisation of CROSS STRESS part of stabilisation */
   if(cross_stress==1)
   {
       /* facMtau*timefac*(u_old, ((u_old *grad ) uold )* grad v) */
       aux = 2 * facMtau * timefac * timefac * fac;      
       eforce[ri*3]   += velint[0] * (conv_old[0] * derxy[0][ri] + conv_old[1] * derxy[1][ri] ) * aux;
       eforce[ri*3+1] += velint[1] * (conv_old[0] * derxy[0][ri] + conv_old[1] * derxy[1][ri] ) * aux;

       
       /* facMtau*timefac*(u_old, u_old * grad v) */
       aux = facMtau * timefac * fac;
       eforce[ri*3]   += conv_c[ri] * velint[0] * aux;
       eforce[ri*3+1] += conv_c[ri] * velint[1] * aux;
       
       
       /* facMtau*timefac*timefac*(u_old, grad p_old * grad v) */
       aux = facMtau * timefac * timefac * fac;
       eforce[ri*3  ]   += velint[0] * (gradp[0] * derxy[0][ri] + gradp[1] * derxy[1][ri]) * aux;
       eforce[ri*3+1]   += velint[1] * (gradp[0] * derxy[0][ri] + gradp[1] * derxy[1][ri]) * aux;


       /* -facMtau*timefac*timefac*((u_old, (div epsilon (u_old)) * grad )v) */
       aux = facMtau * timefac * timefac * 2 * visc * fac;
       eforce[ri*3  ]   -= velint[0] * (visc_old[0] * derxy[0][ri] + visc_old[1] * derxy[1][ri]) * aux;
       eforce[ri*3+1]   -= velint[1] * (visc_old[0] * derxy[0][ri] + visc_old[1] * derxy[1][ri]) * aux;
       
   }

   /* TIME DEPENDENT STABILISATION --- REYNOLDS STRESS PART */
   if(reynolds_stress==1)
   {
       /* facMtau*facMtau*timefac*(sub_vel_trial , (sub_vel_trial * grad )v) */
       aux = facMtau * facMtau * timefac * fac;
       eforce[ri*3  ]   += funct[ri] * sub_vel_trial_wo_facMtau[0] *
	                   (sub_vel_trial_wo_facMtau[0] * derxy[0][ri]
			    +
			    sub_vel_trial_wo_facMtau[0] *  derxy[1][ri]
			   ) * aux;
       eforce[ri*3+1]   += funct[ri] * sub_vel_trial_wo_facMtau[1] *
	                   (sub_vel_trial_wo_facMtau[0] * derxy[0][ri]
			    +
			    sub_vel_trial_wo_facMtau[0] *  derxy[1][ri]
			   ) * aux;
   }
}     /* end row loop (ri) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
}






#endif /*D_FLUID2_TDS*/
#endif /*D_FLUID2*/
/*! @} (documentation module close)*/
