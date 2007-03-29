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
\brief evaluate right hand side for fluid 

<pre>                                                        gammi 11/06

In this routine the Gauss point contributions to the Galerkin part of
the elemental right hand side of a stabilised fluid2 element are
calculated.

The procedure is based on the linearisation discribed in dis. Whiting.

integration schemes:

Gen-Alpha:


for further comments see comment lines within code.


</pre>
\param  DOUBLE        (o)  
\param  DOUBLE        (i)  
\param  DOUBLE        (i)  
\param  DOUBLE        (i)  
\param  DOUBLE        (i)  
\param  DOUBLE        (i)  
\param  DOUBLE        (i)  
\param  INT           (i)  
                            


\return void

------------------------------------------------------------------------*/

void f2_calgalrhs_gen_alpha_tds(
                DOUBLE  *eforce,
		DOUBLE  *velint,
		DOUBLE  *accint,
		DOUBLE   presint,
		DOUBLE  *edeadng,
		DOUBLE  *funct,
		DOUBLE **derxy,
		DOUBLE **derxy2,
		DOUBLE **vderxy,
                DOUBLE **vderxy2,
		DOUBLE   fac,
		DOUBLE   visc,
		int      iel     
              )
{
 int       ri;       /* counter for row index                       */

#ifdef DEBUG
 dstrc_enter("f2_calgalrhs_gen_alpha_tds");
#endif
 
 for (ri=0; ri<iel; ri++)       /* row index    */
 {

     /* (accint,v) */
     eforce[ri*3]   -= funct[ri] * accint[0] * fac;
     eforce[ri*3+1] -= funct[ri] * accint[1] * fac;
     
     /* - ( p_old , div v ) */
     eforce[ri*3  ] += derxy[0][ri] * presint * fac;
     eforce[ri*3+1] += derxy[1][ri] * presint * fac;
     
     /* (2 * nu * epsilon(u_old), epsilon(v)) */
     eforce[ri*3  ] -= visc * fac 
	 * ( (2*vderxy[0][0]) * derxy[0][ri]+
	     (vderxy[0][1]+vderxy[1][0]) * derxy[1][ri]);
     
     eforce[ri*3+1] -= visc * fac 
	 * ( (vderxy[0][1]+vderxy[1][0]) * derxy[0][ri]+
	     (2*vderxy[1][1]) * derxy[1][ri]);

     /* - (f, v) */
/* This expression is already contained in the modified subscale
 * acceleration!!!                                                  */
     
     /* (div u_old, q) */
     eforce[ri*3+2] -=  (vderxy[0][0] + vderxy[1][1]) * funct[ri] * fac;
     
 } /* end loop over row index */
 /*-----------------------------------------------------------------*/

#ifdef DEBUG
 dstrc_exit();
#endif

 return;
}

/*!---------------------------------------------------------------------
\brief evaluate stabilisation part of right hand side for fluid 

<pre>                                                        gammi 11/06

In this routine the Gauss point contributions to the stabilisation part
of the elemental right hand side of a stabilised fluid2 element are
calculated.

The procedure is based on the linearisation discribed in dis. Whiting.

integration schemes:

Gen-Alpha:


for further comments see comment lines within code.


</pre>
\param  DOUBLE        (o)  
\param  DOUBLE        (i)  
\param  DOUBLE        (i)  
\param  DOUBLE        (i)  
\param  DOUBLE        (i)  
\param  DOUBLE        (i)  
\param  DOUBLE        (i)  
\param  INT           (i)  
                            


\return void

------------------------------------------------------------------------*/

void f2_calstabrhs_gen_alpha_tds(
                DOUBLE  *eforce,
		DOUBLE  *velint,
		DOUBLE  *accint,
		DOUBLE   presint,
		DOUBLE  *gradpint,
		DOUBLE  *edeadng,
		DOUBLE  *funct,
		DOUBLE **derxy,
		DOUBLE **derxy2,
		DOUBLE **vderxy,
                DOUBLE **vderxy2,
		DOUBLE   svel[2],
		DOUBLE   sacc[2],
		DOUBLE   spres,
		DOUBLE   fac,
		DOUBLE   visc,
		int      iel     
              )
{
 int       ri;       /* counter for row index                       */

 int       i;

 DOUBLE  viscs2[2][2*MAXNOD]; /* viscous term incluiding 2nd derivatives */

/* time integration */
 FLUID_DYNAMIC   *fdyn;
 
 double  alpha_M;

 
#ifdef DEBUG
 dstrc_enter("f2_calstabrhs_gen_alpha_tds");
#endif


fdyn = alldyn[genprob.numff].fdyn;
alpha_M     = fdyn->alpha_m;

 
 /* Viscous test term  div epsilon(v) */
 for(i=0;i<iel;i++)
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

 }

 
 for (ri=0; ri<iel; ri++)       /* row index    */
 {
     /* -2 * nu * (u_sub , div epsilon(v)) */
     /* viscs2 already contains a - sign */
     eforce[ri*3  ] -= (svel[0]*viscs2[0][2*ri  ]
			+
			svel[1]*viscs2[1][2*ri  ])*2*visc*fac;

     eforce[ri*3+1] -= (svel[0]*viscs2[0][2*ri+1]
			+
			svel[1]*viscs2[1][2*ri+1])*2*visc*fac;


     /* - (p_sub , div v) */
     eforce[ri*3  ] += spres * derxy[0][ri] * fac;
     eforce[ri*3+1] += spres * derxy[1][ri] * fac;

     /* -(u_sub^{n+alpha_F}_i,grad q) */
     eforce[ri*3+2] += (derxy[0][ri] * svel[0]
			+
			derxy[1][ri] * svel[1])* fac;

     /* (sacc^{n+alpha_M},v) */
     eforce[ri*3  ] -= sacc[0] * funct[ri] * fac;
     eforce[ri*3+1] -= sacc[1] * funct[ri] * fac;

     /* this one is bogus, cause it isn't n+alpha_M but n+1 ... */
     
 } /* end loop over row index */

 /*-----------------------------------------------------------------*/
#ifdef DEBUG
 dstrc_exit();
#endif

 return;
}

#endif /*D_FLUID2_TDS*/
#endif /*D_FLUID2*/
/*! @} (documentation module close)*/
