/*!----------------------------------------------------------------------
\file
\brief service routines for fluid time algorithms with adaptive time stepping

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fluid_prototypes.h"
/*#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"*/
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
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;

/*!---------------------------------------------------------------------
\brief calculate fluid acceleration field

<pre>                                                         chfoe 10/03

the fluid acceleration is calculated depending on the actual time
stepping scheme. It is written to sol_increment[5][i]

One-step-Theta and Generalised Alpha:

acc(n+1) = (vel(n+1)-vel(n)) / (Theta * dt(n)) - (1/Theta -1) * acc(n)


BDF2:

	      2*dt(n)+dt(n-1)		     dt(n)+dt(n-1)
acc(n+1) = --------------------- vel(n+1) - --------------- vel(n)
	   dt(n)*[dt(n)+dt(n-1)]	     dt(n)*dt(n-1)

		    dt(n)
	 + ----------------------- vel(n+1)
 	   dt(n-1)*[dt(n)+dt(n-1)]

 NOTE: For the One step Theta scheme the first 10 steps are performed
       by means of the BDF2 acceleration. This is necessary to get
       zero initial field calculations running properly.

 the values at the nodes are
 sol_increment[0][i] .. vel(n-1)
 sol_increment[1][i] .. vel(n)
 sol_increment[2][i] .. lin. combin. of the other values
 sol_increment[3][i] .. vel(n+1)
 sol_increment[4][i] .. acc(n-1)
 sol_increment[5][i] .. acc(n)

</pre>
\param *actfield	FIELD		(i)	the actual field
\param  iop		INT		(i)	flag, which scheme
\return void

------------------------------------------------------------------------*/
void fluid_acceleration(	FIELD 		*actfield,
				INT 	 	iop
			)
{
DOUBLE fact1, fact2, fact3, dta, dtp, sum;
#ifdef DEBUG
dstrc_enter("fluid_acceleration");
#endif

fdyn = alldyn[genprob.numff].fdyn;
switch (iop)
{
case 1:	/* Generalised Alpha time integration 				*/
case 4:	/* One step Theta time integration 				*/
   fact1 = 1.0 / (fdyn->theta*fdyn->dta);
   fact2 =-1.0 / fdyn->theta + 1.0;	/* = -1/Theta + 1		*/
   solserv_sol_zero(actfield,0,node_array_sol_increment,5);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,3,5, fact1);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,1,5,-fact1);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,4,5, fact2);
break;
case 7:	/* 2nd order backward differencing BDF2				*/
   dta = fdyn->dta;
   dtp = fdyn->dtp;
   if (dta*dtp < EPS15)
      dserror("Zero time step size!!!!!");
   sum = dta + dtp;
   fact1 = (2.0 * dta + dtp) / (dta*sum);
   fact2 =-sum / (dta*dtp);
   fact3 = dta / (dtp*sum);
   solserv_sol_zero(actfield,0,node_array_sol_increment,5);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,3,5,fact1);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,1,5,fact2);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,0,5,fact3);
break;
default:
   dserror("Time integration scheme unknown for calculation of acceleration!");
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_acceleration */


/*!---------------------------------------------------------------------
\brief calculate fluid vector to be multiplied with the mass for time rhs

<pre>                                                         chfoe 10/03

This routine prepares the time rhs in mass form. Depending on the time
integration scheme the linear combination of vel(n) and acc(n) is
evaluated which is later multiplied by the mass matrix in order to give
the time rhs.
The linear combination is written to sol_increment[2]


Generalised Alpha:
		         Theta - alpha_m
vel(n) - dt * alpha_f * ---------------- * acc(n)
		            alpha_m

One-step-Theta:

vel(n) + dt*(1-Theta)*acc(n)


BDF2:

for constant time step:

4/3 vel(n) - 1/3 vel(n-1)

for adaptive time step:

/		   2	       \			  2
|	    (dt(n))	       |    		   (dt(n))
|1 + ------------------------- | vel(n) - ------------------------- vel(n-1)
\    dt(n-1)*[2*dt(n)+dt(n-1)] /	  dt(n-1)*[2*dt(n)+dt(n-1)]


 the values at the nodes are
 sol_increment[0][i] .. vel(n-1)
 sol_increment[1][i] .. vel(n)
 sol_increment[2][i] .. lin. combin. of the other values
 sol_increment[3][i] .. vel(n+1)
 sol_increment[4][i] .. acc(n-1)
 sol_increment[5][i] .. acc(n)


</pre>
\param *actfield	FIELD		(i)	the actual field
\param  iop		INT		(i)	flag, which scheme
\return void

------------------------------------------------------------------------*/
void fluid_prep_rhs(FIELD *actfield)
{
DOUBLE 	fact;

#ifdef DEBUG
dstrc_enter("fluid_prep_rhs");
#endif

fdyn = alldyn[genprob.numff].fdyn;

switch (fdyn->iop)
{
case 1:	/* Generalised Alpha time integration 				*/
   fact = -fdyn->dta * fdyn->alpha_f *
          (fdyn->theta - fdyn->alpha_m) / fdyn->alpha_m;
   solserv_sol_copy(actfield,0,node_array_sol_increment,node_array_sol_increment,1,2);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,5,2,fact);
break;
case 4:	/* One step Theta time integration 				*/
   fact = fdyn->dta * (1.0 -fdyn->theta);     /* = dt*(1-Theta)         */
   solserv_sol_copy(actfield,0,node_array_sol_increment,node_array_sol_increment,1,2);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,5,2,fact);
break;
case 7:	/* 2nd order backward differencing BDF2				*/
   fact = DSQR(fdyn->dta)/(fdyn->dtp*(2.0*fdyn->dta+fdyn->dtp));
   solserv_sol_zero(actfield,0,node_array_sol_increment,2);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,1,2,1+fact);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,0,2,-fact);
break;
default:
   dserror("Time integration scheme unknown for mass rhs!");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_prep_rhs */



/*!---------------------------------------------------------------------
\brief calculate predictor for adaptive time stepping

<pre>                                                         chfoe 10/03

This routine evaluates the predicted velocity at the end of the next time
step. The prediction is performed via an explicit partner of the implicit
time stepping scheme used.

One-step-Theta with Theta = 1/2 (Crack Nicolson)
and Generalised Alpha:
with 2nd order Adams-Bashford (AB2):
			     /					           \
		      dt(n) | /	  dt(n)    \	         dt(n)	           |
vel_P(n+1) = vel(n) + ----- || 2 + ------- | * acc(n) - ------- * acc(n-1) |
	 	        2   |\     dt(n-1) /            dt(n-1)            |
		            \					 	  /

One step Theta with Forward Euler: (Theta != 1/2)

vel_P(n+1) = vel(n) + dt(n) * acc(n)

BDF2 with generalised leapfrog:

			     /	 dt(n) 	  \	    /dt(n) \2
vel_P(n+1) = vel(n) + dt(n) | 1 + ------- |*acc(n)-|-------| * [u(n)-u(n-1)]
			    \ 	dt(n-1)   /	   \dt(n-1)/

 the values at the nodes are
 sol_increment[0][i] .. vel(n-1)
 sol_increment[1][i] .. vel(n)
 sol_increment[2][i] .. lin. combin. of the other values
 sol_increment[3][i] .. vel(n+1)
 sol_increment[4][i] .. acc(n-1)
 sol_increment[5][i] .. acc(n)
 sol_increment[6][i] .. predicted vel(n+1)

</pre>
\param *actfield	FIELD		(i)	the actual field
\param  iop		INT		(i)	flag, which scheme
\return void

------------------------------------------------------------------------*/
void fluid_predictor(FIELD *actfield, INT iop)
{
DOUBLE 	fact1, fact2;

#ifdef DEBUG
dstrc_enter("fluid_predictor");
#endif

fdyn = alldyn[genprob.numff].fdyn;

if(fdyn->dtp < EPS15)
   dserror("'zero' previous time step size!");

switch (iop)
{
case 1:
case 4:	/* One step Theta time integration (including TR)		*/
   if (fdyn->theta == 0.5)	/* TR */
   {
      fact1 = fdyn->dta*0.5 * (2.0 + fdyn->dta/fdyn->dtp);
      fact2 =-fdyn->dta*0.5 * fdyn->dta/fdyn->dtp;
      solserv_sol_copy(actfield,0,node_array_sol_increment,node_array_sol_increment,1,6);
      solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,5,6,fact1);
      solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,4,6,fact2);
   }
   else
   {
      solserv_sol_copy(actfield,0,node_array_sol_increment,node_array_sol_increment,1,6);
      solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,5,6,fdyn->dta);
   }
break;
case 7:	/* 2nd order backward differencing BDF2				*/
   fact1 = fdyn->dta*(1.0+fdyn->dta/fdyn->dtp);
   fact2 = DSQR(fdyn->dta/fdyn->dtp);
   solserv_sol_copy(actfield,0,node_array_sol_increment,node_array_sol_increment,1,6);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,5,6, fact1);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,1,6,-fact2);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,0,6, fact2);
break;
default:
   dserror("Time integration scheme unknown for adaptive time stepping!");
}
/*---- copy predicted velocities (no pressure) at (n+1) to sol_field ---*/
solserv_sol_copy(actfield,0,node_array_sol_increment,node_array_sol_increment,6,3);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_predictor */


/*!---------------------------------------------------------------------
\brief calculate local truncation error

<pre>                                                         chfoe 10/03

This routine evaluates the local truncation error of the actual time step
depending on the actual solution (sol_increment[3][i]) and the predicted
solution (sol_increment[6][i]). It acts according to the time integration
method iop


2nd Order Generalised Alpha:
with 2nd order Adams-Bashford:

		6*alpha_f - 2
LTE = ------------------------------- * [vel(n+1) - vel_P(n+1)]
      3*( dt(n-1)/dt(n) + 2*alpha_f )

One-step-Theta with Theta = 1/2 (Crack Nicolson):
with 2nd order Adams-Bashford:

       vel(n+1) - vel_P(n+1)
LTE = -----------------------
      3*( 1 + dt(n-1)/dt(n) )

One-step_Theta general (Theta != 1/2):
with 1st order Forward Euler:

       /      1    \
LTE = |1 - ------- | * [ vel(n+1) - vel_P(n+1) ]
      \    2*Theta /

BDF2 with generalised leapfrog:

	       /    dt(n-1) \2
	      | 1 + ------- |
	      \     dt(n)  /
LTE = -------------------------------------  * [ vel(n+1) - vel_P(n+1) ]
        dt(n-1)     /dt(n-1)\2    /dt(n-1)\3
     1+3------- + 4| -------| + 2| -------|
         dt(n)     \  dt(n) /    \  dt(n) /

 the values at the nodes are
 sol_increment[0][i] .. vel(n-1)
 sol_increment[1][i] .. vel(n)
 sol_increment[2][i] .. lin. combin. of the other values
 sol_increment[3][i] .. vel(n+1)
 sol_increment[4][i] .. acc(n-1)
 sol_increment[5][i] .. acc(n)
 sol_increment[6][i] .. predicted vel(n+1)
 sol_increment[7][i] .. local truncation error (lte)

</pre>
\param *actfield	FIELD		(i)	the actual field
\param  iop		INT		(i)	flag, which scheme
\return void

------------------------------------------------------------------------*/
void fluid_lte(	FIELD	 	*actfield,
		INT 		 iop)
{
DOUBLE 	fact, ratio;

#ifdef DEBUG
dstrc_enter("fluid_lte");
#endif

fdyn = alldyn[genprob.numff].fdyn;

switch (iop)
{
case 1:	/* Generalised Alpha time integration                           */
   fact = (6.0 * fdyn->alpha_f - 2.0) / ( 3.0*( fdyn->dtp/fdyn->dta
           + 2.0*fdyn->alpha_f ) );
   solserv_sol_zero(actfield,0,node_array_sol_increment,7);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,3,7, fact);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,6,7,-fact);
break;
case 4:	/* One step Theta time integration                              */
   if (fdyn->theta == 0.5) 	/* TR! */
   {
      fact = 1.0/3.0*(1.0+fdyn->dtp/fdyn->dta);
      solserv_sol_zero(actfield,0,node_array_sol_increment,7);
      solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,3,7, fact);
      solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,6,7,-fact);
   }
   else
   {
      fact = 1.0 - 1.0/( 2.0*fdyn->theta );
      solserv_sol_zero(actfield,0,node_array_sol_increment,7);
      solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,3,7, fact);
      solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,6,7,-fact);
   }
break;
case 7:	/* 2nd order backward differencing BDF2                         */
   ratio = fdyn->dtp/fdyn->dta;
   fact  = DSQR(1.0+ratio)/
           ( 1.0 + 3.0*ratio + 4.0*DSQR(ratio) + 2.0*DSQR(ratio)*ratio );
   solserv_sol_zero(actfield,0,node_array_sol_increment,7);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,3,7, fact);
   solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,6,7,-fact);
break;
default:
   dserror("Time integration scheme unknown for adaptive time stepping!");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_lte */


/*!---------------------------------------------------------------------
\brief calculate norm of local truncation error

<pre>                                                         chfoe 10/03

This routine evaluates the norm of the local truncation error and sets
the new time step size. It also decides whether or not a time step should
be repeated in order to secure accuracy.
The decision to repeat the timestep bases on the error norm given by
GRESHO/SANI in "Incompressible Flow and the Finite Element Method" p.708.
The maximal velocities in x- and y-direction (and z-direction) are used
to set up the characteristic velocities u0, v0 (w0), respectively.

NOTE: The proposed new time step is actually the proper time step size
      for the recent step. Hence this is repeated in case the new
      proposal is much smaller than the step size which was used.

NOTE: The repeat-timestep/don't-repeat-timestep-decision follows the
      'heuristics' in GRESHO/SANI "Incompressible Flow and the Finite
      Element Method" p. 711 or W.A.Wall (dissertation) p. 98.

</pre>
\param *actpart		PARTITION	(i) 	actual partition
\param *actfield	FIELD		(i)	the actual field
\param *actintra	INTRA		(i)
\param *fdyn		FLUID_DYNAMIC	(i/o)
\param *iststep		INT		(i/o)
\param *repeat		INT		(o)	flag, if to repeat
\param *repeated	INT		(i/o)	flag, if was repeated
\param  itnum		INT		(i)	number of iterations
\return void

------------------------------------------------------------------------*/
void fluid_lte_norm(
			PARTITION 	*actpart,
			INTRA		*actintra,
			INT		*iststep,
			INT		*repeat,
			INT		*repeated,
			INT		 itnum
			)
{

INT      i,j;		/* counters                                     */
INT      nvel;          /* number of free velocity dofs in global dir.  */
#ifdef PARALLEL
INT      get;		/* receive value (int) for MPI-process          */
DOUBLE   recv;          /* receive value (double) for MPI-process       */
#endif
INT      numveldof;

DOUBLE   d_norm = 0.0;  /* norm of LTE                                  */
DOUBLE   sum;
DOUBLE   vel0[3];       /* 'characteristic' velocity in global dir.     */
DOUBLE   getvec[3];     /* vector to receive values                     */
DOUBLE   proposed_dt;   /* the adaptively proposed new delta t          */
DOUBLE   ratio = 0.0;   /* ratio betw. recent step size and new proposal*/

NODE    *actnode;       /* the actual node                              */
GNODE   *actgnode;      /* the corresponding gnode                      */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("fluid_lte_norm");
#endif
/*------------------------------------------- initialise some values ---*/
fdyn = alldyn[genprob.numff].fdyn;

numveldof = fdyn->numdf - 1;
nvel = 0;
sum = ZERO;
vel0[0] = vel0[1] = vel0[2] = ZERO;

/*---------------------- get maximum velocity in x, y, (z)-direction ---*/
/* loop nodal points of this proc. */
for (i=0; i<actpart->pdis[0].numnp; i++)
{
  actnode  = (actpart->pdis[0].node[i]);
  if(actnode->numdf > 4)
    dserror("Too many degrees of freedom!!");
  for (j=0;j<actnode->numdf-1;j++) /* loop all velocity dofs */
    vel0[j] = DMAX(FABS(actnode->sol_increment.a.da[3][j]),vel0[j]);
}

#ifdef PARALLEL
MPI_Reduce(vel0,getvec,3,MPI_DOUBLE,MPI_MAX,0,actintra->MPI_INTRA_COMM);
#else
for (i=0; i<3; i++)  getvec[i] = vel0[i];
#endif

/*---------------- truncate characteristic velocity at lower bound ...
                                       ... to avoid division by zero ---*/
if (par.myrank == 0) for (i=0; i<3; i++) vel0[i]=DMAX(getvec[i],0.0001);

/*------------------------------------------- send vel0 to all procs ---*/
#ifdef PARALLEL
MPI_Bcast(vel0,3,MPI_DOUBLE,0,actintra->MPI_INTRA_COMM);
#endif

/*------------------------- determine norm of local truncation error ---*/

/* loop nodal points of this proc. */
for (i=0; i<actpart->pdis[0].numnp; i++)
{
  actnode  = (actpart->pdis[0].node[i]);
  actgnode = actnode->gnode;
  for (j=0;j<numveldof;j++) /* loop all velocity dofs */
  {
    if (actgnode->dirich->dirich_onoff.a.iv[j]!=0)
      continue; /* do nothing for dbc dofs*/
    sum += DSQR( actnode->sol_increment.a.da[7][j]
         / ( FABS(actnode->sol_increment.a.da[3][j]) + vel0[j] ) );
    nvel++;
  }
}

/*--------------------------------- get error sum and number of dofs ---*/
#ifdef PARALLEL
MPI_Reduce(&sum,&recv,1,MPI_DOUBLE,MPI_SUM,0,actintra->MPI_INTRA_COMM);
sum = recv;
MPI_Reduce(&nvel,&get,1,MPI_INT,MPI_SUM,0,actintra->MPI_INTRA_COMM);
nvel = get;
#endif

/*-------------------- calculate time step size ratio on PROC 0 only ---*/
if (par.myrank == 0)
{
  d_norm = sqrt( sum / nvel );

  /*----------------------------------------- propose new time step ---*/
  if (d_norm < EPS15)	/* if there was almost no error */
      ratio = 1.0;      /* simply go on with the same time step size */
  else
  {
    switch (fdyn->iop)
    {
    case 4: 	/* One step Theta					*/
      {
        if (fdyn->theta == 0.5) /* TR */
          ratio = pow((fdyn->lte/d_norm),0.33333333333333333333333);
	  else
          ratio = pow((fdyn->lte/d_norm),0.49);
      }
    break;
    case 1:	/* Generalised Alpha					*/
    case 7:	/* BDF2							*/
      ratio = pow((fdyn->lte/d_norm),0.33333333333333333333333);
    break;
    }	/* end switch (fdyn->iop) */
  }
}	/* end of par.myrank == 0 */

/*------------------------------------------ send ratio to all procs ---*/
#ifdef PARALLEL
MPI_Bcast(&ratio,1,MPI_DOUBLE,0,actintra->MPI_INTRA_COMM);
#endif

if (*repeated)       /* step has already been repeated                  */
{
   if (ratio < 0.5)   ratio = 0.5;    /* truncate time step size change */
   else if (ratio > 1.5)   ratio = 1.5;
   *repeat = 0; 		      /* don't repeat time step         */
   proposed_dt = fdyn->dta * ratio;
   proposed_dt = DMIN(proposed_dt,fdyn->max_dt);
   proposed_dt = DMAX(proposed_dt,fdyn->min_dt);
   fdyn->dt_prop = proposed_dt;
   *repeated = 0;
   goto end;
}

proposed_dt = fdyn->dta * ratio;    /* set up new proposed dt */

/* upper step size limit reached     */
if (proposed_dt > fdyn->max_dt && ratio < 5.0)
{
  *repeat = 0;
  fdyn->dt_prop = fdyn->max_dt;
  if (par.myrank == 0) printf("Time step size is cut to MAX_DT!\n");
  goto end;
}

/* lower step size limit reached     */
if (proposed_dt < fdyn->min_dt && ratio < 0.1)
{
  *repeat = 0;
  fdyn->dt_prop = fdyn->min_dt;
  if (par.myrank == 0) printf("Time step size is cut to MIN_DT!\n");
  goto end;
}

/*----------------------------------------- act with the proposal ---*/
if (ratio <= 0.1)    /* massive reduction of time step -> warning    */
{
  *repeat = 1; 		     /* repeat time step again       */
  if (par.myrank == 0)
     printf("Warning: massive reduction in time step size, something's wrong!!\n");
  fdyn->acttime -= fdyn->dta;	     /* reset old time values	     */
  fdyn->step--;
  (*iststep)--;
  fdyn->dta = proposed_dt;	     /* set smaller step size	     */
}
else if (ratio > 0.1 && ratio <= 0.8)/* significant reduction	     */
{
  *repeat = 1; 		     /* repeat time step again       */
  if (par.myrank == 0)
    printf("Significant time step size reduction. Step is repeated\n");
  fdyn->acttime -= fdyn->dta;	     /* reset old time values	     */
  fdyn->step--;
  (*iststep)--;
  fdyn->dta = proposed_dt;	     /* set smaller step size	     */
}
else if (ratio > 0.8 && ratio <= 1.5)/* about the same  	     */
{
  *repeat = 0; 		     /* don't repeat time step       */
  proposed_dt = DMIN(proposed_dt,fdyn->max_dt);
  fdyn->dt_prop = proposed_dt;
}
else if (ratio > 1.5 && ratio <= 5.0)/* significant enlargement      */
{
  *repeat = 1; 		     /* repeat time step again       */
  if (par.myrank == 0)
    printf("Significant time step size enlargement. Step is repeated\n");
  fdyn->acttime -= fdyn->dta;	     /* reset old time values	     */
  fdyn->step--;
  (*iststep)--;
  fdyn->dta = proposed_dt;
}
else if (ratio > 5.0)		     /* huge enlargement	     */
{
  if (fdyn->dta == fdyn->max_dt)  /* max timestep already reached */
  {
  *repeat = 0; 		     /* don't repeat time step       */
  proposed_dt = fdyn->max_dt;
  fdyn->dt_prop = proposed_dt;
  }
  else
  {
     *repeat = 1;                   /* repeat time step again       */
  if (par.myrank == 0)
    printf("Huge time step size enlargement. Step is repeated\n");
    fdyn->acttime -= fdyn->dta;     /* reset old time values        */
    fdyn->step--;
    (*iststep)--;
    fdyn->dta = DMIN(5.0 * fdyn->dta,fdyn->max_dt);
  }
}
else
  dserror("Something's wrong in adaptive time stepping.");

end:
/*----------------------------------------------------------------------*/
if (par.myrank == 0)
  if (*repeat == 0)
    plot_lte(fdyn->acttime,fdyn->step,d_norm,fdyn->dta,itnum);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_lte_norm */


#endif
/*! @} (documentation module close)*/
