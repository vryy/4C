/*!---------------------------------------------------------------------
\file
\brief routines writing data to plot file 0.plt

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6127
</pre>

---------------------------------------------------------------------*/
/*!
\addtogroup OUTPUT
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "../ale2/ale2.h"

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;

/*!----------------------------------------------------------------------
\brief plot lift, drag and angular momentum over time

<pre>                                                       chfoe 01/04
Writes lift and drag forces as well as the angular momentum and the actual
time to 0.plt file to serve gnuplot plots.
This routine is usefull for 2D fluid calculations only.

liftdrag[0]     drag force (global x-direction)
liftdrag[1]     lift force (global y-direction)
liftdrag[2]     angular momentum
</pre>

\warning There is nothing special to this routine
\return void
\sa

*-----------------------------------------------------------------------*/
#ifdef D_FLUID
void plot_liftdrag(
    DOUBLE   time,
    DOUBLE  *liftdrag
    )
{

  INT            i;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("plot_liftdrag");
#endif

/*------------------------------------------------------- check proc ---*/
#ifdef PARALLEL
if (par.myrank != 0) goto end;
#endif

fprintf(allfiles.gnu,"%10.5f  ", time);
for (i=0; i<2*(FLUID_NUM_LD+1); i++)
  fprintf(allfiles.gnu,"%8.7f  %8.7f  %8.7f  %8.7f  %8.7f  %8.7f  ",
      liftdrag[i*6+0], liftdrag[i*6+1], liftdrag[i*6+2],
      liftdrag[i*6+3], liftdrag[i*6+4], liftdrag[i*6+5]);
fprintf(allfiles.gnu,"\n");

fflush(allfiles.gnu);

/*----------------------------------------------------------------------*/
#ifdef PARALLEL
end:
#endif

#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of plot_liftdrag */



/*!---------------------------------------------------------------------
\brief plot fluid step size and LTE norm to .plt

<pre>                                                         chfoe 11/03

This routine simply plots fluid dynamic variables to *.plt-file.

</pre>
\param	time		DOUBLE		(i)	actual simulation time
\param	step		INT		(i)	step number
\param	norm		DOUBLE		(i)	actual LTE-norm
\param 	dt		DOUBLE		(i)	actual step size
\param 	itnum		INT		(i)	number of iterations needed
\return void

\warning This routine is automatically used as soon as adaptive time stepping
         is switched on. This could collide with other plotting options.
------------------------------------------------------------------------*/
void plot_lte(	DOUBLE 	time,
			INT 	step,
			DOUBLE 	norm,
			DOUBLE 	dt,
			INT 	itnum)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("plot_lte");
#endif
/*------------------------------------------------------- check proc ---*/
#ifdef PARALLEL
if (par.myrank != 0) goto end;
#endif
/*----------------------------------------------------------------------*/
if (step == 1)
{
   fprintf(allfiles.gnu,"# total time  step nr.    LTE norm      stepsize  itnum\n");
   fprintf(allfiles.gnu,"%12.7f      %i  %14.12f  %12.8f  %i\n",
           0.0, 0, 0.0, dt, 0);
}
fprintf(allfiles.gnu,"%12.7f      %i  %14.12f  %12.7f  %i\n",
        time, step, norm, dt, itnum);
fflush(allfiles.gnu);


#ifdef PARALLEL
end:
#endif

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of plot_lte */

#endif /* D_FLUID */


/*!----------------------------------------------------------------------
\brief  element quality statistics

<pre>                                                             ck 06/03
element quality statistics

</pre>
\param *actfield  FIELD     (i)   actual field
\param  step      INT       (i)   actual time step
\param *actintra  INTRA     (i)   intra-communicator
\param *actpart   PARTITION (i)   actual partition

\warning Works at the moment for ale2-elements only!
\return void
\sa calling:
             called by: dyn_ale();

*----------------------------------------------------------------------*/
#ifdef D_ALE
void plot_ale_quality(FIELD *field,INT step, INTRA *actintra,
                      PARTITION *actpart)
{
INT i;      /* a counter */
INT numel;  /* number of elements in this discretisation */
INT numele_total;

DOUBLE quality;     /* current element quality measure */
DOUBLE square = 0;
DOUBLE min, max;    /* minimal and maximal quality */
DOUBLE stand_degr;  /* standard degression*/
DOUBLE average;
ELEMENT *actele;

#ifdef PARALLEL
DOUBLE recv=ZERO;
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("plot_ale_quality");
#endif

/*----------------------------------------------------------------------*/
numel = actpart->pdis[0].numele;
numele_total=field->dis[0].numele;
average = 0.0;
min = 2.0;
max = -1.0;

for (i=0; i<numel; i++)  /* loop over all elements */
{
   actele = actpart->pdis[0].element[i];
   if (actele->proc!=par.myrank) continue;
   quality = actele->e.ale2->quality;
   average += quality;
   min = DMIN(quality,min);
   max = DMAX(quality,max);
   square += quality * quality;
}

#ifdef PARALLEL
MPI_Reduce(&square,&recv,1,MPI_DOUBLE,MPI_SUM,0,actintra->MPI_INTRA_COMM);
square=recv;
MPI_Reduce(&average,&recv,1,MPI_DOUBLE,MPI_SUM,0,actintra->MPI_INTRA_COMM);
average=recv;
MPI_Reduce(&min,&recv,1,MPI_DOUBLE,MPI_MIN,0,actintra->MPI_INTRA_COMM);
min=recv;
MPI_Reduce(&max,&recv,1,MPI_DOUBLE,MPI_MAX,0,actintra->MPI_INTRA_COMM);
max=recv;
#endif

if (par.myrank==0)
{
   if (numele_total > 1)
      stand_degr = 1.0/(numele_total-1.0) * ( square - 1.0/numele_total*average*average );
   else
      stand_degr = 0.0;
   average = average/numele_total;
/*----------------------------------------------------- gnuplot file ---*/
   if(step == 1)
   fprintf(allfiles.gnu,"# step nr.  av. quality  standard degr.  min quality max quality\n");
   fprintf(allfiles.gnu,"%i  %8.7f  %8.7f  %8.7f  %8.7f\n",
           step-1, average, stand_degr, min, max);
   fflush(allfiles.gnu);
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end plot_ale_quality */
#endif /* ALE */

/*! @} (documentation module close)*/
