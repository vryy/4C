#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | int                   numcurve;                                      |
 | struct _DYNAMIC      *curve;                                         |
 *----------------------------------------------------------------------*/
extern int            numcurve;
extern struct _CURVE *curve;
/*----------------------------------------------------------------------*
 |  init the control of time factors                         m.gee 02/02|
 *----------------------------------------------------------------------*/
void dyn_init_curve(int actcurve,
                   int    nstep,
                   double dt,
                   double maxtime)
{
int      i;
int      numtstep;
int      fdim;
ARRAY    time_a;
ARRAY    rlcurve_a;

#ifdef DEBUG 
dstrc_enter("dyn_init_curve");
#endif
/*----------------------------------------------------------------------*/
switch(curve[actcurve].curvetyp)
{
case curve_polygonal:
   if (curve[actcurve].bystep==0)/*============= curve in absolute time */
   {
      /*---------------- size of array curve[actcurve].time is fdim x 2 */
      fdim     = curve[actcurve].time.fdim;
      /*------------------------------- the number of steps is numtstep */
      numtstep = fdim+1;
      /*---------------------------------------------- allocate storage */
      amdef("time_a",&time_a,numtstep,1,"DV");
      amdef("rlcurve_a",&rlcurve_a,numtstep,1,"DV");
      /*----------------------------set discrete time points and values */
      for (i=0; i<fdim; i++)
      {
         time_a.a.dv[i]    = curve[actcurve].time.a.da[i][0];
         rlcurve_a.a.dv[i] = curve[actcurve].value.a.da[i][0];
      }
      time_a.a.dv[fdim]    = curve[actcurve].time.a.da[fdim-1][1];
      rlcurve_a.a.dv[fdim] = curve[actcurve].value.a.da[fdim-1][1];
      /*----------- check whether the time curve is long enough in time */
 /*     
      if (nstep*dt-EPS14 > time_a.a.dv[fdim]) dserror("Curve not long enough");
      if (maxtime  > time_a.a.dv[fdim])       dserror("Curve not long enough");
 */
      /* copy the arrays time_a and rlcurve_a to the structure curve[actcurve] */
      amdel(&(curve[actcurve].time));
      amdel(&(curve[actcurve].value));
      am_alloc_copy(&time_a,&(curve[actcurve].time));
      am_alloc_copy(&rlcurve_a,&(curve[actcurve].value));
      amdel(&time_a);
      amdel(&rlcurve_a);
   }
   else/*=============================================== curve in steps */
   {
      dserror("Load curve in dynamic only impl. abstime yet");
   }
break;
case curve_explicit: /* there's nothing to initialise */
break;
case curve_none:
   dserror("Type of timecurve unknown");
break;
default:
   dserror("Type of timecurve unknown");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_init_curve */



/*----------------------------------------------------------------------*
 |  get factor at a certain time T                           m.gee 02/02|
 *----------------------------------------------------------------------*/
void dyn_facfromcurve(int actcurve,
                   double T,
                   double *fac)
{
int      i;
int      numtstep;
int      success=0;
double   time1,time2,val1,val2;

#ifdef DEBUG 
dstrc_enter("dyn_facfromcurve");
#endif
/*----------------------------------------------------------------------*/
switch(curve[actcurve].curvetyp)
{
case curve_polygonal:
   if (curve[actcurve].bystep==0)/*============= curve in absolute time */
   {
      numtstep = curve[actcurve].time.fdim;
      for (i=0; i<numtstep-1; i++)
      {
         if (curve[actcurve].time.a.dv[i]  -EPS13 <= T &&
             curve[actcurve].time.a.dv[i+1]+EPS13 >= T)
         {
            time1 = curve[actcurve].time.a.dv[i];
            time2 = curve[actcurve].time.a.dv[i+1];
            val1  = curve[actcurve].value.a.dv[i];
            val2  = curve[actcurve].value.a.dv[i+1];
            
            *fac = val1 + (val2-val1)/(time2-time1) * (T-time1);
            
            success=1;
            break;
         }
      }
      if (!success) 
      dserror("Interpolation went wrong");
   }
   else/*=============================================== curve in steps */
   {
      dserror("Load curve in dynamic only impl. abstime yet");
   }
break;
case curve_explicit:
   *fac = dyn_facexplcurve(actcurve,T);
break;   
case curve_none:
   dserror("Type of timecurve unknown");
break;
default:
   dserror("Type of timecurve unknown");
break;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of dyn_facfromcurve */

/*----------------------------------------------------------------------*
 |  get factor at a certain time T defined by any explict function      |
 |                                                           genk 06/02 |
 *----------------------------------------------------------------------*/
double dyn_facexplcurve(int actcurve,   /* number of actual time curve  */
                        double T        /* actual time                  */
		       )
{
int numex;          /* number of explicit time curve                    */
double val1, fac;
double c1,c2;

#ifdef DEBUG 
dstrc_enter("dyn_facexplcurve");
#endif

numex=curve[actcurve].numex;
c1=curve[actcurve].c1;
c2=curve[actcurve].c2;

/*----------------------------------------------------------------------*/
switch (numex)
{
   
case -1:
   if (T <= c1)
   {
      val1 = T/c1*PI/2;
      fac = sin(val1);      
   }
   else
      fac = ONE;         
   break;
   
case -2:
   if (T < EPS6) 
   {
      fac = ZERO;
   }
   else
   {
      val1 = ONE - ONE/T;
      fac = exp(val1);
   }
   break;
   
default:
   dserror("Number of explicit timecurve (NUMEX) unknown\n");
}   

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return((double)(fac));;
}   /* end of dyn_facexplcurve */
