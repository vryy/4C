/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../pss_full/pss_parser.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _DYNAMIC      *curve;                                         |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
extern struct _CURVE *curve;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 |  init the control of time factors                         m.gee 02/02|
 *----------------------------------------------------------------------*/
void dyn_init_curve(INT actcurve,
                   INT    nstep,
                   DOUBLE dt,
                   DOUBLE maxtime)
{
INT      i;
INT      numtstep;
INT      fdim;
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
case curve_expr:
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
void dyn_facfromcurve(INT actcurve,
                      DOUBLE T,
                      DOUBLE *fac)
{
INT      i;
INT      numtstep;
INT      success=0;
DOUBLE   time1,time2,val1,val2;

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
case curve_expr:
{
  DOUBLE t;
  t = T;

  /* the function expression must be evaluated inside the valid range */
  if (t<curve[actcurve].c1)
    t = curve[actcurve].c1;
  else if (t>curve[actcurve].c2)
    t = curve[actcurve].c2;

  *fac = pss_evaluate_curve(curve[actcurve].funct,t);
  break;
}
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

/*!---------------------------------------------------------------------
\brief get factor at a certain time T defined by any explict function

<pre>                                                         genk 06/02

implemented time functions
numex = -1: f(T) = sin(T/c1*PI/2) for T<c1 else f(T) = 1
numex = -2: f(T) = exp(1-1/T) for T<c1 else f(T) = f(c1)
numex = -3: f(T) = 1-cos(c1*PI*T)
numex = -4: f(T) = exp(-c1*nu*d*d*T)                       beltrami flow
numex = -5: f(T) = exp(-c1*a*a*PI*PI*nu*T)                 kim-moin flow
numex = -6..-9 bitte unten schauen
numex = -10: f(T) =  0.5*Ppeep*(1 - cos(1/Phase*PI*T)) for T < Phase
                    else f(T) = 0.5*(1-Ppeep)*(1 - cos(2*PI*Freq*(T - Phase))) + Ppeep
					lung load curve c1 = Freq, c2 = Ppeep, c3 = Phase
					Maintainer: Robert Metzke (metzke@lnm.mw.tum.de)
tbc

</pre>

\param     actcurve     INT            (i)    number of actual time curve
\param     T            DOUBLE         (i)    actual time
\return DOUBLE

------------------------------------------------------------------------*/
DOUBLE dyn_facexplcurve(INT actcurve,   /* number of actual time curve  */
                        DOUBLE T        /* actual time                  */
		       )
{
INT numex;          /* number of explicit time curve 		     * */
INT i, FileLength=0, h, p, num, size, allocation=1000;
DOUBLE val1, fac=0.0;
DOUBLE c1,c2,c3,c4,c5,c6,c7;    /* function constants                               */
DOUBLE d,visc;      /* parameters for Beltrami-flow                     */
DOUBLE a;           /* parameters for Kim-Moin flow */
DOUBLE t,tinsp,texp,tplateau,tstart,pmaxexp; /* parameters for Lung Curve */
DOUBLE *ArrayLength, *EvenCoeffcient, *OddCoefficient, C;
DOUBLE s0;
static DOUBLE savefac;
CHAR * string=NULL;
FILE *WaveForm;


#ifdef DEBUG
dstrc_enter("dyn_facexplcurve");
#endif

numex=curve[actcurve].numex;
c1=curve[actcurve].c1;
c2=curve[actcurve].c2;
c3=curve[actcurve].c3;
c4=curve[actcurve].c4;
c5=curve[actcurve].c5;
c6=curve[actcurve].c6;
c7=curve[actcurve].c7;

/*----------------------------------------------------------------------*/
switch (numex)
{
case -1: /* f(t)=sin(t:C1*PI:2)_for_t<_C1_else_f(t)=1 */
   if (T <= c1)
   {
      val1 = T/c1*PI/2;
      fac = sin(val1);
   }
   else
      fac = ONE;
break;
case -2: /* f(t)=exp(1-1:t)_for_t<C1_else_const. */
   if (T < EPS6)
   {
      fac = ZERO;
   }
   else if (T<=c1 || c1<EPS6)
   {
      val1 = ONE - ONE/T;
      fac = exp(val1);
      savefac = fac;
   }
   else
      fac = savefac;
break;
case -3: /* f(t)=1-cos(2*PI*C1*t) */
   val1 = TWO*PI*c1*T;
   fac  = ONE - cos(val1);
break;
case -4: /* f(t)=C2*sin(2PI*C1*t) */
   val1 = TWO*c1*PI*T;
   fac  = c2*sin(val1);
break;
case -5: /* f(t)=(sin(PI(t:C1-0.5))+1)*0.5 */
   if (T<=c1)
   {
      val1 = PI*(T/c1-ONE/TWO);
      fac = (sin(val1)+ONE)/TWO;
   }
   else
      fac = ONE;
break;
case -6: /* Beltrami-Flow */
   visc = mat[0].m.fluid->viscosity;
   d = PI/TWO;
   val1 = -c1*visc*d*d*T;
   fac = exp(val1);
break;
/* Maltes "-5" auf "-7" geaendert */
case -7: /* Kim-Moin-Flow */
   visc = mat[0].m.fluid->viscosity;
   a = 2.0;
   val1 = -c1*a*a*PI*PI*visc*T;
   fac = exp(val1);
break;
case -8: /* f(t)=(C2/2PI*C1)*cos(2PI*C1*t) +s0*/
   val1 = TWO*c1*PI;
   s0   = -c2/val1;
   fac = c2/val1*cos(val1*T)+s0;
break;
case -9: /* f(t)=t:2-C1:(2PI)*cos(PI*t:C1-PI:2) */
   if (T<=c1)
   {
      val1 = PI / c1;
      fac = T*0.5 - 0.5/val1 * cos(val1*T-PI*0.5);
   }
   else
      fac = T - c1 * 0.5;
break;
case -10: /* Lung functon shifted sinus with phase */
   if (c3 == 0.0) c3 = 1/c1;
   if (T<=c3) {
	 fac = 0.5*c2*(1 - cos(1/c3*PI*T));
   } else {
	 fac = 0.5*(1-c2)*(1 - cos(2*PI*c1*(T - c3))) + c2;
   }
break;
case -11: /* f(t)=-0.5*cos(PI*(T-c1)/c2)+0.5 */
   	if (T<c1)
    	fac = 0.0;
    else if ((c1 <= T) && (T <= (c1 + c2)))
        fac = -0.5 * cos(PI * ( T - c1 ) / c2) + 0.5;
    else
        fac = 1.0;
break;
case -12: /*Function for alveolar pressure during mechanical ventilation ah/rm 06/07 */
	if (c1<=0) dserror("parameter Frequ must be > 0");
	if (c2<0 || c2>1) dserror("parameter Ppeep must be beetween 0 and 1");
	if (c3<=0) dserror("parameter I/E must be > 0");
	if (c4<0) dserror("parameter IPause/I must be >= 0");
	if (c5<=0) dserror("parameter ShapeIPauseB must be > 0");
	if (c6<c2 || c6>1) dserror("parameter ShapeIPauseC must be beetween Ppeep and 1");
	if (c7<=0) dserror("parameter ShapeEB must be > 0");

	tinsp=c3*(1-c4)/(c1*(1+c3)); /*inspiration time (without end-inspiratory plateau)*/
	texp=1/(c1*(c3+1));     /*expiaration time*/
	tplateau=c4*c3/(c1*(c3+1)); /*period of the rest at the end of inspiration*/
	pmaxexp=(1-c6)*exp(-c5*tplateau)+c6; /*pressure at the beginning of expiration*/
	tstart=2*c2*tinsp/(1-c2); /*time for starting phase*/
	t=T-tstart-(int)((T-tstart)*c1)/c1; /*t runs through one respiration period*/

	if (T<=tstart) /*starting phase, quadratic function*/
	{
		fac=pow(((1-c2)/2/tinsp),2)/c2*pow(T,2);
	} else
	{
		if (t<=tinsp) /*inspiration, linear function*/
		{
			fac = c2+(1-c2)*t/tinsp;
		} else
		{
			if (t<=tplateau+tinsp) /*rest at the end of inspiration, exponential function*/
			{
				fac=(1-c6)*exp(-c5*(t-tinsp))+c6;
			} else /*expiaration, exponential function*/
			{
				fac = (c2-pmaxexp)/(exp(-c7*texp)-1)*exp(-c7*(t-tinsp-tplateau))
				      +pmaxexp-(c2-pmaxexp)/(exp(-c7*texp)-1);
			}
		}
	}
break;
#if 0
case -13: /*Fourier decomposition of respirator data or physiological
	   * blood waveform*/

  ArrayLength = CCAMALLOC(allocation * sizeof(double));
  EvenCoefficient = CCAMALLOC(allocation * sizeof(double));
  OddCoefficient = CCAMALLOC(allocation * sizeof(double));


chdir("$HOME");
if ((WaveForm = fopen("waveform.txt", "r"))==NULL)
       dserror("Error: Unable to open %s for reading","waveform.txt");
chdir("$OLDPWD");


while (getline(&string,&size,WaveForm) !=EOF){

	ArrayLength[FileLength++] = atof(strtok(string,"\n"));
}

 fclose(WaveForm);

 C = (double)FileLength;

 for (p=0; p<=FileLength/2; p++){
   EvenCoefficient[p] = 0;
   OddCoefficient[p] = 0;

   for (num=0; num<=FileLength-1; num++){
     EvenCoefficient[p] = EvenCoefficient[p]+2/C*ArrayLength[num]*cos(2*PI*p*(num+1)/C);
     OddCoefficient[p] = OddCoefficient[p]+2/C*ArrayLength[num]*sin(2*PI*p*(num+1)/C);
   }
 }

EvenCoefficient[FileLength/2] = EvenCoefficient[FileLength/2]/2;
OddCoefficient[FileLength/2] = 0;
fac = EvenCoefficient[0]/2;

  for (h=1; h<=FileLength/2; h++){
    fac = fac+EvenCoefficient[h]*cos(2*PI*h*T/c1)+OddCoefficient[h]*sin(2*PI*h*T/c1);
  }

  CCAFREE(ArrayLength);
  CCAFREE(EvenCoefficient);
  CCAFREE(OddCoefficient);
  break;
#endif
default:
   dserror("Number of explicit timecurve (NUMEX) unknown\n");
} /* end switch(numex) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return((DOUBLE)(fac));;
}   /* end of dyn_facexplcurve */
