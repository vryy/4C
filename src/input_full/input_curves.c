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
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 *----------------------------------------------------------------------*/
extern struct _DYNAMIC *dyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _DYNAMIC      *curve;                                         |
 *----------------------------------------------------------------------*/
INT            numcurve;
struct _CURVE *curve;
/*----------------------------------------------------------------------*
 | input of curves                                        m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_cond_curve()
{
INT  ierr;
INT  i;
#ifdef DEBUG
dstrc_enter("inp_cond_curve");
#endif
/*----------------------------------------------------------------------*/
/*------------------------ count the number of different curves (max=5) */
numcurve=0;
if (frfind("--CURVE1")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numcurve)++;
}

if (frfind("--CURVE2")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numcurve)++;
}

if (frfind("--CURVE3")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numcurve)++;
}

if (frfind("--CURVE4")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numcurve)++;
}

if (frfind("--CURVE5")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numcurve)++;
}


/*------------------------------------------------- allocate the curves */
curve = (CURVE*)CCACALLOC(numcurve,sizeof(CURVE));
for (i=0; i<numcurve; i++) curve[i].Id=0;
/*----------------------------------------------------- read the curves */
inp_read_curve("--CURVE1");
inp_read_curve("--CURVE2");
inp_read_curve("--CURVE3");
inp_read_curve("--CURVE4");
inp_read_curve("--CURVE5");
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inp_cond_curve */







/*----------------------------------------------------------------------*
 | input of curves                                        m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_read_curve(char *string)
{
INT      ierr;
INT      i;
INT      counter=0;
char     buffer[50];
CURVE   *actcurve = NULL;
char    *colpointer;
#ifdef DEBUG
dstrc_enter("inp_read_curve");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------- check whether there is info on this curve */
if (frfind(string)==0) goto end;
frread();
frchk("---",&ierr);
if (ierr==1) goto end;
/*------------------------------------------- find the first free curve */
for (i=0; i<numcurve; i++)
{
   if (curve[i].Id==0)
   {
      actcurve=&(curve[i]);
      break;
   }
}
/*------------------------------- count the lines of info on this curve */
while (ierr==0)
{
   counter++;
   frread();
   frchk("---",&ierr);
}
/* allocate space to store values */
amdef("curve_time" ,&(actcurve->time) ,counter,2,"DA");
amdef("curve_value",&(actcurve->value),counter,2,"DA");
/*------------------------------------------------------- start reading */
if (frfind(string)==0) goto end;
for (i=0; i<counter; i++)
{
   frread();
/*------------------------------------------------------- read curve Id */
   frint("CURVE",&(actcurve->Id),&ierr);
   if (ierr!=1) dserror("cannot read CURVE");
/*--------------------------------------------------- read typ of curve */
   frchk("Polygonal",&ierr);
   if (ierr==1)
   {
     actcurve->curvetyp = curve_polygonal;
     actcurve->bystep=0;
     /* read by abstime */
     colpointer = strstr(allfiles.actplace,"BYABSTIME");
     colpointer += 13;
     actcurve->time.a.da[i][0] = strtod(colpointer,&colpointer);
     actcurve->time.a.da[i][1] = strtod(colpointer,&colpointer);
     colpointer = strstr(allfiles.actplace,"FACTOR");
     colpointer += 6;
     actcurve->value.a.da[i][0] = strtod(colpointer,&colpointer);
     actcurve->value.a.da[i][1] = strtod(colpointer,&colpointer);
     frdouble("T",&(actcurve->T),&ierr);
   }
   frchk("Explicit",&ierr);
   if (ierr==1)
   {
      actcurve->curvetyp = curve_explicit;
      frchar("FUNC",buffer,&ierr);
      if (ierr!=1) dserror("cannot read CURVE\n");
      if (strncmp(buffer,"f(t)=sin(t:C1*PI:2)_for_t<_C1_else_f(t)=1",41)==0)
         actcurve->numex=-1;
      else if (strncmp(buffer,"f(t)=exp(1-1:t)_for_t<C1_else_const.",36)==0)
         actcurve->numex=-2;
      else if (strncmp(buffer,"f(t)=1-cos(2*PI*C1*t)",19)==0)
         actcurve->numex=-3;
      else if (strncmp(buffer,"f(t)=C2*sin(2PI*C1*t)",21)==0)
         actcurve->numex=-4;
      else if (strncmp(buffer,"f(t)=(sin(PI(t:C1-0.5))+1)*0.5",30)==0)
         actcurve->numex=-5;
      else if (strncmp(buffer,"BELTRAMI",8)==0)
         actcurve->numex=-6;
      else if (strncmp(buffer,"KIM-MOIN",8)==0)
         actcurve->numex=-7;
      else if (strncmp(buffer,"f(t)=(C2:2PI*C1)*cos(2PI*C1*t)",30)==0)
         actcurve->numex=-8;
      else if (strncmp(buffer,"f(t)=t:2-C1:(2PI)*cos(PI*t:C1-PI:2)",35)==0)
         actcurve->numex=-9; /* time integral of numex -5 */
      else
         dserror("cannot read function of CURVE\n");
      if (ierr!=1) dserror("cannot read CURVE");
      frdouble("c1",&(actcurve->c1),&ierr);
      if (ierr!=1) dserror("cannot read CURVE");
      frdouble("c2",&(actcurve->c2),&ierr);
      if (ierr!=1) dserror("cannot read CURVE");
   }
} /* end of loop over curve lines */
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inp_read_curve */
