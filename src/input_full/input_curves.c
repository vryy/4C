#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 *----------------------------------------------------------------------*/
#ifdef SUSE73
extern DYNAMIC *dyn;   
#else
extern struct _DYNAMIC *dyn;   
#endif
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | int                   numcurve;                                      |
 | struct _DYNAMIC      *curve;                                         |
 *----------------------------------------------------------------------*/
int            numcurve;
struct _CURVE *curve;
/*----------------------------------------------------------------------*
 | input of curves                                        m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_cond_curve()
{
int  ierr;
int  i;
int  counter;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("inp_cond_curve");
#endif
/*----------------------------------------------------------------------*/
/*------------------------ count the number of different curves (max=5) */
numcurve=0;
frfind("--CURVE1");
frread();
frchk("---",&ierr);
if (ierr==0) (numcurve)++;

frfind("--CURVE2");
frread();
frchk("---",&ierr);
if (ierr==0) (numcurve)++;

frfind("--CURVE3");
frread();
frchk("---",&ierr);
if (ierr==0) (numcurve)++;

frfind("--CURVE4");
frread();
frchk("---",&ierr);
if (ierr==0) (numcurve)++;

frfind("--CURVE5");
frread();
frchk("---",&ierr);
if (ierr==0) (numcurve)++;



/*------------------------------------------------- allocate the curves */
curve = (CURVE*)CALLOC(numcurve,sizeof(CURVE));
if (curve==NULL) dserror("Allocation of CURVE failed");
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
int      ierr;
int      i;
int      counter=0;
char     buffer[50];
CURVE   *actcurve;
char    *colpointer;
#ifdef DEBUG 
dstrc_enter("inp_read_curve");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------- check whether there is info on this curve */
frfind(string);
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
/*-------------------------------------- allocate space to store values */
frfind(string);
frread();
frchar("BYSTEP",buffer,&ierr);
if (strncmp(buffer,"Yes",3)==0 ||
    strncmp(buffer,"YES",3)==0 ||
    strncmp(buffer,"yes",3)==0 )  
{
   actcurve->bystep=1;
   amdef("curve_time" ,&(actcurve->time) ,counter,2,"IA");
   amdef("curve_value",&(actcurve->value),counter,2,"DA");
}
else                              
{
   actcurve->bystep=0;
   amdef("curve_time" ,&(actcurve->time) ,counter,2,"DA");
   amdef("curve_value",&(actcurve->value),counter,2,"DA");
}
/*------------------------------------------------------- start reading */
frfind(string);
for (i=0; i<counter; i++)
{
   frread();
/*------------------------------------------------------- read curve Id */
   frint("CURVE",&(actcurve->Id),&ierr);
   if (ierr!=1) dserror("cannot read CURVE");
/*--------------------------------------------------- read typ of curve */   
   frchk("Polygonal",&ierr);
   if (ierr==1) actcurve->curvetyp = curve_polygonal;
/*-------------------------------------------- read bystep or byabstime */
   if (actcurve->bystep==1)
   {
       colpointer = strstr(allfiles.actplace,"BYSTEP");
       colpointer += 10;
       actcurve->time.a.ia[i][0] = strtol(colpointer,&colpointer,10);
       actcurve->time.a.ia[i][1] = strtol(colpointer,&colpointer,10);
       colpointer = strstr(allfiles.actplace,"FACTOR");
       colpointer += 6;
       actcurve->value.a.da[i][0] = strtod(colpointer,&colpointer);
       actcurve->value.a.da[i][1] = strtod(colpointer,&colpointer);
   }
   else
   {
       colpointer = strstr(allfiles.actplace,"BYABSTIME");
       colpointer += 13;
       actcurve->time.a.da[i][0] = strtod(colpointer,&colpointer);
       actcurve->time.a.da[i][1] = strtod(colpointer,&colpointer);
       colpointer = strstr(allfiles.actplace,"FACTOR");
       colpointer += 6;
       actcurve->value.a.da[i][0] = strtod(colpointer,&colpointer);
       actcurve->value.a.da[i][1] = strtod(colpointer,&colpointer);
   }
   frdouble("T",&(actcurve->T),&ierr);
} /* end of loop over curve lines */
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_read_curve */
