/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
struct _IO_FLAGS        ioflags;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 *----------------------------------------------------------------------*/
extern struct _DYNAMIC *dyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
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

/*---------------------------------------------------------------------*
 | monotoring informations                                  genk 01/03 |
 *---------------------------------------------------------------------*/
extern struct _MONITOR *moni;

/*!---------------------------------------------------------------------
\brief  print out monitoring data 

<pre>                                                         genk 01/03
		     
</pre>

\param *actfield    FIELD   (i)     actual field   
\param  numf        INT     (i)     number of fluid field

\return void                                                                             

------------------------------------------------------------------------*/
void out_monitor(FIELD *actfield, INT numf)
{
INT       i,j;
INT       actpos;
INT       myrank;
INT       numval;
INT       numstep;
INT       numnp;
INT       rest;
INT       full;
INT     **monnodes;
INT     **onoff;
INT       counter;
DOUBLE   *time;
DOUBLE  **val;
MONITOR  *actmoni;
FILE     *mon = allfiles.out_mon;
ARRAY     valdofs_a;
INT      *valdofs;
ARRAY     valIds_a;
INT      *valIds;

#ifdef DEBUG 
dstrc_enter("out_monitor");
#endif
/*----------------------------------------------------------------------*/
myrank = par.myrank;
actmoni = &(moni[numf]);
numstep = actmoni->numstep;
numval = actmoni->numval;
numnp = actmoni->numnp;

if (numstep == 0 || numval == 0) goto end;

if (numval>=3)
{
   rest = numval % 3;
   full = numval/3;
}
else
{
   rest = numval;
   full = 0;
}    

val      = actmoni->val.a.da;
time     = actmoni->time.a.dv;
monnodes = actmoni->monnodes.a.ia;
onoff    = actmoni->onoff.a.ia; 

valdofs = amdef("valdofs",&valdofs_a,numval,1,"IV");
valIds  = amdef("valIds",&valIds_a,numval,1,"IV");
counter=0;
for (i=0;i<numnp;i++)
{
   for (j=0;j<MAXDOFPERNODE;j++)
   {
      if (actmoni->onoff.a.ia[i][j]==-1) continue;
      else 
      {
         valdofs[counter] = j;
	 valIds[counter]=monnodes[i][0];
         counter++;
      }
   }
}

if (counter!=numval)
dserror("Cannot print monitoring data - counter != numval!\n");

/*----------------------------------------------------------------------*/
if (myrank==0)
{

fprintf(mon,"================================================================================\n");
switch (actfield->fieldtyp)
{
case structure:
   fprintf(mon,"MONITORING DATA: STRUCTURE FIELD\n");
break;
case fluid:
   fprintf(mon,"MONITORING DATA: FLUID FIELD\n");
break;
case ale:
   fprintf(mon,"MONITORING DATA: ALE FIELD\n");
break;
default:
break;
}
fprintf(mon,"================================================================================\n");
fprintf(mon,"\n");

actpos=0;
for (i=0;i<full;i++)
{
   fprintf(mon,"================================================================================\n");
   fprintf(mon,"TIME            globalID/DOF           globalID/DOF          globalID/DOF         \n");
   fprintf(mon,"                  %6d/%1d               %6d/%1d              %6d/%1d         \n",
            valIds[actpos],valdofs[actpos],valIds[actpos+1],valdofs[actpos+1],valIds[actpos+2],valdofs[actpos+2]);                  
   fprintf(mon,"================================================================================\n");
   for (j=0;j<numstep;j++)
   {
      fprintf(mon,"%.8lf %18.7E %22.7E %21.7E\n",
      time[j],val[actpos][j],val[actpos+1][j],val[actpos+2][j]);
   }
   fprintf(mon,"________________________________________________________________________________\n\n");
   fprintf(mon,"\n");
   actpos += 3;
}

if (rest==2)
{
   fprintf(mon,"================================================================================\n");
   fprintf(mon,"TIME            globalID/DOF           globalID/DOF                 \n");
   fprintf(mon,"                  %6d/%1d               %6d/%1d                      \n",
           valIds[actpos],valdofs[actpos],valIds[actpos+1],valdofs[actpos+1]);                  
   fprintf(mon,"================================================================================\n");
   for (j=0;j<numstep;j++)
   {
      fprintf(mon,"%.8lf %18.7E %22.7E \n",
      time[j],val[actpos][j],val[actpos+1][j]);
   }
   fprintf(mon,"________________________________________________________________________________\n\n");
   fprintf(mon,"\n");
}
else if (rest==1)
{
   fprintf(mon,"================================================================================\n");
   fprintf(mon,"TIME            globalID/DOF              \n");
   fprintf(mon,"                  %6d/%1d                                \n",
           valIds[actpos],valdofs[actpos]);                  
   fprintf(mon,"================================================================================\n");
   for (j=0;j<numstep;j++)
   {
      fprintf(mon,"%.8lf %18.7E  \n",
      time[j],val[actpos][j]);
   }
   fprintf(mon,"________________________________________________________________________________\n\n");
   fprintf(mon,"\n");

}

amdel(&valdofs_a);
amdel(&valIds_a);
}
if (myrank==0) fflush(mon);
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_monitor */

/*!---------------------------------------------------------------------
\brief  print out area to monitoring file

<pre>                                                         genk 01/03
		     
</pre>

\param *actfield    FIELD   (i)     actual field   
\param  numf        INT     (i)     number of fluid field
\warning this function has to be parallelised 
\return void                                                                             

------------------------------------------------------------------------*/
void out_area(ARRAY totarea_a)
{
INT        i,j;
INT        fdim,sdim;
INT        myrank;
DOUBLE   **totarea;
FILE      *mon = allfiles.out_mon;

#ifdef DEBUG 
dstrc_enter("out_area");
#endif
/*----------------------------------------------------------------------*/
myrank  = par.myrank;
fdim    = totarea_a.fdim;
sdim    = totarea_a.sdim;
totarea = totarea_a.a.da;


/*----------------------------------------------------------------------*/
if (myrank==0)
{
fprintf(mon,"================================================================================\n");
fprintf(mon,"MONITORING DATA: TOTAL AREA FLUID FIELD\n");
fprintf(mon,"================================================================================\n");
fprintf(mon,"\n");

for (i=0;i<fdim;i++)
{
   fprintf(mon,"%6d",i+1);
   for (j=0;j<sdim;j++)
   {
      fprintf(mon,"    %.8lf",totarea[i][j]);
   }
   fprintf(mon,"\n");
}
fprintf(mon,"________________________________________________________________________________\n\n");
fprintf(mon,"\n");
}

if (myrank==0) fflush(mon);

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_monitor */
