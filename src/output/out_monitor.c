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
\param  numf        INT     (i)     number of field
\param  time        DOUBLE  (i)     the actual time
\param  init        INT     (i)     flag

\return void

------------------------------------------------------------------------*/
void out_monitor(FIELD *actfield, INT numf,DOUBLE time,INT init)
{
#ifndef NO_TEXT_OUTPUT
INT       i,j;
INT       myrank;
INT       numval;
INT       numnp;
INT     **monnodes;
INT     **onoff;
INT       counter;
DOUBLE   *val;
ARRAY     valdofs_a;
INT      *valdofs;
ARRAY     valIds_a;
INT      *valIds;
MONITOR  *actmoni;
FILE     *mon = NULL;

#ifdef DEBUG
dstrc_enter("out_monitor");
#endif

/*----------------------------------------------------------------------*/
if (moni==NULL) goto end;
myrank = par.myrank;
actmoni = &(moni[numf]);
if (actmoni==NULL) goto end;

numval = actmoni->numval;
numnp = actmoni->numnp;

if (numnp==0) goto end;

val      = actmoni->val.a.dv;
monnodes = actmoni->monnodes.a.ia;
onoff    = actmoni->onoff.a.ia;

if (myrank>0) goto end;

/*----------------------------------------------------------------------*/
if (init==1)
{

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

   switch (actfield->fieldtyp)
   {
   case structure:
      mon = allfiles.out_smoni;
      dsassert(mon!=NULL,"cannot write to NULL pointer!\n");
      fprintf(mon,"#================================================================================\n");
      fprintf(mon,"# MONITORING DATA: STRUCTURE FIELD\n");
   break;
   case fluid:
      mon = allfiles.out_fmoni;
      dsassert(mon!=NULL,"cannot write to NULL pointer!\n");
      fprintf(mon,"#================================================================================\n");
      fprintf(mon,"# MONITORING DATA: FLUID FIELD\n");
   break;
   case ale:
      mon = allfiles.out_amoni;
      dsassert(mon!=NULL,"cannot write to NULL pointer!\n");
      fprintf(mon,"# =================================================================================\n");
      fprintf(mon,"# MONITORING DATA: ALE FIELD\n");
   break;
   default:
      dserror("fieldtyp unknown!\n");
   }
   fprintf(mon,"#================================================================================\n");
   fprintf(mon,"#\n");
   fprintf(mon,"#==============");
   for (i=0;i<numval;i++) fprintf(mon,"=========================");
   fprintf(mon,"\n");
   fprintf(mon,"# TIME");
   for (i=0;i<numval;i++) fprintf(mon,"           globalID/DOF  ");
   fprintf(mon,"          \n");
   fprintf(mon,"# ");
   for (i=0;i<numval;i++) fprintf(mon,"                 %6d/%1d",valIds[i],valdofs[i]);
   fprintf(mon,"\n");
   fprintf(mon,"#==============");
   for (i=0;i<numval;i++) fprintf(mon,"=========================");
   fprintf(mon,"\n");

   amdel(&valdofs_a);
   amdel(&valIds_a);
   goto end;
}

switch (actfield->fieldtyp)
{
case structure:
   mon = allfiles.out_smoni;
break;
case fluid:
   mon = allfiles.out_fmoni;
break;
case ale:
   mon = allfiles.out_amoni;
break;
default:
break;
}
dsassert(mon!=NULL,"cannot write to NULL pointer!\n");
fprintf(mon,"  %10.6lf   ",time);
for (i=0;i<numval;i++) fprintf(mon," %- 24.7E",val[i]);
fprintf(mon,"\n");
fflush(mon);

end:
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
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
void out_area(ARRAY totarea_a, DOUBLE time, INT itnum, INT init)
{
#ifndef NO_TEXT_OUTPUT
INT        i;
INT        myrank;
DOUBLE    *totarea;
char      *charpointer;
FILE      *mon;

#ifdef DEBUG
dstrc_enter("out_area");
#endif
/*----------------------------------------------------------------------*/
myrank  = par.myrank;
totarea = totarea_a.a.dv;
if (myrank>0) goto end;

if (init==1)
{
   charpointer=allfiles.outputfile_name+strlen(allfiles.outputfile_kenner);
   sprintf(charpointer,"%d",par.myrank);
   charpointer++;
   strncpy(charpointer,".fluidarea.mon",14);
   charpointer+=14;
   sprintf(charpointer,"\0");
   if ( (allfiles.out_monarea=fopen(allfiles.outputfile_name,"w"))==NULL)
   {
      printf("Opening of output file .fluidarea.mon failed\n");
#ifdef PARALLEL
      MPI_Finalize();
#endif
      exit(1);
   }
   mon = allfiles.out_monarea;
   fprintf(mon,"================================================================================\n");
   fprintf(mon,"MONITORING DATA: TOTAL AREA FLUID FIELD\n");
   fprintf(mon,"================================================================================\n");
   fprintf(mon,"\n");
   goto end;
}

/*----------------------------------------------------------------------*/
mon = allfiles.out_monarea;

fprintf(mon,"%10.6lf   ",time);
for (i=0;i<itnum;i++) fprintf(mon," %24.7lf",totarea[i]);
fprintf(mon,"\n");

fflush(mon);

end:

#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_monitor */
