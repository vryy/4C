#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | tracing variables                                                    |
 | defined in pss_ds.c                                                  |
 *----------------------------------------------------------------------*/
#ifdef DEBUG
struct _TRACE         trace;
#endif
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
 extern struct _PAR   par;                      

/*----------------------------------------------------------------------*
 | report an error and stop program                       m.gee 8/00    |
 *----------------------------------------------------------------------*/
void dserror(char string[])
{

int i=0;
char message[300];
char *colptr=NULL;
TRACEROUT *routhis = NULL;
#ifdef DEBUG 
if (trace.trace_on==1)
{
strcpy(message,"PROC ");
colptr = message + strlen(message); 
sprintf(colptr,"%d",par.myrank);
colptr = message + strlen(message);
strcpy(colptr," ERROR IN ");
colptr = message + strlen(message);
strcpy(colptr,trace.actroutine->name);
colptr = message + strlen(message); 
strcpy(colptr," : ");
colptr = message + strlen(message); 
strcpy(colptr,string);

fprintf(allfiles.out_err,"================================================================\n");
fprintf(allfiles.out_err,"%s\n",message); 

printf("=========================================================================\n");
printf("%s\n",message); 



fprintf(allfiles.out_err,"This routine was called by:\n");
printf("This routine was called by:\n");
routhis = trace.actroutine;
for (i=0; i<trace.deepness; i++)
{
routhis = routhis->prev;
fprintf(allfiles.out_err,"%s\n",routhis->name);
printf("%s\n",routhis->name);
}





fprintf(allfiles.out_err,"================================================================\n");
printf("=========================================================================\n");
}
else
{
#endif
fprintf(allfiles.out_err,"================================================================\n");
fprintf(allfiles.out_err,"%s\n",string); 
printf("=========================================================================\n");
printf("%s\n",string); 
fprintf(allfiles.out_err,"================================================================\n");
printf("=========================================================================\n");
#ifdef DEBUG 
}
#endif
fflush(stdout);
fflush(allfiles.out_err);

#ifdef PARALLEL 
MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
exit(EXIT_FAILURE);
#endif

return;
} /* end of dserror */



/*----------------------------------------------------------------------*
 | Initialize bugtracing systems                          m.gee 8/00    |
 *----------------------------------------------------------------------*/
void dsinit()
{
#ifdef DEBUG 
int i=0;
/*================================================tracing of arrays=====*/
trace.num_arrays=0;
/*----------- start with the possibility to trace 1000 int and doubles; */
trace.size_arrays=1000;
/*------------------------allocate space for the pointers to the arrays */
trace.arrays    = (ARRAY**)calloc(trace.size_arrays,sizeof(ARRAY*));
if (trace.arrays == NULL) dserror("Allocation of memory failed\n");
/*-------------------------------------------------------init with NULL */
for (i=0; i<trace.size_arrays; i++) trace.arrays[i]=NULL;
/*==================================================tracing of routines */
/*------------------------------------------------------- init the ring */
trace.routine[0].prev = &(trace.routine[99]);
trace.routine[0].next = &(trace.routine[1]);
trace.routine[0].dsroutcontrol = dsnone;
strcpy(trace.routine[0].name,"xxxxxxxxxxxxxxxxx");
for (i=1; i<99; i++)
{
   trace.routine[i].prev = &(trace.routine[i-1]);
   trace.routine[i].next = &(trace.routine[i+1]);
   strcpy(trace.routine[i].name,"xxxxxxxxxxxxxxxxx");
   trace.routine[i].dsroutcontrol = dsnone;
}
trace.routine[99].prev = &(trace.routine[98]);
trace.routine[99].next = &(trace.routine[0]);
trace.routine[99].dsroutcontrol = dsnone;
strcpy(trace.routine[99].name,"xxxxxxxxxxxxxxxxx");
/*------------------------------------------------- set starting values */
trace.deepness=2;
/*-------------------------------------------------------- this routine */
strcpy(trace.routine[0].name,"main");
trace.routine[0].dsroutcontrol = dsin;
strcpy(trace.routine[1].name,"ntam");
trace.routine[1].dsroutcontrol = dsin;
strcpy(trace.routine[2].name,"ntaini");
trace.routine[2].dsroutcontrol = dsin;

trace.actroutine = &(trace.routine[2]);

/*======================================================tracing of time */
/* nothing implemented yet */

#endif
return;
} /* end of dsinit */


/*----------------------------------------------------------------------*
 | report entrance to routine to the tracing system       m.gee 8/00    |
 *----------------------------------------------------------------------*/
void dstrc_enter(char string[])
{
#ifdef DEBUG 
if (trace.trace_on==1)
{
trace.actroutine = trace.actroutine->next;

strncpy(trace.actroutine->name,string,49);

trace.actroutine->dsroutcontrol=dsin;

trace.deepness++;
}
#endif
return;
} /* end of dstracesize */





/*----------------------------------------------------------------------*
 | report exit to routine to the tracing system           m.gee 8/00    |
 *----------------------------------------------------------------------*/
void dstrc_exit()
{
#ifdef DEBUG 
if (trace.trace_on==1)
{
trace.actroutine->dsroutcontrol=dsout;
trace.actroutine = trace.actroutine->prev;
trace.deepness--;
}
return;
#endif
} /* end of dstracesize */



/*----------------------------------------------------------------------*
 | Check whether current size of tracing is sufficient    m.gee 8/00    |
 | if not, then realloc the tracing system                              |
 *----------------------------------------------------------------------*/
void dstracesize()
{
#ifdef DEBUG 
int i=0;
dstrc_enter("dstracesize");
   if (trace.num_arrays+1 < trace.size_arrays) /* size is large enough */
   {
   }
   else /*----------------------------------tracing is enlarged in steps of 200 */
   {
      trace.arrays=realloc(trace.arrays,(trace.size_arrays+1000)*sizeof(ARRAY*));
      if (trace.arrays==NULL) dserror("Enlargement of tracing failed");
      for (i=trace.size_arrays; i<trace.size_arrays+1000; i++)
      {
         trace.arrays[i]=NULL;
      }
      trace.size_arrays+=1000;
   }


dstrc_exit();
#endif
return;
} /* end of dstracesize */



/*----------------------------------------------------------------------*
 | report a new double array to the bugtracing system     m.gee 8/00    |    
 *----------------------------------------------------------------------*/
void dstracereport(ARRAY *array)
{
#ifdef DEBUG 
int i=0;
dstracesize();
trace.num_arrays++;
while ( trace.arrays[i] != NULL ) i++;
if (i>=trace.size_arrays || trace.num_arrays>=trace.size_arrays)
dserror("Enlargement of ARRAY tracing failed");
trace.arrays[i] = array;
array->place_in_trace=i;
#endif
return;
} /* end of dstracereport */




/*----------------------------------------------------------------------*
 | write a report about all arrays to the .err file       m.gee 8/00    |    
 *----------------------------------------------------------------------*/
void dstrace_to_err()
{
int i=0;
long int byte=0;

#ifdef DEBUG 
if (trace.trace_on==1)
{
fprintf(allfiles.out_err,"===========================================\n");;
fprintf(allfiles.out_err,"dstrace - array report from routine %s\n",trace.actroutine->name);
fprintf(allfiles.out_err,"===========================================\n");;
fprintf(allfiles.out_err,"number of actual arrays: %d\n",trace.num_arrays);
for (i=0; i<trace.size_arrays; i++)
{
   if (trace.arrays[i] != NULL)
   {
      switch(trace.arrays[i]->Typ)
      {
      case DA:
fprintf(allfiles.out_err,"ARRAY NO%4d NAME %9s DIM %4d x%4d TYPE: DOUBLE-ARRAY \n",
        i,
        trace.arrays[i]->name,
        trace.arrays[i]->fdim,
        trace.arrays[i]->sdim);
        byte += trace.arrays[i]->fdim * trace.arrays[i]->sdim * sizeof(double);
        break;
      case DV:
fprintf(allfiles.out_err,"ARRAY NO%4d NAME %9s DIM %4d x%4d TYPE: DOUBLE-VECTOR \n",
        i,
        trace.arrays[i]->name,
        trace.arrays[i]->fdim,
        trace.arrays[i]->sdim);
        byte += trace.arrays[i]->fdim * trace.arrays[i]->sdim * sizeof(double);
        break;
      case IA:
fprintf(allfiles.out_err,"ARRAY NO%4d NAME %9s DIM %4d x%4d TYPE: INTEGER-ARRAY \n",
        i,
        trace.arrays[i]->name,
        trace.arrays[i]->fdim,
        trace.arrays[i]->sdim);
        byte += trace.arrays[i]->fdim * trace.arrays[i]->sdim * sizeof(int);
        break;
      case IV:
fprintf(allfiles.out_err,"ARRAY NO%4d NAME %9s DIM %4d x%4d TYPE: INTEGER-VECTOR \n",
        i,
        trace.arrays[i]->name,
        trace.arrays[i]->fdim,
        trace.arrays[i]->sdim);
        byte += trace.arrays[i]->fdim * trace.arrays[i]->sdim * sizeof(int);
        break;
      default:
fprintf(allfiles.out_err,"ARRAY NO%4d NAME %9s DIM %4d x%4d TYPE: DAMAGED TYPE !!!!!!!!!!\n",
        i,
        trace.arrays[i]->name,
        trace.arrays[i]->fdim,
        trace.arrays[i]->sdim);
      }
   }
} /* end of loop (i=0; i<trace.size_arrays; i++) */
fprintf(allfiles.out_err,"Total number of bytes in ARRAYs defined at the moment: %d\n",byte);
fprintf(allfiles.out_err,"===========================================\n");;
fprintf(allfiles.out_err,"dstrace - array report END                 \n");
fprintf(allfiles.out_err,"===========================================\n");;
fflush(allfiles.out_err);
}
else
{
#endif
fprintf(allfiles.out_err,"bugtracing was switched off - noreport\n");
#ifdef DEBUG 
}
#endif
return;
} /* end of dstrace_to_err */
