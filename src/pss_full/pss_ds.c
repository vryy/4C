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
 |                                                       m.gee 02/02    |
 | global variable  num_byte_allocated                                  |
 | long int num_byte_allocated                                          |
 *----------------------------------------------------------------------*/
#ifdef DEBUG
extern long int num_byte_allocated;
#endif




/*----------------------------------------------------------------------*
 | Initialize bugtracing systems                          m.gee 8/00    |
 *----------------------------------------------------------------------*/
void dsinit()
{
#ifdef DEBUG 
int i=0;
/*====================================================tracing of memory */
/* 
   num_byte_allocated is a global variable which can be seen here and in
   pss_am.c, where it is defined
*/   
num_byte_allocated=0;
/*---------------------------------------------- check for memory sizes */
/*
the routine CALLOC uses unsigned char to allocate internally, so it is
necessary, that unsigned char is exactly one byte in DEBUG mode 
*/
if (sizeof(unsigned char) != 1)
{
   dserror("unsigned char not 1 byte - will have CALLOC problems !!!");
}
/*================================================tracing of arrays=====*/
/*--------------------------------- allocate one initial piece of chain */
trace.arraychain = (TRACEARRAY*)CALLOC(1,sizeof(TRACEARRAY));
if (!trace.arraychain) dserror("Allocation of memory failed");
/*----------- set endchain ptr to this initial piece, rest is automatic */
trace.endarraychain = trace.arraychain;
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
 | report a new double array to the bugtracing system     m.gee 8/00    |    
 *----------------------------------------------------------------------*/
void dsreportarray(void *array, int typ)
{
/*--------------------- count total number of active ARRAYs or ARRAY4Ds */
trace.num_arrays++;
/*--------------------------------- switch for ARRAY (1) or ARRAY4D (2) */
switch (typ)
{
case 1:
   /* set pointer to ARRAY and ptr to tracing */
   trace.endarraychain->arraytyp = array_2d;
   trace.endarraychain->a.a2     = (ARRAY*)array;
   ((ARRAY*)array)->mytracer     = trace.endarraychain;
break;
case 2:
   /* set pointer to ARRAY4D and ptr to tracing */
   trace.endarraychain->arraytyp = array_4d;
   trace.endarraychain->a.a4     = (ARRAY4D*)array;
   ((ARRAY4D*)array)->mytracer   = trace.endarraychain;
break;
default:
   dserror("Unknown type of ARRAY to watch");
break;   
}
/*----------------------- allocate a new piece to the end of the chain */
trace.endarraychain->next = (TRACEARRAY*)CALLOC(1,sizeof(TRACEARRAY));
if (!trace.endarraychain->next) dserror("Allocation of memory failed");
/*---------------------------- set pointer backwards in this new piece */
trace.endarraychain->next->prev = trace.endarraychain;
/*------------------------------------- set endarraychain to new piece */
trace.endarraychain = trace.endarraychain->next;
return;
} /* end of dstracereport */


/*----------------------------------------------------------------------*
 | report a new double array to the bugtracing system     m.gee 8/00    |    
 *----------------------------------------------------------------------*/
void dsdeletearray(void *array, int typ)
{
TRACEARRAY *acttracearray;
/*------------------------------ decrease number of ARRAYs and ARRAY4Ds */
trace.num_arrays--;
/*--------------------------- switch for type: ARRAY (1) or ARRAY4D (2) */
switch (typ)
{
case 1:
/* set pointer acttracearray to the piece of chain belonging to this ARRAY */
acttracearray = ((ARRAY*)array)->mytracer;
/*-------------------------------------------- set ptr in ARRAY to NULL */
((ARRAY*)array)->mytracer=NULL;
break;
case 2:
/* set pointer acttracearray to the piece of chain belonging to this ARRAY4D */
acttracearray = ((ARRAY4D*)array)->mytracer;
/*------------------------------------------ set ptr in ARRAY4D to NULL */
((ARRAY4D*)array)->mytracer=NULL;
break;
default:
   dserror("Unknown type of ARRAY to watch");
break;   
}

/* if acttracearray is NOT the first piece in the chain pointed to by trace.arraychain */
if (acttracearray->prev)
{  
   /*-------------- set forward pointer of previous piece to next piece */ 
   acttracearray->prev->next = acttracearray->next;
   /*------------- set backward pointer of next piece to previous piece */
   acttracearray->next->prev = acttracearray->prev;
}
/* if acttracearray IS the first piece in the chain pointed to by trace.arraychain */
else 
{
   /*-------------------------------- set trace.arraychain to next piece */
   trace.arraychain = acttracearray->next;
   /* set the backward pointer of this next piece to NULL to indicate that 
      it is the first piece in the chain pointed to by trace.arraychain  */
   acttracearray->next->prev = NULL;
}   
/*---------- delete the actual piece which now is taken out of the chain */
FREE(acttracearray);
/*-----------------------------------------------------------------------*/
return;
} /* end of dsdeletearray */





/*----------------------------------------------------------------------*
 | write a report about all arrays to the .err file       m.gee 8/00    |    
 *----------------------------------------------------------------------*/
void dstrace_to_err()
{
int         i=0;
TRACEARRAY *acttracer;

#ifdef DEBUG 
if (trace.trace_on==1)
{
acttracer = trace.arraychain;

fprintf(allfiles.out_err,"===========================================\n");;
fprintf(allfiles.out_err,"dstrace - array report from routine %s\n",trace.actroutine->name);
fprintf(allfiles.out_err,"===========================================\n");;
fprintf(allfiles.out_err,"number of actual arrays: %d\n",trace.num_arrays);
while (acttracer->next != NULL)
{
      /* print for a ARRAY */
   if (acttracer->arraytyp == array_2d)
      switch(acttracer->a.a2->Typ)
      {
      case DA:
fprintf(allfiles.out_err,"ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: DOUBLE-ARRAY \n",
        i,
        acttracer->a.a2->name,
        acttracer->a.a2->fdim,
        acttracer->a.a2->sdim);
        break;
      case DV:
fprintf(allfiles.out_err,"ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: DOUBLE-VECTOR \n",
        i,
        acttracer->a.a2->name,
        acttracer->a.a2->fdim,
        acttracer->a.a2->sdim);
        break;
      case IA:
fprintf(allfiles.out_err,"ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: INTEGER-ARRAY \n",
        i,
        acttracer->a.a2->name,
        acttracer->a.a2->fdim,
        acttracer->a.a2->sdim);
        break;
      case IV:
fprintf(allfiles.out_err,"ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: INTEGER-VECTOR \n",
        i,
        acttracer->a.a2->name,
        acttracer->a.a2->fdim,
        acttracer->a.a2->sdim);
        break;
      default:
fprintf(allfiles.out_err,"ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: DAMAGED TYPE !!!!!!!!!!\n",
        i,
        acttracer->a.a2->name,
        acttracer->a.a2->fdim,
        acttracer->a.a2->sdim);
      }
      
      /* print for a ARRAY4D */ 
   if (acttracer->arraytyp == array_4d)
      switch(acttracer->a.a4->Typ)
      {
      case D3:
fprintf(allfiles.out_err,"ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: DOUBLE 3D ARRAY\n",
        i,
        acttracer->a.a4->name,
        acttracer->a.a4->fdim,
        acttracer->a.a4->sdim,
        acttracer->a.a4->tdim,
        acttracer->a.a4->fodim
        );
      break;
      case D4:
fprintf(allfiles.out_err,"ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: DOUBLE 4D ARRAY\n",
        i,
        acttracer->a.a4->name,
        acttracer->a.a4->fdim,
        acttracer->a.a4->sdim,
        acttracer->a.a4->tdim,
        acttracer->a.a4->fodim
        );
      break;
      case I3:
fprintf(allfiles.out_err,"ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: INTEGER 3D ARRAY\n",
        i,
        acttracer->a.a4->name,
        acttracer->a.a4->fdim,
        acttracer->a.a4->sdim,
        acttracer->a.a4->tdim,
        acttracer->a.a4->fodim
        );
      break;
      case I4:
fprintf(allfiles.out_err,"ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: INTEGER 4D ARRAY\n",
        i,
        acttracer->a.a4->name,
        acttracer->a.a4->fdim,
        acttracer->a.a4->sdim,
        acttracer->a.a4->tdim,
        acttracer->a.a4->fodim
        );
      break;
      default:
fprintf(allfiles.out_err,"ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: DAMAGED TYPE !!!!!!!!!!\n",
        i,
        acttracer->a.a4->name,
        acttracer->a.a4->fdim,
        acttracer->a.a4->sdim,
        acttracer->a.a4->tdim,
        acttracer->a.a4->fodim
        );
      break;
      }

   /* set acttracer to next tracer */
   acttracer = acttracer->next;
   i++;
      
} /* end of loop (i=0; i<trace.num_arrays; i++) */
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




/*----------------------------------------------------------------------*
 | report the amount of actual allocated memory           m.gee 2/02    |
 *----------------------------------------------------------------------*/
void dsmemreport()
{
#ifdef DEBUG 
char   *colptr;
char    message[300];
double  mbyte;
/*----------------------------------------------------------------------*/
mbyte = (double)num_byte_allocated;
mbyte /= 1048576.0;

if (trace.trace_on==1)
{
   strcpy(message,"PROC ");
   colptr = message + strlen(message); 
   sprintf(colptr,"%d",par.myrank);
   colptr = message + strlen(message);
   strcpy(colptr," memory used in ");
   colptr = message + strlen(message);
   strcpy(colptr,trace.actroutine->name);
   colptr = message + strlen(message); 
   strcpy(colptr," : ");
   colptr = message + strlen(message);
   sprintf(colptr,"%.6f MegaByte\n",mbyte); 
}
else
{
   strcpy(message,"PROC ");
   colptr = message + strlen(message); 
   sprintf(colptr,"%d",par.myrank);
   colptr = message + strlen(message);
   strcpy(colptr," memory used : ");
   colptr = message + strlen(message);
   sprintf(colptr,"%.6f MegaByte\n",mbyte); 
}
fprintf(allfiles.out_err,"%s",message);
printf (                 "%s",message);
fflush(allfiles.out_err);
/*----------------------------------------------------------------------*/
#endif
return;
} /* end of dsmemreport */



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
