/*!---------------------------------------------------------------------
\file
\brief contains bug and time tracing routines

---------------------------------------------------------------------*/

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
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;                      
/*!----------------------------------------------------------------------
\brief counter of memory in byte

<pre>                                                         m.gee 02/02    
defined in pss_am.c
</pre>

*----------------------------------------------------------------------*/
#ifdef DEBUG
extern long int num_byte_allocated;
#endif

/*! 
\addtogroup DSSYSTEM 
*//*! @{ */

/*!----------------------------------------------------------------------
\brief the tracing variable

<pre>                                                         m.gee 8/00
defined in pss_ds.c, declared in tracing.h                                                  
</pre>
*----------------------------------------------------------------------*/
#ifdef DEBUG
struct _TRACE         trace;
#endif
/*!----------------------------------------------------------------------
\brief global variables for time tracing
<pre>                                                          genk 05/02
</pre>

*----------------------------------------------------------------------*/
#ifdef PARALLEL
double par_start;
#else
time_t seq_start;
#endif

/*!---------------------------------------------------------------------
\brief Initialize bugtracing systems                                              

<pre>                                                        m.gee 8/00 
-sets num_byte_allocted to zero
-checks for type unsigned char to be exactly one byte
 (Needed for ptr-shifting)
-for ARRAY tracing the start of a linear chained list of unknown length
 is allocated and a ptr is set to the last element in this list.
 With every ARRAY defined, another chain element is attached, with an
 ARRAY deleted, the corresponding chain element is also deleted.
-for routine tracing a chained list closed ring is allocated inside a
 linear vector, which ends are connected. These ring has 100 elements and
 is therefore able to trace routine calls upto a deepness of 100 before it 
 starts overriding itself.
 Everytime a subroutine is entered a ptr is set forward (like on a clock)
 in this chain ring and the name of the routine and status 'dsin' is set.
 On exit of this routine, the ptr is set anticlockwis back to the calling
 routine
</pre>
\return void                                                

------------------------------------------------------------------------*/
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
the routine CCACALLOC uses unsigned char to allocate internally, so it is
necessary, that unsigned char is exactly one byte in DEBUG mode 
*/
if (sizeof(unsigned char) != 1)
{
   dserror("unsigned char not 1 byte - will have CCACALLOC problems !!!");
}
/*================================================tracing of arrays=====*/
/*--------------------------------- allocate one initial piece of chain */
trace.arraychain = (TRACEARRAY*)CCACALLOC(1,sizeof(TRACEARRAY));
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
/*----------------------------------------- initialise CPU-time tracing */
ds_cputime_init();
#endif
return;
} /* end of dsinit */


/*!---------------------------------------------------------------------
\brief report entrance to routine                                              

<pre>                                                        m.gee 8/00 
This routine reports the entry to a subroutine the the ds-system
see dsinit()
</pre>
\param string   char[]   (i)  name of routine                                
\return void                                               
\sa dstrc_exit() , dsinit()                                    

------------------------------------------------------------------------*/
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
} /* end of dstrc_enter */

/*!---------------------------------------------------------------------
\brief report exit to routine                                              

<pre>                                                        m.gee 8/00 
This routine reports the exit of a subroutine the the ds-system
see dsinit()
</pre>
\param string   char[]   (i)  name of routine                                
\return void                                               
\sa dstrc_enter() , dsinit()                                    

------------------------------------------------------------------------*/
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
} /* end of dstrc_exit */





/*!---------------------------------------------------------------------
\brief report a new array to the bugtracing system                                              

<pre>                                                        m.gee 8/00 
This routine reports the creation of an ARRAY to the ds-system
It creates a new chain element in the list of ARRAY-tracing structures
this routine is called by the am-system only !
see dsinit()
</pre>
\param string   char[]   (i)  name of routine                                
\param typ      int      (i)  type of array 1 = ARRAY / 2 = ARRAY4D                                
\return void                                               
\sa dsinit()                                    

------------------------------------------------------------------------*/
void dsreportarray(void *array, int typ)
{
#ifdef DEBUG
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
trace.endarraychain->next = (TRACEARRAY*)CCACALLOC(1,sizeof(TRACEARRAY));
if (!trace.endarraychain->next) dserror("Allocation of memory failed");
/*---------------------------- set pointer backwards in this new piece */
trace.endarraychain->next->prev = trace.endarraychain;
/*------------------------------------- set endarraychain to new piece */
trace.endarraychain = trace.endarraychain->next;
#endif
return;
} /* end of dstracereport */


/*!---------------------------------------------------------------------
\brief report a deletion of an array to the bugtracing system                                              

<pre>                                                        m.gee 8/00 
This routine reports the deletion of an ARRAY to the ds-system
It deleted the correponding chain element in the list of ARRAY-tracing 
structures
this routine is called by the am-system only !
see dsinit()
</pre>
\param string   char[]   (i)  name of routine                                
\param typ      int      (i)  type of array 1 = ARRAY / 2 = ARRAY4D                                
\return void                                               
\sa dsinit() ,  dstracereport()                                  

------------------------------------------------------------------------*/
void dsdeletearray(void *array, int typ)
{
#ifdef DEBUG
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
CCAFREE(acttracearray);
/*-----------------------------------------------------------------------*/
#endif
return;
} /* end of dsdeletearray */





/*!---------------------------------------------------------------------
\brief write a report about all arrays                                              

<pre>                                                        m.gee 8/00 
-write a report about all arrays to the .err file
-does nothing if DEBUG is not defined 
-writes a list of all ARRAY and ARRAY4D structure generated
see dsinit()
</pre>
\return void                                               
\sa dsinit() ,  dstracereport() , dsdeletearray()                                 

------------------------------------------------------------------------*/
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
      case cca_DA:
fprintf(allfiles.out_err,"ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: DOUBLE-ARRAY \n",
        i,
        acttracer->a.a2->name,
        acttracer->a.a2->fdim,
        acttracer->a.a2->sdim);
        break;
      case cca_DV:
fprintf(allfiles.out_err,"ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: DOUBLE-VECTOR \n",
        i,
        acttracer->a.a2->name,
        acttracer->a.a2->fdim,
        acttracer->a.a2->sdim);
        break;
      case cca_IA:
fprintf(allfiles.out_err,"ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: INTEGER-ARRAY \n",
        i,
        acttracer->a.a2->name,
        acttracer->a.a2->fdim,
        acttracer->a.a2->sdim);
        break;
      case cca_IV:
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
      case cca_D3:
fprintf(allfiles.out_err,"ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: DOUBLE 3D ARRAY\n",
        i,
        acttracer->a.a4->name,
        acttracer->a.a4->fdim,
        acttracer->a.a4->sdim,
        acttracer->a.a4->tdim,
        acttracer->a.a4->fodim
        );
      break;
      case cca_D4:
fprintf(allfiles.out_err,"ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: DOUBLE 4D ARRAY\n",
        i,
        acttracer->a.a4->name,
        acttracer->a.a4->fdim,
        acttracer->a.a4->sdim,
        acttracer->a.a4->tdim,
        acttracer->a.a4->fodim
        );
      break;
      case cca_I3:
fprintf(allfiles.out_err,"ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: INTEGER 3D ARRAY\n",
        i,
        acttracer->a.a4->name,
        acttracer->a.a4->fdim,
        acttracer->a.a4->sdim,
        acttracer->a.a4->tdim,
        acttracer->a.a4->fodim
        );
      break;
      case cca_I4:
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




/*!---------------------------------------------------------------------
\brief report the amount of actual allocated memory                                              

<pre>                                                        m.gee 2/02 
-write a report about all memory allocated to the
 .err file. Memory has to be allocated using the CCAMALLOC CCACALLOC CCAREALLOC
 and CCAFREE functions
-does nothing if DEBUG is not defined 
see dsinit()
</pre>
\return void                                               
\sa dsinit()  CCAMALLOC() , CCACALLOC() , CCAREALLOC() , CCAFREE()                              

------------------------------------------------------------------------*/
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


/*!---------------------------------------------------------------------
\brief asserts a boolean criterium                                              

<pre>                                                        m.gee 3/02 
this routine does nothing if the boolean criterium is true           
but aborts the programm, if it is not                                
The routine is empty in an optimezed (not DEBUG) compilation and     
can therefor be excessively used to develop a secure code, without   
making it slow when running as fast-exe                              
</pre>
\param true     int     (i)   boolean value                       
\param string   char[]  (i)   error message, if true==0                    
\return void                                                
\sa dserror()                                    

------------------------------------------------------------------------*/
void dsassert(int true, char string[])
{
#ifdef DEBUG 
/*----------------------------------------------------------------------*/
if (true) return;
else
dserror(string);
/*----------------------------------------------------------------------*/
#endif
return;
} /* end of dsassert */



/*!---------------------------------------------------------------------
\brief report an error and stop program                                              

<pre>                                                        m.gee 8/00 
-report an error and stop program                       
-prints error message string to console and *.err       
-prints call tree, if DEBUG was defined                 
-aborts parallel and sequentiell programm               
</pre>
\param string   char[]  (i)   error message to be printed                    
\return void                                                
\sa dsassert()                                    

------------------------------------------------------------------------*/
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


/*!---------------------------------------------------------------------
\brief routine to initialise the cpu - time                                              

<pre>                                                        genk 05/02 
routine to initialise the cpu - time
</pre>
\return void                                                

------------------------------------------------------------------------*/
void ds_cputime_init()
{
/*#ifdef PARALLEL

#else
#ifdef DEBUG
struct timeval tv;
struct timezone tz;
double sec, usec;
#endif
#endif*/


#ifdef DEBUG
dstrc_enter("ds_cputime_init");
#endif

#ifdef PARALLEL
par_start=MPI_Wtime();
/*#else
#ifdef DEBUG  
gettimeofday(&tv, &tz);   
sec = (double)(tv.tv_sec); 
usec = (double)(tv.tv_usec);
seq_start=sec+0.000001*usec;*/
#else
seq_start=time(NULL);
#endif


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
}


/*!---------------------------------------------------------------------
\brief routine to meassure the cpu - time                                              

<pre>                                                        genk 05/02 
routine to meassure the cpu - time
</pre>
\return void                                                

------------------------------------------------------------------------*/
double ds_cputime()
{
#ifndef SUSE73
#ifdef PARALLEL
double par_end;
/*#else
#ifdef DEBUG
double seq_end;
struct timeval tv;
struct timezone tz;
double sec, usec;*/
#else
time_t seq_end;
#endif
double diff;

#ifdef DEBUG
dstrc_enter("ds_cputime");
#endif

#ifdef PARALLEL
par_end=MPI_Wtime();
diff=par_end-par_start;
/*#else
#ifdef DEBUG
gettimeofday(&tv, &tz);   
sec = (double)(tv.tv_sec);
usec = (double)(tv.tv_usec);
seq_end=sec+0.000001*usec;
diff=seq_end-seq_start;*/
#else
seq_end=time(NULL);
diff = difftime(seq_end,seq_start);
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return ((double)(diff)); 
#endif
}


/*! @} (documentation module close)*/
