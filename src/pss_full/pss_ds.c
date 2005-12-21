/*!---------------------------------------------------------------------
\file
\brief contains bug and time tracing routines

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/
#include <stdarg.h>
#include "../headers/standardtypes.h"
#include "../io/io.h"
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
struct _CCA_TRACE         trace;
#endif
/*!----------------------------------------------------------------------
\brief global variables for time tracing
<pre>                                                          genk 05/02
</pre>

*----------------------------------------------------------------------*/
#ifdef PARALLEL
DOUBLE par_start;
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
INT i=0;
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
if (sizeof(unsigned char) != (unsigned)1)
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
trace.routine[0].name = "xxxxxxxxxxxxxxxxx";
for (i=1; i<99; i++)
{
   trace.routine[i].prev = &(trace.routine[i-1]);
   trace.routine[i].next = &(trace.routine[i+1]);
   trace.routine[i].name = "xxxxxxxxxxxxxxxxx";
   trace.routine[i].dsroutcontrol = dsnone;
}
trace.routine[99].prev = &(trace.routine[98]);
trace.routine[99].next = &(trace.routine[0]);
trace.routine[99].dsroutcontrol = dsnone;
trace.routine[99].name = "xxxxxxxxxxxxxxxxx";

/*------------------------------------------------- set starting values */
trace.deepness=2;

/*-------------------------------------------------------- this routine */
trace.routine[0].name = "main";
trace.routine[0].dsroutcontrol = dsin;
trace.routine[1].name = "ntam";
trace.routine[1].dsroutcontrol = dsin;
trace.routine[2].name = "ntaini";
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
void dstrc_enter(
    char        string[])
{

#ifdef DEBUG
  trace.actroutine = trace.actroutine->next;
  trace.actroutine->name = string;
  trace.actroutine->dsroutcontrol=dsin;
  trace.deepness++;
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
  trace.actroutine->dsroutcontrol=dsout;
  trace.actroutine = trace.actroutine->prev;
  trace.deepness--;
  dsassert(trace.deepness >= 0, "trace stack underflow");
#endif

  return;
} /* end of dstrc_exit */




/*----------------------------------------------------------------------*/
/*!
 \brief print the current function stack.
 */
/*----------------------------------------------------------------------*/
void dstrc_whereami()
{

#ifdef DEBUG
  INT i;
  TRACEROUT *routhis = trace.actroutine;
  for (i=0; i<trace.deepness; i++)
  {
    routhis = routhis->prev;
    fprintf(allfiles.out_err,"%s\n",routhis->name);
    printf("%s\n",routhis->name);
  }
#endif

}


/*!---------------------------------------------------------------------
\brief report a new array to the bugtracing system

<pre>                                                        m.gee 8/00
This routine reports the creation of an ARRAY to the ds-system
It creates a new chain element in the list of ARRAY-tracing structures
this routine is called by the am-system only !
see dsinit()
</pre>
\param string   char[]   (i)  name of routine
\param typ      INT      (i)  type of array 1 = ARRAY / 2 = ARRAY4D
\return void
\sa dsinit()

------------------------------------------------------------------------*/
void dsreportarray(void *array, INT typ)
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
if (!trace.endarraychain->next)
  dserror("Allocation of memory failed");
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
\param typ      INT      (i)  type of array 1 = ARRAY / 2 = ARRAY4D
\return void
\sa dsinit() ,  dstracereport()

------------------------------------------------------------------------*/
void dsdeletearray(void *array, INT typ)
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

#ifdef DEBUG

  INT         i=0;
  TRACEARRAY *acttracer;

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
          fprintf(allfiles.out_err,
              "ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: DOUBLE-ARRAY \n",
              i,
              acttracer->a.a2->name,
              acttracer->a.a2->fdim,
              acttracer->a.a2->sdim);
          break;
        case cca_DV:
          fprintf(allfiles.out_err,
              "ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: DOUBLE-VECTOR \n",
              i,
              acttracer->a.a2->name,
              acttracer->a.a2->fdim,
              acttracer->a.a2->sdim);
          break;
        case cca_IA:
          fprintf(allfiles.out_err,
              "ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: INTEGER-ARRAY \n",
              i,
              acttracer->a.a2->name,
              acttracer->a.a2->fdim,
              acttracer->a.a2->sdim);
          break;
        case cca_IV:
          fprintf(allfiles.out_err,
              "ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: INTEGER-VECTOR \n",
              i,
              acttracer->a.a2->name,
              acttracer->a.a2->fdim,
              acttracer->a.a2->sdim);
          break;
        default:
          fprintf(allfiles.out_err,
              "ARRAY   NO%8d NAME %9s DIM %4d x%4d TYPE: DAMAGED TYPE !!!\n",
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
          fprintf(allfiles.out_err,
              "ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: DOUBLE 3D ARRAY\n",
              i,
              acttracer->a.a4->name,
              acttracer->a.a4->fdim,
              acttracer->a.a4->sdim,
              acttracer->a.a4->tdim,
              acttracer->a.a4->fodim
              );
          break;
        case cca_D4:
          fprintf(allfiles.out_err,
              "ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: DOUBLE 4D ARRAY\n",
              i,
              acttracer->a.a4->name,
              acttracer->a.a4->fdim,
              acttracer->a.a4->sdim,
              acttracer->a.a4->tdim,
              acttracer->a.a4->fodim
              );
          break;
        case cca_I3:
          fprintf(allfiles.out_err,
              "ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: INTEGER 3D ARRAY\n",
              i,
              acttracer->a.a4->name,
              acttracer->a.a4->fdim,
              acttracer->a.a4->sdim,
              acttracer->a.a4->tdim,
              acttracer->a.a4->fodim
              );
          break;
        case cca_I4:
          fprintf(allfiles.out_err,
              "ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: INTEGER 4D ARRAY\n",
              i,
              acttracer->a.a4->name,
              acttracer->a.a4->fdim,
              acttracer->a.a4->sdim,
              acttracer->a.a4->tdim,
              acttracer->a.a4->fodim
              );
          break;
        default:
          fprintf(allfiles.out_err,
              "ARRAY4D NO%8d NAME %9s DIM %4d x%4d x%4d x%4d TYPE: DAMAGED TYPE !!!\n",
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

#else

  fprintf(allfiles.out_err,"bugtracing only in DEBUG - noreport\n");

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
  DOUBLE  mbyte;

  mbyte = (DOUBLE)num_byte_allocated;
  mbyte /= 1048576.0;

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
  fprintf(allfiles.out_err,"%s",message);
  printf (                 "%s",message);
  fflush(allfiles.out_err);

#endif

  return;
} /* end of dsmemreport */




/*
  The latest file position.
 */
static INT latest_line;
static char* latest_file = "{dserror_func call without prototype}";




/*!---------------------------------------------------------------------
\brief asserts a boolean criterium

<pre>                                                        m.gee 3/02
this routine does nothing if the boolean criterium is true
but aborts the programm, if it is not
The routine is empty in an optimezed (not DEBUG) compilation and
can therefor be excessively used to develop a secure code, without
making it slow when running as fast-exe
</pre>
\param true     INT     (i)   boolean value
\param string   char[]  (i)   error message, if true==0
\return void
\sa dserror()

------------------------------------------------------------------------*/
void dsassert_func(
    char*               file,
    INT                 line,
    INT                 test,
    char                string[]
    )
{

#ifdef DEBUG

  if (!test)
  {
    latest_file = file;
    latest_line = line;
    dserror_func(string);
  }

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
void dserror_func(
    char               *string, ...
    )

{

  va_list ap;
  char line[] = "=========================================================================\n";

#ifdef DEBUG
  INT i=0;
  TRACEROUT *routhis = NULL;
#endif


  va_start(ap, string);

#ifdef DEBUG

  printf("\n");
  printf("\n");
  printf(line);
  printf("PROC %d ERROR in %s, line %i:\n",par.myrank,latest_file,latest_line);
  vprintf(string,ap);
  printf("\n");
  printf("\n");

  fprintf(allfiles.out_err,"\n");
  fprintf(allfiles.out_err,"\n");
  fprintf(allfiles.out_err,line);
  fprintf(allfiles.out_err,"PROC %d ERROR in %s, line %i:\n",
      par.myrank,latest_file,latest_line);
  vfprintf(allfiles.out_err,string,ap);
  fprintf(allfiles.out_err,"\n");
  fprintf(allfiles.out_err,"\n");

  routhis = trace.actroutine;

  printf("Calling stack:\n");
  printf("--------------\n");
  printf("%s  <== ERROR is in here!!\n",routhis->name);

  fprintf(allfiles.out_err,"Calling stack:\n");
  fprintf(allfiles.out_err,"--------------\n");
  fprintf(allfiles.out_err,"%s  <== ERROR is in here!!\n",routhis->name);

  for (i=0; i<trace.deepness; i++)
  {
    routhis = routhis->prev;
    printf("%s\n",routhis->name);
    fprintf(allfiles.out_err,"%s\n",routhis->name);
  }

  printf(line);
  printf("\n");

  fprintf(allfiles.out_err,line);
  fprintf(allfiles.out_err,"\n");

#else

  fprintf(allfiles.out_err,"\n");
  fprintf(allfiles.out_err,"\n");
  fprintf(allfiles.out_err,line);
  fprintf(allfiles.out_err,"PROC %d ERROR in %s, line %i:\n",
      par.myrank,latest_file,latest_line);
  vfprintf(allfiles.out_err,string,ap);
  fprintf(allfiles.out_err,"\n");
  fprintf(allfiles.out_err,line);
  fprintf(allfiles.out_err,"\n");
  fprintf(allfiles.out_err,"\n");

  printf("\n");
  printf("\n");
  printf(line);
  printf("PROC %d ERROR in %s, line %i:\n",par.myrank,latest_file,latest_line);
  vprintf(string,ap);
  printf("\n");
  printf(line);
  printf("\n");
  printf("\n");

#endif

  va_end(ap);

  fflush(stdout);
  fflush(allfiles.out_err);

#ifdef BINIO
  io_emergency_close_files();
#endif

#ifdef DSERROR_DUMP
  *((INT*)0x0) = 123456;
#endif


#ifdef PARALLEL
  MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
  exit(EXIT_FAILURE);
#endif

  return;
} /* end of dserror_func */



/*!
  \brief set a file name and line number to be used by the next
  dserror_func

  \return pointer to the dserror_func that uses the latest file
  position info.
 */
void dslatest(char* file, INT line)
{
  latest_file = file;
  latest_line = line;
}





/*!-----------------------------------------------------------------------
\brief get warnings and writes them to the screen at the end

<pre>                                                             ck 07/03
depending on task it initialises all warnings to Zero (= no warning),
collects warnings during running process and reports them to the screen at
the end. A reoccuring warning is written once (per proc).
In parallel processes the warnings are ordered to the fields where they
occure. Warnings are written to the screen and the error files.

task 	= 0 	initialisation to no warnings
	= 1	create warning depending on 'warning'
	= 2	plot warnings to the screen
warning	= 1	when trinangular ale elements are monitored with wrong
		quality meassure
	= 2	aspect ration quality meassure for ale elements is called
		for a distyp it has not been implemented for
	= 3	corner angle quality criterion for ale elements is called
		for a distyp it	has not been implemented for
	= 4	min J quality criterion for ale elements is called for a
		distyp it has not been implemented for
        = 5     Function used for points outside the range 0 < xi < 1
        = 6     Function used for points not on the defined line
        = 7     the value givne for MAXNOD is too large
        = 8     the value givne for MAXELE is too large
        = 9     the value givne for MAXDOFPERNODE is too large
        = 10    the value givne for MAXGAUSS is too large
        = 11    ??
        = 12    ??
        = 13    more than one liftdrag Id to one element found
        = 14    your new warning issue

</pre>

\param task     INT (i)  flag, what to do
\param warning  INT (i)  flag, which warning shall be printed at the end
\return void
\sa

------------------------------------------------------------------------*/
void dswarning(INT task, INT warning)
{
#ifdef PARALLEL
INT     i;                                /* a counter                    */
INTRA  *actintra;
INT     recvbuf;
#endif
static INT called;                      /* flag, if warnings occured    */
static INT ale_quality_min_J_triangles;
static INT ale_quality_ar;
static INT ale_quality_ca;
static INT ale_quality_Je;
static INT funct_range;     /* warning 5 */
static INT funct_line;      /* warning 6 */
static INT dirich_fsi_line; /* warning 7 */
static INT max_nod;         /* warning 8 */
static INT max_ele;         /* warning 9 */
static INT max_dofpernode;  /* warning 10 */
static INT max_gauss;       /* warning 11 */
static INT rescheck_nonode; /* warning 12 */
static INT liftdragId;      /* warning 13 */
/* DEFINE your new warning flag here!!! */

FILE   *err = allfiles.out_err;

#ifdef DEBUG
dstrc_enter("dswarning");
#endif

/*----------------------------------------------------------------------*/
switch (task)
{
  /*------------------------------------------------- initialisation ---*/
  case 0:
    called = 0;
    ale_quality_min_J_triangles = 0;
    ale_quality_ar = 0;
    ale_quality_ca = 0;
    ale_quality_Je = 0;
    funct_range = 0;
    funct_line  = 0;
    dirich_fsi_line=0;
    max_nod        = 0;
    max_ele        = 0 ;
    max_dofpernode = 0;
    max_gauss      = 0;
    rescheck_nonode= 0;
    liftdragId     = 0;
    /* INITIALISE your new warning here!!! */
  break;
  /*------------------------------------------------- create warning ---*/
  case 1:
    called++;
    if (warning == 1)
      ale_quality_min_J_triangles++;
    else if (warning == 2)
      ale_quality_ar++;
    else if (warning == 3)
      ale_quality_ca++;
    else if (warning == 4)
      ale_quality_Je++;
    else if (warning == 5)
      funct_range++;
    else if (warning == 6)
      funct_line++;
    else if (warning == 7)
      dirich_fsi_line++;
    else if (warning == 8)
      max_nod++;
    else if (warning == 9)
      max_ele++;
    else if (warning ==10)
      max_dofpernode++;
    else if (warning ==11)
      max_gauss++;
    else if (warning ==12)
      rescheck_nonode++;
    else if (warning ==13)
      liftdragId++;
    else
      dserror("warning out of range!\n");
    /* COUNT your new warning here!!! */
  break;
  /*------------------------------------------------- write warnings ---*/
  case 2:
  #ifdef PARALLEL
  MPI_Allreduce(&called,&recvbuf,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  called=recvbuf;
  #endif
  if (!called) goto end;  /* supress output if there were no warnings   */
  #ifdef PARALLEL
    for (i=0; i<par.numfld; i++) /* loop all intra communicators        */
    {
      actintra = &(par.intra[i]);
      MPI_Barrier(actintra->MPI_INTRA_COMM);
      /* output result depends on fieldtyp (in parallel case only) */
      switch (actintra->intra_fieldtyp)
      {
      case fluid:
        if (par.myrank==0) printf("\nWARNINGS Fluid Field: \n");
        fprintf(err,"*************************************\n");
        fprintf(err,"WARNINGS Fluid Field: \n");
      break;
      case ale:
        if (par.myrank==0) printf("\nWARNINGS ALE Field: \n");
        fprintf(err,"*************************************\n");
        fprintf(err,"WARNINGS ALE Field: \n");
      break;
      case structure:
        if (par.myrank==0) printf("\nWARNINGS Structure Field: \n");
        fprintf(err,"*************************************\n");
        fprintf(err,"WARNINGS Structure Field: \n");
      break;
      default: dserror("Unknown fieldtyp!");
      }
  #endif
   /*-------------------------------------------------------------------*/
      if (ale_quality_min_J_triangles)
      {
        printf("Warning PROC %i: There is no sense in monitoring lineare triangles with min_J!!!\n                The .plt-file is not useful!\n", par.myrank);
        fprintf(err,"Warning PROC %i: There is no sense in monitoring lineare triangles with min_J!!!\n", par.myrank);
        fprintf(err,"                The .plt-file is not useful!\n");
      }
      if (ale_quality_ar)
      {
        printf("Warning PROC %i: aspect_ratio quality monitoring for one of the distyps not implemented!!\n", par.myrank);
        fprintf(err,"Warning PROC %i: aspect_ratio quality monitoring for one of the distyps not implemented!!\n", par.myrank);
      }
      if (ale_quality_ca)
      {
        printf("Warning PROC %i: corner_angle quality monitoring for one of the distyps not implemented!!\n", par.myrank);
        fprintf(err,"Warning PROC %i: corner_angle quality monitoring for one of the distyps not implemented!!\n", par.myrank);
      }
      if (ale_quality_Je)
      {
        printf("\nWarning PROC %i: min_J quality monitoring for one of the distyps not implemented!!\n", par.myrank);
        fprintf(err,"\n Warning PROC %i: min_J quality monitoring for one of the distyps not implemented!!\n", par.myrank);
      }
      if (funct_range)
      {
        printf("\n WARNING PROC %i: Function used for points outside the range 0<xi<1 !!\n", par.myrank);
        fprintf(err,"\n WARNING PROC %i: Function used for points outside the range 0<xi<1 !!\n", par.myrank);
      }
      if (funct_line)
      {
        printf("\n WARNING PROC %i: Function used for points not on the defined line!!\n", par.myrank);
        fprintf(err,"\n WARNING PROC %i: Function used for points not on the defined line!!\n", par.myrank);
      }
      if (dirich_fsi_line)
      {
        printf("\n %d WARNINGS PROC %i: dirich- and fsi-coupling condition defined on same DLINE\n", dirich_fsi_line,par.myrank);
        fprintf(err,"\n %d WARNINGS PROC %i: dirich- and fsi-coupling condition defined on same DLINE\n", dirich_fsi_line,par.myrank);
      }
      if (max_nod)
      {
        printf("\n WARNING PROC %i: The given value for MAXNOD was too large!! This will reduce the performance!!\n", par.myrank);
        fprintf(err,"\n WARNING PROC %i: The given value for MAXNOD was too large!! This will reduce the performance!!\n", par.myrank);
      }

      if (max_ele)
      {
        printf("\n WARNING PROC %i: The given value for MAXELE was too large!! This will reduce the performance!!\n", par.myrank);
        fprintf(err,"\n WARNING PROC %i: The given value for MAXELE was too large!! This will reduce the performance!!\n", par.myrank);
      }

      if (max_dofpernode)
      {
        printf("\n WARNING PROC %i: The given value for MAXDOFPERNODE was too large!! This will reduce the performance!!\n", par.myrank);
        fprintf(err,"\n WARNING PROC %i: The given value for MAXDOFPERNODE was too large!! This will reduce the performance!!\n", par.myrank);
      }

      if (max_gauss)
      {
        printf("\n WARNING PROC %i: The given value for MAXGAUSS was too large!! This will reduce the performance!!\n", par.myrank);
        fprintf(err,"\n WARNING PROC %i: The given value for MAXGAUSS was too large!! This will reduce the performance!!\n", par.myrank);
      }

      if (rescheck_nonode)
      {
        printf("\n %d WARNINGS PROC %i: The given node for result check was not found in actual partition.\n", rescheck_nonode,par.myrank);
        fprintf(err,"\n %d WARNINGS PROC %i: The given node for result check was not found in actual partition.\n", rescheck_nonode,par.myrank);
      }
      if (liftdragId)
      {
        printf("\n %d WARNINGS PROC %i: More than one liftdrag Id to one node. Liftdrag forces may be incorrect.\n", rescheck_nonode,par.myrank);
        fprintf(err,"\n %d WARNINGS PROC %i: More than one liftdrag Id to one node. Liftdrag forces may be incorrect.\n", rescheck_nonode,par.myrank);
      }
      /* ADD more warning messages here!!! */
  #ifdef PARALLEL
    MPI_Barrier(actintra->MPI_INTRA_COMM);
    } /* end loop over intra groups */
  #endif
   if (par.myrank == 0) printf("\n");
  break;
  /*-------------------------------------------------------- default ---*/
  default:
     printf("in default von dswarning\n");
     dserror("warning task %2i invalid",task);
  break;
} /* end of switch init */
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of dswarning */



/*!---------------------------------------------------------------------
\brief routine to initialise the cpu - time

<pre>                                                        genk 05/02
routine to initialise the cpu - time
</pre>
\return void

------------------------------------------------------------------------*/
void ds_cputime_init()
{

#ifdef DEBUG
  dstrc_enter("ds_cputime_init");
#endif

#ifdef PARALLEL
  par_start=MPI_Wtime();
#endif


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
DOUBLE ds_cputime()
{

#ifdef PARALLEL
  DOUBLE par_end;
#endif

  DOUBLE diff;

#ifdef DEBUG
  dstrc_enter("ds_cputime");
#endif

#ifdef PARALLEL
  par_end = MPI_Wtime();
  diff    = par_end-par_start;
#else
  diff    = perf_time();
#endif

#ifdef DEBUG
  dstrc_exit();
#endif

  return ((DOUBLE)(diff));
}

/*! @} (documentation module close)*/
