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
#include "../io/io.h"

#ifdef CCADISCRET
#include "../drt_lib/drt_init_control.H"
#endif

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
struct _FIELD         *field;



#ifdef D_MLSTRUCT
/*----------------------------------------------------------------------*
 |                                                                      |
 | vector of numfld submesh-FIELDs, defined in global_control.c         |
 *----------------------------------------------------------------------*/
struct _FIELD         *sm_field;
#endif /* D_MLSTRUCT */



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 *----------------------------------------------------------------------*/
struct _GENPROB       genprob;


/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
struct _PARTITION    *partition;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
struct _DESIGN       *design;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | struct _ALLDYNA       *alldyn;                                       |
 *----------------------------------------------------------------------*/
ALLDYNA             *alldyn;


/*----------------------------------------------------------------------*
 |                                                          al 08/02    |
 | pointer to allocate eigensolution variables                          |
 | dedfined in global_control.c                                         |
 | struct _ALLEIG       *alleig;                                        |
 *----------------------------------------------------------------------*/
ALLEIG              *alleig;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
struct _STATIC_VAR   *statvar;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
struct _MATERIAL      *mat;
struct _MULTIMAT      *multimat;


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


/*----------------------------------------------------------------------*
 | main routine                                           m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ntam(
    INT                 argc,
    char               *argv[]
    )

{

  DOUBLE   t0,ti,tc;



  /* init all time counters */
#ifdef PERF
  perf_init_all();
  perf_begin(0);
#endif



  /* init devices, tracing, etc...*/
#ifndef CCADISCRET
  ntaini(argc,argv);
#else
  ntaini_ccadiscret(argc,argv);
#endif

  /* input phase, input of all information */

  t0=ds_cputime();
#ifdef PERF
  perf_begin(1);
#endif

#ifndef CCADISCRET
  ntainp();
#else
  ntainp_ccadiscret();
#endif

#ifdef PERF
  perf_end(1);
#endif
  ti=ds_cputime()-t0;
  if (par.myrank==0)
  {
    printf("\n");
    printf("Total CPU Time for INPUT:       %10.3E sec \n",ti);
    printf("\n");
  }


#ifndef CCADISCRET
  /* close the input file, delete the copy */
  /* -> You cannot read anything from input file beyond this point */
  frend();
#endif


  /*-------------------------------------- check for visualisation mode */
  if (genprob.visual!=0) goto visualisation;
  /*------------------------------- set up field-specific communicators */
#ifdef PARALLEL
  create_communicators();
#endif

#ifdef BINIO

  /* The functions to initialize binary io are quite complex. They make
   * use of some standard ccarat facilities. In particular the error
   * logs must already be available. The parallel version makes use of
   * those communicators as well. */

  if (genprob.restart) {
    /* If we want to restart this we'll have to initialize our binary
     * input. That is we need to read the control file. */
    init_bin_in_main(allfiles.outputfile_kenner);
  }

  /* initialize module for binary output */
  /* This can only be done after the (normal) input is read because
   * here some values are already written. */
  init_bin_out_main(allfiles.outputfile_kenner);

#endif

  /*--------------------------------------------------calculation phase */
  t0=ds_cputime();
#ifdef PERF
  perf_begin(2);
#endif

  ntacal();

#ifdef PERF
  perf_end(2);
#endif
  tc=ds_cputime()-t0;
  if (par.myrank==0)
  {
    printf("\n");
    printf("Total CPU Time for CALCULATION: %10.3E sec \n",tc);
    printf("\n");
  }

#ifdef PERF
  perf_end(0);
  /* print out time counters */
  perf_out();
#endif


  goto endcal;
  /*----------------------------------------------------------------------*/
visualisation:
  if (genprob.visual==2 || genprob.visual==3)
    ntavisual();
  /*----------------------------------------------------------------------*/
endcal:
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ntam */

