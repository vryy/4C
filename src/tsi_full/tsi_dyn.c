/*----------------------------------------------------------------------*/
/*!
\file
\brief Main chooser for different thermal structural interactive analyses

According to the problem type and its KIND variable, different thermal
structure interaction modes are selected. 

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"
#include "../io/io.h"
#include "../solver/solver.h"
#include "tsi_prototypes.h"


/*----------------------------------------------------------------------*/
/*!
\brief File pointers

This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
struct _FILES allfiles;


/*----------------------------------------------------------------------*/
/*!
\brief general problem data

global variable general problem type
global variable GENPROB genprob is defined in global_control.c 

\auther bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
extern struct _GENPROB genprob;


/*----------------------------------------------------------------------*/
/*!
\brief vector of numfld FIELDs

defined in global_control.c

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
extern struct _FIELD *field;


/*----------------------------------------------------------------------*/
/*!
\brief Ranks and communicators of parallelism

This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h

Unsure if we need this here, however, in the long term, it will become
necessary.

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
extern struct _PAR par;


/*----------------------------------------------------------------------*/
/*!
\brief Alldyn dynamic control

pointer to allocate dynamic variables if needed
dedfined in global_control.c

\auther bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
extern ALLDYNA *alldyn;





/*----------------------------------------------------------------------*/
/*!
\brief Thermal structure interaction analyses chooser

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
void tsi_dyn()
{

  /* fields */
  static INT numfld;  /* total number of fields */
  static INT numsf;  /* number (index) of structure field */
  static INT numtf;  /* number (index) of thermal field */
  static FIELD *structfield;  /* pointer to structure field */
  static FIELD *thermfield;  /* pointer to thermal field */
  
  /* discretisations */
  INT s_disnum_calc = 0;
  INT s_disnum_io = 0;
  INT t_disnum_calc = 0;
  INT t_disnum_io = 0;

  /* dynamic control */
  static STRUCT_DYNAMIC *sdyn;
  static THERM_DYNAMIC *tdyn;
  static TSI_DYNAMIC *tsidyn;
  INT itsidyn;

  /* files */
  FILE *out = allfiles.out_out;  /* the OUT file */

#ifdef DEBUG
  dstrc_enter("dyn_tsi");
#endif

  /*--------------------------------------------------------------------*/
  /* check input  of fields */
  /* expected:    2 fields          0. structure
   *                                1. thermal
   */
  numfld = genprob.numfld;  /* total number of fiels; should be 2 */
  dsassert(numfld==2, "Two fields needed for TSI-problem!\n");
  
  /* set structure field */
  numsf = 0;
  if (&(field[numsf]) != NULL)
  {
    structfield = &(field[numsf]);
  }
  else
  {
    dserror("field[%d] is not allocated! Structure field expected.", 
            &numsf);
  }
  dsassert(structfield->fieldtyp==structure, 
           "FIELD 0 has to be structure\n");

  /* set thermal field */
  numtf = 1;
  if (&(field[numsf]) != NULL)
  {
    thermfield = &(field[numtf]);
  }
  else
  {
    dserror("field[%d] is not allocated! Thermal field expected.", 
            &numtf);
  }
  dsassert(thermfield->fieldtyp==thermal, 
           "FIELD 1 has to be thermal\n");


  /*--------------------------------------------------------------------*/
  /* set pointers to dynamic controls */
  sdyn= alldyn[numsf].sdyn;  /* structure */
  tdyn= alldyn[numtf].tdyn;  /* thermal */
  /* We have 2 fields, the first (0) is the structure, the second (1) is
   * the thermal, ie genprob.numfld==2. Each field has a dynamic control 
   * unit stored as alldyn[0] and alldyn[1], respectively. The dynamic
   * control of the TSI coupling is kept at alldyn[2]. */
  itsidyn = genprob.numfld;
  tsidyn = alldyn[itsidyn].tsidyn;  /* TSI */

  /* adjust/initialise dynamic controls */
  tsidyn->time = 0.0;
  sdyn->time = tsidyn->time;
  tdyn->time = tsidyn->time;
  tsidyn->step = 0;
  sdyn->step = tsidyn->step;
  tdyn->step = tsidyn->step;
  


  /*--------------------------------------------------------------------*/
  /* create coupling of structural and thermal field */
  tsi_initcoupling(structfield, s_disnum_calc,
                   thermfield, t_disnum_calc);

  /*--------------------------------------------------------------------*/
  /* select different analyses */
  switch (tsidyn->kind)
  {
    /* */
    case tsi_full:
      /* tsi_dyn_full(); */
      break;
    /* */
    case tsi_therm_stat_struct_dyn: 
      tsi_stat_therm();
      tsi_dyn_struct();
/*      dyn_nln_structural(); */
      break;
    /* */
    case tsi_therm_pred_struct_dyn:
      break;
    /* */
    default:
      dserror("Unknown KIND of TSI");
      break;
  }  /* end of switch (tsidyn->kind) */


  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of dyn_tsi */
