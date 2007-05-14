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
#ifdef BINIO
#include "../io/io.h"
#endif
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
FILES allfiles;


/*----------------------------------------------------------------------*/
/*!
\brief general problem data

global variable general problem type
global variable GENPROB genprob is defined in global_control.c 

\auther bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
extern GENPROB genprob;


/*----------------------------------------------------------------------*/
/*!
\brief vector of numfld FIELDs

defined in global_control.c

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
extern FIELD *field;


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
extern PAR par;


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
  INT numfld;  /* total number of fields */
  INT numsf;  /* number (index) of structure field */
  INT numtf;  /* number (index) of thermal field */
  FIELD *structfield;  /* pointer to structure field */
  FIELD *thermfield;  /* pointer to thermal field */
  
  /* discretisations */
  INT disnum_s;
  INT disnum_t;

  /* dynamic control */
  STRUCT_DYNAMIC *sdyn;
  THERM_DYNAMIC *tdyn;
  TSI_DYNAMIC *tsidyn;

  /* files */
  /*FILE *out = allfiles.out_out;*/  /* the OUT file */ 

#ifdef DEBUG
  dstrc_enter("tsi_dyn");
#endif

  /*--------------------------------------------------------------------*/
  /* check fields */
  tsi_init_chkfld(&numfld, &numsf, &numtf, &structfield, &thermfield);

  /*--------------------------------------------------------------------*/
  /* check discretisation */
  tsi_init_chkdis(structfield, &disnum_s);
  tsi_init_chkdis(thermfield, &disnum_t);

  /*--------------------------------------------------------------------*/
  /* associate dynamic control parameters */
  tsi_init_alldyn(numsf, numtf, &sdyn, &tdyn, &tsidyn);
  
  /*--------------------------------------------------------------------*/
  /* create coupling of structural and thermal field */
  tsi_coupling(structfield, disnum_s,
               thermfield, disnum_t);

  /*--------------------------------------------------------------------*/
  /* set named indices of NODE sol arrays */
  tsi_init_nodsol(structfield, disnum_s,
                  thermfield, disnum_t);

  /*--------------------------------------------------------------------*/
  /* initialise load curves */
  tsi_init_curve(tsidyn);

  /*--------------------------------------------------------------------*/
  /* select different analyses */
  switch (tsidyn->kind)
  {
    /* prescribed TSI */
    case tsi_therm_pred_struct_dyn:
      dserror("Sorry, predefined temperature field is not implemented!");
      break;
    /* semi TSI : static thermal and dynamic structure field */
    case tsi_therm_stat_struct_genalp: 
      tsi_th_stat(disnum_s, disnum_t);
      tsi_st_genalp(disnum_s, disnum_t);
      break;
    /* semi TSI : static thermal and dynamic structure field */
    case tsi_therm_stat_struct_cendif:
      tsi_th_stat(disnum_s, disnum_t);
      tsi_st_cendif(disnum_s, disnum_t);
      break;
    /* semi TSI : static thermal and dynamic structure field */
    case tsi_therm_stat_struct_fehlbg:
      tsi_th_stat(disnum_s, disnum_t);
      tsi_st_fehlbg(disnum_s, disnum_t);
      break;
    /* semi TSI : stationary prescribed thermal and dynamic structure field */
    case tsi_therm_presc_struct_genalp:
      tsi_th_presc_st_genalp(disnum_s, disnum_t);
      break;
    /* semi-coupled TSI : instationary thermal and dynamic structural field */
    case tsi_therm_ost_struct_genalp:
      tsi_th_ost_st_genalp(disnum_s, disnum_t);
      break;
    /* full explicit TIS with Fehlberg */
    case tsi_full_fehlbg:
      dserror("Sorry, full thermo-structure interaction is not implemented!\n");
      break;
    /* default */
    case tsi_none:
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
