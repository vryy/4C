/*======================================================================*/
/*!
\file
\brief Initialise certain variables in thermo-structure-interaction


<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 03/06
*/
#ifdef D_TSI


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
FILES allfiles;


/*----------------------------------------------------------------------*/
/*!
\brief general problem data

global variable general problem type
global variable GENPROB genprob is defined in global_control.c 

\auther bborn
\date 03/06
*/
extern GENPROB genprob;


/*----------------------------------------------------------------------*/
/*!
\brief vector of numfld FIELDs

defined in global_control.c

\author bborn
\date 03/06
*/
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
extern PAR par;


/*----------------------------------------------------------------------*/
/*!
\brief Alldyn dynamic control

pointer to allocate dynamic variables if needed
dedfined in global_control.c

\auther bborn
\date 03/06
*/
extern ALLDYNA *alldyn;

/*----------------------------------------------------------------------*/
/*!
\brief Load curve, load factor function
       number of load curves numcurve
       vector of structures of curves
       defined in input_curves.c
\author bborn
\date 03/06
*/
extern INT numcurve;
extern CURVE* curve;

/*======================================================================*/
/*!
\brief Check fields in thermal-structure-interaction

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
void tsi_init_chkfld(INT *numfld, /*!< total number of fields */
                     INT *numsf,  /*!< number (index) of structure field */
                     INT *numtf,  /*!< number (index) of thermal field */
                     FIELD **structfield,  /*!< structure field */
                     FIELD **thermfield)  /*!< thermal field */

{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_init_chkfld");
#endif

  /*--------------------------------------------------------------------*/
  /* check input  of fields */
  /* expected:    2 fields          0. structure
   *                                1. thermal
   */
  *numfld = genprob.numfld;  /* total number of fields; should be 2 */
  if (*numfld != 2)
  {
    dserror("Two fields are needed for TSI-problem!\n");
  }

  /* set structure field */
  *numsf = 0;
  if (&(field[*numsf]) != NULL)
  {
    *structfield = &(field[*numsf]);
  }
  else
  {
    dserror("field[%d] is not allocated! Structure field expected.", 
            *numsf);
  }

  /* check type of 0th field */
  if ((*structfield)->fieldtyp != structure)
  {
    dserror("FIELD 0 has to be structure\n");
  }

  *numtf = 1;
  if (&(field[*numtf]) != NULL)
  {
    *thermfield = &(field[*numtf]);
  }
  else
  {
    dserror("field[%d] is not allocated! Thermal field expected.", 
            *numtf);
  }

  /* checl type of 1st field */
  if ((*thermfield)->fieldtyp != thermal)
  {
    dserror("FIELD 1 has to be thermal\n");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of dyn_init_chkfld */



/*======================================================================*/
/*!
\brief Check if field consists only of 1 discretisation

\author bborn
\date 10/06
*/
void tsi_init_chkdis(FIELD *actfield,  /* current field */
                     INT *disnum)  /* discretisation index */
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_init_chkfld");
#endif

  if (actfield->ndis == 1)
  {
    *disnum = 0;
  }
  else
  {
    dserror("Only 1 discretisation currently allowed per field!\n");
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of tsi_init_chkdis */



/*======================================================================*/
/*!
\brief Set parameter sets of dynamic control
       Adjust parameters occuring of the 3 different control sets

\author bborn
\date 10/06
*/
void tsi_init_alldyn(INT numsf,  /* index of structure field */
                     INT numtf,  /* index of thermal field */
                     STRUCT_DYNAMIC **sdyn,  /* structure dynamic parameters */
                     THERM_DYNAMIC **tdyn,  /* thermal dynamic parameters */
                     TSI_DYNAMIC **tsidyn)  /* TSI dynamic parameters */
{
  INT itsidyn;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_init_alldyn");
#endif
  
  /* set pointers to dynamic controls */
  *sdyn = alldyn[numsf].sdyn;  /* structure */
  *tdyn = alldyn[numtf].tdyn;  /* thermal */

  /* We have 2 fields, the first (0) is the structure, the second (1) is
   * the thermal, ie genprob.numfld==2. Each field has a dynamic control 
   * unit stored as alldyn[0] and alldyn[1], respectively. The dynamic
   * control of the TSI coupling is kept at alldyn[2]. */
  itsidyn = genprob.numfld;
  *tsidyn = alldyn[itsidyn].tsidyn;  /* TSI */

  /* adjust/initialise dynamic controls */
  (*tsidyn)->time = 0.0;
  (*sdyn)->time = (*tsidyn)->time;
  (*tdyn)->time = (*tsidyn)->time;
  (*tsidyn)->step = 0;
  (*sdyn)->step = (*tsidyn)->step;
  (*tdyn)->step = (*tsidyn)->step;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of tsi_init_alldyn */



/*======================================================================*/
/*!
\brief Set indices of NODE sol arrays in terms of named variables

\author bborn
\date 10/06
*/
void tsi_init_nodsol(FIELD *structfield,
                     INT disnum_s,
                     FIELD *thermfield, 
                     INT disnum_t)
{
  ARRAY_POSITION *ipos_s;
  ARRAY_POSITION_SOL *isol_s;
  ARRAY_POSITION_SOLINC *isolinc_s;
  ARRAY_POSITION_SOLRES *isolres_s;
  ARRAY_POSITION *ipos_t;
  ARRAY_POSITION_SOL *isol_t;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_init_nodsol");
#endif

  /*====================================================================*/
  /* set names for NODE sol array indices of structure */
  ipos_s = &(structfield->dis[disnum_s].ipos);  /* pointer to position */
  /*--------------------------------------------------------------------*/
  /* sol array */
  ipos_s->num = 13;  /* dimension of sol array */
  isol_s = &(ipos_s->isol);  /* pointer to sol names */
  isol_s->disn = 0;  /* displacements d_{n+1}^k at t_{n+1} */
  isol_s->veln = 1;  /* velocities v_{n+1}^k at t_{n+1} */
  isol_s->accn = 2;  /* accelerations a_{n+1}^k at t_{n+1} */
  isol_s->disd = 3;  /* displacement d_n at t_n of presc. DOFs */
  isol_s->disdn = 4; /* displacements d_{n+1} of prescribed DOFs */
  isol_s->disdi = 5;  /* displ. increment d_{n+1} - d_n of presc. DOFs */
  isol_s->veldn = 6;  /* velocities v_{n+1} at t_{n+1} of prescribed DOFs */
  isol_s->accdn = 7;  /* accelerations a_{n+1} at t_{n+1} of presc. DOFs */
  isol_s->dis = 10;  /* displacements d_n at t_n */
  isol_s->vel = 11;  /* velocities v_n at t_n */
  isol_s->acc = 12;  /* accelerations a_n at t_n */
  /*--------------------------------------------------------------------*/
  /* sol_increment array */
  ipos_s->numincr = 3;
  isolinc_s = &(ipos_s->isolinc);
  isolinc_s->disinc = 0;  /* iterative incremental displacements */
  isolinc_s->fint = 1;  /* internal force f_{int;n} at t_n */
  isolinc_s->fintn = 2;  /* internal force at t_{n+1} */
  /*--------------------------------------------------------------------*/
  /* sol_residual array */
  ipos_s->numres = 1;
  isolres_s = &(ipos_s->isolres);
  isolres_s->disres = 0;

  /*====================================================================*/
  /* set named indices of NODE sol for thermal fields */
  ipos_t = &(thermfield->dis[disnum_t].ipos);
  /*--------------------------------------------------------------------*/
  /* sol array */
  ipos_t->num = 1;  /* dimension of sol array */
  isol_t = &(ipos_t->isol);  /* pointer to sol names */
  isol_t->tem = 0;  /* temperature at t_{n} */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of dyn_init_nodsol */



/*======================================================================*/
/*!
\brief Initialise load curves

\author bborn
\date 10/06
*/
void tsi_init_curve(TSI_DYNAMIC* tsidyn)
{
  INT actcurve;  /* curve index */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_init_curve");
#endif

  /*--------------------------------------------------------------------*/
  /* init all applied time curves -*/
  for (actcurve=0; actcurve<numcurve; actcurve++)
  {
    dyn_init_curve(actcurve, tsidyn->nstep, tsidyn->dt, tsidyn->maxtime);
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of tsi_init_curve() */


/*======================================================================*/
#endif  /* end of #ifdef D_TSI */
