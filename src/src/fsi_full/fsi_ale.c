/*!----------------------------------------------------------------------
\file
\brief ale control part of fsi-problems

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
/*!
\addtogroup FSI
*//*! @{ (documentation module open)*/
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../ale3/ale3.h"
#include "fsi_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


void fsi_ale_setup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
  ALE_DYNAMIC *adyn;
  adyn = alldyn[genprob.numaf].adyn;

#ifdef DEBUG
  dstrc_enter("fsi_ale_setup");
#endif

  work->ale_field = actfield;

  switch (adyn->typ)
  {
/*---------------------------------------- purely linear calculation ---*/
  case classic_lin:
    fsi_ale_lin_setup(work,actfield,disnum_calc,disnum_io);
    break;
/*------------------- incremental calculation stiffened with min J_e ---*/
  case min_Je_stiff:
    fsi_ale_nln_setup(work,actfield,disnum_calc,disnum_io);
    break;
/*--------------------------------------------- two step calculation ---*/
/*  calculation in two steps per timestep following Chiandussi et al. in
    'A simple method for automatic update of finite element meshes'
    Commun. Numer. Meth. Engng. 2000; 16: 1-19                          */
  case two_step:
    fsi_ale_2step_setup(work,actfield,disnum_calc,disnum_io);
    break;
/*--------------------------------------------------- spring analogy ---*/
/*  calculation following Farhat et al. in 'Torsional springs for
    two-dimensional dynamic unstructured fluid meshes' Comput. Methods
    Appl. Mech. Engrg. 163 (1998) 231-245 */
  case springs:
    fsi_ale_spring_setup(work,actfield,disnum_calc,disnum_io);
    break;
/*------------------------------------------------ Laplace smoothing ---*/
/*  calculation following Loehner et al. in 'Improved ALE mesh velocities
    for moving bodies' Commun. num. methd. engng. 12 (1996) 599-608 */
  case laplace:
    fsi_ale_laplace_setup(work,actfield,disnum_calc,disnum_io);
    break;
/*-------------------------------------------- Large Amplitude Sloshing */
/*  just a interpolation for large amplitude sloshing                    */
  case LAS:
    dserror("Large Amplitude Sloshing broken");
    fsi_ale_LAS_setup(work,actfield,disnum_calc,disnum_io);
    break;
/*---------------------------------------------------------- default ---*/
  default:
    dserror("unknown ale type %d for fsi", adyn->typ);
    break;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


void fsi_ale_calc(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  )
{
  ALE_DYNAMIC *adyn;
  adyn = alldyn[genprob.numaf].adyn;

#ifdef DEBUG
  dstrc_enter("fsi_ale_calc");
#endif

  switch (adyn->typ)
  {
/*---------------------------------------- purely linear calculation ---*/
  case classic_lin:
    fsi_ale_lin_calc(work,actfield,disnum_calc,disnum_io,structfield,sdisnum);
    break;
/*------------------- incremental calculation stiffened with min J_e ---*/
  case min_Je_stiff:
    fsi_ale_nln_calc(work,actfield,disnum_calc,disnum_io);
    break;
/*--------------------------------------------- two step calculation ---*/
/*  calculation in two steps per timestep following Chiandussi et al. in
    'A simple method for automatic update of finite element meshes'
    Commun. Numer. Meth. Engng. 2000; 16: 1-19                          */
  case two_step:
    fsi_ale_2step_calc(work,actfield,disnum_calc,disnum_io);
    break;
/*--------------------------------------------------- spring analogy ---*/
/*  calculation following Farhat et al. in 'Torsional springs for
    two-dimensional dynamic unstructured fluid meshes' Comput. Methods
    Appl. Mech. Engrg. 163 (1998) 231-245 */
  case springs:
    fsi_ale_spring_calc(work,actfield,disnum_calc,disnum_io);
    break;
/*------------------------------------------------ Laplace smoothing ---*/
/*  calculation following Loehner et al. in 'Improved ALE mesh velocities
    for moving bodies' Commun. num. methd. engng. 12 (1996) 599-608 */
  case laplace:
    fsi_ale_laplace_calc(work,actfield,disnum_calc,disnum_io);
    break;
/*-------------------------------------------- Large Amplitude Sloshing */
/*  just a interpolation for large amplitude sloshing                    */
  case LAS:
    dserror("Large Amplitude Sloshing broken");
    fsi_ale_LAS_calc(work,actfield,disnum_calc,disnum_io);
    break;
/*---------------------------------------------------------- default ---*/
  default:
    dserror("unknown ale type %d for fsi", adyn->typ);
    break;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


void fsi_ale_final(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
  ALE_DYNAMIC *adyn;
  adyn = alldyn[genprob.numaf].adyn;

#ifdef DEBUG
  dstrc_enter("fsi_ale_final");
#endif

  switch (adyn->typ)
  {
/*---------------------------------------- purely linear calculation ---*/
  case classic_lin:
    fsi_ale_lin_final(work,actfield,disnum_calc,disnum_io);
    break;
/*------------------- incremental calculation stiffened with min J_e ---*/
  case min_Je_stiff:
    fsi_ale_nln_final(work,actfield,disnum_calc,disnum_io);
    break;
/*--------------------------------------------- two step calculation ---*/
/*  calculation in two steps per timestep following Chiandussi et al. in
    'A simple method for automatic update of finite element meshes'
    Commun. Numer. Meth. Engng. 2000; 16: 1-19                          */
  case two_step:
    fsi_ale_2step_final(work,actfield,disnum_calc,disnum_io);
    break;
/*--------------------------------------------------- spring analogy ---*/
/*  calculation following Farhat et al. in 'Torsional springs for
    two-dimensional dynamic unstructured fluid meshes' Comput. Methods
    Appl. Mech. Engrg. 163 (1998) 231-245 */
  case springs:
    fsi_ale_spring_final(work,actfield,disnum_calc,disnum_io);
    break;
/*------------------------------------------------ Laplace smoothing ---*/
/*  calculation following Loehner et al. in 'Improved ALE mesh velocities
    for moving bodies' Commun. num. methd. engng. 12 (1996) 599-608 */
  case laplace:
    fsi_ale_laplace_final(work,actfield,disnum_calc,disnum_io);
    break;
/*-------------------------------------------- Large Amplitude Sloshing */
/*  just a interpolation for large amplitude sloshing                    */
  case LAS:
    dserror("Large Amplitude Sloshing broken");
    fsi_ale_LAS_final(work,actfield,disnum_calc,disnum_io);
    break;
/*---------------------------------------------------------- default ---*/
  default:
    dserror("unknown ale type %d for fsi", adyn->typ);
    break;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


void fsi_ale_sd(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  )
{
  ALE_DYNAMIC *adyn;
  adyn = alldyn[genprob.numaf].adyn;

#ifdef DEBUG
  dstrc_enter("fsi_ale_sd");
#endif

  switch (adyn->typ)
  {
/*---------------------------------------- purely linear calculation ---*/
  case classic_lin:
    fsi_ale_lin_sd(work,actfield,disnum_calc,disnum_io,structfield,sdisnum);
    break;
/*------------------- incremental calculation stiffened with min J_e ---*/
  case min_Je_stiff:
    fsi_ale_nln_sd(work,actfield,disnum_calc,disnum_io,structfield,sdisnum);
    break;
/*--------------------------------------------- two step calculation ---*/
/*  calculation in two steps per timestep following Chiandussi et al. in
    'A simple method for automatic update of finite element meshes'
    Commun. Numer. Meth. Engng. 2000; 16: 1-19                          */
  case two_step:
    fsi_ale_2step_sd(work,actfield,disnum_calc,disnum_io,structfield,sdisnum);
    break;
/*--------------------------------------------------- spring analogy ---*/
/*  calculation following Farhat et al. in 'Torsional springs for
    two-dimensional dynamic unstructured fluid meshes' Comput. Methods
    Appl. Mech. Engrg. 163 (1998) 231-245 */
  case springs:
    fsi_ale_spring_sd(work,actfield,disnum_calc,disnum_io,structfield,sdisnum);
    break;
/*------------------------------------------------ Laplace smoothing ---*/
/*  calculation following Loehner et al. in 'Improved ALE mesh velocities
    for moving bodies' Commun. num. methd. engng. 12 (1996) 599-608 */
  case laplace:
    fsi_ale_laplace_sd(work,actfield,disnum_calc,disnum_io,structfield,sdisnum);
    break;
/*-------------------------------------------- Large Amplitude Sloshing */
/*  just a interpolation for large amplitude sloshing                    */
  case LAS:
    dserror("Large Amplitude Sloshing broken");
    fsi_ale_LAS_sd(work,actfield,disnum_calc,disnum_io,structfield,sdisnum);
    break;
/*---------------------------------------------------------- default ---*/
  default:
    dserror("unknown ale type %d for fsi", adyn->typ);
    break;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


void fsi_ale_output(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
  ALE_DYNAMIC *adyn;
  adyn = alldyn[genprob.numaf].adyn;

#ifdef DEBUG
  dstrc_enter("fsi_ale_output");
#endif

  switch (adyn->typ)
  {
/*---------------------------------------- purely linear calculation ---*/
  case classic_lin:
    fsi_ale_lin_output(work,actfield,disnum_calc,disnum_io);
    break;
/*------------------- incremental calculation stiffened with min J_e ---*/
  case min_Je_stiff:
    fsi_ale_nln_output(work,actfield,disnum_calc,disnum_io);
    break;
/*--------------------------------------------- two step calculation ---*/
/*  calculation in two steps per timestep following Chiandussi et al. in
    'A simple method for automatic update of finite element meshes'
    Commun. Numer. Meth. Engng. 2000; 16: 1-19                          */
  case two_step:
    fsi_ale_2step_output(work,actfield,disnum_calc,disnum_io);
    break;
/*--------------------------------------------------- spring analogy ---*/
/*  calculation following Farhat et al. in 'Torsional springs for
    two-dimensional dynamic unstructured fluid meshes' Comput. Methods
    Appl. Mech. Engrg. 163 (1998) 231-245 */
  case springs:
    fsi_ale_spring_output(work,actfield,disnum_calc,disnum_io);
    break;
/*------------------------------------------------ Laplace smoothing ---*/
/*  calculation following Loehner et al. in 'Improved ALE mesh velocities
    for moving bodies' Commun. num. methd. engng. 12 (1996) 599-608 */
  case laplace:
    fsi_ale_laplace_output(work,actfield,disnum_calc,disnum_io);
    break;
/*-------------------------------------------- Large Amplitude Sloshing */
/*  just a interpolation for large amplitude sloshing                    */
  case LAS:
    dserror("Large Amplitude Sloshing broken");
    fsi_ale_LAS_output(work,actfield,disnum_calc,disnum_io);
    break;
/*---------------------------------------------------------- default ---*/
  default:
    dserror("unknown ale type %d for fsi", adyn->typ);
    break;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


void fsi_ale_cleanup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
  ALE_DYNAMIC *adyn;
  adyn = alldyn[genprob.numaf].adyn;

#ifdef DEBUG
  dstrc_enter("fsi_ale_cleanup");
#endif

  switch (adyn->typ)
  {
/*---------------------------------------- purely linear calculation ---*/
  case classic_lin:
    fsi_ale_lin_cleanup(work,actfield,disnum_calc,disnum_io);
    break;
/*------------------- incremental calculation stiffened with min J_e ---*/
  case min_Je_stiff:
    fsi_ale_nln_cleanup(work,actfield,disnum_calc,disnum_io);
    break;
/*--------------------------------------------- two step calculation ---*/
/*  calculation in two steps per timestep following Chiandussi et al. in
    'A simple method for automatic update of finite element meshes'
    Commun. Numer. Meth. Engng. 2000; 16: 1-19                          */
  case two_step:
    fsi_ale_2step_cleanup(work,actfield,disnum_calc,disnum_io);
    break;
/*--------------------------------------------------- spring analogy ---*/
/*  calculation following Farhat et al. in 'Torsional springs for
    two-dimensional dynamic unstructured fluid meshes' Comput. Methods
    Appl. Mech. Engrg. 163 (1998) 231-245 */
  case springs:
    fsi_ale_spring_cleanup(work,actfield,disnum_calc,disnum_io);
    break;
/*------------------------------------------------ Laplace smoothing ---*/
/*  calculation following Loehner et al. in 'Improved ALE mesh velocities
    for moving bodies' Commun. num. methd. engng. 12 (1996) 599-608 */
  case laplace:
    fsi_ale_laplace_cleanup(work,actfield,disnum_calc,disnum_io);
    break;
/*-------------------------------------------- Large Amplitude Sloshing */
/*  just a interpolation for large amplitude sloshing                    */
  case LAS:
    dserror("Large Amplitude Sloshing broken");
    fsi_ale_LAS_cleanup(work,actfield,disnum_calc,disnum_io);
    break;
/*---------------------------------------------------------- default ---*/
  default:
    dserror("unknown ale type %d for fsi", adyn->typ);
    break;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


#endif
/*! @} (documentation module close)*/
#endif
