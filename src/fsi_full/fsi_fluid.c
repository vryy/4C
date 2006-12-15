/*!
\file
\brief fsi fluid.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

\author u.kue
\date 12/06

*/

#ifdef D_FSI

#include "../headers/standardtypes.h"
#include "fsi.h"
#include "fsi_prototypes.h"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void fsi_fluid_setup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  )
{
  FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */
  
#ifdef DEBUG
  dstrc_enter("fsi_fluid_setup");
#endif

  fdyn          = alldyn[genprob.numff].fdyn;
  
  switch (genprob.probtyp)
  {
  case prb_fluid:
  case prb_fsi:
    fsi_fluid_imp_setup(work,actfield,disnum_calc,disnum_io);
    break;
#ifdef D_FLUID_PM
  case prb_pfsi:
    switch (fdyn->dyntyp)
    {
    case dyntyp_pm_cont:
      fsi_fluid_pm_cont_setup(work,actfield,disnum_calc,disnum_io);
      break;
    case dyntyp_pm_discont:
      fsi_fluid_pm_discont_setup(work,actfield,disnum_calc,disnum_io);
      break;
    case dyntyp_pm_cont_laplace:
      fsi_fluid_pm_laplace_setup(work,actfield,disnum_calc,disnum_io);
      break;
    default:
      dserror("unsupported dynamic type %d for fluid projection", fdyn->dyntyp);
    }
    break;
#endif
  default:
    dserror("unsupported problem type %d for fluid in fsi", genprob.probtyp);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void fsi_fluid_calc(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io,
  FIELD          *alefield,
  INT             adisnum_calc
  )
{
  FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */
  
#ifdef DEBUG
  dstrc_enter("fsi_fluid_calc");
#endif

  fdyn          = alldyn[genprob.numff].fdyn;
  
  switch (genprob.probtyp)
  {
  case prb_fluid:
  case prb_fsi:
    fsi_fluid_imp_calc(work,actfield,disnum_calc,disnum_io,alefield,adisnum_calc);
    break;
#ifdef D_FLUID_PM
  case prb_pfsi:
    switch (fdyn->dyntyp)
    {
    case dyntyp_pm_cont:
      fsi_fluid_pm_cont_calc(work,actfield,disnum_calc,disnum_io,alefield,adisnum_calc);
      break;
    case dyntyp_pm_discont:
      fsi_fluid_pm_discont_calc(work,actfield,disnum_calc,disnum_io,alefield,adisnum_calc);
      break;
    case dyntyp_pm_cont_laplace:
      fsi_fluid_pm_laplace_calc(work,actfield,disnum_calc,disnum_io,alefield,adisnum_calc);
      break;
    default:
      dserror("unsupported dynamic type %d for fluid projection", fdyn->dyntyp);
    }
    break;
#endif
  default:
    dserror("unsupported problem type %d for fluid in fsi", genprob.probtyp);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void fsi_fluid_final(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  )
{
  FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */
  
#ifdef DEBUG
  dstrc_enter("fsi_fluid_final");
#endif

  fdyn          = alldyn[genprob.numff].fdyn;
  
  switch (genprob.probtyp)
  {
  case prb_fluid:
  case prb_fsi:
    fsi_fluid_imp_final(work,actfield,disnum_calc,disnum_io);
    break;
#ifdef D_FLUID_PM
  case prb_pfsi:
    switch (fdyn->dyntyp)
    {
    case dyntyp_pm_cont:
      fsi_fluid_pm_cont_final(work,actfield,disnum_calc,disnum_io);
      break;
    case dyntyp_pm_discont:
      fsi_fluid_pm_discont_final(work,actfield,disnum_calc,disnum_io);
      break;
    case dyntyp_pm_cont_laplace:
      fsi_fluid_pm_laplace_final(work,actfield,disnum_calc,disnum_io);
      break;
    default:
      dserror("unsupported dynamic type %d for fluid projection", fdyn->dyntyp);
    }
    break;
#endif
  default:
    dserror("unsupported problem type %d for fluid in fsi", genprob.probtyp);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void fsi_fluid_sd(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  )
{
  FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */
  
#ifdef DEBUG
  dstrc_enter("fsi_fluid_sd");
#endif

  fdyn          = alldyn[genprob.numff].fdyn;
  
  switch (genprob.probtyp)
  {
  case prb_fluid:
  case prb_fsi:
    fsi_fluid_imp_sd(work,actfield,disnum_calc,disnum_io);
    break;
#ifdef D_FLUID_PM
  case prb_pfsi:
    switch (fdyn->dyntyp)
    {
    case dyntyp_pm_cont:
      fsi_fluid_pm_cont_sd(work,actfield,disnum_calc,disnum_io);
      break;
    case dyntyp_pm_discont:
      fsi_fluid_pm_discont_sd(work,actfield,disnum_calc,disnum_io);
      break;
    case dyntyp_pm_cont_laplace:
      fsi_fluid_pm_laplace_sd(work,actfield,disnum_calc,disnum_io);
      break;
    default:
      dserror("unsupported dynamic type %d for fluid projection", fdyn->dyntyp);
    }
    break;
#endif
  default:
    dserror("unsupported problem type %d for fluid in fsi", genprob.probtyp);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void fsi_fluid_output(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  )
{
  FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */
  
#ifdef DEBUG
  dstrc_enter("fsi_fluid_output");
#endif

  fdyn          = alldyn[genprob.numff].fdyn;
  
  switch (genprob.probtyp)
  {
  case prb_fluid:
  case prb_fsi:
    fsi_fluid_imp_output(work,actfield,disnum_calc,disnum_io);
    break;
#ifdef D_FLUID_PM
  case prb_pfsi:
    switch (fdyn->dyntyp)
    {
    case dyntyp_pm_cont:
      fsi_fluid_pm_cont_output(work,actfield,disnum_calc,disnum_io);
      break;
    case dyntyp_pm_discont:
      fsi_fluid_pm_discont_output(work,actfield,disnum_calc,disnum_io);
      break;
    case dyntyp_pm_cont_laplace:
      fsi_fluid_pm_laplace_output(work,actfield,disnum_calc,disnum_io);
      break;
    default:
      dserror("unsupported dynamic type %d for fluid projection", fdyn->dyntyp);
    }
    break;
#endif
  default:
    dserror("unsupported problem type %d for fluid in fsi", genprob.probtyp);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void fsi_fluid_cleanup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  )
{
  FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */
  
#ifdef DEBUG
  dstrc_enter("fsi_fluid_cleanup");
#endif

  fdyn          = alldyn[genprob.numff].fdyn;
  
  switch (genprob.probtyp)
  {
  case prb_fluid:
  case prb_fsi:
    fsi_fluid_imp_cleanup(work,actfield,disnum_calc,disnum_io);
    break;
#ifdef D_FLUID_PM
  case prb_pfsi:
    switch (fdyn->dyntyp)
    {
    case dyntyp_pm_cont:
      fsi_fluid_pm_cont_cleanup(work,actfield,disnum_calc,disnum_io);
      break;
    case dyntyp_pm_discont:
      fsi_fluid_pm_discont_cleanup(work,actfield,disnum_calc,disnum_io);
      break;
    case dyntyp_pm_cont_laplace:
      fsi_fluid_pm_laplace_cleanup(work,actfield,disnum_calc,disnum_io);
      break;
    default:
      dserror("unsupported dynamic type %d for fluid projection", fdyn->dyntyp);
    }
    break;
#endif
  default:
    dserror("unsupported problem type %d for fluid in fsi", genprob.probtyp);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


#endif
