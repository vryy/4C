/*!----------------------------------------------------------------------
\file
\brief main routine fluid3_pro element

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID3_PRO
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3_prototypes.h"
#include "fluid3pro_prototypes.h"
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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             genk 04/02    |
 *----------------------------------------------------------------------*/

/*!---------------------------------------------------------------------
\brief main fluid3_pro control routine

                                                        basol 10/02
</pre>
\param	*action	         CALC_ACTION    (i)
\param  *ele	         ELEMENT	(i)   actual element
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *lmass_global    ARRAY	        (o)   lumped mass matrix
\param  *gradopr_global  ARRAY	        (o)   gradient operator
\param  *etforce_global  ARRAY	        (o)   element time force
\param  *etforce1_global ARRAY	        (o)   element time force part1
\param  *etforce2_global ARRAY	        (o)   element time force part2
\param  *etforce3_global ARRAY	        (o)   element time force part3
\param  *eiforce_global  ARRAY	        (o)   ele iteration force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param	 double	         dt             (i)   incremental time step

\return void

------------------------------------------------------------------------*/
void fluid3_pro(     PARTITION     *actpart,
                     INTRA         *actintra,
		     ELEMENT       *ele,
		     ARRAY         *estif_global,
		     ARRAY         *emass_global,
		     ARRAY         *lmass_global,
		     ARRAY         *gradopr_global,
		     ARRAY         *eforce_global,
		     ARRAY         *edforce_global,
		     ARRAY         *gforce_global,
		     CALC_ACTION   *action,
                     INT           *hasdirich,
                     INT           *hasext
  )
{
/*----------------------------------------------------------------------*/
#ifdef D_FLUID3_PRO
/*----------------------------------------------------------------------*/
  ARRAY_POSITION* ipos;

#ifdef DEBUG
  dstrc_enter("fluid3_pro");
#endif

  ipos = &(field[genprob.numff].dis[0].ipos);

/*------------------------------------------------- switch to do option */
  switch (*action)
  {
/*------------------------------------------------------ initialisation */
  case calc_fluid_init:
    /*---------------------------------------- init the element routines */
    f3_intg(0);
    /*----------- f3_pro element is initialized, some arrays are defined */
    f3pro_calinit(estif_global,emass_global,
		  lmass_global,gradopr_global,eforce_global,
		  edforce_global,gforce_global,ipos);
    break;
/*------------------------------------------- call the element routines */
  case calc_fluid:
    f3pro_calele(ele,estif_global,emass_global,
                 lmass_global,gradopr_global,eforce_global,
                 edforce_global,gforce_global,ipos,hasdirich,hasext);
    break;
/*----------------------------------------------------------------------*/

/* integrate errors for beltrami and kim-moin */
  case calc_fluid_error:
    dserror("Error integration not yet implemented for fluid3pro!!");
    break;

  case calc_fluid_stress:
    /* some time later. maybe. */
    break;

  default:
    dserror("action unknown");
    break;
  } /* end swtich (*action) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
#endif
/*----------------------------------------------------------------------*/
  return;
} /* end of fluid3_pro */
/*! @} (documentation module close)*/
