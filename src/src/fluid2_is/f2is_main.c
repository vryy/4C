/*!----------------------------------------------------------------------
\file
\brief main routine inf-sup fluid2 element

<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET

#include "fluid2_is.h"
#include "../fluid2/fluid2_prototypes.h"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;
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

/*!---------------------------------------------------------------------
\brief main fluid2 control routine

<pre>                                                         genk 03/02
</pre>
\param  *actpart	 PARTITION    (i)
\param	*actintra	 INTRA        (i)
\param	*ele		 ELEMENT      (i)    actual element
\param	*estif_global	 ARRAY        (o)    element stiffness matrix
\param  *emass_global	 ARRAY        (o)    element mass matrix
\param  *eforce_global   ARRAY        (o)    element internal force vector
\param	*edforce_global  ARRAY        (o)    ele dirichl. force vector
\param	*action	         CALC_ACTION  (i)
\param	*hasdirich	 INT          (o)    flag
\param  *hasext          INT          (o)    flag
\param  *container       CONTAINER    (i)    container
\return void

------------------------------------------------------------------------*/
void fluid2_is(PARTITION   *actpart,
	       INTRA       *actintra,
	       ELEMENT     *ele,
	       ELEMENT     *eleke,
	       ARRAY       *estif_global,
	       ARRAY       *emass_global,
	       ARRAY       *eforce_global,
	       ARRAY       *edforce_global,
	       CALC_ACTION *action,
	       INT         *hasdirich,
	       INT         *hasext,
	       CONTAINER   *container
  )
{
/*----------------------------------------------------------------------*/
#ifdef D_FLUID2_IS
/*----------------------------------------------------------------------*/
  static INT              viscstr;
  static FLUID_DYNAMIC   *fdyn;
  ARRAY_POSITION* ipos;

#ifdef DEBUG
  dstrc_enter("fluid2_is");
#endif

  if (container!=NULL)
    ipos = &(field[genprob.numff].dis[container->disnum].ipos);
  else
    ipos = NULL;

/*------------------------------------------------- switch to do option */
  switch (*action)
  {
/*------------------------------------------------------ initialization */
  case calc_fluid_init:
    fdyn   = alldyn[genprob.numff].fdyn;
    viscstr= alldyn[genprob.numff].fdyn->viscstr;
/*------------------------------------------- init the element routines */
    f2_intg(0);
    f2is_calele(NULL,NULL,
		estif_global,emass_global,
		eforce_global,edforce_global,
		NULL,NULL,0,ipos,0,1);
    f2_iedg(NULL,ele,-1,1);
    break;

/*------------------------------------------- call the element routines */
  case calc_fluid:
/*---------------------------------------------------- multi-level FEM? */
    f2is_calele(ele,eleke,
		estif_global,emass_global,
		eforce_global,edforce_global,
		hasdirich,hasext,actintra->intra_rank,ipos,container->is_relax,0);
    break;

  default:
    dserror("action %d unknown",*action);
    break;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
#endif
  return;
}
#endif
