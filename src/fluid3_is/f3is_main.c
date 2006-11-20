/*!----------------------------------------------------------------------
\file
\brief main routine fluid3 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#include "fluid3_is.h"

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
\brief main fluid3 control routine

<pre>                                                         genk 05/02
</pre>
\param  *actpart	 PARTITION    (i)
\param	*actintra	 INTRA        (i)
\param	*ele		 ELEMENT      (i)    actual element
\param	*estif_global	 ARRAY        (o)    element stiffness matrix
\param  *emass_global	 ARRAY        (o)    element mass matrix
\param  *eforce_global   ARRAY        (o)    element force vecotr
\param	*edforce_global  ARRAY        (o)    ele dirichl. force vector
\param	*action	         CALC_ACTION  (i)
\param	*hasdirich	 INT          (o)    flag
\param  *hasext          INT          (o)    flag
\return void

------------------------------------------------------------------------*/
void fluid3_is(PARTITION   *actpart,
	      INTRA       *actintra,
	      ELEMENT     *ele,
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
#ifdef D_FLUID3_IS
  static INT              viscstr;
  static FLUID_DYNAMIC   *fdyn;
  ARRAY_POSITION* ipos;

  INT       i;        /* simply a counter */
  INT       ldflag;
  GVOL     *actgvol;
  GSURF    *actgsurf;
  DSURF    *actdsurf;

#ifdef DEBUG
  dstrc_enter("fluid3is");
#endif

  ipos = &(field[genprob.numff].dis[container->disnum].ipos);

/*------------------------------------------------- switch to do option */
  switch (*action)
  {
/*------------------------------------------------------ initialization */
  case calc_fluid_init:
/* ----------------------------------------- find number of fluid field */
    fdyn   = alldyn[genprob.numff].fdyn;
    viscstr= alldyn[genprob.numff].fdyn->viscstr;

/*------------------------------------------- init the element routines */
    f3_intg(0);
#if 0
    f3is_calele(NULL,estif_global,emass_global,
	      eforce_global,edforce_global,ipos,NULL,NULL,1);
#endif
    break;

/* call the element routines */
  case calc_fluid:
#if 0
      f3is_calele(ele,estif_global,emass_global,
		eforce_global,edforce_global,ipos,hasdirich,hasext,0);
#endif
    break;

  default:
    dserror("action unknown\n");
    break;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
#endif
  return;
} /* end of fluid3 */

