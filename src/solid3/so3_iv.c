/*======================================================================*/
/*!
\file
\brief Iterative update of internal variables

<pre>
Maintainer: Burkhard Bornemann 
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/


/*----------------------------------------------------------------------*/
#ifdef D_SOLID3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

\author mf
\date 10/06
*/
extern GENPROB genprob;



/*======================================================================*/
/*!
\brief Iterative update of element variables

The element stiffness matrix, ie the tangent operator, is determined
for the linear, 3dim elasticity

\param   *container     CONTAINER        (i)  see container.h
\param   *ele           ELEMENT          (i)  pointer to current element
\param   *gpshade       SO3_GPSHAPEDERIV (i)  Gauss point coords 
\param   *mat           MATERIAL         (i)
\param   *eforc_global  ARRAY            (o)  global vector for internal 
                                              forces (initialized!)
\param   *estif_global  ARRAY            (o)  element stiffness matrix
\param   *emass_global  ARRAY            (o)  element mass matrix
\return void

\author mf
\date 10/06
*/
void so3_iv_upd(const CONTAINER* container,
                ELEMENT* ele,
                const MATERIAL* mat)
{

  /*--------------------------------------------------------------------*/
  /* start */
#ifdef DEBUG
  dstrc_enter("so3_iv_upd");
#endif

  /*--------------------------------------------------------------------*/
  /* material internal variables */
  so3_mat_mivupd(container, ele, mat);

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_iv_upd() */


/*======================================================================*/
#endif /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
