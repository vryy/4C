/*!----------------------------------------------------------------------
\file
\brief contains the routine 'wge_call_mat' which calls the material law
       of gradient enhanced wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"

/*!
\addtogroup WALLGE
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief contains the routine 'wge_call_mat' which calls the material law
       of gradient enhanced wall element

\param  *ele             ELEMENT        I: Actual element
\param  *mat             MATERIAL       I: Material of actual element
\param  *bopd            DOUBLE         I: B-operator for displacements
\param  *functe          DOUBLE         I: Ansatz-funct.for equiv.strain
\param **bope            DOUBLE         I: B-operator for equiv.strain
\param   ip              DOUBLE         I: ID of actual Gauss point
\param  *stress          DOUBLE       I/O: stresses at GP(vector)
\param  *eps_vl          DOUBLE       I/O: local equiv. strain at GP(scalar)
\param  *eps_vnl         DOUBLE       I/O: nonlocal equiv. strain at GP(scalar)
\param  *grad_eps_vnl    DOUBLE       I/O: grad of nonl.equi.strain(vector)
\param **D               DOUBLE       I/O: tangent d sig/d eps(matrix)
\param  *E               DOUBLE       I/O: tangent d sig/d eps_v(vector)
\param  *F               DOUBLE       I/O: tangent d eps_v/d eps(vector)
\param   istore          INT            I: flag: istore=1 -> update
\param   newval          INT            I: flag: newval=1 -> just stress

\return void
\sa calling:   wge_mat_damage();
    called by: wgestatic_ke();
*----------------------------------------------------------------------*/
void wge_call_mat(ELEMENT   *ele,     /* actual element                */
                  MATERIAL  *mat,     /* material                      */
                  DOUBLE   **bopd,    /* B-operator for displacements  */
                  DOUBLE    *functe,  /* Ansatz-funct.for equiv.strain */
                  DOUBLE   **bope,    /* B-operator for equiv.strain   */
                  INT        ip,      /* ID of actual Gauss point      */
                  DOUBLE    *stress,  /* stresses at GP                */
                  DOUBLE    *eps_vl,  /* local equiv. strain at GP     */
                  DOUBLE    *eps_vnl, /* nonlocal equiv. strain at GP  */
                  DOUBLE    *grad_eps_vnl, /* grad of nonl.equi.strain */
                  DOUBLE   **D,       /* tangent d sig/d eps           */
                  DOUBLE    *E,       /* tangent d sig/d eps_v         */
                  DOUBLE    *F,       /* tangent d eps_v/d eps         */
                  INT        istore,  /* flag: istore=1 -> update      */
                  INT        newval)  /* flag: newval=1 -> just stress */
{
#ifdef D_WALLGE


#ifdef DEBUG
dstrc_enter("wge_call_mat");
#endif
/*----------------------------------------------------------------------*/
switch(mat->mattyp)
{
case m_damage_ge:
   wge_mat_damage(mat->m.damage_ge->equival,
                  mat->m.damage_ge->damtyp,
                  mat->m.damage_ge->youngs,
                  mat->m.damage_ge->nue,
                  mat->m.damage_ge->kappa_0,
                  mat->m.damage_ge->kappa_m,
                  mat->m.damage_ge->alpha,
                  mat->m.damage_ge->beta,
                  mat->m.damage_ge->k_fac,
                  ele,
                  bopd,
                  functe,
                  bope,
                  ip,
                  stress,
                  eps_vl,
                  eps_vnl,
                  grad_eps_vnl,
                  D,
                  E,
                  F,
                  istore,
                  newval);
break;
/*----------------------------------------------------------------------*/
default:
  dserror(" unknown type of material law");
break;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return;
} /* end of wge_call_mat */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
