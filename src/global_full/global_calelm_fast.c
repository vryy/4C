/*-----------------------------------------------------------------------*/
/*!
\file
\brief contains the routine 'calelm_fast' to loop all sets of fast elements


<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

/*!
\addtogroup GLOBAL
*//*! @{ (documentation module open)*/

#ifdef D_FLUID3_F


#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../fluid3/fluid3.h"
#include "../fluid3_fast/f3f_prototypes.h"

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01  |
  | general problem data                                               |
  | global variable GENPROB genprob is defined in global_control.c     |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
  | enum _CALC_ACTION                                      m.gee 1/02  |
  | command passed from control routine to the element level           |
  | to tell element routines what to do                                |
  | defined globally in global_calelm.c                                |
 *----------------------------------------------------------------------*/
enum _CALC_ACTION calc_action[MAXFIELD];


/*----------------------------------------------------------------------*
  | global dense matrices for element routines             m.gee 7/01  |
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_fast;   /* element stiffness matrix(fortran) */
extern struct _ARRAY emass_fast;   /* element mass matrix (fortran) */
extern struct _ARRAY eforce_fast;  /* element Iteration RHS(fortran) */
extern struct _ARRAY edforce_fast; /* element dirichlet RHS(fortran) */


/*-----------------------------------------------------------------------*/
/*!
  \brief calls all sets of fast elements sets

  This routine loops all sets of fast elements and calculates the matrices of
  all contained elements.

  \param actfield      *FIELD       (i) active field
  \param actsolv       *SOLVAR      (i) active SOLVAR
  \param actpart       *PARTITION   (i) my partition of this field
  \param actintra      *INTRA       (i) my intra-communicator
  \param sysarray1      INT         (i) number of first sparse system matrix
  \param sysarray2      INT         (i) number of secnd system matrix, if present, else -1
  \param container     *CONTAINER   (i) contains variables defined in container.h
  \param action        *CALC_ACTION (i) calculation option passed to element routines

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void calelm_fast(
    FIELD        *actfield,     /* active field */
    SOLVAR       *actsolv,      /* active SOLVAR */
    PARTITION    *actpart,      /* my partition of this field */
    INTRA        *actintra,     /* my intra-communicator */
    INT           sysarray1,    /* number of first sparse system matrix */
    INT           sysarray2,    /* number of secnd system matrix, if present, else -1 */
    CONTAINER    *container,    /* contains variables defined in container.h */
    CALC_ACTION  *action)       /* calculation option passed to element routines */

{
  INT               i,kk;
  ASSEMBLE_ACTION   assemble_action;

  FAST_ELES        *act_fast_eles;


  /*variables for for fluid_fast and shell_fast*/
  INT               hasext_fast[LOOPL];
  INT               hasdirich_fast[LOOPL];
  INT               iel;

#ifdef DEBUG
dstrc_enter("calelm_fast");
#endif



  for(i=0;i<LOOPL;i++)
  {
    hasext_fast[i]=0;
    hasdirich_fast[i]=0;
  }


  /* call elements */
  kk = container->actndis;


  /* loop all vectors of fast elements */
  for (i=0; i<actpart->pdis[kk].num_fele; i++)
  {

    act_fast_eles = &(actpart->pdis[kk].fast_eles[i]);

    switch(act_fast_eles->fast_ele_typ)
    {
      case fele_f3f_hex8_e:
      case fele_f3f_hex8_a:
      case fele_f3f_hex20_e:
      case fele_f3f_hex20_a:
      case fele_f3f_tet4_e:
      case fele_f3f_tet4_a:


    /* perform integration for a group of elements*/
    fluid3_fast(
        actpart,
        actintra,
        act_fast_eles->ele_vec,
        &estif_fast,
        &emass_fast,
        &eforce_fast,
        &edforce_fast,
        action,
        hasdirich_fast,
        hasext_fast,
        container,
        act_fast_eles->aloopl);




    if (*action != calc_fluid_error)
    {
      assemble_action = assemble_one_matrix;

      iel=act_fast_eles->ele_vec[0]->numnp;

      assemble_fast(
          sysarray1,
          actpart,
          actsolv,
          actintra,
          act_fast_eles->ele_vec,
          assemble_action,
          container,
          iel,
          hasext_fast,
          hasdirich_fast,
          act_fast_eles->aloopl);

    }  /* if (*action != fluid_cal_error) */

        break;


      default:
        dserror("unknown typ of fast element");

    } /* switch(act_fast_eles->fast_ele_typ) */

  } /* for (i=0; i<actpart->pdis[kk].num_fele; i++) */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of calelm_fast */


#endif /* ifdef D_FLUID3_F */

/*! @} (documentation module close)*/


