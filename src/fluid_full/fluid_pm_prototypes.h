/*!
\file
\brief fluid_pm_prototypes

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/

#include "../headers/standardtypes.h"
#include "../solver/solver_sparse.h"

/************************************************************************
 | fluid_pm.c                                                           |
 ************************************************************************/
void fluid_pm();


/************************************************************************
 | fluid_pm_service.c                                                   |
 ************************************************************************/

void pm_build_pmat_sparse_mask(FIELD* actfield,
                               PARTITION *actpart,
                               INT disnum,
                               INTRA* actintra,
                               INT numpdof,
                               PARALLEL_SPARSE* pmat);

INT pm_assign_press_dof(FIELD *actfield,
                        PARTITION *actpart,
                        INT disnum,
                        INTRA* actintra);

void pm_fill_gradient_update(PARTITION *actpart,
                             INT disnum,
                             INTRA* actintra,
                             INT numpdof,
                             PARALLEL_SPARSE* grad);

void pm_gradient_mask_mat(FIELD *actfield,
                          PARTITION *actpart,
                          INT disnum,
                          INTRA* actintra,
                          INT numpdof,
                          PARALLEL_SPARSE *grad);

void pm_calelm(FIELD *actfield,
               PARTITION *actpart,
               INT disnum,
               SOLVAR *actsolv,
               INT sysarray,
               INTRA *actintra,
               ARRAY_POSITION *ipos,
               INT numpdof,
               PARALLEL_SPARSE* grad,
               DOUBLE* lmass_vec);

void pm_calprhs(FIELD *actfield,
                PARTITION *actpart,
                INT disnum,
                INTRA *actintra,
                ARRAY_POSITION *ipos,
                INT numpdof,
                DIST_VECTOR* rhs);

void pm_press_update(FIELD *actfield,
                     PARTITION *actpart,
                     INT disnum,
                     INTRA *actintra,
                     ARRAY_POSITION *ipos,
                     INT numpdof,
                     DIST_VECTOR* sol,
                     DOUBLE dta);

void pm_vel_update(FIELD *actfield,
                   PARTITION *actpart,
                   INT disnum,
                   INTRA *actintra,
                   ARRAY_POSITION *ipos,
                   DOUBLE* lmass,
                   DOUBLE* rhs1,
                   DOUBLE* rhs2);

void pm_out_screen_header(INT numeq,
                          INT numeq_total,
                          INTRA *actintra,
                          FILE *out,
                          FLUID_DYNAMIC *fdyn);

/*! @} (documentation module close)*/
