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

void fluid_pm_cont_laplace();

void fluid_pm_cont();

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
			     INT numldof,
                             INT* update);

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
#ifdef PM_TRILINOS
	       TRILINOSMATRIX* grad,
	       TRILINOSMATRIX* lmass
#else
               PARALLEL_SPARSE* grad,
               DOUBLE* lmass_vec
#endif
  );
void pm_calelm_cont(FIELD *actfield,
		    PARTITION *actpart,
		    INT vel_dis,
		    INT press_dis,
		    SOLVAR *actsolv,
		    INT sysarray,
		    INTRA *actintra,
		    ARRAY_POSITION *ipos,
		    TRILINOSMATRIX* grad,
		    TRILINOSMATRIX* lmass
  );
void pm_calelm_laplace(FIELD *actfield,
		       PARTITION *actpart,
		       INT vdisnum,
		       INT pdisnum,
		       SOLVAR *actsolv,
		       INT press_array,
		       INT mass_array,
		       INTRA *actintra,
		       ARRAY_POSITION *ipos
  );

void pm_calprhs(FIELD *actfield,
                PARTITION *actpart,
                INT disnum,
                INTRA *actintra,
                ARRAY_POSITION *ipos,
                INT numpdof,
                DIST_VECTOR* rhs);
void pm_calprhs_cont(FIELD *actfield,
		     PARTITION *actpart,
		     INT disnum,
		     INTRA *actintra,
		     ARRAY_POSITION *ipos,
		     DIST_VECTOR* rhs,
		     DOUBLE* full_rhs,
		     DOUBLE* full_rhs2);
void pm_calvrhs(FIELD *actfield,
		PARTITION *actpart,
		INT disnum,
		INTRA *actintra,
		ARRAY_POSITION *ipos,
		DIST_VECTOR* rhs,
		DOUBLE* full_rhs);

void pm_press_update(FIELD *actfield,
                     PARTITION *actpart,
                     INT disnum,
                     INTRA *actintra,
                     ARRAY_POSITION *ipos,
                     INT numpdof,
                     DIST_VECTOR* sol,
                     DOUBLE thsl);

void pm_vel_update(FIELD *actfield,
                   PARTITION *actpart,
                   INT disnum,
                   INTRA *actintra,
                   ARRAY_POSITION *ipos,
#ifdef PM_TRILINOS
		   TRILINOSMATRIX* lmass,
		   SOLVAR *actsolv,
		   INT actsysarray,
#else
                   DOUBLE* lmass,
#endif
                   DOUBLE* rhs1,
                   DOUBLE* rhs2);

void pm_out_screen_header(INT numeq,
                          INT numeq_total,
                          INTRA *actintra,
                          FILE *out,
                          FLUID_DYNAMIC *fdyn);

/*! @} (documentation module close)*/
