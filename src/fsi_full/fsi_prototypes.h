/*!----------------------------------------------------------------------
\file
\brief fsi_prototypes

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FSI
*//*! @{ (documentation module open)*/

#include "fsi.h"

/* RULE HOW TO ADD NEW FILES AND FUNCTIONS:
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE!!!
*/
/************************************************************************
 | fsi_aitken.c                                                         |
 ************************************************************************/
void fsi_aitken(
    FIELD              *structfield,
    INT                 disnum,
    INT                 itnum,
    INT                 init
    );

void fsi_aitken_force(
  FIELD             *structfield,
  INT                sdisnum,
  FIELD             *fluidfield,
  INT                fdisnum,
  INT                itnum,
  INT                numff,
  INT                init
  );

/************************************************************************
 | fsi_ale.c                                                            |
 ************************************************************************/
void fsi_ale_setup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_calc(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  );

void fsi_ale_final(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_sd(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  );

void fsi_ale_output(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_cleanup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );


/************************************************************************
 | fsi_ale_2step.c                                                            |
 ************************************************************************/
void fsi_ale_2step_setup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_2step_calc(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_2step_final(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_2step_sd(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  );

void fsi_ale_2step_output(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_2step_cleanup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );




/************************************************************************
 | fsi_ale_laplace.c                                                            |
 ************************************************************************/
void fsi_ale_laplace_setup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_laplace_calc(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_laplace_final(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_laplace_sd(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  );

void fsi_ale_laplace_output(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_laplace_cleanup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );


/************************************************************************
 | fsi_ale_LAS.c                                                            |
 ************************************************************************/
void fsi_ale_LAS_setup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_LAS_calc(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_LAS_final(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_LAS_sd(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  );

void fsi_ale_LAS_output(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_LAS_cleanup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );


/************************************************************************
 | fsi_ale_lin.c                                                            |
 ************************************************************************/
void fsi_ale_lin_setup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_lin_calc(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  );

void fsi_ale_lin_final(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_lin_sd(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  );

void fsi_ale_lin_output(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_lin_cleanup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );



/************************************************************************
 | fsi_ale_nln.c                                                            |
 ************************************************************************/
void fsi_ale_nln_setup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_nln_calc(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_nln_final(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_nln_sd(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  );

void fsi_ale_nln_output(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_nln_cleanup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );


/************************************************************************
 | fsi_ale_spring.c                                                            |
 ************************************************************************/
void fsi_ale_spring_setup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_spring_calc(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_spring_final(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_spring_sd(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  );

void fsi_ale_spring_output(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );

void fsi_ale_spring_cleanup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  );


/************************************************************************
 | fsi_coupforce.c                                                      |
 ************************************************************************/
void fsi_cbf(
    PARTDISCRET        *actpdis,
    DOUBLE             *fcouple,
    ARRAY_POSITION     *ipos,
    INT                 numeq_total,
    INT                 init
    );


#ifdef PARALLEL
void fsi_allreduce_coupforce(
    DOUBLE             *fcouple,
    DOUBLE             *recvfcouple,
    INT                 numeq_total,
    INT                 numddof,
    INTRA              *actintra,
    FIELD              *actfield,
    INT                 disnum
    );
#endif


void fsi_load(
    PARTITION          *actpart,
    INT                 disnum,
    FIELD              *fluidfield,
    INT                 fdisnum,
    DOUBLE             *fsiforce,
    INT                 global_numeq
    );


/************************************************************************
 | fsi_coupling.c                                                       |
 ************************************************************************/
void fsi_createfsicoup(void);


void fsi_initcoupling(
    FIELD              *structfield,
    INT                 disnum_s,
    FIELD              *fluidfield,
    INT                 disnum_f,
    FIELD              *alefield,
    INT                 disnum_a
    );


void fsi_struct_intdofs(
    FIELD              *structfield,
    INT                 disnum
    );


/************************************************************************
 | fsi_dyn.c                                                            |
 ************************************************************************/
void dyn_fsi(
    INT                 mctrl
    );


/************************************************************************
 | fsi_energy.c                                                          |
 ************************************************************************/
void fsi_dyneint(
    FIELD              *structfield,
    INT                 disnum,
    INT                 init
    );


void fsi_energycheck( void );


/************************************************************************
 | fsi_fluid.c                                                          |
 ************************************************************************/
void fsi_fluid_setup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_calc(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io,
  FIELD          *alefield,
  INT             adisnum_calc
  );

void fsi_fluid_final(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_sd(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_output(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_cleanup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );



/************************************************************************
 | fsi_fluid_imp.c                                                      |
 ************************************************************************/
void fsi_fluid_imp_setup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_imp_calc(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io,
  FIELD          *alefield,
  INT             adisnum_calc
  );

void fsi_fluid_imp_final(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_imp_sd(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_imp_output(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_imp_cleanup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );



/************************************************************************
 | fsi_fluid_pm_cont.c                                                  |
 ************************************************************************/
void fsi_fluid_pm_cont_setup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_pm_cont_calc(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io,
  FIELD          *alefield,
  INT             adisnum_calc
  );

void fsi_fluid_pm_cont_final(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_pm_cont_sd(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_pm_cont_output(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_pm_cont_cleanup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );



/************************************************************************
 | fsi_fluid_pm_discont.c                                                  |
 ************************************************************************/
void fsi_fluid_pm_discont_setup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_pm_discont_calc(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io,
  FIELD          *alefield,
  INT             adisnum_calc
  );

void fsi_fluid_pm_discont_final(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_pm_discont_sd(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_pm_discont_output(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_pm_discont_cleanup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );



/************************************************************************
 | fsi_fluid_pm_laplace.c                                               |
 ************************************************************************/
void fsi_fluid_pm_laplace_setup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_pm_laplace_calc(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io,
  FIELD          *alefield,
  INT             adisnum_calc
  );

void fsi_fluid_pm_laplace_final(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_pm_laplace_sd(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_pm_laplace_output(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );

void fsi_fluid_pm_laplace_cleanup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  );



/************************************************************************
 | fsi_gradient.c                                                       |
 ************************************************************************/
void fsi_gradient(
  FSI_STRUCT_WORK* struct_work,
  FSI_FLUID_WORK* fluid_work,
  FSI_ALE_WORK* ale_work,
    FIELD              *alefield,
    FIELD              *structfield,
    FIELD              *fluidfield,
    INT                 disnuma_io,
    INT                 disnuma_calc,
    INT                 disnums_io,
    INT                 disnums_calc,
    INT                 disnumf_io,
    INT                 disnumf_calc,
    INT                 numfa,
    INT                 numff,
    INT                 numfs
    );


void fsi_gradient_force(
  FSI_STRUCT_WORK* struct_work,
  FSI_FLUID_WORK* fluid_work,
  FSI_ALE_WORK* ale_work,
  FIELD              *alefield,
  FIELD              *structfield,
  FIELD              *fluidfield,
  INT                 disnuma_io,
  INT                 disnuma_calc,
  INT                 disnums_io,
  INT                 disnums_calc,
  INT                 disnumf_io,
  INT                 disnumf_calc,
  INT                 numfa,
  INT                 numff,
  INT                 numfs
  );


#ifdef FSI_NEWTONCOUPLING

/************************************************************************
 | fsi_newton.c                                                       |
 ************************************************************************/
double *fsi_newton_SolveSystem(
    DOUBLE             *rS,
    FIELD              *structfield,
    FIELD              *fluidfield,
    FIELD              *alefield,
    FSI_STRUCT_WORK    *struct_work,
    FSI_FLUID_WORK     *fluid_work,
    FSI_ALE_WORK       *ale_work,
    INT                s_disnum_calc,
    INT                s_disnum_io,
    INT                f_disnum_calc,
    INT                f_disnum_io,
    INT                a_disnum_calc,
    INT                a_disnum_io,
    INT                itnum
  );

double *fsi_newton_rightSide(
  FIELD           *structfield,
  INT             s_disnum_calc
  );

void fsi_newton_final(
  DOUBLE          *dg,
  FIELD           *structfield,
  INT             s_disnum_calc
  );

#endif

/************************************************************************
 | fsi_mortar.c                                                         |
 ************************************************************************/
#ifdef D_MORTAR
void fsi_initcoupling_intfaces(
    FIELD              *masterfield,
    INT                 m_disnum,
    FIELD              *slavefield,
    INT                 s_disnum,
    INTERFACES         *int_faces
    );


void fsi_init_interfaces(
    FIELD              *masterfield,
    INT                 m_disnum,
    FIELD              *slavefield,
    INT                 s_disnum,
    INTERFACES         *int_faces
    );


void fsi_mortar_coeff(
    FSI_DYNAMIC        *fsidyn,
    INTERFACES         *int_faces
    );


void fsi_detect_intersection(
    DOUBLE              lambdr1_2,
    DOUBLE              lambdr2_1,
    DOUBLE              lambdr2_2,
    DOUBLE              nr1_2,
    DOUBLE              nr2_1,
    DOUBLE              nr2_2,
    DOUBLE             *b1,
    DOUBLE             *b2,
    INT                *intersection
    );


void fsi_calc_disp4ale(
    FSI_DYNAMIC        *fsidyn,
    INTERFACES         *int_faces
    );


void fsi_calc_intforces(
    INTERFACES         *int_faces
    );


void f2_fsiload(
    ELEMENT            *ele
    );


void fsi_put_coupforc2struct(
    FIELD              *masterfield,
    INT                 disnum,
    INTERFACES         *int_faces
    );


void fsiserv_rhs_point_neum(
    DOUBLE             *rhs,
    INT                 dimrhs,
    PARTITION          *actpart
    );
#endif


/************************************************************************
 | fsi_relax_intdisp.c                                                  |
 ************************************************************************/
void fsi_relax_intdisp(
    FIELD              *structfield,
    INT                 disnum
    );

void fsi_relax_intdisp_force(
  FIELD              *structfield,
  INT                 sdisnum,
  FIELD              *fluidfield,
  INT                 fdisnum,
  INT                 numff
  );


/************************************************************************
 | fsi_service.c                                                        |
 ************************************************************************/
void fsi_alecp(
    FIELD              *fluidfield,
    INT                 disnum,
    FIELD              *alefield,
    INT                 adisnum,
    DOUBLE              dt,
    INT                 numdf,
    INT                 phase
    );


void fsi_aleconv(
    FIELD              *fluidfield,
    INT                 disnum,
    INT                 numdf,
    INT                 pos1,
    INT                 pos2
    );


void fsi_fluidstress_result(
    FIELD           *actfield,
    INT              disnum,
    INT              numdf
    );


void fsi_algoout(
    INT                 itnum
    );


void fsi_structpredictor(
    FIELD              *actfield,
    INT                 disnum,
    INT                 init
    );


INT fsi_convcheck(
    FIELD              *structfield,
    INT                 disnum,
    INT                 itnum,
    DOUBLE             *resnorm
    );


INT fsi_convcheck_force(
  FIELD              *structfield,
  INT                 sdisnum,
  FIELD              *fluidfield,
  INT                 fdisnum,
  INT                 itnum,
  INT                 numff
  );


void fsi_init_ale(
    FIELD              *actfield,
    INT                 numr
    );


void fluid_init_pos_ale(
    FIELD              *fluidfield,
    INT                 disnum
    );

/************************************************************************
 | fsi_struct.c                                                         |
 ************************************************************************/
void fsi_struct_setup(
  FSI_STRUCT_WORK    *work,
  FIELD              *actfield,
  INT                 disnum_calc,
  INT                 disnum_io,
  INT                 fsiitnum
  );

void fsi_struct_calc(
  FSI_STRUCT_WORK    *work,
  FIELD              *actfield,
  INT                 disnum_calc,
  INT                 disnum_io,
  INT                 fsiitnum,
  FIELD              *fluidfield,
  INT                 fdisnum_calc
  );

void fsi_struct_final(
  FSI_STRUCT_WORK    *work,
  FIELD              *actfield,
  INT                 disnum_calc,
  INT                 disnum_io,
  INT                 fsiitnum
  );

void fsi_struct_sd(
  FSI_STRUCT_WORK    *work,
  FIELD              *actfield,
  INT                 disnum_calc,
  INT                 disnum_io,
  INT                 fsiitnum
  );

void fsi_struct_output(
  FSI_STRUCT_WORK    *work,
  FIELD              *actfield,
  INT                 disnum_calc,
  INT                 disnum_io,
  INT                 fsiitnum
  );

void fsi_struct_cleanup(
  FSI_STRUCT_WORK    *work,
  FIELD              *actfield,
  INT                 disnum_calc,
  INT                 disnum_io,
  INT                 fsiitnum
  );

/*! @} (documentation module close)*/
