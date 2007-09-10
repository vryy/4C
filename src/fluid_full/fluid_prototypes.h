/*!----------------------------------------------------------------------
\file
\brief fluid_prototypes

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/* RULE HOW TO ADD NEW FILES AND FUNCTIONS:
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE!!!
*/

#ifndef CCADISCRET

/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
/************************************************************************
 | fluid_adapt_service.c                                                |
 ************************************************************************/
void fluid_acceleration(
    FIELD              *actfield,
    INT                 disnum,
    ARRAY_POSITION     *ipos,
    INT                 iop
    );


void fluid_prep_rhs(
    FIELD              *actfield,
    INT                 disnum,
    ARRAY_POSITION     *ipos
    );


void fluid_predictor(
    FIELD              *actfield,
    INT                 disnum,
    ARRAY_POSITION     *ipos,
    INT                 iop
    );


void fluid_lte(
    FIELD              *actfield,
    INT                 disnum,
    ARRAY_POSITION     *ipos,
    INT                 iop
    );


void fluid_lte_norm(
    PARTITION          *actpart,
    INT                 disnum,
    INTRA              *actintra,
    ARRAY_POSITION     *ipos,
    INT                *iststep,
    INT                *repeat,
    INT                *repeated,
    INT                 itnum
    );


/************************************************************************
 | fluid_consistbf.c                                                    |
 ************************************************************************/
void fluid_cbf(
    PARTDISCRET        *actpdis,
    CONTAINER          *container,
    ARRAY_POSITION     *ipos,
    INT                 init
    );


/************************************************************************
 | fluid_curvature.c                                                    |
 ************************************************************************/
void fluid_curvature(
    FIELD              *actfield,
    PARTITION          *actpart,
    INTRA              *actintra,
    CALC_ACTION        *action
    );


/************************************************************************
 | fluid_dirich.c                                                       |
 ************************************************************************/
void fluid_initdirich(
    FIELD              *actfield,
    INT                 disnum,
    ARRAY_POSITION     *ipos
    );


void fluid_setdirich(
    FIELD              *actfield,
    INT                 disnum,
    ARRAY_POSITION     *ipos,
    INT                 pos
    );


void fluid_setdirich_cyl(
    FIELD              *actfield,
    INT                 disnum,
    ARRAY_POSITION     *ipos
    );


void fluid_setdirich_parabolic(
    FIELD              *actfield,
    INT                 disnum,
    ARRAY_POSITION     *ipos
    );


void fluid_setdirich_sd(
    FIELD              *actfield,
    INT                 disnum,
    ARRAY_POSITION     *ipos
    );


void fluid_caldirich(
    ELEMENT            *actele,
    DOUBLE             *dforces,
    DOUBLE            **estif,
    INT                *hasdirich,
    INT                 readfrom
    );


void fluid_reaction_forces(ELEMENT* ele, FLUID_DYNAMIC* fdyn, DOUBLE** estif, DOUBLE* eforce, INT pos);


void fluid_pm_caldirich(
    ELEMENT            *actele,
    DOUBLE             *dforces,
    DOUBLE            **estif,
    DOUBLE            **emass,
    DOUBLE              dt,
    DOUBLE              theta,
    ARRAY_POSITION     *ipos,
    INT                *hasdirich
    );


void fluid_pm_caldirich_parabolic(
    ELEMENT            *actele,
    DOUBLE             *dforces,
    DOUBLE            **estif,
    DOUBLE            **emass,
    DOUBLE              dt,
    DOUBLE              theta,
    INT                *hasdirich
    );


void fluid_pm_caldirich_cyl(
    ELEMENT            *actele,
    DOUBLE             *dforces,
    DOUBLE            **estif,
    DOUBLE            **emass,
    DOUBLE              dt,
    DOUBLE              theta,
    INT                *hasdirich
    );


/************************************************************************
 | fluid_dirich_tu.c                                                     |
 ************************************************************************/
void fluid_setdirich_tu(
    FIELD              *actfield,
    DOUBLE             *lower_limit_kappa,
    DOUBLE             *lower_limit_eps,
    ARRAY_POSITION     *ipos
    );


void fluid_caldirich_tu(
    ELEMENT            *actele,
    DOUBLE             *dforces,
    DOUBLE            **estif,
    ARRAY_POSITION     *ipos,
    INT                *hasdirich
    );


/************************************************************************
 | fluid_dirich_tu_1.c                                                   |
 ************************************************************************/
void fluid_setdirich_tu_1(
    FIELD              *actfield,
    DOUBLE             *lower_limit_kappa,
    DOUBLE             *lower_limit_omega,
    ARRAY_POSITION     *ipos
    );


void fluid_caldirich_tu_1(
    ELEMENT            *actele,
    DOUBLE             *dforces,
    DOUBLE            **estif,
    ARRAY_POSITION     *ipos,
    INT                *hasdirich
    );


/************************************************************************
 | fluid_dyn.c                                                          |
 ************************************************************************/
void dyn_fluid(void);


/************************************************************************
 | fluid_freesurface.c                                                  |
 ************************************************************************/
void fluid_createfreesurf(void);


void fluid_freesurf_setdofs(void);


void fluid_modcoor(void);


void fluid_updfscoor(
    FIELD              *fluidfield,
    FLUID_DYNAMIC      *fdyn,
    DOUBLE              dt,
    ARRAY_POSITION     *ipos,
    INT                 flag
    );


/************************************************************************
 | fluid_heightfunction.c                                               |
 ************************************************************************/
void fluid_heightfunc(
    INT                 hctrl,
    DOUBLE             *grat,
    FIELD              *actfield,
    PARTITION          *actpart,
    INTRA              *actintra,
    CALC_ACTION        *action,
    CONTAINER          *container,
    ARRAY_POSITION     *ipos
    );


/************************************************************************
 | fluid_imp_semimp.c                                                   |
 ************************************************************************/
void fluid_isi(void);


/************************************************************************
 | fluid_incr_acc_gen_alpha.c                                                |
 ************************************************************************/
#ifdef D_FLUID2_TDS
void fluid_incr_acc_gen_alpha(void);
#endif /* D_FLUID2_TDS */


/************************************************************************
 | fluid_imp_semimp_tu.c                                                |
 ************************************************************************/
void fluid_isi_tu(void);


/************************************************************************
 | fluid_imp_semimp_tu_1.c                                              |
 ************************************************************************/
void fluid_isi_tu_1(void);


/************************************************************************
 | fluid_liftdrag.c                                                     |
 ************************************************************************/
void fluid_liftdrag(
    INT                 init,
    CALC_ACTION        *action,
    CONTAINER          *container,
    FIELD              *actfield,
    INT                 disnum,
    SOLVAR             *actsolv,
    PARTITION          *actpart,
    INTRA              *actintra,
    ARRAY_POSITION     *ipos
    );


/************************************************************************
 | fluid_locsys.c                                                       |
 ************************************************************************/
void fluid_locsys(
    FIELD              *actfield,
    INT                 disnum,
    FLUID_DYNAMIC      *fdyn
    );


/************************************************************************
 | fluid_mf.c                                                           |
 ************************************************************************/
void fluid_mf(
    INT                 mctrl
    );


/************************************************************************
 | fluid_mfcoupling.c                                                   |
 ************************************************************************/
void fluid_initmfcoupling(
    FIELD              *fluidfield,
    FIELD              *alefield
    );


/************************************************************************
 | fluid_mlservice.c                                                   |
 ************************************************************************/
void fluid_ml_init(
    FIELD              *actfield
    );


void fluid_smcopy(
    PARTITION          *actpart
    );


void fluid_prgmr(
    DOUBLE            **smmat,
    DOUBLE            **smrhs,
    INT                 numeq,
    INT                 numrhs
    );


void fluid_add_smat(
    DOUBLE            **smat,
    DOUBLE            **semat,
    INT                 numen,
    INT                *slme,
    DOUBLE              fac
    );


void fluid_add_smrhs(
    DOUBLE            **smrhs,
    DOUBLE            **smerhs,
    INT                 numrhs,
    INT                 numen,
    INT                *smlme
    );


void fluid_add_ssrhs(
    DOUBLE             *ssrhs,
    DOUBLE             *sserhs,
    INT                 numen,
    INT                *sslme
    );


void fluid_add_intlhs(
    DOUBLE             *smidiff,
    DOUBLE            **smiediff,
    INT                 numen,
    INT                *smlme
    );


void fluid_add_intrhs(
    DOUBLE             *smirhs,
    DOUBLE             *smierhs,
    INT                 numen,
    INT                *smlme
    );


void fluid_bubint(
    DOUBLE             *bubint,
    DOUBLE             *smfunct,
    DOUBLE            **ebub,
    INT                 smiel,
    INT                 nbub
    );


void fluid_pbubint(
    DOUBLE            **pbubint,
    DOUBLE             *smfunct,
    DOUBLE            **epbub,
    INT                 smiel,
    INT                 iel,
    INT                 nsd
    );


void fluid_bubder(
    DOUBLE            **bubderxy,
    DOUBLE            **smderxy,
    DOUBLE            **ebub,
    INT                 smiel,
    INT                 nbub,
    INT                 nder
    );


void fluid_pbubder(
    DOUBLE           ***pbubderxy,
    DOUBLE            **smderxy,
    DOUBLE            **epbub,
    INT                 smiel,
    INT                 iel,
    INT                 nsd,
    INT                 nder
    );


void fluid_calnofo(
    DOUBLE             *erhs,
    DOUBLE             *funct,
    DOUBLE              fac,
    INT                 iel
    );


void fluid_mlcaldirich(
    ELEMENT            *actele,
    DOUBLE             *dforces,
    DOUBLE            **estif,
    INT                *hasdirich,
    INT                 readfrom
    );


/************************************************************************
 | fluid_normal.c                                                       |
 ************************************************************************/
void fluid_cal_normal(
    FIELD              *actfield,
    INT                 disnum,
    INT                 init,
    CALC_ACTION        *action
    );


/************************************************************************
 | fluid_pseudofsi.c                                                    |
 ************************************************************************/
void fluid_createpseudofsi(void);


/************************************************************************
 | fluid_service.c                                                      |
 ************************************************************************/
void fluid_startproc(
    INT                *nfrastep,
    INT                 init
    );


void fluid_cons( void );


void fluid_tcons( void );


void fluid_icons(
    INT                 itnum
    );


void fluid_init(
    PARTITION          *actpart,
    INTRA              *actintra,
    FIELD              *actfield,
    INT                 disnum_calc,
    INT                 disnum_io,
    CALC_ACTION        *action,
    CONTAINER          *container,
    INT                 numr,
    ARRAY_POSITION     *ipos,
    FLUID_STRESS        str
    );


void fluid_norm(
    FIELD              *actfield,
    INT                 disnum,
    INT                 numeq_total,
    DOUBLE             *vrat,
    DOUBLE             *prat,
    ARRAY_POSITION     *ipos
    );


void fluid_sol_copy(
    FIELD              *actfield,
    INT                 disnum,
    INT                 arrayfrom,
    INT                 arrayto,
    INT                 from,
    INT                 to,
    INT                 numdf
    );


void fluid_transpres(
    FIELD              *actfield,
    INT                 disnum,
    INT                 index,
    INT                 actpos,
    INT                 predof,
    INT                 option
    );


INT fluid_steadycheck(
    FIELD              *actfield,
    INT                 disnum,
    ARRAY_POSITION     *ipos,
    INT                 numeq_total
    );


INT fluid_convcheck(
    DOUBLE              vrat,
    DOUBLE              prat,
    DOUBLE              grat,
    INT                 itnum,
    DOUBLE              te,
    DOUBLE              ts
    );


void fluid_algoout(void);


void fluid_reduceshstr(
    INTRA              *actintra,
    FIELD              *actfield
    );


void fluid_nullshstr(
    INTRA              *actintra,
    PARTITION          *actpart,
    FIELD              *actfield
    );


void fluid_reducestress(
    INTRA              *actintra,
    PARTITION          *actpart,
    FIELD              *actfield,
    INT                 disnum,
    INT                 numdf,
    FLUID_STRESS        str
    );


void fluid_cal_error(
    FIELD              *actfield,
    INT                 disnum,
    ARRAY_POSITION     *ipos,
    INT                 index
    );


void fluid_init_pos_euler(
    ARRAY_POSITION     *ipos
    );


/************************************************************************
 | fluid_service_tu.c                                                   |
 ************************************************************************/
void fluid_init_tu(
    FIELD              *actfield
    );


void fluid_tcons_tu(void);


void fluid_set_check_tu(
    FIELD              *actfield,
    ARRAY_POSITION     *ipos,
    DOUBLE              lower_limit_kappa
    );


void fluid_eddy_pro(
    FIELD              *actfield,
    ARRAY_POSITION     *ipos
    );


INT fluid_convcheck_tu(
    DOUBLE              kapepsrat,
    INT                 itnum1,
    DOUBLE              te,
    DOUBLE              ts
    );


INT fluid_convcheck_test(
    FIELD              *actfield,
    ARRAY_POSITION     *ipos,
    INT                 itnum_check
    );


void fluid_icons_tu(
    INT                 itnum1,
    INT                 itnumke,
    INT                 itnum_n
    );


void fluid_copysol_tu(
    FIELD              *actfield,
    INT                 from,
    INT                 to,
    INT                 flag
    );

void fluid_copysol_test(
    FIELD              *actfield,
    INT                 from,
    INT                 to
    );


void fluid_algoout_tu( void );


void fluid_init_pos_euler_tu(
    ARRAY_POSITION     *ipos
    );


/************************************************************************
 | fluid_service_tu_1.c                                                 |
 ************************************************************************/
void fluid_set_check_tu_1(
    FIELD              *actfield,
    DOUBLE              lower_limit_kappa,
    DOUBLE              lower_limit_omega,
    ARRAY_POSITION     *ipos
    );


/************************************************************************
 | fluid_stationary.c                                                   |
 ************************************************************************/
void fluid_stat(void);


/************************************************************************
 | inp_fluid_start_data.c                                               |
 ************************************************************************/
void inp_fluid_start_data(
    FIELD              *actfield,
    ARRAY_POSITION     *ipos,
    FLUID_DYNAMIC      *fdyn
    );

/*! @} (documentation module close)*/

#endif
