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


/************************************************************************
 | fsi_ale.c                                                            |
 ************************************************************************/
void fsi_ale(
    FIELD              *actfield,
    INT                 disnum_calc,
    INT                 disnum_io,
    INT                 mctrl
    );


/************************************************************************
 | fsi_ale_2step.c                                                            |
 ************************************************************************/
void fsi_ale_2step(
    FIELD              *actfield,
    INT                 disnum_calc,
    INT                 disnum_io,
    INT                 mctrl
    );


/************************************************************************
 | fsi_ale_laplace.c                                                            |
 ************************************************************************/
void fsi_ale_laplace(
    FIELD              *actfield,
    INT                 disnum_calc,
    INT                 disnum_io,
    INT                 mctrl
    );


/************************************************************************
 | fsi_ale_LAS.c                                                            |
 ************************************************************************/
void fsi_ale_LAS(
    FIELD              *actfield,
    INT                 disnum_calc,
    INT                 disnum_io,
    INT                 mctrl
    );


/************************************************************************
 | fsi_ale_lin.c                                                            |
 ************************************************************************/
void fsi_ale_lin(
    FIELD              *actfield,
    INT                 disnum_calc,
    INT                 disnum_io,
    INT                 mctrl
    );


/************************************************************************
 | fsi_ale_nln.c                                                            |
 ************************************************************************/
void fsi_ale_nln(
    FIELD              *actfield,
    INT                 disnum_calc,
    INT                 disnum_io,
    INT                 mctrl
    );


/************************************************************************
 | fsi_ale_spring.c                                                            |
 ************************************************************************/
void fsi_ale_spring(
    FIELD              *actfield,
    INT                 disnum_calc,
    INT                 disnum_io,
    INT                 mctrl
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
void fsi_fluid(
    FIELD              *actfield,
    INT                 disnum_calc,
    INT                 disnum_io,
    INT                 mctrl
    );


/************************************************************************
 | fsi_gradient.c                                                       |
 ************************************************************************/
void fsi_gradient(
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


/************************************************************************
 | fsi_service.c                                                        |
 ************************************************************************/
void fsi_alecp(
    FIELD              *fluidfield,
    INT                 disnum,
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
    INT                 itnum
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
void fsi_struct(
    FIELD              *actfield,
    INT                 disnum_calc,
    INT                 disnum_io,
    INT                 mctrl,
    INT                 fsiitnum
    );


/*! @} (documentation module close)*/
