#ifdef MLPCG

/*!---------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

---------------------------------------------------------------------*/

#ifndef SOLVER_MLPCG_PROTOTYPES_H
#define SOLVER_MLPCG_PROTOTYPES_H


void solver_mlpcg(
    struct _SOLVAR         *actsolv,
    struct _INTRA          *actintra,
    struct _DBCSR          *bdcsr,
    struct _DIST_VECTOR    *sol,
    struct _DIST_VECTOR    *rhs,
    INT                     option);

void mlpcg_solver_create(
    struct _DBCSR       *bdcsr, 
    DIST_VECTOR *sol,
    DIST_VECTOR *rhs,
    MLPCGVARS   *mlpcgvars);

void mlpcg_pcg(
    struct _DBCSR       *bdcsr, 
    DIST_VECTOR *sol,
    DIST_VECTOR *rhs,
    INTRA       *actintra);

void mlpcg_solver_init(
    struct _DBCSR       *bdcsr, 
    DIST_VECTOR *sol,
    DIST_VECTOR *rhs,
    INTRA       *actintra);

void mask_bdcsr(
    FIELD         *actfield, 
    PARTITION     *actpart, 
    SOLVAR        *actsolv,
    INTRA         *actintra, 
    struct _DBCSR         *bdcsr);

void bdcsr_make_csr(
    FIELD         *actfield, 
    PARTITION     *actpart, 
    SOLVAR        *actsolv,
    struct _DBCSR         *bdcsr,
    INT          **dof_connect);

void  bdcsr_nnz_topology(
    FIELD         *actfield, 
    PARTITION     *actpart, 
    SOLVAR        *actsolv,
    INTRA         *actintra,
    struct _DBCSR         *bdcsr,
    INT          **dof_connect);

void bdcsr_update(
    FIELD          *actfield, 
    PARTITION      *actpart, 
    SOLVAR         *actsolv,
    INTRA          *actintra,
    struct _DBCSR          *bdcsr);

void bdcsr_numeq(
    FIELD         *actfield, 
    PARTITION     *actpart, 
    SOLVAR        *actsolv,
    INTRA         *actintra,
    INT           *numeq);

void  add_bdcsr(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _DBCSR         *bdcsr1,
    struct _DBCSR         *bdcsr2);

void mlpcg_precond_create(struct _DBCSR     *bdcsr, 
    MLPCGVARS *mlpcgvars,
    INTRA     *actintra);

void mlpcg_precond_init(struct _DBCSR  *bdcsr,MLPCGVARS *mlpcgvars, INTRA *actintra);
void mlpcg_precond_restrictK(MLLEVEL  *actlev, MLLEVEL *nextlev, INTRA *actintra);
void mlpcg_precond_getdirs(void);
void mlpcg_precond_aggsetdofs(MLLEVEL *actlev,INT numdf, INTRA *actintra);
void mlpcg_precond_agg(MLLEVEL *actlev,INTRA *actintra);
void mlpcg_precond_getneightoagg(INT **neighblock, 
                                 INT  *bpatch[],
                                 INT   nbpatch,
                                 INT **freeblock,
                                 INT   nfreeblock,
                                 INT   numeq,
                                 INT  *update,
                                 INT  *ia,
                                 INT  *ja);
void mlpcg_precond_getfreenblocks(INT  *actblock,
                                  INT **freeblock,
                                  INT   nfreeblock,
                                  INT  *bpatch[],
                                  INT  *nbpatch,
                                  INT   numeq,
                                  INT  *update,
                                  INT  *ia,
                                  INT  *ja);
void mlpcg_matvec(DOUBLE       *y, 
                  struct _DBCSR        *A,
                  DOUBLE       *x,
                  DOUBLE        fac,
                  INT           init,
                  INTRA        *actintra);
void mlpcg_matvec_init(struct _DBCSR       *bdcsr, 
                       INTRA       *actintra);
void mlpcg_matvec_uninit(struct _DBCSR       *bdcsr);
void mlpcg_vecvec(DOUBLE *scalar, DOUBLE *x, DOUBLE *y, const INT dim, INTRA *actintra);
INT mlpcg_getindex(INT dof, INT *update, INT length);
INT mlpcg_getowner(INT dof, INT owner[][2], INT nproc);
void mlpcg_precond_P0(MLLEVEL  *actlev, INTRA *actintra);
void mlpcg_precond_P(MLLEVEL  *actlev, INTRA *actintra);
void mlpcg_precond_oneP_fish(AGG     *actagg,
                             DOUBLE   aggblock[][500],
                             INT      rindex[],
                             INT      cindex[],
                             INT     *nrow,
                             INT     *ncol,
                             struct _DBCSR   *actstiff,
                             INTRA   *actintra);
void mlpcg_smoothP(struct _DBCSR *P, struct _DBCSR *actstiff, INTRA *actintra);
void mlpcg_precond_oneP0_vanek(AGG     *actagg,
                         DOUBLE   aggblock[][500],
                         INT      rindex[],
                         INT      cindex[],
                         INT     *nrow,
                         INT     *ncol,
                         struct _DBCSR   *actstiff);
void mlpcg_precond_oneP_vanek(AGG     *actagg,
                        DOUBLE   aggblock[][500],
                        INT      rindex[],
                        INT      cindex[],
                        INT     *nrow,
                        INT     *ncol,
                        struct _DBCSR   *actstiff,
                        MLLEVEL *prevlevel);
void mlpcg_csr_open(struct _DBCSR*  csr,
                    INT     firstdof,
                    INT     lastdof,
                    INT     numeq_total,
                    INT     nnz_guess,
                    INTRA  *actintra);
void mlpcg_renumberdofs(INT            myrank,
                        INT            nproc,
                        FIELD         *actfield, 
                        PARTDISCRET   *actpdiscret, 
                        INTRA         *actintra,
                        struct _DBCSR         *bdcsr,
			INT            dis);
void mlpcg_extractcollocal(struct _DBCSR *P, INT actcol, DOUBLE *col, 
                           INT *rcol, INT *nrow);
void mlpcg_csr_open(struct _DBCSR*  matrix,
                    INT     firstdof,
                    INT     lastdof,
                    INT     numeq_total,
                    INT     nnz_guess,
                    INTRA  *actintra);
void mlpcg_csr_close(struct _DBCSR*   matrix);
void mlpcg_csr_destroy(struct _DBCSR*   matrix);
void mlpcg_csr_setblock(struct _DBCSR*   matrix,
                        DOUBLE   block[][500],
                        INT     *rindex,
                        INT     *cindex,
                        INT      nrow, 
                        INT      ncol,
                        INTRA   *actintra);
void mlpcg_csr_addblock(struct _DBCSR*   matrix,
                        DOUBLE   block[][500],
                        INT     *rindex,
                        INT     *cindex,
                        INT      nrow, 
                        INT      ncol,
                        INTRA   *actintra);
void mlpcg_csr_addentry(struct _DBCSR*   matrix,
                        DOUBLE  val,
                        INT     rindex,
                        INT     cindex,
                        INTRA   *actintra);
void mlpcg_csr_addrow(struct _DBCSR*   matrix,
                      INT       rownum,
                      DOUBLE   *row,
                      INT     *cindex,
                      INT      ncol,
                      INTRA   *actintra);
void mlpcg_csr_csrtocsc(struct _DBCSR *matrix, INTRA *actintra);
void mlpcg_extractcolcsc(INT      col,
                         INT      numeq,
                         INT     *update,
                         INT     *ia,
                         INT     *ja,
                         DOUBLE  *a,
                         DOUBLE **col_out,
                         INT    **rcol_out,
                         INT     *size,
                         INT     *nrow);
void mlpcg_extractrowcsr(INT      row,
                         INT      numeq,
                         INT     *update,
                         INT     *ia,
                         INT     *ja,
                         DOUBLE  *a,
                         DOUBLE **row_out,
                         INT    **rrow_out,
                         INT     *size,
                         INT     *ncol);
void mlpcg_csr_extractsubblock_dense(struct _DBCSR *from, DOUBLE **A,
                                     INT   *index, INT nindex,
                                     INTRA *actintra);
void mlpcg_csr_setentry(struct _DBCSR*   matrix,
                        DOUBLE  val,
                        INT     rindex,
                        INT     cindex,
                        INTRA   *actintra);
void mlpcg_csr_setentry_overlap(struct _DBCSR*   matrix,
                                DOUBLE  val,
                                INT     rindex,
                                INT     cindex,
                                INTRA   *actintra);
void mlpcg_matvec_asm_overlap(DOUBLE       *y, 
                              struct _DBCSR        *A,
                              DOUBLE       *x,
                              DOUBLE        fac,
                              INT           init,
                              INTRA        *actintra);
void mlpcg_csr_getdinv(DOUBLE *Dinv, struct _DBCSR *csr, INT numeq);
void mlpcg_csr_extractsubblock(struct _DBCSR *from, struct _DBCSR *to,
                               INT    rstart,
                               INT    rend,
                               INT    cstart,
                               INT    cend,
                               INTRA *actintra);
void mlpcg_csr_localnumsf(struct _DBCSR *matrix);
void mlpcg_csr_zero(struct _DBCSR*  matrix,
                    INTRA  *actintra);
void mlpcg_extractcollocal_init(struct _DBCSR    *matrix,
                                INT      *sizes,
                                INT    ***icol,
                                DOUBLE ***dcol);
void mlpcg_extractcollocal_fast(struct _DBCSR *matrix, INT actcol, 
                                DOUBLE **col,INT **rcol, INT *size, INT *nrow,
                                INT *sizes, INT ***icol, DOUBLE ***dcol);
void mlpcg_extractcollocal_uninit(struct _DBCSR    *matrix,
                                  INT      *sizes,
                                  INT    ***icol,
                                  DOUBLE ***dcol);
void mlpcg_precond_PtKP(struct _DBCSR *P, struct _DBCSR *incsr, struct _DBCSR *outcsr, struct _DBCSR *work,
                        AGG *agg, INT nagg, INTRA *actintra);
void mlpcg_precond_presmo(DOUBLE *z, DOUBLE *r, struct _DBCSR *csr, MLLEVEL *lev, INTRA *actintra, INT level);
void mlpcg_precond_postsmo(DOUBLE *z, DOUBLE *r, struct _DBCSR *csr, MLLEVEL *lev, INTRA *actintra, INT level);
void mlpcg_precond_coarsesolv(DOUBLE *z, DOUBLE *r, MLLEVEL *lev, INTRA *actintra);
void mlpcg_precond_prolongz(DOUBLE *zc, DOUBLE *z, struct _DBCSR *P, struct _DBCSR *coarsecsr, 
                            INTRA *actintra);
void mlpcg_precond_restrictr(DOUBLE *rc, DOUBLE *r, struct _DBCSR *P, struct _DBCSR *coarsecsr, 
                            INTRA *actintra);
void mlpcg_precond_check_fcd(struct _DBCSR *matrix, INTRA *actintra);
void mlpcg_precond_checkdirich(struct _DBCSR *matrix, INTRA *actintra);
void mlpcg_precond_smoJacobi(DOUBLE *z, DOUBLE *r, struct _DBCSR *csr, INT nsweep, INTRA *actintra);
void mlpcg_precond_lapacksolve(DOUBLE *z, DOUBLE *r, struct _DBCSR *csr, INTRA *actintra);
void mlpcg_precond_smo_ILUn(DOUBLE *z, DOUBLE *r, struct _DBCSR *csr, INT nsweep, INTRA *actintra);
void mlpcg_precond_smo_ILUn_overlap(DOUBLE *z, DOUBLE *r, struct _DBCSR *csr, INT nsweep, INTRA *actintra);
void mlpcg_precond_gramschmidt(DOUBLE **P, DOUBLE **R,const INT nrow,const INT ncol);
void mlpcg_precond_spoolessolve(DOUBLE *z, DOUBLE *r, struct _DBCSR *csr, INTRA *actintra);
void mlpcg_precond_amgVW(INT    level,
                         ARRAY *z_a,
                         ARRAY *r_a,
                         INTRA *actintra,
                         INT   *gamma);
void mlpcg_printvec(INT          iter,
                   DOUBLE      *z, 
                   struct _DBCSR*       csr,
                   DISCRET     *fielddis,
                   PARTDISCRET *partdis,
                   INTRA       *actintra);
void mlpcg_csr_overlap(struct _DBCSR *csr, struct _DBCSR *ocsr, struct _DBCSR *ilu, INT overlap, INTRA *actintra);
void mlpcg_csr_localnumsf_overlap(struct _DBCSR *matrix);
void mlpcg_precond_P_fish(MLLEVEL  *actlev, INTRA *actintra);

#endif

#endif /* MLPCG */
