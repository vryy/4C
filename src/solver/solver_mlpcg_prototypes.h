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
void mlpcg_precond_create(DBCSR     *bdcsr, 
                          MLPCGVARS *mlpcgvars,
                          INTRA     *actintra);
void mlpcg_precond_init(DBCSR  *bdcsr,MLPCGVARS *mlpcgvars, INTRA *actintra);
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
                  DBCSR        *A,
                  DOUBLE       *x,
                  DOUBLE        fac,
                  INT           init,
                  INTRA        *actintra);
void mlpcg_matvec_init(DBCSR       *bdcsr, 
                       INTRA       *actintra);
void mlpcg_matvec_uninit(DBCSR       *bdcsr);
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
                             DBCSR   *actstiff,
                             INTRA   *actintra);
void mlpcg_smoothP(DBCSR *P, DBCSR *actstiff, INTRA *actintra);
void mlpcg_precond_oneP0_vanek(AGG     *actagg,
                         DOUBLE   aggblock[][500],
                         INT      rindex[],
                         INT      cindex[],
                         INT     *nrow,
                         INT     *ncol,
                         DBCSR   *actstiff);
void mlpcg_precond_oneP_vanek(AGG     *actagg,
                        DOUBLE   aggblock[][500],
                        INT      rindex[],
                        INT      cindex[],
                        INT     *nrow,
                        INT     *ncol,
                        DBCSR   *actstiff,
                        MLLEVEL *prevlevel);
void mlpcg_csr_open(DBCSR*  csr,
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
                        DBCSR         *bdcsr,
			INT            dis);
void mlpcg_extractcollocal(DBCSR *P, INT actcol, DOUBLE *col, 
                           INT *rcol, INT *nrow);
void mlpcg_csr_open(DBCSR*  matrix,
                    INT     firstdof,
                    INT     lastdof,
                    INT     numeq_total,
                    INT     nnz_guess,
                    INTRA  *actintra);
void mlpcg_csr_close(DBCSR*   matrix);
void mlpcg_csr_destroy(DBCSR*   matrix);
void mlpcg_csr_setblock(DBCSR*   matrix,
                        DOUBLE   block[][500],
                        INT     *rindex,
                        INT     *cindex,
                        INT      nrow, 
                        INT      ncol,
                        INTRA   *actintra);
void mlpcg_csr_addblock(DBCSR*   matrix,
                        DOUBLE   block[][500],
                        INT     *rindex,
                        INT     *cindex,
                        INT      nrow, 
                        INT      ncol,
                        INTRA   *actintra);
void mlpcg_csr_addentry(DBCSR*   matrix,
                        DOUBLE  val,
                        INT     rindex,
                        INT     cindex,
                        INTRA   *actintra);
void mlpcg_csr_addrow(DBCSR*   matrix,
                      INT       rownum,
                      DOUBLE   *row,
                      INT     *cindex,
                      INT      ncol,
                      INTRA   *actintra);
void mlpcg_csr_csrtocsc(DBCSR *matrix, INTRA *actintra);
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
void mlpcg_csr_extractsubblock_dense(DBCSR *from, DOUBLE **A,
                                     INT   *index, INT nindex,
                                     INTRA *actintra);
void mlpcg_csr_setentry(DBCSR*   matrix,
                        DOUBLE  val,
                        INT     rindex,
                        INT     cindex,
                        INTRA   *actintra);
void mlpcg_csr_setentry_overlap(DBCSR*   matrix,
                                DOUBLE  val,
                                INT     rindex,
                                INT     cindex,
                                INTRA   *actintra);
void mlpcg_matvec_asm_overlap(DOUBLE       *y, 
                              DBCSR        *A,
                              DOUBLE       *x,
                              DOUBLE        fac,
                              INT           init,
                              INTRA        *actintra);
void mlpcg_csr_getdinv(DOUBLE *Dinv, DBCSR *csr, INT numeq);
void mlpcg_csr_extractsubblock(DBCSR *from, DBCSR *to,
                               INT    rstart,
                               INT    rend,
                               INT    cstart,
                               INT    cend,
                               INTRA *actintra);
void mlpcg_csr_localnumsf(DBCSR *matrix);
void mlpcg_csr_zero(DBCSR*  matrix,
                    INTRA  *actintra);
void mlpcg_extractcollocal_init(DBCSR    *matrix,
                                INT      *sizes,
                                INT    ***icol,
                                DOUBLE ***dcol);
void mlpcg_extractcollocal_fast(DBCSR *matrix, INT actcol, 
                                DOUBLE **col,INT **rcol, INT *size, INT *nrow,
                                INT *sizes, INT ***icol, DOUBLE ***dcol);
void mlpcg_extractcollocal_uninit(DBCSR    *matrix,
                                  INT      *sizes,
                                  INT    ***icol,
                                  DOUBLE ***dcol);
void mlpcg_precond_PtKP(DBCSR *P, DBCSR *incsr, DBCSR *outcsr, DBCSR *work,
                        AGG *agg, INT nagg, INTRA *actintra);
void mlpcg_precond_presmo(DOUBLE *z, DOUBLE *r, DBCSR *csr, MLLEVEL *lev, INTRA *actintra, INT level);
void mlpcg_precond_postsmo(DOUBLE *z, DOUBLE *r, DBCSR *csr, MLLEVEL *lev, INTRA *actintra, INT level);
void mlpcg_precond_coarsesolv(DOUBLE *z, DOUBLE *r, MLLEVEL *lev, INTRA *actintra);
void mlpcg_precond_prolongz(DOUBLE *zc, DOUBLE *z, DBCSR *P, DBCSR *coarsecsr, 
                            INTRA *actintra);
void mlpcg_precond_restrictr(DOUBLE *rc, DOUBLE *r, DBCSR *P, DBCSR *coarsecsr, 
                            INTRA *actintra);
void mlpcg_precond_check_fcd(DBCSR *matrix, INTRA *actintra);
void mlpcg_precond_checkdirich(DBCSR *matrix, INTRA *actintra);
void mlpcg_precond_smoJacobi(DOUBLE *z, DOUBLE *r, DBCSR *csr, INT nsweep, INTRA *actintra);
void mlpcg_precond_lapacksolve(DOUBLE *z, DOUBLE *r, DBCSR *csr, INTRA *actintra);
void mlpcg_precond_smo_ILUn(DOUBLE *z, DOUBLE *r, DBCSR *csr, INT nsweep, INTRA *actintra);
void mlpcg_precond_smo_ILUn_overlap(DOUBLE *z, DOUBLE *r, DBCSR *csr, INT nsweep, INTRA *actintra);
void mlpcg_precond_gramschmidt(DOUBLE **P, DOUBLE **R,const INT nrow,const INT ncol);
void mlpcg_precond_spoolessolve(DOUBLE *z, DOUBLE *r, DBCSR *csr, INTRA *actintra);
void mlpcg_precond_amgVW(INT    level,
                         ARRAY *z_a,
                         ARRAY *r_a,
                         INTRA *actintra,
                         INT   *gamma);
void mlpcg_printvec(INT          iter,
                   DOUBLE      *z, 
                   DBCSR*       csr,
                   DISCRET     *fielddis,
                   PARTDISCRET *partdis,
                   INTRA       *actintra);
void mlpcg_csr_overlap(DBCSR *csr, DBCSR *ocsr, DBCSR *ilu, INT overlap, INTRA *actintra);
void mlpcg_csr_localnumsf_overlap(DBCSR *matrix);
void mlpcg_precond_P_fish(MLLEVEL  *actlev, INTRA *actintra);

