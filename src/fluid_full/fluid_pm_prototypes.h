/*!----------------------------------------------------------------------
\file
\brief fluid_pm_prototypes

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
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
/************************************************************************
 | fluid_pm.c                                                           |
 ************************************************************************/
void fluid_pm(void);

/************************************************************************
 | fluid_pm_calamatrx.c                                                 |
 ************************************************************************/
void fluid_pm_calamatrix(
                          FIELD           *actfield,
			  INTRA           *actintra,
			  PARTITION       *actpart,
			  SOLVAR          *actsolv,
			  SPARSE_ARRAY    *amatrix_oll,
			  FLUID_DYNAMIC   *fdyn,
			  OLL             *gradmatrix_oll,
			  OLL             *rgradmatrix_oll,
			  DOUBLE          *lmass,
			  INT             *numeq_total,
			  INT             *numeq,
			  INT             *numeq_full,
	 		  INT             *numeq_total_full,
			  CALC_ACTION     *action,
			  CONTAINER       *container
			 );
void fluid_pm_dofrenumber(FIELD      *actfield,
#ifdef PARALLEL
                          INT        *saveveldof,
			  INT        *savepredof,
			  INT        *newveldof,
			  INT        *newpredof,
                          INT        *velpointer,
			  INT        *prepointer,
#endif
			  INTRA      *actintra,
			  PARTITION  *actpart,
			  DBCSR      *amatrix_csr,
			  DBCSR      *lumpedmass,
			  INT        *numeq_total_full,
			  INT        *numeq_full,
			  INT         numnp_vel,
			  INT        *firstdof,
			  INT        *lastdof,
			  INT         option
			);
void assemble_fluid_amatrix(
                             CONTAINER *container,
                             ELEMENT   *actvele,
                             ELEMENT   *actpele,
			     INTRA     *actintra
			   );
void fluid_pm_lumpedmass(DBCSR *lumpedmass_csr,INTRA *actintra,
                         DOUBLE *lmass,
#ifdef PARALLEL
			 INT *velpointer, INT *saveveldof,
#endif
			 INT *numeq_total);
void fluid_pm_redcpmat(
#ifdef PARALLEL
		       INT *saverowdof,      INT *rowpointer,
                       INT *savecoldof,      INT *colpointer,
#endif
		       INT numeq_totalr,     INT numeq_totalc,
		       INT numeq_totalr_full,INT numeq_totalc_full,
		       DBCSR *matrix_csr,    OLL *matrix_oll);
void fluid_pm_newdofs(INT            myrank,
                          INT            nproc,
                          FIELD         *actfield,
                          PARTDISCRET   *actpdiscret,
                          INTRA         *actintra,
                          DBCSR         *bdcsr,
			  INT            dis);

/************************************************************************
 | fluid_pm_service.c                                                    |
 ************************************************************************/
void fluid_pm_matlmatmul(
                        SPARSE_ARRAY   *m,
		        SPARSE_TYP     *m_typ,
                        DOUBLE      *lmass
		        );
void fluid_pm_matvecmul(ARRAY *rc_a,     ARRAY *r_a, OLL *P_oll,
                        INTRA *actintra, INT numeq,  INT option);
void fluid_pm_lmatmulvec(
                         DOUBLE   *lmat,
		         DOUBLE   *vec,
		         INT       numeq_total
			);
void fluid_pm_fullvel(FIELD *actfield, INT disnum,
                      DOUBLE *fullvel, INT place);
void fluid_pm_pretovel(FIELD *actfield,INT actpos);
void fluid_init_pos_pm(ARRAY_POSITION *ipos);

/*! @} (documentation module close)*/
