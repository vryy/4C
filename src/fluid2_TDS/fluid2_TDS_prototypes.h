/*!----------------------------------------------------------------------
\file
\brief integration loop for one fluid2 element using USFEM

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/gammi/
            +49-(0)89-289-15235
</pre>


------------------------------------------------------------------------*/
/*!
\addtogroup FLUID2
*/

/*! @{ (documentation module open)*/


/************************************************************************
 | f2_int_tds.c                                                       |
 ************************************************************************/
void f2_int_tds(
    ELEMENT         *ele,
    INT             *hasext,
    DOUBLE         **estif,
    DOUBLE          *eforce,
    DOUBLE         **xyze,
    DOUBLE          *funct,
    DOUBLE         **deriv,
    DOUBLE         **deriv2,
    DOUBLE         **xjm,
    DOUBLE         **derxy,
    DOUBLE         **derxy2,
    DOUBLE         **evelng,
    DOUBLE         **eveln,
    DOUBLE         **evhist,
    DOUBLE         **egridv,
    DOUBLE          *epreng,
    DOUBLE          *epren,
    DOUBLE          *edeadng,
    DOUBLE          *edeadn,
    DOUBLE         **vderxy,
    DOUBLE         **vderxy2,
    DOUBLE         **vderxy_old,
    DOUBLE         **vderxy2_old,
    DOUBLE         **eacc,
    DOUBLE           visc,
    DOUBLE         **wa1,
    DOUBLE         **wa2,
    DOUBLE           estress[3][MAXNOD_F2],
    INT              _relax
    )
    ;

void f2_int_gen_alpha_tds(
                      ELEMENT         *ele,
                      INT             *hasext,
                      DOUBLE         **estif,
	              DOUBLE          *eforce,
	              DOUBLE         **xyze,
	              DOUBLE          *funct,
	              DOUBLE         **deriv,
	              DOUBLE         **deriv2,
	              DOUBLE         **xjm,
	              DOUBLE         **derxy,
	              DOUBLE         **derxy2,
		      DOUBLE	     **eaccng,
	              DOUBLE         **evelng,
		      DOUBLE          *epreng, 
	              DOUBLE          *edeadng,
		      DOUBLE         **vderxy,
                      DOUBLE         **vderxy2,
		      DOUBLE           visc,
	              DOUBLE         **wa1,
	              DOUBLE         **wa2
	             )
    ;

/************************************************************************
 | f2_calmat_tds.c                                                       |
 ************************************************************************/
void f2_calmat_tds(
                DOUBLE **estif,
		DOUBLE  *eforce,
		DOUBLE  *velint,
		DOUBLE   histvec[2],
		DOUBLE   gridvint[2],
		DOUBLE   press,
		DOUBLE **vderxy,
                DOUBLE **vderxy2,
                DOUBLE   gradp[2],
		DOUBLE  *funct,
		DOUBLE **derxy,
		DOUBLE **derxy2,
                DOUBLE  *edeadng,
		DOUBLE   fac,
		DOUBLE   visc,
		INT      iel,
                INT     *hasext,
                INT      isale,
                INT      is_relax,
		DOUBLE   sub_pres,
		DOUBLE   divu_old,
		DOUBLE   sub_vel[2],
		DOUBLE   sub_vel_trial_wo_facMtau[2],
		DOUBLE   old_vel[2],
		DOUBLE   old_acc[2],
		DOUBLE   res_old[2]
              )
    ;

/************************************************************************
 | f2_calmat_incr_acc_gen_alpha_tds.c                                   |
 ************************************************************************/

void f2_calgalmat_gen_alpha_tds(
                DOUBLE **estif,
		DOUBLE  *velint,
		DOUBLE  *funct,
		DOUBLE **derxy,
		DOUBLE **derxy2,
		DOUBLE   fac,
		DOUBLE   visc,
		int      iel
              )
    ;

void f2_calstabmat_gen_alpha_tds(
                DOUBLE **estif,
		DOUBLE  *velint,
		DOUBLE  *funct,
		DOUBLE **derxy,
		DOUBLE **derxy2,
		DOUBLE   fac,
		DOUBLE   visc,
		int      iel
              )
    ;

/************************************************************************
 | f2_update_tds.c                                                       |
 ************************************************************************/
void f2_update_subscale_pres(
    PARTITION      *actpart,
    INTRA          *actintra,
    FIELD          *actfield,
    ARRAY_POSITION *ipos,
    INT             disnum_calc
    );

void f2_update_subscale_vel(
    PARTITION      *actpart,
    INTRA          *actintra,
    FIELD          *actfield,
    ARRAY_POSITION *ipos,
    INT             disnum_calc
    );

void f2_update_subscale_pres_for_inc_gen_alpha(
    PARTITION      *actpart,
    INTRA          *actintra,
    FIELD          *actfield,
    ARRAY_POSITION *ipos,
    INT             disnum_calc
    );

/************************************************************************
 | f2_caltau_tds.c                                                       |
 ************************************************************************/
void f2_get_time_dependent_sub_tau(ELEMENT *ele,
				   DOUBLE  **xyze,
				   DOUBLE   *funct,
				   DOUBLE  **deriv,
				   DOUBLE  **evelng,
				   DOUBLE  **eveln,
				   DOUBLE    visc
    );


/************************************************************************
 | f2_calservice_for_TDS.c                                              |
 ************************************************************************/

void f2_inc_gen_alpha_calset(
	        ELEMENT         *ele,
                DOUBLE         **xyze,
                DOUBLE         **eaccng,
	        DOUBLE         **evelng,
		DOUBLE          *epreng,
		DOUBLE          *edeadng,
                ARRAY_POSITION *ipos,
		double         *visc
	      )
    ;

void fluid_result_incre_for_genalpha(
                          FIELD             *actfield,
                          INT                disnum,
                          INTRA             *actintra,
			  DIST_VECTOR       *sol,
			  DIST_VECTOR       *rhs,
                          ARRAY_POSITION    *ipos,
			  SPARSE_ARRAY      *sysarray,
			  SPARSE_TYP        *sysarray_typ,
			  DOUBLE            *vrat,
			  DOUBLE            *prat,
                          DOUBLE            *grat
		       )
    ;

/************************************************************************
 | f2_calrhs_incr_acc_gen_alpha_tds.c                                   |
 ************************************************************************/
void f2_calgalrhs_gen_alpha_tds(
                DOUBLE  *eforce,
		DOUBLE  *velint,
		DOUBLE  *accint,
		DOUBLE   presint,
		DOUBLE  *edeadng,
		DOUBLE  *funct,
		DOUBLE **derxy,
		DOUBLE **derxy2,
		DOUBLE **vderxy,
                DOUBLE **vderxy2,
		DOUBLE   fac,
		DOUBLE   visc,
		int      iel     
              )
    ;


void f2_calstabrhs_gen_alpha_tds(
                DOUBLE  *eforce,
		DOUBLE  *velint,
		DOUBLE  *accint,
		DOUBLE   presint,
		DOUBLE  *gradpint,
		DOUBLE  *edeadng,
		DOUBLE  *funct,
		DOUBLE **derxy,
		DOUBLE **derxy2,
		DOUBLE **vderxy,
                DOUBLE **vderxy2,
		DOUBLE   fac,
		DOUBLE   visc,
		int      iel     
              )
    ;

/*! @} (documentation module close)*/
