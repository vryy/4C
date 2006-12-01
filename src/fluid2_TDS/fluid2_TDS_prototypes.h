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
    DOUBLE           visc,
    DOUBLE         **wa1,
    DOUBLE         **wa2,
    DOUBLE           estress[3][MAXNOD_F2],
    INT              _relax
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
		DOUBLE   old_vel[2],
		DOUBLE   res_old[2]
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


/*! @} (documentation module close)*/
