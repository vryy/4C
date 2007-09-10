/*!----------------------------------------------------------------------
\file
\brief fluid3_pro prototypes

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/

#ifndef CCADISCRET
#ifdef D_FLUID3_PRO

/************************************************************************
| f3pro_calint.c                                                        |
************************************************************************/
void f3pro_calint(
  ELEMENT         *elev,
  ELEMENT         *elep,
  DOUBLE         **estif,
  DOUBLE         **emass,
  DOUBLE         **gradopr,
  DOUBLE         *etforce,
  DOUBLE         *eiforce,
  DOUBLE         **xyze,
  DOUBLE         *funct,
  DOUBLE         *functpr,
  DOUBLE         **deriv,
  DOUBLE         **derivpr,
  DOUBLE         **xjm,
  DOUBLE         **derxy,
  DOUBLE         **derxypr,
  DOUBLE         **eveln,
  DOUBLE         *epren,
  DOUBLE         *velint,
  DOUBLE         *covint,
  DOUBLE         **vderxy,
  DOUBLE         *pderxy,
  DOUBLE         **wa1,
  DOUBLE         *dirich,
  DOUBLE         **deriv2,
  INT              *dirich_onoff
  );

void f3pro_int_res(
  ELEMENT         *ele,
  INT             *hasext,
  DOUBLE          *force,
  DOUBLE         **xyze,
  DOUBLE          *funct,
  DOUBLE         **deriv,
  DOUBLE         **deriv2,
  DOUBLE          *pfunct,
  DOUBLE         **pderiv,
  DOUBLE         **xjm,
  DOUBLE         **derxy,
  DOUBLE         **derxy2,
  DOUBLE         **pderxy,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE         **egridv,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  DOUBLE         **vderxy,
  DOUBLE         **vderxy2,
  DOUBLE           visc,
  DOUBLE         **wa1,
  DOUBLE         **wa2
  );

void f3pro_calresvec(  DOUBLE  *eforce,
                       DOUBLE  *velint,
                       DOUBLE   histvec[3],
                       DOUBLE   gridvint[3],
                       DOUBLE **vderxy,
                       DOUBLE **vderxy2,
                       DOUBLE  *funct,
                       DOUBLE **derxy,
                       DOUBLE **derxy2,
                       DOUBLE  *edeadng,
                       DOUBLE   press,
                       DOUBLE   gradp[3],
                       DOUBLE   fac,
                       DOUBLE   visc,
                       INT      iel,
                       INT     *hasext,
                       INT      is_ale
  );

/************************************************************************
 | f3pro_calele.c                                                       |
************************************************************************/
void f3pro_calinit(
  ARRAY          *estif_global,
  ARRAY          *emass_global,
  ARRAY          *lmass_global,
  ARRAY          *gradopr_global,
  ARRAY          *eforce_global,
  ARRAY          *edforce_global,
  ARRAY          *gforce_global,
  ARRAY_POSITION *ipos
  );

void f3pro_calele(
  ELEMENT        *ele,
  ARRAY          *estif_global,
  ARRAY          *emass_global,
  ARRAY          *lmass_global,
  ARRAY          *gradopr_global,
  ARRAY          *eforce_global,
  ARRAY          *edforce_global,
  ARRAY          *gforce_global,
  ARRAY_POSITION *ipos,
  INT            *hasdirich,
  INT            *hasext
  );

void f3pro_caleleres(
  ELEMENT        *ele,
  ARRAY          *eforce_global,
  ARRAY_POSITION *ipos,
  INT            *hasdirich,
  INT            *hasext
  );

/************************************************************************
 | f3pro_calgalmat.c                                                    |
************************************************************************/
void f3pro_calkvv(
  DOUBLE         **estif,
  DOUBLE         **derxy,
  DOUBLE           fac,
  DOUBLE           visc,
  DOUBLE           dt,
  INT                iel
  );
void f3pro_lmass(
  DOUBLE         **lmass,
  DOUBLE         **emass,
  INT                iel
  );
void f3pro_gradopr(
  DOUBLE         **gradopr,
  DOUBLE          **derxy,
  DOUBLE         *functpr,
  DOUBLE              fac,
  INT                  ielp,
  INT                   iel
  );

/************************************************************************
 | f3pro_calservice.c                                                   |
************************************************************************/
void f3pro_calset(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **evelnm,
  DOUBLE         **eveln,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE          *eprenm,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  ARRAY_POSITION  *ipos
  );
void f3pro_calseta(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **evelnm,
  DOUBLE         **eveln,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE         **egridv,
  DOUBLE          *eprenm,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  ARRAY_POSITION  *ipos
  );

/************************************************************************
| f3pro_caltimerhs.c                                                    |
************************************************************************/
void f3pro_calgaltfv(
  DOUBLE          *etforce,
  DOUBLE          *eiforce,
  DOUBLE          *velint,
  DOUBLE          *covint,
  DOUBLE         **vderxy,
  DOUBLE          *funct,
  DOUBLE         **derxy,
  DOUBLE           preint,
  DOUBLE           visc,
  DOUBLE           fac,
  DOUBLE           dt,
  INT               iel
  );

/************************************************************************
 | f3pro_inpele.c                                                       |
************************************************************************/
void f3pro_inp(ELEMENT *ele);
void f3pro_dis(ELEMENT* vele, ELEMENT* pele, INT nele, INT nodeshift);

/************************************************************************
 | f3pro_main.c                                                         |
************************************************************************/
void fluid3_pro(     PARTITION     *actpart,
                     INTRA         *actintra,
		     ELEMENT       *ele,
		     ARRAY         *estif_global,
		     ARRAY         *emass_global,
		     ARRAY         *lmass_global,
		     ARRAY         *gradopr_global,
		     ARRAY         *eforce_global,
		     ARRAY         *edforce_global,
		     ARRAY         *gforce_global,
		     CALC_ACTION   *action,
                     INT           *hasdirich,
                     INT           *hasext
  );

void f3pro_phex(DOUBLE* pfunct,
                DOUBLE** pderiv,
                DOUBLE r,
                DOUBLE s,
                DOUBLE t,
                DISMODE dm,
                INT* numpdof);

void f3pro_int_usfem(
  ELEMENT         *ele,
  INT             *hasext,
  DOUBLE         **estif,
  DOUBLE          *eforce,
  DOUBLE          *gforce,
  DOUBLE         **xyze,
  DOUBLE          *funct,
  DOUBLE         **deriv,
  DOUBLE         **deriv2,
  DOUBLE          *pfunct,
  DOUBLE         **pderiv,
  DOUBLE         **xjm,
  DOUBLE         **derxy,
  DOUBLE         **derxy2,
  DOUBLE         **pderxy,
  DOUBLE         **evelng,
  DOUBLE         **eveln,
  DOUBLE         **evelnm,
  DOUBLE         **evhist,
  DOUBLE         **egridv,
  DOUBLE          *eprenm,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  DOUBLE         **vderxy,
  DOUBLE         **vderxy2,
  DOUBLE         **vderxy_n,
  DOUBLE         **vderxy_nm,
  DOUBLE           visc,
  DOUBLE         **wa1,
  DOUBLE         **wa2,
  INT              is_relax
  );

void f3pro_calmat( DOUBLE **estif,
                   DOUBLE  *eforce,
		   DOUBLE  *velint_n,
		   DOUBLE **vderxy_n,
		   DOUBLE  *velint_nm,
		   DOUBLE **vderxy_nm,
                   DOUBLE  *velint,
                   DOUBLE   histvec[3],
                   DOUBLE   gridvint[3],
                   DOUBLE   press_n,
                   DOUBLE   press,
                   DOUBLE **vderxy,
                   DOUBLE **vderxy2,
                   DOUBLE   gradp[3],
                   DOUBLE  *funct,
                   DOUBLE **derxy,
                   DOUBLE **derxy2,
                   DOUBLE  *edeadng,
                   DOUBLE   fac,
                   DOUBLE   visc,
                   INT      iel,
                   INT     *hasext,
                   INT      isale,
                   INT      is_relax
  );

void f3pro_calgradp(
  ELEMENT* ele,
  ARRAY_POSITION *ipos
  );

void f3pro_calprhs(
  ELEMENT* ele,
  ARRAY_POSITION *ipos
  );

void f3pro_calvrhs(ELEMENT* ele, ARRAY_POSITION *ipos);

void f3pro_calvelupdate(
  ELEMENT* ele,
  ARRAY_POSITION *ipos
  );

void f3pro_addnodepressure(ELEMENT* ele, INT k, DOUBLE* pressure);

void f3pro_calpress(ELEMENT* ele,
		    ARRAY* estif_global,
		    ARRAY* emass_global,
		    ARRAY_POSITION* ipos);


#endif
#endif
