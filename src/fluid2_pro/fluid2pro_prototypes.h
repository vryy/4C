/*!----------------------------------------------------------------------
\file
\brief fluid2_pro prototypes

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/

#ifdef D_FLUID2_PRO

/************************************************************************
| f2pro_calint.c                                                        |
************************************************************************/
void f2pro_calint(
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

/************************************************************************
 | f2pro_calele.c                                                       |
************************************************************************/
void f2pro_calinit(
  ARRAY          *estif_global,
  ARRAY          *emass_global,
  ARRAY          *lmass_global,
  ARRAY          *gradopr_global,
  ARRAY          *eforce_global,
  ARRAY          *edforce_global,
  ARRAY          *gforce_global,
  ARRAY_POSITION *ipos
  );

void f2pro_calele(
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

/************************************************************************
 | f2pro_calgalmat.c                                                    |
************************************************************************/
void f2pro_calkvv(
  DOUBLE         **estif,
  DOUBLE         **derxy,
  DOUBLE           fac,
  DOUBLE           visc,
  DOUBLE           dt,
  INT                iel
  );
void f2pro_lmass(
  DOUBLE         **lmass,
  DOUBLE         **emass,
  INT                iel
  );
void f2pro_gradopr(
  DOUBLE         **gradopr,
  DOUBLE          **derxy,
  DOUBLE         *functpr,
  DOUBLE              fac,
  INT                  ielp,
  INT                   iel
  );

/************************************************************************
 | f2pro_calservice.c                                                   |
************************************************************************/
void f2pro_calset(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **eveln,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE          *epren,
  ARRAY_POSITION  *ipos
  );

/************************************************************************
| f2pro_caltimerhs.c                                                    |
************************************************************************/
void f2pro_calgaltfv(
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
 | f2pro_inpele.c                                                       |
************************************************************************/
void f2pro_inp(ELEMENT *ele);

/************************************************************************
 | f2pro_main.c                                                         |
************************************************************************/
void fluid2_pro(     PARTITION     *actpart,
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

void f2pro_prec(DOUBLE* pfunct,
                DOUBLE** pderiv,
                DOUBLE r,
                DOUBLE s,
                DISMODE dm,
                INT* numpdof);

void f2pro_int_usfem(
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

void f2pro_calmat( DOUBLE **estif,
                   DOUBLE  *eforce,
                   DOUBLE  *velint,
                   DOUBLE   histvec[2],
                   DOUBLE   gridvint[2],
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
                   INT     *hasext
  );

void f2pro_calgradp(
  ELEMENT* ele,
  ARRAY_POSITION *ipos
  );

void f2pro_calprhs(
  ELEMENT* ele,
  ARRAY_POSITION *ipos
  );

void f2pro_calvelupdate(
  ELEMENT* ele,
  ARRAY_POSITION *ipos
  );

void f2pro_addnodepressure(ELEMENT* ele, INT k, DOUBLE* pressure);

#if defined(DEBUG) && !defined(PARALLEL)

void f2pro_debugoutpressure(ELEMENT* ele, FILE* f);

#endif

#endif
