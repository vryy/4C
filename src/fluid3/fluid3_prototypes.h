/*!---------------------------------------------------------------------
\file
\brief fluid3 prototypes

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>
#include "../fluid3/fluid3.h"
#include "../fluid2/fluid2.h"

------------------------------------------------------------------------*/
/* RULE HOW TO ADD NEW FILES AND FUNCTIONS:
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE!!!
*/
/************************************************************************
 | f3_calele.c                                                          |
 ************************************************************************/

void f3_calele(
	        ELEMENT        *ele,
                ARRAY          *estif_global,
                ARRAY          *emass_global,
	        ARRAY          *eforce_global,
		ARRAY          *edforce_global,
                ARRAY_POSITION *ipos,
		INT            *hasdirich,
		INT            *hasext,
		INT             init
	       );
void f3_stress(FLUID_STRESS  str,
               INT           viscstr,
	       ELEMENT      *ele,
               ARRAY_POSITION *ipos);
void f3_heightfunc(
                   ELEMENT              *ele,
                   ARRAY                *estif_global,
		   ARRAY                *eiforce_global,
		   CONTAINER            *container,
                   ARRAY_POSITION       *ipos
		   );
void f3_calstab(ELEMENT *ele, ARRAY_POSITION *ipos);
void f3_caleleres(
	           ELEMENT          *ele,
	           ARRAY            *eforce_global,
                   INT              *hasdirich,
                   INT              *hasext,
                   ARRAY_POSITION   *ipos
	       );

/************************************************************************
 | f3_calelesize.c                                                      |
 ************************************************************************/
void f3_calelesize(
                     ELEMENT         *ele,
                     DOUBLE         **xyze,
                     DOUBLE          *funct,
                     DOUBLE         **deriv,
                     DOUBLE         **deriv2,
                     DOUBLE         **derxy,
                     DOUBLE         **xjm,
                     DOUBLE         **evel,
                     DOUBLE         **wa1,
                     INT              cpele
                  );
void f3_calelesize2(
                     ELEMENT         *ele,
                     DOUBLE          *velint,
                     DOUBLE         **derxy,
                     DOUBLE           visc,
                     INT              iel,
                     DIS_TYP          typ
                   );

/************************************************************************
 | f3_calfuncderiv.c                                                    |
 ************************************************************************/
void f3_hex(
               DOUBLE     *funct,
               DOUBLE    **deriv,
               DOUBLE    **deriv2,
               DOUBLE      r,
               DOUBLE      s,
               DOUBLE      t,
               DIS_TYP     typ,
               INT         icode
               );
void f3_rec(
            DOUBLE     *funct,
            DOUBLE    **deriv,
            DOUBLE    **deriv2,
            DOUBLE      r,
            DOUBLE      s,
            DIS_TYP     typ,
            INT         icode
            );
void f3_tet(
            DOUBLE     *funct,
            DOUBLE    **deriv,
            DOUBLE    **deriv2,
            DOUBLE      r,
            DOUBLE      s,
            DOUBLE      t,
            DIS_TYP     typ,
            INT         icode
            );
void f3_tri(
            DOUBLE     *funct,
            DOUBLE    **deriv,
            DOUBLE    **deriv2,
            DOUBLE      r,
            DOUBLE      s,
            DIS_TYP     typ,
            INT         icode
	   );
void f3_jaco(DOUBLE    **xyze,
             DOUBLE    **deriv,
             DOUBLE    **xjm,
             DOUBLE     *det,
             ELEMENT    *ele,
             INT         iel);
void f3_edgejaco( DOUBLE    **xyze,
                  DOUBLE    **deriv,
                  DOUBLE    **xjm,
                  DOUBLE     *det,
                  INT        *iedgnod,
                  INT         iel,
                  ELEMENT    *ele
               );
void f3_tvmr(DOUBLE   **x,
                DOUBLE     akov[3][3],
                DOUBLE    *funct,
                DOUBLE   **deriv,
                INT       *iedgnod,
                INT        iel);
void f3_gder(
               DOUBLE   **derxy,
               DOUBLE   **deriv,
               DOUBLE   **xjm,
               DOUBLE   **xji,
               DOUBLE     det,
               INT        iel
            );
void f3_edgegder(
               DOUBLE   **derxy,
               DOUBLE   **deriv,
               DOUBLE   **xjm,
               DOUBLE     det,
               INT        iel
	    );
void f3_gcoor(
               DOUBLE     *funct,
               ELEMENT    *ele,
               INT         iel,
               DOUBLE     *gcoor
             );
void f3_gder2(
               DOUBLE     **xyze,
               DOUBLE     **xjm,
               DOUBLE     **bm,
               DOUBLE     **xder2,
               DOUBLE     **derxy,
               DOUBLE     **derxy2,
               DOUBLE     **deriv2,
               INT          iel
	     );

/************************************************************************
 | f3_calgalmat.c                                                       |
 ************************************************************************/
void f3_calkvv(
                  ELEMENT         *ele,
                  DOUBLE         **estif,
                  DOUBLE          *velint,
                  DOUBLE          *gridvint,
                  DOUBLE         **vderxy,
                  DOUBLE          *funct,
                  DOUBLE         **derxy,
                  DOUBLE           fac,
                  DOUBLE           visc,
                  INT              iel
              );
void f3_calkvp(
               DOUBLE         **estif,
               DOUBLE          *funct,
               DOUBLE         **derxy,
               DOUBLE           fac,
               INT              iel
              );
void f3_calmvv(
               DOUBLE         **estif,
               DOUBLE          *funct,
               DOUBLE           fac,
               INT              iel
              );
void f3_calkgedge(
                  DOUBLE         **estif,
                  DOUBLE          *funct,
                  DOUBLE           fac,
                  INT             *iedgnod,
                  INT              iel,
                  INT              ngnode
                );

/************************************************************************
 | f3_calheightfunc_sep.c                                               |
 ************************************************************************/
void f3_calint_hfsep(
                     ELEMENT           *ele,
                     DOUBLE            *funct,
                     DOUBLE           **deriv,
                     DOUBLE           **xjm,
                     DOUBLE           **wa2,
                     DOUBLE           **xyze,
                     INT                ngnode,
                     INT                nil,
                     INT               *iedgnod,
                     DOUBLE            *velint,
                     DOUBLE            *vel2int,
                     DOUBLE           **evelng,
                     DOUBLE           **eveln,
                     DOUBLE            *ephing,
                     DOUBLE            *ephin,
                     DOUBLE           **derxy,
                     DIS_TYP            typ,
                     DOUBLE           **estif,
                     DOUBLE            *eiforce
                     );
void f3_stabpar_hfsep(
                        ELEMENT          *ele,
                        DOUBLE          **deriv,
                        DOUBLE          **xjm,
                        DOUBLE          **xyze,
                        DOUBLE           *velint,
                        DOUBLE            phiintn,
                        DOUBLE            phiintng,
                        DOUBLE           *phiderxy,
                        INT              *iedgnod,
                        INT               ngnode,
                        DIS_TYP           typ
                     );
void f3_calmat_vhf_sep(
		       ELEMENT                *ele,
		       DOUBLE                **estif,
		       INT                     ngnode,
		       DOUBLE                 *funct,
		       DOUBLE                **derxy,
		       DOUBLE                 *velint,
		       DOUBLE                  fac
	              );
void f3_caliterhs_vhf_sep(
 	 	           ELEMENT                *ele,
		           DOUBLE                 *eforce,
		           INT                     ngnode,
		           DOUBLE                 *funct,
                           DOUBLE                **derxy,
		           DOUBLE                 *velint,
			   DOUBLE                  fac
	                  );
void f3_caltimerhs_vhf_sep(
                           ELEMENT                *ele,
                           DOUBLE                 *eforce,
                           INT                     ngnode,
                           DOUBLE                 *funct,
                           DOUBLE                **derxy,
                           DOUBLE                 *velint, /* at n+1 */
                           DOUBLE                 *vel2int,/* at n */
                           DOUBLE                  phiint,
                           DOUBLE                 *phiderxy,
                           DOUBLE                  fac
                          );

/************************************************************************
 | f3_calint.c                                                          |
 ************************************************************************/
void f3_calint(
               ELEMENT         *ele,
               DOUBLE         **estif,
               DOUBLE         **emass,
               DOUBLE          *eforce,
               DOUBLE         **xyze,
               DOUBLE          *funct,
               DOUBLE         **deriv,
               DOUBLE         **deriv2,
               DOUBLE         **xjm,
               DOUBLE         **derxy,
               DOUBLE         **derxy2,
               DOUBLE         **evelng,
               DOUBLE         **vderxy,
               DOUBLE         **wa1,
               DOUBLE         **wa2
               );
void f3_calinta(
                  ELEMENT         *ele,
                  DOUBLE         **estif,
                  DOUBLE         **emass,
                  DOUBLE          *eforce,
                  DOUBLE         **xyze,
                  DOUBLE          *funct,
                  DOUBLE         **deriv,
                  DOUBLE         **deriv2,
                  DOUBLE         **xjm,
                  DOUBLE         **derxy,
                  DOUBLE         **derxy2,
                  DOUBLE         **evelng,
                  DOUBLE         **ealecovng,
                  DOUBLE         **egridv,
                  DOUBLE         **vderxy,
                  DOUBLE         **wa1,
                  DOUBLE         **wa2
	      );

/************************************************************************
 | f3_caliterrhs.c                                                      |
 ************************************************************************/
void f3_calgalifv(
                  DOUBLE          *eforce,
                  DOUBLE          *covint,
                  DOUBLE          *funct,
                  DOUBLE           fac,
                  INT              iel
                 );
void f3_calstabifv(
                     STAB_PAR_GLS    *gls,
                     ELEMENT         *ele,
                     DOUBLE          *eforce,
                     DOUBLE          *covint,
                     DOUBLE          *velint,
                     DOUBLE         **derxy,
                     DOUBLE         **derxy2,
                     DOUBLE           fac,
                     DOUBLE           visc,
                     INT              ihoel,
                     INT              iel
                 );
void f3_calstabifp(
                     STAB_PAR_GLS    *gls,
                     DOUBLE          *eforce,
                     DOUBLE          *covint,
                     DOUBLE          **derxy,
                     DOUBLE           fac,
                     INT              iel
                 );

/************************************************************************
 | f3_calmatvec_usfem.c                                                 |
 ************************************************************************/
void f3_calmat( DOUBLE **estif,  
		DOUBLE  *eforce,  
		DOUBLE  *velint,
		DOUBLE   histvec[3], 
		DOUBLE   gridvint[3],
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
                INT      isale
              );
void f3_calresvec(  DOUBLE  *eforce,  
                    DOUBLE  *velint,
                    DOUBLE   histvec[3], 
                    DOUBLE **vderxy, 
                    DOUBLE **vderxy2,
                    DOUBLE  *funct,  
                    DOUBLE **derxy,   
		    DOUBLE **derxy2,
                    DOUBLE  *edeadng,
                    DOUBLE   aleconv[3],
                    DOUBLE  *press,
                    DOUBLE   gradp[3],
                    DOUBLE   fac,    
                    DOUBLE   visc,   
                    INT      iel,
                    INT     *hasext,
                    INT      is_ale
              );

/************************************************************************
 | f3_calservice.c                                                      |
 ************************************************************************/
void f3_calset(
               ELEMENT         *ele,
               DOUBLE         **xyze,
               DOUBLE         **ehist,
               DOUBLE         **evelng,
               DOUBLE          *epren,
               DOUBLE          *edeadng,
               ARRAY_POSITION  *ipos,
               INT             *hasext
	      );
void f3_calseta(
                  ELEMENT         *ele,
                  DOUBLE         **xyze,
                  DOUBLE         **ehist,
                  DOUBLE         **evelng,
                  DOUBLE         **ealecovng,
                  DOUBLE         **egridv,
                  DOUBLE          *epren,
                  DOUBLE          *edeadng,
                  ARRAY_POSITION  *ipos,
                  INT             *hasext
               );
void f3_alecoor(
                  ELEMENT         *ele,
                  DOUBLE         **xyze
	       );
void f3_veci(
             DOUBLE  *vecint,
             DOUBLE  *funct,
             DOUBLE **evec,
             INT      iel
            ) ;
void f3_edgeveli(
                  DOUBLE  *velint,
                  DOUBLE  *funct,
                  DOUBLE **evel,
                  INT      ngnode,
                  INT     *iedgnod
               );
DOUBLE f3_scali(
               DOUBLE  *funct,
               DOUBLE  *epre,
               INT      iel
            );
DOUBLE f3_phii(
                  DOUBLE  *funct,
                  DOUBLE  *ephi,
                  INT     *iedgnod,
                  INT      ngnode
               );
void f3_vder(
               DOUBLE **vderxy,
               DOUBLE **derxy,
               DOUBLE **evel,
               INT      iel
            ) ;
void f3_vder2(
               DOUBLE **vderxy2,
               DOUBLE **derxy2,
               DOUBLE **evel,
               INT      iel
               );
void f3_pder(
               DOUBLE  *pderxy,
               DOUBLE **derxy,
               DOUBLE  *epre,
               INT      iel
            );
void f3_phider(
                  DOUBLE  *phiderxy,
                  DOUBLE **derxy,
                  DOUBLE  *ephi,
                  INT      ngnode,
                  INT     *iedgnod
	        );
void f3_covi(
               DOUBLE **vderxy,
               DOUBLE  *velint,
               DOUBLE  *covint
            );
void f3_permeforce(
                  DOUBLE    *eforce,
                  DOUBLE   **tmp,
                  INT        iel
                  );
void f3_permeforce_ifs(
                        DOUBLE   *eforce,
                        DOUBLE  **tmp,
                        ELEMENT  *ele
	              ) ;
void f3_permestif(
		   DOUBLE         **estif,
		   DOUBLE         **emass,
		   DOUBLE         **tmp,
		   INT              iel
	          );
void f3_permestif_ifs(
                        DOUBLE         **estif,
                        DOUBLE         **emass,
                        DOUBLE         **tmp,
                        ELEMENT         *ele
                     );
void f3_iedg(
               INT     *iegnod,
               ELEMENT *ele,
               INT      surf
	     );

/************************************************************************
 | f3_calstabmat.c                                                      |
 ************************************************************************/
void f3_calstabkvv(
                     ELEMENT         *ele,
                     STAB_PAR_GLS    *gls,
                     DOUBLE         **estif,
                     DOUBLE          *velint,
                     DOUBLE          *vel2int,
                     DOUBLE          *gridvint,
                     DOUBLE         **vderxy,
                     DOUBLE          *funct,
                     DOUBLE         **derxy,
                     DOUBLE         **derxy2,
                     DOUBLE           fac,
                     DOUBLE           visc,
                     INT              iel,
                     INT              ihoel
                  );
void f3_calstabkvp(
                     ELEMENT         *ele,
                     STAB_PAR_GLS    *gls,
                     DOUBLE         **estif,
                     DOUBLE          *velint,
                     DOUBLE          *funct,
                     DOUBLE         **derxy,
                     DOUBLE         **derxy2,
                     DOUBLE           fac,
                     DOUBLE           visc,
                     INT              iel,
                     INT              ihoel
                   );
void f3_calstabkvg(
                     ELEMENT         *ele,
                     STAB_PAR_GLS    *gls,
                     DOUBLE         **estif,
                     DOUBLE         **vderxy,
                     DOUBLE          *funct,
                     DOUBLE         **derxy,
                     DOUBLE         **derxy2,
                     DOUBLE          *alecovint,
                     DOUBLE           fac,
                     DOUBLE           visc,
                     INT              iel,
                     INT              ihoel
                   );
void f3_calstabmvv(
                     ELEMENT         *ele,
                     STAB_PAR_GLS    *gls,
                     DOUBLE         **estif,
                     DOUBLE          *velint,
                     DOUBLE          *funct,
                     DOUBLE         **derxy,
                     DOUBLE         **derxy2,
                     DOUBLE           fac,
                     DOUBLE           visc,
                     INT              iel,
                     INT              ihoel
                   );
void f3_calstabkpv(
                     ELEMENT         *ele,
                     STAB_PAR_GLS    *gls,
                     DOUBLE         **estif,
                     DOUBLE          *velint,
                     DOUBLE          *gridvint,
                     DOUBLE         **vderxy,
                     DOUBLE          *funct,
                     DOUBLE         **derxy,
                     DOUBLE         **derxy2,
                     DOUBLE           fac,
                     DOUBLE           visc,
                     INT              iel,
                     INT              ihoel
                   );
void f3_calstabkpg(
                    STAB_PAR_GLS    *gls,
		    DOUBLE         **estif,
		    DOUBLE          *funct,
		    DOUBLE         **vderxy,
		    DOUBLE         **derxy,
		    DOUBLE           fac,
		    INT              iel
                   );

void f3_calstabkpp(
                     STAB_PAR_GLS    *gls,
                     DOUBLE         **estif,
                     DOUBLE         **derxy,
                     DOUBLE           fac,
                     INT              iel
                   );
void f3_calstabmpv(
                     STAB_PAR_GLS    *gls,
                     DOUBLE         **estif,
                     DOUBLE          *funct,
                     DOUBLE         **derxy,
                     DOUBLE           fac,
                     INT              iel
                   );

/************************************************************************
 | f3_calstabpar.c                                                      |
 ************************************************************************/
void f3_calstabpar(
	            ELEMENT         *ele,
		    DOUBLE          *velint,
		    DOUBLE           visc,
		    INT              iel,
		    DIS_TYP          typ,
		    INT              iflag
                  );

/************************************************************************
 | f3_caltau.c                                                          |
 ************************************************************************/
void f3_caltau(			     
	       ELEMENT         *ele, 
	       DOUBLE         **xyze,
	       DOUBLE          *funct,  
	       DOUBLE         **deriv,  
	       DOUBLE         **derxy,              
	       DOUBLE         **xjm,  
	       DOUBLE         **evelng,
               DOUBLE         **wa1,
               DOUBLE           visc  
              );

/************************************************************************
 | f3_inpele.c                                                          |
 ************************************************************************/
void f3inp(ELEMENT *ele, INT counter);

/************************************************************************
 | f3_intg.c                                                            |
 ************************************************************************/
void f3_intg(
             INT              option);

/************************************************************************
 | f3_int_usfem.c                                                       |
 ************************************************************************/
void f3_int_usfem(
	              ELEMENT         *ele,
                      INT             *hasext,
                      DOUBLE         **estif,
	              DOUBLE          *force,
	              DOUBLE         **xyze,
	              DOUBLE          *funct,
	              DOUBLE         **deriv,
	              DOUBLE         **deriv2,
	              DOUBLE         **xjm,
	              DOUBLE         **derxy,
	              DOUBLE         **derxy2,
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
void f3_int_res(
	        ELEMENT         *ele,
                INT             *hasext,
	        DOUBLE          *force,
	        DOUBLE         **xyze,
	        DOUBLE          *funct,
	        DOUBLE         **deriv,
	        DOUBLE         **deriv2,
	        DOUBLE         **xjm,
	        DOUBLE         **derxy,
	        DOUBLE         **derxy2,
	        DOUBLE         **evelng,
	        DOUBLE         **evhist,
	        DOUBLE         **ealecovng,
	        DOUBLE          *epren,
	        DOUBLE          *edeadng,
	        DOUBLE         **vderxy,
                DOUBLE         **vderxy2,
                DOUBLE           visc,
	        DOUBLE         **wa1,
	        DOUBLE         **wa2
	       );

/************************************************************************
 | f3_main.c                                                            |
 ************************************************************************/
void fluid3(
            PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
            ARRAY       *estif_global,
            ARRAY       *emass_global,
	    ARRAY       *eforce_global,
	    ARRAY       *edforce_global,
            CALC_ACTION *action,
	    INT         *hasdirich,
	    INT         *hasext,
	    CONTAINER   *container
           );

/************************************************************************
 | f3_massrhs.c                                                         |
 ************************************************************************/
void f3_massrhs(ELEMENT *ele,
                DOUBLE **emass,
                DOUBLE **hist,
                DOUBLE  *edeadng,
                DOUBLE *eiforce,
                INT     *hasext);

/*------------------------- the following originates form Neumann... ---*/
void f3_calelestress(
    INT             viscstr,
    ELEMENT        *ele,
    DOUBLE        **evel,
    DOUBLE         *epre,
    DOUBLE         *funct,
    DOUBLE        **deriv,
    DOUBLE        **derxy,
    DOUBLE        **vderxy,
    DOUBLE        **xjm,
    DOUBLE        **wa1,
    DOUBLE        **xyze,
    DOUBLE        **sigmaint,
    ARRAY_POSITION *ipos
    );
DOUBLE f3_rsn (
    INT   node,
    INT   irs
    );
void f3_lgpl (
    INT         i,
    INT         n,
    DOUBLE    *zr,
    DOUBLE      z,
    DOUBLE *value
    );
void f3_hxsm (
    INT nir,
    INT nis,
    INT nit,
    DOUBLE rk,
    DOUBLE sk,
    DOUBLE tk,
    DOUBLE f[6][27],
    DOUBLE *fp,
    DOUBLE *xgr,
    DOUBLE *xgs,
    DOUBLE *xgt
    );
void f3_sext(
    DOUBLE  **nostrs,
    DOUBLE  **gpstress,
    DOUBLE   *xgr,
    DOUBLE   *xgs,
    DOUBLE   *xgt,
    INT       nir,
    INT       nis,
    INT       nit,
    INT       iel
    );
void f3_hex2(
    DOUBLE      funct[MAXNOD_F3],
    DOUBLE      deriv[3][MAXNOD_F3],
    DOUBLE      r,
    DOUBLE      s,
    DOUBLE      t,
    DIS_TYP     typ
    );

/************************************************************************
 | f3_liftdrag.c                                                        |
 ************************************************************************/
void f3_liftdrag(
    ELEMENT       *ele,
    CONTAINER     *container
    );
void f3_jaco2(
    DOUBLE      deriv[3][MAXNOD_F3],
    DOUBLE      xjm[3][3],
    DOUBLE     *det,
    ELEMENT    *ele,
    INT         iel
    );


void f3_out_gid_sol_str(
    FILE       *out,
    FIELD *actfield,
    INT        type,
    INT         init
    );

/************************************************************************
 | f3_mlsubmesh.c                                                       |
 ************************************************************************/
void f3_pdsubmesh(FLUID_ML_SMESH *smesh,
                  INT             xele,
                  INT             yele,
                  INT             zele,
                  INT             order,
		  INT             flag);
void f3_elesubmesh(ELEMENT        *ele,
                   FLUID_ML_SMESH *smesh,
		   INT             flag);

/************************************************************************
 | f3_restart.c                                                         |
 ************************************************************************/
void f3_write_restart(ELEMENT *actele, INT nhandle, long int *handles);
void f3_read_restart(ELEMENT *actele, INT nhandle, long int *handles);
