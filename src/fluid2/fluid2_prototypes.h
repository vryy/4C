/*!----------------------------------------------------------------------
\file
\brief fluid2 prototypes

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
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
/************************************************************************
 | f2_ass_dof_ndis_tu.c                                                 |
 ************************************************************************/
void f2tu_ass_dof(NODE *actnode, INT *counter);

/************************************************************************
 | f2_calcurvature.c                                                    |
 ************************************************************************/
void f2_calq4curv( GLINE                        **actgline,
                     FLUID_FREESURF_CONDITION   **actlinefs,
		     ELEMENT                    *actele,
		     INT                         foundline,
		     INT                         actngline,
		     DOUBLE                    **xyze,
		     DOUBLE                    **deriv,
		     DOUBLE                    **kappa
		   );
void f2_calq8curv( GLINE                        **actgline,
                     FLUID_FREESURF_CONDITION   **actlinefs,
		     ELEMENT                    *actele,
		     INT                         actngline,
		     DOUBLE                    **xyze,
		     DOUBLE                    **kappa
		   );

 /************************************************************************
 | f2_calele.c                                                          |
 ************************************************************************/
void f2_calele(
	        ELEMENT        *ele,
                ELEMENT        *eleke,
                ARRAY          *estif_global,
                ARRAY          *emass_global,
                ARRAY          *eforce_global,
                ARRAY          *edforce_global,
                INT            *hasdirich,
                INT            *hasext,
                INT             imyrank,
                ARRAY_POSITION *ipos,
                INT             is_relax,
                INT             init
	       );
void f2_stress(FLUID_STRESS str, INT viscstr ,
               ELEMENT *ele,
               ARRAY_POSITION *ipos,
               INT is_relax );
void f2_curvature(
                     ELEMENT        *ele,
                     INT             imyrank
		                );
void f2_heightfunc(
                   ELEMENT              *ele,
                   ARRAY                *estif_global,
		   ARRAY                *eforce_global,
                   ARRAY_POSITION       *ipos,
		   CONTAINER            *container
		   );
void f2_calstab(ELEMENT *ele, ARRAY_POSITION *ipos);
void f2_calnormal(ELEMENT *ele);
void f2_caleleres(
	           ELEMENT        *ele,
	           ARRAY          *eforce_global,
                   ARRAY_POSITION *ipos,
                   INT            *hasdirich,
                   INT            *hasext
	       );
void f2_caleleres_relax(ELEMENT        *ele,
			ARRAY          *estif_global,
			ARRAY          *eforce_global,
			ARRAY_POSITION *ipos,
			INT            *hasdirich,
			INT            *hasext);
void f2_calstresspro(ELEMENT   *ele,
                     INT       *hasext,
                     ARRAY     *estif_global,
		     ARRAY     *eforce_global,
                     ARRAY_POSITION *ipos
  );

void f2_calerr(
    ELEMENT          *ele,
    CONTAINER        *container,
    ARRAY_POSITION   *ipos
    );

/************************************************************************
 | f2_calele_tu.c                                                       |
 ************************************************************************/
void f2_calele_tu(
	          ELEMENT        *eleke,
                  ELEMENT        *elev,
                  ARRAY          *estif_global,
                  ARRAY          *emass_global,
	          ARRAY          *etforce_global,
	          ARRAY          *eiforce_global,
                  ARRAY          *edforce_global,
                  ARRAY          *eproforce_global,
                  ARRAY_POSITION *fluid_ipos,
                  INT            *hasdirich,
                  INT            *hasext,
                  INT             init
                );

/************************************************************************
 | f2_calele_tu_1.c                                                     |
 ************************************************************************/
void f2_calele_tu_1(
                     ELEMENT        *eleke,
                     ELEMENT        *elev,
                     ARRAY          *estif_global,
                     ARRAY          *emass_global,
                     ARRAY          *etforce_global,
                     ARRAY          *eiforce_global,
                     ARRAY          *edforce_global,
                     ARRAY          *eproforce_global,
                     ARRAY_POSITION *ipos,
                     INT            *hasdirich,
                     INT            *hasext,
                     INT             init
	                 );
/************************************************************************
 | f2_calelesize.c                                                      |
 ************************************************************************/
void f2_calelesize(
                     ELEMENT         *ele,
                     ELEMENT         *eleke,
                     DOUBLE         **xyze,
                     DOUBLE          *funct,
                     DOUBLE         **deriv,
                     DOUBLE         **deriv2,
                     DOUBLE         **xjm,
                     DOUBLE         **derxy,
                     DOUBLE         **vderxy,
                     DOUBLE         **evel,
                     DOUBLE          *ephin,
                     DOUBLE          *ephing,
                     DOUBLE          *eddy,
                     DOUBLE          *visc,
                     ARRAY_POSITION  *ipos,
                     INT              cpele
                     );
void f2_calelesize2(
                     ELEMENT         *ele,
                     DOUBLE         **xyze,
                     DOUBLE          *funct,
                     DOUBLE          *velint,
                     DOUBLE           visc,
                     INT              iel,
                     DIS_TYP          typ
                    );
void f2_calstrlen(
                  DOUBLE   *strle,
                  DOUBLE  **xyze,
                  DOUBLE   *velint,
                  ELEMENT  *ele,
                  DOUBLE   *gcoor,
                  DIS_TYP   typ
                 );

/************************************************************************
 | f2_calelesize_tu.c                                                   |
 ************************************************************************/
void f2_calelesize_tu(
                        ELEMENT         *ele,
                        ELEMENT         *elev,
                        DOUBLE          *funct,
                        DOUBLE         **deriv,
                        DOUBLE         **deriv2,
                        DOUBLE         **evel,
                        DOUBLE          *eddyg,
                        DOUBLE          *velint,
                        DOUBLE          *velint_dc,
                        DOUBLE          *kapepsn,
                        DOUBLE         **xjm,
                        DOUBLE         **xyze,
                        DOUBLE         **derxy,
                        DOUBLE          *kapepsderxy,
                        DOUBLE         **cutp
                     );

void f2_calstrlen_tu(
		         DOUBLE   *velint_dc,
		         ELEMENT  *ele,
                     DOUBLE   *gcoor,
		         DOUBLE  **cutp,
		         DIS_TYP       typ
                     );

/************************************************************************
 | f2_calelesize_tu_1.c                                                 |
 ************************************************************************/
void f2_calelesize_tu_1(
                        ELEMENT         *ele,
                        ELEMENT         *elev,
                        DOUBLE          *funct,
                        DOUBLE         **deriv,
                        DOUBLE         **deriv2,
                        DOUBLE         **evel,
                        DOUBLE          *eddyg,
                        DOUBLE          *velint,
                        DOUBLE          *velint_dc,
                        DOUBLE          *kapomen,
                        DOUBLE         **xjm,
                        DOUBLE         **xyze,
                        DOUBLE         **derxy,
                        DOUBLE          *kapomederxy,
                        DOUBLE         **cutp
                 );

/************************************************************************
 | f2_calextrhs.c                                                       |
 ************************************************************************/
void f2_calgalexfv(
                     DOUBLE          *eforce,
                     DOUBLE          *funct,
                     DOUBLE          *edeadn,
                     DOUBLE          *edeadng,
                     DOUBLE           fac,
                     INT              iel
              );
void f2_calstabexfv(
                    STAB_PAR_GLS    *gls,
                    ELEMENT         *ele,
                    DOUBLE          *eforce,
                    DOUBLE         **derxy,
                    DOUBLE         **derxy2,
                    DOUBLE          *edead,
                    DOUBLE          *velint,
                    DOUBLE           fac,
                    DOUBLE           visc,
                    INT              iel,
                    INT              ihoel,
                    INT              flag
                   );
void f2_calstabexfp(
                     STAB_PAR_GLS    *gls,
                     DOUBLE          *eforce,
                     DOUBLE         **derxy,
                     DOUBLE          *edead,
                     DOUBLE           fac,
                     INT              iel,
                     INT              flag
                   );
void f2_calneumann(ELEMENT         *ele,
                   INT              imyrank,
                   DOUBLE          *eforce,
                   DOUBLE         **xyze,
                   DOUBLE          *funct,
                   DOUBLE         **deriv,
                   DOUBLE         **xjm,
                   DOUBLE          *edeadn,
                   DOUBLE          *edeadng
  );
/************************************************************************
 | f2_calfuncderiv.c                                                    |
 ************************************************************************/
void f2_rec(
            DOUBLE     *funct,
            DOUBLE    **deriv,
            DOUBLE    **deriv2,
            DOUBLE      r,
            DOUBLE      s,
            DIS_TYP     typ,
            INT         icode
	   );
void f2_recex(
               DOUBLE           *funval,
               DOUBLE           *fpar,
               DOUBLE            r,
               DOUBLE            s,
               DOUBLE           *fval,
               INT               igauss,
               INT               icode
            );
void f2_tri(
            DOUBLE     *funct,
            DOUBLE    **deriv,
            DOUBLE    **deriv2,
            DOUBLE      r,
            DOUBLE      s,
            DIS_TYP     typ,
            INT         icode
	   );
void f2_triex(
               DOUBLE           *funval,
               DOUBLE           *fpar,
               DOUBLE            r,
               DOUBLE            s,
               DOUBLE           *fval,
               INT               igauss,
               INT               icode
            );
void f2_degrectri(
                    DOUBLE     *funct,
                    DOUBLE    **deriv,
                    DOUBLE      r,
                    DIS_TYP     typ,
                    INT         option
		 );
void f2_jaco(
             DOUBLE    **xyze,
             DOUBLE    **deriv,
             DOUBLE    **xjm,
             DOUBLE     *det,
             INT         iel,
             ELEMENT    *ele
	    );
void f2_jaco2(
             DOUBLE    **xyze,
             DOUBLE    **deriv,
             DOUBLE    **xjm,
             DOUBLE     *det,
             INT         iel
	    );
void f2_edgejaco(
                 DOUBLE    **xyze,
                 DOUBLE    **deriv,
                 DOUBLE    **xjm,
                 DOUBLE     *det,
                 INT         iel,
                 INT        *iedgnod
	        );
void f2_gder(
               DOUBLE   **derxy,
               DOUBLE   **deriv,
               DOUBLE   **xjm,
               DOUBLE     det,
               INT        iel
	    );
void f2_edgegder(
                  DOUBLE   **deriv,
                  DOUBLE   **derxy,
                  DOUBLE   **xjm,
                  INT        ngnode
	        );
void f2_gcoor(
               DOUBLE    **xyze,
               DOUBLE     *funct,
               INT         iel,
               DOUBLE     *gcoor
             );
void f2_gder2(
               DOUBLE     **xyze,
               DOUBLE     **xjm,
               DOUBLE     **bi,
               DOUBLE     **xder2,
               DOUBLE     **derxy,
               DOUBLE     **derxy2,
               DOUBLE     **deriv2,
               INT          iel
	     );

/************************************************************************
 | f2_calgalmat.c                                                       |
 ************************************************************************/
void f2_calkvv(
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
void f2_calkvp(
               DOUBLE         **estif,
               DOUBLE          *funct,
               DOUBLE         **derxy,
               DOUBLE           fac,
               INT              iel
              );
void f2_calkvg(
               DOUBLE         **estif,
               DOUBLE         **vderxy,
               DOUBLE          *funct,
               DOUBLE           fac,
               INT              iel
              );
void f2_calmvv(
               DOUBLE         **emass,
               DOUBLE          *funct,
               DOUBLE           fac,
               INT              iel
              );
void f2_calgfskgg(
                  DOUBLE         **estif,
                  DOUBLE          *funct,
                  DOUBLE         **deriv,
                  DOUBLE         **xjm,
                  DOUBLE           det,
                  DOUBLE           fac,
                  INT              iel

                );
void f2_calkgedge(
                  DOUBLE         **estif,
                  DOUBLE          *funct,
                  DOUBLE           fac,
                  INT             *iedgnod,
                  INT              iel,
                  INT              ngnode
                );
void f2_calgfskgedge(
                     DOUBLE         **estif,
                     DOUBLE          *funct,
                     DOUBLE          *vnint,
                     DOUBLE           fac,
                     INT             *iedgnod,
                     INT              iel,
                     INT              ngnode

                  );

/************************************************************************
 | f2_calliftdrag.c                                                     |
 ************************************************************************/
void f2_calliftdrag(
    ELEMENT       *ele,
    CONTAINER     *container);

/************************************************************************
 | f2_calheightfunc_sep.c                                               |
 ************************************************************************/
void f2_stabpar_hfsep(
                        ELEMENT          *ele,
                        DOUBLE          **xyze,
                        DOUBLE          **evelng,
                        DOUBLE           *funct,
                        DOUBLE           *velint,
                        DOUBLE            phiintn,
                        DOUBLE            phiintng,
                        DOUBLE            phiderx,
                        DOUBLE            e1,
                        INT              *iedgnod,
                        INT               ngnode,
                        DIS_TYP           typ
                     );
void f2_calint_hfsep(
                     ELEMENT           *ele,
                     DOUBLE            *funct,
                     DOUBLE           **deriv,
                     DOUBLE           **xjm,
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
void f2_calmat_vhf_sep(
		       ELEMENT                *ele,
		       DOUBLE                **estif,
		       INT                     ngnode,
		       DOUBLE                 *funct,
		       DOUBLE                **derxy,
		       DOUBLE                 *velint,
		       DOUBLE                  fac
	              );
void f2_caliterhs_vhf_sep(
 	 	           ELEMENT                *ele,
		           DOUBLE                 *eforce,
		           INT                     ngnode,
		           DOUBLE                 *funct,
                           DOUBLE                **derxy,
		           DOUBLE                 *velint,
			   DOUBLE                  fac
	                  );
void f2_caltimerhs_vhf_sep(
                           ELEMENT                *ele,
                           DOUBLE                 *eforce,
                           INT                     ngnode,
                           DOUBLE                 *funct,
                           DOUBLE                **derxy,
                           DOUBLE                 *velint, /* at n+1 */
                           DOUBLE                 *vel2int,/* at n */
                           DOUBLE                  phiint,
                           DOUBLE                  phiderx,
                           DOUBLE                  fac
                          );


/************************************************************************
 | f2_calheightfuncmat.c                                                |
 ************************************************************************/
void f2_calmat_vhf(
                        DOUBLE         **emass,
                        DOUBLE         **estif,
                        DOUBLE          *funct,
                        DOUBLE         **derxy,
                        DOUBLE          *velint,
                        DOUBLE           phiderxng,
                        DOUBLE           fac,
                        INT              iel,
                        INT              ngnode,
                        INT             *iedgnod
                         );
void f2_calstabmat_vhf(
                           ELEMENT         *ele,
		           DOUBLE         **emass,
                           DOUBLE         **estif,
		           DOUBLE          *funct,
                           DOUBLE         **derxy,
                           DOUBLE          *velint,
		           DOUBLE           fac,
		           INT              iel,
	  	           INT              ngnode,
                           INT             *iedgnod
                         );
/************************************************************************
 | f2_calheightfuncrhs.c                                                |
 ************************************************************************/
void f2_calrhs_vhf(
                     DOUBLE          *eforce,
		     DOUBLE          *velint,  /* at n+1 */
		     DOUBLE          *vel2int, /* at n */
		     DOUBLE           phiint,
	   	     DOUBLE          *funct,
                     DOUBLE         **derxy,
		     DOUBLE           phiderxng,  /* at n+1 */
		     DOUBLE           phiderxn,   /* at n */
		     DOUBLE           fac,
		     INT              iel,
		     INT              ngnode,
		     INT             *iedgnod
                    );
void f2_calstabtf_vhf(
                     ELEMENT         *ele,
                     DOUBLE          *eforce,
		     DOUBLE          *velint,
		     DOUBLE           phiint,
	   	     DOUBLE          *funct,
                     DOUBLE         **derxy,
		     DOUBLE           phiderx,
		     DOUBLE           fac,
		     INT              iel,
		     INT              ngnode,
		     INT             *iedgnod
                    ) ;
void f2_calgalif_vhf(
                     DOUBLE          *eforce,
	   	     DOUBLE          *funct,
		     DOUBLE          *velint,
		     DOUBLE           phiderx,
		     DOUBLE           fac,
		     INT              iel,
		     INT              ngnode,
		     INT             *iedgnod
                    );

/************************************************************************
 | f2_calgalmat_tu.c                                                    |
 ************************************************************************/
void f2_calkkapeps(
                     DOUBLE         **estif,
                     DOUBLE           kapepsint,
                     DOUBLE          *velint,
                     DOUBLE           eddyint,
                     DOUBLE          *kapepsderxy,
                     DOUBLE          *funct,
                     DOUBLE         **derxy,
                     DOUBLE           fac,
                     DOUBLE           visc,
                     DOUBLE           factor,
                     DOUBLE           sig,
                     INT              iel
              );
void f2_calmkapeps(
                     DOUBLE         **emass,
                     DOUBLE          *funct,
                     DOUBLE           fac,
                     INT              iel
              );

/************************************************************************
 | f2_calgalmat_tu_1.c                                                  |
 ************************************************************************/
void f2_calkkapome(
                     DOUBLE         **estif,
                     DOUBLE           kapomeint,
                     DOUBLE          *velint,
                     DOUBLE           eddyint,
                     DOUBLE          *funct,
                     DOUBLE         **derxy,
                     DOUBLE           fac,
                     DOUBLE           visc,
                     DOUBLE           factor,
                     DOUBLE           sig,
                     INT              iel
              );

void f2_calmkapome(
                     DOUBLE         **emass,
                     DOUBLE          *funct,
                     DOUBLE           fac,
                     INT              iel
              );

/************************************************************************
 | f2_calint.c                                                          |
 ************************************************************************/
void f2_calint(
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
	       DOUBLE         **wa2,
               DOUBLE           visc
	      );
void f2_calinta(
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
                  DOUBLE         **eveln,
                  DOUBLE         **evelng,
                  DOUBLE         **ealecovng,
                  DOUBLE         **egridv,
                  DOUBLE         **vderxy,
                  DOUBLE          *ekappan,
                  DOUBLE          *ekappang,
                  DOUBLE          *ephin,
                  DOUBLE          *ephing,
                  DOUBLE         **evnng,
                  DOUBLE         **wa1,
                  DOUBLE         **wa2
               );

/************************************************************************
 | f2_calint_tu.c                                                       |
 ************************************************************************/
void f2_calint_tu(
                  ELEMENT         *eleke,
                  ELEMENT         *elev,
                  DOUBLE         **estif,
                  DOUBLE         **emass,
                  DOUBLE          *etforce,
                  DOUBLE          *eiforce,
                  DOUBLE          *eproforce,
                  DOUBLE          *funct,
                  DOUBLE         **deriv,
                  DOUBLE         **deriv2,
                  DOUBLE         **xjm,
                  DOUBLE         **xyze,
                  DOUBLE         **derxy,
                  DOUBLE         **derxy2,
                  DOUBLE          *kapepsn,
                  DOUBLE          *kapepsg,
                  DOUBLE          *eddyg,
                  DOUBLE          *eddypro,
                  DOUBLE          *kappa,
                  DOUBLE          *kappan,
                  DOUBLE          *epsilon,
                  DOUBLE          *kapepspro,
                  DOUBLE          *kapepsderxy,
                  DOUBLE          *kapepsderxy2,
                  DOUBLE          *velint,
                  DOUBLE          *velint_dc,
                  DOUBLE         **evel,
                  DOUBLE         **vderxy,
                  DOUBLE         **vderxy2,
                  DOUBLE         **wa1,
                  DOUBLE         **wa2
               );

/************************************************************************
 | f2_calint_tu_1.c                                                     |
 ************************************************************************/
void f2_calint_tu_1(
                     ELEMENT         *ele,
                     ELEMENT         *elev,
                     DOUBLE         **estif,
                     DOUBLE         **emass,
                     DOUBLE          *etforce,
                     DOUBLE          *eiforce,
                     DOUBLE          *eproforce,
                     DOUBLE          *funct,
                     DOUBLE         **deriv,
                     DOUBLE         **deriv2,
                     DOUBLE         **xjm,
                     DOUBLE         **xyze,
                     DOUBLE         **derxy,
                     DOUBLE         **derxy2,
                     DOUBLE          *kapomen,
                     DOUBLE          *kapomeg,
                     DOUBLE          *eddyg,
                     DOUBLE          *eddypro,
                     DOUBLE          *kappan,
                     DOUBLE          *omega,
                     DOUBLE          *kapomepro,
                     DOUBLE          *omegaderxy,
                     DOUBLE          *kapomederxy,
                     DOUBLE          *kapomederxy2,
                     DOUBLE          *velint,
                     DOUBLE          *velint_dc,
                     DOUBLE         **evel,
                     DOUBLE         **vderxy,
                     DOUBLE         **vderxy2,
                     DOUBLE         **wa1,
                     DOUBLE         **wa2
	        );
/************************************************************************
 | f2_caliterrhs.c                                                      |
 ************************************************************************/
void f2_calgalifv(
                  DOUBLE          *eforce,
                  DOUBLE          *covint,
                  DOUBLE          *funct,
                  DOUBLE           fac,
                  INT              iel
                 ) ;
void f2_calstabifv(
                     STAB_PAR_GLS    *gls,
                     ELEMENT         *ele,
                     DOUBLE          *eforce,
                     DOUBLE          *covint,
                     DOUBLE          *velint,
                     DOUBLE          *funct,
                     DOUBLE         **derxy,
                     DOUBLE         **derxy2,
                     DOUBLE           fac,
                     DOUBLE           visc,
                     INT              ihoel,
                     INT              iel
                 );
void f2_calstabifp(
                     STAB_PAR_GLS    *gls,
                     DOUBLE          *eforce,
                     DOUBLE          *covint,
                     DOUBLE         **derxy,
                     DOUBLE           fac,
                     INT              iel
                 );

/************************************************************************
 | f2_caliterrhs_tu.c                                                   |
 ************************************************************************/
void f2_calgalifkapeps(
                        DOUBLE          *eforce,
                        DOUBLE           eddyint,
                        DOUBLE           kapepsint,
                        DOUBLE          *funct,
                        DOUBLE           fac,
                        DOUBLE           factor2,
                        DOUBLE           vderxy_12,
                        DOUBLE           visc,
                        INT              iel
                 );
void f2_calstabifkapeps(
                        ELEMENT         *ele,
                        DOUBLE          *eforce,
                        DOUBLE           kapepsint,
                        DOUBLE          *velint,
                        DOUBLE          *velint_dc,
                        DOUBLE           eddyint,
                        DOUBLE          *funct,
                        DOUBLE         **derxy,
                        DOUBLE           fac,
                        DOUBLE           factor2,
                        DOUBLE           vderxy_12,
                        DOUBLE           visc,
                        INT              iel
                  );

/************************************************************************
 | f2_caliterrhs_tu_1.c                                                 |
 ************************************************************************/
void f2_calgalifkapome(
                        DOUBLE          *eforce,
                        DOUBLE           kapomeint,
                        DOUBLE          *funct,
                        DOUBLE           fac,
                        DOUBLE           factor2,
                        INT              iel
                 );

void f2_calstabifkapome(
                        ELEMENT         *ele,
                        DOUBLE          *eforce,
                        DOUBLE           kapomeint,
                        DOUBLE          *velint,
                        DOUBLE          *velint_dc,
                        DOUBLE          *funct,
                        DOUBLE         **derxy,
                        DOUBLE           fac,
                        DOUBLE           factor2,
                        INT              iel
                  );

/************************************************************************
 | f2_calmatvec_usfem.c                                                 |
 ************************************************************************/
void f2_calmat( DOUBLE **estif,
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
                INT      is_relax
              );
void f2_calresvec(  DOUBLE  *eforce,
                    DOUBLE  *velint,
                    DOUBLE   histvec[2],
                    DOUBLE **vderxy,
                    DOUBLE **vderxy2,
                    DOUBLE  *funct,
                    DOUBLE **derxy,
		    DOUBLE **derxy2,
                    DOUBLE  *edeadng,
                    DOUBLE   aleconv[2],
                    DOUBLE  *press,
                    DOUBLE   gradp[2],
                    DOUBLE   fac,
                    DOUBLE   visc,
                    INT      iel,
                    INT     *hasext,
                    INT      is_ale
              );

/************************************************************************
 | f2_calservice.c                                                      |
 ************************************************************************/
void f2_calset(
	        ELEMENT         *ele,
                DOUBLE         **xyze,
                DOUBLE         **eveln,
	        DOUBLE         **evelng,
	        DOUBLE         **evhist,
	        DOUBLE          *epren,
		DOUBLE          *edeadn,
		DOUBLE          *edeadng,
                ARRAY_POSITION *ipos,
		INT             *hasext
              );
void f2_calseta(
                  ELEMENT         *ele,
                  DOUBLE         **xyze,
                  DOUBLE         **eveln,
                  DOUBLE         **evelng,
                  DOUBLE         **evhist,
                  DOUBLE         **ealecovn,
                  DOUBLE         **ealecovng,
                  DOUBLE         **egridv,
                  DOUBLE          *epren,
                  DOUBLE          *edeadn,
                  DOUBLE          *edeadng,
                  DOUBLE          *ekappan,
                  DOUBLE          *ekappang,
                  DOUBLE          *ephin,
                  DOUBLE          *ephing,
                  DOUBLE         **evnng,
                  DOUBLE         **evnn,
                  ARRAY_POSITION *ipos,
                  INT             *hasext,
                  INT              is_relax
	      );
void f2_alecoor(
                  ELEMENT         *ele,
                  DOUBLE         **xyze
	       );
void f2_alecoor_sd(
                     ELEMENT         *ele,
                     DOUBLE         **xyze
	          );
void f2_veci(
               DOUBLE  *vecint,
               DOUBLE  *funct,
               DOUBLE **evec,
               INT      iel
	    );
void f2_edgeveci(
                  DOUBLE  *vecint,
                  DOUBLE  *funct,
                  DOUBLE **evec,
                  INT      ngnode,
                  INT     *iedgnod
               );
DOUBLE f2_scali(
               DOUBLE  *funct,
               DOUBLE  *escal,
               INT      iel
	     ) ;
DOUBLE f2_edgescali(
                    DOUBLE  *funct,
                    DOUBLE  *escal,
                    INT     *iedgnod,
                    INT      ngnode
	           );
void f2_vder(
               DOUBLE **vderxy,
               DOUBLE **derxy,
               DOUBLE **evel,
               INT      iel
	    );
void f2_vder2(
               DOUBLE **vderxy2,
               DOUBLE **derxy2,
               DOUBLE **evel,
               INT      iel
             );
void f2_pder(
               DOUBLE  *pderxy,
               DOUBLE **derxy,
               DOUBLE  *epre,
               INT      iel
	    );
DOUBLE f2_phider(
                  DOUBLE **derxy,
                  DOUBLE  *ephi,
                  INT      ngnode,
                  INT     *iedgnod
	        );
void f2_covi(
               DOUBLE **vderxy,
               DOUBLE  *velint,
               DOUBLE  *covint
	    );
void f2_permeforce(
                     DOUBLE   *eforce,
                     DOUBLE  **tmp,
                     INT       iel
	          );
void f2_permeforce_ifs(
                        DOUBLE   *eforce,
                        DOUBLE  **tmp,
                        ELEMENT  *ele
	              );
void f2_permestif(
                  DOUBLE         **estif,
                  DOUBLE         **emass,
                  DOUBLE         **tmp,
                  ELEMENT         *ele
                 );
void f2_permestif_ifs(
                        DOUBLE         **estif,
                        DOUBLE         **emass,
                        DOUBLE         **tmp,
                        ELEMENT         *ele
                     );
void f2_iedg(
               INT     *iegnod,
               ELEMENT *ele,
               INT      line,
               INT      init
	     );

/************************************************************************
| f2_calservice_tu.c                                                    |
************************************************************************/
void f2_calset_tu(
                     ELEMENT         *ele,
                     ELEMENT         *elev,
                     DOUBLE          *kapepsn,
                     DOUBLE          *kapepsg,
                     DOUBLE          *kapepspro,
                     DOUBLE          *eddyg,
                     DOUBLE          *eddypro,
                     DOUBLE          *kappa,
                     DOUBLE          *kappan,
                     DOUBLE          *epsilon,
                     DOUBLE         **evel,
                     DOUBLE         **xyze,
                     ARRAY_POSITION *fluid_ipos
	           );

void f2_shearstress(
                     ELEMENT    *ele,
                     ARRAY_POSITION *ipos
                 );

void f2_kapepsi(
                  DOUBLE  *kapepsint,
                  DOUBLE  *funct,
                  DOUBLE  *kapeps,
                  INT      iel
	     );
void f2_eddyi(
               DOUBLE  *eddyint,
               DOUBLE  *funct,
                DOUBLE  *eddy,
               INT      iel
	     );
void f2_kappai_tu(
                  DOUBLE     *kappaint,
                  DOUBLE     *kappanint,
                  DOUBLE     *eps_proint,
                  DOUBLE     *funct,
                  DOUBLE     *kappa,
                  DOUBLE     *kappan,
                  DOUBLE     *eps_pro,
                  INT         iel
	     );

void f2_C_kappa(
                  DOUBLE      kapepsint,
                  DOUBLE     *epsilon,
                  DOUBLE     *funct,
                  DOUBLE      visc,
                  DOUBLE     *C_u,
                  INT         iel
	     );
void f2_C_eps(
                  DOUBLE      kapepsint,
                  DOUBLE      kappaint,
                  DOUBLE      visc,
                  DOUBLE     *C_2,
                  INT         iel
	     );

void f2_v(
            DOUBLE    **vderxy2,
            DOUBLE     *vderxy_12
	   );

void f2_fac_kappa(
                  DOUBLE   C_u,
                  DOUBLE   eddyint,
                  DOUBLE  *factor,
                  DOUBLE  *factor1,
                  DOUBLE  *factor2,
                  DOUBLE  *sig
	           );

void f2_fac_eps(
                  DOUBLE   C_2,
                  DOUBLE   eps_proint,
                  DOUBLE   kappaint,
                  DOUBLE   kappanint,
                  DOUBLE  *factor,
                  DOUBLE  *factor1,
                  DOUBLE  *factor2,
                  DOUBLE  *sig
	           );
void f2_production(
                     DOUBLE  **vderxy,
                     DOUBLE  *production
	           );

void f2_eddyirans(
                     ELEMENT    *eleke,
                     DOUBLE     *eddyint,
                     DOUBLE     *funct,
                     DOUBLE     *eddy,
                     ARRAY_POSITION *fluid_ipos,
                     INT         iel
	     );
void f2_vel_dc(
                  DOUBLE  *velint,
                  DOUBLE  *velint_dc,
                  DOUBLE  *kapepsderxy
	         );

void f2_kapepsder(
                  DOUBLE  *kapepsderxy,
                  DOUBLE **derxy,
                  DOUBLE  *kapeps,
                  INT      iel
	       );
void f2_kapepsder2(
                     DOUBLE  *kapepsderxy2,
                     DOUBLE **derxy2,
                     DOUBLE  *kapepsn,
                     INT      iel
	           );

void f2_estifadd_tu(
                     DOUBLE         **estif,
                     DOUBLE         **emass,
                     DOUBLE         **tmp,
                     INT              iel
	          );

/************************************************************************
| f2_calservice_tu_1.c                                                  |
************************************************************************/
void f2_calset_tu_1(

                     ELEMENT         *ele,
                     ELEMENT         *elev,
                     DOUBLE          *kapomen,
                     DOUBLE          *kapomeg,
                     DOUBLE          *kapomepro,
                     DOUBLE          *eddyg,
                     DOUBLE          *eddypro,
                     DOUBLE          *kappan,
                     DOUBLE          *omega,
                     DOUBLE         **evel,
                     DOUBLE         **xyze,
                     ARRAY_POSITION  *ipos
                   );

void f2_kapomei(
                  DOUBLE  *kapomeint,
                  DOUBLE  *funct,
                  DOUBLE  *kapome,
                  INT      iel
	     );

void f2_kapomeder(
                  DOUBLE  *kapomederxy,
                  DOUBLE **derxy,
                  DOUBLE  *kapome,
                  INT      iel
	       );

void f2_kapomeder2(
                     DOUBLE  *kapomederxy2,
                     DOUBLE **derxy2,
                     DOUBLE  *kapomen,
                     INT      iel
	           );

void f2_xi_kappa(
                  DOUBLE  *kapomederxy,
                  DOUBLE  *omegaderxy,
                  DOUBLE   omegaint,
                  DOUBLE  *xi
	           );

void f2_fac_kappa_1(
                     DOUBLE   xi,
                     DOUBLE   eddyint,
                     DOUBLE   kapomeint,
                     DOUBLE   omegaint,
                     DOUBLE   visc,
                     DOUBLE  *factor,
                     DOUBLE  *factor1,
                     DOUBLE  *factor2,
                     DOUBLE  *sig
	           );

void f2_xi_ome(
                  DOUBLE  **vderxy,
                  DOUBLE   kapomeint,
                  DOUBLE  *xi
	           );

void f2_fac_ome(
                  DOUBLE   xi,
                  DOUBLE   ome_proint,
                  DOUBLE   kappanint,
                  DOUBLE   visc,
                  DOUBLE  *factor,
                  DOUBLE  *factor1,
                  DOUBLE  *factor2,
                  DOUBLE  *sig
	           );

void f2_vel_dc_1(
                  DOUBLE  *velint,
                  DOUBLE  *velint_dc,
                  DOUBLE  *kapomederxy
	         );

void f2_kappain(
                  DOUBLE     *kappanint,
                  DOUBLE     *ome_proint,
                  DOUBLE     *funct,
                  DOUBLE     *kappan,
                  DOUBLE     *kapomepro,
                  INT         iel
	     );

/************************************************************************
 | f2_calstabmat.c                                                      |
 ************************************************************************/
void f2_calstabkvv(
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
void f2_calstabkvp(
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
void f2_calstabkvg(
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
void f2_calstabmvv(
                     ELEMENT         *ele,
                     STAB_PAR_GLS    *gls,
                     DOUBLE         **emass,
                     DOUBLE          *velint,
                     DOUBLE          *funct,
                     DOUBLE         **derxy,
                     DOUBLE         **derxy2,
                     DOUBLE           fac,
                     DOUBLE           visc,
                     INT              iel,
                     INT              ihoel
                   );
void f2_calstabkpv(
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
void f2_calstabkpg(
                    STAB_PAR_GLS    *gls,
		    DOUBLE         **estif,
		    DOUBLE          *funct,
		    DOUBLE         **vderxy,
		    DOUBLE         **derxy,
		    DOUBLE           fac,
		    INT              iel
                   );
void f2_calstabkpp(
                     STAB_PAR_GLS    *gls,
                     DOUBLE         **estif,
                     DOUBLE         **derxy,
                     DOUBLE           fac,
                     INT              iel
                   );
void f2_calstabmpv(
                     STAB_PAR_GLS    *gls,
                     DOUBLE         **emass,
                     DOUBLE          *funct,
                     DOUBLE         **derxy,
                     DOUBLE           fac,
                     INT              iel
                   );

/************************************************************************
 | f2_calstabmat_tu.c                                                   |
 ************************************************************************/
void f2_calstabkkapeps(
                        ELEMENT         *ele,
                        ELEMENT         *elev,
                        DOUBLE         **estif,
                        DOUBLE           kapepsint,
                        DOUBLE          *velint,
                        DOUBLE          *velint_dc,
                        DOUBLE           eddyint,
                        DOUBLE          *kapepsderxy,
                        DOUBLE          *funct,
                        DOUBLE         **derxy,
                        DOUBLE         **derxy2,
                        DOUBLE           fac,
                        DOUBLE           visc,
                        DOUBLE           factor,
                        DOUBLE           sig,
                        INT              iel
                 );
void f2_calstabmkapeps(
                        ELEMENT         *ele,
                        DOUBLE         **emass,
                        DOUBLE          *velint,
                        DOUBLE          *velint_dc,
                        DOUBLE          *funct,
                        DOUBLE         **derxy,
                        DOUBLE           fac,
                       INT              iel
                    );

/************************************************************************
 | f2_calstabmat_tu_1.c                                                 |
 ************************************************************************/
void f2_calstabkkapome(
                        ELEMENT         *ele,
                        ELEMENT         *elev,
                        DOUBLE         **estif,
                        DOUBLE           kapomeint,
                        DOUBLE          *velint,
                        DOUBLE          *velint_dc,
                        DOUBLE           eddyint,
                        DOUBLE          *funct,
                        DOUBLE         **derxy,
                        DOUBLE         **derxy2,
                        DOUBLE           fac,
                        DOUBLE           visc,
                        DOUBLE           factor,
                        DOUBLE           sig,
                        INT              iel
                   );


void f2_calstabmkapome(
                        ELEMENT         *ele,
                        DOUBLE         **emass,
                        DOUBLE          *velint,
                        DOUBLE          *velint_dc,
                        DOUBLE          *funct,
                        DOUBLE         **derxy,
                        DOUBLE           fac,
                        INT              iel
                     );

/************************************************************************
 | f2_calstabpar.c                                                      |
 ************************************************************************/
void f2_calstabpar(
	            ELEMENT         *ele,
		    DOUBLE          *velint,
		    DOUBLE           visc,
		    INT              iel,
		    DIS_TYP          typ,
		    INT              iflag
                  );
void f2_calstabpar_hf(
                        ELEMENT          *ele,
                        INT               ngnode,
                        INT              *iedgnod,
                        DOUBLE          **xyze,
                        DOUBLE            det,
                        DOUBLE           *velint,
                        DOUBLE            phiintng,
                        DOUBLE            phiintn,
                        DOUBLE            phiderx
		     );

/************************************************************************
 | f2_calstabpar_tu.c                                                   |
 ************************************************************************/
void f2_calstabpar_tu(
                        ELEMENT         *ele,
                        ELEMENT         *elev,
                        DOUBLE           eddyint,
                        DOUBLE          *velint,
                        DOUBLE          *velint_dc,
                        DOUBLE           visc
                  );

/************************************************************************
 | f2_calstabpar_tu_1.c                                                  |
 ************************************************************************/
void f2_calstabpar_tu_1(
                        ELEMENT         *ele,
                        ELEMENT         *elev,
                        DOUBLE           eddyint,
                        DOUBLE          *velint,
                        DOUBLE          *velint_dc,
                        DOUBLE           visc
                  );
/************************************************************************
 | f2_calstress.c                                                       |
 ************************************************************************/
void f2_calelestress(
                      INT             viscstr,
       	              ELEMENT        *ele,
		      DOUBLE        **evel,
		      DOUBLE         *epre,
		      DOUBLE         *funct,
		      DOUBLE        **deriv,
		      DOUBLE        **derxy,
		      DOUBLE        **vderxy,
		      DOUBLE        **xjm,
		      DOUBLE        **xyze,
		      DOUBLE        **sigmaint,
                      ARRAY_POSITION *ipos,
		      INT             is_relax
		    );
/************************************************************************
 | f2_calsurfrhs.c                                                      |
 ************************************************************************/
void f2_calsurftenfv(
                        DOUBLE   *eforce,
                        DOUBLE   *funct,
                        DOUBLE   *vn,
                        DOUBLE    sigmaint,
                        DOUBLE    facs,
                        DOUBLE    fac,
                        INT       ngnode,
                        INT      *iedgnod
		    );

/************************************************************************
 | f2_caltau.c                                                          |
 ************************************************************************/
void f2_caltau(
	       ELEMENT         *ele,
	       DOUBLE         **xyze,
	       DOUBLE          *funct,
	       DOUBLE         **deriv,
	       DOUBLE         **xjm,
	       DOUBLE         **evelng,
               DOUBLE           visc
              );

/************************************************************************
 | f2_caltimerhs_tu.c                                                   |
 ************************************************************************/
void f2_calgaltfkapeps(
                  DOUBLE          *eforce,
                  DOUBLE           kapepsint,
                  DOUBLE          *velint,
                  DOUBLE           eddyint,
                  DOUBLE          *funct,
                  DOUBLE         **derxy,
                  DOUBLE         **vderxy,
                  DOUBLE          *kapepsderxy,
                  DOUBLE           visc,
                  DOUBLE           fac,
                  DOUBLE           factor,
                  DOUBLE           factor1,
                  DOUBLE           factor2,
                  DOUBLE           sig,
                  DOUBLE           vderxy_12,
                  DOUBLE           production,
                  INT              iel
                  );
void f2_calstabtfkapeps(
                  ELEMENT         *ele,
                  DOUBLE          *eforce,
                  DOUBLE           kapepsint,
       	          DOUBLE          *velint,
       	          DOUBLE          *velint_dc,
                  DOUBLE           eddyint,
                  DOUBLE         **derxy,
                  DOUBLE          *kapepsderxy2,
                  DOUBLE         **vderxy,
                  DOUBLE          *kapepsderxy,
                  DOUBLE           fac,
                  DOUBLE           visc,
                  DOUBLE           factor,
                  DOUBLE           factor1,
                  DOUBLE           factor2,
                  DOUBLE           sig,
                  DOUBLE           vderxy_12,
                  DOUBLE           production,
                  INT              iel
                  );
/************************************************************************
 | f2_caltimerhs_tu_1.c                                                 |
 ************************************************************************/
void f2_calgaltfkapome(
                  DOUBLE          *eforce,
                  DOUBLE           kapomeint,
                  DOUBLE          *velint,
                  DOUBLE           eddyint,
                  DOUBLE          *funct,
                  DOUBLE         **derxy,
                  DOUBLE         **vderxy,
                  DOUBLE          *kapomederxy,
                  DOUBLE           visc,
                  DOUBLE           fac,
                  DOUBLE           factor,
                  DOUBLE           factor1,
                  DOUBLE           factor2,
                  DOUBLE           sig,
                  DOUBLE           production,
                  INT              iel
                  );

void f2_calstabtfkapome(
                        ELEMENT         *ele,
                        DOUBLE          *eforce,
                        DOUBLE           kapomeint,
                        DOUBLE          *velint,
                        DOUBLE          *velint_dc,
                        DOUBLE           eddyint,
                        DOUBLE         **derxy,
                        DOUBLE          *kapomederxy2,
                        DOUBLE         **vderxy,
                        DOUBLE          *kapomederxy,
                        DOUBLE           visc,
                        DOUBLE           fac,
                        DOUBLE           factor,
                        DOUBLE           factor1,
                        DOUBLE           factor2,
                        DOUBLE           sig,
                        DOUBLE           production,
                        INT              iel
                  );

/************************************************************************
 | f2_caltimerhspro_tu.c                                                |
 ************************************************************************/
void f2_calgalprofkapeps(
                           DOUBLE          *eforce,
                           DOUBLE           eddynint,
                           DOUBLE          *funct,
                           DOUBLE           visc,
                           DOUBLE           fac,
                           DOUBLE           factor1,
                           DOUBLE           production,
                           INT              iel
                  );

void f2_calstabprofkapeps(
                           ELEMENT         *ele,
                           DOUBLE          *eforce,
                           DOUBLE           eddynint,
                           DOUBLE          *funct,
                           DOUBLE           visc,
                           DOUBLE           fac,
                           DOUBLE           factor1,
                           DOUBLE           production,
                           DOUBLE           *velint,
                           DOUBLE           *velint_dc,
                           DOUBLE          **derxy,
                           INT              iel
                  );

/************************************************************************
 | f2_caltimerhspro_tu_1.c                                              |
 ************************************************************************/
void f2_calgalprofkapome(
                           DOUBLE          *eforce,
                           DOUBLE           eddynint,
                           DOUBLE          *funct,
                           DOUBLE           fac,
                           DOUBLE           factor1,
                           DOUBLE           production,
                           INT              iel
                  );

void f2_calstabprofkapome(
                           ELEMENT         *ele,
                           DOUBLE          *eforce,
                           DOUBLE           eddynint,
                           DOUBLE          *funct,
                           DOUBLE           fac,
                           DOUBLE           factor1,
                           DOUBLE           production,
                           DOUBLE           *velint,
                           DOUBLE           *velint_dc,
                           DOUBLE          **derxy,
                           INT              iel
                  );
/************************************************************************
 | f2_caltuvisc.c                                                       |
 ************************************************************************/
DOUBLE f2_calvisc(
                  ELEMENT    *ele,
                  DOUBLE     **vderxy
                );
/************************************************************************
 | f2_calvort.c                                                         |
 ************************************************************************/
void f2_calvort(
                  ELEMENT        *ele,
       	          INT             init
               );

/************************************************************************
 | f2_inpele.c                                                          |
 ************************************************************************/
void f2_inp(ELEMENT *ele, INT counter);
void f2_createale(
    ELEMENT *ele0,
    ELEMENT *ele1,
    INT      numele,
    INT      nodeshift);

/************************************************************************
 | f2_inpele_tu.c                                                       |
 ************************************************************************/
void f2tu_dis(
    ELEMENT *ele0,
    ELEMENT *ele1,
    INT      numele,
    INT      numnode);

/************************************************************************
 | f2_intg.c                                                            |
 ************************************************************************/
void f2_intg(INT                option
	    );
DOUBLE f2_rsn(
               INT            node,
               INT             irs,
               INT             iel
	    );

/************************************************************************
 | f2_int_usfem.c                                                       |
 ************************************************************************/
void f2_int_usfem(
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
	              DOUBLE          *epren,
	              DOUBLE          *edeadng,
	              DOUBLE         **vderxy,
                      DOUBLE         **vderxy2,
                      DOUBLE           visc,
	              DOUBLE         **wa1,
	              DOUBLE         **wa2,
                      DOUBLE           estress[3][MAXNOD_F2],
                      INT              is_relax
	             );
void f2_int_res(
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
	        DOUBLE         **wa2,
                DOUBLE           estress[3][MAXNOD_F2]
	       );
void f2_int_stress_project(
	                   ELEMENT         *ele,
                           INT             *hasext,
                           DOUBLE         **estif,
                           DOUBLE          *force,
                           DOUBLE         **xyze,
                           DOUBLE          *funct,
                           DOUBLE         **deriv,
                           DOUBLE         **xjm,
                           DOUBLE         **derxy,
                           DOUBLE         **evelng,
                           DOUBLE         **vderxy
	                  );
void f2_get_tau(ELEMENT  *ele,
                DOUBLE  **xjm,
                DOUBLE  **xyze,
                DOUBLE   *funct,
                DOUBLE    det,
                DOUBLE   *velint,
                DOUBLE    visc,
                INT       whichtau,
                INT       which_hk);

void f2_int_kim_moin_err(
    ELEMENT         *ele,
    DOUBLE         **xyze,
    DOUBLE          *funct,
    DOUBLE         **deriv,
    DOUBLE         **xjm,
    DOUBLE         **evelng,
    DOUBLE           visc,
    DOUBLE          *epren,
    CONTAINER       *container
    );

/************************************************************************
 | f2_main.c                                                            |
 ************************************************************************/
void fluid2(
            PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
            ELEMENT     *eleke,
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
 | f2_main_tu.c                                                         |
 ************************************************************************/
void fluid2_tu(
               PARTITION   *actpart,
               INTRA       *actintra,
               ELEMENT     *eleke,
               ELEMENT     *elev,
               ARRAY       *estif_global,
               ARRAY       *emass_global,
               ARRAY       *etforce_global,
               ARRAY       *eiforce_global,
               ARRAY       *edforce_global,
               ARRAY       *eproforce_global,
               CALC_ACTION *action,
               INT         *hasdirich,
               INT         *hasext,
               CONTAINER   *container
	   );

/************************************************************************
 | f2_massrhs.c                                            chfoe 09/03  |
 ************************************************************************/
 void f2_massrhs( ELEMENT	 *ele,
                  DOUBLE 	**emass,
                  DOUBLE 	**hist,
                  DOUBLE 	 *edeadng,
                  DOUBLE 	 *eiforce,
                  INT            *hasext);


/************************************************************************
 | f2_mlsubmesh.c                                         gravem 09/03  |
 ************************************************************************/
void f2_pdsubmesh(FLUID_ML_SMESH *smesh,
                  INT             xele,
                  INT             yele,
                  INT             order,
		  INT             flag);
void f2_elesubmesh(ELEMENT        *ele,
                   FLUID_ML_SMESH *smesh,
		   INT             flag);

/************************************************************************
 | f2_restart.c                                                         |
 ************************************************************************/
void f2_write_restart(ELEMENT *actele, INT nhandle, long int *handles);
void f2_read_restart( ELEMENT *actele, INT nhandle, long int *handles);


/*! @} (documentation module close)*/
