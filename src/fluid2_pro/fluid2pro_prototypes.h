/*!----------------------------------------------------------------------
\file
\brief fluid2_pro prototypes

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
 *!---------------------------------------------------------------------
 ************************************************************************
 | f2pro_calbdt.c                                                       |
 ************************************************************************/
void f2pro_calbdt(
                    ELEMENT         *ele,
		    DOUBLE         **estif,
		    DOUBLE          *velint,
		    DOUBLE         **derxy,
		    DOUBLE         **vderxy,
                    DOUBLE          *funct,
                    DOUBLE           fac,
		    DOUBLE          visc,
		    INT               iel
                   );

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
void f2pro_calele(
	        ELEMENT        *elev,
	        ELEMENT        *elep,
                ARRAY          *estif_global,
                ARRAY          *emass_global,
		ARRAY          *lmass_global,
		ARRAY          *gradopr_global,
	        ARRAY          *etforce_global,
		ARRAY          *eiforce_global,
		ARRAY          *edforce_global,
                ARRAY_POSITION *ipos,
		INT             *hasdirich,
		INT             init
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
	        ELEMENT         *elevel,
		ELEMENT         *elepre,
		DOUBLE         **xyze,
                DOUBLE         **eveln,
	        DOUBLE          *epren,
                ARRAY_POSITION  *ipos
	      );
void pro_putdirich_to_dof(FIELD *actfield, INT disnum, DOUBLE scale, INT place);
void pro_putdirich_parabolic_to_dof(FIELD *actfield, INT disnum, DOUBLE scale, INT place);

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
void f2pro_dis(
    ELEMENT *ele0,
    ELEMENT *ele1,
    INT      numele,
    INT      numnode);
void f2pro_ass_dof_q2q1(NODE *actnode, INT *counter);

/************************************************************************
 | f2pro_main.c                                                         |
 ************************************************************************/
void fluid2_pro(     PARTITION     *actpart,
                     INTRA         *actintra,
		     ELEMENT       *elev,
		     ELEMENT       *elep,
		     ARRAY         *estif_global,
		     ARRAY         *emass_global,
		     ARRAY         *lmass_global,
		     ARRAY         *gradopr_global,
		     ARRAY         *etforce_global,
   	             ARRAY         *eiforce_global,
		     ARRAY         *edforce_global,
		     ARRAY         *gforce_global,
		     CALC_ACTION   *action,
                     INT           *hasdirich
	       );
