/*----------------------------------------------------------------------*
 |  f2_inpele.c                                            genk 3/02    |
 *----------------------------------------------------------------------*/
void f2_inp(ELEMENT *ele);
/*----------------------------------------------------------------------*
 | main fluid2  control routine                              genk 03/02 |
 *----------------------------------------------------------------------*/
void fluid2(
            PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
            ARRAY       *estif_global,
            ARRAY       *emass_global, 
	    ARRAY       *etforce_global,
	    ARRAY       *eiforce_global,
	    ARRAY       *edforce_global,
            CALC_ACTION *action,
	    int         *hasdirich
	   );
#ifdef D_FLUID2
/*----------------------------------------------------------------------*
 | integration points                                      genk 6/01    |
 | this routine is a try to organise the integration parameters         |
 | different. ALL paramters are saved in F2_DATA, so that this routine  |
 | has to be (hopefully) called only once!!!                            |
 *----------------------------------------------------------------------*/
void f2_intg(const ELEMENT   *ele,
               F2_DATA         *data,
               int              option);
/*----------------------------------------------------------------------*
 |  control routine for element integration                             |
 | at the moment there's the question if this routine is really         |
 | necessary. it may possible to call directly the integration routine  |
 | 
 |                                                           genk 03/02 |
 *----------------------------------------------------------------------*/
void f2_calele(
                F2_DATA        *data, 
                FLUID_DYN_CALC *dynvar, 
	        ELEMENT        *ele,
                ARRAY          *estif_global,
                ARRAY          *emass_global, 
	        ARRAY          *etforce_global,
	        ARRAY          *eiforce_global,
		ARRAY          *edforce_global,	
		int            *hasdirich,
		int             init
	       );
       
/*----------------------------------------------------------------------*
 |  set all arrays for element calculation                              |
 |  get all the element velocities / pressure at different times        |
 |  NOTE: in contradiction to the old programm the kinematic pressure   |
 |        is stored in the solution history; so one can avoid the       |
 |        transformation in every time step                             |
 |                                                           genk 04/02 |
 *----------------------------------------------------------------------*/
void f2_calset( 
                FLUID_DYN_CALC  *dynvar, 
	        ELEMENT         *ele,
                double         **eveln,
	        double         **evelng,
	        double          *epren
	      );
/*----------------------------------------------------------------------*
 | integration loop for one fluid element                               |
 |                                                                      |
 |                                                           genk 03/02 |
 *----------------------------------------------------------------------*/
void f2_calint(
               F2_DATA         *data, 
	       ELEMENT         *ele,
	       FLUID_DYN_CALC  *dynvar,
               double         **estif,
	       double         **emass,
	       double          *etforce,
	       double          *eiforce,
	       double          *funct,
	       double         **deriv,
	       double         **deriv2,
	       double         **xjm,
	       double         **derxy,
	       double         **derxy2,
	       double         **eveln,
	       double         **evelng,
	       double          *epren,
	       double          *velint,
	       double          *vel2int,
	       double          *covint,
	       double         **vderxy,
	       double          *pderxy,
	       double         **vderxy2,
	       double         **wa1,
	       double         **wa2
	      );
/*----------------------------------------------------------------------*
 | routine for evaluation of shpae functions and their natural first    |
 | and secend derivatives with respect to r/s for                       |
 | R E C T A N G L E S                                     genk 04/02   |
 *----------------------------------------------------------------------*/
void f2_rec(
            double     *funct, 
            double    **deriv, 
            double    **deriv2,
	    double      r, 
            double      s,
            DIS_TYP     typ,
            int         icode
	   );
/*----------------------------------------------------------------------*
 | routine for evaluation of shpae functions and their natural first    |
 | and secend derivatives with respect to r/s for                       |
 | T R I A N G L E S                                       genk 06/02   |
 *----------------------------------------------------------------------*/
void f2_tri(
            double     *funct, 
            double    **deriv, 
            double    **deriv2,
	    double      r, 
            double      s,
            DIS_TYP     typ,
            int         icode
	   );	
/*----------------------------------------------------------------------*
 | calculate operator matrix at point r,s                   genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_jaco(
             double     *funct,
             double    **deriv,
             double    **xjm,
             double     *det,
             ELEMENT    *ele,
             int         iel
	    );	  
/*----------------------------------------------------------------------*
 | routine to calculate the globael derivatives wrt x,y at point r,s    |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_gder(
              double   **derxy,    /* global derivatives wrt. x/y       */
	      double   **deriv,    /* derivatives of shape functions    */
	      double   **xjm,      /* jacobian matrix                   */
	      double     det,      /* jacobian determinant              */
	      int        iel       /* number of nodes in actual element */
	    );	
/*----------------------------------------------------------------------*
 | routine to calculate second global derivatives wrt x/y at point r,s  |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_gder2(
               ELEMENT     *ele,
	       double     **xjm,            
               double     **bi,
	       double     **xder2,
	       double     **derxy,
	       double     **derxy2,
               double     **deriv2,
	       int          iel
	     );
/*---------------------------------------------------------------------*
 | routine to calculate velocities derivativs at integration point     |  
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_vder(
             double **vderxy,   /* velocity derivativs             */ 
             double **derxy,    /* globael derivatives             */
	     double **evel,     /* velocites at element nodes      */
	     int     iel        /* number of nodes in this element */
	    );	         	         	      
/*---------------------------------------------------------------------*
 | routine to 2nd calculate velocity derivatives at integration point  |  
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_vder2(
             double **vderxy2,   /* velocity derivativs             */ 
             double **derxy2,    /* globael derivatives             */
	     double **evel,     /* velocites at element nodes      */
	     int      iel       /* number of nodes in this element */
	    ) ;
/*---------------------------------------------------------------------*
 | routine to calculate velocities at integration point                |  
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_veli(
             double  *velint,   /* velocities at integration point */ 
             double  *funct,    /* shape functions                 */
	     double **evel,     /* velocites at element nodes      */
	     int      iel       /* number of nodes in this element */
	    );
/*---------------------------------------------------------------------*
 | routine to calculate convective velocities at integration point     |
 | u * grad(u)                                                         |
 | e.g. 2D: COVx = Ux*Ux,x + Uy*Ux,y                                   |
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_covi(
             double **vderxy,   /* velocity derivativs                 */ 
             double  *velint,   /* velocity at integration point       */
	     double  *covint    /* convective velocity at int point    */ 
	    ) ;	    
/*---------------------------------------------------------------------*
 | routine to calculate pressure at integration point                |  
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_prei(
             double  *preint,   /* pressure at integration point */ 
             double  *funct,    /* shape functions                 */
	     double  *epre,     /* pressure at element nodes      */
	     int      iel       /* number of nodes in this element */
	    ) ;	    
/*---------------------------------------------------------------------*
 | routine to calculate pressure derivatives at integration point      |  
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_pder(
             double  *pderxy,   /* pressure derivativs             */ 
             double **derxy,    /* globael derivatives             */
	     double  *epre,     /* pressure at element nodes      */
	     int      iel       /* number of nodes in this element */
	    );	    
/*----------------------------------------------------------------------*
 | routine to calculate matrix Kvv                         genk 04/02   |
 |    NOTE: there's only one elestif                                    |
 |          --> Kvv is stored in estif[0..(2*iel-1)][0..(2*iel-1)]      |
 *----------------------------------------------------------------------*/
void f2_calkvv(
                FLUID_DYN_CALC  *dynvar,
		double         **estif,
		double          *velint,
		double         **vderxy,
		double          *funct,
		double         **derxy,
		double           fac,
		double           visc,
		int              iel		 
              );
/*----------------------------------------------------------------------*
 | routine to calculate matrices Kvp and Kpv               genk 04/02   |
 | NOTE: there's only one elestif                                       |
 |       --> Kvp is stored in estif[(0..(2*iel-1)][(2*iel)..(3*iel-1)]  |
 |	 --> Kpv is stored in estif[((2*iel)..(3*iel-1)][0..(2*iel-1)]  |
 *----------------------------------------------------------------------*/

void f2_calkvp(
		double         **estif,
		double          *funct,
		double         **derxy,
		double           fac,
		int              iel                
              );	       	    
/*----------------------------------------------------------------------*
 | routine to calculate matrix Mvv                         genk 04/02   |
 |    NOTE: there's only one elestif                                    |
 |          --> Mvv is stored in estif[0..(2*iel-1)][0..(2*iel-1)]      |
 *----------------------------------------------------------------------*/

void f2_calmvv(
		double         **estif,
		double          *funct,
		double           fac,
		int              iel  
              );	    
/*----------------------------------------------------------------------*
 | routine to calculate elment sizes and stabilisation parameter        |
 | (at center) for one element                             genk 04/02   |
 |									|
 *----------------------------------------------------------------------*/
void f2_calelesize(
	           ELEMENT         *ele,
		   F2_DATA         *data,
		   FLUID_DYN_CALC  *dynvar,
	           double          *funct,
	           double         **deriv,
	           double         **deriv2,	       
		   double         **xjm,
		   double         **evel,	       
		   double          *velint,
		   double         **cutp		
                  );

/*----------------------------------------------------------------------*
 | routine to calculate elment size and stabilisation parameter     |
 | (at center) for one element during integration loop     genk 04/02   |
 *----------------------------------------------------------------------*/
 void f2_calelesize2(
	             ELEMENT         *ele,
	  	     FLUID_DYN_CALC  *dynvar,
	             double          *funct,	       	       
		     double          *velint,
		     double         **cutp,
		     double           visc,
		     int              iel,
		     int              ntyp	
                    );		  	   	    
/*----------------------------------------------------------------------*
 | routine to calculate streamlength                                    |
 | (only approxmation for higher order elements - boundaries are        |
 | assumed to be straight)                                   genk 04/02 |
 *----------------------------------------------------------------------*/
void f2_calstrlen(
                   double   *strle,    /* streamlength  */
		   double   *velint,   /* velocities at integr. point */
		   ELEMENT  *ele,      /* actual element */
		   double   *gcoor,    /* global coordinates of int. point */
		   double  **cutp,     /* cutting points */       
		   int       ntyp      /* flag for element type */
                 );
/*----------------------------------------------------------------------*
 | routine to calculate stability parameter                genk 04/02   |
 *----------------------------------------------------------------------*/
void f2_calstabpar(
	            ELEMENT         *ele,
		    FLUID_DYN_CALC  *dynvar,
		    double          *velint,
		    double           visc,    /* viscosity */
		    int              iel,     /* number of nodes */
		    int              ntyp,
		    int              iflag    /* flag for evaluation */
                  );	
/*----------------------------------------------------------------------*
 | routine to calculate stbilisation matrix Kvv            genk 04/02   |
 |    NOTE: there's only one elestif                                    |
 |          --> Kvv is stored in estif[0..(2*iel-1)][0..(2*iel-1)]      |
 *----------------------------------------------------------------------*/
void f2_calstabkvv(
                    ELEMENT         *ele,
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,
		    double          *velint,
		    double         **vderxy,
		    double          *funct,
		    double         **derxy,
		    double         **derxy2,
		    double           fac,
		    double           visc,
		    int              iel,	 
                    int              ihoel
                   );
/*----------------------------------------------------------------------*
 | routine to calculate stbilisation matrix Kvp            genk 04/02   |
 | NOTE: there's only one elestif                                       |
 |       --> Kvp is stored in estif[(0..(2*iel-1)][(2*iel)..(3*iel-1)]  |
 *----------------------------------------------------------------------*/
void f2_calstabkvp(
                    ELEMENT         *ele,
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,
		    double          *velint,
		    double          *funct,
		    double         **derxy,
		    double         **derxy2,
		    double           fac,
		    double           visc,
		    int              iel,
		    int              ihoel		 
                   );
/*----------------------------------------------------------------------*
 | routine to calculate stbilisation matrix Mvv            genk 04/02   |
 |    NOTE: there's only one elestif                                    |
 |          --> Mvv is stored in estif[0..(2*iel-1)][0..(2*iel-1)]      |
 *----------------------------------------------------------------------*/
void f2_calstabmvv(
                    ELEMENT         *ele,
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,
		    double          *velint,
    		    double          *funct,
		    double         **derxy,
		    double         **derxy2,
		    double           fac,
		    double           visc,
		    int              iel,
		    int              ihoel		 
                   );
/*----------------------------------------------------------------------*
 | routine to calculate stbilisation matrix Kvp            genk 04/02   |
 | NOTE: there's only one elestif                                       |
 |	 --> Kpv is stored in estif[((2*iel)..(3*iel-1)][0..(2*iel-1)]  |
 *----------------------------------------------------------------------*/
void f2_calstabkpv(
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,
		    double          *velint,
		    double         **vderxy,
		    double          *funct,
		    double         **derxy,
		    double         **derxy2,
		    double           fac,
		    double           visc,
		    int              iel,
		    int              ihoel		 
                   );
/*----------------------------------------------------------------------*
 | routine to calculate stbilisation matrix Kpp            genk 04/02   |
 | NOTE: there's only one elestif                                       |
 |	 --> Kpp is stored in                                           |
 |               estif[((2*iel)..(3*iel-1)][((2*iel)..(3*iel-1)]        |
 *----------------------------------------------------------------------*/
void f2_calstabkpp(
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,
		    double         **derxy,
		    double           fac,
		    int              iel		 
                   );
/*----------------------------------------------------------------------*
 | routine to calculate stbilisation matrix Mpv            genk 04/02   |
 | NOTE: there's only one elestif                                       |
 |	 --> Mpv is stored in estif[((2*iel)..(3*iel-1)][0..(2*iel-1)]  |
 *----------------------------------------------------------------------*/
void f2_calstabmpv(
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,
		    double          *funct,
		    double         **derxy,
		    double           fac,
		    int              iel		 
                   );
/*----------------------------------------------------------------------*
 | routine to calculate galerkin part of time forces for vel dofs       |
 | NOTE:							        |
 |    in ONESTEP methods: velint  = vel2int = U(n)		        |
 |    in TWOSTEP methods: velint  = U(n+gamma)			        |
 |    in TWOSTEP methods: vel2int = U(n)			        |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calgaltfv(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,   /* element force vector */
		  double          *velint,
		  double          *vel2int,
		  double          *covint,
		  double          *funct,
		  double         **derxy,
		  double         **vderxy,
		  double           preint,
		  double           visc,
		  double           fac,
		  int              iel
              );
/*----------------------------------------------------------------------*
 | routine to calculate galerkin part of time forces for pre dofs       |
 | NOTE:                                                                |
 |     there's only one full element force vector                       |  
 |     for pre-dofs the pointer eforce points to the entry              |
 |     eforce[2*iel]                                                    |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calgaltfp(
                   FLUID_DYN_CALC  *dynvar, 
                   double          *eforce,   /* element force vector */
		   double          *funct,
		   double         **vderxy,
		   double           fac,
		   int              iel
                  ) ;
/*----------------------------------------------------------------------*
 | routine to calculate stabilisation part of time forces for vel dofs  |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calstabtfv(
                   FLUID_DYN_CALC  *dynvar, 
                   ELEMENT         *ele,
	           double          *eforce,   /* element force vector */
	 	   double          *velint,
		   double          *vel2int,
		   double          *covint,
		   double         **derxy,
		   double         **derxy2,
		   double         **vderxy,
		   double         **vderxy2,
		   double          *pderxy,
		   double           fac,
		   double           visc,
		   int              ihoel,
		   int              iel
                  );
/*----------------------------------------------------------------------*
 | routine to calculate stabilisation part of time forces for pre dofs  |
 | NOTE:                                                                |
 |     there's only one full element force vector                       |  
 |     for pre-dofs the pointer eforce points to the entry              |
 |     eforce[2*iel]                                                    |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calstabtfp(
                   FLUID_DYN_CALC  *dynvar, 
                   double          *eforce,   /* element force vector */ 
     		   double         **derxy,
		   double         **vderxy2,
		   double          *velint,
		   double          *covint,
		   double          *pderxy,
		   double           visc,
		   double           fac,
		   int              ihoel,
		   int              iel
                  );
/*----------------------------------------------------------------------*
 | routine to calculate galerkin part of iteration forces for vel dofs  |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calgalifv(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,   /* element force vector */
		  double          *covint,
		  double          *funct,
		  double           fac,
		  int              iel
                 );
/*----------------------------------------------------------------------*
 | routine to calculate stabilisation part of iteration forces          |
 | for vel dofs                                             genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calstabifv(
                  FLUID_DYN_CALC  *dynvar, 
                  ELEMENT         *ele,
		  double          *eforce,   /* element force vector */
		  double          *covint,
		  double          *velint,
		  double          *funct,
		  double         **derxy,
		  double         **derxy2,
		  double           fac,
		  double           visc,
		  int              ihoel,
		  int              iel
                 );
/*----------------------------------------------------------------------*
 | routine to calculate stabilisation part of iteration forces          |
 | for pre dofs                                                         |
 | NOTE:                                                                |
 |     there's only one full element force vector                       |  
 |     for pre-dofs the pointer eforce points to the entry              |
 |     eforce[2*iel]                                                    |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calstabifp(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,   /* element force vector */
		  double          *covint,
		  double          **derxy,
		  double           fac,
		  int              iel
                 );	
/*---------------------------------------------------------------------*
 | routine to rearrange the element force vector                       |
 | this is necessary since we would like to use the existing assembly  |
 | routines for the RHS                                                |
 | hence a splitting of vel- and pre dof in the element force vector   |
 | is not possible any more!!!!                                        |
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_permeforce(
		   double    *eforce,
		   double   **tmp,
		   int        iel		   		   
	          ) ;	
/*---------------------------------------------------------------------*
 | routine to add galerkin and stabilisation parts of the elment       |
 | stiffness matrix and to rearrange its entries!		       |
 | this is necessary since we would like to use the existing assembly  |
 | routines for the stiffness matrix                                   |
 | hence a splitting of vel- and pre dofs is not possible any more!!!! |  				       |
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_permestif(                  
		   double         **estif,
		   double         **emass,
		   double         **tmp,
		   int              iel,
		   FLUID_DYN_CALC  *dynvar		   		   
	          );	
/*----------------------------------------------------------------------*
 |  routine to calculate the element dirichlet load vector              |
 |                                                  genk 04/02          |
 *----------------------------------------------------------------------*/
void f2_caldirich(
                  ELEMENT   *actele, 
		  double    *dforces,
                  double   **estif, 
		  int       *hasdirich
		 ) ;		  	  	 	 		 			      			      		   		   		   		   		   		   		  	 
#endif
