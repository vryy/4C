/*!----------------------------------------------------------------------
\file
\brief fluid2 prototypes

------------------------------------------------------------------------*/
/* RULE HOW TO ADD NEW FILES AND FUNCTIONS: 
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE!!!
*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
/************************************************************************
 | f2_cacurvature.c                                                     |
 ************************************************************************/
void f2_calq4curv( GLINE                        **actgline, 
                     FLUID_FREESURF_CONDITION   **actlinefs,
		     FLUID_DYN_CALC             *dynvar,
		     ELEMENT                    *actele,
		     int                         foundline,
		     int                         actngline,
		     double                    **xyze,
		     double                    **deriv,
		     double                    **kappa		    
		   );
void f2_calq8curv( GLINE                        **actgline, 
                     FLUID_FREESURF_CONDITION   **actlinefs,
		     FLUID_DYN_CALC             *dynvar,
		     ELEMENT                    *actele,
		     int                         actngline,
		     double                    **xyze,
		     double                    **kappa		    
		   );
		   
 /************************************************************************
 | f2_calele.c                                                          |
 ************************************************************************/
void f2_calele(
                FLUID_DATA     *data, 
                FLUID_DYN_CALC *dynvar, 
	        ELEMENT        *ele,             
                ARRAY          *estif_global,   
                ARRAY          *emass_global,   
	        ARRAY          *etforce_global,       
	        ARRAY          *eiforce_global, 
		ARRAY          *edforce_global,		
		int            *hasdirich,      
                int            *hasext,
                int             imyrank,
		int             init            
	       );
void f2_stress(FLUID_STRESS str, int viscstr ,FLUID_DATA     *data, 
	          ELEMENT        *ele      );
void f2_curvature(
	           FLUID_DATA     *data,
		   FLUID_DYN_CALC *dynvar,  
	           ELEMENT        *ele,
		   int             imyrank
		 );
/************************************************************************
 | f2_calelesize.c                                                      |
 ************************************************************************/
void f2_calelesize(			     
	           ELEMENT         *ele,    
		   FLUID_DATA      *data, 
		   FLUID_DYN_CALC  *dynvar,
	           double         **xyze,
		   double          *funct,  
	           double         **deriv,  
	           double         **deriv2,  		 
		   double         **xjm,    
		   double         **evel,    		  
		   double          *velint, 
		   double         **cutp    
                  );
void f2_calelesize2(			       
	             ELEMENT         *ele,    
	  	     FLUID_DYN_CALC  *dynvar, 
                     double         **xyze,
	             double          *funct,    		   
		     double          *velint, 
		     double         **cutp,   
		     double           visc,   
		     int              iel,    
		     int              ntyp    
                    );
void f2_calstrlen(
                   double   *strle,     
                   double  **xyze,
		   double   *velint,   
		   ELEMENT  *ele,      
                   double   *gcoor,    
		   double  **cutp,             
		   int       ntyp      
                 );

/************************************************************************
 | f2_calelesize.c                                                      |
 ************************************************************************/
void f2_calgalexfv(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,     
		  double          *funct,       
                  double          *edeadn,
		  double          *edeadng,
		  double           fac,      
		  int              iel       
              );
void f2_calstabexfv(
                    FLUID_DYN_CALC  *dynvar, 
                    ELEMENT         *ele,  
                    double          *eforce,     
		    double         **derxy,
		    double         **derxy2,      
                    double          *edead,
	 	    double          *velint,  
		    double           fac,      
                    double           visc,
		    int              iel,
		    int              ihoel,
		    int              flag      
                   );
void f2_calstabexfp(
                    FLUID_DYN_CALC  *dynvar, 
                    double          *eforce,     
		    double         **derxy,       
                    double          *edead,  
		    double           fac,      
		    int              iel,
		    int              flag      
                   );
		   		    
/************************************************************************
 | f2_calfuncderiv.c                                                    |
 ************************************************************************/
void f2_rec(
            double     *funct,     
            double    **deriv,    
            double    **deriv2,   
	    double      r,        
            double      s,        
            DIS_TYP     typ,      
            int         icode     
	   );
void f2_recex(
	      double           *funval,     
	      double           *fpar,
	      double            r,    
	      double            s,
              double           *fval,         
	      int               igauss,
	      int               icode
            );
void f2_tri(
            double     *funct,      
            double    **deriv,    
            double    **deriv2,   
	    double      r,        
            double      s,	  
            DIS_TYP     typ,	  
            int         icode	  
	   );
void f2_triex(
	      double           *funval,     
	      double           *fpar,
	      double            r,    
	      double            s,
              double           *fval,         
	      int               igauss,
	      int               icode
            );
void f2_degrectri(
                    double     *funct, 
                    double    **deriv, 
                    double      r, 
                    DIS_TYP     typ,
                    int         option
		 );
void f2_jaco(
             double    **xyze,
             double     *funct,    
             double    **deriv,   
             double    **xjm,     
             double     *det,          
             int         iel,        
             ELEMENT    *ele
	    );
void f2_jaco2(
             double    **xyze,
             double     *funct,    
             double    **deriv,   
             double    **xjm,     
             double     *det,          
             int         iel,        
             ELEMENT    *ele
	    );
void f2_edgejaco(
                 double    **xyze,
                 double     *funct,    
                 double    **deriv,   
                 double    **xjm,     
                 double     *det,          
                 int         iel,
                 int        *iedgnod
	        );
void f2_gder(
              double   **derxy,     
	      double   **deriv,    
	      double   **xjm,      
	      double     det,      
	      int        iel       
	    );
void f2_gcoor(
              double    **xyze,
	      double     *funct,            
	      int         iel,      
	      double     *gcoor     
             );
void f2_gder2(
               ELEMENT     *ele,
	       double     **xyze,     
	       double     **xjm,      
               double     **bi,     
	       double     **xder2,  
	       double     **derxy,  
	       double     **derxy2, 
               double     **deriv2, 
	       int          iel     
	     );

/************************************************************************
 | f2_calgalmat.c                                                       |
 ************************************************************************/
void f2_calkvv( 
                ELEMENT         *ele,
                FLUID_DYN_CALC  *dynvar,
		double         **estif,   
		double          *velint,
		double          *gridvint,
		double         **vderxy, 
		double          *funct,  
		double         **derxy,  
		double           fac,    
		double           visc,   
		int              iel     
              );
void f2_calkvp(
		double         **estif,   
		double          *funct,  
		double         **derxy,  
		double           fac,    
		int              iel      	     
              );
void f2_calkvg( 
                FLUID_DYN_CALC  *dynvar,
		double         **estif,   
		double         **vderxy, 
		double          *funct,  
		double           fac,    
		int              iel     
              );
void f2_calmvv(
		double         **emass,  
		double          *funct, 
		double           fac,   
		int              iel    
              );
void f2_calkgedge(
	   	   double         **estif,  
		   double          *funct, 
		   double           fac,
		   int             *iedgnod,    		      
		   int              iel,
		   int              ngnode
                );
		
/************************************************************************
 | f2_calint.c                                                          |
 ************************************************************************/
void f2_calint(
               FLUID_DATA      *data,     
	       ELEMENT         *ele,     
	       FLUID_DYN_CALC  *dynvar, 
               int             *hasext,
               double         **estif,   
	       double         **emass,   
	       double          *etforce, 
	       double          *eiforce, 
	       double         **xyze,
	       double          *funct,   
	       double         **deriv,   
	       double         **deriv2,  
	       double         **xjm,     
	       double         **derxy,   
	       double         **derxy2,  
	       double         **eveln,   
	       double         **evelng,  
	       double          *epren,   
	       double          *edeadn,
	       double          *edeadng,	        	       
	       double          *velint,  
	       double          *vel2int, 
	       double          *covint,  
	       double         **vderxy,  
	       double          *pderxy,  
	       double         **vderxy2, 
	       double         **wa1,     
	       double         **wa2      
	      );	      	      	      	     	     	    	   	   
void f2_calinta(
                FLUID_DATA      *data,     
	        ELEMENT         *ele,     
	        FLUID_DYN_CALC  *dynvar, 
                int             *hasext,
                int              imyrank,
                double         **estif,   
	        double         **emass,   
	        double          *etforce, 
	        double          *eiforce, 
	        double         **xyze,
	        double          *funct,   
	        double         **deriv,   
	        double         **deriv2,  
	        double         **xjm,     
	        double         **derxy,   
	        double         **derxy2,  
	        double         **eveln,   
	        double         **evelng, 
		double         **ealecovn,
		double         **ealecovng,
		double         **egridv, 
	        double          *epren,   
	        double          *edeadn,
	        double          *edeadng,	        	       
	        double          *velint,  
	        double          *vel2int,
		double          *alecovint,
		double          *gridvint, 
	        double          *covint,  
	        double         **vderxy,  
	        double          *pderxy,  
	        double         **vderxy2,
		double          *kappan,
		double          *kappang, 
	        double         **wa1,     
	        double         **wa2      
	      );
	      
/************************************************************************
 | f2_caliterrhs.c                                                      |
 ************************************************************************/ 
void f2_calgalifv(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,   
		  double          *covint,  
		  double          *funct,   
		  double           fac,     
		  int              iel      
                 ) ;
void f2_calstabifv(
                  FLUID_DYN_CALC  *dynvar, 
                  ELEMENT         *ele,      
		  double          *eforce,  
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
void f2_calstabifp(
                   FLUID_DYN_CALC  *dynvar, 
                   double          *eforce,   
		   double          *covint,  
		   double         **derxy,   
		   double           fac,     
		   int              iel        
                 );

/************************************************************************
 | f2_calservice.c                                                      |
 ************************************************************************/
void f2_calset( 
                FLUID_DYN_CALC  *dynvar, 
	        ELEMENT         *ele,     
                double         **xyze,
                double         **eveln,    
	        double         **evelng, 
	        double          *epren,
		double          *edeadn,
		double          *edeadng,
		int             *hasext
	      );
void f2_calseta( 
                FLUID_DYN_CALC  *dynvar, 
	        ELEMENT         *ele,     
                double         **xyze,
                double         **eveln,    
	        double         **evelng,
		double         **ealecovn,
		double         **ealecovng,
		double         **egridv, 
	        double          *epren,
		double          *edeadn,
		double          *edeadng,
		double          *ekappan,
		double          *ekappang,
		int             *hasext
	      );
void f2_alecoor( 
                FLUID_DYN_CALC  *dynvar, 
	        ELEMENT         *ele,     
                double         **xyze
	       );
void f2_veli(
             double  *velint,     
             double  *funct,    
	     double **evel,     
	     int      iel       
	    );
void f2_prei(
             double  *preint,     
             double  *funct,    
	     double  *epre,     
	     int      iel       
	    );
double f2_kappai(  
		 double  *funct,    
	         double  *ekappa,  
	         int     *iedgnod,   
	         int      ngnode       
	        );
void f2_vder(
             double **vderxy,     
             double **derxy,    
	     double **evel,     
	     int      iel       
	    );
void f2_vder2(
             double **vderxy2,    
             double **derxy2,    
	     double **evel,      
	     int      iel        
             );
void f2_pder(
             double  *pderxy,     
             double **derxy,    
	     double  *epre,     
	     int      iel        
	    );
void f2_covi(
             double **vderxy,    
             double  *velint,   
	     double  *covint    
	    );
void f2_permeforce( 
		   double   *eforce,  
		   double  **tmp,    
		   int       iel     
	          );
void f2_permeforce_ifs( 
		       double   *eforce,  
		       double  **tmp,    
		       ELEMENT  *ele     
	              );		  
void f2_permestif(                  
		   double         **estif,   
		   double         **emass, 
		   double         **tmp,   
		   ELEMENT         *ele,   
		   FLUID_DYN_CALC  *dynvar		   		    
	          );
void f2_permestif_ifs(                  
		      double         **estif,   
		      double         **emass, 
		      double         **tmp,   
		      ELEMENT         *ele,   
		      FLUID_DYN_CALC  *dynvar		   		    
	             ); 
void f2_iedg(     
                int     *iegnod, 
		ELEMENT *ele, 
		int      line, 
		int      init
	     );
	     
/************************************************************************
 | f2_calstabmat.c                                                      |
 ************************************************************************/
void f2_calstabkvv(			      
                    ELEMENT         *ele,    
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,  
		    double          *velint,
		    double          *vel2int, 
		    double          *gridvint,
		    double         **vderxy, 
		    double          *funct,  
		    double         **derxy,  
		    double         **derxy2, 
		    double           fac,    
		    double           visc,   
		    int              iel,    
                    int              ihoel   
                   );
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
void f2_calstabkvg(			      
                    ELEMENT         *ele,    
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,  
		    double         **vderxy, 
		    double          *funct,  
                    double         **derxy,
		    double         **derxy2, 
                    double          *alecovint,
		    double           fac,    
		    double           visc,   
		    int              iel,    
                    int              ihoel   
                   );
void f2_calstabmvv(
                    ELEMENT         *ele,     
		    FLUID_DYN_CALC  *dynvar,
		    double         **emass,  
		    double          *velint, 
    		    double          *funct,  
		    double         **derxy,  
		    double         **derxy2, 
		    double           fac,    
		    double           visc,   
		    int              iel,    
		    int              ihoel           
                   );
void f2_calstabkpv(
		    ELEMENT         *ele,
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,   
		    double          *velint,
		    double          *gridvint, 
		    double         **vderxy, 
		    double          *funct,  
		    double         **derxy,  
		    double         **derxy2, 
		    double           fac,    
		    double           visc,   
		    int              iel,    
		    int              ihoel          
                   );
void f2_calstabkpg(
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif, 
		    double          *funct,  
		    double         **vderxy, 
		    double         **derxy,  
		    double           fac,    
		    int              iel          
                   );
void f2_calstabkpp(
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,   
		    double         **derxy,  
		    double           fac,    
		    int              iel             
                   );
void f2_calstabmpv(
		    FLUID_DYN_CALC  *dynvar,
		    double         **emass,   
		    double          *funct,  
		    double         **derxy,  
		    double           fac,    
		    int              iel     
                   );

/************************************************************************
 | f2_calstabpar.c                                                      |
 ************************************************************************/ 
void f2_calstabpar(
	            ELEMENT         *ele,      
		    FLUID_DYN_CALC  *dynvar,
		    double          *velint,  
		    double           visc,    
		    int              iel,     
		    int              ntyp,    
		    int              iflag    
                  );
		  
/************************************************************************
 | f2_calstress.c                                                       |
 ************************************************************************/
void f2_calfsistress(
                      int             viscstr,
		      FLUID_DATA     *data, 
       	              ELEMENT        *ele,
		      double        **evel, 
		      double         *epre,
		      double         *funct,
		      double        **deriv,
		      double        **derxy,
		      double        **vderxy,
		      double        **xjm,
		      double        **xyze,
		      double        **sigmaint
		    );
/************************************************************************
 | f2_calsurfrhs.c                                                      |
 ************************************************************************/
void f2_calsurftenfv( 
                     double   *eforce, 
		     double   *funct, 
		     double   *vn, 
		     double    sigmaint,
		     double    facs, 
		     double    fac,
                     int       ngnode,
		     int      *iedgnod     
		    );
		    
/************************************************************************
 | f2_caltimerhs.c                                                      |
 ************************************************************************/
void f2_calgaltfv(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,    
		  double          *vel2int,    
		  double          *covint,   
		  double          *funct,    
		  double         **derxy,    
		  double         **vderxy,   
		  double           preint,   
		  double           visc,     
		  double           fac,      
		  int              iel       
              ) ;
void f2_calgaltfp(
                   FLUID_DYN_CALC  *dynvar,    
                   double          *eforce,   
		   double          *funct,    
		   double         **vderxy,   
		   double           fac,      
		   int              iel       
                  );
void f2_calstabtfv(
                   FLUID_DYN_CALC  *dynvar,
                   ELEMENT         *ele,      
	           double          *eforce,  
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
void f2_calstabtfp(
                   FLUID_DYN_CALC  *dynvar, 
                   double          *eforce,    
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

/************************************************************************
 | f2_calvort.c                                                         |
 ************************************************************************/
void f2_calvort(
                FLUID_DATA     *data, 
                FLUID_DYN_CALC *dynvar, 
	        ELEMENT        *ele,                
       	        int             init            
               );
	       
/************************************************************************
 | f2_inpele.c                                                          |
 ************************************************************************/
void f2_inp(ELEMENT *ele, int counter);

/************************************************************************
 | f2_intg.c                                                            |
 ************************************************************************/
void f2_intg(FLUID_DATA         *data,
             int                option  
	    );
double f2_rsn(
	      int            node,     
	      int             irs,    
	      int             iel       
	    );

/************************************************************************
 | f2_main.c                                                            |
 ************************************************************************/
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
	    int         *hasdirich,
	    int         *hasext,       
	    CONTAINER   *container
            );
/*! @} (documentation module close)*/	    
