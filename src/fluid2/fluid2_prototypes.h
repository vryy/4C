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
                ELEMENT        *eleke, 
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
 | f2_calele_tu.c                                                          |
 ************************************************************************/
void f2_calele_tu(
                FLUID_DATA     *data, 
                FLUID_DYN_CALC *dynvar, 
	          ELEMENT        *eleke,             
                ELEMENT        *elev,          
                ARRAY          *estif_global,   
                ARRAY          *emass_global,   
	          ARRAY          *etforce_global,       
	          ARRAY          *eiforce_global, 
                ARRAY          *edforce_global,		
                ARRAY          *eproforce_global,		
                int            *hasdirich,      
                int            *hasext,
		    int             init            
                );

/************************************************************************
 | f2_calele_tu_1.c                                                          |
 ************************************************************************/
void f2_calele_tu_1(
                FLUID_DATA     *data, 
                FLUID_DYN_CALC *dynvar, 
	          ELEMENT        *eleke,             
                ELEMENT        *elev,          
                ARRAY          *estif_global,   
                ARRAY          *emass_global,   
	          ARRAY          *etforce_global,       
	          ARRAY          *eiforce_global, 
                ARRAY          *edforce_global,		
                ARRAY          *eproforce_global,		
                int            *hasdirich,      
                int            *hasext,
		    int             init            
	         );
/************************************************************************
 | f2_calelesize.c                                                      |
 ************************************************************************/
void f2_calelesize(			     
	           ELEMENT         *ele,    
                 ELEMENT         *eleke,    
		   FLUID_DATA      *data, 
		   FLUID_DYN_CALC  *dynvar,
	           double         **xyze,
		   double          *funct,  
	           double         **deriv,  
	           double         **deriv2,  		 
		   double         **xjm,    
               double         **derxy, 
               double         **vderxy, 
		   double         **evel,    		  
		   double          *velint, 
		   double         **cutp,    
               double          *eddy,
               double          *visc    
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
 | f2_calelesize_tu.c                                                      |
 ************************************************************************/
void f2_calelesize_tu(			     
	           ELEMENT         *ele,    
		     ELEMENT         *elev,    
                 FLUID_DYN_CALC  *dynvar,
		     FLUID_DATA      *data, 
	           double          *funct,  
	           double         **deriv,  
	           double         **deriv2,  
                 double         **evel, 
                 double          *eddyg, 
                 double          *velint, 
                 double          *velint_dc, 
                 double          *kapepsn, 
	           double         **xjm,     
	           double         **xyze,     
	           double         **derxy,   
                 double          *kapepsderxy,  
                 double         **cutp
                 );
		   		    
void f2_calstrlen_tu(
		   double   *velint_dc,   
		   ELEMENT  *ele,      
               double   *gcoor,    
		   double  **cutp,             
		   int       ntyp      
                 );
                 
/************************************************************************
 | f2_calelesize_tu_1.c                                                  |
 ************************************************************************/
void f2_calelesize_tu_1(			     
	           ELEMENT         *ele,    
		     ELEMENT         *elev,    
                 FLUID_DYN_CALC  *dynvar,
		     FLUID_DATA      *data, 
	           double          *funct,  
	           double         **deriv,  
	           double         **deriv2,  		 
                 double         **evel, 
                 double          *eddyg, 
                 double          *velint, 
                 double          *velint_dc, 
                 double          *kapomen, 
	           double         **xjm,     
	           double         **xyze,     
	           double         **derxy,   
                 double          *kapomederxy,  
                 double         **cutp
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
 | f2_calgalmat_tu.c                                                       |
 ************************************************************************/
void f2_calkkapeps(
                FLUID_DYN_CALC  *dynvar,
		    double         **estif,   
		    double           kapepsint, 
	          double          *velint,
                double           eddyint, 
                double          *kapepsderxy, 
                double          *funct,  
		    double         **derxy,  
		    double           fac,    
		    double           visc,   
		    double           factor,
		    double           sig,
                int              iel     
              );
void f2_calmkapeps(
		double         **emass,  
		double          *funct, 
		double           fac,   
		int              iel    
              );

/************************************************************************
 | f2_calgalmat_tu_1.c                                                       |
 ************************************************************************/
void f2_calkkapome(
                FLUID_DYN_CALC  *dynvar,
		    double         **estif,   
		    double           kapomeint, 
	          double          *velint,
                double           eddyint, 
                double          *funct,  
		    double         **derxy,  
		    double           fac,    
		    double           visc,   
		    double           factor,
		    double           sig,
                int              iel     
              );

void f2_calmkapome(
		double         **emass,  
		double          *funct, 
		double           fac,   
		int              iel    
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
	       double         **wa2,      
             double           visc      
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
 | f2_calint_tu.c                                                          |
 ************************************************************************/
void f2_calint_tu(
               FLUID_DATA      *data,     
	         ELEMENT         *eleke,     
	         ELEMENT         *elev, 
               FLUID_DYN_CALC  *dynvar, 
	         double         **estif,   
	         double         **emass,   
	         double          *etforce, 
	         double          *eiforce, 
	         double          *eproforce, 
	         double          *funct,   
	         double         **deriv,   
	         double         **deriv2,  
	         double         **xjm,     
	         double         **xyze,
	         double         **derxy,   
	         double         **derxy2,  
	         double          *kapepsn,   
	         double          *kapepsg,   
               double          *eddyg, 
               double          *eddypro, 
               double          *kappa, 
               double          *kappan, 
               double          *epsilon,    
               double          *kapepspro,    
               double          *kapepsderxy,  
               double          *kapepsderxy2, 
	         double          *velint,  
	         double          *velint_dc,  
               double         **evel,
               double         **vderxy,
               double         **vderxy2,
               double         **wa1,     
	         double         **wa2      
               );	      	      	      	     	     	    	   	   

/************************************************************************
 | f2_calint_tu_1.c                                                          |
 ************************************************************************/
void f2_calint_tu_1(
               FLUID_DATA      *data,     
	         ELEMENT         *ele,     
	         ELEMENT         *elev, 
               FLUID_DYN_CALC  *dynvar, 
               double         **estif,   
	         double         **emass,   
	         double          *etforce, 
	         double          *eiforce, 
	         double          *eproforce, 
	         double          *funct,   
	         double         **deriv,   
	         double         **deriv2,  
	         double         **xjm,     
	         double         **xyze,
	         double         **derxy,   
	         double         **derxy2,  
	         double          *kapomen,   
	         double          *kapomeg,   
               double          *eddyg, 
               double          *eddypro, 
               double          *kappan, 
               double          *omega,    
               double          *kapomepro,    
               double          *omegaderxy,  
               double          *kapomederxy,  
               double          *kapomederxy2, 
	         double          *velint,  
	         double          *velint_dc,  
               double         **evel,
               double         **vderxy,
               double         **vderxy2,
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
 | f2_caliterrhs_tu.c                                                      |
 ************************************************************************/ 
/*!--------------------------------------------------------------------- 
\brief galerkin part of iteration forces for kapeps dofs

<pre>                                                        he  12/02

In this routine the galerkin part of the iteration forces for kapeps dofs
is calculated:

                    / 
          THETA*dt | factor * (kapeps_old)^2 * psi  d_omega
                  /  


LOW-REYNOLD's MODEL only for epsilon:
    
                   / 
   (+)   THETA*dt | 2.0*visc*nue_t* (vderxy2_12)^2 * psi  d_omega
                 /  

      
</pre>
\param  *dynvar      FLUID_DYN_CALC  (i)
\param  *eforce      double	    (i/o)   element force vector
\param   eddyint     double	     (i)    eddy-visc at integr. point
\param   kapepsint   double	     (i)    kapeps at integr. point
\param  *funct       double	     (i)    nat. shape funcs
\param   fac 	   double	     (i)    weighting factor
\param   factor2 	   double	     (i)    factor
\param   vderxy_12   double	     (i)    factor
\param   visc 	   double	     (i)    viscosity
\param   iel	   int           (i)	num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calgalifkapeps(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,
		      double           eddyint,  
                  double           kapepsint,  
                  double          *funct,   
		      double           fac,     
                  double           factor2,  
                  double           vderxy_12,  
                  double           visc,  
                  int              iel      
                 );  
/*!---------------------------------------------------------------------  
\brief stabilisation part of iteration forces for kapeps dofs

<pre>                                                         he  12/02

In this routine the stabilisation part of the iteration forces for kapeps


           /
 THETA*dt | tau_tu * factor * (kapeps_old)^2 * grad(psi) * u  d_omega
         /


LOW-REYNOLD's MODEL only for epsilon:


              / 
(+) THETA*dt | tau_tu * 2.0*visc*nue_t* (vderxy2_12)^2 *  grad(psi) * u  d_omega
            /  

      
</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *ele      ELEMENT	   (i)    actual element
\param   *eforce   double	   (i/o)  element force vector
\param    kapepsint double	   (i)    kapeps at integr. point
\param   *velint   double	   (i)    vel at integr. point
\param   *velint_dc double	   (i)    vel at integr. point for D.C.
\param    eddyint  double	   (i)    eddy-visc. at integr. point
\param   *funct    double	   (i)    nat. shape funcs
\param  **derxy    double	   (i)    global derivative
\param    fac 	 double	   (i)    weighting factor
\param    factor2  double	   (i)    factor
\param    vderxy_12 double	   (i)    factor
\param    visc     double	   (i)    fluid viscosity
\param    iel	   int	   (i)    num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabifkapeps(
                  FLUID_DYN_CALC  *dynvar, 
                  ELEMENT         *ele,      
	            double          *eforce,  
		      double           kapepsint,  
                  double          *velint,  
                  double          *velint_dc,  
		      double           eddyint,  
                  double          *funct,   
		      double         **derxy,   
		      double           fac,     
                  double           factor2, 
                  double           vderxy_12, 
                  double           visc, 
                  int              iel      
                  );  

/************************************************************************
 | f2_caliterrhs_tu_1.c                                                  |
 ************************************************************************/ 
void f2_calgalifkapome(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,
                  double           kapomeint,  
                  double          *funct,   
		      double           fac,     
                  double           factor2,    
                  int              iel      
                 );  

void f2_calstabifkapome(
                  FLUID_DYN_CALC  *dynvar, 
                  ELEMENT         *ele,      
	            double          *eforce,  
		      double           kapomeint,  
                  double          *velint,  
                  double          *velint_dc,  
                  double          *funct,   
		      double         **derxy,   
		      double           fac,     
                  double           factor2,    
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
| f2_calservice_tu.c                                                    |
************************************************************************/
void f2_calset_tu( 
                  FLUID_DYN_CALC  *dynvar, 
                  FLUID_DATA      *data,     
	            ELEMENT         *ele,     
                  ELEMENT         *elev,
                  double          *kapepsn,    
	            double          *kapepsg,
	            double          *kapepspro,
                  double          *eddyg,
                  double          *eddypro,
	            double          *kappa,    
	            double          *kappan,    
	            double          *epsilon,    
	            double         **evel,   
	            double         **xyze
	           );

void f2_shearstress(
	           ELEMENT    *ele,
                 FLUID_DYN_CALC  *dynvar 
                 );

void f2_kapepsi(
             double  *kapepsint,     
             double  *funct,    
	       double  *kapeps,     
             int      iel       
	     ); 
void f2_eddyi(
             double  *eddyint,     
             double  *funct,    
	       double  *eddy,     
             int      iel       
	     ); 
void f2_kappai_tu(	          
             double     *kappaint,     
             double     *kappanint,     
             double     *eps_proint,     
             double     *funct,    
             double     *kappa,    
             double     *kappan,    
             double     *eps_pro,    
             int         iel       
	     ); 

void f2_C_kappa(	          
             double      kapepsint,     
             double     *epsilon,
             double     *funct, 
             double      visc,    
             double     *C_u,
             int         iel       
	     );
void f2_C_eps(	          
             double      kapepsint,     
             double      kappaint,
             double      visc,    
             double     *C_2,
             int         iel       
	     ); 

void f2_v(	          
             double    **vderxy2,
             double     *vderxy_12 
	   ); 
            
void f2_fac_kappa(
                  double   C_u,     
                  double   eddyint,    
	            double  *factor,     
	            double  *factor1,     
	            double  *factor2,     
	            double  *sig     
	           );
                  
void f2_fac_eps(
                  double   C_2,     
                  double   eps_proint,    
                  double   kappaint,    
                  double   kappanint,    
	            double  *factor,     
	            double  *factor1,     
	            double  *factor2,     
	            double  *sig     
	           );
void f2_production(
	            double  **vderxy,     
	            double  *production    
	           );
                  
void f2_eddyirans(	          
             ELEMENT    *eleke,     
             double     *eddyint,     
             double     *funct,    
             double     *eddy,    
             int         iel       
	     ); 
void f2_vel_dc(
		       FLUID_DYN_CALC  *dynvar,
                   double  *velint,    
                   double  *velint_dc,    
	             double  *kapepsderxy      
	         ); 

void f2_kapepsder(
             double  *kapepsderxy,     
             double **derxy,    
	       double  *kapeps,     
             int      iel       
	       ); 
void f2_kapepsder2(
                   double  *kapepsderxy2,    
                   double **derxy2,    
	             double  *kapepsn,      
	             int      iel        
	           ); 

void f2_estifadd_tu(                  
		   double         **estif,   
		   double         **emass, 
		   double         **tmp,   
		   int              iel,   
		   FLUID_DYN_CALC  *dynvar		   		    
	          ); 
                
/************************************************************************
| f2_calservice_tu_1.c                                                  |
************************************************************************/
void f2_calset_tu_1( 
                    FLUID_DYN_CALC  *dynvar, 
                    FLUID_DATA      *data,     
	              ELEMENT         *ele,     
                    ELEMENT         *elev,
                    double          *kapomen,    
	              double          *kapomeg,
	              double          *kapomepro,
                    double          *eddyg,
                    double          *eddypro,
	              double          *kappan,    
	              double          *omega,    
	              double         **evel,
	              double         **xyze
                   );

void f2_kapomei(
             double  *kapomeint,     
             double  *funct,    
	       double  *kapome,     
             int      iel       
	     );
            
void f2_kapomeder(
             double  *kapomederxy,     
             double **derxy,    
	       double  *kapome,    
             int      iel       
	       ); 

void f2_kapomeder2(
                   double  *kapomederxy2,    
                   double **derxy2,    
	             double  *kapomen,      
	             int      iel        
	           ); 

void f2_xi_kappa(
                  double  *kapomederxy,     
                  double  *omegaderxy,    
	            double   omegaint,     
	            double  *xi     
	           ); 

void f2_fac_kappa_1(
                  double   xi,     
                  double   eddyint,    
                  double   kapomeint,    
                  double   omegaint,    
                  double   visc,    
                  double  *factor,    
                  double  *factor1,    
                  double  *factor2,    
	            double  *sig    
	           ); 

void f2_xi_ome(
                  double  **vderxy,     
	            double   kapomeint,     
	            double  *xi     
	           ); 

void f2_fac_ome(
                  double   xi,     
                  double   ome_proint,    
                  double   kappanint,    
                  double   visc,    
                  double  *factor,    
                  double  *factor1,    
                  double  *factor2,    
	            double  *sig    
	           ); 

void f2_vel_dc_1(
		       FLUID_DYN_CALC  *dynvar,
                   double  *velint,    
                   double  *velint_dc,    
	             double  *kapomederxy      
	         ); 

void f2_kappain(	          
             double     *kappanint,     
             double     *ome_proint,     
             double     *funct,    
             double     *kappan,    
             double     *kapomepro,    
             int         iel       
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
 | f2_calstabmat_tu.c                                                      |
 ************************************************************************/
void f2_calstabkkapeps(			      
                ELEMENT         *ele,    
		    ELEMENT         *elev,    
                FLUID_DYN_CALC  *dynvar,
		    double         **estif,  
		    double           kapepsint, 
		    double          *velint, 
		    double          *velint_dc, 
                double           eddyint, 
                double          *kapepsderxy, 
                double          *funct,  
		    double         **derxy,  
		    double         **derxy2, 
		    double           fac,    
		    double           visc,   
		    double           factor, 
		    double           sig, 
                int              iel    
                 );
void f2_calstabmkapeps(
                    ELEMENT         *ele,     
		        FLUID_DYN_CALC  *dynvar,
		        double         **emass,  
    		        double          *velint, 
    		        double          *velint_dc, 
                    double          *funct,  
		        double         **derxy,  
		        double           fac,    
		        int              iel    
                    );
                 
/************************************************************************
 | f2_calstabmat_tu_1.c                                                  |
 ************************************************************************/
void f2_calstabkkapome(			      
                ELEMENT         *ele,    
		    ELEMENT         *elev,    
                FLUID_DYN_CALC  *dynvar,
		    double         **estif,  
		    double           kapomeint, 
		    double          *velint, 
		    double          *velint_dc, 
                double           eddyint, 
                double          *funct,  
		    double         **derxy,  
		    double         **derxy2, 
		    double           fac,    
		    double           visc,   
		    double           factor,
		    double           sig,
                int              iel    
                   );
                   

void f2_calstabmkapome(
                    ELEMENT         *ele,     
		        FLUID_DYN_CALC  *dynvar,
		        double         **emass,  
    		        double          *velint, 
    		        double          *velint_dc, 
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
 | f2_calstabpar_tu.c                                                      |
 ************************************************************************/
void f2_calstabpar_tu(
	            ELEMENT         *ele,      
		      ELEMENT         *elev,
                  FLUID_DYN_CALC  *dynvar,
		      double           eddyint, 
                  double          *velint, 
                  double          *velint_dc, 
                  double           visc    
                  );

/************************************************************************
 | f2_calstabpar_tu_1.c                                                  |
 ************************************************************************/
void f2_calstabpar_tu_1(
	            ELEMENT         *ele,      
		      ELEMENT         *elev,
                  FLUID_DYN_CALC  *dynvar,
		      double           eddyint, 
                  double          *velint, 
                  double          *velint_dc, 
                  double           visc    
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
 | f2_caltimerhs_tu.c                                                      |
 ************************************************************************/
void f2_calgaltfkapeps(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,    
		      double           kapepsint,  
                  double          *velint,   
		      double           eddyint,
                  double          *funct,    
		      double         **derxy,    
		      double         **vderxy,   
		      double          *kapepsderxy,   
                  double           visc,     
		      double           fac,      
                  double           factor,  
                  double           factor1,  
                  double           factor2,  
                  double           sig,  
                  double           vderxy_12,  
                  double           production,  
                  int              iel       
                  );  
void f2_calstabtfkapeps(
                   FLUID_DYN_CALC  *dynvar,
                   ELEMENT         *ele,      
	            double          *eforce,  
	 	      double           kapepsint,  
       	      double          *velint,  
       	      double          *velint_dc,  
		      double           eddyint, 
                  double         **derxy,   
		      double          *kapepsderxy2,   
                  double         **vderxy,  
		      double          *kapepsderxy,
                  double           fac,     
		      double           visc,    
		      double           factor,
                  double           factor1,
                  double           factor2,
                  double           sig,
                  double           vderxy_12,  
                  double           production,  
                  int              iel      
                  );
/************************************************************************
 | f2_caltimerhs_tu_1.c                                                  |
 ************************************************************************/
void f2_calgaltfkapome(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,    
		      double           kapomeint,  
                  double          *velint,   
		      double           eddyint,
                  double          *funct,    
		      double         **derxy,    
		      double         **vderxy,   
		      double          *kapomederxy,   
                  double           visc,     
		      double           fac,      
                  double           factor,  
                  double           factor1,  
                  double           factor2,  
                  double           sig,  
                  double           production,  
                  int              iel       
                  );  

void f2_calstabtfkapome(
                   FLUID_DYN_CALC  *dynvar,
                   ELEMENT         *ele,      
	            double          *eforce,  
	 	      double           kapomeint,  
       	      double          *velint,  
       	      double          *velint_dc,  
		      double           eddyint, 
                  double         **derxy,   
		      double          *kapomederxy2,   
                  double         **vderxy,  
		      double          *kapomederxy,
                  double           visc,     
                  double           fac,     
                  double           factor,
                  double           factor1,
                  double           factor2,
                  double           sig,
                  double           production,  
                  int              iel      
                  );

/************************************************************************
 | f2_caltimerhspro_tu.c                                                      |
 ************************************************************************/
void f2_calgalprofkapeps(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,    
		      double           eddynint,
                  double          *funct,    
                  double           visc,     
		      double           fac,      
                  double           factor1,  
                  double           production,  
                  int              iel       
                  );  

void f2_calstabprofkapeps(
                   FLUID_DYN_CALC  *dynvar,
                   ELEMENT         *ele,      
	            double          *eforce,  
		      double           eddynint, 
                  double          *funct,    
                  double           visc,     
                  double           fac,     
                  double           factor1,
                  double           production,  
                  double           *velint,  
                  double           *velint_dc,  
                  double          **derxy,  
                  int              iel      
                  );
           
/************************************************************************
 | f2_caltimerhspro_tu_1.c                                               |
 ************************************************************************/
void f2_calgalprofkapome(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,    
		      double           eddynint,
                  double          *funct,    
		      double           fac,      
                  double           factor1,  
                  double           production,  
                  int              iel       
                  );
                  
void f2_calstabprofkapome(
                   FLUID_DYN_CALC  *dynvar,
                   ELEMENT         *ele,      
	            double          *eforce,  
		      double           eddynint, 
                  double          *funct,    
                  double           fac,     
                  double           factor1,
                  double           production,  
                  double           *velint,  
                  double           *velint_dc,  
                  double          **derxy,  
                  int              iel      
                  );
/************************************************************************
 | f2_caltuvisc.c                                                        |
 ************************************************************************/
double f2_calvisc(
	           ELEMENT    *ele,
                 double     **vderxy
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
 | f2_inpele_tu.c                                                          |
 ************************************************************************/
void f2tu_dis(ELEMENT *ele0, ELEMENT *ele1);

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
            ELEMENT     *eleke,             
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
/************************************************************************
 | f2_main_tu.c                                                            |
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
	      int         *hasdirich,
	      int         *hasext,
            CONTAINER   *container
	   );
/*! @} (documentation module close)*/	    
