/*!----------------------------------------------------------------------
\file
\brief fluid_prototypes

------------------------------------------------------------------------*/
/* RULE HOW TO ADD NEW FILES AND FUNCTIONS: 
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE!!!
*/   
/*! 
\addtogroup FLUID
*//*! @{ (documentation module open)*/
/************************************************************************
 | fluid_curvature.c                                                       |
 ************************************************************************/
void fluid_curvature(FIELD        *actfield,    
                     PARTITION    *actpart,     
                     INTRA        *actintra,    
                     CALC_ACTION  *action);      

/************************************************************************
 | fluid_dirich.c                                                       |
 ************************************************************************/
void fluid_initdirich(  FIELD          *actfield, 
                        FLUID_DYNAMIC  *fdyn
		     );
void fluid_setdirich(   FIELD           *actfield, 
                        FLUID_DYNAMIC   *fdyn
	            );
void fluid_caldirich(
                        ELEMENT         *actele,  
		        double          *dforces, 
                        double         **estif,   
		        int             *hasdirich
		    ); 

/************************************************************************
 | fluid_dirich_tu.c                                                     |
 ************************************************************************/
void fluid_setdirich_tu(
                       FIELD  *actfield, 
                       FLUID_DYNAMIC *fdyn, 
                       double   *lower_limit_kappa,
                       double   *lower_limit_eps
                       );

void fluid_caldirich_tu( 
                    FLUID_DYN_CALC  *dynvar, 
                     ELEMENT   *actele,  
		        double    *dforces, 
                    double   **estif,   
		       int       *hasdirich
		       );
                        
/************************************************************************
 | fluid_dirich_tu_1.c                                                   |
 ************************************************************************/
void fluid_setdirich_tu_1(
                       FIELD  *actfield, 
                       FLUID_DYNAMIC *fdyn, 
                       double   *lower_limit_kappa,
                       double   *lower_limit_omega
                       );

void fluid_caldirich_tu_1( 
                         FLUID_DYN_CALC  *dynvar, 
                         ELEMENT   *actele,  
		             double    *dforces, 
                         double   **estif,   
		             int       *hasdirich
		             );     


/************************************************************************
 | fluid_dyn.c                                                          |
 ************************************************************************/
void dyn_fluid(void);

/************************************************************************
 | fluid_freesurface.c                                                  |
 ************************************************************************/
void fluid_createfreesurf(void);
void fluid_freesurf_setdofs(void);
void fluid_modcoor(void); 
/************************************************************************
 | fluid_imp_semimp.c                                                   |
 ************************************************************************/
void fluid_isi(FLUID_DYNAMIC *fdyn);

/************************************************************************
 | fluid_imp_semimp_tu.c                                                |
 ************************************************************************/
void fluid_isi_tu(FLUID_DYNAMIC *fdyn);

/************************************************************************
 | fluid_imp_semimp_tu_1.c                                              |
 ************************************************************************/
void fluid_isi_tu_1(FLUID_DYNAMIC *fdyn);

/************************************************************************
 | fluid_mf.c                                                           |
 ************************************************************************/
void fluid_mf(int mctrl);

/************************************************************************
 | fluid_mfcoupling.c                                                   |
 ************************************************************************/
void fluid_initmfcoupling(
                           FIELD         *fluidfield,
			   FIELD         *alefield		
		         );
			   
/************************************************************************
 | fluid_pr_mcorr.c                                                     |
 ************************************************************************/
void fluid_pm(void);

/************************************************************************
 | fluid_service.c                                                      |
 ************************************************************************/
void fluid_startproc(
                          FLUID_DYNAMIC     *fdyn,
		          int               *nfrastep 
		    );
void fluid_tcons(         
                          FLUID_DYNAMIC     *fdyn,
                          FLUID_DYN_CALC    *dynvar
		);
void fluid_icons(         
                          FLUID_DYNAMIC     *fdyn,
                          FLUID_DYN_CALC    *dynvar,
		          int                itnum           
		);
void fluid_init(
		          FIELD             *actfield,  
                          FLUID_DYNAMIC     *fdyn,
		          int                numr,
		          FLUID_STRESS       str	
	       );


void fluid_norm(          
                          FLUID_DYNAMIC     *fdyn, 	     
                          FIELD             *actfield,    
		          int                numeq_total, 
                          double            *vrat,        
		          double            *prat         
	       );
void fluid_sol_copy(       
                          FIELD             *actfield,
			  int                disnum,
			  int                arrayfrom,
			  int                arrayto,  
                          int                from,     
		          int                to,       
		          int                numdf      
		  );
int fluid_steadycheck(    
                          FLUID_DYNAMIC     *fdyn, 	  
                          FIELD             *actfield,   
		          int                numeq_total 
		     );
int fluid_convcheck(      
                          FLUID_DYNAMIC     *fdyn,   
                          double             vrat,  
		          double             prat,
			  double             grat,  
                          int                itnum, 
		          double             te,    
		          double             ts     
		   );
void fluid_algoout(       
                          FLUID_DYNAMIC     *fdyn, 
                          FLUID_DYN_CALC    *dynvar
		  );		   		     		  
void fluid_reduceshstr(INTRA             *actintra,
                         FIELD             *actfield,
                         FLUID_DYN_CALC    *dynvar);
void fluid_nullshstr(INTRA             *actintra,
                        PARTITION         *actpart,
                        FIELD             *actfield);
void fluid_reducestress(  
                          INTRA             *actintra,
                          FIELD             *actfield,
			  int                numdf, 
			  FLUID_STRESS       str
		       );
	       
/************************************************************************
 | fluid_service_tu.c                                                   |
 ************************************************************************/
void fluid_init_tu(
		        FIELD  *actfield,  
                    FLUID_DYNAMIC *fdyn	
	            );      
void fluid_tcons_tu(FLUID_DYNAMIC *fdyn,
                    FLUID_DYN_CALC *dynvar);

void fluid_set_check_tu(
                       FIELD  *actfield, 
                       FLUID_DYNAMIC *fdyn, 
                       double lower_limit_kappa);

void fluid_eddy_pro(FIELD         *actfield 
                  );

int fluid_convcheck_tu(FLUID_DYNAMIC *fdyn,   
                       double         kapepsrat,  
                       int            itnum1, 
		           double         te,    
		           double         ts     
		          );
int fluid_convcheck_test(
                     FLUID_DYNAMIC *fdyn, 
                     FIELD         *actfield, 
                     int            itnum_check 
		         );

void fluid_icons_tu(FLUID_DYNAMIC *fdyn,
                    FLUID_DYN_CALC *dynvar,
		        int itnum1,
		        int itnumke,
                    int itnum_n
		        );

void fluid_copysol_tu(FLUID_DYNAMIC *fdyn, 
                       FIELD         *actfield,  
                       int            from,     
		           int            to,       
		           int            flag      
		          );
                       
void fluid_copysol_test(FLUID_DYNAMIC *fdyn, 
                   FIELD         *actfield,  
                   int            from,     
                   int            to     
		       );

void fluid_algoout_tu(FLUID_DYNAMIC  *fdyn, 
                      FLUID_DYN_CALC *dynvar
		         );
                   
/************************************************************************
 | fluid_service_tu_1.c                                                 |
 ************************************************************************/
void fluid_set_check_tu_1(
                       FIELD  *actfield, 
                       FLUID_DYNAMIC *fdyn, 
                       double lower_limit_kappa,
                       double lower_limit_omega
                       );

/************************************************************************
 | fluid_stationary.c                                                   |
 ************************************************************************/
void fluid_stat(void);


/************************************************************************
 | inp_fluid_start_data.c                                               |
 ************************************************************************/
void inp_fluid_start_data(void);

/*! @} (documentation module close)*/
