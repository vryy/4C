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
                        FLUID_DYNAMIC   *fdyn,
			INT              pos
	            );
void fluid_setdirich_cyl(FIELD  *actfield, 
                         FLUID_DYNAMIC *fdyn
			 );
void fluid_setdirich_parabolic(FIELD  *actfield, 
                         FLUID_DYNAMIC *fdyn
			 );			 		    
void fluid_setdirich_sd(   FIELD           *actfield, 
                           FLUID_DYNAMIC   *fdyn
	               );
void fluid_caldirich(
                        ELEMENT         *actele,  
		        DOUBLE          *dforces, 
                        DOUBLE         **estif,   
		        INT             *hasdirich,
			INT		 readfrom
		    ); 
void fluid_pm_caldirich(
                     ELEMENT   *actele,  
		     DOUBLE   *dforces, 
                     DOUBLE   **estif,
		     DOUBLE   **emass,
		     DOUBLE   dt,   
		     DOUBLE   theta,
		     INT       *hasdirich
		    );
void fluid_pm_caldirich_parabolic(
                     ELEMENT   *actele,  
		     DOUBLE   *dforces, 
                     DOUBLE   **estif,
		     DOUBLE   **emass,
		     DOUBLE   dt,   
                     DOUBLE   theta,
		     INT       *hasdirich
		    );	
void fluid_pm_caldirich_cyl(
                     ELEMENT   *actele,  
		     DOUBLE   *dforces, 
                     DOUBLE   **estif,
		     DOUBLE   **emass,
		     DOUBLE   dt,   
                     DOUBLE   theta,
		     INT       *hasdirich
		    );		 
		       
/************************************************************************
 | fluid_dirich_tu.c                                                     |
 ************************************************************************/
void fluid_setdirich_tu(
                       FIELD  *actfield, 
                       FLUID_DYNAMIC *fdyn, 
                       DOUBLE   *lower_limit_kappa,
                       DOUBLE   *lower_limit_eps
                       );

void fluid_caldirich_tu( 
                    FLUID_DYN_CALC  *dynvar, 
                     ELEMENT   *actele,  
		        DOUBLE    *dforces, 
                    DOUBLE   **estif,   
		       INT       *hasdirich
		       );
                        
/************************************************************************
 | fluid_dirich_tu_1.c                                                   |
 ************************************************************************/
void fluid_setdirich_tu_1(
                       FIELD  *actfield, 
                       FLUID_DYNAMIC *fdyn, 
                       DOUBLE   *lower_limit_kappa,
                       DOUBLE   *lower_limit_omega
                       );

void fluid_caldirich_tu_1( 
                         FLUID_DYN_CALC  *dynvar, 
                         ELEMENT   *actele,  
		             DOUBLE    *dforces, 
                         DOUBLE   **estif,   
		             INT       *hasdirich
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
void fluid_mf(INT mctrl);

/************************************************************************
 | fluid_mfcoupling.c                                                   |
 ************************************************************************/
void fluid_initmfcoupling(
                           FIELD         *fluidfield,
			   FIELD         *alefield		
		         );

/************************************************************************
 | fluid_mlservice.c                                                   |
 ************************************************************************/
void fluid_ml_init(FIELD         *actfield,  
                   FLUID_DYNAMIC *fdyn);      
void fluid_smcopy(PARTITION       *actpart,      
                  FLUID_DYNAMIC   *fdyn);      
void fluid_prgmr(DOUBLE    **smmat,
                 DOUBLE    **smrhs,
		 INT         numeq,
                 INT         numrhs);      
void fluid_add_smat(DOUBLE  **smat,
		    DOUBLE  **semat,
		    INT       numen,
		    INT      *slme,  
                    DOUBLE    fac);      
void fluid_add_smrhs(DOUBLE  **smrhs,
		     DOUBLE  **smerhs,
                     INT       numrhs,
		     INT       numen,
		     INT      *smlme);	  
void fluid_add_ssrhs(DOUBLE   *ssrhs,
		     DOUBLE   *sserhs,
		     INT       numen,
		     INT      *sslme);	  
void fluid_add_intlhs(DOUBLE   *smidiff,
		      DOUBLE  **smiediff,
		      INT       numen,
		      INT      *smlme);      
void fluid_add_intrhs(DOUBLE   *smirhs,
		      DOUBLE   *smierhs,
		      INT       numen,
		      INT      *smlme);      
void fluid_bubint(DOUBLE  *bubint,     
                  DOUBLE  *smfunct,    
	          DOUBLE **ebub,
		  INT      smiel,     
	          INT      nbub); 
void fluid_pbubint(DOUBLE **pbubint,     
                   DOUBLE  *smfunct,    
	           DOUBLE **epbub,
		   INT      smiel,     
		   INT      iel,     
		   INT      nsd); 
void fluid_bubder(DOUBLE **bubderxy,     
                  DOUBLE **smderxy,    
	          DOUBLE **ebub,
		  INT      smiel,     
	          INT      nbub,
		  INT      nder); 
void fluid_pbubder(DOUBLE ***pbubderxy,     
                   DOUBLE  **smderxy,    
	           DOUBLE  **epbub,
		   INT       smiel,     
		   INT       iel,     
	           INT       nsd,
		   INT       nder); 
void fluid_calnofo(DOUBLE       *erhs,  
		   DOUBLE  	*funct,  
		   DOUBLE  	 fac,	 
		   INT		 iel);
void fluid_mlcaldirich(
                       ELEMENT   *actele,  
		       DOUBLE    *dforces, 
                       DOUBLE   **estif,   
		       INT       *hasdirich
		      );     

/************************************************************************
 | fluid_service.c                                                      |
 ************************************************************************/
void fluid_startproc(
                          FLUID_DYNAMIC     *fdyn,
		          INT               *nfrastep 
		    );
void fluid_cons(         
                          FLUID_DYNAMIC     *fdyn,
                          FLUID_DYN_CALC    *dynvar
		);
void fluid_tcons(         
                          FLUID_DYNAMIC     *fdyn,
                          FLUID_DYN_CALC    *dynvar
		);
void fluid_icons(         
                          FLUID_DYNAMIC     *fdyn,
                          FLUID_DYN_CALC    *dynvar,
		          INT                itnum           
		);
void fluid_init(
                          PARTITION	    *actpart,
                          INTRA	            *actintra,
			  FIELD             *actfield,  
                          FLUID_DYNAMIC     *fdyn,
                          CALC_ACTION       *action,
			  CONTAINER         *container,
		          INT                numr,
		          FLUID_STRESS       str	
	       );
void fluid_norm(          
                          FLUID_DYNAMIC     *fdyn, 	     
                          FIELD             *actfield,    
		          INT                numeq_total, 
                          DOUBLE            *vrat,        
		          DOUBLE            *prat         
	       );
void fluid_sol_copy(       
                          FIELD             *actfield,
			  INT                disnum,
			  INT                arrayfrom,
			  INT                arrayto,  
                          INT                from,     
		          INT                to,       
		          INT                numdf      
		  );
void fluid_transpres(       
                          FIELD             *actfield,
			  INT                disnum,
			  INT                index,
			  INT                actpos,
                          INT                predof,
			  INT                option     
		    );		  
INT fluid_steadycheck(    
                          FLUID_DYNAMIC     *fdyn, 	  
                          FIELD             *actfield,   
		          INT                numeq_total 
		     );
INT fluid_convcheck(      
                          FLUID_DYNAMIC     *fdyn,   
                          DOUBLE             vrat,  
		          DOUBLE             prat,
			  DOUBLE             grat,  
                          INT                itnum, 
		          DOUBLE             te,    
		          DOUBLE             ts     
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
			  INT                numdf, 
			  FLUID_STRESS       str
		       );
void fluid_acceleration(	FIELD 		*actfield, 
				INT 	 	iop, 
				FLUID_DYN_CALC 	*dynvar,
				FLUID_DYNAMIC	*fdyn
			);
void fluid_prep_rhs(FIELD 		*actfield, 
		    FLUID_DYN_CALC 	*dynvar,
		    FLUID_DYNAMIC	*fdyn);
void fluid_predictor(FIELD *actfield, INT iop, FLUID_DYN_CALC *dynvar);
void fluid_lte(	FIELD	 	*actfield, 
		INT 		 iop, 
		FLUID_DYN_CALC 	*dynvar,
		FLUID_DYNAMIC	*fdyn );
void fluid_lte_norm(	
			PARTITION 	*actpart,
			INTRA		*actintra,
			FLUID_DYN_CALC	*dynvar,
                        FLUID_DYNAMIC	*fdyn,
			INT		*iststep,
			INT		*repeat,
			INT		*repeated
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
                       DOUBLE lower_limit_kappa);

void fluid_eddy_pro(FIELD         *actfield 
                  );

INT fluid_convcheck_tu(FLUID_DYNAMIC *fdyn,   
                       DOUBLE         kapepsrat,  
                       INT            itnum1, 
		           DOUBLE         te,    
		           DOUBLE         ts     
		          );
INT fluid_convcheck_test(
                     FLUID_DYNAMIC *fdyn, 
                     FIELD         *actfield, 
                     INT            itnum_check 
		         );

void fluid_icons_tu(FLUID_DYNAMIC *fdyn,
                    FLUID_DYN_CALC *dynvar,
		        INT itnum1,
		        INT itnumke,
                    INT itnum_n
		        );

void fluid_copysol_tu(FLUID_DYNAMIC *fdyn, 
                       FIELD         *actfield,  
                       INT            from,     
		           INT            to,       
		           INT            flag      
		          );
                       
void fluid_copysol_test(FLUID_DYNAMIC *fdyn, 
                   FIELD         *actfield,  
                   INT            from,     
                   INT            to     
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
                       DOUBLE lower_limit_kappa,
                       DOUBLE lower_limit_omega
                       );

/************************************************************************
 | fluid_stationary.c                                                   |
 ************************************************************************/
void fluid_stat(void);


/************************************************************************
 | inp_fluid_start_data.c                                               |
 ************************************************************************/
void inp_fluid_start_data( FIELD   *actfield,
                           FLUID_DYNAMIC *fdyn
			 )
;

/*! @} (documentation module close)*/
