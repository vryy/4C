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
void fluid_reducestress(  
                          INTRA             *actintra,
                          FIELD             *actfield,
			  int                numdf, 
			  FLUID_STRESS       str
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
