/*!----------------------------------------------------------------------
\file
\brief fsi_prototypes

------------------------------------------------------------------------*/
/*! 
\addtogroup FSI
*//*! @{ (documentation module open)*/

/* RULE HOW TO ADD NEW FILES AND FUNCTIONS: 
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE!!!
*/   
/************************************************************************
 | fsi_aitken.c                                                         |
 ************************************************************************/
void fsi_aitken(  
                 FIELD          *structfield, 
                 FSI_DYNAMIC    *fsidyn, 
		 int             itnum
	       );
 
/************************************************************************
 | fsi_ale.c                                                            |
 ************************************************************************/
void fsi_ale(
               FSI_DYNAMIC      *fsidyn,
               STRUCT_DYNAMIC   *sdyn,
               FIELD            *actfield,
               int               mctrl,
               int               numfa
	    );

/************************************************************************
 | fsi_coupling.c                                                       |
 ************************************************************************/
void fsi_createfsicoup(void);
void fsi_initcoupling( 
                          FIELD       *structfield,
                          FIELD       *fluidfield, 
		          FIELD       *alefield
		      );
void fsi_struct_intdofs(
                          FIELD       *structfield, 
			  FSI_DYNAMIC *fsidyn
		       );
		       
/************************************************************************
 | fsi_dyn.c                                                            |
 ************************************************************************/
void dyn_fsi(int mctrl);

/************************************************************************
 | fsi_energy.c                                                          |
 ************************************************************************/
void fsi_dyneint( 
                       FIELD          *structfield, 
                       FSI_DYNAMIC    *fsidyn,
                       STRUCT_DYNAMIC *sdyn,
		       FLUID_DYNAMIC  *fdyn,
		       int             init
		);
void fsi_energycheck(
                        FSI_DYNAMIC   *fsidyn
		    );
		 
/************************************************************************
 | fsi_fluid.c                                                          |
 ************************************************************************/ 
void fsi_fluid(
                       FSI_DYNAMIC    *fsidyn,
		       FLUID_DYNAMIC  *fdyn, 
		       FIELD          *actfield, 
		       int             mctrl,
		       int             numff
	      );

/************************************************************************
 | fsi_relax_intdisp.c                                                  |
 ************************************************************************/
void fsi_relax_intdisp(
                            FIELD          *structfield, 
			    FSI_DYNAMIC    *fsidyn
		      );
 
/************************************************************************
 | fsi_service.c                                                        |
 ************************************************************************/
void fsi_alecp(
		             FIELD           *fluidfield,  
                             FLUID_DYN_CALC  *dynvar,
		             int               numdf,
		             int               phase
	       );
void fsi_aleconv(
		              FIELD           *fluidfield,
		              int              numdf,  
                              int              pos1, 
		              int              pos2
	        );
void fsi_copysol( 
                              FIELD           *actfield,  
                              int              from,	 
		              int              to,	 
		              int              flag	 
		  );
void fsi_fluidstress_result(  
                              FIELD           *actfield, 
                              int              numdf
			   );
void fsi_algoout(             
                              FSI_DYNAMIC      *fsidyn, 
			      int               itnum
	        );
void fsi_structpredictor(
                              FSI_DYNAMIC      *fsidyn, 
                              FIELD            *actfield, 
                              int               init
		         );
int fsi_convcheck(            FIELD            *structfield, 
                              FSI_DYNAMIC      *fsidyn, 
			      int               itnum
		 );
		 			 
/************************************************************************
 | fsi_struct.c                                                         |
 ************************************************************************/
void fsi_struct(   
                   FSI_DYNAMIC       *fsidyn,
                   STRUCT_DYNAMIC    *sdyn, 
		   FIELD             *actfield, 
		   int                mctrl, 
		   int                numfs,
		   int                fsiitnum
	       );

/*! @} (documentation module close)*/
