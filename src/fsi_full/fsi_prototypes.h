/*!----------------------------------------------------------------------
\file
\brief fsi_prototypes

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

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
		 INT             itnum,
                 INT             init
	       );
 
/************************************************************************
 | fsi_ale.c                                                            |
 ************************************************************************/
void fsi_ale(
               FSI_DYNAMIC      *fsidyn,
               ALE_DYNAMIC      *sdyn,
               FIELD            *actfield,
               INT               mctrl,
               INT               numfa
	    );

/************************************************************************
 | fsi_ale_2step.c                                                            |
 ************************************************************************/
void fsi_ale_2step(
                    FSI_DYNAMIC      *fsidyn,
                    ALE_DYNAMIC      *sdyn,
                    FIELD            *actfield,
                    INT               mctrl,
                    INT               numfa
	           );

/************************************************************************
 | fsi_ale_laplace.c                                                            |
 ************************************************************************/
void fsi_ale_laplace(
                     FSI_DYNAMIC      *fsidyn,
                     ALE_DYNAMIC      *sdyn,
                     FIELD            *actfield,
                     INT               mctrl,
                     INT               numfa
	            );

/************************************************************************
 | fsi_ale_lin.c                                                            |
 ************************************************************************/
void fsi_ale_lin(
                  FSI_DYNAMIC      *fsidyn,
                  ALE_DYNAMIC      *sdyn,
                  FIELD            *actfield,
                  INT               mctrl,
                  INT               numfa
	       );

/************************************************************************
 | fsi_ale_nln.c                                                            |
 ************************************************************************/
void fsi_ale_nln(
                   FSI_DYNAMIC      *fsidyn,
                   ALE_DYNAMIC      *sdyn,
                   FIELD            *actfield,
                   INT               mctrl,
                   INT               numfa
	         );

/************************************************************************
 | fsi_ale_spring.c                                                            |
 ************************************************************************/
void fsi_ale_spring(
                    FSI_DYNAMIC      *fsidyn,
                    ALE_DYNAMIC      *sdyn,
                    FIELD            *actfield,
                    INT               mctrl,
                    INT               numfa
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
void dyn_fsi(INT mctrl);

/************************************************************************
 | fsi_energy.c                                                          |
 ************************************************************************/
void fsi_dyneint( 
                       FIELD          *structfield, 
                       FSI_DYNAMIC    *fsidyn,
                       STRUCT_DYNAMIC *sdyn,
		       FLUID_DYNAMIC  *fdyn,
		       INT             init
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
		       INT             mctrl,
		       INT             numff
	      );
/************************************************************************
 | fsi_gradient.c                                                       |
 ************************************************************************/
void fsi_gradient(  
                  FIELD          *alefield, 
                  FIELD          *structfield, 
                  FIELD          *fluidfield, 
                  FSI_DYNAMIC    *fsidyn, 
		  ALE_DYNAMIC    *adyn,
		  FLUID_DYNAMIC  *fdyn,
		  STRUCT_DYNAMIC *sdyn,
		  INT             numfa,
		  INT             numff,
 		  INT             numfs
	         );
void fsi_omega_sg(FIELD          *actfield,
                  FSI_DYNAMIC    *fsidyn  );
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
		             INT               numdf,
		             INT               phase
	       );
void fsi_aleconv(
		              FIELD           *fluidfield,
		              INT              numdf,  
                              INT              pos1, 
		              INT              pos2
	        );
void fsi_copysol( 
                              FIELD           *actfield,  
                              INT              from,	 
		              INT              to,	 
		              INT              flag	 
		  );
void fsi_fluidstress_result(  
                              FIELD           *actfield, 
                              INT              numdf
			   );
void fsi_algoout(             
                              FSI_DYNAMIC      *fsidyn, 
			      INT               itnum
	        );
void fsi_structpredictor(
                              FSI_DYNAMIC      *fsidyn, 
                              FIELD            *actfield, 
                              INT               init
		         );
INT fsi_convcheck(            FIELD            *structfield, 
                              FSI_DYNAMIC      *fsidyn, 
			      INT               itnum
		 );
void fsi_init_ale(FIELD *actfield,INT numr);

/************************************************************************
 | fsi_struct.c                                                         |
 ************************************************************************/
void fsi_struct(   
                   FSI_DYNAMIC       *fsidyn,
                   STRUCT_DYNAMIC    *sdyn, 
		   FIELD             *actfield, 
		   INT                mctrl, 
		   INT                numfs,
		   INT                fsiitnum
	       );

/*! @} (documentation module close)*/
