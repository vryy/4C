/*!----------------------------------------------------------------------
\file
\brief ssi_prototypes

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup SSI
*//*! @{ (documentation module open)*/

/* RULE HOW TO ADD NEW FILES AND FUNCTIONS: 
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE!!!
*/   
/************************************************************************
 | ssi_aitken.c                                                         |
 ************************************************************************/ 
void ssi_aitken(  
                 FIELD          *masterfield, 
                 SSI_DYNAMIC    *ssidyn, 
		 INT             itnum,
                 INT             init
	       ); 
/************************************************************************
 | ssi_coupling.c                                                       |
 ************************************************************************/
void ssi_createssicoup(void);
void ssi_initcoupling( 
                          FIELD       *structfield,
                          FIELD       *fluidfield, 
                          INTERFACES  *interfaces
		      );
void ssi_master_intdofs(
                          FIELD       *masterfield, 
			  SSI_DYNAMIC *ssidyn
		       );

/************************************************************************
 | ssi_dyn.c                                                            |
 ************************************************************************/
void dyn_ssi(void);

/************************************************************************
 | ssi_service.c                                                         |
 ************************************************************************/
void ssiserv_put_disp2slave(
                             FIELD *actfield,
                             CONTAINER *container,
                             DOUBLE *rldfac,
                             SSI_DYNAMIC *ssidyn 
                           );
void ssiserv_put_true_disp2slave(
                                 FIELD *actfield,
                                 CONTAINER *container
                                );

void calc_ssi_coupforce_mod( ELEMENT *actele, 
                             ARRAY *estif_global, 
                             ARRAY *emass_global, 
                             CONTAINER *container
                       );
void init_sol_a_tozero(      FIELD *actfield
                       );
void init_sol_mf_a_tozero(   FIELD *actfield
                         );
void ssiserv_rhs_point_neum( DOUBLE *rhs, 
                             INT dimrhs, 
                             PARTITION *actpart
                           );     
INT ssi_convcheck(           FIELD *structfield, 
                             SSI_DYNAMIC *ssidyn, 
                             INT itnum
                 );
void ssiserv_update_coupforc(FIELD *actfield, 
                             CONTAINER *container
                            );
void ssiserv_erase_coupforc( FIELD *actfield, 
                             CONTAINER *container
                           );
void ssi_detect_nmb_nods_int(FIELD *masterfield, 
                             FIELD *slavefield,
                             INTERFACES *int_faces
                            );
void ssi_mortar_coeff(       FIELD *masterfield, 
                             FIELD *slavefield, 
                             INT step,
                             SSI_DYNAMIC *ssidyn, 
                             INTERFACES *int_faces
                     );
void ssi_detect_intersection(DOUBLE lambdr1_2, 
                             DOUBLE lambdr2_1, 
                             DOUBLE lambdr2_2, 
                             DOUBLE nr1_2, 
                             DOUBLE nr2_1, 
                             DOUBLE nr2_2, 
                             DOUBLE *b1, 
                             DOUBLE *b2, 
                             INT *intersection
                            );
void ssi_put_coupforc2mst(   FIELD *masterfield, 
                             FIELD *slavefield, 
                             INTERFACES *int_faces
                         );
void ssi_interpolation_disp( FIELD *field_to, 
                             FIELD *field_from, 
                             INT place,
                             SSI_DYNAMIC *ssidyn
                      );
void ssi_interpolation_frc(  FIELD *field_to, 
                             FIELD *field_from, 
                             INT place,
                             SSI_DYNAMIC *ssidyn
                      );
void ssi_calc_disp4slv(      FIELD *masterfield, 
                             FIELD *slavefield, 
                             SSI_DYNAMIC *ssidyn,
                             INTERFACES *int_faces 
                      );
void ssi_dbg_plot_frcs(      FIELD *masterfield, 
                             FIELD *slavefield, 
                             DOUBLE time, 
                             INT itnum, 
                             SSI_DYNAMIC *ssidyn,
                             INTERFACES *int_faces
                      );
void ssi_dbg_plot_frc(       FIELD *masterfield, 
                             FIELD *slavefield, 
                             DOUBLE time, 
                             INT itnum,
                             DOUBLE step, 
                             DOUBLE nstep, 
                             INTERFACES *int_faces
                      );
void ssi_dbg_plot_disp(      FIELD *masterfield, 
                             FIELD *slavefield, 
                             DOUBLE step, 
                             DOUBLE nstep, 
                             INT time, 
                             INT maxtime,
                             INTERFACES *int_faces
                      );
void ssi_alloc_dirich(       FIELD *slavefield
                     );
void ssi_dbg_plot_disps(     FIELD *masterfield, 
                             FIELD *slavefield, 
                             DOUBLE step, 
                             INT itnum, 
                             SSI_DYNAMIC *ssidyn,
                             INTERFACES *int_faces
                       );
/************************************************************************
 | ssi_struct.c                                                         |
 ************************************************************************/
void ssi_struct(   
                   SSI_DYNAMIC       *ssidyn,
                   STRUCT_DYNAMIC    *sdyn, 
		   FIELD             *actfield, 
		   INT                mctrl, 
		   INT                ssiitnum,
		   enum _SSI_COUPTYP  ssi_couptyp
	       );
void ssi_struct_decide(INT mctrl, INT numaf); 

/*! @} (documentation module close)*/
