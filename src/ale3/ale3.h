/*!----------------------------------------------------------------------
\file
\brief headerfile for 3D ale element, containing structures and prototypes

*----------------------------------------------------------------------*/
#ifdef D_ALE
/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief 3D ale element

<pre>                                                              mn 06/02
This structure contains all specific information for a 3D ale element
</pre>

*----------------------------------------------------------------------*/
typedef struct _ALE3
{

INT                nGP[3];         /*!< number of gaussian points in rst direction */
INT                jacobi;         /*!< flag whether to use the Jacobean matrix    */
struct _ELEMENT   *my_fluid;       /*!< pointer to fluid element associated to me  */

} ALE3;


/*!----------------------------------------------------------------------
\brief 3D ale element data

<pre>                                                              mn 06/02
This structure contains the coordinates and weights for numerical integration
of a 3D ale element
</pre>

*----------------------------------------------------------------------*/
typedef struct _ALE3_DATA
{
DOUBLE        xgpr[4];         /*!< natural coordinates r of gaussian points */
DOUBLE        wgtr[4];         /*!< weights at natural coordinates r */

DOUBLE        xgps[4];         /*!< natural coordinates s of gaussian points */
DOUBLE        wgts[4];         /*!< weights at natural coordinates r */

DOUBLE        xgpt[4];         /*!< natural coordinates t of gaussian points */
DOUBLE        wgtt[4];         /*!< weights at natural coordinates r */
} ALE3_DATA;


/*!----------------------------------------------------------------------
\brief  number of time curves for ale elements 
                                                      
<pre>                                                              mn 08/02
number of time curves for ale elements
</pre>

*----------------------------------------------------------------------*/
#define ALENUMTIMECURVE (5)

/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL ALE ROUTINES                           mn 06/02   |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | ale3_bop.c                                                 mn 06/02  |
 *----------------------------------------------------------------------*/
void ale3_bop(DOUBLE **b, DOUBLE **deriv, DOUBLE **xjm, DOUBLE det, INT iel);
/*----------------------------------------------------------------------*
 | ale_call_stiff.c                                          mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_keku(DOUBLE **s, DOUBLE **bs, DOUBLE **d,
	      DOUBLE fac, INT nd, INT neps);
/*----------------------------------------------------------------------*
 | ale_dirich.c                                              mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_setdirich(FIELD  *actfield, ALE_DYNAMIC *adyn,INT actpos);
void ale_setdirich_increment(FIELD *actfield, ALE_DYNAMIC *adyn);
void ale_caldirich(ELEMENT *actele, DOUBLE *fullvec, INT dim,
		   ARRAY *estif_global);
void ale_caldirich_increment(ELEMENT *actele, DOUBLE *fullvec,
		             INT dim, ARRAY *estif_global, INT place);
void ale_setdirich_increment_fsi(FIELD        *actfield, 
                                 ALE_DYNAMIC  *adyn, 
				 INT           actpos);
INT check_ale_dirich(ELEMENT *actele);
/*----------------------------------------------------------------------*
 | ale_dyn_control.c                                         mn 06/02  |
 *----------------------------------------------------------------------*/
void dyn_ale(void);
void dyn_ale_lin(void);
void dyn_ale_nln(void);
void dyn_ale_2step(void);
void dyn_ale_spring(void);
void dyn_ale_laplace(void);
/*----------------------------------------------------------------------*
 | ale_service.c                                            genk 03/03  |
 *----------------------------------------------------------------------*/
void ale_init(
		    FIELD *actfield
             );
/*----------------------------------------------------------------------*
 | ale3_funct_deriv.c                                         mn 06/02  |
 *----------------------------------------------------------------------*/
void ale3_funct_deriv(DOUBLE *funct, DOUBLE **deriv, DOUBLE r, DOUBLE s,
		     DOUBLE t, DIS_TYP typ, INT option);
/*----------------------------------------------------------------------*
 | ale3_inpele.c                                              mn 06/02  |
 *----------------------------------------------------------------------*/
void ale3inp(ELEMENT *ele);
void fluid_to_ale(const FIELD *fluidfield, const FIELD *alefield);
void find_compatible_ele(const ELEMENT *ele1, const ELEMENT *ele2, INT *ierr);
/*----------------------------------------------------------------------*
 | ale3_intg.c                                                mn 06/02  |
 *----------------------------------------------------------------------*/
void ale3_intg(const ELEMENT *ele, ALE3_DATA  *data);
/*----------------------------------------------------------------------*
 | ale3_jaco.c                                                mn 06/02  |
 *----------------------------------------------------------------------*/
void ale3_jaco(DOUBLE **deriv, DOUBLE **xjm, DOUBLE *det, ELEMENT *ele,
	      INT iel);
/*----------------------------------------------------------------------*
 | ale3_main.c                                                mn 06/02  |
 *----------------------------------------------------------------------*/
void ale3(PARTITION *actpart, INTRA *actintra, ELEMENT *ele,
	  ARRAY *estif_global, CALC_ACTION *action, CONTAINER *container);
/*----------------------------------------------------------------------*
 | ale3_mat_linel.c                                           mn 06/02  |
 *----------------------------------------------------------------------*/
void ale3_mat_linel(STVENANT *mat, DOUBLE **d);
/*----------------------------------------------------------------------*
 | ale3_static_ke.c                                           mn 06/02  |
 *----------------------------------------------------------------------*/
void ale3_static_ke(ELEMENT *ele, ALE3_DATA *data, MATERIAL *mat,
		   ARRAY *estif_global, INT init);
/*----------------------------------------------------------------------*
 | ale3_hourglass.c                                           mn 06/02  |
 *----------------------------------------------------------------------*/
 void ale3_hourglass(ELEMENT *ele, DOUBLE **s, DOUBLE vol); 
/*! @} (documentation module close)*/
#endif
