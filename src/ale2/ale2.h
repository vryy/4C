/*!----------------------------------------------------------------------
\file
\brief headerfile for 2D ale element, containing structures and prototypes

*----------------------------------------------------------------------*/
#ifdef D_ALE

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief 2D ale element

<pre>                                                              mn 06/02
This structure contains all specific information for a 2D ale element
</pre>

*----------------------------------------------------------------------*/
typedef struct _ALE2
{

INT                nGP[2];         /*!< number of gaussian points in rs direction */
INT                jacobi;         /*!< flag whether to use the Jacobean matrix    */
struct _ELEMENT   *my_fluid;       /*!< pointer to fluid element associated to me  */

} ALE2;


/*!----------------------------------------------------------------------
\brief 2D ale element data

<pre>                                                              mn 06/02
This structure contains the coordinates and weights for numerical
integration of a 2D ale element
</pre>

*----------------------------------------------------------------------*/
typedef struct _ALE2_DATA
{
DOUBLE        xgpr[3];         /*!< natural coordinates r of gaussian points */
DOUBLE        wgtr[3];         /*!< weights at natural coordinates r */

DOUBLE        xgps[3];         /*!< natural coordinates s of gaussian points */
DOUBLE        wgts[3];         /*!< weights at natural coordinates s */
} ALE2_DATA;

/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL ALE ROUTINES                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | ale2_bop.c                                                  mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_bop(DOUBLE **b, DOUBLE **deriv, DOUBLE **xjm, DOUBLE det, INT iel);
/*----------------------------------------------------------------------*
 | ale2_funct_deriv.c                                          mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_funct_deriv(DOUBLE *funct, DOUBLE **deriv, DOUBLE r, DOUBLE s,
		     DIS_TYP typ, INT option);
/*----------------------------------------------------------------------*
 | ale2_inpele.c                                               mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2inp(ELEMENT *ele);
/*----------------------------------------------------------------------*
 | ale2_INTg.c                                                 mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_intg(const ELEMENT *ele, ALE2_DATA  *data);
/*----------------------------------------------------------------------*
 | ale2_jaco.c                                                 mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_jaco(DOUBLE **deriv, DOUBLE **xjm, DOUBLE *det, ELEMENT *ele,
	      INT iel);
/*----------------------------------------------------------------------*
 | ale2_main.c                                                 mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2(PARTITION *actpart, INTRA *actintra, ELEMENT *ele,
	  ARRAY *estif_global, CALC_ACTION *action);
/*----------------------------------------------------------------------*
 | ale2_mat_linel.c                                            mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_mat_linel(STVENANT *mat, DOUBLE **d);
/*----------------------------------------------------------------------*
 | ale2_static_ke.c                                           mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_static_ke(ELEMENT *ele, ALE2_DATA *data, MATERIAL *mat,
		   ARRAY *estif_global, INT init);
/*----------------------------------------------------------------------*
 | ale2_hourglass.c                                           mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_hourglass(ELEMENT *ele, DOUBLE **s);
#endif
/*! @} (documentation module close)*/

