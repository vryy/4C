/*!----------------------------------------------------------------------
\file
\brief headerfile for 2D ale element, containing structures and prototypes

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6127
</pre>

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
DOUBLE             quality;        /*!< my element quality */

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
 | ale2_intg.c                                                 mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_intg(const ELEMENT *ele, ALE2_DATA  *data);
/*----------------------------------------------------------------------*
 | ale2_jaco.c                                                 mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_jaco(DOUBLE **deriv, DOUBLE **xjm, DOUBLE *det, DOUBLE **xyz,
	      INT iel);
void ale2_min_jaco(enum _DIS_TYP distyp, DOUBLE **xyz, DOUBLE *min_detF);
/*----------------------------------------------------------------------*
 | ale2_main.c                                                 mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2(PARTITION *actpart, INTRA *actintra, ELEMENT *ele,
	  ARRAY *estif_global, CALC_ACTION *action, CONTAINER *container);
/*----------------------------------------------------------------------*
 | ale2_mat_linel.c                                            mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_mat_linel(STVENANT *mat, DOUBLE **d);
/*----------------------------------------------------------------------*
 | ale2_static_ke.c                                           mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_static_ke(ELEMENT *ele, ALE2_DATA *data, MATERIAL *mat,
		   ARRAY *estif_global, INT init);
void ale2_static_ke_stiff(ELEMENT *ele,
                          ALE2_DATA      *data,
                          MATERIAL       *mat,
                          ARRAY          *estif_global,
                          INT             init,
		          INT             quality);
void ale2_static_ke_prestress(ELEMENT    *ele,
                              ALE2_DATA  *data,
                              MATERIAL   *mat,
                              ARRAY      *estif_global,
                              INT         init,
		              DOUBLE     *rhs,
       			      INT         total_dim,
		              INT         quality);
void ale2_static_ke_step2(ELEMENT   *ele,
                          ALE2_DATA  *data,
                          MATERIAL  *mat,
                          ARRAY     *estif_global,
                          INT	     init,
                          INT	     quality,
			  DOUBLE    *min_stiff,
			  DOUBLE    *max_stiff,
			  DOUBLE    *min,
			  DOUBLE    *max);
void ale2_static_ke_spring(ELEMENT   *ele,
                           ARRAY     *estif_global,
			   INT        quality,
			   INT        init);
void ale2_static_ke_laplace(ELEMENT     *ele,
                            ALE2_DATA   *data,
                            ARRAY       *estif_global,
                            INT          init,
		            INT          quality);
/*----------------------------------------------------------------------*
 | ale2_hourglass.c                                           mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_hourglass(ELEMENT *ele, DOUBLE **s);
/*----------------------------------------------------------------------*
 | ale2_elm_geometry.c                                        ck 06/03  |
 *----------------------------------------------------------------------*/
DOUBLE ale2_el_area(DOUBLE **xyz);
DOUBLE ale2_area_tria(DOUBLE **xyz, INT i, INT j, INT k);
void edge_geometry(INT      i,
                   INT      j,
         	   DOUBLE **xyz,
	           DOUBLE  *length,
	           DOUBLE  *sin_alpha,
	           DOUBLE  *cos_alpha);
void ale2_torsional(INT      i,
                    INT      j,
		    INT      k,
		    DOUBLE **xyz,
		    DOUBLE **k_torsion,
		    INT      init);
void ale2_tors_spring_quad4(DOUBLE **estif, DOUBLE **xyz, INT init);
void ale2_tors_spring_tri3(DOUBLE **estif, DOUBLE **xyz, INT init);
void ale2_deriv_xy(DOUBLE    **deriv_xy,
                   DOUBLE    **deriv,
                   DOUBLE    **xjm,
                   DOUBLE      det,
                   INT         iel);
/*----------------------------------------------------------------------*
 | ale2_mesh_quality.c                                        ck 06/03  |
 *----------------------------------------------------------------------*/
DOUBLE ale2_corner_angle(DOUBLE **xyz);
DOUBLE ale2_corner_angle_tria(DOUBLE **xyz);
DOUBLE ale2_aspect_ratio(DOUBLE **xyz);
DOUBLE ale2_aspect_ratio_tria(DOUBLE **xyz);
void write_element_quality(ELEMENT  *ele,
                           INT       quality,
			   DOUBLE  **xyz,
			   DOUBLE    min_detF);
void ale_quality(FIELD *field,INT step,
                 INTRA  *actintra, PARTITION    *actpart);
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/

