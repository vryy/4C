#ifdef D_ALE
/*----------------------------------------------------------------------*
 | ALE element                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _ALE2
{

int                nGP[2];         /* number of gaussian points in rs direction */
int                jacobi;         /* flag whether to use the Jacobean matrix    */
struct _ELEMENT   *my_fluid;       /* pointer to fluid element associated to me  */

} ALE2;


/*----------------------------------------------------------------------*
 | ALE data                                               m.gee 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _ALE2_DATA
{
double        xgpr[3];         /* natural coordinates r of gaussian points */
double        wgtr[3];         /* weights at natural coordinates r */

double        xgps[3];         /* natural coordinates s of gaussian points */
double        wgts[3];         /* weights at natural coordinates s */
} ALE2_DATA;

/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL ALE ROUTINES                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | ale2_bop.c                                                  mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_bop(double **b, double **deriv, double **xjm, double det, int iel);
/*----------------------------------------------------------------------*
 | ale2_funct_deriv.c                                          mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_funct_deriv(double *funct, double **deriv, double r, double s,
		     DIS_TYP typ, int option);
/*----------------------------------------------------------------------*
 | ale2_inpele.c                                               mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2inp(ELEMENT *ele);
/*----------------------------------------------------------------------*
 | ale2_intg.c                                                 mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_intg(const ELEMENT *ele, ALE2_DATA  *data, int option);
/*----------------------------------------------------------------------*
 | ale2_jaco.c                                                 mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_jaco(double **deriv, double **xjm, double *det, ELEMENT *ele,
	      int iel);
/*----------------------------------------------------------------------*
 | ale2_main.c                                                 mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2(PARTITION *actpart, INTRA *actintra, ELEMENT *ele,
	  ARRAY *estif_global, CALC_ACTION *action);
/*----------------------------------------------------------------------*
 | ale2_mat_linel.c                                            mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_mat_linel(STVENANT *mat, double **d);
/*----------------------------------------------------------------------*
 | ale2_static_ke.c                                           mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_static_ke(ELEMENT *ele, ALE2_DATA *data, MATERIAL *mat,
		   ARRAY *estif_global, int init);
/*----------------------------------------------------------------------*
 | ale2_hourglass.c                                           mn 06/02  |
 *----------------------------------------------------------------------*/
void ale2_hourglass(ELEMENT *ele, double **s);
#endif
