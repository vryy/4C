#ifdef D_ALE
/*----------------------------------------------------------------------*
 | ALE element                                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _ALE3
{

int                nGP[3];         /* number of gaussian points in rst direction */
int                jacobi;         /* flag whether to use the Jacobean matrix    */
struct _ELEMENT   *my_fluid;       /* pointer to fluid element associated to me  */

} ALE3;


/*----------------------------------------------------------------------*
 | ALE data                                               m.gee 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _ALE3_DATA
{
double        xgpr[3];         /* natural coordinates r of gaussian points */
double        wgtr[3];         /* weights at natural coordinates r */

double        xgps[3];         /* natural coordinates s of gaussian points */
double        wgts[3];         /* weights at natural coordinates r */

double        xgpt[3];         /* natural coordinates t of gaussian points */
double        wgtt[3];         /* weights at natural coordinates r */
} ALE3_DATA;

/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL ALE ROUTINES                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | ale3_bop.c                                                 mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_bop(double **b, double **deriv, double **xjm, double det, int iel);
/*----------------------------------------------------------------------*
 | ale3_calelm.c                                              mn 06/02  |
 *----------------------------------------------------------------------*/
/*void ale_calelm(FIELD *actfield, SOLVAR *actsolv, PARTITION *actpart, 
                INTRA *actintra, int sysarray1, int sysarray2, CALC_ACTION  *action);*/
/*----------------------------------------------------------------------*
 | ale3_call_stiff.c                                          mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_keku(double **s, double **bs, double **d,
	      double fac, int nd, int neps);
/*----------------------------------------------------------------------*
 | ale3_dirich.c                                              mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_setdirich(FIELD  *actfield, STRUCT_DYNAMIC *sdyn);
void ale_caldirich(ELEMENT *actele, double *fullvec, int dim,
		   ARRAY *estif_global);
/*----------------------------------------------------------------------*
 | ale3_dyn_control.c                                         mn 06/02  |
 *----------------------------------------------------------------------*/
void dyn_ale();
/*----------------------------------------------------------------------*
 | ale3_funct_deriv.c                                         mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_funct_deriv(double *funct, double **deriv, double r, double s,
		     double t, DIS_TYP typ, int option);
/*----------------------------------------------------------------------*
 | ale3_inpele.c                                              mn 06/02  |
 *----------------------------------------------------------------------*/
void ale3inp(ELEMENT *ele);
void fluid_to_ale(const FIELD *fluidfield, const FIELD *alefield);
void find_compatible_ele(const ELEMENT *ele1, const ELEMENT *ele2, int *ierr);
/*----------------------------------------------------------------------*
 | ale3_intg.c                                                mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_intg(const ELEMENT *ele, ALE3_DATA  *data, int option);
/*----------------------------------------------------------------------*
 | ale3_jaco.c                                                mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_jaco(double **deriv, double **xjm, double *det, ELEMENT *ele,
	      int iel);
/*----------------------------------------------------------------------*
 | ale3_main.c                                                mn 06/02  |
 *----------------------------------------------------------------------*/
void ale3(PARTITION *actpart, INTRA *actintra, ELEMENT *ele,
	  ARRAY *estif_global, CALC_ACTION *action);
/*----------------------------------------------------------------------*
 | ale3_mat_linel.c                                           mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_mat_linel(STVENANT *mat, double **d);
/*----------------------------------------------------------------------*
 | ale3_rhs.c                                                 mn 06/02  |
 *----------------------------------------------------------------------*/
/*void ale_rhs(FIELD *actfield, SOLVAR *actsolv, PARTITION *actpart,
             INTRA *actintra, int sysarray1, int sysarray2, double *dirich,
             int global_numeq, int kstep, CALC_ACTION *action);*/
/*----------------------------------------------------------------------*
 | ale3_static_ke.c                                           mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_static_ke(ELEMENT *ele, ALE3_DATA *data, MATERIAL *mat,
		   ARRAY *estif_global, int init);
/*----------------------------------------------------------------------*
 | ale3_hourglass.c                                           mn 06/02  |
 *----------------------------------------------------------------------*/
 void ale3_hourglass(ELEMENT *ele, double **s, double vol); 
#endif
