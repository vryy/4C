/*----------------------------------------------------------------------*
 | shell8                                                 m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _SHELL8
{
/*----------------------------- general variables needed by the element */
struct _ARRAY    a3ref;
struct _ARRAY    thick_node;
double           thick;
double           sdc;
int              nGP[3];
int              nGP_tri;
/*---------------------------------------------- array of forces at GPs */
enum
    {
    s8_xyz,
    s8_rst,
    s8_rst_ortho
    }            forcetyp;   
struct _ARRAY4D  forces;
/*-------------------------------------------- array of internal forces */
struct _ARRAY    intforce;
/*-------------------------------------------- variables needed for eas */
int              eas[5];
int              nhyb; 
struct _ARRAY    alfa;
struct _ARRAY    Dtildinv;
struct _ARRAY    Lt;
struct _ARRAY    Rtilde;


/*-------------------------------------------- variables needed for ans */
int              ans;

} SHELL8;

/*----------------------------------------------------------------------*
 | shell8 data                                            m.gee 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _S8_DATA
{
double        xgpr[3];
double        wgtr[3];

double        xgps[3];
double        wgts[3];

double        xgpt[3];
double        wgtt[3];
} S8_DATA;


/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL SHELL8 ROUTINES                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  s8_a3.c                                              m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8a3ref_extern(double   *funct,
                       double  **deriv,
                       double   *thick,
                       double  **a3ref,
                       ELEMENT  *ele);
void s8averdir(double **dir_list, int numa3, double *a3);
/*----------------------------------------------------------------------*
 |  s8_btdb.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_BtDB(double **estif, double **bop, double **D, int iel,
                int numdf, double weight, double **work);
/*----------------------------------------------------------------------*
 |  s8_eas.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_eas(const int    nhyb,
               const double e1,
               const double e2,
               const int    iel, 
               const int   *eas, 
               double     **P);
void s8_transeas(double      **P, 
                    double      **transP,
                    double      **T,
                    double      **akovr,
                    double      **akonr0,
                    double        detr,
                    double        detr0,
                    int           nhyb);
/*----------------------------------------------------------------------*
 |  s8_eps.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_eps(double *strain,double **gmkovc, double **gmkovr);
/*----------------------------------------------------------------------*
 |  s8_funcderiv.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_funct_deriv(double     *funct, 
                       double    **deriv, 
                       double      r, 
                       double      s,
                       ELEMENT_TYP typ,
                       int         option);
/*----------------------------------------------------------------------*
 |  s8_init.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8init(FIELD *actfield);
/*----------------------------------------------------------------------*
 |  s8_inpele.c                                          m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8inp(ELEMENT *ele);
/*----------------------------------------------------------------------*
 |  s8_intg.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8intg(const ELEMENT   *ele,
               S8_DATA         *data,
               int              option);
/*----------------------------------------------------------------------*
 |  s8_jaco.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8jaco(double    *funct,
               double   **deriv,
               double   **x,
               double   **xjm,
               double    *hte,
               double   **a3ref,
               double     e3,
               int        iel,
               double    *det,
               double    *deta,
               int        init);
/*----------------------------------------------------------------------*
 |  s8_load1.c                                           m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8eleload(ELEMENT  *ele,
                  S8_DATA  *data,
                  MATERIAL *mat,
                  double   *global_vec,
                  int       global_numeq,
                  int       init);
void s8loadGP(ELEMENT    *ele,
                double    **eload,
                double      hhi,
                double      wgt,
                double    **xjm,
                double     *funct,
                double    **deriv,
                int         iel,
                double      xi,
                double      yi,
                double      zi);
/*----------------------------------------------------------------------*
 |  s8_loccoordnode.c                                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
double s8_local_coord_node(int node, int flag, ELEMENT_TYP typ);
/*----------------------------------------------------------------------*
 |  s8_main.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void shell8(      FIELD      *actfield,
                  PARTITION  *actpart,
                  INTRA      *actintra,
                  ELEMENT    *ele,
                  ARRAY      *estif_global,
                  ARRAY      *emass_global,
                  double     *global_vec,
                  int         global_numeq,
                  int         kstep,
                  CALC_ACTION *action);
/*----------------------------------------------------------------------*
 |  s8_mat_linel.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_mat_linel(STVENANT *mat, double **g, double **CC);
void s8_mat_stress1(double *stress, double *strain, double **C);
/*----------------------------------------------------------------------*
 |  s8_mtr.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8mtr(double   **x,
              double   **a3,
              double     e3,
              double   **gkov,
              double   **gkon,
              double   **gmkov,
              double   **gmkon,
              double    *det,
              double    *funct,
              double   **deriv,
              double    *hte,
              int        iel,
              double     condfac,
              char       string[]);
/*----------------------------------------------------------------------*
 |  s8_static_ke.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8static_ke(ELEMENT   *ele,                        
                    S8_DATA   *data,                    
                    MATERIAL  *mat,                     
                    ARRAY     *estif_global,            
                    double    *force,                   
                    int        iforce,                  
                    int        kstep,                   
                    int        init);                   
/*----------------------------------------------------------------------*
 |  s8_static_keug.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8static_keug(ELEMENT   *ele,                         /* the element structure */
                    S8_DATA   *data,                       /* element integration data */
                    MATERIAL  *mat,                        /* the material structure */
                    ARRAY     *estif_global,               /* element stiffness matrix (NOT initialized!) */
                    ARRAY     *emass_global,               /* element mass matrix      (NOT initialized!) */
                    double    *force,                      /* global vector for internal forces (initialized!) */
                    int        iforce,                     /* size of force */
                    int        kstep,                      /* actual step in nonlinear analysis */
                    int        init);                      /* init=1 -> init phase / init=0 -> calc. phase / init=-1 -> uninit phase */
/*----------------------------------------------------------------------*
 |  s8_tmat.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tmat(ELEMENT    *ele,
                MATERIAL   *mat,
                double     *stress,
                double     *strain,
                double    **C,
                double    **gmkovc,
                double    **gmkonc,
                double    **gmkovr,
                double    **gmkonr,
                double    **gkovc,
                double    **gkonc,
                double    **gkovr,
                double    **gkonr,
                double      detc,
                double      detr,
                double      e3,
                int         option);
/*----------------------------------------------------------------------*
 |  s8_tmtr.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tmtr(double   **x,
                double   **a3,
                double     e3,
                double   **gkov,
                double   **gkon,
                double   **gmkov,
                double   **gmkon,
                double    *det,
                double    *funct,
                double   **deriv,
                int        iel,
                double     condfac,
                int        flag);
/*----------------------------------------------------------------------*
 |  s8_tvbo.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tvbo(double      e1,
                double      e2,
                double    **bop,
                double     *funct,
                double    **deriv,
                int         iel,
                int         numdf,
                double    **akov,
                double    **a3kvp);
/*----------------------------------------------------------------------*
 |  s8_tvhe_linear.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tvhe_linear(double **gmkovr,
                       double **gmkovc,
                       double **gmkonr,
                       double **gmkonc,
                       double **gkovr,
                       double **gkovc,
                       double  *detr,
                       double  *detc);
/*----------------------------------------------------------------------*
 |  s8_tvma.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tvma(double **D, double **C, double *stress, double *stress_r,
                double e3, double fact, double condfac);
/*----------------------------------------------------------------------*
 |  s8_tvmr.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tvmr(double   **x,
                double   **a3,
                double   **akov,
                double   **akon,
                double   **amkov,
                double   **amkon,
                double    *det,
                double    *funct,
                double   **deriv,
                int        iel,
                double   **a3kvp,
                int        flag);
/*----------------------------------------------------------------------*
 |  s8_tforce.c                                          m.gee 12/01    |
 *----------------------------------------------------------------------*/
void s8_tforce(double **force,
             int      ngauss,
             double **akov,
             double **akon,
             ELEMENT *ele);
void s8_tettr(double x[3][3], double **a, double b[3][3]);
