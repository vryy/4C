/*!----------------------------------------------------------------------
\file
\brief headerfile for multilayer shell element, containing structures and
       prototypes

*----------------------------------------------------------------------*/
#ifdef D_SHELL9

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | shell9             modified from shell8              by     sh 10/02 |
 *----------------------------------------------------------------------*/
typedef struct _SHELL9
{
/*----------------------------- general variables needed by the element */
struct _ARRAY    a3ref;
struct _ARRAY    thick_node;
double           thick;
double           sdc;
int              nGP[3];
int              nGP_tri;
int              num_klay;      /* number of kinematic layers */
double          *klayhgt;       /* hgt of a kinematic layer in % of total thickness of the shell */
struct _KINLAY  *kinlay;        /* one kinematic layer -> struct is defined in materials.h */
int              numdf;         /* number of dofs per node of the s9 element */
/*---------------------------------------------- array of forces at GPs */
enum
    {
    s9_xyz,
    s9_rst,
    s9_rst_ortho
    }            forcetyp;    /*-> forces and stresses in same coordinate system*/
struct _ARRAY4D  forces;      /*stress resultants at GPs -> mid surface */
struct _ARRAY4D  stresses;    /*physical stresses at GPs*/
/*struct _ARRAY    energy;*/
/*-------------------------------------------- array of internal forces */
struct _ARRAY    intforce;
/*-------------------------------------------- variables needed for eas */
int              oldkstep;                
int              eas[5];
int              nhyb; 
struct _ARRAY    alfa;
struct _ARRAY    Dtildinv[MAXKLAY_SHELL9];
struct _ARRAY    Lt[MAXKLAY_SHELL9];
struct _ARRAY    Rtilde[MAXKLAY_SHELL9];
/*-------------------------------------------- variables needed for ans */
int              ans;
} SHELL9;


/*----------------------------------------------------------------------*
 | shell9 data                                            m.gee 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _S9_DATA
{
double        xgpr[3];
double        wgtr[3];

double        xgps[3];
double        wgts[3];

double        xgpt[3];
double        wgtt[3];
} S9_DATA;


/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL SHELL9 ROUTINES      (shell8: m.gee 11/01)        |
 |                modified for shell9                      by sh 10/02  |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  s9_a3.c                                              m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9a3(ELEMENT   *ele);
void s9a3ref_extern(double   *funct,
                    double  **deriv,
                    double  **a3ref, 
                    ELEMENT  *ele);
void s9averdir(double **dir_list, int numa3, double *a3, double h2);
/*----------------------------------------------------------------------*
 |  s9_ans.c                                            m.gee 03/02    |
 *----------------------------------------------------------------------*/
void s9_ans_colloqpoints(int       nsansq,     int       iel,
                         int       ans,        DIS_TYP   distyp,
                         double    xr1[],      double    xs1[],
                         double    xr2[],      double    xs2[],
                         double   *funct1q[],  double  **deriv1q[], 
                         double   *funct2q[],  double  **deriv2q[],
                         double  **xrefe,      double ***a3r,
                         double  **xcure,      double ***a3c,
                         double ***akovr1q[],  double ***akonr1q[],
                         double ***amkovr1q[], double ***amkonr1q[],
                         double ***a3kvpr1q[],
                         double ***akovc1q[] , double ***akonc1q[],
                         double ***amkovc1q[], double ***amkonc1q[],
                         double ***a3kvpc1q[],
                         double ***akovr2q[] , double ***akonr2q[],
                         double ***amkovr2q[], double ***amkonr2q[],
                         double ***a3kvpr2q[],
                         double ***akovc2q[] , double ***akonc2q[],
                         double ***amkovc2q[], double ***amkonc2q[],
                         double ***a3kvpc2q[],
                         double  **akovh,
                         double  **akonh,
                         double  **amkovh,
                         double  **amkonh,
                         int       num_klay,   int       numdf);
void s9_ans_colloqcoords(double xqr1[], double xqs1[],
                         double xqr2[], double xqs2[],
                         int    iel,    int    ans);
void s9_ansq_funct(double frq[], double fsq[], double r, double s,
                   int    iel,   int    nsansq);
void s9_ans_bbar_q(double  **bop,        double    frq[],      double fsq[],
                   double   *funct1q[],  double   *funct2q[],
                   double  **deriv1q[],  double  **deriv2q[],
                   double ***akovc1q[],  double ***akovc2q[],
                   double ***a3kvpc1q[], double ***a3kvpc2q[],
                   int       iel,        int       numdf, 
                   int       num_klay,   int       klay,
                   double    condfac,    int       nsansq);
void s9_ans_tvhe_q(double  **gmkovr,     double  **gmkovc,    
                   double  **gmkonr,     double  **gmkonc,
                   double  **gkovr,      double  **gkovc,     
                   double ***amkovc,     double ***amkovr,
                   double ***akovc,      double ***akovr,     
                   double ***a3kvpc,     double ***a3kvpr,
                   double   *detr,       double   *detc,
                   double ***amkovr1q[], double ***amkovc1q[], 
                   double ***akovr1q[] , double ***akovc1q[] ,
                   double ***a3kvpr1q[], double ***a3kvpc1q[],
                   double ***amkovr2q[], double ***amkovc2q[], 
                   double ***akovr2q[] , double ***akovc2q[] ,
                   double ***a3kvpr2q[], double ***a3kvpc2q[],
                   double    frq[],      double    fsq[], 
                   double    e3,         int       nansq, 
                   int       iel,       
                   double    h,           /* total thickness of this element */
                   double   *klayhgt,     /* hight of kin layer in % of total thickness of shell */
                   double   *mlayhgt,     /* hight of mat layer in % of adjacent kin layer */
                   int       num_klay,    /* number of kin layers to this element */
                   int       num_mlay,    /* number of mat layers to this kin layer */
                   int       klay,        /* actual kin layer */
                   int       mlay,        /* actual mat layer of this kin layer */
                   double    condfac);
void s9_ans_tvkg(double **estif,     double  *stress_r,
                 double  *funct,     double **deriv,
                 int      numdf,     int      iel,         
                 double   weight,
                 double   frq[],     double   fsq[],
                 double  *funct1q[], double  *funct2q[],
                 double **deriv1q[], double **deriv2q[],
                 int      ansq,      int       nsansq,
                 int      klay,                       /* actual kin layer */
                 int      num_klay);                  /* number of kin layers to this element */ 
/*----------------------------------------------------------------------*
 |  s9_btdb.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_BtDB(double **estif, double **bop,    double **D,   int iel,
             int      numdf, double   weight, double **work);
/*----------------------------------------------------------------------*
 |  s9_cal_dyn.c                                         m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s9_tmas(double *funct, double *thick, double **emass, int iel, int numdf,
             double facv,   double facw,   double facvw);
/*----------------------------------------------------------------------*
 |  s9_call_mat.c                                           sh 02/03    |
 *----------------------------------------------------------------------*/
void s9_call_mat(ELEMENT    *ele,
                 MULTIMAT   *multimat,    /* material of actual material layer */
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
                 int         rot_axis,
                 double      phi);
void s9_getdensity(MATERIAL   *mat, double *density);
/*----------------------------------------------------------------------*
 |  s9_eas.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_eas(const int    nhyb,
            const double e1,
            const double e2,
            const int    iel, 
            const int   *eas, 
            double     **P);
void s9_transeas(double      **P, 
                 double      **transP,
                 double      **T,
                 double     ***akovr,
                 double      **akonr0,
                 double        detr,
                 double        detr0,
                 int           nhyb,
                 int           klay);
/*----------------------------------------------------------------------*
 |  s9_eps.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_eps(double *strain,double **gmkovc, double **gmkovr);
/*----------------------------------------------------------------------*
 |  s9_funcderiv.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_funct_deriv(double     *funct, 
                    double    **deriv, 
                    double      r, 
                    double      s,
                    DIS_TYP     typ,
                    int         option);
/*----------------------------------------------------------------------*
 |  s9_init.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9init(FIELD *actfield);
/*----------------------------------------------------------------------*
 |  s9_inpele.c                                          m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9inp(ELEMENT *ele);
/*----------------------------------------------------------------------*
 |  s9_intforce.c                                        m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s9_intforce(double *intforce, double *stress_r, double **bop,
                 int iel, int numdf, int nstress_r, double weight);
/*----------------------------------------------------------------------*
 |  s9_intg.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9intg(const ELEMENT   *ele,
            S9_DATA         *data,
            int              option);
/*----------------------------------------------------------------------*
 |  s9_jaco.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9jaco(double    *funct,
            double   **deriv,
            double   **x,
            double   **xjm,
            double    *hte,
            double  ***a3r,
            double     e3,
            int        iel,
            double    *deta,
            int        init,
            int        num_klay,
            double    *klayhgt,     /* hight of kin layer in % of total thickness of shell */
            double    *mlayhgt);    /* hight of mat layer in % of adjacent kin layer */
/*----------------------------------------------------------------------*
 |  s9_load1.c                                           m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9eleload(ELEMENT  *ele,
               S9_DATA  *data,
               double	*loadvec,
               int	 init);
void s9loadGP(ELEMENT    *ele,
              double    **eload1,
              double      wgt,
              double    **xjm,
              double     *funct,
              double    **deriv,
              int         iel,
              double      xi,
              double      yi,
              double      zi,
              double      shift);    /*shell shifter if top,bot */
/*----------------------------------------------------------------------*
 |  s9_load_surf.c                                          sh 12/02    |
 *----------------------------------------------------------------------*/
void s9_surf(double **eload, int nsurf, int numklay, int iel);
void s9_surf_P(double *pload, int nsurf, int numklay);
void s9_surf_onoff(int *pload, int nsurf, int numklay);
/*----------------------------------------------------------------------*
 |  s9_loccoordnode.c                                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
double s9_local_coord_node(int node, int flag, enum _DIS_TYP typ);
/*----------------------------------------------------------------------*
 |  s9_main.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void shell9(FIELD       *actfield,
            PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
            ARRAY       *estif_global,
            ARRAY       *emass_global,
            ARRAY       *intforce_global,
            CALC_ACTION *action,
            CONTAINER   *container);    /* contains variables defined in container.h */
/*----------------------------------------------------------------------*
 |  s9_mat_linel.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_mat_linel(STVENANT *mat, double **g, double **CC);
void s9_mat_stress1(double *stress, double *strain, double **C);
/*----------------------------------------------------------------------*
 |  s9_mat_linel3D.c                                        sh 02/03    |
 *----------------------------------------------------------------------*/
void s9_mat_linel3D(double youngs,
                    double possionratio,
                    double **d);
/*----------------------------------------------------------------------*
 |  s9_mat_orth3D.c                                         sh 02/03    |
 *----------------------------------------------------------------------*/
void s9_mat_orth3D(ORTHOTROPIC *mat, double **d);
/*----------------------------------------------------------------------*
 |  s9_math.c                                               sh 02/03    |
 *----------------------------------------------------------------------*/
void s9_matrix_sort(double **M, char flag1[], char flag2[]);
void s9tcuca (double *str, double **gkov, char flag1[],char flag2[]);
void s9T4sym (double **A, double **C);
void s9shear_cor(double **d, double xsi); 
void s9_Tmaca (double **akov, double phi, int rot_axis, double **C);
void s9_Tcacu (double **akov, double **C);
/*----------------------------------------------------------------------*
 |  s9_out.c                                                sh 12/02    |
 *----------------------------------------------------------------------*/
void s9out_dof(NODE *actnode,int j);
void s9out_nodal_dis(NODE *actnode,int place);
void s9_out_gid_allcoords(FILE *out,int type);
void s9_out_gid_eletop(FILE *out,int type,ELEMENT *actele);
void s9_out_gid_sol_dis(FILE *out,FIELD *actfield,int place);
void s9_out_gid_sol_str(FILE *out,FIELD *actfield,int place);
/*----------------------------------------------------------------------*
 |  s9_restart.c                                         m.gee 05/02    |
 *----------------------------------------------------------------------*/
void s9_write_restart(ELEMENT *actele, int nhandle, long int *handles);
void s9_read_restart(ELEMENT *actele, int nhandle, long int *handles);
/*----------------------------------------------------------------------*
 |  s9_static_keug.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9static_keug(ELEMENT   *ele,                        /* the element structure */
                   S9_DATA   *data,                       /* element integration data */
                   MATERIAL  *mat,                        /* the material structure */
                   ARRAY     *estif_global,               /* element stiffness matrix (NOT initialized!) */
                   ARRAY     *emass_global,               /* element mass matrix (NOT initialized!) */
                   int        kintyp,                     /* typ of kinematic formulation */
                   double    *force,                      /* global vector for internal forces (initialized!) */
                   int        kstep,                      /* actual step in nonlinear analysis */
                   int        init);                      /* init=1 -> init phase / init=0 -> calc. phase / init=-1 -> uninit phase */
/*----------------------------------------------------------------------*
 |  s9_stress.c                                          m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s9_stress(ELEMENT      *ele,
               S9_DATA      *data,
               MATERIAL     *mat,
               int           kintyp,   /* typ of kinematic formulation */
               int           kstep,
               int           init);
void s9_stress_reduce(FIELD     *actfield,
                      PARTITION *actpart,
                      INTRA     *actintra,
                      int        kstep);
/*----------------------------------------------------------------------*
 |  s9_tforce.c                                          m.gee 12/01    |
 *----------------------------------------------------------------------*/
void s9_tforce(double   **force,
               int        ngauss,
               double  ***akov,
               int        klay,    /* actual kin layer */
               ELEMENT   *ele);
void s9_tstress(double   **force,
                double    *stress,
                int        ngauss,
                double   **akov,
                ELEMENT   *ele);
void s9_tettr(double x[3][3], double **a, double b[3][3]);
/*----------------------------------------------------------------------*
 |  s9_tfte.c         m.gee 02/02     modified for shell9   by sh02/03  |
 *----------------------------------------------------------------------*/
void s9_tfte(double  **force,
             int       ngauss,
             double   *stress,
             double  **gkov,
             double ***akon,
             double ***amkov,
             double    h,
             double    e3,
             double   *klayhgt,   /* hight of kin layer in % of total thickness of shell */
             double   *mlayhgt,   /* hight of mat layer in % of adjacent kin layer */
             int       klay,      /* actual kin layer */
             int       mlay,      /* actual mat layer of this kin layer */
             double    fact,
             double    detsm,
             double    detsr);
/*----------------------------------------------------------------------*
 |  s9_tmtr.c                                               sh 10/02    |
 *----------------------------------------------------------------------*/
void s9_tmtr(double   **x,
             double  ***a3,
             double     e3,
             double   **gkov,
             double   **gkon,
             double   **gmkov,
             double   **gmkon,
             double    *det,
             double    *funct,
             double   **deriv,
             int        iel,
             double  ***akov,
             double  ***a3kvp,
             double     h,           /* total thickness of this element */
             double    *klayhgt,     /* hight of kin layer in % of total thickness of shell */
             double    *mlayhgt,     /* hight of mat layer in % of adjacent kin layer */
             int        num_klay,    /* number of kin layers to this element */  
             int        num_mlay,    /* number of mat layers to this kin layer */
             int        klay,        /* actual kin layer */
             int        mlay,        /* actual mat layer of this kin layer */
             double     condfac);
/*----------------------------------------------------------------------*
 |  s9_tvbo.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_tvbo(double      e1,
             double      e2,
             double    **bop,
             double     *funct,
             double    **deriv,
             int         iel,
             int         numdf,
             double   ***akov,
             double   ***a3kvp,
             int         num_klay,    /* number of kin layers to this element */  
             int         klay,        /* actual kin layer */
             double      condfac,
             int         nsansq);
/*----------------------------------------------------------------------*
 |  s9_tvhe.c          m.gee 02/02     modified for shell9     sh 10/02 |
 *----------------------------------------------------------------------*/
void s9_tvhe(double  **gmkovr,
             double  **gmkovc,
             double  **gmkonr,
             double  **gmkonc,
             double  **gkovr,
             double  **gkovc,
             double   *detr,
             double   *detc,
             double ***amkovc,
             double ***amkovr,
             double ***akovc,
             double ***akovr,
             double ***a3kvpc,
             double ***a3kvpr,
             double    e3,
             int       kintyp,    /* typ of kinematic formulation */
             double    h,         /* total thickness of this element */
             double   *klayhgt,   /* hight of kin layer in % of total thickness of shell */
             double   *mlayhgt,   /* hight of mat layer in % of adjacent kin layer */
             int       num_klay,  /* number of kin layers to this element */
             int       num_mlay,  /* number of mat layers to this kin layer */
             int       klay,      /* actual kin layer */
             int       mlay,      /* actual mat layer of this kin layer */
             double    condfac);
void s9_tvhe_lin(double  **gmkovr,
                 double  **gmkovc,
                 double  **gmkonr,
                 double  **gmkonc,
                 double  **gkovr,
                 double  **gkovc,
                 double   *detr,
                 double   *detc);
/*----------------------------------------------------------------------*
 |  s9_tvma.c       m.gee 11/01    modified for shell9   by   sh 10/02  |
 *----------------------------------------------------------------------*/
void s9_tvma(double   **D,
             double   **C,
             double    *stress,
             double    *stress_r,
             double     e3,
             double     fact,
             double     h,           /* total thickness of this element */
             double    *klayhgt,     /* hight of kin layer in % of total thickness of shell */
             double    *mlayhgt,     /* hight of mat layer in % of adjacent kin layer */
             int        num_klay,    /* number of kin layers to this element */  
             int        num_mlay,    /* number of mat layers to this kin layer */
             int        klay,        /* actual kin layer */
             int        mlay,        /* actual mat layer of this kin layer */
             double     condfac);
/*----------------------------------------------------------------------*
 |  s9_tvmr.c                                                sh 10/02   |
 |                                                                      |
 | callculates all the metrics at the middle surface of each kinematic  |
 | layer                                                                |
 *----------------------------------------------------------------------*/
void s9_tvmr(double    **x,
             double   ***a3,
             double   ***akov,
             double   ***akon,
             double   ***amkov,
             double   ***amkon,
             double    **akovh,
             double    **akonh,
             double    **amkovh,
             double    **amkonh,
             double     *det,
             double     *funct,
             double    **deriv,
             int         iel,
             double   ***a3kvp,
             int         num_klay);
/*----------------------------------------------------------------------*
 |  s9_vthv.c          m.gee 02/02     modified for shell9     sh 02/03 |
 *----------------------------------------------------------------------*/
void s9_vthv(double **gmkovc,
             double **gmkonc,
             double  *epsh,
             double  *detc,
             double   e3,
             double   h,         /* total thickness of this element */
             double  *klayhgt,   /* hight of kin layer in % of total thickness of shell */
             double  *mlayhgt,   /* hight of mat layer in % of adjacent kin layer */
             int      num_klay,  /* number of kin layers to this element */
             int      num_mlay,  /* number of mat layers to this kin layer */
             int      klay,      /* actual kin layer */
             int      mlay,      /* actual mat layer of this kin layer */
             double   condfac);
/*----------------------------------------------------------------------*
 |  s9_tvkg.c         m.gee 02/02     modified for shell9   by sh02/03  |
 *----------------------------------------------------------------------*/
void s9_tvkg(double **estif,
             double  *stress_r,
             double  *funct,
             double **deriv,
             int      numdf,
             int      iel,
             double   weight,
             int      klay,        /* actual kin layer */
             int      num_klay);   /* number of kin layers to this element */  
/*----------------------------------------------------------------------*
 |  s9_xint.c                                            m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s9_xint(double *result, double *values, double *funct, int iel);
double s9con(double    e3,       /* zeta_kl */
             int       numlay,   /* number of kinematic layers */
             int       ilay,     /* kin layer where point is located */
             int       jlay,     /* kin layer over which is looped */   
             double    condfac);
double s9notr(int       numlay,   /* number of kinematic layers */
              int       ilay,     /* kin layer where point is located */
              int       jlay);    /* kin layer over which is looped */   
double s9ksi(int       numlay,   /* number of kinematic layers */
             int       ilay,     /* kin layer where point is located */
             int       jlay,     /* kin layer over which is looped */  
             double    condfac); 
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
