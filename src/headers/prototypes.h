/*----------------------------------------------------------------------*
 |  main_ccarat.c                                        m.gee 11/01    |
 *----------------------------------------------------------------------*/
void main(int argc, char *argv[]);
/*----------------------------------------------------------------------*
 |  global_ass_dof.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void assign_dof();
/*----------------------------------------------------------------------*
 |  global_cal_control.c                                 m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntacal();
/*----------------------------------------------------------------------*
 |  cal_dyn_control.c                                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
void caldyn();
/*----------------------------------------------------------------------*
 | cal_static_control.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calsta();
void stalin();
/*----------------------------------------------------------------------*
 | global_control.c                                      m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntam(int argc, char *argv[]);
/*----------------------------------------------------------------------*
 | global_init_control.c                                 m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntaini(int argc, char *argv[]);
/*----------------------------------------------------------------------*
 | global_inp_control.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntainp();
/*----------------------------------------------------------------------*
 | global_mask_matrices.c                                   m.gee 11/01    |
 *----------------------------------------------------------------------*/
void mask_global_matrices();
/*----------------------------------------------------------------------*
 |  machine_hpux.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntadev(int argc, char *argv[]);
/*----------------------------------------------------------------------*
 |  map_node_find.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void iscouple_find_node_comp(NODE  *actnode, 
                                FIELD *searchfield, 
                                NODE **partnernode,
                                int    coupleID,
                                int    dof);
void cheque_distance(double *x1, double *x2, double tol, int *ierr);
void find_assign_coupset(FIELD *actfield, 
                            int    coupleID, 
                            int   *counter);
/*----------------------------------------------------------------------*
 |  out_gid_sol.c                                        m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_sol_init();
void out_gid_domains(FIELD *actfield);
void out_gid_sol(char string[], FIELD *actfield, INTRA  *actintra, int step,
                 int place);
/*----------------------------------------------------------------------*
 |  out_gid_msh.c                                        m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_msh();
void out_gid_allcoords(FILE *out);
/*----------------------------------------------------------------------*
 |  input_cond_couple.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inp_couple(FIELD *field);
/*----------------------------------------------------------------------*
 |  input_cond_fluid.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inp_cond_nodal_fluid(FIELD *field);
void inp_cond_nodal_ale(FIELD *field);
/*----------------------------------------------------------------------*
 |  input_cond_struct.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inp_cond_nodal_struct(FIELD *field);
void inp_cond_ele_struct(FIELD *field);
/*----------------------------------------------------------------------*
 |  input_conditions.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inp_conditions();
/*----------------------------------------------------------------------*
 |  input_control_global.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inpctr();
void inpctrprob();
void inpctrdyn();
void inpctrstat();
/*----------------------------------------------------------------------*
 |  input_ctr_head.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inpctrhed();
void inptrace();
/*----------------------------------------------------------------------*
 |  input_curves.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inp_cond_curve();
void inp_read_curve(char *string);
/*----------------------------------------------------------------------*
 |  input_design.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inpdesign();
void inp_dnode();
void read_1_dnode(DNODE *dnode, int readId);
void inp_dline();
void read_1_dline(DLINE *dline, int readId);
void inp_dsurface();
void read_1_dsurf(DSURF *dsurf, int readId);
void inp_dvolume();
void read_1_dvol(DVOL *dvol, int readId);
/*----------------------------------------------------------------------*
 |  input_design_top.c                                  m.gee 11/01    |
 *---------------------------------------------------------------------*/
void inpdesign_topology_design();
void inpdesign_topology_fe();
/*----------------------------------------------------------------------*
 |  input_material.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inp_material();
/*----------------------------------------------------------------------*
 |  input_mesh.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inpfield();
void inp_assign_nodes(FIELD *field);
void inpnodes();
void inp_struct_field(FIELD *structfield);
void inp_fluid_field(FIELD *fluidfield);
void inp_ale_field(FIELD *alefield);
/*----------------------------------------------------------------------*
 |  input_topology.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inp_topology(FIELD *field);
/*----------------------------------------------------------------------*
 |  math1.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void math_array_copy(double **from, int n, int m, double **to);
void math_inv3(double **a, double *det);
void math_tran(double **a, int n);
void math_unvc(double *enorm,double *vec, int n);
void math_matvecdense(double  *r,
                         double **A,
                         double  *b,
                         int      ni,
                         int      nk,
                         int      init,
                         double   factor);
void math_mattrnvecdense(double  *r,
                         double **A,
                         double  *b,
                         int      ni,
                         int      nk,
                         int      init,
                         double   factor);
void math_matmatdense(double **R,
                         double **A,
                         double **B,
                         int      ni,
                         int      nk,
                         int      nj,
                         int      init,
                         double   factor);
void math_mattrnmatdense(double **R,
                            double **A,
                            double **B,
                            int      ni,
                            int      nk,
                            int      nj,
                            int      init,
                            double   factor);
void math_matmattrndense(double **R,
                            double **A,
                            double **B,
                            int      ni,
                            int      nk,
                            int      nj,
                            int      init,
                            double   factor);
void math_sym_inv(double **A, int dim);
void math_sppr(double *spat, double *a, double *b, double *c);
/*----------------------------------------------------------------------*
 |  sort_find.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
int cmp_int(const void *a, const void *b );
double cmp_double(const void *a, const void *b );
void mg_sort(int list[], int N, int list2[], double list3[]);
int quick_find(int key, int list[], int length, int shift, int bins[]);
void init_quick_find(int list[], int length, int *shift, int *bins);
int find_index(int key, int list[], int length);
/*----------------------------------------------------------------------*
 |  par_assignmesh.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void part_assignfield();
/*----------------------------------------------------------------------*
 |  par_initmetis.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void part_fields();
/*----------------------------------------------------------------------*
 |  par_make_comm.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void create_communicators();
/*----------------------------------------------------------------------*
 |  pss_am.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void* amdef(char *namstr,ARRAY *a,int fdim, int sdim, char typstr[]);
void* amredef(ARRAY *a,int newfdim, int newsdim, char newtypstr[]);
void  amdel(ARRAY *array);
void  amzero(ARRAY *array);
void  aminit(ARRAY *array, void *value);
void* am_alloc_copy(ARRAY *array_from, ARRAY *array_to);
void* amcopy(ARRAY *array_from, ARRAY *array_to);
void* am4def(char    *namstr,
             ARRAY4D *a,
             int      fdim, 
             int      sdim,
             int      tdim,
             int      fodim, 
             char     typstr[]);
void  am4del(ARRAY4D *array);
void  am4zero(ARRAY4D *array);
void  am4init(ARRAY4D *array, void *value);
void* am4_alloc_copy(ARRAY4D *array_from, ARRAY4D *array_to);
void* am4copy(ARRAY4D *array_from, ARRAY4D *array_to);
void* am4redef(ARRAY4D *array, 
               int newfdim, 
               int newsdim, 
               int newtdim,
               int newfodim);
/*----------------------------------------------------------------------*
 |  pss_ds.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void dserror(char string[]);
void dsinit();
void dstrc_enter(char string[]);
void dstrc_exit();
void dstracesize();
void dstracereport(ARRAY *array);
void dstrace_to_err();
/*----------------------------------------------------------------------*
 |  pss_fr.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void frinit();
void frrewind();
void frfind(char string[]);
void frread();
void frint_n(char string[],int *var,int num, int *ierr);
void frint(char string[],int *var, int *ierr);
void frdouble_n(char string[],double *var,int num, int *ierr);
void frdouble(char string[],double *var, int *ierr);
void frchar(char string[],char *var, int *ierr);
void frchk(char string[], int *ierr);
void frend();
/*----------------------------------------------------------------------*
 |  pss_pss.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void pss_write(char          *name, 
                  int            fdim, 
                  int            sdim,
                  int            byte,
                  const void    *startaddress,
                  int           *handle, 
                  int           *ierr);
void pss_write_array(const ARRAY *array, 
                        int         *handle, 
                        int         *ierr);
void pss_read_name(char       *name, 
                      int        *fdim, 
                      int        *sdim,
                      int        *byte,
                      void       *ziel,
                      int        *handle, 
                      int        *ierr);
void pss_read_name_handle(char       *name, 
                             int        *fdim, 
                             int        *sdim,
                             int        *byte,
                             void       *ziel, 
                             int        *handle, 
                             int        *ierr);
void pss_read_array_name(char       *name, 
                            ARRAY      *array,
                            int        *handle,
                            int        *ierr);
void pss_read_array_name_handle(char       *name, 
                                   ARRAY      *array,
                                   int        *handle,
                                   int        *ierr);
void pss_read_array_handle(ARRAY      *array,
                              int        *handle,
                              int        *ierr);
void pss_chck(char       *name,
                 int        *handle, 
                 int        *ierr);
void pss_chck_handle(char       *name,
                        int        *handle, 
                        int        *ierr);
void pss_getdims_name(char       *name, 
                         int        *fdim,
                         int        *sdim,
                         int        *byte,
                         int        *handle,
                         int        *ierr);
void pss_getdims_name_handle(char       *name, 
                                int        *fdim,
                                int        *sdim,
                                int        *byte,
                                int        *handle,
                                int        *ierr);
void pss_status_to_err();
