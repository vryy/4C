/*----------------------------------------------------------------------*
 |  routine to control linear static structural analysis       al 6/01  |
 *----------------------------------------------------------------------*/
void opt_stalin(CALSTA_EXEC stalact);
/*----------------------------------------------------------------------*
 |  routine to control nonlinear static execution          al 11/01     |
 *----------------------------------------------------------------------*/
void opt_stanln(CALSTA_EXEC stalact) ;
/*----------------------------------------------------------------------*
 |  routine to control dynamic eigenvalue analysis              al 08/02|
 *----------------------------------------------------------------------*/
void calfrq(INT init); 
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | control program for update of fem-arrays for act. variable           |
 *----------------------------------------------------------------------*/
void optupd(int init);
/*----------------------------------------------------------------------*
 | control output of optimization data                      al 05/01    |
 *----------------------------------------------------------------------*/
void opt_g_out(OPT_GR_OUT gract);
/*----------------------------------------------------------------------*
 | output of fe-mesh, loads, dirichlet conditions...        al 05/01    |
 *----------------------------------------------------------------------*/
void og_write_mesh(int nmesh);
/*----------------------------------------------------------------------*
 | output of element density in case of topoopt             al 05/01    |
 *----------------------------------------------------------------------*/
void og_write_eledens(int numdataw);
/*----------------------------------------------------------------------*/
void og_write_displacements(int kstep);
/*----------------------------------------------------------------------*
 | initialize execution stage of optimization           a.lipka 5/01    |
 *----------------------------------------------------------------------*/
void opcini(void);
/*----------------------------------------------------------------------*
 | main control routine for oc-methods in fsd               al 05/01    |
 *----------------------------------------------------------------------*/
void optfsd(void);
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | direct update of optimization variables                              |
 *----------------------------------------------------------------------*/
void updvar(void);
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | control program for calculation of objective functions               |
 *----------------------------------------------------------------------*/
void optobj(double *objective);
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | control program for evaluation of equality constraints               |
 *----------------------------------------------------------------------*/
void opteqc(double *constraint,int init);
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | evaluate values of objective function and constraints                |
 *----------------------------------------------------------------------*/
void func(int *m, double *f, double *g);
/*----------------------------------------------------------------------*
 |  routine to control static execution                  m.gee 6/01     |
 *----------------------------------------------------------------------*/
void opt_calsta(CALSTA_EXEC stalact);
/*----------------------------------------------------------------------*
 | variational sensitivity analysis                         al 05/01    |
 *----------------------------------------------------------------------*/
void optvsa(double *grdobj, double *grdcon,int init);
/*----------------------------------------------------------------------*
 | variational sensitivity analysis                         al 05/01    |
 *----------------------------------------------------------------------*/
void optsmo(double *vvar, int init);
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | create objective function with eigen frequencies or buckling         |
 *----------------------------------------------------------------------*/
void objeig(double *objctval);
/*----------------------------------------------------------------------*
 | prototypes for fortran routines                               al 9/01|
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | control program for eigenvalue analysis                  al 05/01    |
 *----------------------------------------------------------------------*/
void sspace(double *A, 
            double *ACOP,
            double *B,
            int *MAXA,
            double *R,
            double *EIGV,
            double *TT,
            double *W,
            double *AR,
            double *BR,
            double *VEC,
            double *D,
            double *RTOLV,
            int *IND,
            int *NN,
            int *NNM,
            int *NWK,
            int *NWM,
            int *NC,
            int *NSTA,
            int *ISOL,
            int *NNZ,
            double *EIGFOU,
            double *RFOU,
            int *NROOT,
            int *NITEM,
            int *NSMAX,
            int *IFSS,
            int *ISUB,
            double *TOLEIG,
            double *TOLJAC,
            double *SHIFT,
            int *ISTLDL,
            int *INIT,
            int *IPRINT);
void fortranpow (double *V,double *R,double *RE);
void fsdoc  (double *var,double *df,double *dg,double *etai,double *etha,
             double *xdgo,double *resu,double *resl,double *varup,
             double *varlo,int *numvar,double *beta,double *accit,
             double *delta,int *iprint );
