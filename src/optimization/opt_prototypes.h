/*----------------------------------------------------------------------*
<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0771 - 685-6575
</pre>

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
void optupd(INT init);
/*----------------------------------------------------------------------*
 | control output of optimization data                      al 05/01    |
 *----------------------------------------------------------------------*/
void opt_g_out(OPT_GR_OUT gract);
/*----------------------------------------------------------------------*
 | output of fe-mesh, loads, dirichlet conditions...        al 05/01    |
 *----------------------------------------------------------------------*/
void og_write_mesh(INT nmesh);
/*----------------------------------------------------------------------*
 | output of element density in case of topoopt             al 05/01    |
 *----------------------------------------------------------------------*/
void og_write_eledens(INT numdataw);
/*----------------------------------------------------------------------*/
void og_write_displacements(INT kstep);
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
void optobj(DOUBLE *objective);
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | control program for evaluation of equality constraints               |
 *----------------------------------------------------------------------*/
void opteqc(DOUBLE *constraint,INT init);
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | evaluate values of objective function and constraints                |
 *----------------------------------------------------------------------*/
void func(INT *m, DOUBLE *f, DOUBLE *g);
/*----------------------------------------------------------------------*
 |  routine to control static execution                  m.gee 6/01     |
 *----------------------------------------------------------------------*/
void opt_calsta(CALSTA_EXEC stalact);
/*----------------------------------------------------------------------*
 | variational sensitivity analysis                         al 05/01    |
 *----------------------------------------------------------------------*/
void optvsa(DOUBLE *grdobj, DOUBLE *grdcon,INT init);
/*----------------------------------------------------------------------*
 | variational sensitivity analysis                         al 05/01    |
 *----------------------------------------------------------------------*/
void optsmo(DOUBLE *vvar, INT init);
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | create objective function with eigen frequencies or buckling         |
 *----------------------------------------------------------------------*/
void objeig(DOUBLE *objctval);
/*----------------------------------------------------------------------*
 | prototypes for fortran routines                               al 9/01|
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | control program for eigenvalue analysis                  al 05/01    |
 *----------------------------------------------------------------------*/
void sspace(DOUBLE *A, 
            DOUBLE *ACOP,
            DOUBLE *B,
            INT *MAXA,
            DOUBLE *R,
            DOUBLE *EIGV,
            DOUBLE *TT,
            DOUBLE *W,
            DOUBLE *AR,
            DOUBLE *BR,
            DOUBLE *VEC,
            DOUBLE *D,
            DOUBLE *RTOLV,
            INT *IND,
            INT *NN,
            INT *NNM,
            INT *NWK,
            INT *NWM,
            INT *NC,
            INT *NSTA,
            INT *ISOL,
            INT *NNZ,
            DOUBLE *EIGFOU,
            DOUBLE *RFOU,
            INT *NROOT,
            INT *NITEM,
            INT *NSMAX,
            INT *IFSS,
            INT *ISUB,
            DOUBLE *TOLEIG,
            DOUBLE *TOLJAC,
            DOUBLE *SHIFT,
            INT *ISTLDL,
            INT *INIT,
            INT *IPRINT);
void fortranpow (DOUBLE *V,DOUBLE *R,DOUBLE *RE);
void fsdoc  (DOUBLE *var,DOUBLE *df,DOUBLE *dg,DOUBLE *etai,DOUBLE *etha,
             DOUBLE *xdgo,DOUBLE *resu,DOUBLE *resl,DOUBLE *varup,
             DOUBLE *varlo,INT *numvar,DOUBLE *beta,DOUBLE *accit,
             DOUBLE *delta,INT *iprint );
