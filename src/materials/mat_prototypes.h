/*!----------------------------------------------------------------------
\file
\brief headerfile for all 3D formulated material laws containing all the 
       prototypes; sorting [11,22,33,12,23,13]


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_MAT

/*! 
\addtogroup MAT 
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL 3D formulated MATERIAL ROUTINES         sh 04/03  |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  mat_el_iso.c                                              sh 04/03  |
 |                                                                      |
 |  linear elastic isotropic material (St.Venant-Kirchhoff              |
 *----------------------------------------------------------------------*/
void mat_el_iso(DOUBLE   youngs,
                DOUBLE   possionratio,
                DOUBLE **d);
void mat_el_iso_inv(DOUBLE   youngs,
                    DOUBLE   possionratio,
                    DOUBLE   dinv[36]);
/*----------------------------------------------------------------------*
 |  mat_el_orth.c                                             sh 04/03  |
 |                                                                      |
 |  linear elastic orthotropic material                                 |
 *----------------------------------------------------------------------*/
void mat_el_orth(DOUBLE    emod1,
                 DOUBLE    emod2,
                 DOUBLE    emod3,
                 DOUBLE    xnue12,
                 DOUBLE    xnue23,
                 DOUBLE    xnue13,
                 DOUBLE    gmod12,
                 DOUBLE    gmod23,
                 DOUBLE    gmod13,
                 DOUBLE  **d);
/*----------------------------------------------------------------------*
 |  mat_pl_hoff.c                                             sh 04/03  |
 |                                                                      |
 |  anisotropic plasticity model, based on the hoffman yield criterion  |
 |  see Dis. Hoermann !                                                 |
 *----------------------------------------------------------------------*/
void mat_pl_hoff_main(
                 PL_HOFF   *mat,      /*!< material properties          */
                 DOUBLE     stress[6],/*!< vector of stresses [11,22,33,12,23,13]  */
                 DOUBLE   **C,        /*!< constitutive matrix          */
                 INT       *iupd,     /*!< controls update of new stresses to wa */
             /*zusaetzliche Parameter*/                       
                 INT       *yip,      /*!< from WA*/
                 DOUBLE    *dhard,    /*!< from WA*/
                 DOUBLE     strain[6],/*!< actual strains from displacements*/
                 DOUBLE     sig[6],   /*!< stresses from WA*/
                 DOUBLE     eps[6],   /*!< strains from WA*/
                 DOUBLE     dkappa[6],    
                 DOUBLE     gamma[6],
                 DOUBLE     rkappa[9]);    
void mat_pl_hoff_homa(PL_HOFF   *mat,       /*!< material properties          */
                      DOUBLE     P[6][6],   /*!< */
                      DOUBLE     Q[6],      /*!< */
                      DOUBLE     rkappa[9]);/*!< */
void mat_pl_hoff_mapl(DOUBLE     stress[6],   /*!< current stresses          */
                      DOUBLE   **C,           /*!< material matrix to be calculated */
                      DOUBLE     P[6][6],     /*!< coupling matrix P */
                      DOUBLE     Q[6],        /*!< coupling vector Q */
                      DOUBLE     dlam,        /*!< plastic multiplier */
                      DOUBLE    *dhard,       /*!< hardening value from WA*/
                      DOUBLE     gamma[6]);   /*!< from WA*/ 
DOUBLE mat_pl_hoff_yilcr(DOUBLE   stress[6], 
                         DOUBLE   P[6][6],
                         DOUBLE   Q[6],
                         DOUBLE   uniax);
void mat_pl_hoff_radi(PL_HOFF   *mat,
                      DOUBLE     sigtr[6], 
                      DOUBLE   **C,
                      DOUBLE     P[6][6],
                      DOUBLE     Q[6],
                      DOUBLE     rkappal[9],
                      DOUBLE    *dhard,
                      DOUBLE    *dlam,
                      DOUBLE     signew[6],
                      DOUBLE     P_hat[6][6],
                      DOUBLE     Q_hat[6],
                      DOUBLE     gamma[6],
                      DOUBLE     dkappa[6]);
void mat_pl_hoff_hard(DOUBLE     DGDS[6],
                      DOUBLE    *dlam,
                      DOUBLE     rkappal[9],
                      DOUBLE     rkappa[9],
                      DOUBLE     dkappa_hat[6],
                      INT        iter,
                      DOUBLE     RM[54]);
void mat_pl_hoff_serv(DOUBLE     DKDL[6],
                      PL_HOFF   *mat,
                      DOUBLE     rkappa[9],
                      DOUBLE     P[6][6],
                      DOUBLE     Q[6],
                      DOUBLE     sig[6],
                      DOUBLE     eta[9],
                      DOUBLE    *dlam,
                      DOUBLE     DKDS[6][6]);
/*----------------------------------------------------------------------*
 |  mat_pl_mises_lin.c                                        sh 08/02  |
 |                                                                      |
 |  plasticity model based on the von Mises yield criterion             |
 |  -> combined linear iso./kin. hardening law                          |
 *----------------------------------------------------------------------*/
void mat_pl_mises_lin_main(DOUBLE   ym,       
                           DOUBLE   pv,       
                           DOUBLE   sigy,     
                           DOUBLE   hard,     
                           DOUBLE   gf,       
                           DOUBLE   betah,    
                           DOUBLE  *stress,   
                           DOUBLE   strain[6],
                           DOUBLE **d,        
                           INT     *iupd,     
                           INT     *yip,      
                           DOUBLE  *epstn,    
                           DOUBLE   sig[6],   
                           DOUBLE   eps[6],   
                           DOUBLE   qn[6],    
                           DOUBLE   dia);     
void mat_pl_mises_lin_yilcr(DOUBLE  E, 
                            DOUBLE  Eh,
                            DOUBLE  betah,
                            DOUBLE  sigy,
                            DOUBLE  epstn,
                            INT     isoft,
                            DOUBLE  dia,
                            DOUBLE *tau,
                            DOUBLE *ft);
void mat_pl_mises_lin_radi(DOUBLE  e, 
                           DOUBLE  eh,
                           DOUBLE  betah,
                           DOUBLE  sigy,
                           DOUBLE  vnu,
                           DOUBLE  dia,
                           DOUBLE *sigma,
                           DOUBLE *qn,
                           INT     isoft,
                           DOUBLE *epstn,
                           DOUBLE *dlam);
void mat_pl_mises_lin_mapl(DOUBLE   e, 
                           DOUBLE   eh,
                           DOUBLE   sigy,
                           DOUBLE   vnu,
                           DOUBLE   dia,
                           DOUBLE  *tau,   /*Praediktorspannungen*/
                           INT      isoft,
                           DOUBLE  *epstn,
                           DOUBLE   dlam,
                           DOUBLE **d);
/*----------------------------------------------------------------------*
 |  mat_pl_dp_lin.c                                           sh 09/03  |
 |                                                                      |
 |  plasticity model based on the 'Drucker Prager' yield criterion      |
 |  -> combined linear iso./kin. hardening law                          |
 |  for theory see Dis. Menrath                                         |
 *----------------------------------------------------------------------*/
void mat_pl_dp_lin_main(DOUBLE   ym,       
                        DOUBLE   pv,       
                        DOUBLE   sigy,     
                        DOUBLE   eh,       
                        DOUBLE   gf,       
                        DOUBLE   betah,    
                        DOUBLE   phi,      
                        DOUBLE   stress[6],
                        DOUBLE   strain[6],
                        DOUBLE **d,        
                        INT     *iupd,     
                        INT     *yip,      
                        DOUBLE  *epstn,    
                        DOUBLE   sig[6],   
                        DOUBLE   eps[6],   
                        DOUBLE   qn[6],    
                        DOUBLE   dia);      
void mat_pl_dp_lin_radi(DOUBLE  hards, 
                        DOUBLE  betah,
                        DOUBLE *epstn,
                        DOUBLE  G,
                        DOUBLE  K,
                        DOUBLE  alpha,
                        DOUBLE  sigma[6],
                        DOUBLE  qn[6],
                        DOUBLE *dlam,
                        DOUBLE  ft_tr,
                        DOUBLE  n[6]);
void mat_pl_dp_lin_mapl(DOUBLE   hards,
                        DOUBLE   betah,
                        DOUBLE   alpha,
                        DOUBLE   dlam,
                        DOUBLE   e,
                        DOUBLE   vnu,
                        DOUBLE **d,
                        DOUBLE   norm,
                        DOUBLE   n[6]);
void mat_pl_dp_lin_radi_apex(DOUBLE  hards, 
                             DOUBLE  betah,
                             DOUBLE *epstn,
                             DOUBLE  G,
                             DOUBLE  K,
                             DOUBLE  alpha1,
                             DOUBLE  sigma[6],
                             DOUBLE  qn[6],
                             DOUBLE *dlam,
                             DOUBLE  ft1,
                             DOUBLE  ft2,
                             DOUBLE  n[6]);
void mat_pl_dp_lin_mapl_apex(DOUBLE   hards,
                             DOUBLE   betah,
                             DOUBLE   alpha1,
                             DOUBLE   dlam,
                             DOUBLE   e,
                             DOUBLE   vnu,
                             DOUBLE **d,
                             DOUBLE   norm, 
                             DOUBLE   n[6]); 
void mat_pl_dp_lin_preval(DOUBLE   sig[6], 
                          DOUBLE  *norm,
                          DOUBLE   n[6]);
/*----------------------------------------------------------------------*
 |  mat_pl_epc.c                                              sh 10/03  |
 |                                                                      |
 |  elasto plastic concrete material model, containing two              |
 |  'Drucker-Prager' cones (+ one for the apex region == inverted cone) |
 |  and a spherical cap. This model is described in detail in           | 
 |  Dis. Menrath (& Dis. Haufe)                                         |
 *----------------------------------------------------------------------*/
void mat_pl_epc_main(DOUBLE   Ec,        
                     DOUBLE   vc,        
                     DOUBLE   ftm,       
                     DOUBLE   fcm,       
                     DOUBLE   gt,        
                     DOUBLE   gc, 
                     DOUBLE   gamma1,    
                     DOUBLE   gamma2,    
                     DOUBLE   gamma3,    
                     DOUBLE   gamma4,    
                     DOUBLE   dia,       
                     DOUBLE   stress[6], 
                     DOUBLE   strain[6], 
                     DOUBLE **d,         
                     INT     *iupd,      
                     INT     *yip,       
                     DOUBLE  *kappa_t,   
                     DOUBLE  *kappa_c,   
                     DOUBLE   sig[6],    
                     DOUBLE   eps[6]);    
void mat_pl_epc_radi_dp(DOUBLE *kappa_t,
                        DOUBLE *kappa_c,
                        DOUBLE  beta[2],
                        DOUBLE  k_bar[3],
                        DOUBLE  kappa_tu,
                        DOUBLE  kappa_cu,
                        DOUBLE  kappa_e,
                        DOUBLE  gamma3,
                        DOUBLE  ftm,
                        DOUBLE  fcm,
                        DOUBLE  G,
                        DOUBLE  K,
                        DOUBLE  alpha,
                        DOUBLE  sigma[6],
                        DOUBLE  dlam[2],
                        DOUBLE  ft_tr,
                        DOUBLE  n[6],
                        INT     yip);
void mat_pl_epc_mapl_dp(DOUBLE   kappa_t,
                        DOUBLE   kappa_c,
                        DOUBLE   beta[2],
                        DOUBLE   kappa_tu,
                        DOUBLE   kappa_cu,
                        DOUBLE   kappa_e,
                        DOUBLE   gamma3,
                        DOUBLE   ftm,
                        DOUBLE   fcm,
                        DOUBLE   alpha,
                        DOUBLE   dlam,
                        DOUBLE   e,
                        DOUBLE   vnu,
                        DOUBLE **d,
                        DOUBLE   norm,
                        DOUBLE   n[6],
                        INT      yip);
void mat_pl_epc_radi_ms(DOUBLE   *kappa_t, 
                        DOUBLE   *kappa_c, 
                        DOUBLE    beta[2], 
                        DOUBLE    k_bar[3],
                        DOUBLE    kappa_tu,
                        DOUBLE    kappa_cu,
                        DOUBLE    kappa_e, 
                        DOUBLE    gamma3,  
                        DOUBLE    ftm,     
                        DOUBLE    fcm,     
                        DOUBLE    G,       
                        DOUBLE    K,       
                        DOUBLE    alpha[3],
                        DOUBLE    sigma[6],
                        DOUBLE    dlam[2], 
                        DOUBLE    ft1,     
                        DOUBLE    ft2,     
                        DOUBLE    n[6],    
                        INT       flag);    
void mat_pl_epc_mapl_ms(DOUBLE   kappa_t,
                        DOUBLE   kappa_c,
                        DOUBLE   beta[2],
                        DOUBLE   kappa_tu,
                        DOUBLE   kappa_cu,
                        DOUBLE   kappa_e,
                        DOUBLE   gamma3,
                        DOUBLE   ftm,
                        DOUBLE   fcm,
                        DOUBLE   alpha[3],
                        DOUBLE   dlam[2],
                        DOUBLE   e,
                        DOUBLE   vnu,
                        DOUBLE **d,
                        DOUBLE   norm,
                        DOUBLE   n[6],
                        INT      yip);
void mat_pl_epc_radi_cap(DOUBLE *kappa_c,
                         DOUBLE  kappa_cu,
                         DOUBLE  kappa_e,
                         DOUBLE  gamma2,
                         DOUBLE  gamma3,
                         DOUBLE  fcm,
                         DOUBLE  I1_tr,
                         DOUBLE  dev_tr,
                         DOUBLE  G,
                         DOUBLE  K,
                         DOUBLE  alpha2,
                         DOUBLE  sigma[6],
                         DOUBLE *dlam);
void mat_pl_epc_mapl_cap(DOUBLE   kappa_c,
                         DOUBLE   kappa_cu,
                         DOUBLE   kappa_e,
                         DOUBLE   gamma2,
                         DOUBLE   gamma3,
                         DOUBLE   fcm,
                         DOUBLE   alpha2,
                         DOUBLE   dlam,
                         DOUBLE   e,
                         DOUBLE   vnu,
                         DOUBLE   R,
                         DOUBLE   L,
                         DOUBLE   sigma[6],
                         DOUBLE **d);
void mat_pl_epc_dia(DOUBLE    Ec, 
                    DOUBLE    gt, 
                    DOUBLE    gc, 
                    DOUBLE    ftm,
                    DOUBLE    fcm,
                    DOUBLE   *dia);
void mat_pl_epc_hard(DOUBLE    beta[2], 
                     DOUBLE    kappa_tu,
                     DOUBLE    kappa_e, 
                     DOUBLE    kappa_cu,
                     DOUBLE    kappa_c, 
                     DOUBLE    kappa_t,
                     DOUBLE    ftm,
                     DOUBLE    fcm,
                     DOUBLE    gamma3,
                     DOUBLE    k_bar[3],
                     DOUBLE    hard[2],
                     DOUBLE    dhdl[2]);
void mat_pl_epc_preval(DOUBLE    sig[6],
                       DOUBLE   *dev,
                       DOUBLE    n[6],
                       DOUBLE   *I1,
                       DOUBLE    gamma2,  
                       DOUBLE    gamma3, 
                       DOUBLE    ftm,    
                       DOUBLE    fcm,    
                       DOUBLE    beta[2],
                       DOUBLE    alpha[3],
                       DOUBLE    k_bar[3],
                       DOUBLE    kappa_t, 
                       DOUBLE    kappa_tu,
                       DOUBLE    kappa_c, 
                       DOUBLE    kappa_e, 
                       DOUBLE    kappa_cu,
                       DOUBLE   *L,       
                       DOUBLE   *R,       
                       DOUBLE   *I1_0);
/*----------------------------------------------------------------------*
 | prototypes for fortran routines                               al 9/01|
 *----------------------------------------------------------------------*/
void c1inv6(DOUBLE *A,DOUBLE *B,INT *IRC);
void c1ab (DOUBLE *A,DOUBLE *B,DOUBLE *R,INT *NZA,INT *NSA,INT *NSB,
           DOUBLE *XNULL);                                 
/*----------------------------------------------------------------------*/
#endif /*D_MAT*/
/*! @} (documentation module close)*/
