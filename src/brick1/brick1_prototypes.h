/*!----------------------------------------------------------------------
\file
\brief headerfile for 3D hex element, containing prototypes

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | evaluates linear/nonlinear strains from               al    9/01     |
 | displacement derivatives                                             |
 *----------------------------------------------------------------------*/
void c1_eps( double   *disd,                /* displacement derivatives */
             double   *eps,                 /* strain vector            */
             int      iform);        /* index for nonlinear formulation */
/*----------------------------------------------------------------------*
 | evaluates element forces                              al    9/01     |
 *----------------------------------------------------------------------*/
void c1fi( double  *F,   /*  force vector integral (stress-resultants)  */
           double   fac, /*  multiplier for numerical integration       */
           double **bop, /*  b-operator matrix                          */
           int      nd,  /*  total number degrees of freedom of element */
           double  *fie);/*  internal force vector                      */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   radial return for elements with mises material model               |
 |                                                                      |
 *----------------------------------------------------------------------*/
void c1rad1(double e,        /* young's modulus                         */
            double eh,       /* hardening modulus                       */
            double uniax,    /* yield stresse                           */
            double vnu,      /* poisson's ratio                         */
            double sig2,
            double *dev,     /* elastic predicor projected onto yield   */
            double *epstn,   /* equivalent uniaxial plastic strain      */
            double *dlam);   /* increment of plastic multiplier         */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   forms the elasto-plastic consistent tangent material tensor        |
 *----------------------------------------------------------------------*/
void c1matp1(double e,       /* young's modulus                         */
             double fhard,   /* hardening modulus                       */
             double uniax,   /* yield stresse                           */
             double vnu,     /* poisson's ratio                         */
             double sig2,
             double *tau,    /* current stresses (local)                */
             double epstn,   /* equivalent uniaxial plastic strain      */ 
             double dlam,    /* increment of plastic multiplier         */
             double **cc);   /* material matrix to be calculated        */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   radial return for elements with mises material model               |
 |                                                                      |
 *----------------------------------------------------------------------*/
void c1rad (double e,        /* young's modulus                         */
            double fhard,    /* hardening modulus                       */
            double uniax,    /* yield stresse                           */
            double vnu,      /* poisson's ratio                         */
            double *sigma,   /* elastic predicor projected onto yield   */
            double *epstn,   /* equivalent uniaxial plastic strain      */
            double *dlam);   /* increment of plastic multiplier         */
void c1matp (double e,       /* young's modulus                         */
             double fhard,   /* hardening modulus                       */
             double vnu,     /* poisson's ratio                         */
             double *tau,    /* current stresses (local)                */
             double epstn,   /* equivalent uniaxial plastic strain      */ 
             double dlam,    /* increment of plastic multiplier         */
             double **d);    /* material matrix to be calculated        */
/*----------------------------------------------------------------------*
 | transformation of local stress vector to global axes         al 9/01 |
 *----------------------------------------------------------------------*/
void c1trss2global(double *s,       /* stress vector to be transformed  */
                   double g[6][6]); /* transformation matrix            */
/*----------------------------------------------------------------------*
 | transformation of global stress vector to local axes         al 9/01 |
 *----------------------------------------------------------------------*/
void c1trss2local(double *s,       /* stress vector to be transformed   */
                  double gi[6][6]);/* inverse of transformation matrix  */
/*----------------------------------------------------------------------*
 | transformation of local material-matrix to global axes       al 9/01 |
 *----------------------------------------------------------------------*/
void c1gld(double **d,     /* material matrix                           */
           double g[6][6]);/* transformation matrix                     */
/*----------------------------------------------------------------------*
 | program for evaluation of material transformation matricies  al 9/01 |
 *----------------------------------------------------------------------*/
void c1tram(double **xjm,   /* jacobian matrix r,s,t-direction          */
            double g[6][6], /* transformation matrix s(glob)=g*s(loc)   */
            double gi[6][6]);/* inverse of g          s(loc) =gi*s(glob)*/
/*----------------------------------------------------------------------*
 | integration points                                        al 6/01    |
 -----------------------------------------------------------------------|
 | hex-element                                                          |
 | coordinates and weighting factors of gauss-integration-points for    |
 | numerical integration                                                |
 *----------------------------------------------------------------------*/
void c1intg(ELEMENT         *ele,       /* actual element */
            C1_DATA         *data,      /* element data: wa, stress...  */
            int              option);   /* flag: calc.force, init...    */
/*----------------------------------------------------------------------*
 | update of strain paramenters                          al    9/01     |
 *----------------------------------------------------------------------*/
void c1upenh(  
             ELEMENT    *ele,
             double      edis[60],  
             double      ehdis[3][10],
             int         l1,
             int         l3
            );
/*----------------------------------------------------------------------*
 | evaluation of residual                                al    9/01     |
 *----------------------------------------------------------------------*/
void c1res( 
           double  *F,    /*  force vector integral (stress-resultants)*/
           double   fac,  /*  multiplier for numerical integration     */
           double **bop9, /*  b-operator matrix                        */
           int      iel,  /*  number nodes of element                  */
           double  *fieh, /*  internal force vector                    */
           int      l3);  /*                                           */
/*----------------------------------------------------------------------*
 | modify stiffness matrix and internal forces              al 9/01     |
 | for enhanced strain elements                                         |
 | update displacement parameters                                       |
 *----------------------------------------------------------------------*/
/*      subroutine c1rkefi (estif9,fieh,estif,fie,ns,lwah,nch,iel,nel,l1,
     *                    sbb,sbbi,sab,sba,sbai,saa,fiehi,wah)*/
void c1rkefi(
            ELEMENT    *ele,
            double **estif9, /* element stiffness-matrix          */
            double  **estif, /* modified element stiffness-matrix */
            double    *fieh, /* element residuals                 */
            double     *fie, /*                                   */
            int          l1
           ); 
/*----------------------------------------------------------------------*
 | include initial displacements to derivative operator     al 9/01     |
 | (eas-part)                                                           |
 *----------------------------------------------------------------------*/
void c1bdish(
            double    **bop,     /* b-operator  matrix */
            double      bn[3][10],/* bn-operator matrix */
            double      disd[9],
            int         iel,
            int         l3 
           ); 
/*----------------------------------------------------------------------*
 |program for evaluation of bop with extended deformations  al 9/01     |
 |(eas - kleine verz.)                                                  |
 *----------------------------------------------------------------------*/
void c1bop9(
            double    **bop9,/* b-operator matrix modified */
            double      bn1[3][10],
            double      fi[6][6],
            double      disd1[9],
            double      ehdis[3][10],
            double      det0,
            double      det1,
            double      e1,
            double      e2,
            double      e3,
            int         iel,
            int         l1,
            int         l3 
           ); 
/*----------------------------------------------------------------------*
 | determine the transformation matrix for orthogonal    al    9/01     |
 | systems and its inverse                                              |
 | jacobian matrix at point 0,0,0                                       |
 *----------------------------------------------------------------------*/
void c1t0( double   fi[6][6],    
           double   ff[6][6],    
           double **xjm0); 
/*----------------------------------------------------------------------*
 | shape functions and derivatives                               al 9/01|
 *----------------------------------------------------------------------*/
void c1_funct_deriv(double     *funct, 
                    double    **deriv, 
                    double      r, 
                    double      s,
                    double      t,
                    int         typ,
                    int         option);
/*----------------------------------------------------------------------*
 | calculate operator matrix at point r,s,t                  al 9/01    |
 *----------------------------------------------------------------------*/
void c1_jaco(double  **deriv,
             double    **xjm,
             double     *det,
             double    *xyze,
             int         iel);
/*----------------------------------------------------------------------*
 | calculate operator matrix at point r,s                    al 9/01    |
 *----------------------------------------------------------------------*/
void c1_bop(double    **b,
            double    **bn,
            double    **deriv,
            double    **xjm,
            double      det,
            int         iel);
/*----------------------------------------------------------------------*
 | include initial displacements to derivative operator      al 9/01    |
 *----------------------------------------------------------------------*/
void c1_bdis( double   **bop, /* b-operator matrix                      */
              double   *disd, /* displacement derivatives               */
              int        iel);/* number of nodes at actual element      */
/*----------------------------------------------------------------------*
 | compute displacement derivatives                          al 9/01    |
 *----------------------------------------------------------------------*/
void c1_disd(double    **bop, /* b-operator matrix                      */
             double    *edis, /* element displacements                  */ 
             double    *disd, /* displacement derivatives               */
             int         iel);/* number of element nodes                */
/*----------------------------------------------------------------------*
 |                                                          al 01/02    |
 |         control program for formulation of material law              |
 |                 select proper material law                           |
 |              and evaluation of element stresses                      |
 *----------------------------------------------------------------------*/
void c1_call_mat(ELEMENT   *ele,
                 MATERIAL  *mat, 
                 double **bop,
                 double **xjm,
                 int ip,       
                 double *stress,
                 double *strain,
                 double **d,
                 double *disd,
                 double g[6][6], /* transformation matrix s(glob)=g*s(loc)  */
                 double gi[6][6],/* inverse of g          s(loc) =gi*s(glob)*/
                 int istore, /* controls storing of new stresses to wa */
                 int newval);/* controls evaluation of new stresses    */
/*----------------------------------------------------------------------*
 |                                                          al 01/02    |
 |         control program for formulation of material law              |
 |                 select proper material law                           |
 |              and evaluation of element stresses                      |
 *----------------------------------------------------------------------*/
void c1_call_matd(ELEMENT   *ele,
                 MATERIAL  *mat, 
                 double **bop,
                 double **xjm,
                 int ip,       
                 double *stress,
                 double *strain,
                 double **d,
                 double *disd,
                 double g[6][6], /* transformation matrix s(glob)=g*s(loc)  */
                 double gi[6][6],/* inverse of g          s(loc) =gi*s(glob)*/
                 int istore,/* controls storing of new stresses to wa */
                 int newval);/* controls evaluation of new stresses    */
/*----------------------------------------------------------------------*
 | program for evaluation of principal stresses at given gauss  al 9/01 |
 | point for isoparametric brick element                                |
 |----------------------------------------------------------------------|
 | stress  (1) = sig-rr                                                 |
 | stress  (2) = sig-ss                                                 |
 | stress  (3) = sig-tt                                                 |
 | stress  (4) = sig-rs                                                 |
 | stress  (5) = sig-st                                                 |
 | stress  (6) = sig-tr                                                 |
 | stress  (7) = sig-  i   )                                            |
 | stress  (8) = sig- ii   ) --> principal stresses                     |
 | stress  (9) = sig-iii   )                                            |
 | stress (10) = alpha  (r,i)  )                                        |
 | stress (11) = alpha  (s,i)  )                                        |
 | stress (12) = alpha  (t,i)  )                                        |
 | stress (13) = alpha (r,ii)  )                                        |
 | stress (14) = alpha (s,ii)  ) -->   directions                       |
 | stress (15) = alpha (t,ii)  )                                        |
 | stress (16) = alpha(r,iii)  )                                        |
 | stress (17) = alpha(s,iii)  )                                        |
 | stress (18) = alpha(t,iii)  )                                        |
 *----------------------------------------------------------------------*/
void c1pstr(double    *srst,/* stresses at given gauss point            */
            double    *s123 /* principal stresses and direction at g.p. */
            );
/*----------------------------------------------------------------------*
 | calculate element stresses for postprocessing             al 9/01    |
 |----------------------------------------------------------------------|
 | stress  (1) = sig-rr                                                 |
 | stress  (2) = sig-ss                                                 |
 | stress  (3) = sig-tt                                                 |
 | stress  (4) = sig-rs                                                 |
 | stress  (5) = sig-st                                                 |
 | stress  (6) = sig-tr                                                 |
 | stress  (7) = sig-  i   )                                            |
 | stress  (8) = sig- ii   ) --> principal stresses                     |
 | stress  (9) = sig-iii   )                                            |
 | stress (10) = alpha  (r,i)  )                                        |
 | stress (11) = alpha  (s,i)  )                                        |
 | stress (12) = alpha  (t,i)  )                                        |
 | stress (13) = alpha (r,ii)  )                                        |
 | stress (14) = alpha (s,ii)  ) -->   directions                       |
 | stress (15) = alpha (t,ii)  )                                        |
 | stress (16) = alpha(r,iii)  )                                        |
 | stress (17) = alpha(s,iii)  )                                        |
 | stress (18) = alpha(t,iii)  )                                        |
 *----------------------------------------------------------------------*/
void c1_cstr(double    *srst,     /* element vector with stress resultants */
             double    *s123,
             double    *pstrs
            );
/*----------------------------------------------------------------------*
 | program to determine global coordinates referring to the     al 9/01 |
 | natural ones                                                         |
 *----------------------------------------------------------------------*/
void c1gcor( 
             double     *funct, /* value of form functions              */
             double     *xyze,  /* element coordinates                  */
             int         iel,   /* number of nodes                      */
             double     *gpcod  /* global coordinates of actual point   */
            );
/*----------------------------------------------------------------------*/
double c1rsn (
             int node,
             int irs,
             int iel  /* number of nodes */
             );
/*----------------------------------------------------------------------*
 |program to evaluate nth order legendre polynomial of degree 'n' at 'z'|
 *----------------------------------------------------------------------*/
void c1lgpl (
              int i,
              int n,
              double *zr,  
              double z,
              double *value
              );
void c1hxsm (
              int nir,
              int nis,
              int nit,  
              double rk,
              double sk,
              double tk,
              double f[8][27],
              double *fp,
              double *xgr,
              double *xgs,
              double *xgt
              );
/*----------------------------------------------------------------------*
 | calculate element stresses for postprocessing             al 9/01    |
 |----------------------------------------------------------------------|
 | stress  (1) = sig-rr                                                 |
 | stress  (2) = sig-ss                                                 |
 | stress  (3) = sig-tt                                                 |
 | stress  (4) = sig-rs                                                 |
 | stress  (5) = sig-st                                                 |
 | stress  (6) = sig-tr                                                 |
 | stress  (7) = sig-  i   )                                            |
 | stress  (8) = sig- ii   ) --> principal stresses                     |
 | stress  (9) = sig-iii   )                                            |
 | stress (10) = alpha  (r,i)  )                                        |
 | stress (11) = alpha  (s,i)  )                                        |
 | stress (12) = alpha  (t,i)  )                                        |
 | stress (13) = alpha (r,ii)  )                                        |
 | stress (14) = alpha (s,ii)  ) -->   directions                       |
 | stress (15) = alpha (t,ii)  )                                        |
 | stress (16) = alpha(r,iii)  )                                        |
 | stress (17) = alpha(s,iii)  )                                        |
 | stress (18) = alpha(t,iii)  )                                        |
 *----------------------------------------------------------------------*/
void c1_nstr(double    *srst,     /* element vector with stress resultants */
             double    *s123,
             double    *pstrs
            );
/*----------------------------------------------------------------------*
 | extrapolation of stress from gauss points to nodal points    al 9/01 |
 | for isoparametric brick elements                                     |
 *----------------------------------------------------------------------*/
void c1_sext(
            double **nostrs,
            double *funct,
            double **deriv,
            double **xjm,
            double *xyze,
            double **gpstress,
            double *xgr,
            double *xgs,
            double *xgt,
            int nir,
            int nis, 
            int nit,
            int iel /* number of nodes */
            );
/*----------------------------------------------------------------------*
 | usual stiffness matrix total lagrangian formulation           al 9/01|
 |----------------------------------------------------------------------|
 | S       -->  ELEMENT STIFFNESS MATRIX                                |
 | BS      -->  DERIVATIVE OPERATOR                                     |
 |              CONTENTS DEPENDS ON CALLING PROGRAM                     |
 |              IF T.L BS INCLUDES INITIAL DISPLACEMENT PARTS           |
 | D       -->  CONSTITUTIVE MATRIX                                     |
 | FAC     -->  INTEGRATION FACTOR                                      |
 | ND      -->  TOTAL NUMBER DEGREES OF FREEDOM OF ELEMENT              |
 | NEPS    -->  ACTUAL NUMBER OF STRAIN COMPONENTS                      |
 |              =3 PLANE STRESS/PLANE STRAIN                            |
 |              =4 AXISYMMETRIC                                         |
 |                                                                      |
 *----------------------------------------------------------------------*/
void c1_keku(double  **s, 
             double  **bs, 
             double  **d, 
             double    fac, 
             int       nd,
             int       neps);
/*----------------------------------------------------------------------*
 | evaluation of geometric stiffness-matrix (initial stress)     al 9/01|
 *----------------------------------------------------------------------*/
void c1vkg(  double  **s,   /* element stiffness-matrix                 */
             double  **bn,  /* b-operator matrix                        */
             double   *f,   /* force vector integral (stress-resultants)*/
             double    fac, /* multiplier for numerical integration     */
             int       iel);/* number of nodes at actual element         */
/*----------------------------------------------------------------------*
 | reorder stiffness matrix for 'gid'  hex20                     al 9/01|
 *----------------------------------------------------------------------*/
void c1kgid( double  **si,  /* element stiffness-matrix - convent.      */
             double  **so); /* element stiffness-matrix - gid.          */
/*----------------------------------------------------------------------*
 | reorder element force vector for 'gid'  hex20                 al 9/01|
 *----------------------------------------------------------------------*/
void c1fgid( double  *fi,  /* element stiffness-matrix - convent.       */
             double  *fo); /* element stiffness-matrix - gid.           */
/*----------------------------------------------------------------------*
 | program to establish local material law for brick element    al 9/01 |
 | stress-strain law for isotropic material.                            |
 *----------------------------------------------------------------------*/
void c1_mat_linel(double youngs,
                  double possionratio,
                  double **d);
/*----------------------------------------------------------------------*
 | program to establish local material law for brick element    al 9/01 |
 | stress-strain law for isotropic porous material.                     |
 *----------------------------------------------------------------------*/
void c1_mat_stvpor(MATERIAL  *mat,
                   double *matdata,
                   double **d);
/*----------------------------------------------------------------------*
 | program to establish local material law for brick element    al 9/01 |
 | stress-strain law for isotropic porous material.                     |
 *----------------------------------------------------------------------*/
void c1_matd_stvpor(MATERIAL  *mat,
                   double *matdata,
                   double **d);
/*----------------------------------------------------------------------*
 |  evaluate stresses for elastic material                      al 9/01 |
 *----------------------------------------------------------------------*/
void c1mefm(double *strain, /* global strains                           */
            double **d,     /* material matrices                        */
            double *stress);/* forces/moments/additional terms          */
/*----------------------------------------------------------------------*
 | constitutive matrix - forces - linear elastic- von Mises - 3D al 9/01|
 *----------------------------------------------------------------------*/
void c1_mat_plast_mises(double ym,      /* young's modulus              */
                        double pv,      /* poisson's ratio              */
                        double alfat,   /* temperature expansion factor */
                        double uniax,   /* yield stresse                */
                        double fhard,   /* hardening modulus            */
                        double gf,      /* fracture energy              */
                        ELEMENT   *ele, /* actual element               */
                        double **bop,   /* derivative operator          */
                        int ip,         /* integration point Id         */
                        double *stress, /*ele stress (-resultant) vector*/      
                        double **d,     /* material matrix              */
                        double  *disd,  /* displacement derivatives     */
                        double g[6][6], /* transformation matrix        */
                        double gi[6][6],/* inverse of g                 */
                        int istore,     /* controls storing of stresses */
                        int newval);    /* controls eval. of stresses   */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   evaluation of the elastic constitutive matrix (large strains)      |
 *----------------------------------------------------------------------*/
void c1mate(
            double detf,   /*   */
            double rmu,    /*   */
            double rk,     /*   */
            double bmu,    /*   */
            double sig2,   /*   */
            double *devn,  /*   */
            double **d);   /*   */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   radial return for elements with mises material model               |
 |   topic : radial return (large def. model)                           |
 *----------------------------------------------------------------------*/
void c1radg(
            double fhard,    /* hardening modulus                       */
            double uniax,    /* yield stresse                           */
            double bmu,      /*                        */
            double sig2,     /*                                         */
            double *dhard,   /* hardening modulus                       */
            double *dev,     /* elastic predicor projected onto yield   */
            double *epstn,   /* equivalent uniaxial plastic strain      */
            double *dlam);   /* increment of plastic multiplier         */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   forms the elasto-plastic consistent tangent material tensor        |
 |   (large strains)                                                    |
 *----------------------------------------------------------------------*/
void c1matpg(
             double dlam,    /* increment of plastic multiplier         */
             double detf,
             double rmu,
             double rk,
             double bmu,
             double sig2,
             double hard,
             double *devn,
             double **d);    /* material matrix to be calculated        */
/*----------------------------------------------------------------------*
 |   push forward(deformations)/pull-back(stresses)       al    9/01    |
 *----------------------------------------------------------------------*/
void c1pushf(
            double *be,   
            double *bet,  
            double *fn
            );  
/*----------------------------------------------------------------------*
 |                                                              al 9/01 |
 | topic :  calculate  e l a s t o p l a s t i c  stresses and          |
 |          stress increments  -  s t r e s s  p r o j e c t i o n      |
 |          for large deformation model                                 |
 |                                                                      |
 *----------------------------------------------------------------------*/
void c1elpag(
             double ym,      /* young's modulus              */
             double pv,      /* poisson's ratio              */
             double uniax,   /* yield stresse                */
             double fhard,   /* hardening modulus            */
             double *stress, /**/      
             double *sig,    /**/      
             double *fn,     /**/      
             double *fni,    /**/      
             double detf,    /**/      
             double **d,     /* material matrix              */
             double *epstn,  /**/
             int    *iupd,   /* controls storing of stresses */
             int    *yip);   /*    */
/*----------------------------------------------------------------------*
 |                                                              al 9/01 |
 | constitutive matrix - forces - plastic large strain - von Mises - 3D |
 *----------------------------------------------------------------------*/
void c1_mat_plast_mises_ls(
                        double ym,      /* young's modulus              */
                        double pv,      /* poisson's ratio              */
                        double alfat,   /* temperature expansion factor */
                        double uniax,   /* yield stresse                */
                        double fhard,   /* hardening modulus            */
                        ELEMENT   *ele, /* actual element               */
                        double **bop,   /* derivative operator          */
                        int ip,         /* integration point Id         */
                        double *stress, /*ele stress (-resultant) vector*/      
                        double **d,     /* material matrix              */
                        double  *disd,  /* displacement derivatives     */
                        double g[6][6], /* transformation matrix        */
                        double gi[6][6],/* inverse of g                 */
                        int istore,     /* controls storing of stresses */
                        int newval);    /* controls eval. of stresses   */
/*----------------------------------------------------------------------*
 | initialize the element                                    al 6/01    |
 *----------------------------------------------------------------------*/
void c1init(PARTITION *actpart,MATERIAL    *mat );
/*----------------------------------------------------------------------*
 | integration of linear stiffness ke for BRICK1 element     al 9/01    |
 *----------------------------------------------------------------------*/
void c1_cint(
             ELEMENT   *ele, 
             C1_DATA   *data, 
             MATERIAL  *mat,
             ARRAY     *estif_global, 
             double    *force,  /* global vector for internal forces (initialized!) */
             int        init
             );
/*----------------------------------------------------------------------*
 | evaluates element forces                              al    9/01     |
 *----------------------------------------------------------------------*/
void c1fi( double  *F,   /*  force vector integral (stress-resultants)  */
           double   fac, /*  multiplier for numerical integration       */
           double **bop, /*  b-operator matrix                          */
           int      nd,  /*  total number degrees of freedom of element */
           double  *fie);/*  internal force vector                      */
/*----------------------------------------------------------------------*
 | program for evaluation of material transformation matricies  al 9/01 |
 *----------------------------------------------------------------------*/
void c1tram( double **xjm,   /* jacobian matrix r,s,t-direction         */
             double g[6][6], /* transformation matrix s(glob)=g*s(loc)  */
             double gi[6][6]);/* inverse of g          s(loc) =gi*s(glob)*/
/*----------------------------------------------------------------------*
 | transformation of local material-matrix to global axes       al 9/01 |
 *----------------------------------------------------------------------*/
void c1gld(double **d,     /* material matrix                           */
           double g[6][6]);/* transformation matrix                     */
/*----------------------------------------------------------------------*
 | transformation of global stress vector to local axes         al 9/01 |
 *----------------------------------------------------------------------*/
void c1trss2local(double *s,       /* stress vector to be transformed   */
                  double gi[6][6]);/* inverse of transformation matrix  */
/*----------------------------------------------------------------------*
 | transformation of local stress vector to global axes         al 9/01 |
 *----------------------------------------------------------------------*/
void c1trss2global(double *s,       /* stress vector to be transformed  */
                   double g[6][6]); /* transformation matrix            */
/*----------------------------------------------------------------------*
 | calculate global coordinates referring to the natural ones   al 9/01 |
 *----------------------------------------------------------------------*/
void  c1gcor( 
             double     *funct, /* value of form functions              */
             double     *xyze,  /* element coordinates                  */
             int         iel,   /* number of nodes                      */
             double     *gpcod  /* global coordinates of actual point   */
            );

/*----------------------------------------------------------------------*
 | calculate element loads (edge, surface, volume)              al 9/01 |
 *----------------------------------------------------------------------*/
void c1_eleload(
             ELEMENT   *ele, 
             C1_DATA   *data, 
             MATERIAL  *mat,
             double    *loadvec,  
             int        init
             );
/*----------------------------------------------------------------------*
 | calculate constitutive matrix                                al 9/01 |
 *----------------------------------------------------------------------*/
void c1_mat_mfoc(  MATERIAL  *mat,
                   double *matdata,
                   double **d);
/*----------------------------------------------------------------------*
 | calculate constitutive matrix                                al 9/01 |
 *----------------------------------------------------------------------*/
void c1_mat_mfcc(  MATERIAL  *mat,
                   double *matdata,
                   double **d);
/*----------------------------------------------------------------------*
 | calculate derivative of constitutive matrix                  al 9/01 |
 *----------------------------------------------------------------------*/
void c1_matd_mfoc( MATERIAL  *mat,
                   double *matdata,
                   double **d);
/*----------------------------------------------------------------------*
 | calculate derivative of constitutive matrix                  al 9/01 |
 *----------------------------------------------------------------------*/
void c1_matd_mfcc( MATERIAL  *mat,
                   double *matdata,
                   double **d);
/*----------------------------------------------------------------------*
 | calculate constitutive matrix                                al 9/01 |
 *----------------------------------------------------------------------*/
void c1_mat_elorth(double   emod1 ,
                   double   emod2 ,
                   double   emod3 ,
                   double   xnue23,
                   double   xnue13,
                   double   xnue12,
                   double   gmod12,
                   double   gmod23,
                   double   gmod13,
                   double      **c);
/*----------------------------------------------------------------------*
 | - delamination hashin plasticity based - 3D                   al 9/01|
 *----------------------------------------------------------------------*/
void c1_mat_plast_hashdel
                        (
                        MATERIAL  *mat, /* material data                */
                        ELEMENT   *ele, /* actual element               */
                        double **bop,   /* derivative operator          */
                        int       ip,   /* integration point            */
                        double *stress, /*ele stress (-resultant) vector*/      
                        double **c,     /* material matrix              */
                        double  *disd,  /* displacement derivatives     */
                        double g[6][6], /* transformation matrix        */
                        double gi[6][6],/* inverse of g                 */
                        int istore,     /* controls storing of stresses */
                        int newval);    /* controls eval. of stresses   */
/*----------------------------------------------------------------------*
 | prototypes for fortran routines                               al 9/01|
 *----------------------------------------------------------------------*/
void fortranpow (double *V,double *R,double *RE);
void MXMAB (double *A,double *B,double *R,int *NZA,int *NSA,int *NSB);
void mxmab (double *A,double *B,double *R,int *NZA,int *NSA,int *NSB);
void MXMABT (double *A,int *M,int *N,int *L,double *B,double *R);
void mxmabt (double *A,int *M,int *N,int *L,double *B,double *R);
void mxmatb (double *A,double *B,double *R,int *NZA,int *NSA,
             int *NSB,double *XNULL);
void c1jacb (double *A,double *V);
void solveq (double *A,double *B,int *NN,int *LL);
void c1inv6(double *A,double *B,int *IRC);
void c1invf (double *FN,double *FNI,double *DETF);
void c1ab (double *A,double *B,double *R,int *NZA,int *NSA,int *NSB,
           double *XNULL);                                 


#endif
/*! @} (documentation module close)*/
