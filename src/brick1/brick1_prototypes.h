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
void c1_eps( DOUBLE   *disd,                /* displacement derivatives */
             DOUBLE   *eps,                 /* strain vector            */
             INT      iform);        /* index for nonlinear formulation */
/*----------------------------------------------------------------------*
 | evaluates element forces                              al    9/01     |
 *----------------------------------------------------------------------*/
void c1fi( DOUBLE  *F,   /*  force vector integral (stress-resultants)  */
           DOUBLE   fac, /*  multiplier for numerical integration       */
           DOUBLE **bop, /*  b-operator matrix                          */
           INT      nd,  /*  total number degrees of freedom of element */
           DOUBLE  *fie);/*  internal force vector                      */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   radial return for elements with mises material model               |
 |                                                                      |
 *----------------------------------------------------------------------*/
void c1rad1(DOUBLE e,        /* young's modulus                         */
            DOUBLE eh,       /* hardening modulus                       */
            DOUBLE uniax,    /* yield stresse                           */
            DOUBLE vnu,      /* poisson's ratio                         */
            DOUBLE sig2,
            DOUBLE *dev,     /* elastic predicor projected onto yield   */
            DOUBLE *epstn,   /* equivalent uniaxial plastic strain      */
            DOUBLE *dlam);   /* increment of plastic multiplier         */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   forms the elasto-plastic consistent tangent material tensor        |
 *----------------------------------------------------------------------*/
void c1matp1(DOUBLE e,       /* young's modulus                         */
             DOUBLE fhard,   /* hardening modulus                       */
             DOUBLE vnu,     /* poisson's ratio                         */
             DOUBLE sig2,
             DOUBLE *tau,    /* current stresses (local)                */
             DOUBLE epstn,   /* equivalent uniaxial plastic strain      */ 
             DOUBLE dlam,    /* increment of plastic multiplier         */
             DOUBLE **cc);   /* material matrix to be calculated        */
/*----------------------------------------------------------------------*
 | transformation of local stress vector to global axes         al 9/01 |
 *----------------------------------------------------------------------*/
void c1trss2global(DOUBLE *s,       /* stress vector to be transformed  */
                   DOUBLE g[6][6]); /* transformation matrix            */
/*----------------------------------------------------------------------*
 | transformation of global stress vector to local axes         al 9/01 |
 *----------------------------------------------------------------------*/
void c1trss2local(DOUBLE *s,       /* stress vector to be transformed   */
                  DOUBLE gi[6][6]);/* inverse of transformation matrix  */
/*----------------------------------------------------------------------*
 | transformation of local material-matrix to global axes       al 9/01 |
 *----------------------------------------------------------------------*/
void c1gld(DOUBLE **d,     /* material matrix                           */
           DOUBLE g[6][6]);/* transformation matrix                     */
/*----------------------------------------------------------------------*
 | program for evaluation of material transformation matricies  al 9/01 |
 *----------------------------------------------------------------------*/
void c1tram(DOUBLE **xjm,   /* jacobian matrix r,s,t-direction          */
            DOUBLE g[6][6], /* transformation matrix s(glob)=g*s(loc)   */
            DOUBLE gi[6][6]);/* inverse of g          s(loc) =gi*s(glob)*/
/*----------------------------------------------------------------------*
 | integration points                                        al 6/01    |
 -----------------------------------------------------------------------|
 | hex-element                                                          |
 | coordinates and weighting factors of gauss-integration-points for    |
 | numerical integration                                                |
 *----------------------------------------------------------------------*/
void c1intg(ELEMENT         *ele,       /* actual element */
            C1_DATA         *data);     /* element data: wa, stress...  */
/*----------------------------------------------------------------------*
 | update of strain paramenters                          al    9/01     |
 *----------------------------------------------------------------------*/
void c1upenh(  
             ELEMENT    *ele,
             DOUBLE      edis[60],  
             DOUBLE      ehdis[3][10],
             INT         l1,
             INT         l3
            );
/*----------------------------------------------------------------------*
 | evaluation of residual                                al    9/01     |
 *----------------------------------------------------------------------*/
void c1res( 
           DOUBLE  *F,    /*  force vector integral (stress-resultants)*/
           DOUBLE   fac,  /*  multiplier for numerical integration     */
           DOUBLE **bop9, /*  b-operator matrix                        */
           INT      iel,  /*  number nodes of element                  */
           DOUBLE  *fieh, /*  internal force vector                    */
           INT      l3);  /*                                           */
/*----------------------------------------------------------------------*
 | modify stiffness matrix and internal forces              al 9/01     |
 | for enhanced strain elements                                         |
 | update displacement parameters                                       |
 *----------------------------------------------------------------------*/
/*      subroutine c1rkefi (estif9,fieh,estif,fie,ns,lwah,nch,iel,nel,l1,
     *                    sbb,sbbi,sab,sba,sbai,saa,fiehi,wah)*/
void c1rkefi(
            ELEMENT    *ele,
            DOUBLE **estif9, /* element stiffness-matrix          */
            DOUBLE  **estif, /* modified element stiffness-matrix */
            DOUBLE    *fieh, /* element residuals                 */
            DOUBLE     *fie, /*                                   */
            INT          l1
           ); 
/*----------------------------------------------------------------------*
 | include initial displacements to derivative operator     al 9/01     |
 | (eas-part)                                                           |
 *----------------------------------------------------------------------*/
void c1bdish(
            DOUBLE    **bop,     /* b-operator  matrix */
            DOUBLE      bn[3][10],/* bn-operator matrix */
            DOUBLE      disd[9],
            INT         iel,
            INT         l3 
           ); 
/*----------------------------------------------------------------------*
 |program for evaluation of bop with extended deformations  al 9/01     |
 |(eas - kleine verz.)                                                  |
 *----------------------------------------------------------------------*/
void c1bop9(
            DOUBLE    **bop9,/* b-operator matrix modified */
            DOUBLE      bn1[3][10],
            DOUBLE      fi[6][6],
            DOUBLE      disd1[9],
            DOUBLE      ehdis[3][10],
            DOUBLE      det0,
            DOUBLE      det1,
            DOUBLE      e1,
            DOUBLE      e2,
            DOUBLE      e3,
            INT         iel,
            INT         l1,
            INT         l3 
           ); 
/*----------------------------------------------------------------------*
 | determine the transformation matrix for orthogonal    al    9/01     |
 | systems and its inverse                                              |
 | jacobian matrix at point 0,0,0                                       |
 *----------------------------------------------------------------------*/
void c1t0( DOUBLE   fi[6][6],    
           DOUBLE   ff[6][6],    
           DOUBLE **xjm0); 
/*----------------------------------------------------------------------*
 | shape functions and derivatives                               al 9/01|
 *----------------------------------------------------------------------*/
void c1_funct_deriv(DOUBLE     *funct, 
                    DOUBLE    **deriv, 
                    DOUBLE      r, 
                    DOUBLE      s,
                    DOUBLE      t,
                    INT         typ,
                    INT         option);
/*----------------------------------------------------------------------*
 | calculate operator matrix at point r,s,t                  al 9/01    |
 *----------------------------------------------------------------------*/
void c1_jaco(DOUBLE  **deriv,
             DOUBLE    **xjm,
             DOUBLE     *det,
             DOUBLE    *xyze,
             INT         iel);
/*----------------------------------------------------------------------*
 | calculate operator matrix at point r,s                    al 9/01    |
 *----------------------------------------------------------------------*/
void c1_bop(DOUBLE    **b,
            DOUBLE    **bn,
            DOUBLE    **deriv,
            DOUBLE    **xjm,
            DOUBLE      det,
            INT         iel);
/*----------------------------------------------------------------------*
 | include initial displacements to derivative operator      al 9/01    |
 *----------------------------------------------------------------------*/
void c1_bdis( DOUBLE   **bop, /* b-operator matrix                      */
              DOUBLE   *disd, /* displacement derivatives               */
              INT        iel);/* number of nodes at actual element      */
/*----------------------------------------------------------------------*
 | compute displacement derivatives                          al 9/01    |
 *----------------------------------------------------------------------*/
void c1_disd(DOUBLE    **bop, /* b-operator matrix                      */
             DOUBLE    *edis, /* element displacements                  */ 
             DOUBLE    *disd, /* displacement derivatives               */
             INT         iel);/* number of element nodes                */
/*----------------------------------------------------------------------*
 |                                                          al 01/02    |
 |         control program for formulation of material law              |
 |                 select proper material law                           |
 |              and evaluation of element stresses                      |
 *----------------------------------------------------------------------*/
void c1_call_mat(ELEMENT   *ele,
                 MATERIAL  *mat, 
                 INT ip,       
                 DOUBLE *stress,
                 DOUBLE *strain,
                 DOUBLE **d,
                 DOUBLE *disd,
                 DOUBLE g[6][6], 
                 DOUBLE gi[6][6],
                 INT istore,
                 INT newval);
/*----------------------------------------------------------------------*
 |                                                          al 01/02    |
 |         control program for formulation of material law              |
 |                 select proper material law                           |
 |              and evaluation of element stresses                      |
 *----------------------------------------------------------------------*/
void c1_call_matd(ELEMENT   *ele,
                 MATERIAL  *mat, 
                 DOUBLE *stress,
                 DOUBLE *strain,
                 DOUBLE **d,
                 DOUBLE g[6][6]);
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
void c1pstr(DOUBLE    *srst,/* stresses at given gauss point            */
            DOUBLE    *s123 /* principal stresses and direction at g.p. */
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
void c1_cstr(DOUBLE    *srst,     /* element vector with stress resultants */
             DOUBLE    *s123,
             DOUBLE    *pstrs
            );
/*----------------------------------------------------------------------*
 | program to determine global coordinates referring to the     al 9/01 |
 | natural ones                                                         |
 *----------------------------------------------------------------------*/
void c1gcor( 
             DOUBLE     *funct, /* value of form functions              */
             DOUBLE     *xyze,  /* element coordinates                  */
             INT         iel,   /* number of nodes                      */
             DOUBLE     *gpcod  /* global coordinates of actual point   */
            );
/*----------------------------------------------------------------------*/
DOUBLE c1rsn (
             INT node,
             INT irs
             );
/*----------------------------------------------------------------------*
 |program to evaluate nth order legendre polynomial of degree 'n' at 'z'|
 *----------------------------------------------------------------------*/
void c1lgpl (
              INT i,
              INT n,
              DOUBLE *zr,  
              DOUBLE z,
              DOUBLE *value
              );
void c1hxsm (
              INT nir,
              INT nis,
              INT nit,  
              DOUBLE rk,
              DOUBLE sk,
              DOUBLE tk,
              DOUBLE f[8][27],
              DOUBLE *fp,
              DOUBLE *xgr,
              DOUBLE *xgs,
              DOUBLE *xgt
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
void c1_nstr(DOUBLE    *srst,     /* element vector with stress resultants */
             DOUBLE    *s123,
             DOUBLE    *pstrs
            );
/*----------------------------------------------------------------------*
 | extrapolation of stress from gauss points to nodal points    al 9/01 |
 | for isoparametric brick elements                                     |
 *----------------------------------------------------------------------*/
void c1_sext(
            DOUBLE nostrs[20][26],
            DOUBLE *funct,
            DOUBLE **deriv,
            DOUBLE **xjm,
            DOUBLE *xyze,
            DOUBLE gpstress[27][26],
            DOUBLE *xgr,
            DOUBLE *xgs,
            DOUBLE *xgt,
            INT nir,
            INT nis, 
            INT nit,
            INT iel /* number of nodes */
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
void c1_keku(DOUBLE  **s, 
             DOUBLE  **bs, 
             DOUBLE  **d, 
             DOUBLE    fac, 
             INT       nd,
             INT       neps);
/*----------------------------------------------------------------------*
 | evaluation of geometric stiffness-matrix (initial stress)     al 9/01|
 *----------------------------------------------------------------------*/
void c1vkg(  DOUBLE  **s,   /* element stiffness-matrix                 */
             DOUBLE  **bn,  /* b-operator matrix                        */
             DOUBLE   *f,   /* force vector integral (stress-resultants)*/
             DOUBLE    fac, /* multiplier for numerical integration     */
             INT       iel);/* number of nodes at actual element         */
/*----------------------------------------------------------------------*
 | reorder stiffness matrix for 'gid'  hex20                     al 9/01|
 *----------------------------------------------------------------------*/
void c1kgid( DOUBLE  **si,  /* element stiffness-matrix - convent.      */
             DOUBLE  **so); /* element stiffness-matrix - gid.          */
/*----------------------------------------------------------------------*
 | reorder element force vector for 'gid'  hex20                 al 9/01|
 *----------------------------------------------------------------------*/
void c1fgid( DOUBLE  *fi,  /* element stiffness-matrix - convent.       */
             DOUBLE  *fo); /* element stiffness-matrix - gid.           */
/*----------------------------------------------------------------------*
 | program to establish local material law for brick element    al 9/01 |
 | stress-strain law for isotropic material.                            |
 *----------------------------------------------------------------------*/
void c1_mat_linel(DOUBLE youngs,
                  DOUBLE possionratio,
                  DOUBLE **d);
/*----------------------------------------------------------------------*
 | program to establish local material law for brick element    al 9/01 |
 | stress-strain law for isotropic porous material.                     |
 *----------------------------------------------------------------------*/
void c1_mat_stvpor(MATERIAL  *mat,
                   DOUBLE *matdata,
                   DOUBLE **d);
/*----------------------------------------------------------------------*
 | program to establish local material law for brick element    al 9/01 |
 | stress-strain law for isotropic porous material.                     |
 *----------------------------------------------------------------------*/
void c1_matd_stvpor(MATERIAL  *mat,
                   DOUBLE *matdata,
                   DOUBLE **d);
/*----------------------------------------------------------------------*
 |  evaluate stresses for elastic material                      al 9/01 |
 *----------------------------------------------------------------------*/
void c1mefm(DOUBLE *strain, /* global strains                           */
            DOUBLE **d,     /* material matrices                        */
            DOUBLE *stress);/* forces/moments/additional terms          */
/*----------------------------------------------------------------------*
 | constitutive matrix - forces - linear elastic- von Mises - 3D al 9/01|
 *----------------------------------------------------------------------*/
void c1_mat_plast_mises(DOUBLE ym,      /* young's modulus              */
                        DOUBLE pv,      /* poisson's ratio              */
                        DOUBLE uniax,   /* yield stresse                */
                        DOUBLE fhard,   /* hardening modulus            */
                        ELEMENT   *ele, /* actual element               */
                        INT ip,         /* integration point Id         */
                        DOUBLE *stress, /*ele stress (-resultant) vector*/      
                        DOUBLE **d,     /* material matrix              */
                        DOUBLE  *disd,  /* displacement derivatives     */
                        DOUBLE g[6][6], /* transformation matrix        */
                        DOUBLE gi[6][6],/* inverse of g                 */
                        INT istore,     /* controls storing of stresses */
                        INT newval);    /* controls eval. of stresses   */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   evaluation of the elastic constitutive matrix (large strains)      |
 *----------------------------------------------------------------------*/
void c1mate(
            DOUBLE detf,   /*   */
            DOUBLE rk,     /*   */
            DOUBLE bmu,    /*   */
            DOUBLE sig2,   /*   */
            DOUBLE *devn,  /*   */
            DOUBLE **d);   /*   */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   radial return for elements with mises material model               |
 |   topic : radial return (large def. model)                           |
 *----------------------------------------------------------------------*/
void c1radg(
            DOUBLE fhard,    /* hardening modulus                       */
            DOUBLE uniax,    /* yield stresse                           */
            DOUBLE bmu,      /*                        */
            DOUBLE sig2,     /*                                         */
            DOUBLE *dhard,   /* hardening modulus                       */
            DOUBLE *dev,     /* elastic predicor projected onto yield   */
            DOUBLE *epstn,   /* equivalent uniaxial plastic strain      */
            DOUBLE *dlam);   /* increment of plastic multiplier         */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   forms the elasto-plastic consistent tangent material tensor        |
 |   (large strains)                                                    |
 *----------------------------------------------------------------------*/
void c1matpg(
             DOUBLE dlam,    /* increment of plastic multiplier         */
             DOUBLE detf,
             DOUBLE rk,
             DOUBLE bmu,
             DOUBLE sig2,
             DOUBLE hard,
             DOUBLE *devn,
             DOUBLE **d);    /* material matrix to be calculated        */
/*----------------------------------------------------------------------*
 |   push forward(deformations)/pull-back(stresses)       al    9/01    |
 *----------------------------------------------------------------------*/
void c1pushf(
            DOUBLE *be,   
            DOUBLE *bet,  
            DOUBLE *fn
            );  
void c1elpag(
             DOUBLE ym,      
             DOUBLE pv,      
             DOUBLE uniax,   
             DOUBLE fhard,   
             DOUBLE *stress, 
             DOUBLE *sig,    
             DOUBLE *fn,     
             DOUBLE *fni,    
             DOUBLE detf,    
             DOUBLE **d,     
             DOUBLE *epstn,  
             INT    *iupd,   
             INT    *yip);    
/*----------------------------------------------------------------------*
 |                                                              al 9/01 |
 | constitutive matrix - forces - plastic large strain - von Mises - 3D |
 *----------------------------------------------------------------------*/
void c1_mat_plast_mises_ls(
                        DOUBLE ym,      /* young's modulus              */
                        DOUBLE pv,      /* poisson's ratio              */
                        DOUBLE uniax,   /* yield stresse                */
                        DOUBLE fhard,   /* hardening modulus            */
                        ELEMENT   *ele, /* actual element               */
                        INT ip,         /* integration point Id         */
                        DOUBLE *stress, /*ele stress (-resultant) vector*/      
                        DOUBLE **d,     /* material matrix              */
                        DOUBLE  *disd,  /* displacement derivatives     */
                        INT istore,     /* controls storing of stresses */
                        INT newval);    /* controls eval. of stresses   */
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
             DOUBLE    *force,  /* global vector for internal forces (initialized!) */
             INT        init
             );
/*----------------------------------------------------------------------*
 | evaluates element forces                              al    9/01     |
 *----------------------------------------------------------------------*/
void c1fi( DOUBLE  *F,   /*  force vector integral (stress-resultants)  */
           DOUBLE   fac, /*  multiplier for numerical integration       */
           DOUBLE **bop, /*  b-operator matrix                          */
           INT      nd,  /*  total number degrees of freedom of element */
           DOUBLE  *fie);/*  internal force vector                      */
/*----------------------------------------------------------------------*
 | program for evaluation of material transformation matricies  al 9/01 |
 *----------------------------------------------------------------------*/
void c1tram( DOUBLE **xjm,   /* jacobian matrix r,s,t-direction         */
             DOUBLE g[6][6], /* transformation matrix s(glob)=g*s(loc)  */
             DOUBLE gi[6][6]);/* inverse of g          s(loc) =gi*s(glob)*/
/*----------------------------------------------------------------------*
 | transformation of local material-matrix to global axes       al 9/01 |
 *----------------------------------------------------------------------*/
void c1gld(DOUBLE **d,     /* material matrix                           */
           DOUBLE g[6][6]);/* transformation matrix                     */
/*----------------------------------------------------------------------*
 | transformation of global stress vector to local axes         al 9/01 |
 *----------------------------------------------------------------------*/
void c1trss2local(DOUBLE *s,       /* stress vector to be transformed   */
                  DOUBLE gi[6][6]);/* inverse of transformation matrix  */
/*----------------------------------------------------------------------*
 | transformation of local stress vector to global axes         al 9/01 |
 *----------------------------------------------------------------------*/
void c1trss2global(DOUBLE *s,       /* stress vector to be transformed  */
                   DOUBLE g[6][6]); /* transformation matrix            */
/*----------------------------------------------------------------------*
 | calculate global coordinates referring to the natural ones   al 9/01 |
 *----------------------------------------------------------------------*/
void  c1gcor( 
             DOUBLE     *funct, /* value of form functions              */
             DOUBLE     *xyze,  /* element coordinates                  */
             INT         iel,   /* number of nodes                      */
             DOUBLE     *gpcod  /* global coordinates of actual point   */
            );

void c1_lint(
             INT        ngr,       INT ngs,        INT ngt, 
             DOUBLE    *xgp,   DOUBLE *wgx,    DOUBLE *ygp,
             DOUBLE    *wgy,   DOUBLE *zgp,    DOUBLE *wgz,
             DOUBLE   *xyze, DOUBLE *funct, DOUBLE **deriv, DOUBLE **xjm,
             INT        iel,    INT ngnode,       INT *shn,  RSTF rstgeo,
             INT    *lonoff,  DOUBLE *lval, 
             DOUBLE **eload
             );
/*----------------------------------------------------------------------*
 | calculate element loads (edge, surface, volume)              al 9/01 |
 *----------------------------------------------------------------------*/
void c1_eleload(
             ELEMENT   *ele, 
             C1_DATA   *data, 
             DOUBLE    *loadvec,  
             INT        init
             );
/*----------------------------------------------------------------------*
 | calculate constitutive matrix                                al 9/01 |
 *----------------------------------------------------------------------*/
void c1_mat_mfoc(  MATERIAL  *mat,
                   DOUBLE *matdata,
                   DOUBLE **d);
/*----------------------------------------------------------------------*
 | calculate constitutive matrix                                al 9/01 |
 *----------------------------------------------------------------------*/
void c1_mat_mfcc(  MATERIAL  *mat,
                   DOUBLE *matdata,
                   DOUBLE **d);
/*----------------------------------------------------------------------*
 | calculate derivative of constitutive matrix                  al 9/01 |
 *----------------------------------------------------------------------*/
void c1_matd_mfoc( MATERIAL  *mat,
                   DOUBLE *matdata,
                   DOUBLE **d);
/*----------------------------------------------------------------------*
 | calculate derivative of constitutive matrix                  al 9/01 |
 *----------------------------------------------------------------------*/
void c1_matd_mfcc( MATERIAL  *mat,
                   DOUBLE *matdata,
                   DOUBLE **d);
/*----------------------------------------------------------------------*
 | calculate constitutive matrix                                al 9/01 |
 *----------------------------------------------------------------------*/
void c1_mat_elorth(DOUBLE   emod1 ,
                   DOUBLE   emod2 ,
                   DOUBLE   emod3 ,
                   DOUBLE   xnue23,
                   DOUBLE   xnue13,
                   DOUBLE   xnue12,
                   DOUBLE   gmod12,
                   DOUBLE   gmod23,
                   DOUBLE   gmod13,
                   DOUBLE      **c);
/*----------------------------------------------------------------------*
 | - delamination hashin plasticity based - 3D                   al 9/01|
 *----------------------------------------------------------------------*/
void c1_mat_plast_hashdel
                        (
                        MATERIAL  *mat, /* material data                */
                        ELEMENT   *ele, /* actual element               */
                        DOUBLE **bop,   /* derivative operator          */
                        INT       ip,   /* integration point            */
                        DOUBLE *stress, /*ele stress (-resultant) vector*/      
                        DOUBLE **c,     /* material matrix              */
                        DOUBLE  *disd,  /* displacement derivatives     */
                        DOUBLE g[6][6], /* transformation matrix        */
                        DOUBLE gi[6][6],/* inverse of g                 */
                        INT istore,     /* controls storing of stresses */
                        INT newval);    /* controls eval. of stresses   */
/*----------------------------------------------------------------------*
 |  integration routine for BRICK1 element                     al 6/01  |
 *----------------------------------------------------------------------*/
void c1_oint(
             ELEMENT   *ele, 
             C1_DATA   *data, 
             MATERIAL  *mat,
             DOUBLE    *retval,  /* return value */
             INT        init
             );
/*----------------------------------------------------------------------*
 | get density out of material law                          al 08/02    |
 *----------------------------------------------------------------------*/
void c1_getdensity(MATERIAL   *mat, DOUBLE *density);
/*----------------------------------------------------------------------*
 | prototypes for fortran routines                               al 9/01|
 *----------------------------------------------------------------------*/
void fortranpow (DOUBLE *V,DOUBLE *R,DOUBLE *RE);
void MXMAB (DOUBLE *A,DOUBLE *B,DOUBLE *R,INT *NZA,INT *NSA,INT *NSB);
void mxmab (DOUBLE *A,DOUBLE *B,DOUBLE *R,INT *NZA,INT *NSA,INT *NSB);
void MXMABT (DOUBLE *A,INT *M,INT *N,INT *L,DOUBLE *B,DOUBLE *R);
void mxmabt (DOUBLE *A,INT *M,INT *N,INT *L,DOUBLE *B,DOUBLE *R);
void mxmatb (DOUBLE *A,DOUBLE *B,DOUBLE *R,INT *NZA,INT *NSA,
             INT *NSB,DOUBLE *XNULL);
void c1jacb (DOUBLE *A,DOUBLE *V);
void solveq (DOUBLE *A,DOUBLE *B,INT *NN,INT *LL);
void c1inv6(DOUBLE *A,DOUBLE *B,INT *IRC);
void c1invf (DOUBLE *FN,DOUBLE *FNI,DOUBLE *DETF);
void c1ab (DOUBLE *A,DOUBLE *B,DOUBLE *R,INT *NZA,INT *NSA,INT *NSB,
           DOUBLE *XNULL);                                 


#endif
/*! @} (documentation module close)*/
