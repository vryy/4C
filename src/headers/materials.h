/*----------------------------------------------------------------------*
 | materials                                              m.gee 4/01    |
 | structure to hold all types of material laws                         |
 *----------------------------------------------------------------------*/
typedef struct _MATERIAL
{
     INT                       Id;           /* Id of the material */

     enum _MATERIAL_TYP        mattyp;       /* type of material */

     union
     {
     struct _STVENANT         *stvenant;     /* St. Venant-Kirchhoff material */
     struct _PL_MISES         *pl_mises;     /* von Mises material */
     struct _PL_FOAM          *pl_foam;      /* foam material - large strains */
     struct _PL_MISES_LS      *pl_mises_ls;  /* von Mises material - large strains*/
     struct _PL_DP            *pl_dp;        /* Drucker Prager material */
     struct _PL_EPC           *pl_epc;       /* elastoplastic concrete material */
     struct _STVENPOR         *stvenpor;     /* porous St. Ven.-Kirch. material */
     struct _PL_POR_MISES     *pl_por_mises; /* porous von Mises material */
     struct _NEO_HOOKE        *neohooke;     /* Neo-Hooke material */
     struct _COMPOGDEN        *compogden;    /* compressible ogden hyperelastic material */
     struct _FLUID            *fluid;        /* fluid material */
     struct _PL_HASH          *pl_hash;      /* elpl. hashin delamination material */
     struct _EL_ORTH          *el_orth;      /* elastic orthotropic material */
     struct _MFOC             *mfoc;         /* metal foam, open cell  */
     struct _MFCC             *mfcc;         /* metal foam, closed cell  */
     struct _NHMFCC           *nhmfcc;       /* foam, closed cell, based on modified Neo Hook */
     struct _MULTI_LAYER      *multi_layer;  /* multi layer material*/
     struct _ORTHOTROPIC      *orthotropic;  /* linear elastic, orthotropic material*/
     }                         m;            /* union pointer to material specific structure */

} MATERIAL;





/*----------------------------------------------------------------------*
 | St. Venant-Kirchhoff material                          m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _STVENANT
{
     DOUBLE                    youngs;         /* Young's modulus */
     DOUBLE                    possionratio;   /* Possion ratio */
     DOUBLE                    density;        /* material specific weight */
} STVENANT;


/*----------------------------------------------------------------------*
 | porous St. Venant-Kirchhoff material (optimization)       al 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _STVENPOR
{
     DOUBLE                    youngs;         /* Young's modulus */
     DOUBLE                    possionratio;   /* Possion ratio */
     DOUBLE                    density;        /* material specific weight */
     DOUBLE                    refdens;        /* reference density */
     DOUBLE                    exponent;       /* material parameter */
} STVENPOR;


/*----------------------------------------------------------------------*
 | Neo Hooke material                                     m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _NEO_HOOKE
{
     DOUBLE                    youngs;         /* Young's modulus */
     DOUBLE                    possionratio;   /* Possion ratio */
     DOUBLE                    density;        /* material specific weight */
} NEO_HOOKE;

/*----------------------------------------------------------------------*
 | compressible ogden material                            m.gee 6/03    |
 *----------------------------------------------------------------------*/
typedef struct _COMPOGDEN
{
     INT                       init;           /* init flag */
     DOUBLE                    nue;            /* Possion ratio */
     DOUBLE                    beta;           /* the unphysical material constant called beta */
     DOUBLE                    alfap[3];       /* three parameters alfap */
     DOUBLE                    mup[3];         /* three parameters nuep */  
     DOUBLE                    density;        /* material specific weight */
     DOUBLE                    lambda;         /* 1. lame constant */
     DOUBLE                    kappa;          /* bulkmodulus */
} COMPOGDEN;


/*----------------------------------------------------------------------*
 | fluid material                                         m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _FLUID
{
     DOUBLE                    viscosity;
     DOUBLE                    density;
     DOUBLE                    gamma;     /* surface tension coeficient */
} FLUID;
/*----------------------------------------------------------------------*
 | plastic mises material                              a.lipka 17/05    |
 *----------------------------------------------------------------------*/
typedef struct _PL_MISES
{
     DOUBLE                    youngs;        /* Young's modulus */
     DOUBLE                    possionratio;  /* Possion ratio */
     DOUBLE                    ALFAT;
     DOUBLE                    Sigy;
     DOUBLE                    Hard;
     DOUBLE                    GF;
     DOUBLE                    betah;
} PL_MISES;

/*----------------------------------------------------------------------*
 | plastic mises material including large strains      a.lipka 17/05    |
 *----------------------------------------------------------------------*/
typedef struct _PL_MISES_LS
{
     DOUBLE                    youngs;        /* Young's modulus */
     DOUBLE                    possionratio;  /* Possion ratio */
     DOUBLE                    ALFAT;
     DOUBLE                    Sigy;
     DOUBLE                    Hard;
     DOUBLE                    GF;
} PL_MISES_LS;

/*----------------------------------------------------------------------*
 | plastic foam  material                              a.lipka 02/12    |
 *----------------------------------------------------------------------*/
typedef struct _PL_FOAM
{
     DOUBLE                    youngs;        /* Young's modulus */
     DOUBLE                    possionratio;  /* Possion ratio */
     DOUBLE                    ALFAT;
     DOUBLE                    Sigy;
     DOUBLE                    Hard;
     DOUBLE                    GF;
} PL_FOAM;
/*----------------------------------------------------------------------*
 | plastic drucker prager material                     a.lipka 17/05    |
 *----------------------------------------------------------------------*/
typedef struct _PL_DP
{
     DOUBLE                    youngs;        /* Young's modulus */
     DOUBLE                    possionratio;  /* Possion ratio */
     DOUBLE                    ALFAT;
     DOUBLE                    Sigy;
     DOUBLE                    Hard;
     DOUBLE                    PHI;
} PL_DP;
/*----------------------------------------------------------------------*
 | elastoplastic concrete material                     a.lipka 17/05    |
 *----------------------------------------------------------------------*/
typedef struct _PL_EPC
{
     DOUBLE                    dens;
     /* concrete */
     DOUBLE                    youngs;       /* Young's modulus */
     DOUBLE                    possionratio; /* Possion ratio */
     DOUBLE                    alfat;
     DOUBLE                    sigy;
     DOUBLE                    phi;
     DOUBLE                    xsi;
     DOUBLE                    ftm;        /* tensile strength */
     DOUBLE                    gt;         /* tensile fracture energy */
     DOUBLE                    fcm;        /* compressive strength */
     DOUBLE                    gc;         /* compressive fracture energy */
     DOUBLE                    gamma1;     /* fitting factor yield function 1 */
     DOUBLE                    gamma2;     /* symm. biaxial compression stressfactor */
     DOUBLE                    dfac;       /* dammage factor: 0.0 plastic - 1.0 full damaged */
     /* tension stiffening */
     INT                       nstiff;     /* ==1 in consideration of tension stiffening */
     /* rebars */
     INT                       maxreb;     /* number of*/
     INT                      *rebar;      /* Id */
     DOUBLE                   *reb_area;   /* area   */
     DOUBLE                   *reb_ang;    /* angel  */
     DOUBLE                   *reb_so;     /* minimum bond length  */
     DOUBLE                   *reb_ds;     /* diameter  */
     DOUBLE                   *reb_rgamma; /* =4: deformed bars =2: plane par */
     DOUBLE                   *reb_dens; 
     DOUBLE                   *reb_alfat; 
     DOUBLE                   *reb_emod; 
     DOUBLE                   *reb_rebnue; 
     DOUBLE                   *reb_sigy; 
     DOUBLE                   *reb_hard;
} PL_EPC;
/*----------------------------------------------------------------------*
 | plastic mises porous material                       a.lipka 17/05    |
 *----------------------------------------------------------------------*/
typedef struct _PL_POR_MISES
{
     DOUBLE                    youngs;
     DOUBLE                    DP_YM;
     DOUBLE                    possionratio;
     DOUBLE                    ALFAT;
     DOUBLE                    Sigy;
     DOUBLE                    DP_Sigy;
     DOUBLE                    Hard;
     DOUBLE                    DP_Hard;
} PL_POR_MISES;
/*----------------------------------------------------------------------*
 | delamination material                               a.lipka 17/05    |
 *----------------------------------------------------------------------*/
typedef struct _PL_HASH
{
     DOUBLE                    emod1;
     DOUBLE                    emod2;
     DOUBLE                    emod3;
     DOUBLE                    xnue23;
     DOUBLE                    xnue13;
     DOUBLE                    xnue12;
     DOUBLE                    gmod12;
     DOUBLE                    gmod23;
     DOUBLE                    gmod13;
     DOUBLE                    s33;
     DOUBLE                    sinf33;
     DOUBLE                    s23;
     DOUBLE                    s13;
     DOUBLE                    gamma;
     DOUBLE                    gc;
     DOUBLE                    deltat;
     DOUBLE                    eta_i;
     DOUBLE                    c1hgt;
     DOUBLE                    c1layhgt;
     INT                       ivisco;
} PL_HASH;

/*----------------------------------------------------------------------*
 | elastic orthotropic material                              al 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _EL_ORTH
{
     DOUBLE                    emod1;
     DOUBLE                    emod2;
     DOUBLE                    emod3;
     DOUBLE                    xnue23;
     DOUBLE                    xnue13;
     DOUBLE                    xnue12;
     DOUBLE                    gmod12;
     DOUBLE                    gmod23;
     DOUBLE                    gmod13;
} EL_ORTH;

/*----------------------------------------------------------------------*
 | open cell metal foam material (optimization)              al 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _MFOC
{
     DOUBLE                    es;             /* Young's modulus (cell) */
     DOUBLE                    pr;             /* Possion ratio */
     DOUBLE                    dens;           /* density foam  */
     DOUBLE                    denss;          /* density (bulk) */
     DOUBLE                    denmin;         /* min. dens. foam (opti.)*/
     DOUBLE                    denmax;         /* max. dens. foam (opti.)*/
     DOUBLE                    refdens;        /* reference density */
     DOUBLE                    oce;            /* exponent  */
     DOUBLE                    ocf;            /* factor    */
} MFOC;

/*----------------------------------------------------------------------*
 | closed cell metal foam material (optimization)            al 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _MFCC
{
     DOUBLE                    es;             /* Young's modulus (cell) */
     DOUBLE                    pr;             /* Possion ratio */
     DOUBLE                    dens;           /* density foam  */
     DOUBLE                    denss;          /* density (bulk) */
     DOUBLE                    denmin;         /* min. dens. foam (opti.)*/
     DOUBLE                    denmax;         /* max. dens. foam (opti.)*/
     DOUBLE                    refdens;        /* reference density */
     DOUBLE                    cce;            /* exponent  */
     DOUBLE                    ccf;            /* factor    */
} MFCC;
/*----------------------------------------------------------------------*
 | foam, closed cell, based on modified Neo Hook             al 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _NHMFCC
{
     DOUBLE                    es;             /* Young's modulus (cell) */
     DOUBLE                    pr;             /* Possion ratio */
     DOUBLE                    dens;           /* density foam  */
     DOUBLE                    denss;          /* density (bulk) */
     DOUBLE                    denmin;         /* min. dens. foam (opti.)*/
     DOUBLE                    denmax;         /* max. dens. foam (opti.)*/
     DOUBLE                    refdens;        /* reference density */
     DOUBLE                    cce;            /* exponent  */
     DOUBLE                    ccf;            /* factor    */
} NHMFCC;
/*----------------------------------------------------------------------*
 | multi layer material  -> shell9                          sh 10/02    |
 | material, that can have differnt materials in different layers       |
 | in here is just the definition of the cross section                  |
 *----------------------------------------------------------------------*/
typedef struct _MULTI_LAYER
{
     INT                       num_klay;      /* number of kinematic layers */
     DOUBLE                   *klayhgt;       /* hgt of a kinematic layer in % of total thickness of the shell */
     struct _KINLAY           *kinlay;        /* one kinematic layer */

} MULTI_LAYER;
/*----------------------------------------------------------------------*
 | for multi layer material  -> shell9                      sh 10/02    |
 | information about one kinematic layer                                |
 *----------------------------------------------------------------------*/
typedef struct _KINLAY
{
     INT                       num_mlay;     /* number of material layer to this kinematic layer*/
     DOUBLE                   *mlayhgt;      /* hgt of a material layer in % of thickness of the adjacent kinematic layer*/ 
     INT                      *mmatID;       /* ID of multilayer material in every material layer */
     DOUBLE                   *phi;          /* rotation of the material in one material layer */
     INT                      *rot;          /* axis of rotation of the material: x=1, y=2, z=3*/
} KINLAY;
/*----------------------------------------------------------------------*
 | multilayer materials                                     sh 10/02    |
 | structure to hold all types of material laws                         |
 | is equivalent to the struct MATERIAL but is used for shell9 only     |
 *----------------------------------------------------------------------*/
typedef struct _MULTIMAT
{
     INT                       Id;           /* Id of the material */

     enum _MATERIAL_TYP        mattyp;       /* type of material */

     union
     {
     struct _STVENANT         *stvenant;     /* St. Venant-Kirchhoff material */
     struct _NEO_HOOKE        *neohooke;     /* Neo-Hooke material */
     struct _ORTHOTROPIC      *orthotropic;  /* linear elastic, orthotropic material*/
     }                         m;            /* union pointer to material specific structure */

} MULTIMAT;
/*----------------------------------------------------------------------*
 | linear elastic, orthotropic material                     sh 02/03    |
 *----------------------------------------------------------------------*/
typedef struct _ORTHOTROPIC
{
     DOUBLE                    emod1;         
     DOUBLE                    emod2;         
     DOUBLE                    emod3;         
     DOUBLE                    gmod12;         
     DOUBLE                    gmod13;         
     DOUBLE                    gmod23;         
     DOUBLE                    xnue12;         
     DOUBLE                    xnue13;         
     DOUBLE                    xnue23;         
} ORTHOTROPIC;
