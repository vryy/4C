/*----------------------------------------------------------------------*
 | materials                                              m.gee 4/01    |
 | structure to hold all types of material laws                         |
 *----------------------------------------------------------------------*/
typedef struct _MATERIAL
{
     int                       Id;           /* Id of the material */

     enum _MATERIAL_TYP        mattyp;       /* type of material */

     union
     {
     struct _STVENANT         *stvenant;     /* St. Venant-Kirchhoff material */
     struct _PL_MISES         *pl_mises;     /* von Mises material */
     struct _PL_MISES_LS      *pl_mises_ls;  /* von Mises material - large strains*/
     struct _PL_DP            *pl_dp;        /* Drucker Prager material */
     struct _PL_EPC           *pl_epc;       /* elastoplastic concrete material */
     struct _STVENPOR         *stvenpor;     /* porous St. Ven.-Kirch. material */
     struct _PL_POR_MISES     *pl_por_mises; /* porous von Mises material */
     struct _NEO_HOOKE        *neohooke;     /* Neo-Hooke material */
     struct _FLUID            *fluid;        /* fluid material */
     struct _PL_HASH          *pl_hash;      /* elpl. hashin delamination material */
     struct _EL_ORTH          *el_orth;      /* elastic orthotropic material */
     struct _MFOC             *mfoc;         /* metal foam, open cell  */
     struct _MFCC             *mfcc;         /* metal foam, closed cell  */
     }                         m;            /* union pointer to material specific structure */

} MATERIAL;





/*----------------------------------------------------------------------*
 | St. Venant-Kirchhoff material                          m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _STVENANT
{
     double                    youngs;         /* Young's modulus */
     double                    possionratio;   /* Possion ratio */
     double                    density;        /* material specific weight */
} STVENANT;


/*----------------------------------------------------------------------*
 | porous St. Venant-Kirchhoff material (optimization)       al 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _STVENPOR
{
     double                    youngs;         /* Young's modulus */
     double                    possionratio;   /* Possion ratio */
     double                    density;        /* material specific weight */
     double                    refdens;        /* reference density */
     double                    exponent;       /* material parameter */
} STVENPOR;


/*----------------------------------------------------------------------*
 | Neo Hooke material                                     m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _NEO_HOOKE
{
     double                    youngs;         /* Young's modulus */
     double                    possionratio;   /* Possion ratio */
     double                    density;        /* material specific weight */
} NEO_HOOKE;



/*----------------------------------------------------------------------*
 | fluid material                                         m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _FLUID
{
     double                    viscosity;
     double                    density;
} FLUID;
/*----------------------------------------------------------------------*
 | plastic mises material                              a.lipka 17/05    |
 *----------------------------------------------------------------------*/
typedef struct _PL_MISES
{
     double                    youngs;        /* Young's modulus */
     double                    possionratio;  /* Possion ratio */
     double                    ALFAT;
     double                    Sigy;
     double                    Hard;
     double                    GF;
} PL_MISES;

/*----------------------------------------------------------------------*
 | plastic mises material including large strains      a.lipka 17/05    |
 *----------------------------------------------------------------------*/
typedef struct _PL_MISES_LS
{
     double                    youngs;        /* Young's modulus */
     double                    possionratio;  /* Possion ratio */
     double                    ALFAT;
     double                    Sigy;
     double                    Hard;
     double                    GF;
} PL_MISES_LS;

/*----------------------------------------------------------------------*
 | plastic drucker prager material                     a.lipka 17/05    |
 *----------------------------------------------------------------------*/
typedef struct _PL_DP
{
     double                    youngs;        /* Young's modulus */
     double                    possionratio;  /* Possion ratio */
     double                    ALFAT;
     double                    Sigy;
     double                    Hard;
     double                    PHI;
} PL_DP;
/*----------------------------------------------------------------------*
 | elastoplastic concrete material                     a.lipka 17/05    |
 *----------------------------------------------------------------------*/
typedef struct _PL_EPC
{
     double                    dens;
     /* concrete */
     double                    youngs;       /* Young's modulus */
     double                    possionratio; /* Possion ratio */
     double                    alfat;
     double                    sigy;
     double                    phi;
     double                    xsi;
     double                    ftm;        /* tensile strength */
     double                    gt;         /* tensile fracture energy */
     double                    fcm;        /* compressive strength */
     double                    gc;         /* compressive fracture energy */
     double                    gamma1;     /* fitting factor yield function 1 */
     double                    gamma2;     /* symm. biaxial compression stressfactor */
     double                    dfac;       /* dammage factor: 0.0 plastic - 1.0 full damaged */
     /* tension stiffening */
     int                       nstiff;     /* ==1 in consideration of tension stiffening */
     /* rebars */
     int                       maxreb;     /* number of*/
     int                      *rebar;      /* Id */
     double                   *reb_area;   /* area   */
     double                   *reb_ang;    /* angel  */
     double                   *reb_so;     /* minimum bond length  */
     double                   *reb_ds;     /* diameter  */
     double                   *reb_rgamma; /* =4: deformed bars =2: plane par */
     double                   *reb_dens; 
     double                   *reb_alfat; 
     double                   *reb_emod; 
     double                   *reb_rebnue; 
     double                   *reb_sigy; 
     double                   *reb_hard;
} PL_EPC;
/*----------------------------------------------------------------------*
 | plastic mises porous material                       a.lipka 17/05    |
 *----------------------------------------------------------------------*/
typedef struct _PL_POR_MISES
{
     double                    youngs;
     double                    DP_YM;
     double                    possionratio;
     double                    ALFAT;
     double                    Sigy;
     double                    DP_Sigy;
     double                    Hard;
     double                    DP_Hard;
} PL_POR_MISES;
/*----------------------------------------------------------------------*
 | delamination material                               a.lipka 17/05    |
 *----------------------------------------------------------------------*/
typedef struct _PL_HASH
{
     double                    emod1;
     double                    emod2;
     double                    emod3;
     double                    xnue23;
     double                    xnue13;
     double                    xnue12;
     double                    gmod12;
     double                    gmod23;
     double                    gmod13;
     double                    s33;
     double                    sinf33;
     double                    s23;
     double                    s13;
     double                    gamma;
     double                    gc;
     double                    deltat;
     double                    eta_i;
     double                    c1hgt;
     double                    c1layhgt;
     int                       ivisco;
} PL_HASH;

/*----------------------------------------------------------------------*
 | elastic orthotropic material                              al 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _EL_ORTH
{
     double                    emod1;
     double                    emod2;
     double                    emod3;
     double                    xnue23;
     double                    xnue13;
     double                    xnue12;
     double                    gmod12;
     double                    gmod23;
     double                    gmod13;
} EL_ORTH;

/*----------------------------------------------------------------------*
 | open cell metal foam material (optimization)              al 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _MFOC
{
     double                    es;             /* Young's modulus (cell) */
     double                    pr;             /* Possion ratio */
     double                    dens;           /* density foam  */
     double                    denss;          /* density (bulk) */
     double                    denmin;         /* min. dens. foam (opti.)*/
     double                    denmax;         /* max. dens. foam (opti.)*/
     double                    refdens;        /* reference density */
     double                    oce;            /* exponent  */
     double                    ocf;            /* factor    */
} MFOC;

/*----------------------------------------------------------------------*
 | closed cell metal foam material (optimization)            al 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _MFCC
{
     double                    es;             /* Young's modulus (cell) */
     double                    pr;             /* Possion ratio */
     double                    dens;           /* density foam  */
     double                    denss;          /* density (bulk) */
     double                    denmin;         /* min. dens. foam (opti.)*/
     double                    denmax;         /* max. dens. foam (opti.)*/
     double                    refdens;        /* reference density */
     double                    cce;            /* exponent  */
     double                    ccf;            /* factor    */
} MFCC;

