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
