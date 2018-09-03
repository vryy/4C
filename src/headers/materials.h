/*!---------------------------------------------------------------------
\file materials.h
\brief C-style material definitions
\level 3
\maintainer Martin Kronbichler
---------------------------------------------------------------------*/
#if defined(D_SHELL8)
typedef enum _MATERIAL_TYP
{
  m_stvenant,  /* St.Venant Kirchhoff material */
  m_neohooke,  /* Neo-Hooke material */
  m_compogden, /* compressible Ogden material (with shell8) */
  m_viscohyper /* compressible viscous Ogden material (with shell8) */
} MATERIAL_TYP;

/*----------------------------------------------------------------------*
 | materials                                              m.gee 4/01    |
 | structure to hold all types of material laws                         |
 *----------------------------------------------------------------------*/
typedef struct _MATERIAL
{
  INT Id; /* Id of the material */

  enum _MATERIAL_TYP mattyp; /* type of material */

  union {
    struct _STVENANT *stvenant;     /* St. Venant-Kirchhoff material */
    struct _NEO_HOOKE *neohooke;    /* Neo-Hooke material */
    struct _COMPOGDEN *compogden;   /* compressible ogden hyperelastic material */
    struct _VISCOHYPER *viscohyper; /* viscoelastic compressible ogden hyperelastic material */
  } m;                              /* union pointer to material specific structure */

} MATERIAL;

/*----------------------------------------------------------------------*
 | St. Venant-Kirchhoff material                          m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _STVENANT
{
  DOUBLE youngs;       /* Young's modulus */
  DOUBLE possionratio; /* Possion ratio */
  DOUBLE density;      /* material specific weight */
  DOUBLE thermexpans;  /* coefficient of thermal expansion */
} STVENANT;


/*----------------------------------------------------------------------*
 | Neo Hooke material                                     m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _NEO_HOOKE
{
  DOUBLE youngs;       /* Young's modulus */
  DOUBLE possionratio; /* Possion ratio */
  DOUBLE density;      /* material specific weight */
} NEO_HOOKE;


/*----------------------------------------------------------------------*
 | compressible ogden material                            m.gee 6/03    |
 *----------------------------------------------------------------------*/
typedef struct _COMPOGDEN
{
  INT init;        /* init flag */
  DOUBLE nue;      /* Possion ratio */
  DOUBLE beta;     /* the unphysical material constant called beta */
  DOUBLE alfap[3]; /* three parameters alfap */
  DOUBLE mup[3];   /* three parameters nuep */
  DOUBLE density;  /* material specific weight */
  DOUBLE lambda;   /* 1. lame constant */
  DOUBLE kappa;    /* bulkmodulus */
#if 1
  DOUBLE l[3];
#endif
} COMPOGDEN;

/*----------------------------------------------------------------------*
 | viscoelastic compressible ogden material               m.gee 9/03    |
 *----------------------------------------------------------------------*/
typedef struct _VISCOHYPER
{
  INT init;        /* init flag */
  DOUBLE nue;      /* Possion ratio */
  DOUBLE beta;     /* the unphysical material constant called beta */
  DOUBLE alfap[3]; /* three parameters alfap */
  DOUBLE mup[3];   /* three parameters nuep */
  DOUBLE density;  /* material specific weight */
  DOUBLE lambda;   /* 1. lame constant */
  DOUBLE kappa;    /* bulkmodulus */
  INT nmaxw;       /* number of maxwell elements in the material (1-4) */
  DOUBLE tau[4];   /* relaxation times of hte maxwell elements */
  DOUBLE betas[4]; /* strain energy factors of the springs of the maxwell elements */
} VISCOHYPER;

#endif /* !defined(D_SHELL8) */
