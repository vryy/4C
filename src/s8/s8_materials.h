#ifndef S8_MATERIALS_H
#define S8_MATERIALS_H

/*---------------------------------------------------------------------*/
/*! \file

\brief C-style material definitions

\level 3


*/
/*---------------------------------------------------------------------*/
typedef enum _MATERIAL_TYP
{
  m_stvenant,  /* St.Venant Kirchhoff material */
  m_neohooke,  /* Neo-Hooke material */
  m_viscohyper /* compressible viscous Ogden material (with shell8) */
} MATERIAL_TYP;

/*----------------------------------------------------------------------*
 | materials                                              m.gee 4/01    |
 | structure to hold all types of material laws                         |
 *----------------------------------------------------------------------*/
typedef struct _MATERIAL
{
  int Id; /* Id of the material */

  enum _MATERIAL_TYP mattyp; /* type of material */

  union
  {
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
  double youngs;       /* Young's modulus */
  double possionratio; /* Possion ratio */
  double density;      /* material specific weight */
  double thermexpans;  /* coefficient of thermal expansion */
} STVENANT;


/*----------------------------------------------------------------------*
 | Neo Hooke material                                     m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _NEO_HOOKE
{
  double youngs;       /* Young's modulus */
  double possionratio; /* Possion ratio */
  double density;      /* material specific weight */
} NEO_HOOKE;


/*----------------------------------------------------------------------*
 | compressible ogden material                            m.gee 6/03    |
 *----------------------------------------------------------------------*/
typedef struct _COMPOGDEN
{
  int init;        /* init flag */
  double nue;      /* Possion ratio */
  double beta;     /* the unphysical material constant called beta */
  double alfap[3]; /* three parameters alfap */
  double mup[3];   /* three parameters nuep */
  double density;  /* material specific weight */
  double lambda;   /* 1. lame constant */
  double kappa;    /* bulkmodulus */
  double l[3];
} COMPOGDEN;

/*----------------------------------------------------------------------*
 | viscoelastic compressible ogden material               m.gee 9/03    |
 *----------------------------------------------------------------------*/
typedef struct _VISCOHYPER
{
  int init;        /* init flag */
  double nue;      /* Possion ratio */
  double beta;     /* the unphysical material constant called beta */
  double alfap[3]; /* three parameters alfap */
  double mup[3];   /* three parameters nuep */
  double density;  /* material specific weight */
  double lambda;   /* 1. lame constant */
  double kappa;    /* bulkmodulus */
  int nmaxw;       /* number of maxwell elements in the material (1-4) */
  double tau[4];   /* relaxation times of hte maxwell elements */
  double betas[4]; /* strain energy factors of the springs of the maxwell elements */
} VISCOHYPER;

#endif
