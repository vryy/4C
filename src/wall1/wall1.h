/*!----------------------------------------------------------------------
\file
\brief headerfile for wall element, containing structures and prototypes

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_WALL1

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief 2D ale element data

<pre>                                                              al 06/01
This structure contains the coordinates and weights for numerical
integration of a wall element
</pre>

 *----------------------------------------------------------------------*
 | wall1 data                                                al 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _W1_DATA
{
DOUBLE        xgrr[13];         /*!< natural coordinates r of gaussian points */
DOUBLE        wgtr[13];         /*!< weights at natural coordinates r */

DOUBLE        xgss[13];         /*!< natural coordinates s of gaussian points */
DOUBLE        wgts[13];         /*!< weights at natural coordinates s */
DOUBLE        txgr[MAXTINTP][MAXTINTC];   /*!< coordinates in r for TRIS and TETS */
DOUBLE        txgs[MAXTINTP][MAXTINTC];   /*!< coordinates in s for TRIS and TETS*/
DOUBLE        twgt[MAXTINTP][MAXTINTC];   /*!< weights for TRIS and TETS*/
DOUBLE        qxg[MAXQINTP][MAXQINTC];    /*!< coordinates for QUADS and HEX */
DOUBLE        qwgt[MAXQINTP][MAXQINTC];   /*!< weights for QUADS and HEX */
} W1_DATA;
/*----------------------------------------------------------------------*
 | working array                                             al 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _W1_IMODE_WA
{
DOUBLE      knninv[4][4]; /* inverse of stiffness of incomp modes   */
DOUBLE      knc[8][4];    /* mixed stiffness compat. - incompatible */
DOUBLE      fintn[4];     /* internal forces of incomp modes        */
DOUBLE      alpha[4];     /* INT.dof for inc modes of last loadst k */
} W1_IMODE_WA;

typedef struct _W1_IP_WA
{
DOUBLE      sig[4]; /* global stresses                              */
DOUBLE      t_d[2]; /* global discrete stresses                     */
DOUBLE      sig_esz[4]; /* global stresses for plain stress         */
DOUBLE      eps[4]; /* global strains                               */
DOUBLE      eps_esz[4]; /* global strains for plain stress          */
DOUBLE      d4_esz[4];  /* tangentstiffness (x,3) for plain stress  */
DOUBLE         *qn; /* 'backstress' vec                             */
DOUBLE       epstn; /* equivalent strain                            */
INT          yip;   /* stress state: 1=elastic 2=plastic            */
DOUBLE       kap;   /* history variable for damage-model            */
DOUBLE       dam;   /* damage variable for output                   */
/* concrete */
DOUBLE       *sigc; /* condensed global stresses                    */
DOUBLE       *grad; /* actual gradient                              */
DOUBLE       *dlam; /* uniaxial plastic strain  (- dlam)            */
DOUBLE       *sigi; /* stresses from last iteration step            */
DOUBLE       *epsi; /* strains from last iteration step             */
DOUBLE       *di  ; /* components d41,d42,d43,d44 of the            */
                    /* constitutive tensor from last iteration step */
/* epc3D   ...*/  /*SH 12/03*/
DOUBLE   kappa_c;   /*!< internal variable -> controls evolution in compression*/
DOUBLE   kappa_t;   /*!< internal variable -> controls evolution in tension*/
                    /* constitutive tensor from last iteration step */
/* rebar */
DOUBLE       *rsig;   /* stress */
DOUBLE       *reps;   /* strain */
DOUBLE       *repstn; /* equivalent strain */
INT          *ryip;   /* flag */
/* tension stiffening */
DOUBLE       stcappae;  /* max. elastic strain of the concrete     */
DOUBLE       stcappaut; /* fracture strain of the concrete         */
DOUBLE       stalpha;   /* angle of the principal concrete stress  */
DOUBLE       stthick ;  /* thickness ot the structure              */
DOUBLE       stfbd;     /* tension stiffening stress               */
/* incompatible mode method */
DOUBLE       D[3];      /* material tangente:D0=C11,D1=C22,D2=C12=C21 */
/* damgage */
DOUBLE       kappa;     /* max. equival. strain ever reached*/
DOUBLE       aequistrain;     /* actual equival. strain */
DOUBLE       damage;    /* damage value*/
} W1_IP_WA;
typedef struct _W1_ELE_WA
{
W1_IP_WA      *ipwa;
W1_IMODE_WA   *imodewa;
DOUBLE         dia;
DOUBLE        *matdata; /* element material data, actual density ... */
INT           *optdata; /* optimization variable number ... */
} W1_ELE_WA;

/*----------------------------------------------------------------------*
 | type of 2D problem                                       al 01/02    |
 *----------------------------------------------------------------------*/
typedef enum _WALL_TYPE
{
                       plane_strain,
                       plane_stress,
                       rotat_symmet
} WALL_TYPE;
/*----------------------------------------------------------------------*
 | type of Kinematik                                       ah 06/02    |
 *----------------------------------------------------------------------*/
typedef enum _KINEMATIK_TYPE
{
                       geo_lin,
                       total_lagr,
                       updated_lagr
} KINEMATIK_TYPE;
/*----------------------------------------------------------------------*
 | type of model (displacem. or mixed) for QUAD4            ah 08/02    |
 *----------------------------------------------------------------------*/
typedef enum _MODEL_TYPE
{
                       displ_model,   /* normal QUAD4, diplacem. model  */
                       incomp_mode    /* Q1 with incompatible modes     */
} MODEL_TYPE;
/*----------------------------------------------------------------------*
 | wall                                                     al 01/02    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | result struct: stresses                                   fh 6/02    |
 *----------------------------------------------------------------------*/
typedef struct _WALL1
{
enum _WALL_TYPE    wtype;             /* type of 2D problem, see above */
enum _KINEMATIK_TYPE  kintype;	      /* type of Kinematik             */
enum _MODEL_TYPE modeltype;	      /* type of model for QUAD4       */

INT           nGP[4];

DOUBLE        thick;

#ifdef GEMM
DOUBLE strain_energy;
DOUBLE kinetic_energy;
DOUBLE angular_momentum;
DOUBLE linmom[2];
#endif

W1_ELE_WA     *elewa;                         /* element working array */
#ifdef D_MLSTRUCT
  INT            isinomegaprime;     /* flag, if element is in Omega' */
  INT            firstinomegaprime;  /* flag, if element is first time in Omega' */
  DOUBLE         translation[2];     /* translation vector from submesh origin to element origin */
  struct  _SM_ELEMENT_DATA    *sm_eledata;
  struct  _SM_NODAL_DATA      *sm_nodaldata;
  struct  _DBCSR              *stiff_mi_ma_csr; /* stored stiffness for rueckrechnung to d' */
  union   _SPARSE_ARRAY       *stiff_mi_mi_ccf; /* stored stiffness for rueckrechnung to d' */
  DOUBLE   fint_mi[3000];                       /* stored int.forces for rueckrechnung to d' */
#endif /* D_MLSTRUCT */

/*---------------------------------------------- array of forces ------*/
enum
    {
    w1_xy,
    w1_rs
    }            stresstyp;
#define WALL1_STRESSTYPE { "w1_xy", "w1_rs", NULL }

struct _ARRAY4D  stress_GP;
struct _ARRAY4D  stress_ND;
/*---------------------------------------------------------------------*/
/*------------------------------------------History for  EM Int. Scheme*/
#ifdef GEMM
struct _ARRAY4D b_bar_history;
struct _ARRAY4D PK_history;
#endif
#ifdef D_SSI
   enum _SSI_COUPTYP         ssi_couptyp;
#endif
#ifdef D_FSI
   enum _FSI_COUPTYP         fsi_couptyp;
#endif
#ifdef D_TSI
  enum _TSI_COUPTYP tsi_couptyp;  /* TSI coupling type */
  struct _ELEMENT *therm_ele;      /* pointer to conforming thermal2 element */
#endif
} WALL1;
/*----------------------------------------------------------------------*
 |  prototypes of main and input routines                   al 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | main wall1  control routine                               al 6/01    |
 *----------------------------------------------------------------------*/
void wall1(PARTITION   *actpart,
           INTRA       *actintra,
           ELEMENT     *ele,
           ARRAY       *estif_global,
           ARRAY       *emass_global,
           ARRAY       *intforce_global,
           CALC_ACTION *action,
           CONTAINER   *container);   /* contains variables defined in container.h */

/*----------------------------------------------------------------------*
 | read wall element                                         al 9/01    |
 *----------------------------------------------------------------------*/
void w1inp(ELEMENT *ele);
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
#endif
