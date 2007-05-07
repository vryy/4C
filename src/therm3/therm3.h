/*======================================================================*/
/*!
\file
\brief headerfile for 3dim thermal element (THERM3), containing 
structures and prototypes

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/
/*======================================================================*/

#ifdef D_THERM3

/*!
\addtogroup THERM3
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "../output/gid.h"

/*======================================================================*/
/* define constants of THERM3 */
#define NDIM_THERM3      (3)    /* 3dim problem */

#ifndef MAXNOD_THERM3
#ifdef MAXNOD
#define MAXNOD_THERM3    (MAXNOD)
#else
#define MAXNOD_THERM3    (27)   /* maximum of element nodes : max hex27 */
#endif
#endif

#ifndef MAXEDG_THERM3
#define MAXEDG_THERM3    (12)   /* maximal number element edges */
#endif

#ifndef MAXSID_THERM3
#define MAXSID_THERM3    (6)    /* maximal number element sides */
#endif

#ifndef MAXNE_THERM3
#define MAXNE_THERM3     (3)    /* maximal number of nodes on an edge */
#endif

#ifndef MAXNS_THERM3
#define MAXNS_THERM3     (9)    /* maximal number of nodes on a side */
#endif

#define DIMSID_THERM3    (2)    /* dimension of element sides */

#define NUMDOF_THERM3    (1)    /* number of thermal DOFs at each node :
                                 * temperature */

#define NUMTMGR_THERM3   (3)    /* number of temperature gradients */

#define NUMHFLX_THERM3   (3)    /* number of heat fluxes : q_x, q_y, q_z */

#ifndef GLINTC_THERM3
#define GLINTC_THERM3    (6)    /* line Gauss integration cases */
#endif

#ifndef GLMAXP_THERM3
#define GLMAXP_THERM3    (6)    /* line domain max. number of Gauss points */
#endif

#ifndef GTINTC_THERM3 
#define GTINTC_THERM3    (3)    /* triangle domain Gauss integration cases */
#endif

#ifndef GTMAXP_THERM3
#define GTMAXP_THERM3    (5)    /* max. number of Gauss points tet domain */
#endif

#ifndef GSINTC_THERM3
#define GSINTC_THERM3    (5)    /* tet sides Gauss integration cases */
#endif

#ifndef GSMAXP_THERM3
#define GSMAXP_THERM3    (6)    /* tet sides max. number of Gauss points */
#endif

#ifdef DEBUG
#define TEST_THERM3             /* compile test routines if ... */
#undef TEST_THERM3              /* ... not undefined */
#else
#endif

/*======================================================================*/
/* global declarations, variables etc */

/*----------------------------------------------------------------------*/
/*!
\brief Gauss points and weights

The constants are defined in headers/define_sizes.h

\author bborn
\date 09/06
*/
typedef struct _TH3_DATA
{
  /*--------------------------------------------------------------------*/
  /* Gauss coordinates and weights */
  /* hexahedron domain, sides and edges --> line [-1,+1] */
  DOUBLE ghlc[GLINTC_THERM3][GLMAXP_THERM3];  /* coordinates */
  DOUBLE ghlw[GLINTC_THERM3][GLMAXP_THERM3];  /* weights */
  /* tetrahedron domain [T.J.R. Hughes, "The FEM", Dover 2000] */
  DOUBLE gtdc[GTINTC_THERM3][GTMAXP_THERM3][NDIM_THERM3];  /* coordinates 
                                                            * in r,s,t */
  DOUBLE gtdw[GTINTC_THERM3][GTMAXP_THERM3];  /* weights */
  /* tetrahedron sides */
  DOUBLE gtsc[GSINTC_THERM3][GSMAXP_THERM3][DIMSID_THERM3];  /* coordinates 
                                                                in side */
  DOUBLE gtsw[GSINTC_THERM3][GSMAXP_THERM3];  /* weights */
  /* triangle edges --> line [0,+1] */
  DOUBLE gtlc[GLINTC_THERM3][GLMAXP_THERM3];  /* coordinates */
  DOUBLE gtlw[GLINTC_THERM3][GLMAXP_THERM3];  /* weights */
  /*--------------------------------------------------------------------*/
  /* numbering of element nodes, edges, sides in parameter space */
  /* parameter coordinates of nodes */
  DOUBLE nodhrst[MAXNOD_THERM3][NDIM_THERM3];  /* hexahedron */
  DOUBLE nodtrst[MAXNOD_THERM3][NDIM_THERM3];  /* tetrahedron */
  /* nodes on sides (surfaces) */
  INT nodsidh[MAXSID_THERM3][MAXNS_THERM3];  /* hexahedra */
  INT nodsidt[MAXSID_THERM3][MAXNS_THERM3];  /* tetrahedra */
  /* nodes on edges */
  INT nodedghl[MAXEDG_THERM3][MAXNE_THERM3];  /* linear hex8 */       
  INT nodedghq[MAXEDG_THERM3][MAXNE_THERM3];  /* quadratic hex20,27 */
  INT nodedgtl[MAXEDG_THERM3][MAXNE_THERM3];  /* linear tet4 */       
  INT nodedgtq[MAXEDG_THERM3][MAXNE_THERM3];  /* quadratic tet10 */
  /*--------------------------------------------------------------------*/
  /* anchor and span vectors for sides and edges in param. space */
  /* sides hex */
  DOUBLE ancsidh[MAXSID_THERM3][NDIM_THERM3];  /* anchors hex */
  DOUBLE redsidh[MAXSID_THERM3][DIMSID_THERM3][NDIM_THERM3];  /* dim red 
                                                               * matrix */
  /* sides tet */
  DOUBLE ancsidt[MAXSID_THERM3][NDIM_THERM3];  /* anchors tet */
  DOUBLE redsidt[MAXSID_THERM3][DIMSID_THERM3][NDIM_THERM3];  /* dim red 
                                                               * matrix */
  /* edges hex */
  DOUBLE ancedgh[MAXEDG_THERM3][NDIM_THERM3];  /* anchors hex */
  DOUBLE rededgh[MAXEDG_THERM3][NDIM_THERM3];  /* dimension reduction
                                                * matrix multiplied
                                                * on Jacobi matrix */
  /* edges tet */
  DOUBLE ancedgt[MAXEDG_THERM3][NDIM_THERM3];  /* anchors tet */
  DOUBLE rededgt[MAXEDG_THERM3][NDIM_THERM3];  /* dimension reduction
                                                * matrix multiplied
                                                * on Jacobi matrix */
} TH3_DATA;


/*----------------------------------------------------------------------*/
/*!
\brief Type of kinematics

\author bborn
\date 09/06
*/
typedef enum _TH3_KINEMATICS
{
  th3_geo_lin,  /* geometrically linear */
  th3_total_lagr,  /* geometrically non-linear => total Lagrangian */
  th3_updated_lagr  /* geometrically non-linear => updated Lagrangian */
} TH3_KINEMATICS;


/*----------------------------------------------------------------------*/
/*!
\brief Type of heat flux output

\author bborn
\date 09/06
*/
typedef enum _TH3_HEATFLUXOUT
{
  th3_hflux_none,
  th3_hflux_gpxyz,  /* globally xyz-oriented at Gauss points */
  th3_hflux_gp123,  /* modulus and its direction */
  th3_hflux_gprst,  /* locally rs-oriented at Gauss points */
  th3_hflux_gpxyz123,  /* all of the above */
  th3_hflux_ndxyz,
  th3_hflux_nd123,
  th3_hflux_ndrst,
  th3_hflux_ndxyz123
} TH3_HEATFLUXOUT;

/*----------------------------------------------------------------------*/
/*!
\brief Definition of THERM3 type holding THERM3 element properties

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
typedef struct _THERM3
{
  TH3_KINEMATICS kintype;  /* type of kinematics */
  TH3_HEATFLUXOUT hfluxtype;  /* output type of heat flux */

  /* number of Gauss points as obtained at read-in
   * hexahedra:
   *    gpnum[0] : in r-direction : 1,2,3,4,5,6 : read-in
   *    gpnum[1] : in s-direction : 1,2,3,4,5,6 : read-in/set
   *    gpnum[2] : in t-direction : 1,2,3,4,5,6 : read-in/set
   * tetrahedra:
   *    gpnum[0] : total number of GPs in domain: 1,4,5 : read-in
   *    gpnum[1] : total number of GPs on sides: 1,3,4,6 : read-in/set
   *    gpnum[2] : total number of GPs on edges: 1,2,3,4,5,6 : read-in/set */
  INT gpnum[NDIM_THERM3];
  /* Gauss integration case
   * hexahedra:
   *    intc[i] = gpnum[0] - 1
   * tetrahedra:
   *    ... */
  INT gpintc[NDIM_THERM3];
  /* total number of GPs in domain */
  INT gptot;

  /* heat flux vector at Gauss points or nodes
   *    hflux_gp_xyz : at Gauss points in global XYZ-components
   *    hflux_gp_rst : at Gauss points in parameter rst-compoonents
   *    hflux_gp_123 : at Gauss points modulus and angles vs XYZ-direct. */
  ARRAY4D hflux_gp_xyz;
  ARRAY4D hflux_gp_rst;
  ARRAY4D hflux_gp_123;
  ARRAY4D hflux_nd_xyz;
  ARRAY4D hflux_nd_rst;
  ARRAY4D hflux_nd_123;

#ifdef D_TSI
  TSI_COUPTYP tsi_couptyp;
  ELEMENT *struct_ele;
#endif
} THERM3;



/*======================================================================*/
/* Declarations of functions in therm2 directory
 * Order: Firstly, alphabetically list file names, secondly, list
 *        functions in file according to appearance */


/*----------------------------------------------------------------------*/
/* file th3_bop.c */
void th3_bop(INT        enod,
             DOUBLE     deriv[MAXNOD_THERM3][NDIM_THERM3],
             DOUBLE     xji[NDIM_THERM3][NDIM_THERM3],
             DOUBLE     bop[NUMTMGR_THERM3][NUMDOF_THERM3*MAXNOD_THERM3]);

/*---------------------------------------------------------------------*/
/* file th3_cfg.c */
void th3_cfg_chkdef();
void th3_cfg_init(TH3_DATA *data);
void th3_cfg_noderst(ELEMENT *ele,
                     TH3_DATA *data,
                     INT inode,
                     DOUBLE *rst);
#ifdef TEST_THERM3
void th3_cfg_test(TH3_DATA *data);
#endif

/*----------------------------------------------------------------------*/
/* file th3_gid.c */
void th3_gid_init(ELEMENT *actele,
                  GIDSET *actgid);
void th3_gid_msh(FIELD *actfield,
                 INT ldis,
                 GIDSET *actgid,
                 FILE *out);
void th3_gid_gpset(INT jdis,
                   GIDSET *actgid,
                   FILE *out);
void th3_gid_dom(FIELD *actfield,
                 INT disnum,
                 GIDSET *actgid,
                 FILE *out);
void th3_gid_hflux(char resstring[],
                   FIELD *actfield,
                   INT disnum,
                   INT step,
                   GIDSET *actgid,
                   FILE *out);

/*----------------------------------------------------------------------*/
/* file th3_hflux.c */
void th3_hflux_init(PARTITION *actpart);
void th3_hflux_final(PARTITION *actpart);
void th3_hflux_cal(CONTAINER *cont,
                   ELEMENT *ele,
                   TH3_DATA *data,
                   MATERIAL *mat);
void th3_hflux_rst(DOUBLE xjm[NDIM_THERM3][NDIM_THERM3],
                   DOUBLE *hflux,
                   INT igp,
                   DOUBLE **hfluxrst);
void th3_hflux_modang(DOUBLE *hflux,
                      INT igp,
                      DOUBLE **hflux123);
void th3_hflux_extrpol(ELEMENT *ele,
                       TH3_DATA *data,
                       INT ngauss,
                       DOUBLE **hfluxgp,
                       DOUBLE *rst,
                       DOUBLE *hfluxnd);

/*----------------------------------------------------------------------*/
/* file th3_inp.c */
void th3_inp(ELEMENT *ele);

/*----------------------------------------------------------------------*/
/* file th3_intg.c */
void th3_intg_eleinp(ELEMENT *actele,
                     INT *ierr);
void th3_intg_init(TH3_DATA *data);

/*----------------------------------------------------------------------*/
/* file th3_metr.c */
void th3_metr_jaco(ELEMENT *ele,
                   INT      enod,
                   DOUBLE   deriv[MAXNOD_THERM3][NDIM_THERM3],
                   INT      flag,
                   DOUBLE   xjm[NDIM_THERM3][NDIM_THERM3],
                   DOUBLE  *det,
                   DOUBLE   xji[NDIM_THERM3][NDIM_THERM3]);
void th3_metr_surf(ELEMENT *ele, 
                   INT      nelenod, 
                   DOUBLE   deriv[MAXNOD_THERM3][NDIM_THERM3], 
                   DOUBLE   sidredm[DIMSID_THERM3][NDIM_THERM3],
                   DOUBLE  *metr);
void th3_metr_line(ELEMENT *ele, 
                   INT      nelenod, 
                   DOUBLE   deriv[MAXNOD_THERM3][NDIM_THERM3], 
                   DOUBLE   linredv[NDIM_THERM3],
                   DOUBLE  *metr);

/*----------------------------------------------------------------------*/
/* file th2_lin.c */
void th3_lin_tang(CONTAINER* cont,
                  ELEMENT* ele,
                  TH3_DATA* data,
                  MATERIAL* mat,
                  ARRAY* estif_global,
                  ARRAY* emass_global,
                  DOUBLE* force);
void th3_lin_bcb(INT neledof,
                 DOUBLE bop[NDIM_THERM3][NUMDOF_THERM3*MAXNOD_THERM3],
                 DOUBLE cmat[NUMHFLX_THERM3][NUMTMGR_THERM3],
                 DOUBLE fac,
                 DOUBLE** tmat);
void th3_lin_fint(INT neledof,
                  DOUBLE bop[NDIM_THERM3][NUMDOF_THERM3*MAXNOD_THERM3],
                  DOUBLE hflux[NUMHFLX_THERM3],
                  DOUBLE fac,
                  DOUBLE* intfor);

/*----------------------------------------------------------------------*/
/* file th3_load.c */
void th3_load_heat(ELEMENT* ele,
                   TH3_DATA* data,
                   INT imyrank,
                   DOUBLE* loadvec);
void th3_load_vol(ELEMENT *ele,
                  INT nelenod,
                  DOUBLE shape[MAXNOD_THERM3],
                  DOUBLE fac,
                  DOUBLE eload[NUMDOF_THERM3][MAXNOD_THERM3]);
void th3_load_surf(ELEMENT *ele,
                   INT nelenod,
                   GSURF *gsurf,
                   DOUBLE shape[MAXNOD_THERM3],
                   DOUBLE fac,
                   DOUBLE eload[NUMDOF_THERM3][MAXNOD_THERM3]);
void th3_load_line(ELEMENT *ele,
                   INT nelenod,
                   GLINE *gline,
                   DOUBLE shape[MAXNOD_THERM3],
                   DOUBLE fac,
                   DOUBLE eload[NUMDOF_THERM3][MAXNOD_THERM3]);



/*----------------------------------------------------------------------*/
/* file th3_main.c */
void therm3(PARTITION *actpart,
            INTRA *actintra,
            ELEMENT *ele,
            ARRAY *estif_global,
            ARRAY *emass_global,
            ARRAY *intforce_global,
            CALC_ACTION *action,
            CONTAINER *container);   /* contains variables defined 
                                      * in container.h */

/*----------------------------------------------------------------------*/
/* file th3_mat.c */
void th3_mat_sel(CONTAINER *cont,
                 ELEMENT *ele,
                 MATERIAL *mat,
                 DOUBLE bop[NDIM_THERM3][NUMDOF_THERM3*MAXNOD_THERM3],
                 INT ip,
                 DOUBLE heatflux[NUMHFLX_THERM3],
                 DOUBLE cmat[NUMHFLX_THERM3][NUMTMGR_THERM3]);

/*----------------------------------------------------------------------*/
/* file th3_matlin.c */
void th3_matlin_iso(CONTAINER *cont,
                    DOUBLE con,
                    ELEMENT *ele,
                    DOUBLE bop[NDIM_THERM3][NUMDOF_THERM3*MAXNOD_THERM3],
                    DOUBLE heatflux[NUMHFLX_THERM3],
                    DOUBLE cmat[NUMHFLX_THERM3][NUMTMGR_THERM3]);
void th3_matlin_gen(DOUBLE **con,
                    ELEMENT *ele,
                    DOUBLE bop[NDIM_THERM3][NUMDOF_THERM3*MAXNOD_THERM3],
                    DOUBLE heatflux[NUMHFLX_THERM3],
                    DOUBLE cmat[NUMHFLX_THERM3][NUMTMGR_THERM3]);

/*----------------------------------------------------------------------*/
/* file th3_out.c */
void th3_out_hflux(ELEMENT *actele,
                   FILE *out);

/*----------------------------------------------------------------------*/
/* file th3_shape.c */
void th3_shape_deriv(DIS_TYP     typ,
                     DOUBLE      r,
                     DOUBLE      s,
                     DOUBLE      t,
                     INT         option,
                     DOUBLE      shape[MAXNOD_THERM3],
                     DOUBLE      deriv[MAXNOD_THERM3][NDIM_THERM3]);

/*----------------------------------------------------------------------*/
/* file th2_temper.c */
void th3_temper_init();
void th3_temper_final();
void th3_temper_cal(const CONTAINER *container,
                    const ELEMENT *ele,
                    const DOUBLE r,
                    const DOUBLE s,
                    const DOUBLE t,
                    DOUBLE *tem);

/*----------------------------------------------------------------------*/
#endif /*end of #ifdef D_THERM3 */
/*! @} (documentation module close) */
