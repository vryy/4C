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

/*======================================================================*/
/* defines constants of THERM3 */
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
#define GTMAXP_THERM3    (5)    /* max. number of Gauss points */
#endif

#ifndef GSINTC_THERM3
#define GSINTC_THERM3    (5)    /* triangle domain Gauss integration cases */
#endif

#ifndef GSMAXP_THERM3
#define GSMAXP_THERM3    (6)    /* max. number of Gauss points */
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
  DOUBLE gtdcr[GTINTC_THERM3][GTMAXP_THERM3];  /* coordinates in r */
  DOUBLE gtdcs[GTINTC_THERM3][GTMAXP_THERM3];   /* coordinates in s */
  DOUBLE gtdct[GTINTC_THERM3][GTMAXP_THERM3];  /* coordinates in t */
  DOUBLE gtdw[GTINTC_THERM3][GTMAXP_THERM3];  /* weights */
  /* tetrahedron sides */
  DOUBLE gtscr[GSINTC_THERM3][GSMAXP_THERM3];  /* coordinates in r */
  DOUBLE gtscs[GSINTC_THERM3][GSMAXP_THERM3];  /* coordinates in s */
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
  /* sides */
  DOUBLE ancsidh[MAXSID_THERM3][NDIM_THERM3];  /* anchors hex */
  INT dirsidh[MAXSID_THERM3][NDIM_THERM3];  /* base vector directions */
  /* edges */
  DOUBLE ancedgh[MAXEDG_THERM3][NDIM_THERM3];  /* anchors hex */
  INT diredgh[MAXEDG_THERM3][NDIM_THERM3];  /* base vector direction */
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
   *       nGP[0] : in r-direction : 1,2,3,4,5,6 : read-in
   *       nGP[1] : in s-direction : 1,2,3,4,5,6 : read-in/set
   *       nGP[2] : in t-direction : 1,2,3,4,5,6 : read-in/set
   * tetrahedra:
   *       nGP[0] : total number of GPs in domain: 1,4,5 : read-in
   *       nGP[1] : total number of GPs on sides: 1,3,4,6 : read-in/set
   *       nGP[2] : total number of GPs on edges: 1,2,3,4,5,6 : read-in/set */
  INT gpnum[NDIM_THERM3];
  INT gpintc[NDIM_THERM3];

  ARRAY4D hflux_gp;
  ARRAY4D hflux_nd;

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
             DOUBLE   **deriv,
             DOUBLE   **xji,
             DOUBLE   **bop);

/*---------------------------------------------------------------------*/
/* file th3_cfg.c */
void th3_cfg_chkdef();
void th3_cfg_init(TH3_DATA *data);
void th3_cfg_noderst(ELEMENT *ele,
		     TH3_DATA *data,
		     INT inode,
		     DOUBLE *rst);

/*----------------------------------------------------------------------*/
/* file th2_hflux.c */
/* void th2_hflux_init(PARTITION *actpart); */
/* void th2_hflux_final(PARTITION *actpart); */
/* void th2_hflux_cal(ELEMENT *ele, */
/*                    THERM2_DATA *data, */
/*                    MATERIAL *mat, */
/*                    INT kstep); */
/* void th2_hflux_steep(DOUBLE *hflux, */
/*                      DOUBLE *hfluxmod, */
/*                      DOUBLE *hfluxang, */
/*                      DOUBLE *dum); */
/* void th2_hflux_extrpol(ELEMENT *ele, */
/*                        THERM2_DATA *data, */
/*                        INT ngauss, */
/*                        DOUBLE *hfluxgp, */
/*                        DOUBLE *rs, */
/*                        DOUBLE *hfluxnd); */

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
void th3_jaco(ELEMENT *ele,
              INT      enod,
              DOUBLE **deriv,
              INT      flag,
              DOUBLE **xjm,
              DOUBLE  *det,
              DOUBLE **xji);

/*----------------------------------------------------------------------*/
/* file th2_lin.c */
void th3_lin_tang(ELEMENT *ele,
                  TH3_DATA *data,
                  MATERIAL *mat,
                  ARRAY *estif_global,
                  ARRAY *emass_global,
                  DOUBLE *force);
void th3_lin_bcb(INT       neledof,
                 DOUBLE  **bop,
                 DOUBLE  **cmat,
                 DOUBLE    fac,
                 DOUBLE  **tmat);
void th3_lin_fint(INT      neledof,
                  DOUBLE **bop,
                  DOUBLE  *hflux,
                  DOUBLE   fac,
                  DOUBLE  *intfor);

/*----------------------------------------------------------------------*/
/* file th2_load.c */
/* void th2_load_init(); */
/* void th2_load_final(); */
/* void th2_load_heat(ELEMENT *ele, */
/*                    THERM2_DATA *data, */
/*                    DOUBLE *loadvec, */
/*                    INT imyrank); */
/* void th2_load_heatsurf(ELEMENT *ele, */
/*                        DOUBLE **eload, */
/*                        DOUBLE *funct, */
/*                        DOUBLE fac, */
/*                        INT iel); */

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
void th3_mat_sel(ELEMENT *ele,
                 MATERIAL *mat,
                 DOUBLE **bop,
                 INT ip,
                 DOUBLE *heatflux,
                 DOUBLE **cmat);

/*----------------------------------------------------------------------*/
/* file th3_matlin.c */
void th2_matlin_iso(DOUBLE con,
                    ELEMENT *ele,
                    DOUBLE **bop,
                    DOUBLE *heatflux,
                    DOUBLE **cmat);
void th3_matlin_gen(DOUBLE **con,
                    ELEMENT *ele,
                    DOUBLE **bop,
                    DOUBLE *heatflux,
                    DOUBLE **cmat);

/*----------------------------------------------------------------------*/
/* file th3_shape.c */
void th3_shape_deriv(DIS_TYP     typ,
                     DOUBLE      r,
                     DOUBLE      s,
                     DOUBLE      t,
                     INT         option,
                     DOUBLE     *shape,
                     DOUBLE    **deriv);

/*----------------------------------------------------------------------*/
/* file th2_temper.c */
/* void th2_temper_init(); */
/* void th2_temper_final(); */
/* void th2_temper_cal(ELEMENT *ele, */
/*                      DOUBLE r, */
/*                      DOUBLE s, */
/*                      DOUBLE *tem); */

/*----------------------------------------------------------------------*/
#endif /*end of #ifdef D_THERM3 */
/*! @} (documentation module close) */
