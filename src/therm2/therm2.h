/*======================================================================*/
/*!
\file
\brief headerfile for planar thermal element (therm2), containing 
structures and prototypes

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/
/*======================================================================*/

#ifdef D_THERM2

/*!
\addtogroup THERM2
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"


/*======================================================================*/
/* global declarations, variables etc */

/*----------------------------------------------------------------------*/
/*!
\brief Gauss points and weights

The constants are defined in headers/define_sizes.h

\author bborn
\date 03/06
*/
typedef struct _THERM2_DATA
{
  /* quadrilateral domain and edges --> line [-1,+1] */
  DOUBLE gqlc[GLMAXP_THERM2][GLINTC_THERM2];  /* coordinates on line [-1,+1] */
  DOUBLE gqlw[GLMAXP_THERM2][GLINTC_THERM2];  /* weights on line [-1,+1] */
  /* triangle domain */
  DOUBLE gtdcr[GTMAXP_THERM2][GTINTC_THERM2];  /* coordinates in r */
  DOUBLE gtdcs[GTMAXP_THERM2][GTINTC_THERM2];  /* coordinates in s */
  DOUBLE gtdw[GTMAXP_THERM2][GTINTC_THERM2];  /* weights */
  /* triangle edges --> line [0,+1] */
  DOUBLE gtlc[GLMAXP_THERM2][GLINTC_THERM2];  /* coordinates on line [0,+1] */
  DOUBLE gtlw[GLMAXP_THERM2][GLINTC_THERM2];  /* weights on line [0,+1] */
} THERM2_DATA;


/*----------------------------------------------------------------------*/
/*!
\brief Type of plane states, ie planar heat flux or planar temperature
gradient

IMPORTANT: The magnitudes of types is globally defined in 
define_sizes.h  :  #define NUMPLST_THERM2   (2)
If types are added here, then this number MUST be updated accordingly!

\author bborn
\date 03/06
*/
typedef enum _TH2_PLANESTATES
{
                       th2_plane_tempgrad,
                       th2_plane_heatflux
} TH2_PLANESTATES;


/*----------------------------------------------------------------------*/
/*!
\brief Type of kinematics

\author bborn
\date 03/06
*/
typedef enum _TH2_KINEMATIK_TYPE
{
  th2_geo_lin,  /* geometrically linear */
  th2_total_lagr,  /* geometrically non-linear => total Lagrangian */
  th2_updated_lagr  /* geometrically non-linear => updated Lagrangian */
} TH2_KINEMATIK_TYPE;


/*----------------------------------------------------------------------*/
/*!
\brief Type of heat flux output

\author bborn
\date 03/06
*/
typedef enum _TH2_HEATFLUXOUT_TYPE
{
  th2_hflux_none,
  th2_hflux_gpxy,  /* globally xy-oriented at Gauss points */
  th2_hflux_gp12,  /* modulus and its direction */
  th2_hflux_gprs,  /* locally rs-oriented at Gauss points */
  th2_hflux_gpxy12,  /* all of the above */
  th2_hflux_ndxy,
  th2_hflux_nd12,
  th2_hflux_ndrs,
  th2_hflux_ndxy12
} TH2_HEATFLUXOUT_TYPE;

/*----------------------------------------------------------------------*/
/*!
\brief Definition of THERM2 type holding THERM2 element properties

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
typedef struct _THERM2
{
  enum _TH2_PLANESTATES planestate;  /* type of plane state */
  enum _TH2_KINEMATIK_TYPE kintype;  /* type of Kinematik */
  enum _TH2_HEATFLUXOUT_TYPE hfluxtype;  /* output type of heat flux */

  INT nGP[2];  /* number of Gauss points as obtained during read-in
                * quads:
                *       nGP[0] : in r-direction 
                *       nGP[1] : in s-direction
                * tris:
                *       nGP[0] : total number of GPs: 1,3,4,6,7,9,12,13
                *       nGP[1] : set 0 or 1         :   x   x x
                */
  INT gptotal;  /* total number of GPs */
  INT gpintc;  /* GP integration case; important to triangles */
  INT gpned;  /* amount of GPs on line/edge, important for triangles */
    
  DOUBLE thick;  /* thickness */

  struct _ARRAY4D hflux_gp;
  struct _ARRAY4D hflux_nd;

#ifdef D_TSI
  TSI_COUPTYP tsi_couptyp;
  ELEMENT *struct_ele;
#endif
} THERM2;



/*======================================================================*/
/* Declarations of functions in therm2 directory
 * Order: Firstly, alphabetically list file names, secondly, list
 *        functions in file according to appearance */


/*----------------------------------------------------------------------*/
/* file th2_bop.c */
void th2_bop(DOUBLE    **bop,
             DOUBLE    **deriv,
             DOUBLE    **xjm,
             DOUBLE      det,
             INT         iel);

/*---------------------------------------------------------------------*/
/* file th2_cfg.c */
void th2_cfg_init();
void th2_cfg_iedg(INT *iegnod, 
                  ELEMENT *ele, 
                  INT line);
void th2_cfg_noders(ELEMENT *ele,
                    INT inode,
                    DOUBLE *rs);

/*----------------------------------------------------------------------*/
/* file th2_hflux.c */
void th2_hflux_init(PARTITION *actpart);
void th2_hflux_final(PARTITION *actpart);
void th2_hflux_cal(ELEMENT *ele,
                   THERM2_DATA *data,
                   MATERIAL *mat,
                   INT kstep);
void th2_hflux_steep(DOUBLE *hflux,
                     DOUBLE *hfluxmod,
                     DOUBLE *hfluxang,
                     DOUBLE *dum);
void th2_hflux_extrpol(ELEMENT *ele,
                       THERM2_DATA *data,
                       INT ngauss,
                       DOUBLE *hfluxgp,
                       DOUBLE *rs,
                       DOUBLE *hfluxnd);

/*----------------------------------------------------------------------*/
/* file th2_inp.c */
void th2_inp(ELEMENT *ele);


/*----------------------------------------------------------------------*/
/* file th2_intg.c */
void th2_intg_eleinp(ELEMENT *actele, 
                     INT *ierr);
void th2_intg_init(THERM2_DATA *data);

/*----------------------------------------------------------------------*/
/* file th2_jaco.c */
void th2_jaco(DOUBLE **deriv,
              DOUBLE **xjm, 
              DOUBLE *det,
              ELEMENT *ele,
              INT iel);

/*----------------------------------------------------------------------*/
/* file th2_lin.c */
void th2_lin_init();
void th2_lin_final();
void th2_lin_stiff(ELEMENT *ele,
                   THERM2_DATA *data,
                   MATERIAL *mat,
                   ARRAY *estif_global,
                   ARRAY *emass_global,
                   DOUBLE *force);
void th2_lin_bcb(DOUBLE  **s,
                 DOUBLE  **bs,
                 DOUBLE  **d,
                 DOUBLE    fac,
                 INT       neledof,
                 INT       ntmgr);
void th2_lin_fint(DOUBLE  *hflux,
                  DOUBLE   fac,
                  DOUBLE **bop,
                  INT      neledof,
                  DOUBLE  *fie);

/*----------------------------------------------------------------------*/
/* file th2_load.c */
void th2_load_init();
void th2_load_final();
void th2_load_heat(ELEMENT *ele,
                   THERM2_DATA *data,
                   DOUBLE *loadvec,
                   INT imyrank);
void th2_load_heatsurf(ELEMENT *ele,
                       DOUBLE **eload,
                       DOUBLE *funct,
                       DOUBLE fac,
                       INT iel);

/*----------------------------------------------------------------------*/
/* file th2_main.c */
void therm2(PARTITION *actpart,
	    INTRA *actintra,
	    ELEMENT *ele,
	    ARRAY *estif_global,
	    ARRAY *emass_global,
	    ARRAY *intforce_global,
	    CALC_ACTION *action,
	    CONTAINER *container);   /* contains variables defined 
				      * in container.h */

/*----------------------------------------------------------------------*/
/* file th2_mat.c */
void th2_mat_sel(ELEMENT   *ele,
                 MATERIAL  *mat,
                 DOUBLE **bop,
                 INT ip,
                 DOUBLE *heatflux,
                 DOUBLE **cmat);

/*----------------------------------------------------------------------*/
/* file th2_matlin.c */
void th2_matlin_iso(DOUBLE con,
                    ELEMENT *ele,
                    DOUBLE **bop,
                    DOUBLE *heatflux,
                    DOUBLE **cmat);
void th2_matlin_gen(DOUBLE **con,
                    ELEMENT *ele,
                    DOUBLE **bop,
                    DOUBLE* heatflux,
                    DOUBLE **cmat);

/*----------------------------------------------------------------------*/
/* file th2_shape.c */
void th2_shape_deriv(DOUBLE     *shape,
                     DOUBLE    **deriv,
                     DOUBLE      r,
                     DOUBLE      s,
                     DIS_TYP     typ,
                     INT         option);

/*----------------------------------------------------------------------*/
/* file th2_temper.c */
void th2_temper_init();
void th2_temper_final();
void th2_temper_cal(ELEMENT *ele,
                     DOUBLE r,
                     DOUBLE s,
                     DOUBLE *tem);

/*----------------------------------------------------------------------*/
#endif /*end of #ifdef D_THERM2 */
/*! @} (documentation module close) */
