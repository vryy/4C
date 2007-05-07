/*======================================================================*/
/*!
\file
\brief headerfile for 3dim thermal element (SOLID3), containing 
structures and prototypes

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089-289-15240
</pre>
*/
/*======================================================================*/

#ifdef D_SOLID3

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "../output/gid.h"

/*======================================================================*/
/* defines constants of SOLID3 */
#define NDIM_SOLID3      (3)    /* 3dim problem */

#ifndef MAXNOD_SOLID3
#ifdef MAXNOD
#define MAXNOD_SOLID3    (MAXNOD)
#else
#define MAXNOD_SOLID3    (27)   /* maximum of element nodes : max hex27 */
#endif
#endif

#ifndef MAXEDG_SOLID3
#define MAXEDG_SOLID3    (12)   /* maximal number element edges */
#endif

#ifndef MAXSID_SOLID3
#define MAXSID_SOLID3    (6)    /* maximal number element sides */
#endif

#ifndef MAXNE_SOLID3
#define MAXNE_SOLID3     (3)    /* maximal number of nodes on an edge */
#endif

#ifndef MAXNS_SOLID3
#define MAXNS_SOLID3     (9)    /* maximal number of nodes on a side */
#endif

#define DIMSID_SOLID3    (2)    /* dimension of element faces */

#define NUMDOF_SOLID3    (3)    /* number of structural DOFs at each node :
                                 * displacement ux uy uz */
#define NUMDFGR_SOLID3   (9)    /* number of deformation gradient components
                                 * (deformation gradient is non-symmetric */
#define NUMSTR_SOLID3    (6)    /* number of strains and stresses :
				 * GL strain: E_XX,E_YY,E_ZZ,2E_XY,2E_YZ,2E_ZX
                                 * 2PK stress: S_XX,S_YY,S_ZZ,S_XY,S_YZ,S_ZX */

#define MAXDOF_SOLID3    (MAXNOD_SOLID3*NUMDOF_SOLID3)  /* maximal element
                                                         * DOFs */

#ifndef MAXGAUSS_SOLID3
#ifdef MAXGAUSS
#define MAXGAUSS_SOLID3  (MAXGAUSS)
#else
#define MAXGAUSS_SOLID3  (27)   /* maximum of total Gauss points in domain */
#endif
#endif

#ifndef GLINTC_SOLID3
#define GLINTC_SOLID3    (6)    /* line domain Gauss integration cases */
#endif

#ifndef GLMAXP_SOLID3
#define GLMAXP_SOLID3    (6)    /* line domain max. number of Gauss points */
#endif
#ifndef GTINTC_SOLID3
#define GTINTC_SOLID3    (3)    /* tetrahedron domain Gauss integration cases*/
#endif

#ifndef GTMAXP_SOLID3
#define GTMAXP_SOLID3    (5)    /* tet domain max. number of Gauss points */
#endif

#ifndef GSINTC_SOLID3
#define GSINTC_SOLID3    (5)    /* triangle domain Gauss integration cases */
#endif

#ifndef GSMAXP_SOLID3
#define GSMAXP_SOLID3    (6)    /* triangle max. number of Gauss points */
#endif

#ifdef DEBUG
/*#define TEST_SOLID3*/             /* compile test routines if ... */
#endif

/* debug: */
/*#define TESTROBIN_SOLID3*/    /* Special debug areas are activated during
                                 * development of the Robinson material.
                                 * This CPP definition should be turned
                                 * off as soon as possible */

/* quick hack */
/*#define TANGLINELOAD_SOLID3*/ /* quick hack to apply tangential 
                                 * 'line load' on a cylindrical 
                                 * cantilever beam (length 120) subjected 
                                 * to a torsional torque at its tip. 
                                 * (bborn/mgit 04/07) */

/*======================================================================*/
/* global declarations, variables etc */

/*----------------------------------------------------------------------*/
/*!
\brief Gauss points and weights. This is static data.

Most of the preprocessorily constants are defined on the top of the file.
A few constants are defined in headers/define_sizes.h.

\author mf
\date 10/06
*/
typedef struct _SO3_DATA
{
  /*--------------------------------------------------------------------*/
  /* Gauss coordinates and weights */
  /* hexahedron domain, sides and edges --> line [-1,+1] */
  DOUBLE ghlc[GLINTC_SOLID3][GLMAXP_SOLID3];  /* coordinates */
  DOUBLE ghlw[GLINTC_SOLID3][GLMAXP_SOLID3];  /* weights */
  /* tetrahedron domain [T.J.R. Hughes, "The FEM", Dover 2000] */
  DOUBLE gtdc[GTINTC_SOLID3][GTMAXP_SOLID3][NDIM_SOLID3];  /* coordinates 
                                                            * in r,s,t */
  DOUBLE gtdw[GTINTC_SOLID3][GTMAXP_SOLID3];  /* weights */
  /* tetrahedron sides */
  DOUBLE gtsc[GSINTC_SOLID3][GSMAXP_SOLID3][DIMSID_SOLID3];  /* coordinates 
                                                              * in side */
  DOUBLE gtsw[GSINTC_SOLID3][GSMAXP_SOLID3];  /* weights */
  /* triangle edges --> line [0,+1] */
  DOUBLE gtlc[GLINTC_SOLID3][GLMAXP_SOLID3];  /* coordinates */
  DOUBLE gtlw[GLINTC_SOLID3][GLMAXP_SOLID3];  /* weights */
  /*--------------------------------------------------------------------*/
  /* numbering of element nodes, edges, sides in parameter space */
  /* parameter coordinates of nodes */
  DOUBLE nodhrst[MAXNOD_SOLID3][NDIM_SOLID3];  /* hexahedron */
  DOUBLE nodtrst[MAXNOD_SOLID3][NDIM_SOLID3];  /* tetrahedron */
  /* nodes on sides (surfaces) */
  INT nodsidh[MAXSID_SOLID3][MAXNS_SOLID3];  /* hexahedra */
  INT nodsidt[MAXSID_SOLID3][MAXNS_SOLID3];  /* tetrahedra */
  /* nodes on edges */
  INT nodedghl[MAXEDG_SOLID3][MAXNE_SOLID3];  /* linear hex8 */       
  INT nodedghq[MAXEDG_SOLID3][MAXNE_SOLID3];  /* quadratic hex20,27 */
  INT nodedgtl[MAXEDG_SOLID3][MAXNE_SOLID3];  /* linear tet4 */       
  INT nodedgtq[MAXEDG_SOLID3][MAXNE_SOLID3];  /* quadratic tet10 */
  /*--------------------------------------------------------------------*/
  /* anchor and span vectors for sides and edges in param. space */
  /* sides hex */
  DOUBLE ancsidh[MAXSID_SOLID3][NDIM_SOLID3];  /* anchors hex */
  DOUBLE redsidh[MAXSID_SOLID3][DIMSID_SOLID3][NDIM_SOLID3];  /* dim red 
                                                               * matrix */
  /* sides tet */
  DOUBLE ancsidt[MAXSID_SOLID3][NDIM_SOLID3];  /* anchors tet */
  DOUBLE redsidt[MAXSID_SOLID3][DIMSID_SOLID3][NDIM_SOLID3];  /* dim red 
                                                               * matrix */
  /* edges hex */
  DOUBLE ancedgh[MAXEDG_SOLID3][NDIM_SOLID3];  /* anchors hex */
  DOUBLE rededgh[MAXEDG_SOLID3][NDIM_SOLID3];  /* dimension reduction
                                                * matrix multiplied
                                                * on Jacobi matrix */
  /* edges tet */
  DOUBLE ancedgt[MAXEDG_SOLID3][NDIM_SOLID3];  /* anchors tet */
  DOUBLE rededgt[MAXEDG_SOLID3][NDIM_SOLID3];  /* dimension reduction
                                                * matrix multiplied
                                                * on Jacobi matrix */
  /*--------------------------------------------------------------------*/
  /* base vectors of surface plus outward normal are common 
   * coordinate system (Rechtssystem) if ==0 else ==1 not */
  INT nrmsidh[MAXEDG_SOLID3];  /* hex */
  INT nrmsidt[MAXEDG_SOLID3];  /* tet */
} SO3_DATA;



/*----------------------------------------------------------------------*/
/*!
\brief All Gauss point coordinates, shape functions and their 
       parametric derivates evaluated
\author bborn
\date 12/06
*/
typedef struct _SO3_GPSHAPEDERIV
{
  DIS_TYP distyp;  /* discretisation */
  INT gpintc[NDIM_SOLID3];  /* Gauss integration case */
  INT gptot;  /* total number of Gauss points in domain */
  DOUBLE gpco[MAXGAUSS_SOLID3][NDIM_SOLID3];
  DOUBLE gpwg[MAXGAUSS_SOLID3];
  DOUBLE gpshape[MAXGAUSS_SOLID3][MAXNOD_SOLID3];
  DOUBLE gpderiv[MAXGAUSS_SOLID3][MAXNOD_SOLID3][NDIM_SOLID3];
} SO3_GPSHAPEDERIV;


/*----------------------------------------------------------------------*/
/*!
\brief The (mostly tensorial) variables of this datum type describe
       the geometric transformations between the 3 frames:
       parameter (r,s,t), material (X,Y,Z) and spatial (x,y,z).
       Moreover a few (mostly tensorial) variables defined in 
       these frames are contained.
       These data are handed down to the material routines.
\author bborn
\date 12/06
*/
typedef struct _SO3_GEODEFSTR
{
  /*--------------------------------------------------------------------*/
  /* geometry:
   *     These variables refer to the transformation between 
   *     parameter frame (r,s,t coordinates) and material frame
   *     (X,Y,Z coordinates)
   */
  /*--------------------------------------------------------------------*/
  /* (total) Gauss point index */
  INT igp;
  /* Gauss point co-ordinate */
  DOUBLE gpc[NDIM_SOLID3];
  /* (FE-) Jacobi matrix J (isoparametric)
   *         [ J_11  J_12  J_13 ]   [ X_{,r}  Y_{,r}  Z_{,r} ]
   *     J = [ J_21  J_22  J_23 ] = [ X_{,s}  Y_{,s}  Z_{,s} ]
   *         [ J_31  J_32  J_33 ]   [ X_{,t}  Y_{,t}  Z_{,t} ]
   * HINT: This is the classical ordering found in FE text books.
   *       This ordering transposes the ordering fashion
   *       of J viewed as the `parametric deformation gradient'
   *       (being responsible for the `deformation' of paramter
   *       space to material space) */
  DOUBLE xjm[NDIM_SOLID3][NDIM_SOLID3];
  /* Jacobi determinant det(J) */
  DOUBLE xjdet;
  /* inverted Jacobi matrix J^{-1}*/
  DOUBLE xji[NDIM_SOLID3][NDIM_SOLID3];
  /* rotational component of  J = R . U   ie  R = J . U^{-1}  */
  /* HINT: This matrix adheres to the ordering of xjm */
  DOUBLE xrm[NDIM_SOLID3][NDIM_SOLID3];
  /* rotational component of J prepared such that a symmetric 2-tensor
   * A living in (rst) denoted vectorially Av can be computed 
   * by matrix-vector product
   *     (A)_{XYZ} = (xrm) . (A)_{rst} . (xrm^T)
   * becomes
   *     (Av)_{XYZ} = xrvm . (Av)_{rst} */
  DOUBLE xrvm[NUMSTR_SOLID3][NUMSTR_SOLID3];
  /* rotational component of J^{-1} prepared such that a symmetric 2-tensor
   * A living in (XYZ) denoted vectorially Av can be computed 
   * by matrix-vector product 
   *     (A)_{rst} = (xrm^T) . (A)_{XYZ} . (xrm)
   * becomes
   *     (Av)_{rst} = xrvi . (Av)_{XYZ} */
  DOUBLE xrvi[NUMSTR_SOLID3][NUMSTR_SOLID3];
  /*--------------------------------------------------------------------*/
  /* deformation: 
   *     These variables refer to the transformation between
   *     material frame (X,Y,Z coordinates) and spatial frame
   *     (x,y,z coordinates)
   */
  /* material deformation gradient (commonly denoted F)
   *         [ F_11  F_12  F_13 ]
   *     F = [ F_21  F_22  F_23 ]
   *         [ F_31  F_32  F_33 ] */
  DOUBLE defgrd[NDIM_SOLID3][NDIM_SOLID3];
  /* material deformation gradient in vectorial notion Fv
   *     Fv^T = [ F_11  F_22  F_33  F_12  F_21  F_23  F_32  F_31  F_13 ] */
  /*  DOUBLE defgrdv[NUMDFGR_SOLID3]; */
  /* material displacement gradient in vectorial notion Hv
   *     Hv^T = [ u_,X  v_,Y  w_,Z  u_,Y  v_,X  v_,Z  w_,Y  w_,X  u_,Z ] */
  DOUBLE disgrdv[NUMDFGR_SOLID3];
  /*--------------------------------------------------------------------*/
  /* strain:
   *     Strains referred to different frames */
  /* Green-Lagrange strain tensor in vector notion Ev
   *     Ev^T = [ E_11  E_22  E_33  2*E_12  2*E_23  2*E_31 ] */
  DOUBLE stnglv[NUMSTR_SOLID3];
  /* Linear (engineering strain) in vector notion epsv
   *     epsv^T = [ eps_11  eps_22  eps_33  2*eps_12  2*eps_23  2*eps_31 ] */
  DOUBLE stnengv[NUMSTR_SOLID3];
  /*--------------------------------------------------------------------*/
  /* stress:
   *     Stresses referred to different frames */
  /* 2nd Piola-Kirchhoff stress vector Sv
   *     Sv^T = [ S_11  S_22  S_33  S_12  S_23  S_31 ] */
  /* DOUBLE sts2pkv[NUMSTR_SOLID3]; */
  /* 1st Piola-Kirchhoff stress vector Pv
   *     Pv^T = [ P_11  P_22  P_33  P_12  P_21  P_23  P_32  P_31  P_13 ] */
  /* DOUBLE sts1pkv[NUMDFGR_SOLID3]; */
  /*--------------------------------------------------------------------*/
  /* B-operators */
  /* nodal B-operator */
  DOUBLE bopn[MAXNOD_SOLID3][NDIM_SOLID3];
  /* B-operator (geometrically linear/non-linear) */
  DOUBLE bop[NUMSTR_SOLID3][MAXDOF_SOLID3];
} SO3_GEODEFSTR;


/*----------------------------------------------------------------------*/
/*!
\brief Type of kinematics

\author mf
\date 10/06
*/
typedef enum _SO3_KINEMATICS
{
  so3_geo_lin,      /* geometrically linear */
  so3_total_lagr,   /* geometrically non-linear => total Lagrangian */
  so3_updated_lagr  /* geometrically non-linear => updated Lagrangian */
} SO3_KINEMATICS;


/*----------------------------------------------------------------------*/
/*!
\brief Type of stress output

\author mf
\date 10/06
*/
#define SOLID3_STRESSTYPE { "so3_stress_none",   \
      "so3_stress_gpxyz",                        \
      "so3_stress_gp123",                        \
      "so3_stress_gprst",                        \
      "so3_stress_gpxyz123",                     \
      "so3_stress_ndxyz",                        \
      "so3_stress_nd123",                        \
      "so3_stress_ndrst",                        \
      "so3_stress_ndxyz123",                     \
      "so3_stress_ndeqv",                        \
      NULL }
typedef enum _SO3_STRESSOUT
{
  so3_stress_none,
  so3_stress_gpxyz,  /* globally xyz-oriented at Gauss points */
  so3_stress_gp123,  /* in principal directions */
  so3_stress_gprst,  /* locally rs-oriented at Gauss points */
  so3_stress_gpxyz123,  /* all of the above */
  so3_stress_ndxyz,     /* the same for nodes*/
  so3_stress_nd123,
  so3_stress_ndrst,
  so3_stress_ndxyz123,
  so3_stress_ndeqv
} SO3_STRESSOUT;

/*----------------------------------------------------------------------*/
/*!
\brief Mode/Status in Robinson's visco-plastic material
\author bborn
\date 05/07
*/
typedef enum _SO3_MAT_ROBINSON_STATE
{
  so3_mat_robinson_state_vague = 0,
  so3_mat_robinson_state_elastic,
  so3_mat_robinson_state_inelastic,
} SO3_MAT_ROBINSON_STATE;


/*----------------------------------------------------------------------*/
/*!
\brief Material internal variables (MIV) used for the
       visco-plastic Robinson material
       These fields are allocated at every Gauss point.
       Not every field is allocated, it depends on the used
       time integration of the MIV flow rules.
       ==> so3_mat_robinson.c
           ==> so3_mat_robinson_fb4.c  (Fehlberg4)
           ==> so3_mat_robinson_be.c  (Backward Euler)

\note Mief, Mief, Miiiieeeef!

\author bborn
\date 03/07
*/
typedef struct _SO3_MIV_ROBINSON
{
  /* visco-plsatic strain state (mode) at each <g>
   *    ==0: undetermined;  ==1: 'elastic';  ==2: 'viscous' */
  ARRAY vscstns;
  /* visco-plastic strain vector Ev^<g> at t_{n} for every Gauss point g
   *    Ev^<g>T = [ E_11  E_22  E_33  2*E_12  2*E_23  2*E_31 ]^<g> */
  ARRAY vicstn;
  /* visco-plastic strain vector Ev^<g> at t_{n+1} for every Gauss point g
   *    Ev^<g>T = [ E_11  E_22  E_33  2*E_12  2*E_23  2*E_31 ]^<g> */
  ARRAY vicstnn;
  /* back stress state (mode) at each <g>
   *    ==0: undetermined;  ==1: 'elastic';  ==2: 'viscous' */
  ARRAY bckstss;
  /* back stress vector Av^<g> at t_n for every Gauss point g
   *    Av^<g>T = [ A_11  A_22  A_33  A_12  A_23  A_31 ]^<g> */
  ARRAY bacsts;
  /* back stress vector Av^<g> at t_{n+1} for every Gauss point g
   *    Av^<g>T = [ A_11  A_22  A_33  A_12  A_23  A_31 ]^<g> */
  ARRAY bacstsn;
  /* visco-plastic strain rate vector dEv^<g><i> 
   * for every intermediate time stages t_{n+c_i} (with i=1,...,s)
   * and every Gauss point g
   *    dEv^<g><i>T = [ dE_11  dE_22  dE_33  2*dE_12  2*dE_23  2*dE_31 ] */
  ARRAY4D dvicstn;
  /* back stress rate vector dAv^<g><i> 
   * for every intermediate time stages t_{n+c_i} (with i=1,...,s)
   * and every Gauss point g
   *    dAv^<g><i>T = [ dA_11  dA_22  dA_33  dA_12  dA_23  dA_31 ] */
  ARRAY4D dbacsts;
  /* update vector for MIV iterative increments
   * (stored as Fortranesque vector 2*NUMSTR_SOLID3x1)
   *             [ kvv  kva ]^{-1}   [ res^v  ]
   *    kvarva = [ kav  kaa ]      . [ res^al ] */
  ARRAY kvarva;
  /* update matrix for MIV iterative increments
   * (stored as Fortranesque vector (2*NUMSTR_SOLID3*NUMSTR_SOLID3)x1)
   *              [ kvv  kva ]^{-1}   [ kve ]
   *    kvakvae = [ kav  kaa ]      . [ kae ] */
  ARRAY kvakvae;
} SO3_MIV_ROBINSON;

/*----------------------------------------------------------------------*/
/*!
\brief Definition of SOLID3 type holding SOLID3 element properties

\author mf
\date 10/06
*/
typedef struct _SOLID3
{
  SO3_KINEMATICS kintype;  /* type of kinematics */
  SO3_STRESSOUT stresstype;  /* output type of stress */

  /* number of Gauss points as obtained at read-in
   * hexahedra:
   *    gpnum[0] : in r-direction : 1,2,3,4,5,6 : read-in
   *    gpnum[1] : in s-direction : 1,2,3,4,5,6 : read-in/set
   *    gpnum[2] : in t-direction : 1,2,3,4,5,6 : read-in/set
   * tetrahedra:
   *    gpnum[0] : total number of GPs in domain: 1,4,5 : read-in
   *    gpnum[1] : total number of GPs on sides: 1,3,4,6 : read-in/set
   *    gpnum[2] : total number of GPs on edges: 1,2,3,4,5,6 : read-in/set */
  INT gpnum[NDIM_SOLID3];
  INT gpintc[NDIM_SOLID3];
  /* total number of GPs in domain */
  INT gptot;

  /* coordinates of Gauss points in material frame */
  ARRAY gpco_xyz;

  /* stress vector at Gauss points or nodes
   *    stress_gpxyz : at Gauss points in global XYZ-components
   *    stress_gprst : at Gauss points in parameter rst-compoonents
   *    stress_gp123 : at Gauss points modulus and angles vs XYZ-direct. */
  ARRAY stress_gpxyz;
  ARRAY stress_gprst;
  ARRAY stress_gp123;
  ARRAY stress_ndxyz;
  ARRAY stress_ndrst;
  ARRAY stress_nd123;
  /* ARRAY stress_ndeqv; */

  /* thermo-structure-interaction */
#ifdef D_TSI
  TSI_COUPTYP tsi_couptyp;  /* TSI coupling type */
  ELEMENT *therm_ele;  /* pointer to conforming THERM3 element */
#endif

  /* Robinson material, only with TSI */
#ifdef D_TSI
  SO3_MIV_ROBINSON* miv_rob;
#endif

} SOLID3;






/*======================================================================*/
/* Declarations of functions in SOLID3 directory
 * Order: Firstly, alphabetically list file names, secondly, list
 *        functions in file according to appearance */


/*----------------------------------------------------------------------*/
/* file so3_bop.c */
void so3_bop_lin(INT enod,
                 DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                 DOUBLE xji[NUMDOF_SOLID3][NUMDOF_SOLID3],
                 DOUBLE boplin[NUMSTR_SOLID3][MAXDOF_SOLID3]);
void so3_bop(ELEMENT *ele,
             INT enod,
             DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
             DOUBLE xji[NUMDOF_SOLID3][NUMDOF_SOLID3],
             DOUBLE defgrd[NUMDOF_SOLID3][NUMDOF_SOLID3],
             DOUBLE bopn[MAXDOF_SOLID3][NUMDOF_SOLID3],
             DOUBLE bop[NUMSTR_SOLID3][MAXDOF_SOLID3]);

/*---------------------------------------------------------------------*/
/* file so3_cfg.c */
void so3_cfg_chkdef();
void so3_cfg_init(SO3_DATA *data);
void so3_cfg_noderst(ELEMENT *ele,
                     SO3_DATA *data,
                     INT inode,
                     DOUBLE *rst);
#ifdef TEST_SOLID3
void so3_cfg_test(SO3_DATA *data);
#endif

/*----------------------------------------------------------------------*/
/* file so3_def.c */
void so3_def_grad(INT enod,
                  DOUBLE edis[MAXNOD_SOLID3][NDIM_SOLID3],
                  DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                  DOUBLE xji[NDIM_SOLID3][NDIM_SOLID3],
                  DOUBLE disgrdv[NUMDFGR_SOLID3],
                  DOUBLE defgrd[NDIM_SOLID3][NDIM_SOLID3]);

/*----------------------------------------------------------------------*/
/* file so3_gid.c */
void so3_gid_init(ELEMENT *actele,
                  GIDSET *actgid);
void so3_gid_msh(FIELD *actfield,
                 INT ldis,
                 GIDSET *actgid,
                 FILE *out);
void so3_gid_msh_val(FIELD *actfield,
                     CHAR *gidname,
                     INT ldis,
                     GIDSET *actgid,
                     CHAR *disname,
                     CHAR *volname,
                     INT nelenod,
                     FILE *out);
void so3_gid_gpset(INT jdis,
                   GIDSET *actgid,
                   FILE *out);
void so3_gid_gpset_int(GIDSET *actgid,
                       INT jdis,
                       CHAR *gpname,
                       CHAR *gidname,
                       CHAR *volname,
                       INT ngauss,
                       FILE *out);
void so3_gid_dom(FIELD *actfield,
                 INT disnum,
                 GIDSET *actgid,
                 FILE *out);
void so3_gid_dom_val(FIELD *actfield,
                     INT disnum,
                     CHAR *gidname,
                     GIDSET *actgid,
                     INT nelenod,
                     INT ngauss,
                     FILE *out);
void so3_gid_stress(CHAR resstring[],
                    FIELD *actfield,
                    INT disnum,
                    INT step,
                    GIDSET *actgid,
                    FILE *out);
void so3_gid_stress_gp(FIELD *actfield,
                       INT disnum,
                       ELEMENT *actele,
                       GIDSET *actgid,
                       SO3_STRESSOUT stresstype,
                       CHAR *resultname,
                       INT step,
                       CHAR *resulttype,
                       CHAR *resultplace,
                       CHAR *gpset,
                       INT ncomp,
                       CHAR *componentnames[],
                       INT nelenod,
                       INT *gperm,
                       INT ngauss,
                       FILE *out);

/*----------------------------------------------------------------------*/
/* file so3_inp.c */
void so3_inp(ELEMENT *ele);

/*----------------------------------------------------------------------*/
/* file so3_int.c */
void so3_int_fintstifmass(CONTAINER *container,
                          ELEMENT *ele,
                          SO3_GPSHAPEDERIV *gpshade,
                          MATERIAL *mat,
                          ARRAY *eforc_global,
                          ARRAY *estif_global,
                          ARRAY *emass_global);
void so3_int_fintcont(INT neledof,
                      DOUBLE bop[NUMSTR_SOLID3][MAXDOF_SOLID3],
                      DOUBLE stress[NUMSTR_SOLID3],
                      DOUBLE fac,
                      DOUBLE *intfor);
void so3_int_stiffeu(INT neledof,
                     DOUBLE bop[NUMSTR_SOLID3][MAXDOF_SOLID3],
                     DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3],
                     DOUBLE fac,
                     DOUBLE **estif);
void so3_int_stiffgeo(INT enod,
                      DOUBLE bopn[MAXDOF_SOLID3][NUMDOF_SOLID3],
                      DOUBLE stress[NUMSTR_SOLID3],
                      DOUBLE fac,
                      DOUBLE **estif);
void so3_int_mass(MATERIAL *mat,
                  INT nnod,
                  DOUBLE shape[MAXNOD_SOLID3],
                  DOUBLE fac,
                  DOUBLE **emass);

/*----------------------------------------------------------------------*/
/* file so3_intg.c */
void so3_intg_eleinp(ELEMENT *actele,
                     INT *ierr);
void so3_intg_init(SO3_DATA *data);

/*----------------------------------------------------------------------*/
/* file so3_iv.c */
void so3_iv_updreq(ELEMENT* ele,
                   MATERIAL* actmat,
                   INT* updreq);
void so3_iv_upditer(CONTAINER* container,
                    ELEMENT* ele,
                    SO3_GPSHAPEDERIV* gpshade,
                    MATERIAL* mat);
void so3_iv_updincr(const CONTAINER* container,
                    ELEMENT* ele,
                    const MATERIAL* mat);

/*----------------------------------------------------------------------*/
/* file so3_load.c */
void so3_load(ELEMENT *ele,
              SO3_DATA *data,
              SO3_GPSHAPEDERIV *gpshade,
              INT imyrank,
              ARRAY *eforc_global);

/*----------------------------------------------------------------------*/
/* file so3_load_line.c */
void so3_load_line_int(ELEMENT* ele,
                       SO3_DATA* data,
                       DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                       DOUBLE exs[MAXNOD_SOLID3][NDIM_SOLID3],
                       const DOUBLE timen,
                       INT ngline,
                       GLINE* gsurf[MAXEDG_SOLID3],
                       DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3]);
void so3_load_line_valh(ELEMENT* ele,
                        SO3_DATA* data,
                        INT igline,
                        GLINE* gline,
                        DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE exs[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE gpc[NDIM_SOLID3],
                        DOUBLE shape[MAXNOD_SOLID3],
                        DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE fac,
                        DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3]);
void so3_load_line_valt(ELEMENT* ele,
                        SO3_DATA* data,
                        INT igline,
                        GLINE* gline,
                        DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE exs[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE gpc[NDIM_SOLID3],
                        DOUBLE shape[MAXNOD_SOLID3],
                        DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE fac,
                        DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3]);

/*----------------------------------------------------------------------*/
/* file so3_load_surf.c */
void so3_load_surf_int(ELEMENT* ele,
                       SO3_DATA* data,
                       DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                       DOUBLE exc[MAXNOD_SOLID3][NDIM_SOLID3],
                       const DOUBLE timen,
                       INT ngsurf,
                       GSURF* gsurf[MAXSID_SOLID3],
                       DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3]);
void so3_load_surf_valh(ELEMENT* ele,
                        SO3_DATA* data,
                        INT igsurf,
                        GSURF* gsurf,
                        DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE exs[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE gpc[NDIM_SOLID3],
                        DOUBLE shape[MAXNOD_SOLID3],
                        DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE fac,
                        DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3]);
void so3_load_surf_valt(ELEMENT* ele,
                        SO3_DATA* data,
                        INT igsurf,
                        GSURF* gsurf,
                        DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE exs[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE gpc[NDIM_SOLID3],
                        DOUBLE shape[MAXNOD_SOLID3],
                        DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                        DOUBLE fac,
                        DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3]);

/*----------------------------------------------------------------------*/
/* file so3_load_vol.c */
void so3_load_vol_int(ELEMENT* ele,
                      SO3_GPSHAPEDERIV *gpshade,
                      DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                      const DOUBLE timen,
                      GVOL* gvol,
                      DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3]);
void so3_load_vol_val(ELEMENT* ele,
                      INT nelenod,
                      DOUBLE shape[MAXNOD_SOLID3],
                      DOUBLE fac,
                      DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3]);

/*----------------------------------------------------------------------*/
/* file so3_main.c */
void solid3(PARTITION *actpart,
	    INTRA *actintra,
	    ELEMENT *ele,
	    ARRAY *estif_global,
	    ARRAY *emass_global,
	    ARRAY *intforce_global,
	    CALC_ACTION *action,
	    CONTAINER *container);

/*----------------------------------------------------------------------*/
/* file so3_mat.c */
void so3_mat_init(PARTITION* part,
                  const MATERIAL* mat);
void so3_mat_final(PARTITION* part,  /*!< partition */
                   const MATERIAL* mat);  /*!< material */
void so3_mat_sel(CONTAINER* container,
                 ELEMENT* ele,
                 MATERIAL* mat,
                 SO3_GPSHAPEDERIV* gpshade,
                 INT ip,
                 SO3_GEODEFSTR* gds,
                 DOUBLE stress[NUMSTR_SOLID3],
                 DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mat_stress(CONTAINER *container,
                    ELEMENT *ele,
                    MATERIAL *mat,
                    INT ip,
                    SO3_GEODEFSTR *gds,
                    DOUBLE stress[NUMSTR_SOLID3],
                    DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mat_density(MATERIAL *mat, 
                     DOUBLE *density);
void so3_mat_mivupdreq(ELEMENT* ele,
                       MATERIAL* mat,
                       INT* updreq);
void so3_mat_mivupditer(const CONTAINER* container,
                        ELEMENT* ele,
                        const MATERIAL* mat,
                        const INT ip,
                        const SO3_GEODEFSTR* gds,
                        const DOUBLE epsii[NUMSTR_SOLID3]);
void so3_mat_mivupdincr(const CONTAINER* container,
                        ELEMENT* ele,
                        const MATERIAL* mat);

/*----------------------------------------------------------------------*/
/* file so3_mat_robinson.c */
#ifdef D_TSI
void so3_mat_robinson_init(ELEMENT* ele);
void so3_mat_robinson_final(ELEMENT* ele);
void so3_mat_robinson_prmbytmpr(const MAT_PARAM_MULT prm,
                                const DOUBLE tmpr,
                                DOUBLE* prmbytempr);
void so3_mat_robinson_elmat(const VP_ROBINSON* mat_robin,
                            const DOUBLE tem,
                            DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mat_robinson_sel(const CONTAINER* container,
                          ELEMENT* ele,
                          const VP_ROBINSON* mat_robin,
                          SO3_GPSHAPEDERIV* gpshade,
                          const INT ip,
                          SO3_GEODEFSTR* gds,
                          DOUBLE stress[NUMSTR_SOLID3],
                          DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mat_robinson_mivupdreq(ELEMENT* ele,
                                VP_ROBINSON* mat_robin,
                                INT* updreq);
void so3_mat_robinson_mivupditer(const CONTAINER* container,
                                 ELEMENT* ele,
                                 const VP_ROBINSON* mat_robin,
                                 const INT ip,
                                 const DOUBLE epsii[NUMSTR_SOLID3]);
void so3_mat_robinson_mivupdincr(const CONTAINER* container,
                                 ELEMENT* ele,
                                 const VP_ROBINSON* mat_robin);
void so3_mat_robinson_stress(const CONTAINER* container,
                             const ELEMENT* ele,
                             const VP_ROBINSON* mat_robin,
                             const INT ip,
                             const SO3_GEODEFSTR* gds,
                             DOUBLE stress[NUMSTR_SOLID3],
                             DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3]);
#endif

/*----------------------------------------------------------------------*/
/* file so3_mat_robinson_be.c */
#ifdef D_TSI
void so3_mat_robinson_be_init(SOLID3* actso3);
void so3_mat_robinson_be_final(SOLID3* actso3);
void so3_mat_robinson_be_sel(const CONTAINER* container,
                             ELEMENT* ele,
                             const VP_ROBINSON* mat_robin,
                             SO3_GPSHAPEDERIV* gpshade,
                             const INT ip,
                             SO3_GEODEFSTR* gds,
                             DOUBLE stress[NUMSTR_SOLID3],
                             DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mat_robinson_be_rvscstn(ELEMENT* ele,
                                 const VP_ROBINSON* mat_robin,
                                 const DOUBLE dt,
                                 const DOUBLE tmpr,
                                 const DOUBLE vscstn[NUMSTR_SOLID3],
                                 const DOUBLE vscstnn[NUMSTR_SOLID3],
                                 const DOUBLE devstsn[NUMSTR_SOLID3],
                                 const DOUBLE ovrstsn[NUMSTR_SOLID3],
                                 INT* vscstns,
                                 DOUBLE vscstnr[NUMSTR_SOLID3],
                                 DOUBLE kve[NUMSTR_SOLID3][NUMSTR_SOLID3],
                                 DOUBLE kvv[NUMSTR_SOLID3][NUMSTR_SOLID3],
                                 DOUBLE kva[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mat_robinson_be_rbcksts(ELEMENT* ele,
                                 const VP_ROBINSON* mat_robin,
                                 const DOUBLE dt,
                                 const DOUBLE tmpr,
                                 const DOUBLE vscstn[NUMSTR_SOLID3],
                                 const DOUBLE vscstnn[NUMSTR_SOLID3],
                                 const DOUBLE devstsn[NUMSTR_SOLID3],
                                 const DOUBLE bacsts[NUMSTR_SOLID3],
                                 const DOUBLE bacstsn[NUMSTR_SOLID3],
                                 INT* bckstss,
                                 DOUBLE bckstsr[NUMSTR_SOLID3],
                                 DOUBLE kae[NUMSTR_SOLID3][NUMSTR_SOLID3],
                                 DOUBLE kav[NUMSTR_SOLID3][NUMSTR_SOLID3],
                                 DOUBLE kaa[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mat_robinson_be_red(DOUBLE stress[NUMSTR_SOLID3], 
                             DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             DOUBLE kev[NUMSTR_SOLID3][NUMSTR_SOLID3], 
                             DOUBLE kea[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             const DOUBLE vscstnr[NUMSTR_SOLID3],
                             DOUBLE kve[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             DOUBLE kvv[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             DOUBLE kva[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             const DOUBLE bckstsr[NUMSTR_SOLID3],
                             DOUBLE kae[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             DOUBLE kav[NUMSTR_SOLID3][NUMSTR_SOLID3], 
                             DOUBLE kaa[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             DOUBLE* kvarva,
                             DOUBLE* kvakvae);
void so3_mat_robinson_be_mivupdreq(ELEMENT* ele,
                                   VP_ROBINSON* mat_robin,
                                   INT* updreq);
void so3_mat_robinson_be_mivupditer(const CONTAINER* container,
                                    ELEMENT* ele,
                                    const VP_ROBINSON* mat_robin,
                                    const INT ip,
                                    const DOUBLE epsii[NUMSTR_SOLID3]);
void so3_mat_robinson_be_mivupdincr(ELEMENT* ele,
                                    const VP_ROBINSON* mat_robin);
void so3_mat_robinson_be_stress(const CONTAINER* container,
                                const ELEMENT* ele,
                                const VP_ROBINSON* mat_robin,
                                const INT ip,
                                const SO3_GEODEFSTR* gds,
                                DOUBLE stress[NUMSTR_SOLID3],
                                DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3]);
#endif

/*----------------------------------------------------------------------*/
/* file so3_mat_robinson_fb4.c */
#ifdef D_TSI
void so3_mat_robinson_fb4_init(SOLID3* actso3);
void so3_mat_robinson_fb4_final(SOLID3* actso3);
void so3_mat_robinson_fb4_sel(const CONTAINER* container,
                              const ELEMENT* ele,
                              const VP_ROBINSON* mat_robin,
                              const INT ip,
                              const SO3_GEODEFSTR* gds,
                              DOUBLE stress[NUMSTR_SOLID3],
                              DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mat_robinson_fb4_stnvscrat(const VP_ROBINSON* mat_robin,
                                    const DOUBLE tmpr,
                                    const DOUBLE stsdev[NUMSTR_SOLID3],
                                    const DOUBLE stsovr[NUMSTR_SOLID3],
                                    DOUBLE dstnvsc[NUMSTR_SOLID3]);
void so3_mat_robinson_fb4_stsbckrat(const VP_ROBINSON* mat_robin,
                                    const DOUBLE tmpr,
                                    const DOUBLE stsdev[NUMSTR_SOLID3],
                                    const DOUBLE stsbck[NUMSTR_SOLID3],
                                    const DOUBLE dstnvsc[NUMSTR_SOLID3],
                                    DOUBLE dstsbck[NUMSTR_SOLID3]);
void so3_mat_robinson_fb4_mivupdincr(ELEMENT* ele,
                                 const VP_ROBINSON* mat_robin);
void so3_mat_robinson_fb4_stress(const CONTAINER* container,
                                 const ELEMENT* ele,
                                 const VP_ROBINSON* mat_robin,
                                 const INT ip,
                                 const SO3_GEODEFSTR* gds,
                                 DOUBLE stress[NUMSTR_SOLID3],
                                 DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3]);
#endif

/*----------------------------------------------------------------------*/
/* file so3_mat_stvenant.c */
void so3_mat_stvenant_sel(CONTAINER* container,
                          ELEMENT* ele,
                          MATERIAL* mat,
                          INT ip,
                          SO3_GEODEFSTR* gds,
                          DOUBLE stress[NUMSTR_SOLID3],
                          DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3]);

/*----------------------------------------------------------------------*/
/* file so3_metr.c */
void so3_metr_jaco(ELEMENT *ele,
                   INT enod,
                   DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                   DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                   INT flag,
                   DOUBLE xjm[NDIM_SOLID3][NDIM_SOLID3],
                   DOUBLE *det,
                   DOUBLE xji[NDIM_SOLID3][NDIM_SOLID3]);
void so3_metr_surf(ELEMENT *ele, 
                   INT nelenod, 
                   DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                   DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                   DOUBLE sidredm[DIMSID_SOLID3][NDIM_SOLID3],
                   DOUBLE gamtt[DIMSID_SOLID3][NDIM_SOLID3],
                   DOUBLE *metr);
void so3_metr_line(ELEMENT *ele, 
                   INT nelenod, 
                   DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                   DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3], 
                   DOUBLE linredv[NDIM_SOLID3],
                   DOUBLE gamtt[NDIM_SOLID3],
                   DOUBLE *metr);
void so3_metr_rot(DOUBLE xjm[NDIM_SOLID3][NDIM_SOLID3],
                  DOUBLE xrm[NDIM_SOLID3][NDIM_SOLID3],
                  DOUBLE xrvm[NUMSTR_SOLID3][NUMSTR_SOLID3],
                  DOUBLE xrvi[NUMSTR_SOLID3][NUMSTR_SOLID3]);

/*----------------------------------------------------------------------*/
/* file so3_out.c */
void so3_out_stress(ELEMENT *actele,
                    FILE *out);

/*----------------------------------------------------------------------*/
/* file so3_shape.c */
void so3_shape_deriv(DIS_TYP typ,
                     DOUBLE r,
                     DOUBLE s,
                     DOUBLE t,
                     INT option,
                     DOUBLE shape[MAXNOD_SOLID3],
                     DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3]);
void so3_shape_gpshade_init(SO3_GPSHAPEDERIV* so3_gpshade);
void so3_shape_gpshade(ELEMENT *ele,
                       SO3_DATA *data,
                       SO3_GPSHAPEDERIV *so3_gpshade);
#ifdef TEST_SOLID3
void so3_shape_test(SO3_DATA *data);
#endif
#ifdef TEST_SOLID3
void so3_shape_test_shp(SO3_DATA *data,
                        DIS_TYP dis,
                        CHAR *text, 
                        INT nelenod,
                        FILE *filetest);
#endif
#ifdef TEST_SOLID3
void so3_shape_test_drv(SO3_DATA *data,
                        DIS_TYP dis,
                        CHAR *text, 
                        INT nelenod,
                        FILE *filetest);
#endif

/*----------------------------------------------------------------------*/
/* file so3_strain.c */
void so3_strain_lin(ELEMENT *ele,
                    DOUBLE disgrdv[NUMDFGR_SOLID3],
                    DOUBLE strain[NUMSTR_SOLID3]);
void so3_strain_gl(ELEMENT *ele,
                   DOUBLE disgrdv[NUMDFGR_SOLID3],
                   DOUBLE strain[NUMSTR_SOLID3]);

/*----------------------------------------------------------------------*/
/* file so3_stress.c */
void so3_stress_init(PARTITION *actpart);
void so3_stress_final(PARTITION *actpart);
void so3_stress(CONTAINER *cont,
                ELEMENT *ele,
                SO3_DATA *data,
                SO3_GPSHAPEDERIV *gpshade,
                INT imyrank,
                MATERIAL *mat);
void so3_stress_rst(DOUBLE xrvi[NUMSTR_SOLID3][NUMSTR_SOLID3],
                    DOUBLE stress[NUMSTR_SOLID3],
                    DOUBLE *stressrst);
void so3_stress_123(DOUBLE stress[NUMSTR_SOLID3],
                    DOUBLE *stress123);
void so3_stress_extrpol(ELEMENT *ele,
                        SO3_DATA *data,
                        INT ngauss,
                        DOUBLE **stressgp,
                        DOUBLE rst[NDIM_SOLID3],
                        DOUBLE stressnd[NUMSTR_SOLID3]);

/*----------------------------------------------------------------------*/
/* file so3_tns3.c */
void so3_tns3_zero(DOUBLE ot[NDIM_SOLID3][NDIM_SOLID3]);
void so3_tns3_id(DOUBLE it[NDIM_SOLID3][NDIM_SOLID3]);
void so3_tns3_tr(DOUBLE at[NDIM_SOLID3][NDIM_SOLID3],
                 DOUBLE *tr);
void so3_tns3_det(DOUBLE at[NDIM_SOLID3][NDIM_SOLID3],
                  DOUBLE *det);
void so3_tns3_inva(DOUBLE at[NDIM_SOLID3][NDIM_SOLID3],
                   DOUBLE *ai,
                   DOUBLE *aii,
                   DOUBLE *aiii);
void so3_tns3_dotprod(DOUBLE at[NDIM_SOLID3][NDIM_SOLID3],
                      DOUBLE bt[NDIM_SOLID3][NDIM_SOLID3],
                      DOUBLE ct[NDIM_SOLID3][NDIM_SOLID3]);
void so3_tns3_dotprod_tl(DOUBLE at[NDIM_SOLID3][NDIM_SOLID3],
                         DOUBLE bt[NDIM_SOLID3][NDIM_SOLID3],
                         DOUBLE ct[NDIM_SOLID3][NDIM_SOLID3]);
void so3_tns3_dotprod_tr(DOUBLE at[NDIM_SOLID3][NDIM_SOLID3],
                         DOUBLE bt[NDIM_SOLID3][NDIM_SOLID3],
                         DOUBLE ct[NDIM_SOLID3][NDIM_SOLID3]);
void so3_tns3_plrdcmp(DOUBLE ft[NDIM_SOLID3][NDIM_SOLID3],
                      DOUBLE rt[NDIM_SOLID3][NDIM_SOLID3],
                      DOUBLE ut[NDIM_SOLID3][NDIM_SOLID3],
                      DOUBLE vt[NDIM_SOLID3][NDIM_SOLID3]);
void so3_tns3_tsym2v(DOUBLE at[NDIM_SOLID3][NDIM_SOLID3],
                     DOUBLE av[NUMSTR_SOLID3]);
void so3_tns3_v2tsym(DOUBLE av[NUMSTR_SOLID3],
                     DOUBLE at[NDIM_SOLID3][NDIM_SOLID3]);
void so3_tns3_spcdcmp(DOUBLE at[NDIM_SOLID3][NDIM_SOLID3],
                      DOUBLE ew[NDIM_SOLID3],
                      INT *err);
void so3_tns3_symspcdcmp_jit(DOUBLE at[NDIM_SOLID3][NDIM_SOLID3],
                             DOUBLE itertol,
                             INT itermax,
                             DOUBLE ew[NDIM_SOLID3],
                             DOUBLE ev[NDIM_SOLID3][NDIM_SOLID3],
                             INT* err);
void so3_tns3_norm2(const DOUBLE av[NDIM_SOLID3],
                    DOUBLE* norm);
void so3_tns3_unitvct(DOUBLE av[NDIM_SOLID3]);
void so3_tns3_crsprd(const DOUBLE av[NDIM_SOLID3],
                     const DOUBLE bv[NDIM_SOLID3],
                     DOUBLE cv[NDIM_SOLID3]);
void so3_tns3_unrm(const DOUBLE av[NDIM_SOLID3],
                   const DOUBLE bv[NDIM_SOLID3],
                   const INT swp,
                   DOUBLE cv[NDIM_SOLID3]);

/*----------------------------------------------------------------------*/
/* file so3_tsi.c */
#ifdef D_TSI
void so3_tsi_temper(const CONTAINER *container,
                    const ELEMENT *ele,
                    const DOUBLE r, 
                    const DOUBLE s,
                    const DOUBLE t,
                    DOUBLE *temper);
#endif

/*----------------------------------------------------------------------*/
/* file so3_mv6.c */
void so3_mv6_v_zero(DOUBLE iv[NUMSTR_SOLID3]);
void so3_mv6_v_id(DOUBLE iv[NUMSTR_SOLID3]);
void so3_mv6_v_ass(const DOUBLE av[NUMSTR_SOLID3],
                  DOUBLE bv[NUMSTR_SOLID3]);
void so3_mv6_v_assscl(const DOUBLE scl,
                     const DOUBLE av[NUMSTR_SOLID3],
                     DOUBLE bv[NUMSTR_SOLID3]);
void so3_mv6_v_updscl(const DOUBLE scl,
                     const DOUBLE av[NUMSTR_SOLID3],
                     DOUBLE bv[NUMSTR_SOLID3]);
void so3_mv6_v_tr(const DOUBLE av[NUMSTR_SOLID3],
                 DOUBLE *tr);
void so3_mv6_v_det(const DOUBLE av[NUMSTR_SOLID3],
                  DOUBLE *det);
void so3_mv6_v_dev(const DOUBLE av[NUMSTR_SOLID3],
                  DOUBLE adev[NUMSTR_SOLID3]);
void so3_mv6_v_sub(const DOUBLE av[NUMSTR_SOLID3],
                   const DOUBLE bv[NUMSTR_SOLID3],
                   DOUBLE cv[NUMSTR_SOLID3]);
void so3_mv6_v_dblctr(const DOUBLE av[NUMSTR_SOLID3],
                      const DOUBLE bv[NUMSTR_SOLID3],
                      DOUBLE* prd);
void so3_mv6_v_assmvp(DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                      const DOUBLE bv[NUMSTR_SOLID3],
                      DOUBLE cv[NUMSTR_SOLID3]);
void so3_mv6_v2_assscl(const DOUBLE scl,
                       const DOUBLE av[NUMSTR_SOLID3],
                       DOUBLE bv[NUMSTR_SOLID3]);
void so3_mv6_v05_assscl(const DOUBLE scl,
                        const DOUBLE av[NUMSTR_SOLID3],
                        DOUBLE bv[NUMSTR_SOLID3]);
void so3_mv6_v05_updvtov05(DOUBLE av[NUMSTR_SOLID3]);
void so3_mv6_m_zero(DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_id(DOUBLE id[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_idscl(const DOUBLE scale,
                     DOUBLE id[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_one(DOUBLE om[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_updonescl(const DOUBLE scl,
                         DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_ass(const DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                   DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_assscl(const DOUBLE scl,
                      DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                      DOUBLE mv[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_updscl(const DOUBLE scl,
                      DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                      DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_sub(DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                   DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3],
                   DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_assdyd(const DOUBLE av[NUMSTR_SOLID3],
                      const DOUBLE bv[NUMSTR_SOLID3],
                      DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_upddydscl(const DOUBLE scl,
                         const DOUBLE av[NUMSTR_SOLID3],
                         const DOUBLE bv[NUMSTR_SOLID3],
                         DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_mprd(DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                    DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3],
                    DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_mprdscl(const DOUBLE scl,
                       DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                       DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3],
                       DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_updmprd(DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                       DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3],
                       DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m_updmprdscl(const DOUBLE scl,
                          DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                          DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3],
                          DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m2_idscl(const DOUBLE scale,
                      DOUBLE im[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m2_updmtom2(DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3]);
void so3_mv6_m05_idscl(const DOUBLE scale,
                       DOUBLE im[NUMSTR_SOLID3][NUMSTR_SOLID3]);

/*----------------------------------------------------------------------*/
#endif /*end of #ifdef D_SOLID3 */
/*! @} (documentation module close) */
