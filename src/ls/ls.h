/*!----------------------------------------------------------------------
\file
\brief ls.h

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_LS

#ifndef LS_H
#define LS_H

/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/



/*!----------------------------------------------------------------------
\brief 2D level set interface data

<pre>                                                            irhan 05/04
This structure contains data related to interface
</pre>

*----------------------------------------------------------------------*/
typedef struct _LS_INT_DATA
{
  DOUBLE               p[2][2];      /* start and end  points of interface */
  DOUBLE               pd[2];        /* intersection with diagonal */
  INT                  edge[2];      /* start and end edges */
  INT                  polycon[2][3];/* polygon connectivity */
  INT                  is_diagcut;   /* is diagonal cut */
  INT                  reconstruct;  /* flag for reconstruction */
  struct _ELEMENT     *nbor;         /* element neighbor to eedge */
  struct _ELEMENT     *nbor_s;       /* element neighbor to sedge */
} LS_INT_DATA;



/*!----------------------------------------------------------------------
\brief 2D level set polygon data

<pre>                                                            irhan 05/04
This structure contains data related to polygonization
</pre>
*----------------------------------------------------------------------*/
typedef struct _LS_POLY_DATA
{
  INT                  ind[2];               /* index for refined integration */
  INT                  polygonmat[2][7];     /* material properties for each subpolygon */
  DOUBLE               polygonwgt[2][7];     /* weight corresponding to each subpolygon */
  DOUBLE               polygonGP[2][2][7];   /* Gauss points corresponding to each subpolygon */
} LS_POLY_DATA;



/*!----------------------------------------------------------------------
\brief 2D level set element

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

<pre>                                                          irhan 05/04
This structure contains all specific information for a 2D level set element
</pre>

*----------------------------------------------------------------------*/
typedef struct _LS2
{
  INT                  ntyp;         /* flag for element type */
  INT                  nGP[2];       /* # of GP in r and s directions */
  INT                  is_elcut;     /* flag to control element activation */
  INT                  is_sedge_set; /* flag to control element edge activation */
  INT                  is_elsearch;  /* flag to confine element search */
  INT                  nlayer;       /* corresponding layer number */
  INT                  nenode;       /* number of enriched nodes */
  INT                  enode[4];     /* local node numbers of enriched nodes */
  INT                  intcnt;       /* index to interface number involved */
  INT                  prncnt;       /* counter used in printing */
  INT                  rstcnt;       /* counter used in resetting */

  struct _LS_INT_DATA  intdata[2];   /* struct to store interface data */
  struct _LS_POLY_DATA polydata[2];  /* struct to store polygon data */
  struct _ELEMENT      *my_fluid;    /* pointer to fluid element associated */
} LS2;


/*!----------------------------------------------------------------------
\brief 2D level set dynamic data

<pre>                                                            irhan 05/04
This structure contains the necessary data to control time history analysis
for level set problem
</pre>

*----------------------------------------------------------------------*/
typedef struct _LS_DYNAMIC
{
  INT                     step;        /* the actual step */
  INT 	                  nstep;       /* number of time steps */
  INT 	                  nfstep;      /* number of fluid only time steps */
  INT 	                  nfreinit;    /* frequency of reinitialization */
  INT 	                  nfreprn;     /* frequency of print */
  INT  	                  itemax;      /* maximum num. of itn. at a certain time step */
  INT                     iop;         /* time integration scheme */
  INT	                  ite;         /* iteration scheme */
  INT                     itchk;       /* convergence check during nonlin. iteration */
  INT                     stchk;       /* flag for steady state check */
  INT                     itnrm;       /* norm for conv. check d. nonlin. iteration */
  INT 	                  init;	       /* initialization of level set profile */
  INT 	                  nif;         /* flag for evaluation of time rhs */
  INT 	                  nii;         /* flag for evaluation of iteration rhs */
  DOUBLE                  time;	       /* actual time */
  DOUBLE                  dt;	       /* time increment */
  DOUBLE                  maxtime;     /* maximum time */
  DOUBLE                  ittol;       /* tolerance to check convergence */
  DOUBLE                  sttol;       /* tolerance to check convergence */
  DOUBLE                  theta;       /* integration parameter for one step theta scheme */
  struct _LS_GEN_DATA    *lsdata;      /* struct to store data specific to levels set formulation */
} LS_DYNAMIC;



/*!----------------------------------------------------------------------
\brief 2D level set general data

<pre>                                                            irhan 05/04
This structure contains data specific to level set problem
</pre>

*----------------------------------------------------------------------*/
typedef struct _LS_GEN_DATA
{
  INT        numcirc;                /* number of circular interfaces */
  INT        numline;                /* number of line interfaces */
  INT        is_sharp;               /* flag for initial level set profile */
  INT        setvel;                 /* flag to set velocity field (by user, or by fluid) */
  INT        flag_vel;               /* if setvel==byuser => flag for description of velocity profile */

  DOUBLE     xc1;                    /* data for circular interfaces */
  DOUBLE     yc1;
  DOUBLE     rad1;
  DOUBLE     xc2;
  DOUBLE     yc2;
  DOUBLE     rad2;
  DOUBLE     xc3;
  DOUBLE     yc3;
  DOUBLE     rad3;

  DOUBLE     xs1;                    /* data for first line interfaces */
  DOUBLE     ys1;
  DOUBLE     xe1;
  DOUBLE     ye1;
  DOUBLE     xs2;
  DOUBLE     ys2;
  DOUBLE     xe2;
  DOUBLE     ye2;
  DOUBLE     xs3;
  DOUBLE     ys3;
  DOUBLE     xe3;
  DOUBLE     ye3;

  INT 	     isstab;              /* flag to control stabilization */
  INT 	     localization;        /* flag to control localization (localization not implemented yet!) */
  INT 	     numlayer;            /* flag to control localization region */

  INT 	     reinitialization;    /* flag to control reinitialization */
  INT 	     numreinit;           /* number of reinitialization */
  INT 	     reinitflag;          /* initialization flag */
  DOUBLE     rdt;                 /* dt for reinitialization */
  INT        algo;                /* treatment of velocity (imp. or exp.) */
  DOUBLE     epsilon;             /* smoothing parameter for sign function */
  INT 	     anchor;              /*  */

  INT        boundary_on_off;     /* treatment of boundary */
  INT        reconstruct_on_off;  /* flag used in reconstruction phase */
  INT        print_on_off;        /* flag used to control print out */

  INT        probdescr;           /* ( = 0 ) => no description */
                                  /* ( = 1 ) => bubble problem */
                                  /* ( = 2 ) => breaking dam   */
  INT        searchband;          /* maximum number of layer. it is used for element search in interface construction */
} LS_GEN_DATA;



/*!----------------------------------------------------------------------
\brief 2D level set update

<pre>                                                            irhan 05/04
This structure contains interface information
</pre>

*----------------------------------------------------------------------*/
typedef struct _LS_UPDATE
{
  INT 	      *nael;              /* number of active elements */
  INT 	       nand;              /* number of active nodes */
  ELEMENT     *first_in_list;     /* first element in the list */
  INT          numinterface;      /* number of interfaces */
  INT         *typeinterface;     /* type of interface (open or closed) */
} LS_UPDATE;



/*!----------------------------------------------------------------------
\brief 2D level set element integration data

<pre>                                                            irhan 05/04
This structure contains the coordinates and weights for numerical integration
of a 2D level set element
</pre>

*----------------------------------------------------------------------*/
typedef struct _LS2_INTG_DATA
{
  /* for quadrilateral element */
  DOUBLE     wgtq[3][3];
  DOUBLE     xgq[3][3];
  /* for triangular element */
  DOUBLE     wgtt[3][3];
  DOUBLE     xgtr[3][3];
  DOUBLE     xgts[3][3];
} LS2_INTG_DATA;



/*!----------------------------------------------------------------------
\brief 2D level set lsflag

<pre>                                                            irhan 05/04
This enumerator contains execution flags used in level set formulation
</pre>

*----------------------------------------------------------------------*/
typedef enum _LSFLAG
{
  ls_initphase,
  ls_updtphase,
  ls_finaphase
} LSFLAG;



/*!----------------------------------------------------------------------
\brief 2D level set frontlsflag

<pre>                                                            irhan 05/04
This enumerator contains execution flags used in level set interface update
</pre>

*----------------------------------------------------------------------*/
typedef enum _FRONTLSFLAG
{
  front_ls_set,
  front_ls_init,
  front_ls_initmat,
  front_ls_updt,
  front_ls_updtmat,
  front_ls_activate,
  front_ls_localize,
  front_ls_polygonize,
  front_ls_write,
  front_ls_write_fld_soln,
  front_ls_finalize,
  front_ls_init_pressure_breaking_dam
} FRONTLSFLAG;
/*! @} (documentation module close)*/

#endif
#endif
