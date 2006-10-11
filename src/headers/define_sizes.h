/*!---------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/

#include "compile_settings.h"

/*----------------------------------------------------------------------*
 *                                                                      *
 * General definitions of max. sizes                                    *
 *                                                                      *
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | maximum number of nodes to an element                                |
 *----------------------------------------------------------------------*/
#ifndef MAXNOD
#define MAXNOD           (4)
#endif

/*----------------------------------------------------------------------*
 | maximum number of elements to a node                                 |
 *----------------------------------------------------------------------*/
#ifndef MAXELE
#define MAXELE           (4)
#endif


/*----------------------------------------------------------------------*
 | maximum number of dofs to a node -> can be more than 6 for shell9    |
 *----------------------------------------------------------------------*/
#ifndef MAXDOFPERNODE
#define MAXDOFPERNODE    (3)
#endif

/*----------------------------------------------------------------------*
 | maximum number of gaussian points in an element                      |
 *----------------------------------------------------------------------*/
#ifndef MAXGAUSS
#define MAXGAUSS         (4)
#endif

/*----------------------------------------------------------------------*
 | maximum number of fields                                             |
 *----------------------------------------------------------------------*/
#ifndef MAXFIELD
#define MAXFIELD         (3)
#endif

/*----------------------------------------------------------------------*
 | maximum number of fields                                             |
 *----------------------------------------------------------------------*/
#ifndef MAXDIS
#define MAXDIS           (2)
#endif

/*----------------------------------------------------------------------*
 | maximum numberof fields                                              |
 *----------------------------------------------------------------------*/
#ifndef MAXTIMECURVE
#define MAXTIMECURVE     (5)
#endif

/*----------------------------------------------------------------------*
 | maximum numberof restart records per element in dynamic and nonlinear|
 | static restarting from and to pss file                               |
 *----------------------------------------------------------------------*/
#ifndef MAXRECORDPERELE
#define MAXRECORDPERELE  (6)
#endif


/*
 * The maximum number of matrices per field. Needed for aztec's reuse
 * feature.
 */
#ifndef MAXNUMMATRICES
#define MAXNUMMATRICES  (10)
#endif

/*----------------------------------------------------------------------*
 *                                                                      *
 * element specific definitions of max. sizes                           *
 *                                                                      *
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | axishell                                                             |
 *----------------------------------------------------------------------*/
#ifndef MAXNOD_AXISHELL
#define MAXNOD_AXISHELL  (2)
#endif

/*----------------------------------------------------------------------*
 | beam3                                                                |
 *----------------------------------------------------------------------*/
#ifndef MAXNOD_BEAM3
#define MAXNOD_BEAM3     (3)
#endif

/*----------------------------------------------------------------------*
 | brick8                                                               |
 *----------------------------------------------------------------------*/
#ifndef MAXNOD_BRICK1
#define MAXNOD_BRICK1    (20)
#endif
#ifndef NUMDOF_BRICK1
#define NUMDOF_BRICK1    (3)
#endif

/*----------------------------------------------------------------------*
 | fluid                                                                |
 *----------------------------------------------------------------------*/
#ifndef MAXQINTC
#define MAXQINTC         (6)  /* max. number of implemented
                               * integration cases for QUADS */
#endif
#ifndef MAXQINTP
#define MAXQINTP         (6)  /* max. number of integration parameters
                               * for QUADS */
#endif
#ifndef MAXTINTC
#define MAXTINTC         (11) /* max. number of implemented
                               * integration cases for TRIS */
#endif
#ifndef MAXTINTP
#define MAXTINTP         (13) /* max. number of integration parameters
                               * for TRIS */
#endif
#ifndef FLUID_NUM_LD
#define FLUID_NUM_LD     (4)  /* number of lift and drag conditions */
#endif

#ifndef FLUID_NUM_BNODE
#define FLUID_NUM_BNODE  (10)  /* max. number of fluid boundary nodes with
                                * lift & drag condition or FSI coupling */
#endif

/*----------------------------------------------------------------------*
 | fluid2                                                               |
 *----------------------------------------------------------------------*/
#ifndef NUM_F2_VELDOF
#define NUM_F2_VELDOF    (2)  /* number of velocity dofs per node */
#endif
#ifndef NUMDOF_FLUID2
#define NUMDOF_FLUID2    (3)  /* number of dofs per node */
#endif
#ifndef MAXNOD_F2
#define MAXNOD_F2        (9)  /* max. number of nodes per element */
#endif

/*----------------------------------------------------------------------*
 | fluid3                                                               |
 *----------------------------------------------------------------------*/
#ifndef NUM_F3_VELDOF
#define NUM_F3_VELDOF    (3)  /* number of velocity dofs per node */
#endif
#ifndef NUMDOF_FLUID3
#define NUMDOF_FLUID3    (4)  /* number of dofs per node */
#endif
#ifndef MAXNOD_F3
#define MAXNOD_F3        (27) /* max. number of nodes per element */
#endif

/*----------------------------------------------------------------------*
 | ale3                                                               |
 *----------------------------------------------------------------------*/
#ifndef MAXNOD_ALE3
#define MAXNOD_ALE3      (8) /* max. number of nodes per element */
#endif

/*----------------------------------------------------------------------*
 | shell8                                                               |
 *----------------------------------------------------------------------*/
#ifndef MAXNOD_SHELL8
#define MAXNOD_SHELL8    (9)
#endif
#ifndef NUMDOF_SHELL8
#define NUMDOF_SHELL8    (6)
#endif
#ifndef MAXHYB_SHELL8
#define MAXHYB_SHELL8    (45)
#endif

/*----------------------------------------------------------------------*
 | shell9                                                               |
 |                                                                      |
 |  set MAXDOFPERNODE = 3+3*klay                                        |
 |  set MAXGAUSS      = lr*ls*lt*numlay (numlay=total number of layers) |
 |                                                                      |
 *----------------------------------------------------------------------*/
#ifndef MAXNOD_SHELL9
#define MAXNOD_SHELL9    (9)
#endif
#ifndef MAXLAY_SHELL9
#define MAXLAY_SHELL9    (30)                   /* max. total number
                                                 * of layers */
#endif
#ifndef MAXKLAY_SHELL9
#define MAXKLAY_SHELL9   (10)                   /* max. nr of kin lay
                                                 * */
#endif
#ifndef NUMDOF_SHELL9
#define NUMDOF_SHELL9    (3+3*MAXKLAY_SHELL9)   /* numdf = 3 + 3*klay
                                                 * */
#endif
#ifndef MAXHYB_SHELL9
#define MAXHYB_SHELL9    (45)
#endif
#ifndef A3FAC_SHELL9
#define A3FAC_SHELL9     (1.0) /* makes it possible to change the norm
                                * of a3L */
#endif
                               /* A3FAC_SHELL9 = 0.5 : |a3L| = 0.5 * hL -> as in shell8 */
                               /* A3FAC_SHELL9 = 1.0 : |a3L| = 1.0 * hL -> like in Dis. Braun */
#ifndef MAXNODESTRESS_SHELL9
#define MAXNODESTRESS_SHELL9 (72) /* Numnod*lt*numlay with Numnod=4
                                   * for Quad4   */
#endif
                                  /*                       Numnod=9 for Quad8/9 */

/*----------------------------------------------------------------------*
 | wall1                                                                |
 *----------------------------------------------------------------------*/
#ifndef MAXNOD_WALL1
#define MAXNOD_WALL1     (9)
#endif

/*----------------------------------------------------------------------*
 | therm2                                                   bborn 03/06 |
 *----------------------------------------------------------------------*/
#ifndef NDIM_THERM2
#define NDIM_THERM2      (2)    /* planar problem also genprob.ndim */
#endif
#ifndef MAXNOD_THERM2
#define MAXNOD_THERM2    (9)    /* maximum of element nodes :
                                 * max quad9 */
#endif
#ifndef NUMDOF_THERM2
#define NUMDOF_THERM2    (1)    /* number of thermal DOFs at each node :
                                 * temperature */
#endif
#ifndef NUMPLST_THERM2
#define NUMPLST_THERM2   (2)    /* number of plane states :
                                 * plane temperature gradient,
                                 * plane heat flux */
#endif
#ifndef NUMTMGR_THERM2
#define NUMTMGR_THERM2   (2)    /* number of temperature gradients */
#endif
#ifndef NUMHFLX_THERM2
#define NUMHFLX_THERM2   (3)    /* number of heat fluxes : q_x, q_y, q_z */
#endif
#ifndef GLINTC_THERM2
#define GLINTC_THERM2    (6)    /* line domain Gauss integration cases */
#endif
#ifndef GLMAXP_THERM2
#define GLMAXP_THERM2    (6)    /* line domain max. number of
                                 * Gauss points */
#endif
#ifndef GTINTC_THERM2
#define GTINTC_THERM2    (11)   /* triangle domain Gauss
                                 * integration cases */
#endif
#ifndef GTMAXP_THERM2
#define GTMAXP_THERM2    (13)   /* max. number of Gauss points */
#endif

/*----------------------------------------------------------------------*
 | hex20 elements                                                      |
 *----------------------------------------------------------------------*/
#ifndef NODESHIFT_HEX20
#define NODESHIFT_HEX20  (10000) /* must be larger than genprob.nnode */
#endif

#ifndef ELESHIFT_HEX20
#define ELESHIFT_HEX20  (10000) /* must be larger than genprob.nele */
#endif


/* --------------------------------------------------------------------- *
 * subdivide
 * --------------------------------------------------------------------- */
#ifndef MAX_DIVIDE
#define MAX_DIVIDE  (10)
#endif


#ifndef MAX_WARNINGS
#define MAX_WARNINGS (20)
#endif
