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

/*----------------------------------------------------------------------*
 *                                                                      *
 * General definitions of max. sizes                                    *
 *                                                                      *
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | maximum number of nodes to an element                                |
 *----------------------------------------------------------------------*/
#define MAXNOD           (8)

/*----------------------------------------------------------------------*
 | maximum number of elements to a node                                 |
 *----------------------------------------------------------------------*/
#define MAXELE           (8)

/*----------------------------------------------------------------------*
 | maximum number of dofs to a node -> can be more than 6 for shell9    |
 *----------------------------------------------------------------------*/
#define MAXDOFPERNODE    (5)

/*----------------------------------------------------------------------*
 | maximum number of gaussian points in an element                      |
 *----------------------------------------------------------------------*/
#define MAXGAUSS         (8)  

/*----------------------------------------------------------------------*
 | maximum number of fields                                             |
 *----------------------------------------------------------------------*/
#define MAXFIELD         (3)

/*----------------------------------------------------------------------*
 | maximum numberof fields                                              |
 *----------------------------------------------------------------------*/
#define MAXTIMECURVE     (5)

/*----------------------------------------------------------------------*
 | maximum numberof restart records per element in dynamic and nonlinear|
 | static restarting from and to pss file                               |
 *----------------------------------------------------------------------*/
#define MAXRECORDPERELE  (6)



/*----------------------------------------------------------------------*
 *                                                                      *
 * element specific definitions of max. sizes                           *
 *                                                                      *
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | axishell                                                             |
 *----------------------------------------------------------------------*/
#define MAXNOD_AXISHELL  (2)

/*----------------------------------------------------------------------*
 | beam3                                                                |
 *----------------------------------------------------------------------*/
#define MAXNOD_BEAM3     (3)

/*----------------------------------------------------------------------*
 | brick8                                                               |
 *----------------------------------------------------------------------*/
#define MAXNOD_BRICK1    (20)
#define NUMDOF_BRICK1    (3)

/*----------------------------------------------------------------------*
 | fluid                                                                |
 *----------------------------------------------------------------------*/
#define MAXQINTC         (6)  /* max. number of implemented integration cases for QUADS */
#define MAXQINTP         (6)  /* max. number of integration parameters  for QUADS */ 
#define MAXTINTC         (11) /* max. number of implemented integration cases for TRIS */
#define MAXTINTP         (13) /* max. number of integration parameters  for TRIS */
#define FLUID_NUM_LD     (4)  /* number of lift and drag conditions */

/*----------------------------------------------------------------------*
 | fluid2                                                               |
 *----------------------------------------------------------------------*/
#define NUM_F2_VELDOF    (2)  /* number of velocity dofs per node */
#define NUMDOF_FLUID2    (3)  /* number of dofs per node */
#define MAXNOD_F2        (9)  /* max. number of nodes per element */

/*----------------------------------------------------------------------*
 | fluid3                                                               |
 *----------------------------------------------------------------------*/
#define NUM_F3_VELDOF    (3)  /* number of velocity dofs per node */
#define MAXNOD_F3        (27) /* max. number of nodes per element */

/*----------------------------------------------------------------------*
 | shell8                                                               |
 *----------------------------------------------------------------------*/
#define MAXNOD_SHELL8    (9)
#define NUMDOF_SHELL8    (6)
#define MAXHYB_SHELL8    (45)

/*----------------------------------------------------------------------*
 | shell9                                                               |
 |                                                                      |
 |  set MAXDOFPERNODE = 3+3*klay                                        |
 |  set MAXGAUSS      = lr*ls*lt*numlay (numlay=total number of layers) |
 |                                                                      |
 *----------------------------------------------------------------------*/
#define MAXNOD_SHELL9    (9)
#define MAXLAY_SHELL9    (30)                   /* max. total number of layers */
#define MAXKLAY_SHELL9   (10)                   /* max. nr of kin lay */
#define NUMDOF_SHELL9    (3+3*MAXKLAY_SHELL9)   /* numdf = 3 + 3*klay */
#define MAXHYB_SHELL9    (45)
#define A3FAC_SHELL9     (1.0) /* makes it possible to change the norm of a3L */
                               /* A3FAC_SHELL9 = 0.5 : |a3L| = 0.5 * hL -> as in shell8 */
                               /* A3FAC_SHELL9 = 1.0 : |a3L| = 1.0 * hL -> like in Dis. Braun */
#define MAXNODESTRESS_SHELL9 (72) /* Numnod*lt*numlay with Numnod=4 for Quad4   */
                                  /*                       Numnod=9 for Quad8/9 */

/*----------------------------------------------------------------------*
 | wall1                                                                |
 *----------------------------------------------------------------------*/
#define MAXNOD_WALL1     (9)

/*----------------------------------------------------------------------*
 | ls2                                                                  |
 *----------------------------------------------------------------------*/
#define MAXNOD_LS2 (4)
