/*-----------------------------------------------------------------------*/
/*! \file
\brief A set of preprocessor defines for mortar methods

\level 1

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------*/
#ifndef MORTAR_DEFINES_H
#define MORTAR_DEFINES_H

/************************************************************************/
/* Mortar algorithm parameters                                          */
/************************************************************************/

// MORTAR INTEGRATION
#define MORTARINTTOL 0.0 /* tolerance for assembling gp-values*/

// MORTAR PROJECTION (2D/3D)
#define MORTARMAXITER 10      /* max. no. iterations for local Newton */
#define MORTARCONVTOL 1.0e-12 /* convergence tolerance for local Newton */

// MORTAR PROJECTION AND INTEGRATION (2D)
#define MORTARPROJTOL 0.05   /* projection tolerance for overlap */
#define MORTARPROJLIM 1.0e-8 /* exact projection limit (no tolerance!) */

// MORTAR PROJECTION AND INTEGRATION (3D)
#define MORTARCLIPTOL 1.0e-8 /* tolerance for polygon clipping */
#define MORTARINTLIM 1.0e-12 /* min(area-%) cell/slave for integration */

/************************************************************************/
/* Mortar debugging options                                             */
/************************************************************************/

// GMSH DEBUGGING OPTIONS
//#define MORTARGMSH1          /* gmsh output of interface in each time step */
//#define MORTARGMSH2          /* gmsh output of interface in each iteration */
//#define MORTARGMSH3          /* gmsh output of interface at t = 0.0 */
//#define MORTARGMSHCELLS      /* gmsh output of intcells in each iteration */
//#define MORTARGMSHTN         /* gmsh output of all treenodes */
//#define MORTARGMSHCTN        /* gmsh output of coupling treenodes */

#endif  // #ifndef MORTAR_DEFINES_H
