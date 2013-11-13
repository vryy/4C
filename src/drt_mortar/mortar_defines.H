/*!----------------------------------------------------------------------
\file mortar_defines.H

\brief A set of preprocessor defines for mortar methods

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef MORTAR_DEFINES_H
#define MORTAR_DEFINES_H


/************************************************************************/
/* Mortar algorithm options                                             */
/************************************************************************/

//#define MORTARTRAFO           /* flag for basis trafo for quadratic FE */

/************************************************************************/
/* Mortar algorithm parameters                                          */
/************************************************************************/

// MORTAR PROJECTION (2D/3D)
#define MORTARMAXITER 10       /* max. no. iterations for local Newton */
#define MORTARCONVTOL 1.0e-12  /* convergence tolerance for local Newton */

// MORTAR PROJECTION AND INTEGRATION (2D)
#define MORTARPROJTOL 0.05     /* projection tolerance for overlap */
#define MORTARPROJLIM 1.0e-8   /* exact projection limit (no tolerance!) */

// MORTAR PROJECTION AND INTEGRATION (3D)
#define MORTARCLIPTOL 1.0e-8   /* tolerance for polygon clipping */
#define MORTARINTLIM 1.0e-12   /* min(area-%) cell/slave for integration */

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
