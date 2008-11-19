/*!----------------------------------------------------------------------
\file contactdefines.H

\brief A set of preprocessor defines for contact methods

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
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONTACTDEFINES_H
#define CONTACTDEFINES_H

/************************************************************************/
/* Contact algorithm options                                            */
/************************************************************************/

//#define CONTACTWNORMAL         /* flag for weighted nodal normal calc. */
//#define CONTACTBOUNDMOD        /* flag for modification of boundary nodes */
//#define CONTACTONEMORTARLOOP   /* flag for computing D and M in one loop */

/************************************************************************/
/* Debugging options                                                    */
/************************************************************************/

//#define CONTACTFDNORMAL        /* flag for FD check of normal derivative */
//#define CONTACTFDMORTARD       /* flag for FD check of mortar matrix D */
//#define CONTACTFDMORTARM       /* flag for FD check of mortar matrix M */
//#define CONTACTFDGAP           /* flag for FD check of weighted gap g */
//#define CONTACTFDTANGLM        /* flag for FD check of tangential LM */
//#define CONTACTGMSH1           /* flag for gmsh output in each time step */
//#define CONTACTGMSH2           /* flag for gmsh output in each iteration */

/************************************************************************/
/* Contact algorithm parameters                                         */
/************************************************************************/

#define CONTACTPROJTOL 0.05      /* projection tol. for overlap decision */
#define CONTACTPROJLIM 1.0e-8    /* no problems with exact node projections */
#define CONTACTMAXITER 10        /* max. no. iterations for local Newton */
#define CONTACTCONVTOL 1.0e-12   /* convergence tolerance for local Newton */
#define CONTACTCRITDIST 0.5      /* critical distance for potential contact */
#define CONTACTCLIPTOL 1.0e-6    /* tolerance for polygon clipping in 3D */
#define CONTACTDEPTHOCTREE 5     /* max. depth for contact search octree */

namespace CONTACT {
  // static solution routine
  void contact_stru_static_drt();
}

#endif  // #ifndef CONTACTDEFINES_H
#endif  // #ifdef CCADISCRET
