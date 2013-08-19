/*!----------------------------------------------------------------------
\file meshtying_defines.H

\brief A set of preprocessor defines for meshtying methods

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
#ifndef MESHTYING_DEFINES_H
#define MESHTYING_DEFINES_H

/************************************************************************/
/* Meshtying algorithm options                                          */
/************************************************************************/

#define MESHTYINGTWOCON      /* flag for applying 2 static condensations */
#define MESHTYINGUCONSTR     /* flag for computing constraints b(u,\mu) */
#define DELTALM 		     /* flag for incremental lagrange multipliers */

/************************************************************************/
/* Debugging options                                                    */
/************************************************************************/

// GENERAL DEBUGGING OPTIONS
#define MESHTYINGASOUTPUT    /* flag for detailed active set output */

#endif  // #ifndef MESHTYING_DEFINES_H
