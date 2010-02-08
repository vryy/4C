/*!----------------------------------------------------------------------
\file contact_defines.H

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
#ifndef CONTACT_DEFINES_H
#define CONTACT_DEFINES_H

/************************************************************************/
/* Contact algorithm options                                            */
/************************************************************************/

//#define CONTACTFRICTIONLESSFIRST  /* flag for frictionless contact in first contact step */
#define CONTACTCOMPHUEBER           /* flag for enhanced complementarity function mu(zn-cn.g) */
//#define CONTACTPSEUDO2D           /* flag for pseudo-3D in 2D situations (u_z=0)*/
//#define CONTACTPETROVGALERKINFRIC /* flag for PetrovGalerkin approach for LM with friction */

/************************************************************************/
/* Debugging options                                                    */
/************************************************************************/

// GENERAL DEBUGGING OPTIONS
//#define CONTACTTIME               /* flag for contact time measurement */
//#define CONTACTEIG                /* flag for eigenvalue+condition analysis */
//#define CONTACTFDNORMAL           /* flag for FD check of normal derivative */
//#define CONTACTFDMORTARD          /* flag for FD check of mortar matrix D */
//#define CONTACTFDMORTARM          /* flag for FD check of mortar matrix M */

// LAGRANGE DEBUGGING OPTIONS
//#define CONTACTFDGAP              /* flag for FD check of weighted gap g */
//#define CONTACTFDTANGLM           /* flag for FD check of tangential LM */

// PENALTY DEBUGGING OPTIONS
//#define CONTACTFDPENALTYKTEFF     /* flag for FD check of global assembled contact tangent stiffness matrix */
//#define CONTACTFDPENALTYKC1       /* flag for FD check of global assembled LinD and LinM part of contact tangent stiffness matrix */
//#define CONTACTFDPENALTYDERIVTRAC /* flag for FD check of nodal derivatives of lagrange multipliers */

// FRICTIONAL DEBUGGING OPTIONS
//#define CONTACTFDSTICK            /* flag for FD check of stick condition */
//#define CONTACTFDSLIP             /* flag for FD check of slip condition */

// GMSH DEBUGGING OPTIONS
//#define CONTACTGMSH1              /* flag for gmsh output in each time step */
//#define CONTACTGMSH2              /* flag for gmsh output in each iteration */
//#define CONTACTGMSH3              /* flag for gmsh output at t = 0.0 */

#endif  // #ifndef CONTACT_DEFINES_H
#endif  // #ifdef CCADISCRET
