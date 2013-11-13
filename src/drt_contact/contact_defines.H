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
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef CONTACT_DEFINES_H
#define CONTACT_DEFINES_H


/************************************************************************/
/* Contact algorithm options                                            */
/************************************************************************/

//#define CONTACTPSEUDO2D           /* flag for pseudo-3D in 2D situations (u_z=0)*/
//#define CONSISTENTSTICK		        /* flag for consistent stick branch of the complementarity function */
//#define CONSISTENTSLIP			      /* flag for consistent slip branch of the complementarity function */
//#define OBJECTVARSLIPINCREMENT    /* flag for object-variant slip increment created in integrator*/
#define SIMPLETANGENT               /* flag for simple tangent formulation*/
//#define TANGENTCHECK              /* Check tangents from solution step n-1 to n --> enforcing same orientation */
//#define CONTACTCONSTRAINTXYZ

/************************************************************************/
/* Debugging options                                                    */
/************************************************************************/

// GENERAL DEBUGGING OPTIONS
//#define CONTACTEXPORT             /* flag for export of inttime, active set... to textfiles. USE WITH CONTACTASOUTPUT AND CONTACTFORCEOUTPUT */
//#define CONTACTASOUTPUT           /* flag for detailed active set output */
//#define CONTACTFORCEOUTPUT        /* flag for output of contact forces */
//#define CONTACTTIME               /* flag for contact time measurement */
//#define CONTACTSTATUS             /* flag for contact segmentation measurement */
//#define CONTACTEIG                /* flag for eigenvalue+condition analysis */
//#define CONTACTFORCEREFCONFIG     /* flag for evaluation of forces in ref. config. */
//#define CONTACTFDNORMAL           /* flag for FD check of normal derivative */
//#define CONTACTFDMORTARD          /* flag for FD check of mortar matrix D */
//#define CONTACTFDMORTARM          /* flag for FD check of mortar matrix M */
//#define CONTACTFDGAP              /* flag for FD check of weighted gap g */
//#define CONTACTFDSLIPINCR         /* flag for FD check of weighted obj.-var. slip */

// LAGRANGE DEBUGGING OPTIONS
//#define CONTACTFDTANGLM           /* flag for FD check of (frictionless) tangential LM */
//#define CONTACTFDSTICK            /* flag for FD check of (frictional) stick condition */
//#define CONTACTFDSLIP             /* flag for FD check of (frictional) slip condition */

// PENALTY DEBUGGING OPTIONS
//#define CONTACTFDPENALTYTRAC      /* flag for FD check of penalty tractions */

// IMPLICIT WEAR DEBUGGING OPTIONS
//#define WEARIMPLICITFD            /* flag for FD-check w.r.t. disp. */
//#define WEARIMPLICITFDLM          /* flag for FD-check w.r.t. lm-mult. */

// DISCRETE WEAR DEBUGGING OPTIONS
//#define CONTACTFDT_D              /* flag for FD check of mortar matrix T lin for wear cond */
//#define CONTACTFDE_D              /* flag for FD check of mortar matrix E lin for wear cond */
//#define CONTACTFDMORTART          /* flag for FD check of mortar matrix T */
//#define CONTACTFDMORTARE          /* flag for FD check of mortar matrix E */

/************************************************************************/
/* Printing and patient specific options                                */
/************************************************************************/

// PRINTING OF SPECIFIC INFORMATION
//#define MASTERNODESINCONTACT           /* flag for printing master nodes in contact - CONTACTFORCEOUTPUT also has to be switched on*/
//#define PRINTSTRUCTDEFORMEDNODECOORDS  /* flag for printing all nodes and xspatial in current configuration - only works for 1 processor*/

// PRESTRESSING CUSTOMIZATIONS
//#define RESTARTAFTERPRESTRESSING       /* flag for allowing restart with contact switched on after prestressing - NOTE: no contact has yet occured */

#endif  // #ifndef CONTACT_DEFINES_H
