/*!----------------------------------------------------------------------
\file beam3contact_defines.H

\brief A set of preprocessor defines for contact methods

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef BEAM3CONTACT_DEFINES_H
#define BEAM3CONTACT_DEFINES_H

/************************************************************************/
/* Beam contact algorithm options                                       */
/************************************************************************/
#define RELCONSTRTOL              /* flag for relative evaluation of the constraint tolerance */
//#define MEASURETIME             /* flag for time measurement of octree setup and intermediate steps */
//#define OCTREEDEBUG             /* flag for debug options in octree search routine*/

/************************************************************************/
/* Beam contact algorithm parameters                                    */
/************************************************************************/
#define BEAMCONTACTMAXITER 50     /* max number of steps for local Newton */
#define BEAMCONTACTTOL 1.0e-12    /* convergence tolerance for local Newton */
#define NDIM 3                    /* number of spatial dimensions */
#define XIETATOL 1.0e-9           /* tolerance for valid vaulues of beam parameters Xi and Eta
                                     (e.g. Xi \in [-1 - XIETATOL; 1 + XIETATOL]) */
#define NORMTOL  1.0e-9           /* tolerance when the to contact points are to close and one
                                     beam has to be shifted artificially */
#define SHIFTVALUE 1.0e-7         /* value, one contact point will be shiftet, if contact points are to close */
#define COLINEARTOL 1.0e-12       /* tolerance when two elements are almost colinear */
#define NORMALTOL 0.1             /* tolerance how small the scalar product of normal_old and normal is allowed to be */
#define MANIPULATERADIUS 1        /* with this parameter the geometrical radius used for contact evaluation can be
                                     chosen different from the physical radius appearing in the bending stiffness */
#define NEIGHBORTOL 2             /* parameter that allows to also evaluate the contact normal of neighbor elements
                                     needed for new gap function. It has proven useful that NEIGHBORTOL \in [1 ; 3] */
#define GAPTOL 1.0e-12            /* tolerance for the allowable difference between old and new gap definition
                                     which should be identical in the converged state */

/************************************************************************/
/* Debugging options                                                    */
/************************************************************************/
//#define GMSHTIMESTEPS           /* write gmsh output in each time step */
//#define GMSHNEWTONSTEPS         /* write gmsh and console output in each Newton step */
//#define REACTIONFORCES          /* output of reaction forces and moments */
//#define BEAMCONTACTFDCHECKS     /* perform finite difference checks */

#endif /* BEAM3CONTACT_DEFINES_H */
#endif /* CCADISCRET */
