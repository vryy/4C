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
#ifndef BEAM3CONTACT_DEFINES_H
#define BEAM3CONTACT_DEFINES_H

#include "Sacado.hpp"

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
//default: 10
#define BEAMCONTACTTOL 1.0e-10    /* convergence tolerance for local Newton */
//default: 1.0e-12
#define XIETATOL 1.0e-7           /* tolerance for valid values of beam parameters Xi and Eta
//default: 1.0e-9                                     (e.g. Xi \in [-1 - XIETATOL; 1 + XIETATOL]) */
#define NORMTOL  1.0e-12           /* tolerance when the to contact points are to close and one
//default: 1.0e-15                                    beam has to be shifted artificially */
#define SHIFTVALUE 1.0e-7         /* value, one contact point will be shiftet, if contact points are to close */
//default: 1.0e-7
#define COLINEARTOL 1.0e-12       /* tolerance when two elements are almost colinear */
//default: 1.0e-12
#define DETERMINANTTOL 1.0e-12    /* tolerance how small the determinant appearing in the contact point linearization is allowed to be */
//default: 1.0e-12
#define NORMALTOL 0.00001         /* tolerance how small the scalar product of normal_old and normal is allowed to be */
//default: 0.001 / 0.00001
#define MANIPULATERADIUS 1.0        /* with this parameter the geometrical radius used for contact evaluation can be
//default: 1                         chosen different from the physical radius appearing in the bending stiffness */
#define MANIPULATERADIUSVIS 1.0     /* with this parameter the geometrical radius used for gmsh visualization can be
//default: 1                         chosen different from the physical radius appearing in the bending stiffness */
#define NEIGHBORTOL 3            /* parameter that allows to also evaluate the contact normal of neighbor elements
//default: 3                         needed for new gap function. It has proven useful that NEIGHBORTOL ca. 3,
                                     since the sum of |xi1| and |xi2| has to be larger than three for a contact pair
                                     consisting of element i and element i+2 of the same physical beam! */
#define PARALLELTOL 1.0e-10        //Tolerance necessary for start values eta1_0 and eta2_0 of closest point projection
//default: 1.0e-10
#define GAPTOL 1.0e-12            /* tolerance for the allowable difference between old and new gap definition
//default: 1.0e-12                   which should be identical in the converged state */

#define DAMPTOL 1.0e-8              //Tolerance for minimal difference of the damping regularization parameter

#define NEIGHBORNORMALTOL 5         //tolerance necessary in order to calcualte normal_old_ out of the neighbor element pair!

#define MAXDELTADFAC       0.8       //Maximal allowed nodal displacement per time step as multiple bounding box increment

#define MAXDELTAXIETA       2.0     //Maximal allowed parameter step xi_-xi_old per time step

#define PARALLEL_DEACTIVATION_VAL 0.9  //Value of scalar product x1_xi*x2_xi, at which the beams are considered as parallel

#define PARALLEL_ACTIVATION_VAL 0.6  //Value of above scalar product x1_xi*x2_xi, at which the beams are again considered as not parallel

//#define AUTOMATICDIFF           //Decide weather automatic differentiation via Sacado is used or not
//default: Off
#ifdef AUTOMATICDIFF
  typedef Sacado::Fad::DFad<double> TYPE;
  #define FADCHECKS //Decide weather FAD checks via Sacado are printed as output or not
#else
  typedef double TYPE;
#endif

#define CHECKBOUNDARYCONTACT      //Detection of beam ends in order to open the contact, when one beam reaches the end of the second beam
//default: true

//#define MAXFORCE           50      //prescribe maximum penalty force (cut off value in combination with secant stiffness)
//#define BASICSTIFFWEIGHT    1.1    //Additional weighting of the basis stiffness term

#define ITERMAX              15      //Newton iteration at which PTC or ALGORITHMICDAMP is applied
//#define ALGORITHMICDAMP      50      //Apply algorithmicdamp proportional to iterative gap change (tends to zero in converged configuration)
//#define ALGDAMPREGFAC1       0.001
//#define ALGDAMPREGFAC2      -0.03
//#define ALGDAMPBASICSTIFFFAC 1.0
//#define ALGCOMPLETESTIFF
//#define BEAMCONTACTPTC       10        //Apply PTC terms to contact elements
//#define BEAMCONTACTPTCROT    0.0

/************************************************************************/
/* Debugging options                                                    */
/************************************************************************/
#define N_CIRCUMFERENTIAL 10      /* number of gmsh visualization devisions of beam element in circumferential direction*/
#define N_Axial 10                /* number of gmsh visualization devisions of beam element in axial direction*/


#define OUTPUTEVERY 1          //Write output only every OUTPUTEVERY^th step
//#define GMSHTIMESTEPS           /* write gmsh output in each time step */
//#define GMSHNEWTONSTEPS           /* write gmsh and console output in each Newton step */
//#define REACTIONFORCES          /* output of reaction forces and moments */
//#define BEAMCONTACTFDCHECKS     /* perform finite difference checks */

#endif /* BEAM3CONTACT_DEFINES_H */
