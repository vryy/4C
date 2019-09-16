/*----------------------------------------------------------------------------*/
/*! \file

\brief A set of preprocessor defines for beam contact methods

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/
#ifndef BEAM3CONTACT_DEFINES_H
#define BEAM3CONTACT_DEFINES_H

#include "Sacado.hpp"

/************************************************************************/
/* Beam contact algorithm options                                       */
/************************************************************************/
#define RELCONSTRTOL /* flag for relative evaluation of the constraint tolerance */
//#define MEASURETIME             /* flag for time measurement of octree setup and intermediate
// steps */ #define OCTREEDEBUG             /* flag for debug options in octree search routine*/

/************************************************************************/
/* Beam contact algorithm parameters                                    */
/************************************************************************/
#define BEAMCONTACTMAXITER 50 /* max number of steps for local Newton */
// default: 50
#define BEAMCONTACTTOL 1.0e-10 /* convergence tolerance for local Newton */
// default: 1.0e-10
#define RELBEAMCONTACTTOL 1.0e-6 /* apply relative convergence tolerance for local Newton */
// default: 1.0e-8
#define CRITICALRELNORM                                                                      \
  1.0e-6 /*we get an error message, if the relative norm of the local Newton is decreased by \
            CRITICALRELNORM, but no convergence is reached nevertheless*/
// default: 1.0e-8
#define XIETAITERATIVEDISPTOL                                                                     \
  1.0e-10 /* second convergence tolerance for local Newton: norm of iterative displacements in xi \
             and eta */
// default: 1.0e-10
#define XIETARESOLUTIONFAC                                                                   \
  10 /*this quantity defines the smallest possible difference in parameter space between two \
        different closest points as multiple of XIETAITERATIVEDISPTOL*/
// default: 10
#define XIETASAFETYFAC /*safety factor for closest point projection*/
// default: 10
#define XIETATOL                                                                                  \
  1.0e-7 /* tolerance for valid values of beam parameters Xi and Eta                              \
//default: 1.0e-9                                     (e.g. Xi \in [-1 - XIETATOL; 1 + XIETATOL]) \
*/
#define NORMTOL                                                        \
  1.0e-12 /* tolerance when the to contact points are to close and one \
//default: 1.0e-15                                    beam has to be shifted artificially */
#define SHIFTVALUE \
  1.0e-7 /* value, one contact point will be shiftet, if contact points are to close */
// default: 1.0e-7
#define COLINEARTOL 1.0e-12 /* tolerance when two elements are almost colinear */
// default: 1.0e-12
#define DETERMINANTTOL                                                                           \
  1.0e-12 /* tolerance how small the determinant appearing in the contact point linearization is \
             allowed to be */
// default: 1.0e-12
#define NORMALTOL \
  0.00001 /* tolerance how small the scalar product of normal_old and normal is allowed to be */
// default: 0.001 / 0.00001

#define MANIPULATERADIUS                                                                        \
  1.0 /* with this parameter the geometrical radius used for contact evaluation can be          \
//default: 1                         chosen different from the physical radius appearing in the \
bending stiffness */
#define MANIPULATERADIUSVIS                                                                     \
  1.0 /* with this parameter the geometrical radius used for gmsh visualization can be          \
//default: 1                         chosen different from the physical radius appearing in the \
bending stiffness */
#define NEIGHBORTOL                                                                                \
  2 /* parameter that allows to also evaluate the contact normal of neighbor elements              \
//default: 3                         needed for new gap function. It has proven useful that        \
NEIGHBORTOL ca. 3, since the sum of |xi1| and |xi2| has to be larger than three for a contact pair \
        consisting of element i and element i+2 of the same physical beam! */
#define PARALLELTOL \
  1.0e-10  // Tolerance necessary for start values eta1_0 and eta2_0 of closest point projection
// default: 1.0e-10
#define GAPTOL                                                                         \
  1.0e-12 /* tolerance for the allowable difference between old and new gap definition \
//default: 1.0e-12                   which should be identical in the converged state */

#define DAMPTOL 1.0e-8  // Tolerance for minimal difference of the damping regularization parameter

#define NEIGHBORNORMALTOL \
  5  // tolerance necessary in order to calcualte normal_old_ out of the neighbor element pair!

#define MAXDELTADFAC \
  0.8  // Maximal allowed nodal displacement per time step as multiple bounding box increment

#define MAXDELTAXIETA 2.0  // Maximal allowed parameter step xi_-xi_old per time step

#define PARALLEL_DEACTIVATION_VAL \
  0.9  // Value of scalar product x1_xi*x2_xi, at which the beams are considered as parallel

#define PARALLEL_ACTIVATION_VAL \
  0.6  // Value of above scalar product x1_xi*x2_xi, at which the beams are again considered as not
       // parallel

//#define AUTOMATICDIFF           //Decide whether automatic differentiation via Sacado is used or
// not default: Off
#ifdef AUTOMATICDIFF
// TODO
dserror(
    "check functionality of FAD based linearization in case of beam3k and beam3r with Hermite cl "
    "intpol");
typedef Sacado::Fad::DFad<double> TYPE;
//#define FADCHECKS //Decide whether FAD checks via Sacado are printed as output or not
#else
typedef double TYPE;
#endif

#define CHECKBOUNDARYCONTACT  // Detection of beam ends in order to open the contact, when one beam
                              // reaches the end of the second beam
// default: true

//#define MAXFORCE           50      //prescribe maximum penalty force (cut off value in combination
// with secant stiffness) #define BASICSTIFFWEIGHT    1.1    //Additional weighting of the basis
// stiffness term

#define ITERMAX 15  // Newton iteration at which PTC or ALGORITHMICDAMP is applied
//#define ALGORITHMICDAMP      50      //Apply algorithmicdamp proportional to iterative gap change
//(tends to zero in converged configuration) #define ALGDAMPREGFAC1       0.001 #define
// ALGDAMPREGFAC2      -0.03 #define ALGDAMPBASICSTIFFFAC 1.0 #define ALGCOMPLETESTIFF #define
// BEAMCONTACTPTC       10        //Apply PTC terms to contact elements #define BEAMCONTACTPTCROT
// 0.0


#define MAXNUMSEG 256  // default: 256         //maximal number of segments
#define ANGLETOL \
  0.1 / 180 * M_PI  // within the search algorithm, for angles<ANGLETOL two segments are assumed to
                    // be parallel
#define UNIQUECPPANGLE \
  1.0 / 180.0 *        \
      M_PI  // angle above which it is assumed that a unique cpp can be found by the local Newton

#define BEAMCONTACTGAUSSRULE intrule_line_5point  // Gauss rule applied in each segment
// default: 5 GP

//#define CONSISTENTTRANSITION 2                    //apply moment based transition: 1 -> variant 1,
// 2 -> variant 2 #define ENDPOINTSEGMENTATION                       //build integration segments at
// the physical beam endpoints
#define RELSEGMENTTOL \
  0.0000000001  // default:0.0000000001              //minimal (relative) size of a cut integration
                // interval
//#define CONVERGENCECONTROLE                      //check, if all local Newton loops converge

//#define CHANGEENDPOINTPROJECTION                 //If CHANGEENDPOINTPROJECTION is defined, the
// endpoint projection necessary for ENDPOINTSEGMENTATION is done as an orthogonal projection and
// not
// as an closest point projection

// TODO
//#define ONLYLEFTENDPOINTCONTACT                  //hack necessary for example "dynamic failure of
// rope" -> endpoint contact is only considered at one end of the ropes

#define CPP_APPROX  // Use approximation based on point-to-line projection in case of non-converent
                    // cpp

#define MAXCROSSSECTIONTOCURVATURE \
  0.005  // maximal allowed ratio of cross-section radius to curvature radius

#define INITSEG1 1  // minimal number of segments of element1 (necessary for debugging)
#define INITSEG2 1  // minimal number of segments of element2 (necessary for debugging)
//#define NOSEGMENTATION                             //switch segmentation of -> one element=one
// segment

// TODO: Change back to default
#define MAXPENETRATIONSAFETYFAC \
  0.8  // 0.8         //maximal allowed penetration is MAXPENETRATIONSAFETYFAC*(R1_+R2_)
// default: 0.4

#if defined(ENDPOINTPENALTY) and defined(CONSISTENTTRANSITION)
dserror("So far, CONSISTENTTRANSITION does not work in combination with ENDPOINTPENALTY!");
#endif

#if defined(ENDPOINTSEGMENTATION) and !defined(AUTOMATICDIFF)
dserror("ENDPOINTSEGMENTATION only works in combination with AUTOMATICDIFF!");
#endif

#if defined(CONSISTENTTRANSITION) and !defined(AUTOMATICDIFF)
dserror("CONSISTENTTRANSITION only works in combination with AUTOMATICDIFF!")
#endif

// For the CONSISTENTTRANSITION case, the variations of the contact point coordinates are necessary.
// For the small-angle contact formulation, this is only implemented in a way that considers
// variable eta while xi is fixed. However, in case of ENDPOINTSEGMENTATION both, xi and eta, are
// variable!
#if defined(CONSISTENTTRANSITION) and defined(ENDPOINTSEGMENTATION)
    dserror("CONSISTENTTRANSITION does not work in combination with ENDPOINTSEGMENTATION!")
#endif


/************************************************************************/
/* Debugging options                                                    */
/************************************************************************/
#define N_CIRCUMFERENTIAL \
  10 /* number of gmsh visualization devisions of beam element in circumferential direction*/
#define N_Axial 10 /* number of gmsh visualization devisions of beam element in axial direction*/

//#define DISTINGUISHCONTACTCOLOR
//#define CONTACTPAIRSPECIFICOUTPUT
//#define OUTPUTALLPROCS
//#define OUTPUTEVERY 0.000001        //Write gmsh output in time intervals of size OUTPUTEVERY
//#define OUTPUTEVERYNEWTONSTEP      //Write output also after each Newton step
//#define GMSHTIMESTEPS              /*write gmsh output in each time step */
//#define GMSHNEWTONSTEPS            /*write gmsh and console output in each Newton step */
//#define REACTIONFORCES             /* output of reaction forces and moments */
//#define BEAMCONTACTFDCHECKS        /* perform finite difference checks */

//#define PRINTGAPFILE //print file with gap/gap error
//#define PRINTNUMCONTACTSFILE//print number of active contact points over time
//#define PRINTGAPSOVERLENGTHFILE //print all GP gaps and forces over element length

//#define FDCHECK                //perform finite difference check of contact linearization



/************************************************************************/
/* Beam to solid contact options                                        */
/************************************************************************/

#define LINPENALTY
//#define QUADPENALTY
//#define ARBITPENALTY 2
#ifdef ARBITPENALTY
#define G0 3.0e-2  // G0 > 0
#endif
//#define GAP0 // GAP0 > 0                        // Decide if complete stiffness is used
#define GAUSSRULE \
  intrule_line_6point  // Define gauss rule for integrating residual and stiffness terms

//#define AUTOMATICDIFFBTS                        // Decide if automatic differentiation via Sacado
// is used default: Off
#ifdef AUTOMATICDIFFBTS
        typedef Sacado::Fad::DFad<double> TYPEBTS;
//  #define FADCHECKSTIFFNESS                     // Decide if FAD checks for contact stiffness are
//  printed as output
#ifndef FADCHECKSTIFFNESS
//    #define FADCHECKLINCONTACTINTERVALBORDER    // Decide if FAD checks for linearization of
//    contact interval borders are printed as output #define FADCHECKLINORTHOGONALITYCONDITION   //
//    Decide if FAD checks for linearization of orthogonality conditions are printed as output
//    #define FADCHECKLINGAUSSPOINT               // Decide if FAD checks for linearizations at
//    gauss points are printed as output
#endif
#else
typedef double TYPEBTS;
#endif

//#define FDCHECKSTIFFNESS                        // Decide if differentiation via finite difference
// for calculating contact stiffness is used default: Off

#ifndef AUTOMATICDIFFBTS
//#define BTSOLGMSH                           // gmsh output of solids
//#define GMSHDEBUG                             // Draw gap as 3D text at each Gauss point in Gmsh
//#define GMSHFORCE                             // Create extra Gmsh file with contact forces acting
// on solid element #define SAVEFORCE                             // Create text file with contact
// forces at all contact interval borders and Gauss points
#endif
// -----------------------------------------------------------------

#endif /* BEAM3CONTACT_DEFINES_H */
