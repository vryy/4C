/*!-----------------------------------------------------------------------------------------------*
\file scatra_reinit_defines.H

  \brief preprocessor define flags for reinitialization via solving an advection equation
         only used for level set problems

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef SCATRA_REINIT_DEFINES_H
#define SCATRA_REINIT_DEFINES_H



#define NONLINEAR                // nonlinear incremental formulation, typically switched on
#define PHIN_INSTEAD_OF_PHI_0    // for sign for transport use the solution of the last pseudo time step

// evaluate l1-norm for steady state check within a transition region
#define L1_NORM_TRANSITION_REGION // transition region is set to 3*h


//---------------------------
// Characteristic Galerkin method for reinitialization
//---------------------------
//#define FIXPOINTLIKE           // fixpoint like iteration
//#define HIGHER_ORDER_NEW       // additional terms for higher order elements, see Lin2005
//#define MIDPOINT               // penalty midpoint for intersection penalty method


//----------------------------
// LinearAdvection with SUPG stabilization
#define REINIT_LINEAR_ADVECTION_PHIGRADIENT    // original phi-gradients, typically switched on
//#define REINIT_LINEAR_ADVECTION_RECONSTRUCTED_NORMALS // reconstructed phi-gradients

//#define DONT_EVALUATE_SMALL_GRADIENTS // do not assemble elements with too small gradients (only mass matrix)
//----------------------------



#endif // #ifdef SCATRA_REINIT_DEFINES_H
#endif // #ifdef CCADISCRET
