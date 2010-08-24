/*!-----------------------------------------------------------------------------------------------*
\file combust_defines.H

  \brief preprocessor define flags for combustion problems

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef COMBUST_DEFINES_H
#define COMBUST_DEFINES_H

//-----------------------------
// combustion algorithm options
//-----------------------------
//#define COMBUST_DECOUPLEDXFEM   // turn on/off boundary integrals for premixed combustion system matrix
//#define COMBUST_NITSCHE
#define COMBUST_STRESS_BASED
//#define TWOPHASEFLOW_NITSCHE
//#define COMBUST_PROJCELLEDGES // turn on/off old method to find facing flame front patches
//------------------------------------------
// combustion domain integral terms standard
//------------------------------------------
#define COMBUST_DOM_INERTIA    // turn on/off inertia term
#define COMBUST_DOM_CONVECTION // turn on/off concevtion term
#define COMBUST_DOM_DIFFUSION  // turn on/off diffusion term
#define COMBUST_DOM_PRESSURE   // turn on/off pressure term
#define COMBUST_DOM_CONTINUITY // turn on/off continuity term
#define COMBUST_DOM_RHS        // turn on/off RHS term

#define COMBUST_SUPG_INERTIA    // turn on/off SUPG inertia term
#define COMBUST_SUPG_CONVECTION // turn on/off SUPG concevtion term
#define COMBUST_SUPG_DIFFUSION  // turn on/off SUPG diffusion term
#define COMBUST_SUPG_PRESSURE   // turn on/off SUPG pressure term
#define COMBUST_SUPG_BODYFORCE
#define COMBUST_SUPG_RHS

#define COMBUST_PSPG_INERTIA    // turn on/off PSPG inertia term
#define COMBUST_PSPG_CONVECTION // turn on/off PSPG concevtion term
#define COMBUST_PSPG_DIFFUSION  // turn on/off PSPG diffusion term
#define COMBUST_PSPG_PRESSURE   // turn on/off PSPG pressure term
#define COMBUST_PSPG_RHS        // turn on/off PSPG RHS term

#define COMBUST_CSTAB // turn on/off CSTAB inertia term

//-------------------------------------------
// combustion Nitsche boundary integral terms
//-------------------------------------------
#define COMBUST_NITSCHE_VISCCONSISTMEAN   // turn on/off
#define COMBUST_NITSCHE_PRESCONSISTMEAN   // turn on/off
//#define COMBUST_NITSCHE_VISCCONSISTJUMP   // turn on/off
#define COMBUST_NITSCHE_VISCADJCONSIST    // turn on/off
//#define COMBUST_NITSCHE_PRESCONSISTJUMP   // turn on/off
#define COMBUST_NITSCHE_PRESADJCONSIST    // turn on/off pressure adjoint consistency term
#define COMBUST_NITSCHE_FLUXJUMP          // turn on/off pressure/flux jump RHS term
#define COMBUST_NITSCHE_NITSCHEVEL        // turn on/off Nitsche velocity terms
//#define COMBUST_NITSCHE_NITSCHEPRES       // turn on/off Nitsche pressure terms

//----------------------------------------------------------
// combustion domain integral terms stress-based formulation
//----------------------------------------------------------
#define COMBUST_STRESS_BASED_TAU_COUPLING
#define COMBUST_STRESS_BASED_PRESSURE_COUPLING
#define COMBUST_STRESS_BASED_INCOMPR
//--------------------------------------------------------------------
// combustion stress-based double one-sided formulation boundary terms
//--------------------------------------------------------------------
//#define COMBUST_STRESS_BASED_DOUBLE_ONESIDED
//------------------------------------------------
// combustion stress-based boundary integral terms
//------------------------------------------------
//#define COMBUST_STRESS_BASED_CONSISTVISCJUMP    // turn on/off
//#define COMBUST_STRESS_BASED_CONSISTPRESJUMP    // turn on/off
#define COMBUST_STRESS_BASED_FLUXJUMP           // turn on/off
#define COMBUST_STRESS_BASED_CONSISTVISCMEAN    // turn on/off
#define COMBUST_STRESS_BASED_CONSISTPRESMEAN    // turn on/off
#define COMBUST_STRESS_BASED_ADJCONSISTVISCJUMP // turn on/off
#define COMBUST_STRESS_BASED_ADJCONSISTPRESJUMP // turn on/off
//#define COMBUST_STRESS_BASED_NEUMANN_INFLOW     // turn on/off
//#define COMBUST_STRESS_BASED_ADVECTIVEFLUX      // turn on/off

#define COMBUST_STRESS_BASED_DOUBLE_ONESIDED_CONSIST_TAU_PLUS
#define COMBUST_STRESS_BASED_DOUBLE_ONESIDED_CONSIST_TAU_MINUS
#define COMBUST_STRESS_BASED_DOUBLE_ONESIDED_CONSIST_PRES_PLUS
#define COMBUST_STRESS_BASED_DOUBLE_ONESIDED_CONSIST_PRES_MINUS
#define COMBUST_STRESS_BASED_DOUBLE_ONESIDED_TAU_DIRICH_PLUS
#define COMBUST_STRESS_BASED_DOUBLE_ONESIDED_TAU_DIRICH_MINUS
#define COMBUST_STRESS_BASED_DOUBLE_ONESIDED_PRES_DIRICH_PLUS
#define COMBUST_STRESS_BASED_DOUBLE_ONESIDED_PRES_DIRICH_MINUS

//-----------------------------
// combustion test case options
//-----------------------------
//#define COMBUST_TESTCOUETTEFLOW          // modify Dirichlet boundary conditions manually for Couette flow test case
//#define COMBUST_TESTCOUETTEFLOWDECOUPLED // velocity jump for combust_stat_inclined_couette_decoupled.dat
//#define COMBUST_XVELFREE // free Dirichlet boundary conditions x-component manually
//#define COMBUST_YVELFREE // free Dirichlet boundary conditions y-component manually
//#define COMBUST_PRESFREE

//----------------------
// combustion tolerances
//----------------------
//#define COMBUSTCONVTOL 1.0e-12 // just an example
#define COMBUST_GMSH_NORMALFIELD
//#define COMBUST_GMSH_SMOOTHED_NORMALFIELD


//---------------------------------
// TwoPhaseFlow jump-conditions
//---------------------------------
//#define TPF_COMPARE_to_ANALYT_Sol_Nitsche

//---------------------------velocity jump ------------------------
//#define TPF_NORMAL_VEL_JUMP_CONST 	                   // else zero jump in velocity
//#define TPF_NORMAL_VEL_JUMP_PREMIXED_COMBUSTION          // else zero jump in velocity

//---------------------------normal tension jump ------------------
//#define TPF_NORMAL_TENSION_JUMP_CONST                    // else zero jump in normal tension
//#define TPF_NORMAL_TENSION_JUMP_PREMIXED_COMBUSTION      // else zero jump in normal tension
//#define TPF_SURFACETENSION // else zero jump in normal tension


//---------------------------surface tension options---------------
// => -----choose surface tension approximation type -----------------
//#define TPF_SURFACETENSION_NOT_CONNECTED_INTERFACE

//#define TPF_LAPLACE_BELTRAMI
//#define TPF_LAPLACE_BELTRAMI_SMOOTHED 	//Vers.2 of Laplace-Beltrami surface tension for improved calculation of curvature
//#define TPF_SURFACETENSION_CONST_CURVATURE
//#define TPF_SURFACETENSION_DIVERGENCE

//#define TPF_PHI_RECONSTRUCT



//---------------------------------
// TwoPhaseFlow domain integral terms
//---------------------------------

#define TPF_DOM_INERTIA    // turn on/off inertia term
#define TPF_DOM_CONVECTION // turn on/off concevtion term
#define TPF_DOM_DIFFUSION  // turn on/off diffusion term
#define TPF_DOM_PRESSURE   // turn on/off pressure term
#define TPF_DOM_CONTINUITY // turn on/off continuity term
#define TPF_DOM_RHS        // turn on/off RHS term

#define TPF_SUPG_INERTIA    // turn on/off SUPG inertia term
#define TPF_SUPG_CONVECTION // turn on/off SUPG concevtion term
#define TPF_SUPG_DIFFUSION  // turn on/off SUPG diffusion term
#define TPF_SUPG_PRESSURE   // turn on/off SUPG pressure term
#define TPF_SUPG_BODYFORCE
#define TPF_SUPG_RHS
//
#define TPF_PSPG_INERTIA    // turn on/off PSPG inertia term
#define TPF_PSPG_CONVECTION // turn on/off PSPG concevtion term
#define TPF_PSPG_DIFFUSION  // turn on/off PSPG diffusion term
#define TPF_PSPG_PRESSURE   // turn on/off PSPG pressure term
#define TPF_PSPG_RHS        // turn on/off PSPG RHS term

#define TPF_CSTAB // turn on/off CSTAB inertia term

//-----------------------------------
// TwoPhaseFlow boundary integral terms
//-----------------------------------
#define TPF_VISCCONSISTMEAN    // turn on/off
//#define TPF_VISCCONSISTJUMP    // turn on/off
#define TPF_VISCADJCONSIST     // turn on/off
#define TPF_PRESCONSISTMEAN    // turn on/off
//#define TPF_PRESCONSISTJUMP    // turn on/off
#define TPF_PRESADJCONSIST     // turn on/off pressure adjoint consistency term
//#define TPF_PRESJUMPRHS        // turn on/off pressure/flux jump RHS term
#define TPF_NITSCHEVEL         // turn on/off Nitsche velocity terms
//#define TPF_NITSCHEPRES      // turn on/off Nitsche pressure terms

// #define TPF_NITSCHE_LPS			// turn on/off Nitsche pressure_LPS_stabilization at the interface
//#define TPF_NITSCHE_LPS_DOMAIN
//-----------------------------
// TwoPhaseFlow test case options
//-----------------------------
//#define COMBUST_TESTCOUETTEFLOW// modify Dirichlet boundary conditions manually for Couette flow test case
//#define COMBUST_XVELFREE      // free Dirichlet boundary conditions manually for Couette flow test case

//----------------------
// TwoPhaseFlow tolerances
//----------------------
//#define COMBUSTCONVTOL 1.0e-12 // just an example

#define TRAFO_OWN






#endif // #ifndef COMBUST_DEFINES_H
#endif // #ifdef CCADISCRET
