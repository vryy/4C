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
#define COMBUST_NITSCHE_VISCCONSISTMEAN // turn on/off
#define COMBUST_NITSCHE_PRESCONSISTMEAN // turn on/off
//#define COMBUST_NITSCHE_VISCCONSISTJUMP // turn on/off
#define COMBUST_NITSCHE_VISCADJCONSIST  // turn on/off
//#define COMBUST_NITSCHE_PRESCONSISTJUMP // turn on/off
#define COMBUST_NITSCHE_PRESADJCONSIST  // turn on/off pressure adjoint consistency term
#define COMBUST_NITSCHE_PRESJUMPRHS     // turn on/off pressure/flux jump RHS term
#define COMBUST_NITSCHE_NITSCHEVEL      // turn on/off Nitsche velocity terms
//#define COMBUST_NITSCHE_NITSCHEPRES     // turn on/off Nitsche pressure terms

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
#define COMBUST_STRESS_BASED_CONSISTVISCMEAN    // turn on/off
#define COMBUST_STRESS_BASED_CONSISTPRESMEAN    // turn on/off
#define COMBUST_STRESS_BASED_ADJCONSISTVISCJUMP // turn on/off
#define COMBUST_STRESS_BASED_ADJCONSISTPRESJUMP // turn on/off

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
//#define COMBUST_TESTCOUETTEFLOW // modify Dirichlet boundary conditions manually for Couette flow test case
//#define COMBUST_XVELFREE // free Dirichlet boundary conditions x-component manually
//#define COMBUST_YVELFREE // free Dirichlet boundary conditions y-component manually

//----------------------
// combustion tolerances
//----------------------
//#define COMBUSTCONVTOL 1.0e-12 // just an example
#define COMBUST_GMSH_NORMALFIELD

#endif // #ifndef COMBUST_DEFINES_H
#endif // #ifdef CCADISCRET
