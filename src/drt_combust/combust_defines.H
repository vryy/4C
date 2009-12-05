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
#define COMBUST_DECOUPLEDXFEM   // turn on/off boundary integrals for premixed combustion system matrix
//#define COMBUST_NITSCHE
//#define COMBUST_AXEL
//---------------------------------
// combustion domain integral terms
//---------------------------------
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

//-----------------------------------
// combustion boundary integral terms
//-----------------------------------
//#define COMBUST_VISCCONSIST    // turn on/off
//#define COMBUST_VISCADJCONSIST // turn on/off
//#define COMBUST_PRESCONSIST    // turn on/off
//#define COMBUST_PRESADJCONSIST // turn on/off pressure adjoint consistency term
//#define COMBUST_PRESJUMPRHS    // turn on/off pressure/flux jump RHS term
//#define COMBUST_NITSCHEVEL     // turn on/off Nitsche velocity terms
//#define COMBUST_NITSCHEPRES    // turn on/off Nitsche pressure terms

//-----------------------------
// combustion test case options
//-----------------------------
//#define COMBUST_TESTCOUETTEFLOW // modify Dirichlet boundary conditions manually for Couette flow test case
#define COMBUST_XVELFREE      // free Dirichlet boundary conditions manually for Couette flow test case

//----------------------
// combustion tolerances
//----------------------
//#define COMBUSTCONVTOL 1.0e-12 // just an example


#endif // #ifndef COMBUST_DEFINES_H
#endif // #ifdef CCADISCRET
