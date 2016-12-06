/*!-----------------------------------------------------------------------------------------------*
\file combust_defines.H

  \brief preprocessor define flags for combustion problems

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef COMBUST_DEFINES_H
#define COMBUST_DEFINES_H

//---------------------------
// general combustion options
//---------------------------
//#define COMBUST_2D // turn on/off modifications for 2D (pseudo 3D) problems
#define COMBUST_CUT
//-----------------------------
// combustion algorithm options
//-----------------------------
//#define COMBUST_DECOUPLEDXFEM   // turn on/off boundary integrals for premixed combustion system matrix
//#define COMBUST_PROJCELLEDGES // turn on/off old method to find facing flame front patches
//#define COMBUST_SXFEM
//#define COMBUST_SETJUMP
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

#define COMBUST_GRAD_DIV // turn on/off GRAD_DIV inertia term

//-------------------------------------------
// combustion Nitsche boundary integral terms
//-------------------------------------------
#define COMBUST_NITSCHE_VISCCONSISTMEAN   // turn on/off viscous consistency term
#define COMBUST_NITSCHE_PRESCONSISTMEAN   // turn on/off pressure consistency term
#define COMBUST_NITSCHE_VISCADJCONSIST    // turn on/off viscous adjoint consistency term
#define COMBUST_NITSCHE_PRESADJCONSIST    // turn on/off pressure adjoint consistency term
#define COMBUST_NITSCHE_FLUXJUMP          // turn on/off flux jump RHS term
#define COMBUST_NITSCHE_NITSCHEVEL        // turn on/off Nitsche velocity stabilization term
#define COMBUST_NITSCHE_CONVFLUX          // turn on/off convective inflow flux terms
#define COMBUST_NITSCHE_CONVSTAB          // turn on/off Nitsche convective stabilization terms
//#define COMBUST_NITSCHE_CONVPENALTY       // turn on/off Nitsche normal and tangential penalty/stabilization terms

//-----------------------------
// combustion test case options
//-----------------------------
//#define COMBUST_TESTCOUETTEFLOW          // modify Dirichlet boundary conditions manually for Couette flow test case
//#define COMBUST_TESTCOUETTEFLOWDECOUPLED // velocity jump for combust_stat_inclined_couette_decoupled.dat
//#define COLLAPSE_FLAME
//#define FLAME_VORTEX
//#define DL_INSTAB
//#define ORACLES
//#define COMBUST_TWO_FLAME_FRONTS
//#define COMBUST_BUNSENBURNER
//#define COMBUST_XVELFREE // free Dirichlet boundary conditions x-component manually
//#define COMBUST_YVELFREE // free Dirichlet boundary conditions y-component manually
//#define COMBUST_PRESFREE

//--------------------------
// combustion output options
//--------------------------
//#define GMSH_AVERAGE_JUMP
#define GMSH_WRITE_JUMPS
//#define COMBUST_GMSH_NORMALFIELD
//#define GMSH_REF_FIELDS


#endif // #ifndef COMBUST_DEFINES_H
