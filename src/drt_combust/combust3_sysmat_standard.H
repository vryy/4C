/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat_standard.H

\brief system matrix standard FEM terms
       for premixed combustion and two-phase flow problems / 3D combustion element

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef COMBUST3_SYSMAT_STANDARD_H
#define COMBUST3_SYSMAT_STANDARD_H

#include <Teuchos_TimeMonitor.hpp>

#include "combust3_local_assembler.H"
#include "combust3_utils.H"
#include "combust_defines.H"
#include "../drt_xfem/xfem_element_utils.H"


using namespace XFEM::PHYSICS;

namespace COMBUST
{

/*!
  \brief build standard domain integrals (matrix entries) for combustion or two-phase flow problem

  some remarks on the following sysmat entries:
  - this is an incremental formulation (details see e.g. Diss. Wall and BACI guide)
  - in comparison to standard fluid sysmat entries, here entries are scaled by the density
  - scaling has to be performed with quantities (density, viscosity) from corresponding time step (n, n+1)
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          size_t NUMDOF,
          size_t shpVecSize,
          size_t shpVecSizeStress>
void BuildStandardDomainIntegrals(
    LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>& assembler,
    const XFEM::ApproxFunc<2,shpVecSize>& shp,         // shape functions and their derivatives
    const double& fac,                                 // factor for integration in space
    const double& timefac,                             // factor for integration in time
    const double& timefacfac,                          // overall integration factor
    const double& densaf,                              // density at n+1 (OST), alpha_f (GenAlpha)
    const double& dynvisc,                             // dynamic viscosity at n+1
    const LINALG::Matrix<3,1>& velnp,                  // velocity at Gaussian point (GP) at n+1
    const double&              pres,                   // pressure at GP at n+1
    const LINALG::Matrix<3,1>& gradp,                  // pressure gradient at GP at n+1
    const LINALG::Matrix<3,3>& vderxy,                 // velocity derivative at GP at n+1
    const LINALG::Matrix<3,1>& rhsint,                 // time RHS (at n) and body force (at n+1) at GP
    const LINALG::Matrix<3,1>& res_old,                // residual at GP at n
    const LINALG::Matrix<3,1>& visc_old,               // viscous term at GP at n
    const LINALG::Matrix<shpVecSize,1>& enr_conv_c,    // convective term at GP at n
    const COMBUST::EnrViscs2<shpVecSize>& enr_viscs2,  // viscous term (second derivatves) at GP at n
    const bool instationary,                           // turn on/off instationary problem
    const bool newton,                                 // turn on/off linearization of convective term
    const bool pstab,                                  // turn on/off PSPG stabilization
    const bool supg,                                   // turn on/off SUPG stabilization
    const bool cstab,                                  // turn on/off BV   stabilization
    const double& tau_stab_Mu,                         // SUPG stabilization parameter
    const double& tau_stab_Mp,                         // PSPG stabilization parameter
    const double& tau_stab_C                           // BV   stabilization parameter
)
{
  const double densfac = densaf * fac;
  //----------------------------------------------------------------------
  //                            GALERKIN PART

#ifdef COMBUST_DOM_INERTIA
  if (instationary)
  {
    // inertia term (contribution to mass matrix)
    /*
                           /        \
                          |          |
                          |  v , Du  |
                          |          |
                           \        /
     */
    assembler.template Matrix<Velx,Velx>(shp.d0, densfac, shp.d0);
    assembler.template Matrix<Vely,Vely>(shp.d0, densfac, shp.d0);
    assembler.template Matrix<Velz,Velz>(shp.d0, densfac, shp.d0);

    assembler.template Vector<Velx>(shp.d0, -densfac*velnp(0));
    assembler.template Vector<Vely>(shp.d0, -densfac*velnp(1));
    assembler.template Vector<Velz>(shp.d0, -densfac*velnp(2));
  }
#endif

#ifdef COMBUST_DOM_CONVECTION
  // convection term, convective part
  /*
               /                       \
              |      / n+1       \      |
              | v , | u   o nabla | Du  |
              |      \ (i)       /      |
               \                       /
   */
  const double denstimefacfac = densaf * timefacfac;
  assembler.template Matrix<Velx,Velx>(shp.d0, denstimefacfac, enr_conv_c);
  assembler.template Matrix<Vely,Vely>(shp.d0, denstimefacfac, enr_conv_c);
  assembler.template Matrix<Velz,Velz>(shp.d0, denstimefacfac, enr_conv_c);

  assembler.template Vector<Velx>(shp.d0, -denstimefacfac*(velnp(0)*vderxy(0,0)
                                                          +velnp(1)*vderxy(0,1)
                                                          +velnp(2)*vderxy(0,2)));
  assembler.template Vector<Vely>(shp.d0, -denstimefacfac*(velnp(0)*vderxy(1,0)
                                                          +velnp(1)*vderxy(1,1)
                                                          +velnp(2)*vderxy(1,2)));
  assembler.template Vector<Velz>(shp.d0, -denstimefacfac*(velnp(0)*vderxy(2,0)
                                                          +velnp(1)*vderxy(2,1)
                                                          +velnp(2)*vderxy(2,2)));

  if (newton)
  {
    // convection term, reactive part
    /*
             /                         \
            |      /          \   n+1   |
            | v , | Du o nabla | u      |
            |      \          /   (i)   |
             \                         /
     */
    assembler.template Matrix<Velx,Velx>(shp.d0, denstimefacfac*vderxy(0,0), shp.d0);
    assembler.template Matrix<Velx,Vely>(shp.d0, denstimefacfac*vderxy(0,1), shp.d0);
    assembler.template Matrix<Velx,Velz>(shp.d0, denstimefacfac*vderxy(0,2), shp.d0);
    assembler.template Matrix<Vely,Velx>(shp.d0, denstimefacfac*vderxy(1,0), shp.d0);
    assembler.template Matrix<Vely,Vely>(shp.d0, denstimefacfac*vderxy(1,1), shp.d0);
    assembler.template Matrix<Vely,Velz>(shp.d0, denstimefacfac*vderxy(1,2), shp.d0);
    assembler.template Matrix<Velz,Velx>(shp.d0, denstimefacfac*vderxy(2,0), shp.d0);
    assembler.template Matrix<Velz,Vely>(shp.d0, denstimefacfac*vderxy(2,1), shp.d0);
    assembler.template Matrix<Velz,Velz>(shp.d0, denstimefacfac*vderxy(2,2), shp.d0);
  }
#endif

#ifdef COMBUST_DOM_DIFFUSION
  // viscous term
  /*
                /                                 \
               |       / \                  /  \   |
               |  eps | v | , 2\mu epsilon | Du |  |
               |       \ /                  \  /   |
                \                                 /
   */
  assembler.template Matrix<Velx,Velx>(shp.dx, 2.0*dynvisc*timefacfac, shp.dx);
  assembler.template Matrix<Velx,Velx>(shp.dy,     dynvisc*timefacfac, shp.dy);
  assembler.template Matrix<Velx,Vely>(shp.dy,     dynvisc*timefacfac, shp.dx);
  assembler.template Matrix<Velx,Velx>(shp.dz,     dynvisc*timefacfac, shp.dz);
  assembler.template Matrix<Velx,Velz>(shp.dz,     dynvisc*timefacfac, shp.dx);

  assembler.template Matrix<Vely,Vely>(shp.dx,     dynvisc*timefacfac, shp.dx);
  assembler.template Matrix<Vely,Velx>(shp.dx,     dynvisc*timefacfac, shp.dy);
  assembler.template Matrix<Vely,Vely>(shp.dy, 2.0*dynvisc*timefacfac, shp.dy);
  assembler.template Matrix<Vely,Vely>(shp.dz,     dynvisc*timefacfac, shp.dz);
  assembler.template Matrix<Vely,Velz>(shp.dz,     dynvisc*timefacfac, shp.dy);

  assembler.template Matrix<Velz,Velz>(shp.dx,     dynvisc*timefacfac, shp.dx);
  assembler.template Matrix<Velz,Velx>(shp.dx,     dynvisc*timefacfac, shp.dz);
  assembler.template Matrix<Velz,Velz>(shp.dy,     dynvisc*timefacfac, shp.dy);
  assembler.template Matrix<Velz,Vely>(shp.dy,     dynvisc*timefacfac, shp.dz);
  assembler.template Matrix<Velz,Velz>(shp.dz, 2.0*dynvisc*timefacfac, shp.dz);

  assembler.template Vector<Velx>(shp.dx,     -dynvisc*timefacfac*(vderxy(0, 0) + vderxy(0, 0)));
  assembler.template Vector<Velx>(shp.dy,     -dynvisc*timefacfac*(vderxy(0, 1) + vderxy(1, 0)));
  assembler.template Vector<Velx>(shp.dz,     -dynvisc*timefacfac*(vderxy(0, 2) + vderxy(2, 0)));

  assembler.template Vector<Vely>(shp.dx,     -dynvisc*timefacfac*(vderxy(1, 0) + vderxy(0, 1)));
  assembler.template Vector<Vely>(shp.dy,     -dynvisc*timefacfac*(vderxy(1, 1) + vderxy(1, 1)));
  assembler.template Vector<Vely>(shp.dz,     -dynvisc*timefacfac*(vderxy(1, 2) + vderxy(2, 1)));

  assembler.template Vector<Velz>(shp.dx,     -dynvisc*timefacfac*(vderxy(2, 0) + vderxy(0, 2)));
  assembler.template Vector<Velz>(shp.dy,     -dynvisc*timefacfac*(vderxy(2, 1) + vderxy(1, 2)));
  assembler.template Vector<Velz>(shp.dz,     -dynvisc*timefacfac*(vderxy(2, 2) + vderxy(2, 2)));
#endif

#ifdef COMBUST_DOM_PRESSURE
  // pressure term
  /*
                  /                \
                 |                  |
               - |  nabla o v , Dp  |
                 |                  |
                  \                /
   */
  assembler.template Matrix<Velx,Pres>(shp.dx, -timefacfac, shp.d0);
  assembler.template Matrix<Vely,Pres>(shp.dy, -timefacfac, shp.d0);
  assembler.template Matrix<Velz,Pres>(shp.dz, -timefacfac, shp.d0);

  assembler.template Vector<Velx>(shp.dx, timefacfac*pres);
  assembler.template Vector<Vely>(shp.dy, timefacfac*pres);
  assembler.template Vector<Velz>(shp.dz, timefacfac*pres);
#endif

#ifdef COMBUST_DOM_CONTINUITY
  // cout << "/!\\ warning === sign of continuity term has been changed to make formulation symmetric" << endl;
  // solenoidality term - continuity equation
  /*
                 /              \
                |                |
                | q , nabla o Du |
                |                |
                 \              /
   */
  assembler.template Matrix<Pres,Velx>(shp.d0, timefacfac, shp.dx);
  assembler.template Matrix<Pres,Vely>(shp.d0, timefacfac, shp.dy);
  assembler.template Matrix<Pres,Velz>(shp.d0, timefacfac, shp.dz);

  const double trace_gamma = (vderxy(0,0) + vderxy(1,1) + vderxy(2,2));
  assembler.template Vector<Pres>(shp.d0, -timefacfac*trace_gamma);
#endif

#ifdef COMBUST_DOM_RHS
  // source term of the right hand side
  /*
                  /    \
                 |      |
                 | v, f |
                 |      |
                  \    /
   */
  assembler.template Vector<Velx>(shp.d0, fac*rhsint(0));
  assembler.template Vector<Vely>(shp.d0, fac*rhsint(1));
  assembler.template Vector<Velz>(shp.d0, fac*rhsint(2));
#endif

  //----------------------------------------------------------------------
  //                 PRESSURE STABILISATION PART
  if(pstab)
  {
    const double tauMp = tau_stab_Mp * fac;
    const double denstauMp = densaf * tauMp;
#ifdef COMBUST_PSPG_INERTIA
    if (instationary)
    {
      /* pressure stabilisation: inertia */
      /*
                      /              \
                     |                |
                     |  Du , nabla q  |
                     |                |
                      \              /
       */
      assembler.template Matrix<Pres,Velx>(shp.dx, denstauMp, shp.d0);
      assembler.template Matrix<Pres,Vely>(shp.dy, denstauMp, shp.d0);
      assembler.template Matrix<Pres,Velz>(shp.dz, denstauMp, shp.d0);
    }
#endif

    const double denstimetauMp = timefac * tau_stab_Mp * densfac;
#ifdef COMBUST_PSPG_CONVECTION
    /* pressure stabilisation: convection, convective part */
    /*
                /                             \
               |             / n+1       \     |
               | nabla q ,  | u   o nabla | Du |
               |             \ i         /     |
                \                             /
     */
    assembler.template Matrix<Pres,Velx>(shp.dx, denstimetauMp, enr_conv_c);
    assembler.template Matrix<Pres,Vely>(shp.dy, denstimetauMp, enr_conv_c);
    assembler.template Matrix<Pres,Velz>(shp.dz, denstimetauMp, enr_conv_c);

    if (newton)
    {
      /*  pressure stabilisation: convection, reactive part
                /                             \
               |           /          \   n+1  |
               | grad q , | Du o nabla | u     |
               |           \          /   (i)  |
                \                             /
       */
      assembler.template Matrix<Pres,Velx>(shp.dx, denstimetauMp*vderxy(0,0), shp.d0);
      assembler.template Matrix<Pres,Velx>(shp.dy, denstimetauMp*vderxy(1,0), shp.d0);
      assembler.template Matrix<Pres,Velx>(shp.dz, denstimetauMp*vderxy(2,0), shp.d0);

      assembler.template Matrix<Pres,Vely>(shp.dx, denstimetauMp*vderxy(0,1), shp.d0);
      assembler.template Matrix<Pres,Vely>(shp.dy, denstimetauMp*vderxy(1,1), shp.d0);
      assembler.template Matrix<Pres,Vely>(shp.dz, denstimetauMp*vderxy(2,1), shp.d0);

      assembler.template Matrix<Pres,Velz>(shp.dx, denstimetauMp*vderxy(0,2), shp.d0);
      assembler.template Matrix<Pres,Velz>(shp.dy, denstimetauMp*vderxy(1,2), shp.d0);
      assembler.template Matrix<Pres,Velz>(shp.dz, denstimetauMp*vderxy(2,2), shp.d0);
    }
#endif

    const double timetauMp = timefac * tauMp;
#ifdef COMBUST_PSPG_DIFFUSION
    /* pressure stabilisation: viscosity (-L_visc_u) */
    /*
                 /                             \
                |                         /  \  |
              - |  nabla q , nabla o tau | Du | |
                |                         \  /  |
                 \                             /
     */
    assembler.template Matrix<Pres,Velx>(shp.dx, -2.0*dynvisc*timetauMp, enr_viscs2.xx);
    assembler.template Matrix<Pres,Vely>(shp.dx, -2.0*dynvisc*timetauMp, enr_viscs2.xy);
    assembler.template Matrix<Pres,Velz>(shp.dx, -2.0*dynvisc*timetauMp, enr_viscs2.xz);

    assembler.template Matrix<Pres,Velx>(shp.dy, -2.0*dynvisc*timetauMp, enr_viscs2.xy);
    assembler.template Matrix<Pres,Vely>(shp.dy, -2.0*dynvisc*timetauMp, enr_viscs2.yy);
    assembler.template Matrix<Pres,Velz>(shp.dy, -2.0*dynvisc*timetauMp, enr_viscs2.yz);

    assembler.template Matrix<Pres,Velx>(shp.dz, -2.0*dynvisc*timetauMp, enr_viscs2.xz);
    assembler.template Matrix<Pres,Vely>(shp.dz, -2.0*dynvisc*timetauMp, enr_viscs2.yz);
    assembler.template Matrix<Pres,Velz>(shp.dz, -2.0*dynvisc*timetauMp, enr_viscs2.zz);
#endif

#ifdef COMBUST_PSPG_PRESSURE
    /* pressure stabilisation: pressure( L_pres_p) */
    /*
                /                    \
               |                      |
               |  nabla q , nabla Dp  |
               |                      |
                \                    /
     */
    assembler.template Matrix<Pres,Pres>(shp.dx, timetauMp, shp.dx);
    assembler.template Matrix<Pres,Pres>(shp.dy, timetauMp, shp.dy);
    assembler.template Matrix<Pres,Pres>(shp.dz, timetauMp, shp.dz);
#endif

#ifdef COMBUST_PSPG_RHS
    // pressure stabilization
    assembler.template Vector<Pres>(shp.dx, -tauMp*res_old(0));
    assembler.template Vector<Pres>(shp.dy, -tauMp*res_old(1));
    assembler.template Vector<Pres>(shp.dz, -tauMp*res_old(2));
#endif
  }

  //----------------------------------------------------------------------
  //                     SUPG STABILISATION PART
  if(supg)
  {
    const double tauMu = tau_stab_Mu * fac;
    const double denstauMu = densaf * tauMu;
    const double densdenstauMu = densaf * denstauMu;
#ifdef COMBUST_SUPG_INERTIA
    if (instationary)
    {
      /* supg stabilisation: inertia  */
      /*
                    /                       \
                   |   / n+1       \         |
                   |  | u   o nabla | v, Du  |
                   |   \ (i)       /         |
                    \                       /
       */
      assembler.template Matrix<Velx,Velx>(enr_conv_c, densdenstauMu, shp.d0);
      assembler.template Matrix<Vely,Vely>(enr_conv_c, densdenstauMu, shp.d0);
      assembler.template Matrix<Velz,Velz>(enr_conv_c, densdenstauMu, shp.d0);

      if (newton)
      {
        /* supg stabilisation: inertia, linearisation of testfunction  */
        /*
                         /                           \
                        |   n+1      /          \     |
                        |  u      , | Du o nabla | v  |
                        |   (i)      \          /     |
                         \                           /
         */
        assembler.template Matrix<Velx,Velx>(shp.dx, densdenstauMu*velnp(0), shp.d0);
        assembler.template Matrix<Velx,Vely>(shp.dy, densdenstauMu*velnp(0), shp.d0);
        assembler.template Matrix<Velx,Velz>(shp.dz, densdenstauMu*velnp(0), shp.d0);
        assembler.template Matrix<Vely,Velx>(shp.dx, densdenstauMu*velnp(1), shp.d0);
        assembler.template Matrix<Vely,Vely>(shp.dy, densdenstauMu*velnp(1), shp.d0);
        assembler.template Matrix<Vely,Velz>(shp.dz, densdenstauMu*velnp(1), shp.d0);
        assembler.template Matrix<Velz,Velx>(shp.dx, densdenstauMu*velnp(2), shp.d0);
        assembler.template Matrix<Velz,Vely>(shp.dy, densdenstauMu*velnp(2), shp.d0);
        assembler.template Matrix<Velz,Velz>(shp.dz, densdenstauMu*velnp(2), shp.d0);
      }
    }
#endif

    const double denstimetauMu  = timefac * denstauMu;
    const double densdenstimetauMu  = timefac * densdenstauMu;
#ifdef COMBUST_SUPG_CONVECTION
    /* supg stabilisation: convective part ( L_conv_u) */
    /*
           /                                          \
          |  / n+1        \        / n+1        \      |
          | | u    o nabla | v ,  | u    o nabla | Du  |
          |  \ (i)        /        \ (i)        /      |
           \                                          /
     */
    assembler.template Matrix<Velx,Velx>(enr_conv_c, densdenstimetauMu, enr_conv_c);
    assembler.template Matrix<Vely,Vely>(enr_conv_c, densdenstimetauMu, enr_conv_c);
    assembler.template Matrix<Velz,Velz>(enr_conv_c, densdenstimetauMu, enr_conv_c);
#endif

#ifdef COMBUST_SUPG_PRESSURE
    /* supg stabilisation: pressure part  ( L_pres_p) */
    /*
                /                             \
               |   / n+1       \               |
               |  | u   o nabla | v , nabla Dp |
               |   \ (i)       /               |
                \                             /
     */
    assembler.template Matrix<Velx,Pres>(enr_conv_c, denstimetauMu, shp.dx);
    assembler.template Matrix<Vely,Pres>(enr_conv_c, denstimetauMu, shp.dy);
    assembler.template Matrix<Velz,Pres>(enr_conv_c, denstimetauMu, shp.dz);
#endif

#ifdef COMBUST_SUPG_DIFFUSION
    /* supg stabilisation: viscous part  (-L_visc_u) */
    /*
            /                                        \
           |               /  \    / n+1        \     |
         - |  nabla o eps | Du |, | u    o nabla | v  |
           |               \  /    \ (i)        /     |
            \                                        /
     */
    assembler.template Matrix<Velx,Velx>(enr_conv_c, -2.0*dynvisc*denstimetauMu, enr_viscs2.xx);
    assembler.template Matrix<Velx,Vely>(enr_conv_c, -2.0*dynvisc*denstimetauMu, enr_viscs2.xy);
    assembler.template Matrix<Velx,Velz>(enr_conv_c, -2.0*dynvisc*denstimetauMu, enr_viscs2.xz);

    assembler.template Matrix<Vely,Velx>(enr_conv_c, -2.0*dynvisc*denstimetauMu, enr_viscs2.yx);
    assembler.template Matrix<Vely,Vely>(enr_conv_c, -2.0*dynvisc*denstimetauMu, enr_viscs2.yy);
    assembler.template Matrix<Vely,Velz>(enr_conv_c, -2.0*dynvisc*denstimetauMu, enr_viscs2.yz);

    assembler.template Matrix<Velz,Velx>(enr_conv_c, -2.0*dynvisc*denstimetauMu, enr_viscs2.zx);
    assembler.template Matrix<Velz,Vely>(enr_conv_c, -2.0*dynvisc*denstimetauMu, enr_viscs2.zy);
    assembler.template Matrix<Velz,Velz>(enr_conv_c, -2.0*dynvisc*denstimetauMu, enr_viscs2.zz);
#endif

    if (newton)
    {
#ifdef COMBUST_SUPG_CONVECTION
      /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
      /*
                     /                                           \
                    |    /          \   n+1    / n+1        \     |
                    |   | Du o nabla | u    , | u    o nabla | v  |  \\ change order in documentation!
                    |    \          /   (i)    \ (i)        /     |
                     \                                           /
       */
      assembler.template Matrix<Velx,Velx>(enr_conv_c, densdenstimetauMu*vderxy(0,0), shp.d0);
      assembler.template Matrix<Velx,Vely>(enr_conv_c, densdenstimetauMu*vderxy(0,1), shp.d0);
      assembler.template Matrix<Velx,Velz>(enr_conv_c, densdenstimetauMu*vderxy(0,2), shp.d0);

      assembler.template Matrix<Vely,Velx>(enr_conv_c, densdenstimetauMu*vderxy(1,0), shp.d0);
      assembler.template Matrix<Vely,Vely>(enr_conv_c, densdenstimetauMu*vderxy(1,1), shp.d0);
      assembler.template Matrix<Vely,Velz>(enr_conv_c, densdenstimetauMu*vderxy(1,2), shp.d0);

      assembler.template Matrix<Velz,Velx>(enr_conv_c, densdenstimetauMu*vderxy(2,0), shp.d0);
      assembler.template Matrix<Velz,Vely>(enr_conv_c, densdenstimetauMu*vderxy(2,1), shp.d0);
      assembler.template Matrix<Velz,Velz>(enr_conv_c, densdenstimetauMu*vderxy(2,2), shp.d0);

      /*
                   /                                           \
                  |    / n+1        \   n+1    /          \     |
                  |   | u    o nabla | u    , | Du o nabla | v  |
                  |    \ (i)        /   (i)    \          /     |
                   \                                           /
       */
      const double con0 = densdenstimetauMu*(velnp(0)*vderxy(0,0) + velnp(1)*vderxy(0,1) + velnp(2)*vderxy(0,2));
      assembler.template Matrix<Velx,Velx>(shp.dx, con0, shp.d0);
      assembler.template Matrix<Velx,Vely>(shp.dy, con0, shp.d0);
      assembler.template Matrix<Velx,Velz>(shp.dz, con0, shp.d0);

      const double con1 = densdenstimetauMu*(velnp(0)*vderxy(1,0) + velnp(1)*vderxy(1,1) + velnp(2)*vderxy(1,2));
      assembler.template Matrix<Vely,Velx>(shp.dx, con1, shp.d0);
      assembler.template Matrix<Vely,Vely>(shp.dy, con1, shp.d0);
      assembler.template Matrix<Vely,Velz>(shp.dz, con1, shp.d0);

      const double con2 = densdenstimetauMu*(velnp(0)*vderxy(2,0) + velnp(1)*vderxy(2,1) + velnp(2)*vderxy(2,2));
      assembler.template Matrix<Velz,Velx>(shp.dx, con2, shp.d0);
      assembler.template Matrix<Velz,Vely>(shp.dy, con2, shp.d0);
      assembler.template Matrix<Velz,Velz>(shp.dz, con2, shp.d0);
#endif

#ifdef COMBUST_SUPG_PRESSURE
      /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
      /*
                          /                               \
                         |         n+1    /          \     |
                         |  nabla p    , | Du o nabla | v  |
                         |         (i)    \          /     |
                          \                               /
       */
      assembler.template Matrix<Velx,Velx>(shp.dx, denstimetauMu*gradp(0), shp.d0);
      assembler.template Matrix<Velx,Vely>(shp.dy, denstimetauMu*gradp(0), shp.d0);
      assembler.template Matrix<Velx,Velz>(shp.dz, denstimetauMu*gradp(0), shp.d0);

      assembler.template Matrix<Vely,Velx>(shp.dx, denstimetauMu*gradp(1), shp.d0);
      assembler.template Matrix<Vely,Vely>(shp.dy, denstimetauMu*gradp(1), shp.d0);
      assembler.template Matrix<Vely,Velz>(shp.dz, denstimetauMu*gradp(1), shp.d0);

      assembler.template Matrix<Velz,Velx>(shp.dx, denstimetauMu*gradp(2), shp.d0);
      assembler.template Matrix<Velz,Vely>(shp.dy, denstimetauMu*gradp(2), shp.d0);
      assembler.template Matrix<Velz,Velz>(shp.dz, denstimetauMu*gradp(2), shp.d0);
#endif

#ifdef COMBUST_SUPG_DIFFUSION
      /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
      /*
                    /                                         \
                   |               / n+1 \    /          \     |
                 - |  nabla o eps | u     |, | Du o nabla | v  |
                   |               \ (i) /    \          /     |
                    \                                         /
       */
      assembler.template Matrix<Velx,Velx>(shp.dx, -2.0*dynvisc*denstimetauMu*visc_old(0), shp.d0);
      assembler.template Matrix<Velx,Vely>(shp.dy, -2.0*dynvisc*denstimetauMu*visc_old(0), shp.d0);
      assembler.template Matrix<Velx,Velz>(shp.dz, -2.0*dynvisc*denstimetauMu*visc_old(0), shp.d0);

      assembler.template Matrix<Vely,Velx>(shp.dx, -2.0*dynvisc*denstimetauMu*visc_old(1), shp.d0);
      assembler.template Matrix<Vely,Vely>(shp.dy, -2.0*dynvisc*denstimetauMu*visc_old(1), shp.d0);
      assembler.template Matrix<Vely,Velz>(shp.dz, -2.0*dynvisc*denstimetauMu*visc_old(1), shp.d0);

      assembler.template Matrix<Velz,Velx>(shp.dx, -2.0*dynvisc*denstimetauMu*visc_old(2), shp.d0);
      assembler.template Matrix<Velz,Vely>(shp.dy, -2.0*dynvisc*denstimetauMu*visc_old(2), shp.d0);
      assembler.template Matrix<Velz,Velz>(shp.dz, -2.0*dynvisc*denstimetauMu*visc_old(2), shp.d0);
#endif

#ifdef COMBUST_SUPG_BODYFORCE
      /* supg stabilisation: bodyforce part, linearisation of test function */
      /*
                        /                             \
                       |              /          \     |
                     - |  rhsint   , | Du o nabla | v  |
                       |              \          /     |
                        \                             /
       */
      assembler.template Matrix<Velx,Velx>(shp.dx, -denstauMu*rhsint(0), shp.d0);
      assembler.template Matrix<Velx,Vely>(shp.dy, -denstauMu*rhsint(0), shp.d0);
      assembler.template Matrix<Velx,Velz>(shp.dz, -denstauMu*rhsint(0), shp.d0);

      assembler.template Matrix<Vely,Velx>(shp.dx, -denstauMu*rhsint(1), shp.d0);
      assembler.template Matrix<Vely,Vely>(shp.dy, -denstauMu*rhsint(1), shp.d0);
      assembler.template Matrix<Vely,Velz>(shp.dz, -denstauMu*rhsint(1), shp.d0);

      assembler.template Matrix<Velz,Velx>(shp.dx, -denstauMu*rhsint(2), shp.d0);
      assembler.template Matrix<Velz,Vely>(shp.dy, -denstauMu*rhsint(2), shp.d0);
      assembler.template Matrix<Velz,Velz>(shp.dz, -denstauMu*rhsint(2), shp.d0);
#endif
    } // if newton

#ifdef COMBUST_SUPG_RHS
    // supg stabilisation
    assembler.template Vector<Velx>(enr_conv_c, -denstauMu*res_old(0));
    assembler.template Vector<Vely>(enr_conv_c, -denstauMu*res_old(1));
    assembler.template Vector<Velz>(enr_conv_c, -denstauMu*res_old(2));
#endif
  }


  //----------------------------------------------------------------------
  //                     STABILISATION, CONTINUITY PART
  if(cstab)
  {
#ifdef COMBUST_CSTAB
    const double timefac_tau_C = timefac * tau_stab_C * fac;
    const double timefac_tau_C_divunp = timefac_tau_C * (vderxy(0,0)+vderxy(1,1)+vderxy(2,2));
    /* continuity stabilisation on left hand side */
    /*
               /                        \
              |                          |
              | nabla o Du  , nabla o v  |
              |                          |
               \                        /
     */
    assembler.template Matrix<Velx,Velx>(shp.dx, timefac_tau_C, shp.dx);
    assembler.template Matrix<Velx,Vely>(shp.dx, timefac_tau_C, shp.dy);
    assembler.template Matrix<Velx,Velz>(shp.dx, timefac_tau_C, shp.dz);

    assembler.template Matrix<Vely,Velx>(shp.dy, timefac_tau_C, shp.dx);
    assembler.template Matrix<Vely,Vely>(shp.dy, timefac_tau_C, shp.dy);
    assembler.template Matrix<Vely,Velz>(shp.dy, timefac_tau_C, shp.dz);

    assembler.template Matrix<Velz,Velx>(shp.dz, timefac_tau_C, shp.dx);
    assembler.template Matrix<Velz,Vely>(shp.dz, timefac_tau_C, shp.dy);
    assembler.template Matrix<Velz,Velz>(shp.dz, timefac_tau_C, shp.dz);

    assembler.template Vector<Velx>(shp.dx, -timefac_tau_C_divunp);
    assembler.template Vector<Vely>(shp.dy, -timefac_tau_C_divunp);
    assembler.template Vector<Velz>(shp.dz, -timefac_tau_C_divunp);
#endif
  } // endif cstab
}

} // namespace COMBUST

#endif
#endif
#endif
