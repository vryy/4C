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
          size_t shpVecSize>
void BuildStandardDomainIntegrals(
    COMBUST::LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>& assembler,
    const XFEM::ApproxFunc<2,shpVecSize>& shpvel,      // shape functions and their derivatives (velocity)
    const XFEM::ApproxFunc<2,shpVecSize>& shppres,     // shape functions and their derivatives (pressure)
    const double& fac,                                 // factor for integration in space
    const double& timefac,                             // factor for integration in time
    const double& timefacfac,                          // overall integration factor
    const double& alphaM,
    const double& densaf,                              // density at n+1 (OST), alpha_f (GenAlpha)
    const double& dynvisc,                             // dynamic viscosity at n+1
    const LINALG::Matrix<3,1>& velnp,                  // velocity at Gaussian point (GP) at n+1
    const double&              pres,                   // pressure at GP at n+1
    const LINALG::Matrix<3,1>& accam,                  // acceleration at n+alphaM
    const LINALG::Matrix<3,1>& gradp,                  // pressure gradient at GP at n+1
    const LINALG::Matrix<3,3>& vderxy,                 // velocity derivative at GP at n+1
    const LINALG::Matrix<3,1>& rhsint,                 // time RHS (at n) and body force (at n+1) at GP
    const LINALG::Matrix<3,1>& res_old,                // residual at GP at n
    const LINALG::Matrix<3,1>& visc_old,               // viscous term at GP at n
    const LINALG::Matrix<shpVecSize,1>& enr_conv_c,    // convective term at GP at n
    const LINALG::Matrix<shpVecSize,1>& enr_conv_c_sg, // convective subgrid term at GP at n
    const COMBUST::EnrViscs2<shpVecSize>& enr_viscs2,  // viscous term (second derivatves) at GP at n
    const bool instationary,                           // turn on/off instationary problem
    const bool genalpha,                               // switch to generalized alpha terms
    const bool newton,                                 // turn on/off linearization of convective term
    const bool pstab,                                  // turn on/off PSPG stabilization
    const bool supg,                                   // turn on/off SUPG stabilization
    const bool cstab,                                  // turn on/off BV   stabilization
    const double& tau_stab_Mu,                         // SUPG stabilization parameter
    const double& tau_stab_Mp,                         // PSPG stabilization parameter
    const double& tau_stab_C                           // BV   stabilization parameter
)
{
  const double rhstimefacfac = timefacfac/alphaM;
  const double densfac = densaf * fac;
  //----------------------------------------------------------------------
  //                            GALERKIN PART

#ifdef COMBUST_DOM_INERTIA
  if (instationary)
  {
    // inertia term (contribution to mass matrix)
      /*
                                   /           \
                                  |             |
                                  | roh Du , v  |
                                  |             |
                                   \           /
       */
    assembler.template Matrix<Velx,Velx>(shpvel.d0, densfac, shpvel.d0);
    assembler.template Matrix<Vely,Vely>(shpvel.d0, densfac, shpvel.d0);
    assembler.template Matrix<Velz,Velz>(shpvel.d0, densfac, shpvel.d0);

    if (genalpha)
    {
      const bool densrhstimefacfac = densaf * rhstimefacfac;
      // TODO const bool densrhstimefacfac = densam * rhstimefacfac;

      assembler.template Vector<Velx>(shpvel.d0, -densrhstimefacfac*accam(0));
      assembler.template Vector<Vely>(shpvel.d0, -densrhstimefacfac*accam(1));
      assembler.template Vector<Velz>(shpvel.d0, -densrhstimefacfac*accam(2));
    }
    else
    {
      assembler.template Vector<Velx>(shpvel.d0, -densfac*velnp(0));
      assembler.template Vector<Vely>(shpvel.d0, -densfac*velnp(1));
      assembler.template Vector<Velz>(shpvel.d0, -densfac*velnp(2));
    }
  }
#endif

#ifdef COMBUST_DOM_CONVECTION
  // convection term, convective part
  /* convection, convective part */
  /*
                     /                                                      \
                    |         / n+1       \                                  |
                    | v , roh| u   o nabla | Du - rho*(tau_M*R_M o nabla) Du |
                    |         \ (i)       /          subgrid velocity        |
                     \                            cross-stress              /
   */
  const double denstimefacfac = densaf * timefacfac;
  const double rhsdenstimefacfac = densaf * rhstimefacfac;
  assembler.template Matrix<Velx,Velx>(shpvel.d0, denstimefacfac, enr_conv_c_sg);
  assembler.template Matrix<Vely,Vely>(shpvel.d0, denstimefacfac, enr_conv_c_sg);
  assembler.template Matrix<Velz,Velz>(shpvel.d0, denstimefacfac, enr_conv_c_sg);

  assembler.template Vector<Velx>(shpvel.d0, -rhsdenstimefacfac*(velnp(0)*vderxy(0,0)
                                                          +velnp(1)*vderxy(0,1)
                                                          +velnp(2)*vderxy(0,2)));
  assembler.template Vector<Vely>(shpvel.d0, -rhsdenstimefacfac*(velnp(0)*vderxy(1,0)
                                                          +velnp(1)*vderxy(1,1)
                                                          +velnp(2)*vderxy(1,2)));
  assembler.template Vector<Velz>(shpvel.d0, -rhsdenstimefacfac*(velnp(0)*vderxy(2,0)
                                                          +velnp(1)*vderxy(2,1)
                                                          +velnp(2)*vderxy(2,2)));

  if (newton)
  {
    // convection term, reactive part
      /*
                     /                             \
                    |         /          \   n+1   |
                    | v ,roh | Du o nabla | u      |
                    |         \          /   (i)   |
                     \                             /
       */
    assembler.template Matrix<Velx,Velx>(shpvel.d0, denstimefacfac*vderxy(0,0), shpvel.d0);
    assembler.template Matrix<Velx,Vely>(shpvel.d0, denstimefacfac*vderxy(0,1), shpvel.d0);
    assembler.template Matrix<Velx,Velz>(shpvel.d0, denstimefacfac*vderxy(0,2), shpvel.d0);
    assembler.template Matrix<Vely,Velx>(shpvel.d0, denstimefacfac*vderxy(1,0), shpvel.d0);
    assembler.template Matrix<Vely,Vely>(shpvel.d0, denstimefacfac*vderxy(1,1), shpvel.d0);
    assembler.template Matrix<Vely,Velz>(shpvel.d0, denstimefacfac*vderxy(1,2), shpvel.d0);
    assembler.template Matrix<Velz,Velx>(shpvel.d0, denstimefacfac*vderxy(2,0), shpvel.d0);
    assembler.template Matrix<Velz,Vely>(shpvel.d0, denstimefacfac*vderxy(2,1), shpvel.d0);
    assembler.template Matrix<Velz,Velz>(shpvel.d0, denstimefacfac*vderxy(2,2), shpvel.d0);
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
  assembler.template Matrix<Velx,Velx>(shpvel.dx, 2.0*dynvisc*timefacfac, shpvel.dx);
  assembler.template Matrix<Velx,Velx>(shpvel.dy,     dynvisc*timefacfac, shpvel.dy);
  assembler.template Matrix<Velx,Vely>(shpvel.dy,     dynvisc*timefacfac, shpvel.dx);
  assembler.template Matrix<Velx,Velx>(shpvel.dz,     dynvisc*timefacfac, shpvel.dz);
  assembler.template Matrix<Velx,Velz>(shpvel.dz,     dynvisc*timefacfac, shpvel.dx);

  assembler.template Matrix<Vely,Vely>(shpvel.dx,     dynvisc*timefacfac, shpvel.dx);
  assembler.template Matrix<Vely,Velx>(shpvel.dx,     dynvisc*timefacfac, shpvel.dy);
  assembler.template Matrix<Vely,Vely>(shpvel.dy, 2.0*dynvisc*timefacfac, shpvel.dy);
  assembler.template Matrix<Vely,Vely>(shpvel.dz,     dynvisc*timefacfac, shpvel.dz);
  assembler.template Matrix<Vely,Velz>(shpvel.dz,     dynvisc*timefacfac, shpvel.dy);

  assembler.template Matrix<Velz,Velz>(shpvel.dx,     dynvisc*timefacfac, shpvel.dx);
  assembler.template Matrix<Velz,Velx>(shpvel.dx,     dynvisc*timefacfac, shpvel.dz);
  assembler.template Matrix<Velz,Velz>(shpvel.dy,     dynvisc*timefacfac, shpvel.dy);
  assembler.template Matrix<Velz,Vely>(shpvel.dy,     dynvisc*timefacfac, shpvel.dz);
  assembler.template Matrix<Velz,Velz>(shpvel.dz, 2.0*dynvisc*timefacfac, shpvel.dz);

  assembler.template Vector<Velx>(shpvel.dx,     -dynvisc*rhstimefacfac*(vderxy(0, 0) + vderxy(0, 0)));
  assembler.template Vector<Velx>(shpvel.dy,     -dynvisc*rhstimefacfac*(vderxy(0, 1) + vderxy(1, 0)));
  assembler.template Vector<Velx>(shpvel.dz,     -dynvisc*rhstimefacfac*(vderxy(0, 2) + vderxy(2, 0)));

  assembler.template Vector<Vely>(shpvel.dx,     -dynvisc*rhstimefacfac*(vderxy(1, 0) + vderxy(0, 1)));
  assembler.template Vector<Vely>(shpvel.dy,     -dynvisc*rhstimefacfac*(vderxy(1, 1) + vderxy(1, 1)));
  assembler.template Vector<Vely>(shpvel.dz,     -dynvisc*rhstimefacfac*(vderxy(1, 2) + vderxy(2, 1)));

  assembler.template Vector<Velz>(shpvel.dx,     -dynvisc*rhstimefacfac*(vderxy(2, 0) + vderxy(0, 2)));
  assembler.template Vector<Velz>(shpvel.dy,     -dynvisc*rhstimefacfac*(vderxy(2, 1) + vderxy(1, 2)));
  assembler.template Vector<Velz>(shpvel.dz,     -dynvisc*rhstimefacfac*(vderxy(2, 2) + vderxy(2, 2)));
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
  assembler.template Matrix<Velx,Pres>(shpvel.dx, -timefacfac, shppres.d0);
  assembler.template Matrix<Vely,Pres>(shpvel.dy, -timefacfac, shppres.d0);
  assembler.template Matrix<Velz,Pres>(shpvel.dz, -timefacfac, shppres.d0);

  assembler.template Vector<Velx>(shpvel.dx, rhstimefacfac*pres);
  assembler.template Vector<Vely>(shpvel.dy, rhstimefacfac*pres);
  assembler.template Vector<Velz>(shpvel.dz, rhstimefacfac*pres);
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
  assembler.template Matrix<Pres,Velx>(shppres.d0, timefacfac, shpvel.dx);
  assembler.template Matrix<Pres,Vely>(shppres.d0, timefacfac, shpvel.dy);
  assembler.template Matrix<Pres,Velz>(shppres.d0, timefacfac, shpvel.dz);

  const double trace_gamma = (vderxy(0,0) + vderxy(1,1) + vderxy(2,2));
  assembler.template Vector<Pres>(shppres.d0, -rhstimefacfac*trace_gamma);
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
  if (genalpha)
  {
    const double rhsfac = fac/alphaM;
    assembler.template Vector<Velx>(shpvel.d0, rhsfac*rhsint(0));
    assembler.template Vector<Vely>(shpvel.d0, rhsfac*rhsint(1));
    assembler.template Vector<Velz>(shpvel.d0, rhsfac*rhsint(2));
  }
  else
  {
    assembler.template Vector<Velx>(shpvel.d0, fac*rhsint(0));
    assembler.template Vector<Vely>(shpvel.d0, fac*rhsint(1));
    assembler.template Vector<Velz>(shpvel.d0, fac*rhsint(2));
  }
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
                              /                 \
                             |                   |
                             | roh Du , nabla q  |
                             |                   |
                              \                 /
      */
      assembler.template Matrix<Pres,Velx>(shppres.dx, denstauMp, shpvel.d0);
      assembler.template Matrix<Pres,Vely>(shppres.dy, denstauMp, shpvel.d0);
      assembler.template Matrix<Pres,Velz>(shppres.dz, denstauMp, shpvel.d0);
    }
#endif

    const double denstimetauMp = timefac * tau_stab_Mp * densfac;
#ifdef COMBUST_PSPG_CONVECTION
    /* pressure stabilisation: convection, convective part */
    /*
                      /                                \
                     |                / n+1       \     |
                     | nabla q ,roh  | u   o nabla | Du |
                     |                \ i         /     |
                      \                                /
     */
    assembler.template Matrix<Pres,Velx>(shppres.dx, denstimetauMp, enr_conv_c);
    assembler.template Matrix<Pres,Vely>(shppres.dy, denstimetauMp, enr_conv_c);
    assembler.template Matrix<Pres,Velz>(shppres.dz, denstimetauMp, enr_conv_c);

    if (newton)
    {
        /*  pressure stabilisation: convection, reactive part
                        /                                \
                       |              /          \   n+1  |
                       | grad q , roh| Du o nabla | u     |
                       |              \          /   (i)  |
                        \                                /
         */
      assembler.template Matrix<Pres,Velx>(shppres.dx, denstimetauMp*vderxy(0,0), shpvel.d0);
      assembler.template Matrix<Pres,Velx>(shppres.dy, denstimetauMp*vderxy(1,0), shpvel.d0);
      assembler.template Matrix<Pres,Velx>(shppres.dz, denstimetauMp*vderxy(2,0), shpvel.d0);

      assembler.template Matrix<Pres,Vely>(shppres.dx, denstimetauMp*vderxy(0,1), shpvel.d0);
      assembler.template Matrix<Pres,Vely>(shppres.dy, denstimetauMp*vderxy(1,1), shpvel.d0);
      assembler.template Matrix<Pres,Vely>(shppres.dz, denstimetauMp*vderxy(2,1), shpvel.d0);

      assembler.template Matrix<Pres,Velz>(shppres.dx, denstimetauMp*vderxy(0,2), shpvel.d0);
      assembler.template Matrix<Pres,Velz>(shppres.dy, denstimetauMp*vderxy(1,2), shpvel.d0);
      assembler.template Matrix<Pres,Velz>(shppres.dz, denstimetauMp*vderxy(2,2), shpvel.d0);
    }
#endif

    const double timetauMp = timefac * tauMp;
#ifdef COMBUST_PSPG_DIFFUSION
    /* pressure stabilisation: viscosity (-L_visc_u) */
    /*
                       /                                      \
                      |                                  /  \  |
                    - |  nabla q , nabla o 2/mu epsilon | Du | |
                      |                                  \  /  |
                       \                                      /
     */
    assembler.template Matrix<Pres,Velx>(shppres.dx, -2.0*dynvisc*timetauMp, enr_viscs2.xx);
    assembler.template Matrix<Pres,Vely>(shppres.dx, -2.0*dynvisc*timetauMp, enr_viscs2.xy);
    assembler.template Matrix<Pres,Velz>(shppres.dx, -2.0*dynvisc*timetauMp, enr_viscs2.xz);

    assembler.template Matrix<Pres,Velx>(shppres.dy, -2.0*dynvisc*timetauMp, enr_viscs2.yx);
    assembler.template Matrix<Pres,Vely>(shppres.dy, -2.0*dynvisc*timetauMp, enr_viscs2.yy);
    assembler.template Matrix<Pres,Velz>(shppres.dy, -2.0*dynvisc*timetauMp, enr_viscs2.yz);

    assembler.template Matrix<Pres,Velx>(shppres.dz, -2.0*dynvisc*timetauMp, enr_viscs2.zx);
    assembler.template Matrix<Pres,Vely>(shppres.dz, -2.0*dynvisc*timetauMp, enr_viscs2.zy);
    assembler.template Matrix<Pres,Velz>(shppres.dz, -2.0*dynvisc*timetauMp, enr_viscs2.zz);
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
    assembler.template Matrix<Pres,Pres>(shppres.dx, timetauMp, shppres.dx);
    assembler.template Matrix<Pres,Pres>(shppres.dy, timetauMp, shppres.dy);
    assembler.template Matrix<Pres,Pres>(shppres.dz, timetauMp, shppres.dz);
#endif

#ifdef COMBUST_PSPG_RHS
    // pressure stabilization
    assembler.template Vector<Pres>(shppres.dx, -tauMp*res_old(0));
    assembler.template Vector<Pres>(shppres.dy, -tauMp*res_old(1));
    assembler.template Vector<Pres>(shppres.dz, -tauMp*res_old(2));
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
                             /                                                        \
                            |              / n+1       \                               |
                            | roh Du , roh| u   o nabla | v - rho*(tau_M*R_M o nabla)v |
                            |              \ (i)       /           subgrid velocity    |
                             \                                    reynolds-stress     /
      */
      assembler.template Matrix<Velx,Velx>(enr_conv_c_sg, densdenstauMu, shpvel.d0);
      assembler.template Matrix<Vely,Vely>(enr_conv_c_sg, densdenstauMu, shpvel.d0);
      assembler.template Matrix<Velz,Velz>(enr_conv_c_sg, densdenstauMu, shpvel.d0);

      if (newton)
      {
        /* supg stabilisation: inertia, linearisation of testfunction  */
          /*
                                 /                                  \
                                |      n+1          /          \     |
                                | roh u      , roh | Du o nabla | v  |
                                |      (i)          \          /     |
                                 \                                   /

           */
        assembler.template Matrix<Velx,Velx>(shpvel.dx, densdenstauMu*velnp(0), shpvel.d0);
        assembler.template Matrix<Velx,Vely>(shpvel.dy, densdenstauMu*velnp(0), shpvel.d0);
        assembler.template Matrix<Velx,Velz>(shpvel.dz, densdenstauMu*velnp(0), shpvel.d0);
        assembler.template Matrix<Vely,Velx>(shpvel.dx, densdenstauMu*velnp(1), shpvel.d0);
        assembler.template Matrix<Vely,Vely>(shpvel.dy, densdenstauMu*velnp(1), shpvel.d0);
        assembler.template Matrix<Vely,Velz>(shpvel.dz, densdenstauMu*velnp(1), shpvel.d0);
        assembler.template Matrix<Velz,Velx>(shpvel.dx, densdenstauMu*velnp(2), shpvel.d0);
        assembler.template Matrix<Velz,Vely>(shpvel.dy, densdenstauMu*velnp(2), shpvel.d0);
        assembler.template Matrix<Velz,Velz>(shpvel.dz, densdenstauMu*velnp(2), shpvel.d0);
      }
    }
#endif

    const double denstimetauMu  = timefac * denstauMu;
    const double densdenstimetauMu  = timefac * densdenstauMu;
#ifdef COMBUST_SUPG_CONVECTION
    /* supg stabilisation: convective part ( L_conv_u) */
    /*
                 /                                                                           \
                |     / n+1        \                                      / n+1        \      |
                | roh| u    o nabla | v - rho*(tau_M*R_M o nabla)v , roh | u    o nabla | Du  |
                |     \ (i)        /         subgrid velocity             \ (i)        /      |
                 \                           reynolds-stress                                 /
     */
    assembler.template Matrix<Velx,Velx>(enr_conv_c_sg, densdenstimetauMu, enr_conv_c);
    assembler.template Matrix<Vely,Vely>(enr_conv_c_sg, densdenstimetauMu, enr_conv_c);
    assembler.template Matrix<Velz,Velz>(enr_conv_c_sg, densdenstimetauMu, enr_conv_c);
#endif

#ifdef COMBUST_SUPG_PRESSURE
    /* supg stabilisation: pressure part  ( L_pres_p) */
    /*
                      /                                                           \
                     |      / n+1       \                                          |
                     | roh | u   o nabla | v - rho*(tau_M*R_M o nabla)v , nabla Dp |
                     |      \ (i)       /         subgrid velocity                 |
                      \                            reynolds-stress                /
     */
    assembler.template Matrix<Velx,Pres>(enr_conv_c_sg, denstimetauMu, shppres.dx);
    assembler.template Matrix<Vely,Pres>(enr_conv_c_sg, denstimetauMu, shppres.dy);
    assembler.template Matrix<Velz,Pres>(enr_conv_c_sg, denstimetauMu, shppres.dz);
#endif

#ifdef COMBUST_SUPG_DIFFUSION
    /* supg stabilisation: viscous part  (-L_visc_u) */
    /*
                  /                                                                     \
                 |               /  \       / n+1        \                               |
               - |  nabla o eps | Du |, roh| u    o nabla | v - rho*(tau_M*R_M o nabla)v |
                 |               \  /       \ (i)        /      subgrid velocity         |
                  \                                             reynolds-stress         /
     */
    assembler.template Matrix<Velx,Velx>(enr_conv_c_sg, -2.0*dynvisc*denstimetauMu, enr_viscs2.xx);
    assembler.template Matrix<Velx,Vely>(enr_conv_c_sg, -2.0*dynvisc*denstimetauMu, enr_viscs2.xy);
    assembler.template Matrix<Velx,Velz>(enr_conv_c_sg, -2.0*dynvisc*denstimetauMu, enr_viscs2.xz);

    assembler.template Matrix<Vely,Velx>(enr_conv_c_sg, -2.0*dynvisc*denstimetauMu, enr_viscs2.yx);
    assembler.template Matrix<Vely,Vely>(enr_conv_c_sg, -2.0*dynvisc*denstimetauMu, enr_viscs2.yy);
    assembler.template Matrix<Vely,Velz>(enr_conv_c_sg, -2.0*dynvisc*denstimetauMu, enr_viscs2.yz);

    assembler.template Matrix<Velz,Velx>(enr_conv_c_sg, -2.0*dynvisc*denstimetauMu, enr_viscs2.zx);
    assembler.template Matrix<Velz,Vely>(enr_conv_c_sg, -2.0*dynvisc*denstimetauMu, enr_viscs2.zy);
    assembler.template Matrix<Velz,Velz>(enr_conv_c_sg, -2.0*dynvisc*denstimetauMu, enr_viscs2.zz);
#endif

    if (newton)
    {
#ifdef COMBUST_SUPG_CONVECTION
      /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
      /*
                     /                                           \
                    |    /          \   n+1    / n+1        \     |
                    |roh| Du o nabla | u  ,roh| u    o nabla | v  |  \\ change order in documentation!
                    |    \          /   (i)    \ (i)        /     |
                     \                                           /
       */
      assembler.template Matrix<Velx,Velx>(enr_conv_c, densdenstimetauMu*vderxy(0,0), shpvel.d0);
      assembler.template Matrix<Velx,Vely>(enr_conv_c, densdenstimetauMu*vderxy(0,1), shpvel.d0);
      assembler.template Matrix<Velx,Velz>(enr_conv_c, densdenstimetauMu*vderxy(0,2), shpvel.d0);

      assembler.template Matrix<Vely,Velx>(enr_conv_c, densdenstimetauMu*vderxy(1,0), shpvel.d0);
      assembler.template Matrix<Vely,Vely>(enr_conv_c, densdenstimetauMu*vderxy(1,1), shpvel.d0);
      assembler.template Matrix<Vely,Velz>(enr_conv_c, densdenstimetauMu*vderxy(1,2), shpvel.d0);

      assembler.template Matrix<Velz,Velx>(enr_conv_c, densdenstimetauMu*vderxy(2,0), shpvel.d0);
      assembler.template Matrix<Velz,Vely>(enr_conv_c, densdenstimetauMu*vderxy(2,1), shpvel.d0);
      assembler.template Matrix<Velz,Velz>(enr_conv_c, densdenstimetauMu*vderxy(2,2), shpvel.d0);

      /*
                   /                                           \
                  |    / n+1        \   n+1    /          \     |
                  |roh| u    o nabla | u  ,roh| Du o nabla | v  |
                  |    \ (i)        /   (i)    \          /     |
                   \                                           /
       */
      const double con0 = densdenstimetauMu*(velnp(0)*vderxy(0,0) + velnp(1)*vderxy(0,1) + velnp(2)*vderxy(0,2));
      assembler.template Matrix<Velx,Velx>(shpvel.dx, con0, shpvel.d0);
      assembler.template Matrix<Velx,Vely>(shpvel.dy, con0, shpvel.d0);
      assembler.template Matrix<Velx,Velz>(shpvel.dz, con0, shpvel.d0);

      const double con1 = densdenstimetauMu*(velnp(0)*vderxy(1,0) + velnp(1)*vderxy(1,1) + velnp(2)*vderxy(1,2));
      assembler.template Matrix<Vely,Velx>(shpvel.dx, con1, shpvel.d0);
      assembler.template Matrix<Vely,Vely>(shpvel.dy, con1, shpvel.d0);
      assembler.template Matrix<Vely,Velz>(shpvel.dz, con1, shpvel.d0);

      const double con2 = densdenstimetauMu*(velnp(0)*vderxy(2,0) + velnp(1)*vderxy(2,1) + velnp(2)*vderxy(2,2));
      assembler.template Matrix<Velz,Velx>(shpvel.dx, con2, shpvel.d0);
      assembler.template Matrix<Velz,Vely>(shpvel.dy, con2, shpvel.d0);
      assembler.template Matrix<Velz,Velz>(shpvel.dz, con2, shpvel.d0);
#endif

#ifdef COMBUST_SUPG_PRESSURE
      /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
      /*
                                /                                  \
                               |         n+1       /          \     |
                               |  nabla p    ,roh | Du o nabla | v  |
                               |         (i)       \          /     |
                                \                                  /
       */
      assembler.template Matrix<Velx,Velx>(shpvel.dx, denstimetauMu*gradp(0), shpvel.d0);
      assembler.template Matrix<Velx,Vely>(shpvel.dy, denstimetauMu*gradp(0), shpvel.d0);
      assembler.template Matrix<Velx,Velz>(shpvel.dz, denstimetauMu*gradp(0), shpvel.d0);

      assembler.template Matrix<Vely,Velx>(shpvel.dx, denstimetauMu*gradp(1), shpvel.d0);
      assembler.template Matrix<Vely,Vely>(shpvel.dy, denstimetauMu*gradp(1), shpvel.d0);
      assembler.template Matrix<Vely,Velz>(shpvel.dz, denstimetauMu*gradp(1), shpvel.d0);

      assembler.template Matrix<Velz,Velx>(shpvel.dx, denstimetauMu*gradp(2), shpvel.d0);
      assembler.template Matrix<Velz,Vely>(shpvel.dy, denstimetauMu*gradp(2), shpvel.d0);
      assembler.template Matrix<Velz,Velz>(shpvel.dz, denstimetauMu*gradp(2), shpvel.d0);
#endif

#ifdef COMBUST_SUPG_DIFFUSION
      /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
      /*
                          /                                            \
                         |               / n+1 \       /          \     |
                       - |  nabla o eps | u     |,roh | Du o nabla | v  |
                         |               \ (i) /       \          /     |
                          \                                            /
       */
      assembler.template Matrix<Velx,Velx>(shpvel.dx, -2.0*dynvisc*denstimetauMu*visc_old(0), shpvel.d0);
      assembler.template Matrix<Velx,Vely>(shpvel.dy, -2.0*dynvisc*denstimetauMu*visc_old(0), shpvel.d0);
      assembler.template Matrix<Velx,Velz>(shpvel.dz, -2.0*dynvisc*denstimetauMu*visc_old(0), shpvel.d0);

      assembler.template Matrix<Vely,Velx>(shpvel.dx, -2.0*dynvisc*denstimetauMu*visc_old(1), shpvel.d0);
      assembler.template Matrix<Vely,Vely>(shpvel.dy, -2.0*dynvisc*denstimetauMu*visc_old(1), shpvel.d0);
      assembler.template Matrix<Vely,Velz>(shpvel.dz, -2.0*dynvisc*denstimetauMu*visc_old(1), shpvel.d0);

      assembler.template Matrix<Velz,Velx>(shpvel.dx, -2.0*dynvisc*denstimetauMu*visc_old(2), shpvel.d0);
      assembler.template Matrix<Velz,Vely>(shpvel.dy, -2.0*dynvisc*denstimetauMu*visc_old(2), shpvel.d0);
      assembler.template Matrix<Velz,Velz>(shpvel.dz, -2.0*dynvisc*denstimetauMu*visc_old(2), shpvel.d0);
#endif

#ifdef COMBUST_SUPG_BODYFORCE
      /* supg stabilisation: bodyforce part, linearisation of test function */
      /*
                              /                                \
                             |                 /          \     |
                           - |  rhsint   , roh| Du o nabla | v  |
                             |                 \          /     |
                              \                                /

       */
      assembler.template Matrix<Velx,Velx>(shpvel.dx, -denstauMu*rhsint(0), shpvel.d0);
      assembler.template Matrix<Velx,Vely>(shpvel.dy, -denstauMu*rhsint(0), shpvel.d0);
      assembler.template Matrix<Velx,Velz>(shpvel.dz, -denstauMu*rhsint(0), shpvel.d0);

      assembler.template Matrix<Vely,Velx>(shpvel.dx, -denstauMu*rhsint(1), shpvel.d0);
      assembler.template Matrix<Vely,Vely>(shpvel.dy, -denstauMu*rhsint(1), shpvel.d0);
      assembler.template Matrix<Vely,Velz>(shpvel.dz, -denstauMu*rhsint(1), shpvel.d0);

      assembler.template Matrix<Velz,Velx>(shpvel.dx, -denstauMu*rhsint(2), shpvel.d0);
      assembler.template Matrix<Velz,Vely>(shpvel.dy, -denstauMu*rhsint(2), shpvel.d0);
      assembler.template Matrix<Velz,Velz>(shpvel.dz, -denstauMu*rhsint(2), shpvel.d0);
#endif
    } // if newton

#ifdef COMBUST_SUPG_RHS
    // supg stabilisation
    assembler.template Vector<Velx>(enr_conv_c_sg, -denstauMu*res_old(0));
    assembler.template Vector<Vely>(enr_conv_c_sg, -denstauMu*res_old(1));
    assembler.template Vector<Velz>(enr_conv_c_sg, -denstauMu*res_old(2));
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
    assembler.template Matrix<Velx,Velx>(shpvel.dx, timefac_tau_C, shpvel.dx);
    assembler.template Matrix<Velx,Vely>(shpvel.dx, timefac_tau_C, shpvel.dy);
    assembler.template Matrix<Velx,Velz>(shpvel.dx, timefac_tau_C, shpvel.dz);

    assembler.template Matrix<Vely,Velx>(shpvel.dy, timefac_tau_C, shpvel.dx);
    assembler.template Matrix<Vely,Vely>(shpvel.dy, timefac_tau_C, shpvel.dy);
    assembler.template Matrix<Vely,Velz>(shpvel.dy, timefac_tau_C, shpvel.dz);

    assembler.template Matrix<Velz,Velx>(shpvel.dz, timefac_tau_C, shpvel.dx);
    assembler.template Matrix<Velz,Vely>(shpvel.dz, timefac_tau_C, shpvel.dy);
    assembler.template Matrix<Velz,Velz>(shpvel.dz, timefac_tau_C, shpvel.dz);

    assembler.template Vector<Velx>(shpvel.dx, -timefac_tau_C_divunp);
    assembler.template Vector<Vely>(shpvel.dy, -timefac_tau_C_divunp);
    assembler.template Vector<Velz>(shpvel.dz, -timefac_tau_C_divunp);
#endif
  } // endif cstab
}


template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          size_t NUMDOF,
          size_t shpVecSizeVel,
          size_t shpVecSizePres>
void BuildStandardDomainIntegralsNormal(
    COMBUST::LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>& assembler,
    const XFEM::ApproxFuncNormalVector<2,shpVecSizeVel>& shp, // shape functions and derivatives
    const XFEM::ApproxFunc<1,shpVecSizePres>& shp_pres,    // shape functions and derivatives
    const double& fac,                                     // factor for integration in space
    const double& timefac,                                 // factor for integration in time
    const double& timefacfac,                              // overall integration factor
    const double& densaf,                                  // density at n+1 (OST), alpha_f (GenAlpha)
    const double& dynvisc,                                 // dynamic viscosity at n+1
    const LINALG::Matrix<3,1>& velnp,                      // velocity at Gaussian point (GP) at n+1
    const double&              pres,                       // pressure at GP at n+1
    const LINALG::Matrix<3,1>& accam,                      // acceleration at n+alphaM
    const LINALG::Matrix<3,1>& gradp,                      // pressure gradient at GP at n+1
    const LINALG::Matrix<3,3>& vderxy,                     // velocity derivative at GP at n+1
    const LINALG::Matrix<3,1>& rhsint,                     // time RHS (at n) and body force (at n+1) at GP
    const LINALG::Matrix<3,1>& res_old,                    // residual at GP at n
    const LINALG::Matrix<3,1>& visc_old,                   // viscous term at GP at n
    const XFEM::ApproxFuncNormalVector<0,shpVecSizeVel>& enr_conv_c, // convective term at GP at n
    const XFEM::EnrViscs2N<shpVecSizeVel>& enr_viscs2,     // viscous term (second derivatves) at GP at n
    const bool instationary,                               // turn on/off instationary problem
    const bool newton,                                     // turn on/off linearization of convective term
    const bool pstab,                                      // turn on/off PSPG stabilization
    const bool supg,                                       // turn on/off SUPG stabilization
    const bool cstab,                                      // turn on/off BV   stabilization
    const double& tau_stab_Mu,                             // SUPG stabilization parameter
    const double& tau_stab_Mp,                             // PSPG stabilization parameter
    const double& tau_stab_C                               // BV   stabilization parameter
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
    assembler.template MatrixNN<Velx,Velx>(shp.velx.d0, densfac, shp.velx.d0);
    assembler.template MatrixNN<Vely,Vely>(shp.vely.d0, densfac, shp.vely.d0);
    assembler.template MatrixNN<Velz,Velz>(shp.velz.d0, densfac, shp.velz.d0);

    assembler.template VectorN<Velx>(shp.velx.d0, -densfac*velnp(0));
    assembler.template VectorN<Vely>(shp.vely.d0, -densfac*velnp(1));
    assembler.template VectorN<Velz>(shp.velz.d0, -densfac*velnp(2));
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
  assembler.template MatrixNN<Velx,Velx>(shp.velx.d0, denstimefacfac, enr_conv_c.velx.d0);
  assembler.template MatrixNN<Vely,Vely>(shp.vely.d0, denstimefacfac, enr_conv_c.vely.d0);
  assembler.template MatrixNN<Velz,Velz>(shp.velz.d0, denstimefacfac, enr_conv_c.velz.d0);

  assembler.template VectorN<Velx>(shp.velx.d0, -denstimefacfac*(velnp(0)*vderxy(0,0)
                                                                +velnp(1)*vderxy(0,1)
                                                                +velnp(2)*vderxy(0,2)));
  assembler.template VectorN<Vely>(shp.vely.d0, -denstimefacfac*(velnp(0)*vderxy(1,0)
                                                                +velnp(1)*vderxy(1,1)
                                                                +velnp(2)*vderxy(1,2)));
  assembler.template VectorN<Velz>(shp.velz.d0, -denstimefacfac*(velnp(0)*vderxy(2,0)
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
    assembler.template MatrixNN<Velx,Velx>(shp.velx.d0, denstimefacfac*vderxy(0,0), shp.velx.d0);
    assembler.template MatrixNN<Velx,Vely>(shp.velx.d0, denstimefacfac*vderxy(0,1), shp.vely.d0);
    assembler.template MatrixNN<Velx,Velz>(shp.velx.d0, denstimefacfac*vderxy(0,2), shp.velz.d0);
    assembler.template MatrixNN<Vely,Velx>(shp.vely.d0, denstimefacfac*vderxy(1,0), shp.velx.d0);
    assembler.template MatrixNN<Vely,Vely>(shp.vely.d0, denstimefacfac*vderxy(1,1), shp.vely.d0);
    assembler.template MatrixNN<Vely,Velz>(shp.vely.d0, denstimefacfac*vderxy(1,2), shp.velz.d0);
    assembler.template MatrixNN<Velz,Velx>(shp.velz.d0, denstimefacfac*vderxy(2,0), shp.velx.d0);
    assembler.template MatrixNN<Velz,Vely>(shp.velz.d0, denstimefacfac*vderxy(2,1), shp.vely.d0);
    assembler.template MatrixNN<Velz,Velz>(shp.velz.d0, denstimefacfac*vderxy(2,2), shp.velz.d0);
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
  assembler.template MatrixNN<Velx,Velx>(shp.velx.dx, 2.0*dynvisc*timefacfac, shp.velx.dx);
  assembler.template MatrixNN<Velx,Velx>(shp.velx.dy,     dynvisc*timefacfac, shp.velx.dy);
  assembler.template MatrixNN<Velx,Vely>(shp.velx.dy,     dynvisc*timefacfac, shp.vely.dx);
  assembler.template MatrixNN<Velx,Velx>(shp.velx.dz,     dynvisc*timefacfac, shp.velx.dz);
  assembler.template MatrixNN<Velx,Velz>(shp.velx.dz,     dynvisc*timefacfac, shp.velz.dx);

  assembler.template MatrixNN<Vely,Vely>(shp.vely.dx,     dynvisc*timefacfac, shp.vely.dx);
  assembler.template MatrixNN<Vely,Velx>(shp.vely.dx,     dynvisc*timefacfac, shp.velx.dy);
  assembler.template MatrixNN<Vely,Vely>(shp.vely.dy, 2.0*dynvisc*timefacfac, shp.vely.dy);
  assembler.template MatrixNN<Vely,Vely>(shp.vely.dz,     dynvisc*timefacfac, shp.vely.dz);
  assembler.template MatrixNN<Vely,Velz>(shp.vely.dz,     dynvisc*timefacfac, shp.velz.dy);

  assembler.template MatrixNN<Velz,Velz>(shp.velz.dx,     dynvisc*timefacfac, shp.velz.dx);
  assembler.template MatrixNN<Velz,Velx>(shp.velz.dx,     dynvisc*timefacfac, shp.velx.dz);
  assembler.template MatrixNN<Velz,Velz>(shp.velz.dy,     dynvisc*timefacfac, shp.velz.dy);
  assembler.template MatrixNN<Velz,Vely>(shp.velz.dy,     dynvisc*timefacfac, shp.vely.dz);
  assembler.template MatrixNN<Velz,Velz>(shp.velz.dz, 2.0*dynvisc*timefacfac, shp.velz.dz);

  assembler.template VectorN<Velx>(shp.velx.dx,     -dynvisc*timefacfac*(vderxy(0, 0) + vderxy(0, 0)));
  assembler.template VectorN<Velx>(shp.velx.dy,     -dynvisc*timefacfac*(vderxy(0, 1) + vderxy(1, 0)));
  assembler.template VectorN<Velx>(shp.velx.dz,     -dynvisc*timefacfac*(vderxy(0, 2) + vderxy(2, 0)));

  assembler.template VectorN<Vely>(shp.vely.dx,     -dynvisc*timefacfac*(vderxy(1, 0) + vderxy(0, 1)));
  assembler.template VectorN<Vely>(shp.vely.dy,     -dynvisc*timefacfac*(vderxy(1, 1) + vderxy(1, 1)));
  assembler.template VectorN<Vely>(shp.vely.dz,     -dynvisc*timefacfac*(vderxy(1, 2) + vderxy(2, 1)));

  assembler.template VectorN<Velz>(shp.velz.dx,     -dynvisc*timefacfac*(vderxy(2, 0) + vderxy(0, 2)));
  assembler.template VectorN<Velz>(shp.velz.dy,     -dynvisc*timefacfac*(vderxy(2, 1) + vderxy(1, 2)));
  assembler.template VectorN<Velz>(shp.velz.dz,     -dynvisc*timefacfac*(vderxy(2, 2) + vderxy(2, 2)));
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
  assembler.template MatrixNS<Velx,Pres,shpVecSizePres>(shp.velx.dx, -timefacfac, shp_pres.d0);
  assembler.template MatrixNS<Vely,Pres,shpVecSizePres>(shp.vely.dy, -timefacfac, shp_pres.d0);
  assembler.template MatrixNS<Velz,Pres,shpVecSizePres>(shp.velz.dz, -timefacfac, shp_pres.d0);

  assembler.template VectorN<Velx>(shp.velx.dx, timefacfac*pres);
  assembler.template VectorN<Vely>(shp.vely.dy, timefacfac*pres);
  assembler.template VectorN<Velz>(shp.velz.dz, timefacfac*pres);
#endif

#ifdef COMBUST_DOM_CONTINUITY
  //cout << "/!\\ warning === sign of continuity term has been changed to make formulation symmetric" << endl;
  // solenoidality term - continuity equation
  /*
                 /              \
                |                |
                | q , nabla o Du |
                |                |
                 \              /
   */
  assembler.template MatrixSN<Pres,Velx,shpVecSizePres>(shp_pres.d0, timefacfac, shp.velx.dx);
  assembler.template MatrixSN<Pres,Vely,shpVecSizePres>(shp_pres.d0, timefacfac, shp.vely.dy);
  assembler.template MatrixSN<Pres,Velz,shpVecSizePres>(shp_pres.d0, timefacfac, shp.velz.dz);

  const double trace_gamma = (vderxy(0,0) + vderxy(1,1) + vderxy(2,2));
  assembler.template Vector<Pres>(shp_pres.d0, -timefacfac*trace_gamma);
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
  assembler.template VectorN<Velx>(shp.velx.d0, fac*rhsint(0));
  assembler.template VectorN<Vely>(shp.vely.d0, fac*rhsint(1));
  assembler.template VectorN<Velz>(shp.velz.d0, fac*rhsint(2));
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
      assembler.template MatrixSN<Pres,Velx,shpVecSizePres>(shp_pres.dx, denstauMp, shp.velx.d0);
      assembler.template MatrixSN<Pres,Vely,shpVecSizePres>(shp_pres.dy, denstauMp, shp.vely.d0);
      assembler.template MatrixSN<Pres,Velz,shpVecSizePres>(shp_pres.dz, denstauMp, shp.velz.d0);
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
    assembler.template MatrixSN<Pres,Velx>(shp_pres.dx, denstimetauMp, enr_conv_c.velx.d0);
    assembler.template MatrixSN<Pres,Vely>(shp_pres.dy, denstimetauMp, enr_conv_c.vely.d0);
    assembler.template MatrixSN<Pres,Velz>(shp_pres.dz, denstimetauMp, enr_conv_c.velz.d0);

    if (newton)
    {
      /*  pressure stabilisation: convection, reactive part
                /                             \
               |           /          \   n+1  |
               | grad q , | Du o nabla | u     |
               |           \          /   (i)  |
                \                             /
       */
      assembler.template MatrixSN<Pres,Velx,shpVecSizePres>(shp_pres.dx, denstimetauMp*vderxy(0,0), shp.velx.d0);
      assembler.template MatrixSN<Pres,Velx,shpVecSizePres>(shp_pres.dy, denstimetauMp*vderxy(1,0), shp.velx.d0);
      assembler.template MatrixSN<Pres,Velx,shpVecSizePres>(shp_pres.dz, denstimetauMp*vderxy(2,0), shp.velx.d0);

      assembler.template MatrixSN<Pres,Vely,shpVecSizePres>(shp_pres.dx, denstimetauMp*vderxy(0,1), shp.vely.d0);
      assembler.template MatrixSN<Pres,Vely,shpVecSizePres>(shp_pres.dy, denstimetauMp*vderxy(1,1), shp.vely.d0);
      assembler.template MatrixSN<Pres,Vely,shpVecSizePres>(shp_pres.dz, denstimetauMp*vderxy(2,1), shp.vely.d0);

      assembler.template MatrixSN<Pres,Velz,shpVecSizePres>(shp_pres.dx, denstimetauMp*vderxy(0,2), shp.velz.d0);
      assembler.template MatrixSN<Pres,Velz,shpVecSizePres>(shp_pres.dy, denstimetauMp*vderxy(1,2), shp.velz.d0);
      assembler.template MatrixSN<Pres,Velz,shpVecSizePres>(shp_pres.dz, denstimetauMp*vderxy(2,2), shp.velz.d0);
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
    assembler.template MatrixSN<Pres,Velx>(shp_pres.dx, -2.0*dynvisc*timetauMp, enr_viscs2.xx);
    assembler.template MatrixSN<Pres,Vely>(shp_pres.dx, -2.0*dynvisc*timetauMp, enr_viscs2.xy);
    assembler.template MatrixSN<Pres,Velz>(shp_pres.dx, -2.0*dynvisc*timetauMp, enr_viscs2.xz);

    assembler.template MatrixSN<Pres,Velx>(shp_pres.dy, -2.0*dynvisc*timetauMp, enr_viscs2.yx);
    assembler.template MatrixSN<Pres,Vely>(shp_pres.dy, -2.0*dynvisc*timetauMp, enr_viscs2.yy);
    assembler.template MatrixSN<Pres,Velz>(shp_pres.dy, -2.0*dynvisc*timetauMp, enr_viscs2.yz);

    assembler.template MatrixSN<Pres,Velx>(shp_pres.dz, -2.0*dynvisc*timetauMp, enr_viscs2.zx);
    assembler.template MatrixSN<Pres,Vely>(shp_pres.dz, -2.0*dynvisc*timetauMp, enr_viscs2.zy);
    assembler.template MatrixSN<Pres,Velz>(shp_pres.dz, -2.0*dynvisc*timetauMp, enr_viscs2.zz);
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
    assembler.template Matrix<Pres,Pres>(shp_pres.dx, timetauMp, shp_pres.dx);
    assembler.template Matrix<Pres,Pres>(shp_pres.dy, timetauMp, shp_pres.dy);
    assembler.template Matrix<Pres,Pres>(shp_pres.dz, timetauMp, shp_pres.dz);
#endif

#ifdef COMBUST_PSPG_RHS
    // pressure stabilization
    assembler.template Vector<Pres>(shp_pres.dx, -tauMp*res_old(0));
    assembler.template Vector<Pres>(shp_pres.dy, -tauMp*res_old(1));
    assembler.template Vector<Pres>(shp_pres.dz, -tauMp*res_old(2));
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
      assembler.template MatrixNN<Velx,Velx>(enr_conv_c.velx.d0, densdenstauMu, shp.velx.d0);
      assembler.template MatrixNN<Vely,Vely>(enr_conv_c.vely.d0, densdenstauMu, shp.vely.d0);
      assembler.template MatrixNN<Velz,Velz>(enr_conv_c.velz.d0, densdenstauMu, shp.velz.d0);

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
        assembler.template MatrixNN<Velx,Velx>(shp.velx.dx, densdenstauMu*velnp(0), shp.velx.d0);
        assembler.template MatrixNN<Velx,Vely>(shp.velx.dy, densdenstauMu*velnp(0), shp.vely.d0);
        assembler.template MatrixNN<Velx,Velz>(shp.velx.dz, densdenstauMu*velnp(0), shp.velz.d0);
        assembler.template MatrixNN<Vely,Velx>(shp.vely.dx, densdenstauMu*velnp(1), shp.velx.d0);
        assembler.template MatrixNN<Vely,Vely>(shp.vely.dy, densdenstauMu*velnp(1), shp.vely.d0);
        assembler.template MatrixNN<Vely,Velz>(shp.vely.dz, densdenstauMu*velnp(1), shp.velz.d0);
        assembler.template MatrixNN<Velz,Velx>(shp.velz.dx, densdenstauMu*velnp(2), shp.velx.d0);
        assembler.template MatrixNN<Velz,Vely>(shp.velz.dy, densdenstauMu*velnp(2), shp.vely.d0);
        assembler.template MatrixNN<Velz,Velz>(shp.velz.dz, densdenstauMu*velnp(2), shp.velz.d0);
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
    assembler.template MatrixNN<Velx,Velx>(enr_conv_c.velx.d0, densdenstimetauMu, enr_conv_c.velx.d0);
    assembler.template MatrixNN<Vely,Vely>(enr_conv_c.vely.d0, densdenstimetauMu, enr_conv_c.vely.d0);
    assembler.template MatrixNN<Velz,Velz>(enr_conv_c.velz.d0, densdenstimetauMu, enr_conv_c.velz.d0);
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
    assembler.template MatrixNS<Velx,Pres>(enr_conv_c.velx.d0, denstimetauMu, shp_pres.dx);
    assembler.template MatrixNS<Vely,Pres>(enr_conv_c.vely.d0, denstimetauMu, shp_pres.dy);
    assembler.template MatrixNS<Velz,Pres>(enr_conv_c.velz.d0, denstimetauMu, shp_pres.dz);
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
    assembler.template MatrixNN<Velx,Velx>(enr_conv_c.velx.d0, -2.0*dynvisc*denstimetauMu, enr_viscs2.xx);
    assembler.template MatrixNN<Velx,Vely>(enr_conv_c.velx.d0, -2.0*dynvisc*denstimetauMu, enr_viscs2.xy);
    assembler.template MatrixNN<Velx,Velz>(enr_conv_c.velx.d0, -2.0*dynvisc*denstimetauMu, enr_viscs2.xz);

    assembler.template MatrixNN<Vely,Velx>(enr_conv_c.vely.d0, -2.0*dynvisc*denstimetauMu, enr_viscs2.yx);
    assembler.template MatrixNN<Vely,Vely>(enr_conv_c.vely.d0, -2.0*dynvisc*denstimetauMu, enr_viscs2.yy);
    assembler.template MatrixNN<Vely,Velz>(enr_conv_c.vely.d0, -2.0*dynvisc*denstimetauMu, enr_viscs2.yz);

    assembler.template MatrixNN<Velz,Velx>(enr_conv_c.velz.d0, -2.0*dynvisc*denstimetauMu, enr_viscs2.zx);
    assembler.template MatrixNN<Velz,Vely>(enr_conv_c.velz.d0, -2.0*dynvisc*denstimetauMu, enr_viscs2.zy);
    assembler.template MatrixNN<Velz,Velz>(enr_conv_c.velz.d0, -2.0*dynvisc*denstimetauMu, enr_viscs2.zz);
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
      assembler.template MatrixNN<Velx,Velx>(enr_conv_c.velx.d0, densdenstimetauMu*vderxy(0,0), shp.velx.d0);
      assembler.template MatrixNN<Velx,Vely>(enr_conv_c.velx.d0, densdenstimetauMu*vderxy(0,1), shp.vely.d0);
      assembler.template MatrixNN<Velx,Velz>(enr_conv_c.velx.d0, densdenstimetauMu*vderxy(0,2), shp.velz.d0);

      assembler.template MatrixNN<Vely,Velx>(enr_conv_c.vely.d0, densdenstimetauMu*vderxy(1,0), shp.velx.d0);
      assembler.template MatrixNN<Vely,Vely>(enr_conv_c.vely.d0, densdenstimetauMu*vderxy(1,1), shp.vely.d0);
      assembler.template MatrixNN<Vely,Velz>(enr_conv_c.vely.d0, densdenstimetauMu*vderxy(1,2), shp.velz.d0);

      assembler.template MatrixNN<Velz,Velx>(enr_conv_c.velz.d0, densdenstimetauMu*vderxy(2,0), shp.velx.d0);
      assembler.template MatrixNN<Velz,Vely>(enr_conv_c.velz.d0, densdenstimetauMu*vderxy(2,1), shp.vely.d0);
      assembler.template MatrixNN<Velz,Velz>(enr_conv_c.velz.d0, densdenstimetauMu*vderxy(2,2), shp.velz.d0);

      /*
                   /                                           \
                  |    / n+1        \   n+1    /          \     |
                  |   | u    o nabla | u    , | Du o nabla | v  |
                  |    \ (i)        /   (i)    \          /     |
                   \                                           /
       */
      const double con0 = densdenstimetauMu*(velnp(0)*vderxy(0,0) + velnp(1)*vderxy(0,1) + velnp(2)*vderxy(0,2));
      assembler.template MatrixNN<Velx,Velx>(shp.velx.dx, con0, shp.velx.d0);
      assembler.template MatrixNN<Velx,Vely>(shp.velx.dy, con0, shp.vely.d0);
      assembler.template MatrixNN<Velx,Velz>(shp.velx.dz, con0, shp.velz.d0);

      const double con1 = densdenstimetauMu*(velnp(0)*vderxy(1,0) + velnp(1)*vderxy(1,1) + velnp(2)*vderxy(1,2));
      assembler.template MatrixNN<Vely,Velx>(shp.vely.dx, con1, shp.velx.d0);
      assembler.template MatrixNN<Vely,Vely>(shp.vely.dy, con1, shp.vely.d0);
      assembler.template MatrixNN<Vely,Velz>(shp.vely.dz, con1, shp.velz.d0);

      const double con2 = densdenstimetauMu*(velnp(0)*vderxy(2,0) + velnp(1)*vderxy(2,1) + velnp(2)*vderxy(2,2));
      assembler.template MatrixNN<Velz,Velx>(shp.velz.dx, con2, shp.velx.d0);
      assembler.template MatrixNN<Velz,Vely>(shp.velz.dy, con2, shp.vely.d0);
      assembler.template MatrixNN<Velz,Velz>(shp.velz.dz, con2, shp.velz.d0);
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
      assembler.template MatrixNN<Velx,Velx>(shp.velx.dx, denstimetauMu*gradp(0), shp.velx.d0);
      assembler.template MatrixNN<Velx,Vely>(shp.velx.dy, denstimetauMu*gradp(0), shp.vely.d0);
      assembler.template MatrixNN<Velx,Velz>(shp.velx.dz, denstimetauMu*gradp(0), shp.velz.d0);

      assembler.template MatrixNN<Vely,Velx>(shp.vely.dx, denstimetauMu*gradp(1), shp.velx.d0);
      assembler.template MatrixNN<Vely,Vely>(shp.vely.dy, denstimetauMu*gradp(1), shp.vely.d0);
      assembler.template MatrixNN<Vely,Velz>(shp.vely.dz, denstimetauMu*gradp(1), shp.velz.d0);

      assembler.template MatrixNN<Velz,Velx>(shp.velz.dx, denstimetauMu*gradp(2), shp.velx.d0);
      assembler.template MatrixNN<Velz,Vely>(shp.velz.dy, denstimetauMu*gradp(2), shp.vely.d0);
      assembler.template MatrixNN<Velz,Velz>(shp.velz.dz, denstimetauMu*gradp(2), shp.velz.d0);
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
      assembler.template MatrixNN<Velx,Velx>(shp.velx.dx, -2.0*dynvisc*denstimetauMu*visc_old(0), shp.velx.d0);
      assembler.template MatrixNN<Velx,Vely>(shp.velx.dy, -2.0*dynvisc*denstimetauMu*visc_old(0), shp.vely.d0);
      assembler.template MatrixNN<Velx,Velz>(shp.velx.dz, -2.0*dynvisc*denstimetauMu*visc_old(0), shp.velz.d0);

      assembler.template MatrixNN<Vely,Velx>(shp.vely.dx, -2.0*dynvisc*denstimetauMu*visc_old(1), shp.velx.d0);
      assembler.template MatrixNN<Vely,Vely>(shp.vely.dy, -2.0*dynvisc*denstimetauMu*visc_old(1), shp.vely.d0);
      assembler.template MatrixNN<Vely,Velz>(shp.vely.dz, -2.0*dynvisc*denstimetauMu*visc_old(1), shp.velz.d0);

      assembler.template MatrixNN<Velz,Velx>(shp.velz.dx, -2.0*dynvisc*denstimetauMu*visc_old(2), shp.velx.d0);
      assembler.template MatrixNN<Velz,Vely>(shp.velz.dy, -2.0*dynvisc*denstimetauMu*visc_old(2), shp.vely.d0);
      assembler.template MatrixNN<Velz,Velz>(shp.velz.dz, -2.0*dynvisc*denstimetauMu*visc_old(2), shp.velz.d0);
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
      assembler.template MatrixNN<Velx,Velx>(shp.velx.dx, -denstauMu*rhsint(0), shp.velx.d0);
      assembler.template MatrixNN<Velx,Vely>(shp.velx.dy, -denstauMu*rhsint(0), shp.vely.d0);
      assembler.template MatrixNN<Velx,Velz>(shp.velx.dz, -denstauMu*rhsint(0), shp.velz.d0);

      assembler.template MatrixNN<Vely,Velx>(shp.vely.dx, -denstauMu*rhsint(1), shp.velx.d0);
      assembler.template MatrixNN<Vely,Vely>(shp.vely.dy, -denstauMu*rhsint(1), shp.vely.d0);
      assembler.template MatrixNN<Vely,Velz>(shp.vely.dz, -denstauMu*rhsint(1), shp.velz.d0);

      assembler.template MatrixNN<Velz,Velx>(shp.velz.dx, -denstauMu*rhsint(2), shp.velx.d0);
      assembler.template MatrixNN<Velz,Vely>(shp.velz.dy, -denstauMu*rhsint(2), shp.vely.d0);
      assembler.template MatrixNN<Velz,Velz>(shp.velz.dz, -denstauMu*rhsint(2), shp.velz.d0);
#endif
    } // if newton

#ifdef COMBUST_SUPG_RHS
    // supg stabilisation
    assembler.template VectorN<Velx>(enr_conv_c.velx.d0, -denstauMu*res_old(0));
    assembler.template VectorN<Vely>(enr_conv_c.vely.d0, -denstauMu*res_old(1));
    assembler.template VectorN<Velz>(enr_conv_c.velz.d0, -denstauMu*res_old(2));
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
    assembler.template MatrixNN<Velx,Velx>(shp.velx.dx, timefac_tau_C, shp.velx.dx);
    assembler.template MatrixNN<Velx,Vely>(shp.velx.dx, timefac_tau_C, shp.vely.dy);
    assembler.template MatrixNN<Velx,Velz>(shp.velx.dx, timefac_tau_C, shp.velz.dz);

    assembler.template MatrixNN<Vely,Velx>(shp.vely.dy, timefac_tau_C, shp.velx.dx);
    assembler.template MatrixNN<Vely,Vely>(shp.vely.dy, timefac_tau_C, shp.vely.dy);
    assembler.template MatrixNN<Vely,Velz>(shp.vely.dy, timefac_tau_C, shp.velz.dz);

    assembler.template MatrixNN<Velz,Velx>(shp.velz.dz, timefac_tau_C, shp.velx.dx);
    assembler.template MatrixNN<Velz,Vely>(shp.velz.dz, timefac_tau_C, shp.vely.dy);
    assembler.template MatrixNN<Velz,Velz>(shp.velz.dz, timefac_tau_C, shp.velz.dz);

    assembler.template VectorN<Velx>(shp.velx.dx, -timefac_tau_C_divunp);
    assembler.template VectorN<Vely>(shp.vely.dy, -timefac_tau_C_divunp);
    assembler.template VectorN<Velz>(shp.velz.dz, -timefac_tau_C_divunp);
#endif
  } // endif cstab
}

} // namespace COMBUST

#endif
#endif
#endif
