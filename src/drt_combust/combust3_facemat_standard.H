/*----------------------------------------------------------------------*/
/*!
\file combust3_facemat_standard.H

\brief patch matrix standard FEM terms
       for premixed combustion and two-phase flow problems / 3D combustion element

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_FACEMAT_STANDARD_H
#define COMBUST3_FACEMAT_STANDARD_H


#include "combust3_local_assembler.H"
#include "../drt_inpar/inpar_fluid.H"


namespace COMBUST
{

//template <DRT::Element::DiscretizationType DISTYPE,
//          XFEM::AssemblyType ASSTYPE,
//          size_t NUMDOF,
//          size_t shpVecSize>
template <DRT::Element::DiscretizationType M_DISTYPE,
          DRT::Element::DiscretizationType S_DISTYPE,
          XFEM::AssemblyType M_ASSTYPE,
          XFEM::AssemblyType S_ASSTYPE,
          size_t NUMDOF,
          size_t m_shpVecSize,
          size_t s_shpVecSize>
void BuildStandardFaceIntegrals(
      COMBUST::LocalPatchAssembler<M_DISTYPE, S_DISTYPE, M_ASSTYPE, S_ASSTYPE, NUMDOF>& assembler,
    const XFEM::ApproxFunc<2,m_shpVecSize>& m_shpvel,    // shape functions and their derivatives (velocity) (master side)
    const XFEM::ApproxFunc<2,s_shpVecSize>& s_shpvel,    // shape functions and their derivatives (velocity) (slave side)
    const XFEM::ApproxFunc<2,m_shpVecSize>& m_shppres,   // shape functions and their derivatives (pressure) (master side)
    const XFEM::ApproxFunc<2,s_shpVecSize>& s_shppres,   // shape functions and their derivatives (pressure) (slave side)
    const double& fac,                                 // factor for integration in space
    const double& timefac,                             // factor for integration in time
    const double& rhsfac,                              // overall integration factor for the rhs
    const bool genalpha,                               // switch to generalized alpha terms
    const bool instationary,                           // turn on/off instationary problem
    const double& densaf,                              // density at n+1 (OST), alpha_f (GenAlpha)
    const double& dynvisc,                             // dynamic viscosity at n+1
    const LINALG::Matrix<3,1>& velnp,                  // velocity at Gaussian point (GP) at n+1
    const LINALG::Matrix<3,3>& m_vderxy,               // velocity derivative at GP at n+1 (master side)
    const LINALG::Matrix<3,3>& s_vderxy,               // velocity derivative at GP at n+1 (slave side)
    const double&              pres,                   // pressure at GP at n+1
    const LINALG::Matrix<3,1>& m_gradp,                // pressure gradient at GP at n+1 (master side)
    const LINALG::Matrix<3,1>& s_gradp,                // pressure gradient at GP at n+1 (slave side)
    const LINALG::Matrix<3,1>& grad_pres_jump,         // pressure jump at interface at GP
    const double& vdiv_jump,                           // divergence jump at interface at GP
    const double& normal_vel,                          // projection of velocity in normal direction of face
    LINALG::Matrix<3,3>& grad_vel_jump,            // gradient jump of velocity at interface at GP
    const INPAR::FLUID::EOS_Pres        pres_stab,       ///< face term to evaluate: pressure stab
    const INPAR::FLUID::EOS_Conv_Stream conv_stream_stab,///< face term to evaluate: conv stab stream
    const INPAR::FLUID::EOS_Conv_Cross  conv_cross_stab, ///< face term to evaluate: conv stab cross
    const INPAR::FLUID::EOS_Div         conti_stab,      ///< face term to evaluate: conti stab
    const double& tau_conv,                            // convection-stabilization parameter
    const double& tau_div,                             // continuity-stabilization parameter
    const double& tau_p                                // pressure-stabilization parameter
)
{
  // TODO: reicht es ja nein abzufragen oder muss ich hier zwischen den optionen unterscheiden
  // TODO: Zeitintegration

  if (conv_stream_stab != INPAR::FLUID::EOS_CONV_STREAM_none)
  {
    /*
    //
    //
    //             /                                              \
    //            |   lin                                          |
    //  + tau_u * | | P  ( u )*n | * |[ grad Du ]| : |[ grad v ]|  |
    //            |                                                |
    //             \                                              / surface
    //
    */

//    std::cout << "Eval conv" << std::endl;
    const double conv_fac = tau_conv * normal_vel * fac * timefac;
//    std::cout << normal_vel << std::endl;
//    std::cout << conv_fac << std::endl;

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac, m_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac, m_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac, m_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac, m_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac, m_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac, m_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac, m_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac, m_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac, m_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dx, -conv_fac, s_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dy, -conv_fac, s_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dz, -conv_fac, s_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dx, -conv_fac, s_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dy, -conv_fac, s_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dz, -conv_fac, s_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dx, -conv_fac, s_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dy, -conv_fac, s_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dz, -conv_fac, s_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dx, -conv_fac, m_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dy, -conv_fac, m_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dz, -conv_fac, m_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dx, -conv_fac, m_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dy, -conv_fac, m_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dz, -conv_fac, m_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dx, -conv_fac, m_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dy, -conv_fac, m_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dz, -conv_fac, m_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dx, conv_fac, s_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dy, conv_fac, s_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dz, conv_fac, s_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dx, conv_fac, s_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dy, conv_fac, s_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dz, conv_fac, s_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dx, conv_fac, s_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dy, conv_fac, s_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dz, conv_fac, s_shpvel.dz);

//    assembler.PrintMat();

    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac * grad_vel_jump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac * grad_vel_jump(0,1));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac * grad_vel_jump(0,2));

    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac * grad_vel_jump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac * grad_vel_jump(1,1));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac * grad_vel_jump(1,2));

    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac * grad_vel_jump(2,0));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac * grad_vel_jump(2,1));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac * grad_vel_jump(2,2));

    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dx, -conv_fac * grad_vel_jump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dy, -conv_fac * grad_vel_jump(0,1));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dz, -conv_fac * grad_vel_jump(0,2));

    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dx, -conv_fac * grad_vel_jump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dy, -conv_fac * grad_vel_jump(1,1));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dz, -conv_fac * grad_vel_jump(1,2));

    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dx, -conv_fac * grad_vel_jump(2,0));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dy, -conv_fac * grad_vel_jump(2,1));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dz, -conv_fac * grad_vel_jump(2,2));

    // TODO: verguengung mit grad u nicht grad u ^ T: ist wohl falsch
    //    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac, m_shpvel.dx);
    //    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac, m_shpvel.dx);
    //    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac, m_shpvel.dx);
    //
    //    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac, m_shpvel.dy);
    //    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac, m_shpvel.dy);
    //    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac, m_shpvel.dy);
    //
    //    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac, m_shpvel.dz);
    //    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac, m_shpvel.dz);
    //    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac, m_shpvel.dz);
    //
    //    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dx, -conv_fac, s_shpvel.dx);
    //    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dy, -conv_fac, s_shpvel.dx);
    //    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dz, -conv_fac, s_shpvel.dx);
    //
    //    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dx, -conv_fac, s_shpvel.dy);
    //    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dy, -conv_fac, s_shpvel.dy);
    //    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dz, -conv_fac, s_shpvel.dy);
    //
    //    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dx, -conv_fac, s_shpvel.dz);
    //    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dy, -conv_fac, s_shpvel.dz);
    //    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dz, -conv_fac, s_shpvel.dz);
    //
    //    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dx, -conv_fac, m_shpvel.dx);
    //    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dy, -conv_fac, m_shpvel.dx);
    //    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dz, -conv_fac, m_shpvel.dx);
    //
    //    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dx, -conv_fac, m_shpvel.dy);
    //    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dy, -conv_fac, m_shpvel.dy);
    //    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dz, -conv_fac, m_shpvel.dy);
    //
    //    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dx, -conv_fac, m_shpvel.dz);
    //    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dy, -conv_fac, m_shpvel.dz);
    //    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dz, -conv_fac, m_shpvel.dz);
    //
    //    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dx, conv_fac, s_shpvel.dx);
    //    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dy, conv_fac, s_shpvel.dx);
    //    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dz, conv_fac, s_shpvel.dx);
    //
    //    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dx, conv_fac, s_shpvel.dy);
    //    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dy, conv_fac, s_shpvel.dy);
    //    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dz, conv_fac, s_shpvel.dy);
    //
    //    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dx, conv_fac, s_shpvel.dz);
    //    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dy, conv_fac, s_shpvel.dz);
    //    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dz, conv_fac, s_shpvel.dz);

    //    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac * grad_vel_jump(0,0));
    //    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac * grad_vel_jump(1,0));
    //    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac * grad_vel_jump(2,0));
    //
    //    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac * grad_vel_jump(0,1));
    //    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac * grad_vel_jump(1,1));
    //    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac * grad_vel_jump(2,1));
    //
    //    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac * grad_vel_jump(0,2));
    //    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac * grad_vel_jump(1,2));
    //    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac * grad_vel_jump(2,2));

    //    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dx, -conv_fac * grad_vel_jump(0,0));
    //    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dy, -conv_fac * grad_vel_jump(1,0));
    //    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dz, -conv_fac * grad_vel_jump(2,0));
    //
    //    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dx, -conv_fac * grad_vel_jump(0,1));
    //    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dy, -conv_fac * grad_vel_jump(1,1));
    //    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dz, -conv_fac * grad_vel_jump(2,1));
    //
    //    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dx, -conv_fac * grad_vel_jump(0,2));
    //    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dy, -conv_fac * grad_vel_jump(1,2));
    //    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dz, -conv_fac * grad_vel_jump(2,2));


//    assembler.PrintRhs();
  }

  if (pres_stab != INPAR::FLUID::EOS_PRES_none)
  {
    /*
    //
    //             /                             \
    //            |                               |
    //  + tau_p * |  |[ grad q ]| , |[ grad p ]|  |
    //            |                               |
    //             \                             / surface
    //
     */


//      std::cout << "Eval pres" << std::endl;
    const double pres_fac = tau_p * fac * timefac;
//    std::cout << pres_fac << std::endl;

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shppres.dx, pres_fac, m_shppres.dx);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shppres.dy, pres_fac, m_shppres.dy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shppres.dz, pres_fac, m_shppres.dz);

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shppres.dx, -pres_fac, s_shppres.dx);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shppres.dy, -pres_fac, s_shppres.dy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shppres.dz, -pres_fac, s_shppres.dz);

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shppres.dx, -pres_fac, m_shppres.dx);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shppres.dy, -pres_fac, m_shppres.dy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shppres.dz, -pres_fac, m_shppres.dz);

//    std::cout << "slave-salve" << std::endl;
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shppres.dx, pres_fac, s_shppres.dx);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shppres.dy, pres_fac, s_shppres.dy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shppres.dz, pres_fac, s_shppres.dz);

//    assembler.PrintMat();

    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::master>(m_shppres.dx, pres_fac*grad_pres_jump(0));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::master>(m_shppres.dy, pres_fac*grad_pres_jump(1));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::master>(m_shppres.dz, pres_fac*grad_pres_jump(2));

    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::slave>(s_shppres.dx, -pres_fac*grad_pres_jump(0));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::slave>(s_shppres.dy, -pres_fac*grad_pres_jump(1));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::slave>(s_shppres.dz, -pres_fac*grad_pres_jump(2));
//    assembler.PrintRhs();
  }

  if (conti_stab != INPAR::FLUID::EOS_DIV_none)
  {
    /*
    //
    //
    //               /                               \
    //              |                                 |
    //  + tau_div * |   |[ div(u) ]| ,  |[ div(v) ]|  |
    //              |                                 |
    //               \                               / surface
    //
    */

//      std::cout << "Eval div" << std::endl;
    const double div_fac = tau_div * fac * timefac;
//    std::cout << div_fac << std::endl;

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dx, div_fac, m_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dx, div_fac, m_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dx, div_fac, m_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dy, div_fac, m_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dy, div_fac, m_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dy, div_fac, m_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dz, div_fac, m_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dz, div_fac, m_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dz, div_fac, m_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dx, -div_fac, s_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dx, -div_fac, s_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dx, -div_fac, s_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dy, -div_fac, s_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dy, -div_fac, s_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dy, -div_fac, s_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dz, -div_fac, s_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dz, -div_fac, s_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dz, -div_fac, s_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dx, -div_fac, m_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dx, -div_fac, m_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dx, -div_fac, m_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dy, -div_fac, m_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dy, -div_fac, m_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dy, -div_fac, m_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dz, -div_fac, m_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dz, -div_fac, m_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dz, -div_fac, m_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dx, div_fac, s_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dx, div_fac, s_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dx, div_fac, s_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dy, div_fac, s_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dy, div_fac, s_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dy, div_fac, s_shpvel.dz);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dz, div_fac, s_shpvel.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dz, div_fac, s_shpvel.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dz, div_fac, s_shpvel.dz);

//    assembler.PrintMat();

    const double div_fac_rhs = div_fac * vdiv_jump;

    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dx, div_fac_rhs);
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dy, div_fac_rhs);
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dz, div_fac_rhs);

    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dx, -div_fac_rhs);
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dy, -div_fac_rhs);
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dz, -div_fac_rhs);
//    assembler.PrintRhs();
  }


  return;
}

}// end namespace COMBUST
#endif
