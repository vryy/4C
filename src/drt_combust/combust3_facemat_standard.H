/*----------------------------------------------------------------------*/
/*!
\file combust3_facemat_standard.H

\brief patch matrix standard FEM terms
       for premixed combustion and two-phase flow problems / 3D combustion element

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
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
    const double fac,                                 // factor for integration in space
    const double timefac,                             // factor for integration in time
    const double rhstimefac,                          // time integration factor for the rhs
    const bool genalpha,                               // switch to generalized alpha terms
    const bool instationary,                           // turn on/off instationary problem
    const double densaf,                              // density at n+1 (OST), alpha_f (GenAlpha)
    const double dynvisc,                             // dynamic viscosity at n+1
    const LINALG::Matrix<3,1>& velnp,                  // velocity at Gaussian point (GP) at n+1
    const LINALG::Matrix<3,3>& m_vderxy,               // velocity derivative at GP at n+1 (master side)
    const LINALG::Matrix<3,3>& s_vderxy,               // velocity derivative at GP at n+1 (slave side)
    const double              pres,                   // pressure at GP at n+1
    const LINALG::Matrix<3,1>& m_gradp,                // pressure gradient at GP at n+1 (master side)
    const LINALG::Matrix<3,1>& s_gradp,                // pressure gradient at GP at n+1 (slave side)
    const LINALG::Matrix<3,1>& grad_pres_jump,         // pressure jump at interface at GP
    const double vdiv_jump,                           // divergence jump at interface at GP
    LINALG::Matrix<3,3>& grad_vel_jump,            // gradient jump of velocity at interface at GP
    LINALG::Matrix<3,6>& der2_vel_jump,             // 2nd derivative jump of velocity at interface at GP
    LINALG::Matrix<6,1>& der2_pres_jump,             // 2nd derivative jump of pressure at interface at GP
    const INPAR::FLUID::EOS_Pres        pres_stab,       /// face term to evaluate: pressure stab
    const INPAR::FLUID::EOS_Conv_Stream conv_stream_stab,/// face term to evaluate: conv stab stream
    const INPAR::FLUID::EOS_Conv_Cross  conv_cross_stab, /// face term to evaluate: conv stab cross
    const INPAR::FLUID::EOS_Div         conti_stab,      /// face term to evaluate: conti stab
    const double tau_conv,                            // convection-stabilization parameter
    const double tau_div,                             // continuity-stabilization parameter
    const double tau_p,                               // pressure-stabilization parameter
    const bool add_ghost_penalties,                   // face term to evaluate: ghost penalty
    const double tau_ghost_first_order,               // ghost penalty stabilization parameter: first order terms
    const double tau_ghost_second_order,              // ghost penalty stabilization parameter: second order terms
    const double tau_pres_second_order,               // ghost penalty stabilization parameter: second order term pressure
    const bool higher_order_ele                       // higher order ghost penalties
)
{
  if (conv_stream_stab == INPAR::FLUID::EOS_CONV_STREAM_std_eos or
      conti_stab == INPAR::FLUID::EOS_DIV_vel_jump_std_eos or
      add_ghost_penalties)
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

    double tau_all = 0.0;
    if (not add_ghost_penalties)
    {
      if (conv_stream_stab == INPAR::FLUID::EOS_CONV_STREAM_std_eos and
          conti_stab == INPAR::FLUID::EOS_DIV_vel_jump_std_eos)
        tau_all = tau_conv + tau_div;
      else if (conv_stream_stab != INPAR::FLUID::EOS_CONV_STREAM_std_eos and
               conti_stab == INPAR::FLUID::EOS_DIV_vel_jump_std_eos)
        tau_all = tau_div;
      else if (conv_stream_stab == INPAR::FLUID::EOS_CONV_STREAM_std_eos and
               conti_stab != INPAR::FLUID::EOS_DIV_vel_jump_std_eos)
        tau_all = tau_conv;
    }

    if (add_ghost_penalties)
    {
      tau_all = tau_conv + tau_ghost_first_order;
      if (conti_stab == INPAR::FLUID::EOS_DIV_vel_jump_std_eos)
        tau_all += tau_div;
    }

    const double conv_fac = tau_all * fac * timefac;
    const double conv_fac_rhs = tau_all * fac * rhstimefac;

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

    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac_rhs * grad_vel_jump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac_rhs * grad_vel_jump(0,1));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac_rhs * grad_vel_jump(0,2));

    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac_rhs * grad_vel_jump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac_rhs * grad_vel_jump(1,1));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac_rhs * grad_vel_jump(1,2));

    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dx, conv_fac_rhs * grad_vel_jump(2,0));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dy, conv_fac_rhs * grad_vel_jump(2,1));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dz, conv_fac_rhs * grad_vel_jump(2,2));

    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dx, -conv_fac_rhs * grad_vel_jump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dy, -conv_fac_rhs * grad_vel_jump(0,1));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dz, -conv_fac_rhs * grad_vel_jump(0,2));

    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dx, -conv_fac_rhs * grad_vel_jump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dy, -conv_fac_rhs * grad_vel_jump(1,1));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dz, -conv_fac_rhs * grad_vel_jump(1,2));

    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dx, -conv_fac_rhs * grad_vel_jump(2,0));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dy, -conv_fac_rhs * grad_vel_jump(2,1));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dz, -conv_fac_rhs * grad_vel_jump(2,2));
  }

  if (pres_stab == INPAR::FLUID::EOS_PRES_std_eos or
      add_ghost_penalties)
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

    const double pres_fac = tau_p * fac * timefac;
    const double pres_fac_rhs = tau_p * fac * rhstimefac;

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shppres.dx, pres_fac, m_shppres.dx);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shppres.dy, pres_fac, m_shppres.dy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shppres.dz, pres_fac, m_shppres.dz);

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shppres.dx, -pres_fac, s_shppres.dx);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shppres.dy, -pres_fac, s_shppres.dy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shppres.dz, -pres_fac, s_shppres.dz);

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shppres.dx, -pres_fac, m_shppres.dx);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shppres.dy, -pres_fac, m_shppres.dy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shppres.dz, -pres_fac, m_shppres.dz);

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shppres.dx, pres_fac, s_shppres.dx);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shppres.dy, pres_fac, s_shppres.dy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shppres.dz, pres_fac, s_shppres.dz);

    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::master>(m_shppres.dx, pres_fac_rhs*grad_pres_jump(0));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::master>(m_shppres.dy, pres_fac_rhs*grad_pres_jump(1));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::master>(m_shppres.dz, pres_fac_rhs*grad_pres_jump(2));

    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::slave>(s_shppres.dx, -pres_fac_rhs*grad_pres_jump(0));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::slave>(s_shppres.dy, -pres_fac_rhs*grad_pres_jump(1));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::slave>(s_shppres.dz, -pres_fac_rhs*grad_pres_jump(2));
  }

  if (conti_stab == INPAR::FLUID::EOS_DIV_div_jump_std_eos)
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

    const double div_fac = tau_div * fac * timefac;
    const double div_fac_rhs = tau_div * fac * rhstimefac * vdiv_jump;

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

    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dx, div_fac_rhs);
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dy, div_fac_rhs);
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dz, div_fac_rhs);

    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dx, -div_fac_rhs);
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dy, -div_fac_rhs);
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dz, -div_fac_rhs);
  }

  // higher order ghost penalties for hex elements
  if (add_ghost_penalties and higher_order_ele)
  {
    const double vel_fac = tau_ghost_second_order * fac * timefac;
    const double vel_fac_rhs = tau_ghost_second_order * fac * rhstimefac;

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dxdx, vel_fac, m_shpvel.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dydy, vel_fac, m_shpvel.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dzdz, vel_fac, m_shpvel.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dxdy, vel_fac, m_shpvel.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dxdz, vel_fac, m_shpvel.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dydz, vel_fac, m_shpvel.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dxdx, vel_fac, m_shpvel.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dydy, vel_fac, m_shpvel.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dzdz, vel_fac, m_shpvel.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dxdy, vel_fac, m_shpvel.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dxdz, vel_fac, m_shpvel.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dydz, vel_fac, m_shpvel.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dxdx, vel_fac, m_shpvel.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dydy, vel_fac, m_shpvel.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dzdz, vel_fac, m_shpvel.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dxdy, vel_fac, m_shpvel.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dxdz, vel_fac, m_shpvel.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shpvel.dydz, vel_fac, m_shpvel.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dxdx, -vel_fac, s_shpvel.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dydy, -vel_fac, s_shpvel.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dzdz, -vel_fac, s_shpvel.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dxdy, -vel_fac, s_shpvel.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dxdz, -vel_fac, s_shpvel.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dydz, -vel_fac, s_shpvel.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dxdx, -vel_fac, s_shpvel.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dydy, -vel_fac, s_shpvel.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dzdz, -vel_fac, s_shpvel.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dxdy, -vel_fac, s_shpvel.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dxdz, -vel_fac, s_shpvel.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dydz, -vel_fac, s_shpvel.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dxdx, -vel_fac, s_shpvel.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dydy, -vel_fac, s_shpvel.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dzdz, -vel_fac, s_shpvel.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dxdy, -vel_fac, s_shpvel.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dxdz, -vel_fac, s_shpvel.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shpvel.dydz, -vel_fac, s_shpvel.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dxdx, -vel_fac, m_shpvel.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dydy, -vel_fac, m_shpvel.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dzdz, -vel_fac, m_shpvel.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dxdy, -vel_fac, m_shpvel.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dxdz, -vel_fac, m_shpvel.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dydz, -vel_fac, m_shpvel.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dxdx, -vel_fac, m_shpvel.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dydy, -vel_fac, m_shpvel.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dzdz, -vel_fac, m_shpvel.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dxdy, -vel_fac, m_shpvel.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dxdz, -vel_fac, m_shpvel.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dydz, -vel_fac, m_shpvel.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dxdx, -vel_fac, m_shpvel.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dydy, -vel_fac, m_shpvel.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dzdz, -vel_fac, m_shpvel.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dxdy, -vel_fac, m_shpvel.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dxdz, -vel_fac, m_shpvel.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shpvel.dydz, -vel_fac, m_shpvel.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dxdx, vel_fac, s_shpvel.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dydy, vel_fac, s_shpvel.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dzdz, vel_fac, s_shpvel.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dxdy, vel_fac, s_shpvel.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dxdz, vel_fac, s_shpvel.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dydz, vel_fac, s_shpvel.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dxdx, vel_fac, s_shpvel.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dydy, vel_fac, s_shpvel.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dzdz, vel_fac, s_shpvel.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dxdy, vel_fac, s_shpvel.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dxdz, vel_fac, s_shpvel.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dydz, vel_fac, s_shpvel.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dxdx, vel_fac, s_shpvel.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dydy, vel_fac, s_shpvel.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dzdz, vel_fac, s_shpvel.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dxdy, vel_fac, s_shpvel.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dxdz, vel_fac, s_shpvel.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shpvel.dydz, vel_fac, s_shpvel.dydz);

    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dxdx, vel_fac_rhs * der2_vel_jump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dydy, vel_fac_rhs * der2_vel_jump(0,1));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dzdz, vel_fac_rhs * der2_vel_jump(0,2));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dxdy, vel_fac_rhs * der2_vel_jump(0,3));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dxdz, vel_fac_rhs * der2_vel_jump(0,4));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::master>(m_shpvel.dydz, vel_fac_rhs * der2_vel_jump(0,5));

    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dxdx, vel_fac_rhs * der2_vel_jump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dydy, vel_fac_rhs * der2_vel_jump(1,1));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dzdz, vel_fac_rhs * der2_vel_jump(1,2));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dxdy, vel_fac_rhs * der2_vel_jump(1,3));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dxdz, vel_fac_rhs * der2_vel_jump(1,4));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::master>(m_shpvel.dydz, vel_fac_rhs * der2_vel_jump(1,5));

    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dxdx, vel_fac_rhs * der2_vel_jump(2,0));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dydy, vel_fac_rhs * der2_vel_jump(2,1));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dzdz, vel_fac_rhs * der2_vel_jump(2,2));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dxdy, vel_fac_rhs * der2_vel_jump(2,3));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dxdz, vel_fac_rhs * der2_vel_jump(2,4));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::master>(m_shpvel.dydz, vel_fac_rhs * der2_vel_jump(2,5));

    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dxdx, -vel_fac_rhs * der2_vel_jump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dydy, -vel_fac_rhs * der2_vel_jump(0,1));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dzdz, -vel_fac_rhs * der2_vel_jump(0,2));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dxdy, -vel_fac_rhs * der2_vel_jump(0,3));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dxdz, -vel_fac_rhs * der2_vel_jump(0,4));
    assembler.template Vector<XFEM::PHYSICS::Velx,COMBUST::UTILS::slave>(s_shpvel.dydz, -vel_fac_rhs * der2_vel_jump(0,5));

    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dxdx, -vel_fac_rhs * der2_vel_jump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dydy, -vel_fac_rhs * der2_vel_jump(1,1));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dzdz, -vel_fac_rhs * der2_vel_jump(1,2));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dxdy, -vel_fac_rhs * der2_vel_jump(1,3));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dxdz, -vel_fac_rhs * der2_vel_jump(1,4));
    assembler.template Vector<XFEM::PHYSICS::Vely,COMBUST::UTILS::slave>(s_shpvel.dydz, -vel_fac_rhs * der2_vel_jump(1,5));

    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dxdx, -vel_fac_rhs * der2_vel_jump(2,0));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dydy, -vel_fac_rhs * der2_vel_jump(2,1));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dzdz, -vel_fac_rhs * der2_vel_jump(2,2));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dxdy, -vel_fac_rhs * der2_vel_jump(2,3));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dxdz, -vel_fac_rhs * der2_vel_jump(2,4));
    assembler.template Vector<XFEM::PHYSICS::Velz,COMBUST::UTILS::slave>(s_shpvel.dydz, -vel_fac_rhs * der2_vel_jump(2,5));

    const double pres_fac = tau_pres_second_order * fac * timefac;
    const double pres_fac_rhs = tau_pres_second_order * fac * rhstimefac;

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shppres.dxdx, pres_fac, m_shppres.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shppres.dydy, pres_fac, m_shppres.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shppres.dzdz, pres_fac, m_shppres.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shppres.dxdy, pres_fac, m_shppres.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shppres.dxdz, pres_fac, m_shppres.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::master>(m_shppres.dydz, pres_fac, m_shppres.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shppres.dxdx, -pres_fac, s_shppres.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shppres.dydy, -pres_fac, s_shppres.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shppres.dzdz, -pres_fac, s_shppres.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shppres.dxdy, -pres_fac, s_shppres.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shppres.dxdz, -pres_fac, s_shppres.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::master,COMBUST::UTILS::slave>(m_shppres.dydz, -pres_fac, s_shppres.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shppres.dxdx, -pres_fac, m_shppres.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shppres.dydy, -pres_fac, m_shppres.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shppres.dzdz, -pres_fac, m_shppres.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shppres.dxdy, -pres_fac, m_shppres.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shppres.dxdz, -pres_fac, m_shppres.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::master>(s_shppres.dydz, -pres_fac, m_shppres.dydz);

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shppres.dxdx, pres_fac, s_shppres.dxdx);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shppres.dydy, pres_fac, s_shppres.dydy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shppres.dzdz, pres_fac, s_shppres.dzdz);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shppres.dxdy, pres_fac, s_shppres.dxdy);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shppres.dxdz, pres_fac, s_shppres.dxdz);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres,COMBUST::UTILS::slave,COMBUST::UTILS::slave>(s_shppres.dydz, pres_fac, s_shppres.dydz);

    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::master>(m_shppres.dxdx, pres_fac_rhs*der2_pres_jump(0));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::master>(m_shppres.dydy, pres_fac_rhs*der2_pres_jump(1));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::master>(m_shppres.dzdz, pres_fac_rhs*der2_pres_jump(2));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::master>(m_shppres.dxdy, pres_fac_rhs*der2_pres_jump(3));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::master>(m_shppres.dxdz, pres_fac_rhs*der2_pres_jump(4));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::master>(m_shppres.dydz, pres_fac_rhs*der2_pres_jump(5));

    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::slave>(s_shppres.dxdx, -pres_fac_rhs*der2_pres_jump(0));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::slave>(s_shppres.dydy, -pres_fac_rhs*der2_pres_jump(1));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::slave>(s_shppres.dzdz, -pres_fac_rhs*der2_pres_jump(2));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::slave>(s_shppres.dxdy, -pres_fac_rhs*der2_pres_jump(3));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::slave>(s_shppres.dxdz, -pres_fac_rhs*der2_pres_jump(4));
    assembler.template Vector<XFEM::PHYSICS::Pres,COMBUST::UTILS::slave>(s_shppres.dydz, -pres_fac_rhs*der2_pres_jump(5));
  }


  return;
}

}// end namespace COMBUST
#endif
