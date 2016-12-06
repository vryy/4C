/*----------------------------------------------------------------------*/
/*!
\file combust3_stabilization.H

\brief stabilization related, templated functions

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 */
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_STABILIZATION_H
#define COMBUST3_STABILIZATION_H

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_element.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_xfem/enrichment_utils.H"
#include "combust3.H"


namespace COMBUST
{
namespace UTILS
{
  //! lines perpendicular to surface of Hex element
  //! 6 surfaces
  //! 4 lines each
  const int get_lines_perpendicular_surf_hex_ele[6][4]={
                                        { 4,  5,  6,  7},
                                        { 1,  3,  9, 11},
                                        { 0,  2,  8, 10},
                                        { 1,  3,  9, 11},
                                        { 0,  2,  8, 10},
                                        { 4,  5,  6,  7}};

  //! for each distype there should be an m_k needed for stabilization computation
  template<DRT::Element::DiscretizationType DISTYPE>
  double MK()
  {
    switch (DISTYPE)
    {
      case DRT::Element::tet4:
      case DRT::Element::pyramid5:
      case DRT::Element::hex8:
      case DRT::Element::wedge6:
        return 0.333333333333333333333; // 1/3
        break;
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      case DRT::Element::tet10:
      case DRT::Element::wedge15:
        return 0.083333333333333333333; // 1/12
        break;
      default:
        dserror("type unknown!\n");
    }
    return -1.0;
  }


  //! get one point gauss rule to calculate tau at element center
  template<DRT::Element::DiscretizationType DISTYPE>
  DRT::UTILS::GaussRule3D getIntegrationRuleForStabilization()
  {
    switch (DISTYPE)
    {
      case DRT::Element::hex8:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
        return DRT::UTILS::intrule_hex_1point;
        break;
      case DRT::Element::tet4:
      case DRT::Element::tet10:
        return DRT::UTILS::intrule_tet_1point;
        break;
      case DRT::Element::wedge6:
      case DRT::Element::wedge15:
        return DRT::UTILS::intrule_wedge_1point;
        break;
      case DRT::Element::pyramid5:
        return DRT::UTILS::intrule_pyramid_1point;
        break;
      default:
        dserror("invalid discretization type for fluid3");
        return DRT::UTILS::intrule3D_undefined;
    }
  }


  //! calculate characteristic element size for stabilization
  template <class V>
  double Streamlength(
      const V&     shp_dx,
      const V&     shp_dy,
      const V&     shp_dz,
      const LINALG::Matrix<3,1>& gpvelnp,  /// velocity at Gaussian point
      const double               vel_norm, /// velocity norm
      const int                  numparamvelx
  )
  {
    // normed velocity at Gauss point
    LINALG::Matrix<3,1> velino;
    if (vel_norm >= 1e-6)
    {
      velino.Update(1.0/vel_norm,gpvelnp);
    }
    else
    {
      //std::cout << "/!\\ set x-velocity for streamlength computation" << std::endl;
      velino.Clear();
      velino(0) = 1.0;
    }

    // get streamlength
    //const double val = sum(abs(velino(j)*derxy(j,i)));
    LINALG::SerialDenseVector velinoder(numparamvelx,true);
    for (int iparam = 0; iparam < numparamvelx; ++iparam)
    {
      velinoder(iparam) += velino(0)*shp_dx(iparam);
      velinoder(iparam) += velino(1)*shp_dy(iparam);
      velinoder(iparam) += velino(2)*shp_dz(iparam);
    }
    // remark: "const double val = velinoder.Norm1();" ;
    double val = 0.0;
    for (int iparam = 0; iparam < numparamvelx; ++iparam)
      val += fabs(velinoder(iparam));

    if (val == 0.0) dserror("infinite stream length detected!");
    const double strle = 2.0/val;

    return strle;
  }

  //! calculate characteristic element size for stabilization
  // to be used for XFEM problems with normal enrichment strategy
  template <class V>
  double Streamlength(
      const V&                   shp,
      const LINALG::Matrix<3,1>& gpvelnp,  /// velocity at Gaussian point
      const double               vel_norm, /// velocity norm
      const int                  numparamvel
  )
  {
#ifdef DEBUG
  if (numparamvel!=8) dserror("This will only work for hex8 elements");
#endif
    // normed velocity at Gauss point
    LINALG::Matrix<3,1> velino(true);
    if (vel_norm >= 1e-6)
    {
      velino.Update(1.0/vel_norm,gpvelnp);
    }
    else
    {
      //std::cout << "/!\\ set x-velocity for streamlength computation" << std::endl;
      velino.Clear();
      velino(0) = 1.0;
    }

    // get streamlength
    //const double val = sum(abs(velino(j)*derxy(j,i)));
    // remark: extend vector from 8 to 16 to bring in normal shape functions correctly
    LINALG::SerialDenseVector velinoder(numparamvel+8,true);
    for (int iparam = 0; iparam < numparamvel; ++iparam)
    {
      velinoder(iparam) += velino(0)*shp.velx.dx.s(iparam);
      velinoder(iparam) += velino(1)*shp.vely.dy.s(iparam);
      velinoder(iparam) += velino(2)*shp.velz.dz.s(iparam);

      // store values belonging to normal shape functions
      velinoder(iparam+8) += velino(0)*shp.velx.dx.n(iparam);
      velinoder(iparam+8) += velino(1)*shp.vely.dy.n(iparam);
      velinoder(iparam+8) += velino(2)*shp.velz.dz.n(iparam);
    }
    // remark: "const double val = velinoder.Norm1();" ;
    double val = 0.0;
    for (int iparam = 0; iparam < numparamvel+8; ++iparam)
      val += fabs(velinoder(iparam));

    if (val == 0.0) dserror("infinite stream length detected!");
    const double strle = 2.0/val;

    return strle;
  }


  //! calculate volume based on numerical integration
  template<DRT::Element::DiscretizationType DISTYPE, class M2>
  double VolumeViaNumIntegration(
      const M2&  xyze     ///< element nodal position array
  )
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // use one point integration rule to calculate hk at element center
    const DRT::UTILS::GaussRule3D integrationrule_stabili = getIntegrationRuleForStabilization<DISTYPE>();

    // integration points
    const DRT::UTILS::IntegrationPoints3D intpoints(integrationrule_stabili);

    // shape functions and derivs at element center
    LINALG::Matrix<3,1> e;
    e(0) = intpoints.qxg[0][0];
    e(1) = intpoints.qxg[0][1];
    e(2) = intpoints.qxg[0][2];
    const double wquad = intpoints.qwgt[0];

    LINALG::Matrix<3,numnode> deriv;
    DRT::UTILS::shape_function_3D_deriv1(deriv, e(0), e(1), e(2), DISTYPE);

    // get Jacobian matrix and determinant
    // xjm_ = deriv_(i,k)*xyze(j,k);
    LINALG::Matrix<3,3> xjm;
    xjm.MultiplyNT(deriv,xyze);

    const double vol = wquad * xjm.Determinant();

    return vol;
  }

  //! calculate characteristic element size for stabilization
  template<DRT::Element::DiscretizationType DISTYPE, class M1>
  double HK(
      const M1&  xyze     ///< element nodal position array
  )
  {
    const double vol = VolumeViaNumIntegration<DISTYPE>(xyze);

    // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
    const double hk = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

    return hk;
  }


  template<DRT::Element::DiscretizationType M_DISTYPE,
           DRT::Element::DiscretizationType S_DISTYPE,
           int dim, class M1, class M2>
  double computeHKfacestab(
     const DRT::ELEMENTS::Combust3*          master,
     const DRT::ELEMENTS::Combust3*          slave,
     const DRT::ELEMENTS::Combust3IntFace*   ele,
     M1&                                     m_xyze,
     M2&                                     s_xyze,
     const INPAR::FLUID::EOS_ElementLength&  eos_element_length
     )
  {
    double hk = 0.0;

    // do some checks first
    if (M_DISTYPE != S_DISTYPE)
     dserror("Same element type for master and slave expected");

    // max distance to opposite surface
    if(eos_element_length == INPAR::FLUID::EOS_he_max_dist_to_opp_surf)
    {
      // get all lines perpendicular to face with respect to master element
      int face_id_wrt_master = ele->FaceMasterNumber();
      std::vector<int> lines_master;

      // get all lines perpendicular to face with respect to slave element
      int face_id_wrt_slave = ele->FaceSlaveNumber();
      std::vector<int> lines_slave;

      switch (M_DISTYPE)
      {
        case DRT::Element::hex8:
        {
          // we have 4 lines for hex elements
          lines_master.resize(4);
          lines_slave.resize(4);
          for (int i=0; i<4; i++)
          {
            lines_master[i] = get_lines_perpendicular_surf_hex_ele[face_id_wrt_master][i];
            lines_slave[i] = get_lines_perpendicular_surf_hex_ele[face_id_wrt_slave][i];
          }
          break;
        }
        default:
        {
          dserror("Unknown element type!");
          break;
        }
      }

      // get nodes of lines
      std::vector<std::vector<int> > eleNodeNumbering_hex_lines = DRT::UTILS::getEleNodeNumberingLines(M_DISTYPE);

      // do the master element first
      for (std::size_t iline=0; iline<lines_master.size(); iline++)
      {
        const int line_id = lines_master[iline];
        const int node_1 = eleNodeNumbering_hex_lines[line_id][0];
        const int node_2 = eleNodeNumbering_hex_lines[line_id][1];

        LINALG::Matrix<dim,1> line_vec (true);
        for (int idim=0; idim<dim; idim++)
            line_vec(idim,0) = m_xyze(idim,node_2) - m_xyze(idim,node_1);

        const double distance = line_vec.Norm2();

        // take the longest distance to the other surface
        hk = std::max(hk,distance);
      }

      // repeat for the slave element
      for (std::size_t iline=0; iline<lines_slave.size(); iline++)
      {
        const int line_id = lines_slave[iline];
        const int node_1 = eleNodeNumbering_hex_lines[line_id][0];
        const int node_2 = eleNodeNumbering_hex_lines[line_id][1];

        LINALG::Matrix<dim,1> line_vec (true);
        for (int idim=0; idim<dim; idim++)
            line_vec(idim,0) = s_xyze(idim,node_2) -s_xyze(idim,node_1);

        const double distance = line_vec.Norm2();

        // take the longest distance to the other surface
        hk = std::max(hk,distance);
      }
    }
    // largest surface
    else if(eos_element_length == INPAR::FLUID::EOS_he_surf_with_max_diameter)
    {
      dserror("EOS_he_surf_with_max_diameter not considered for combustion problems");
    }
    else if (eos_element_length == INPAR::FLUID::EOS_he_vol_eq_diameter)
    {
      // master
      hk = HK<M_DISTYPE>(m_xyze);
      // slave
      hk = std::max(hk,HK<S_DISTYPE>(s_xyze));
    }
    else
      dserror("Unknown reference length for face stabilization in combust!");

    return hk;
  }


  void computeStabilizationParams(
      const LINALG::Matrix<3,1>& gpvelnp,  /// velocity at Gaussian point
      const LINALG::Matrix<3,3>& xji,      /// inverse of transposed Jacobian matrix
      const bool   instationary,
      const double dynvisc,                /// dyamic viscosity
      const double dens,                   /// density
      const double vel_norm,
      const double strle,
      const double hk,
      const double mk,
      const double timefac,
      const double dt,
      const enum INPAR::FLUID::TauType tautype,
      double& tau_stab_Mu,
      double& tau_stab_Mp,
      double& tau_stab_C
      );


  void computeStabilizationParamsEdgeBased(
      const double dynvisc,                /// dyamic viscosity
      const double dens,                   /// density
      const double vel_norm,
      const double normal_vel,             /// projection of velocity in normal direction of face
      const double hk,
      const enum INPAR::FLUID::EOS_TauType tautype,
      double& tau_conv,
      double& tau_div,
      double& tau_p,
      const bool instationary,
      const double timefac,
      const bool add_ghost_penalties,
      double& tau_ghost_first_order,
      double& tau_ghost_second_order,
      double& tau_pres_second_order
      );


} // namespace UTILS
} // namespace COMBUST

#endif
