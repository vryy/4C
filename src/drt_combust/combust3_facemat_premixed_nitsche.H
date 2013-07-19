/*----------------------------------------------------------------------*/
/*!
\file combust3_facemat_premixed_nitsche.H

\brief system matrix formulation for two-phase flow problem / 3D combustion element

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

// note: implementation has be carefully checked against implementation in
//       standard fluid environment -> without interface the same results are
//       obtained for all quantities evaluated in the preparation, for the
//       element matrices for the master, the salve and the coupling as well as
//       for the final patch matrix assembled in the finish step of the assembler (10.04.2013)

#ifndef COMBUST3_FACEMAT_PREMIXED_NITSCHE_H
#define COMBUST3_FACEMAT_PREMIXED_NITSCHE_H

#include "combust3.H"
#include "combust3_facemat_standard.H"
#include "combust3_stabilization.H"
#include "combust_time_integration_element.H"
#include "../drt_geometry/position_array.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_xfem/xfem_element_utils.H"

namespace COMBUST
{

template <DRT::Element::DiscretizationType DISTYPE,
          DRT::Element::DiscretizationType M_DISTYPE,
          DRT::Element::DiscretizationType S_DISTYPE,
          XFEM::AssemblyType M_ASSTYPE,
          XFEM::AssemblyType S_ASSTYPE,
          int NUMDOF,
          class M1, class M2, class V1, class V2>
void FaceMatNitsche(
     const DRT::ELEMENTS::Combust3IntFace*      ele,     ///< the element whose matrix is calculated
     COMBUST::LocalPatchAssembler<M_DISTYPE, S_DISTYPE, M_ASSTYPE, S_ASSTYPE, NUMDOF>& assembler,
     const M1&                                  evelnp_m,
     const M1&                                  eveln_m,
     const M1&                                  evelnm_m,
     const M1&                                  eaccn_m,
     const M1&                                  eaccam_m,
     const V1&                                  eprenp_m,
     const M2&                                  evelnp_s,
     const M2&                                  eveln_s,
     const M2&                                  evelnm_s,
     const M2&                                  eaccn_s,
     const M2&                                  eaccam_s,
     const V2&                                  eprenp_s,
     Teuchos::RCP<const MAT::Material>          material,      ///< fluid material
     const INPAR::FLUID::TimeIntegrationScheme  timealgo,      ///< time discretization type
     const double                               time,          ///< current time
     const double                               dt,            ///< delta t (time step size)
     const double                               theta,         ///< factor for one step theta scheme
     const double                               ga_alphaF,
     const double                               ga_alphaM,
     const double                               ga_gamma,
     const INPAR::FLUID::EOS_Pres               pres_stab,       ///< face term to evaluate: pressure stab
     const INPAR::FLUID::EOS_Conv_Stream        conv_stream_stab,///< face term to evaluate: conv stab stream
     const INPAR::FLUID::EOS_Conv_Cross         conv_cross_stab, ///< face term to evaluate: conv stab cross
     const INPAR::FLUID::EOS_Div                conti_stab,      ///< face term to evaluate: conti stab
     const INPAR::FLUID::EOS_ElementLength      hk_def,          ///< definition of characteristic element length
     const INPAR::FLUID::EOS_TauType            tau_def          ///< definition of stab parameter
)
{
  // hier kommt keine Parameter Liste an, denn die werden uebergeben
  // state groessen fur master und slave werden ebenfalls uebergeben

//  std::cout << "Jetzt bin ich in der mat" << std::endl;
//  std::cout << "evelnp_m " << evelnp_m << std::endl;
//  std::cout << "evelnp_s " << evelnp_s << std::endl;
//
//  std::cout << "eprenp_m " << eprenp_m << std::endl;
//  std::cout << "eprenp_s " << eprenp_s << std::endl;


  //--------------------------------------------------
  // get adjacent elements
  //--------------------------------------------------
  // get master and slave element
  const DRT::ELEMENTS::Combust3* master_ele = ele->ParentMasterElement();
  const DRT::ELEMENTS::Combust3* slave_ele = ele->ParentSlaveElement();

  //--------------------------------------------------
  // get coordinates of face and adjacent elements
  //--------------------------------------------------
  // space dimension for 3d fluid element
  // TODO: stimmt das auch fuer face
  const size_t nsd = 3;
  const size_t face_nsd = 2;

  // number of nodes for element, i.e., face
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // get node coordinates of the current element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  // get node coordinates of master and slave element
  if (M_DISTYPE != S_DISTYPE)
   dserror("Same element type for master and slave expected");
  const size_t ms_numnode = DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement;

  static LINALG::Matrix<nsd,ms_numnode> m_xyze;
  GEO::fillInitialPositionArray<M_DISTYPE>(master_ele,m_xyze);
  static LINALG::Matrix<nsd,ms_numnode> s_xyze;
  GEO::fillInitialPositionArray<S_DISTYPE>(slave_ele,s_xyze);

//  std::cout << "m_xyze " << m_xyze << std::endl;
//  std::cout << "s_xyze " << s_xyze << std::endl;
//  std::cout << "xyze " << xyze << std::endl;

  // flag for higher order elements
  //TODO: klaeren: kommen die ueberhaupt vor?
  const bool higher_order_ele = false; //DRT::UTILS::secondDerivativesZero<DISTYPE>();

  // get number of parameters (dofs) for each field on master and slave side
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t m_numparamvelx = XFEM::NumParam<ms_numnode,M_ASSTYPE>::get(*(master_ele->GetEleDofManager()), XFEM::PHYSICS::Velx);
//  std::cout << "m_numparamvelx " << m_numparamvelx << std::endl;
  const size_t m_numparampres = XFEM::NumParam<ms_numnode,M_ASSTYPE>::get(*(master_ele->GetEleDofManager()), XFEM::PHYSICS::Pres);
  dsassert(m_numparamvelx == m_numparampres, "assumption violation");
  const size_t s_numparamvelx = XFEM::NumParam<ms_numnode,S_ASSTYPE>::get(*(slave_ele->GetEleDofManager()), XFEM::PHYSICS::Velx);
  const size_t s_numparampres = XFEM::NumParam<ms_numnode,S_ASSTYPE>::get(*(slave_ele->GetEleDofManager()), XFEM::PHYSICS::Pres);
  dsassert(s_numparamvelx == s_numparampres, "assumption violation");

  //--------------------------------------------------
  // get material parameters
  //--------------------------------------------------
  // TODO: fuer den augenblick sehe ich mal beide vor
  //       aber eigentlich braucht man hier nur eine, weil
  //       - beide element im selben Gebiet liegen
  //       - fuer geschnitte Elemente zweimal, also getrennt fuer beide seiten integriert werden muss
  // density \rho^{n+1}
  double densaf  = -1.0;
  // dynamic viscosity \mu^{n+1}
  double dynvisc = -1.0;
  // set density and dynamic viscosity
  {
    double densplus = -1.0;
    double densminus = -1.0;
    double dynviscplus = -1.0;
    double dynviscminus = -1.0;
    COMBUST::GetMaterialParams(material, true, densplus, dynviscplus);
    COMBUST::GetMaterialParams(material, false, densminus, dynviscminus);

    // blend smoothly to minus (gas bubble) domain
    // TODO: klaeren
    //COMBUST::BlendMaterial(ele, &(*cell), time, densaf, dynvisc, densplus, dynviscplus, densminus, dynviscminus);

    if ((master_ele->Bisected() or master_ele->Touched()) and (slave_ele->Bisected() or slave_ele->Touched()))
      dserror("Cut elements not yet considered");

    densaf = densminus;
    dynvisc = dynviscplus;
  }

  //----------------------------------------------
  // get time integration parameters
  //----------------------------------------------
  // instationary formulation
  bool instationary = true;
  if (timealgo == INPAR::FLUID::timeint_stationary) instationary = false;
  // generalized alpha time integration scheme
  bool genalpha = false;
  if (timealgo == INPAR::FLUID::timeint_afgenalpha) genalpha = true;
  // time integration constant
  const double timefac = COMBUST::TIMEINT::ComputeTimeFac(timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma);
  // TODO: passt das?
  // integration factors and coefficients of single terms
  double rhstimefac = timefac;
  if (genalpha)
    rhstimefac /= ga_alphaF; // (fac * gamma * dt ) / alpha_M

  //----------------------------------------------
  // get Gaussian points
  //----------------------------------------------
  // TODO: quadratische Funktionen der Kinks beruecksichtigen
  const DRT::UTILS::GaussRule2D gaussrule = XFEM::getOptimalFaceGaussrule(DISTYPE);
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  //---------------------------------------------
  // transform coordinates of of Gaussian points
  //---------------------------------------------

  // get coordinates of nodes of face in master element coordinates (xi)
  const LINALG::Matrix<nsd,numnode> face_xi_master;

  // get coordinates of nodes of face in slave element coordinates (xi)
  const LINALG::Matrix<nsd,numnode> face_xi_slave;

  // transform the face's Gaussian points to both parent elements

  // number of Gaussian points
  int numgp = intpoints.nquad;

  // TODO: checken was das hier genau macht
  // local coordinates of the face's gausspoints w.r.t master and slave element
  Epetra_SerialDenseMatrix master_xi_points(numgp,nsd);
  Epetra_SerialDenseMatrix slave_xi_points(numgp,nsd);
  Epetra_SerialDenseMatrix face_xi_points_master(numgp,face_nsd);
  Epetra_SerialDenseMatrix face_xi_points_slave(numgp,face_nsd);

  // local coordinates of the face nodes w.r.t slave side
  LINALG::Matrix<face_nsd, numnode> local_slave_coordiantes_trafo(true);

  std::vector<int> localtrafomap = ele->GetLocalTrafoMap();

  for(std::size_t i=0; i< numnode; i++)
  {
    for(std::size_t isd= 0; isd< face_nsd; isd++)
    {
      switch(DISTYPE)
      {
      case DRT::Element::quad4:
      {
        local_slave_coordiantes_trafo(isd,localtrafomap[i]) = DRT::UTILS::eleNodeNumbering_quad9_nodes_reference[i][isd];
        break;
      }
      default:
      {
        dserror("ele type not supported %d", DISTYPE);
        break;
      }
      }
    }
  }

  // coordinates of all integration points as with local coordinates w.r.t the respective local side
  // of the respective parent element
  for(int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    LINALG::Matrix<face_nsd,1> face_xi_points_master_linalg(true);
    LINALG::Matrix<face_nsd,1> face_xi_points_slave_linalg(true);

    // Gaussian point in face's element's local coordinates w.r.t master element
    for (std::size_t idim=0;idim<face_nsd;idim++)
    {
      face_xi_points_master(iquad,idim) = intpoints.qxg[iquad][idim];
      face_xi_points_master_linalg(idim) = intpoints.qxg[iquad][idim];
    }

    // transform the local coordinates from the local coordinate system of the face w.r.t master face
    // to the local coordinate system of the face w.r.t slave face
//    DRT::UTILS::shape_function<distype>(face_xi_points_master_linalg,funct_);
    static LINALG::Matrix<numnode,1> funct;
    DRT::UTILS::shape_function_2D(funct,face_xi_points_master_linalg(0,0),face_xi_points_master_linalg(1,0),DISTYPE);

    face_xi_points_slave_linalg.Multiply(local_slave_coordiantes_trafo,funct);

    for (std::size_t idim=0;idim<face_nsd;idim++)
    {
      face_xi_points_slave(iquad,idim) = face_xi_points_slave_linalg(idim);
    }
  }

  // get the local gp coordinates w.r.t master element
  DRT::UTILS::BoundaryGPToParentGP3(master_xi_points,
      face_xi_points_master,
      M_DISTYPE,
      DISTYPE,
      ele->SurfaceMasterNumber());

  // get the local gp coordinates w.r.t slave element
  DRT::UTILS::BoundaryGPToParentGP3(slave_xi_points,
      face_xi_points_slave,
      S_DISTYPE,
      DISTYPE,
      ele->SurfaceSlaveNumber());

  //--------------------------------------------
  // prepare evaluation of enrichment function
  //--------------------------------------------
  // TODO: xfem klaeren
  // TODO: das muss wohl noch ueergeben werden
  const bool indomainplus =  false;
  LINALG::Matrix<ms_numnode,1> m_ephi(true);
  LINALG::Matrix<ms_numnode,1> s_ephi(true);
  for (std::size_t inode=0; inode<ms_numnode; inode++)
  {
    m_ephi(inode,0) = -1;
    s_ephi(inode,0) = -1;
  }
  // create dummy cell for master and slave
  LINALG::SerialDenseMatrix xietazetae(nsd,ms_numnode);
  LINALG::SerialDenseMatrix m_xyze_serial_dens_mat(nsd,ms_numnode);
  LINALG::SerialDenseMatrix s_xyze_serial_dens_mat(nsd,ms_numnode);
  for (std::size_t inode=0; inode<ms_numnode; inode++)
  {
    for (std::size_t idim=0; idim<nsd; idim++)
    {
      xietazetae(idim,inode) =  DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[inode][idim];
      m_xyze_serial_dens_mat(idim,inode) = m_xyze(idim,inode);
      s_xyze_serial_dens_mat(idim,inode) = s_xyze(idim,inode);
    }
  }
  GEO::DomainIntCell master_cell(M_DISTYPE,xietazetae,m_xyze_serial_dens_mat,indomainplus);
  GEO::DomainIntCell slave_cell(S_DISTYPE,xietazetae,s_xyze_serial_dens_mat,indomainplus);

  const XFEM::ElementEnrichmentValues master_enrvals(*master_ele,*(master_ele->GetEleDofManager()),master_cell,m_ephi);
  const XFEM::ElementEnrichmentValues slave_enrvals(*slave_ele,*(slave_ele->GetEleDofManager()),slave_cell,s_ephi);



  //----------------------------------------------
  // setup of stabilization parameter
  //----------------------------------------------
  // get characteristic element length
  const double hk = COMBUST::UTILS::computeHKfacestab<M_DISTYPE,S_DISTYPE,nsd>(master_ele, slave_ele, ele, m_xyze, s_xyze, hk_def);
//  std::cout << "hk  " << hk << std::endl;

  // get characteristic velocity
  double vel_norm = 0.0;
  // get the L_inf-norm of the parent's element velocity for stabilization
  // -> maximum value of all nodal values of the patch
  for(std::size_t idim=0; idim<nsd; idim++)
  {
    for(std::size_t inode=0; inode<ms_numnode; inode++)
    {
      if (fabs(evelnp_m(idim,inode)) > vel_norm)
        vel_norm = fabs(evelnp_m(idim,inode));
    }
  }
  for(std::size_t idim=0; idim<nsd; idim++)
  {
    for(std::size_t inode=0; inode<ms_numnode; inode++)
    {
      if (fabs(evelnp_s(idim,inode)) > vel_norm)
        vel_norm = fabs(evelnp_s(idim,inode));
    }
  }

//  std::cout << "vel_norm  " << vel_norm << std::endl;

  //--------------------------------------
  // integration loop over Gaussian points
  //--------------------------------------
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {

//    std::cout<< "gp" << std::endl;


    //--------------------------------------
    // evaluate shape functions w.r.t. master
    //--------------------------------------

    // shape functions and derivs of master element at gausspoint
    static LINALG::Matrix<ms_numnode,1> m_funct;
    static LINALG::Matrix<nsd,ms_numnode> m_deriv;
    DRT::UTILS::shape_function_3D(m_funct,master_xi_points(iquad,0),master_xi_points(iquad,1),master_xi_points(iquad,2),M_DISTYPE);
    DRT::UTILS::shape_function_3D_deriv1(m_deriv,master_xi_points(iquad,0),master_xi_points(iquad,1),master_xi_points(iquad,2),M_DISTYPE);

//    std::cout << "gp 1: master " << master_xi_points(iquad,0) << std::endl;
//    std::cout << "gp 2: master " << master_xi_points(iquad,1) << std::endl;
//    std::cout << "gp 3: master " << master_xi_points(iquad,2) << std::endl;
//
//    std::cout << "gp 1: slave " << slave_xi_points(iquad,0) << std::endl;
//    std::cout << "gp 2: slave " << slave_xi_points(iquad,1) << std::endl;
//    std::cout << "gp 3: slave " << slave_xi_points(iquad,2) << std::endl;

    // get transposed of the jacobian matrix d x / d \xi
    // xjm(i,j) = deriv(i,k)*xyze(j,k)
    static LINALG::Matrix<nsd,nsd> m_xjm;
    m_xjm.MultiplyNT(m_deriv,m_xyze);
    const double m_det = m_xjm.Determinant();
    if (m_det < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", master_ele->Id(), m_det);
    }
    // inverse of jacobian
    static LINALG::Matrix<nsd,nsd> m_xji;
    m_xji.Invert(m_xjm);

    // compute global derivates
    static LINALG::Matrix<3,ms_numnode> m_derxy;
    // derxy(i,j) = xji(i,k) * deriv(k,j)
    m_derxy.Multiply(m_xji,m_deriv);

    //--------------------------------------
    // evaluate shape functions w.r.t. slave
    //--------------------------------------

    // shape functions and derivs of slave element at gausspoint
    static LINALG::Matrix<ms_numnode,1> s_funct;
    static LINALG::Matrix<nsd,ms_numnode> s_deriv;
    DRT::UTILS::shape_function_3D(s_funct,slave_xi_points(iquad,0),slave_xi_points(iquad,1),slave_xi_points(iquad,2),S_DISTYPE);
    DRT::UTILS::shape_function_3D_deriv1(s_deriv,slave_xi_points(iquad,0),slave_xi_points(iquad,1),slave_xi_points(iquad,2),S_DISTYPE);

    // get transposed of the jacobian matrix d x / d \xi
    // xjm(i,j) = deriv(i,k)*xyze(j,k)
    static LINALG::Matrix<nsd,nsd>s_xjm;
    s_xjm.MultiplyNT(s_deriv,s_xyze);
    const double s_det = s_xjm.Determinant();
    if (s_det < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", slave_ele->Id(), s_det);
    }
    // inverse of jacobian
    static LINALG::Matrix<nsd,nsd> s_xji;
    s_xji.Invert(s_xjm);

    // compute global derivates
    static LINALG::Matrix<3,ms_numnode> s_derxy;
    // derxy(i,j) = xji(i,k) * deriv(k,j)
    s_derxy.Multiply(s_xji,s_deriv);

    //--------------------------------------
    // evaluate shape functions w.r.t. face
    //--------------------------------------

    // shape functions and derivs of face element at gausspoint
    static LINALG::Matrix<numnode,1> funct;
    static LINALG::Matrix<face_nsd,numnode> deriv;
    DRT::UTILS::shape_function_2D(funct,intpoints.qxg[iquad][0],intpoints.qxg[iquad][1],DISTYPE);
    DRT::UTILS::shape_function_2D_deriv1(deriv,intpoints.qxg[iquad][0],intpoints.qxg[iquad][1],DISTYPE);

    // jacobian(eta->x)*jacobian(eta->x)^T
    static LINALG::Matrix<face_nsd,face_nsd> metrictensor;
    // determinant of metric tensor
    double sqrtdetmetrictensor = 0.0;
    // outward pointing normal vector on face
    static LINALG::Matrix<nsd,1> normal;
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<DISTYPE>(xyze,deriv,metrictensor,sqrtdetmetrictensor,&normal);

    const double fac = intpoints.qwgt[iquad] * sqrtdetmetrictensor;

    //--------------------------------------
    // second order derivatives
    //--------------------------------------
    // compute second global derivative for master and slave
    static LINALG::Matrix<6,ms_numnode> m_derxy2;
    static LINALG::Matrix<6,ms_numnode> s_derxy2;
    if (higher_order_ele)
    {
      static LINALG::Matrix<6,ms_numnode> deriv2;
      DRT::UTILS::shape_function_3D_deriv2(deriv2,master_xi_points(iquad,0),master_xi_points(iquad,1),master_xi_points(iquad,2),M_DISTYPE);
      DRT::UTILS::gder2<M_DISTYPE>(m_xjm,m_derxy,deriv2,m_xyze,m_derxy2);
      deriv2.Clear();
      DRT::UTILS::shape_function_3D_deriv2(deriv2,slave_xi_points(iquad,0),slave_xi_points(iquad,1),slave_xi_points(iquad,2),M_DISTYPE);
      DRT::UTILS::gder2<M_DISTYPE>(s_xjm,s_derxy,deriv2,s_xyze,s_derxy2);
    }
    else
    {
      m_derxy2.Clear();
      s_derxy2.Clear();
    }

    //--------------------------------------------------------------------------------
    // rearrange (enriched) shape functions and derivatives as approximation functions
    //--------------------------------------------------------------------------------
    // distinguish between master and slave to used implementation not only for
    // xfem stability but also for edge-based stabilization
    const size_t m_shpVecSize       = COMBUST::SizeFac<M_ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement;
    const size_t s_shpVecSize       = COMBUST::SizeFac<S_ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement;

    //-------------------------------------------------
    // evaluate enriched shape functions w.r.t. master
    //-------------------------------------------------

    // (enriched) shape functions = approximation functions (P = N * \Psi)
    static XFEM::ApproxFunc<2,m_shpVecSize> m_shp;

    if (M_ASSTYPE == XFEM::xfem_assembly)
    {
      // temporary arrays holding enriched shape functions (N * \Psi)
      static LINALG::Matrix<m_shpVecSize,1> m_enr_funct;
      m_enr_funct.Clear();
      static LINALG::Matrix<3,m_shpVecSize> m_enr_derxy;
      m_enr_derxy.Clear();
      static LINALG::Matrix<6,m_shpVecSize> m_enr_derxy2;
      m_enr_derxy2.Clear();

      // shape functions and derivatives for nodal parameters (dofs)
      master_enrvals.ComputeModifiedEnrichedNodalShapefunction(
          XFEM::PHYSICS::Velx,
          m_funct,
          m_derxy,
          m_derxy2,
          m_enr_funct,
          m_enr_derxy,
          m_enr_derxy2);

      // fill approximation functions for XFEM
      for (size_t iparam = 0; iparam != m_numparamvelx; ++iparam)
      {
        m_shp.d0(iparam) = m_enr_funct(iparam);
        m_shp.dx(iparam) = m_enr_derxy(0,iparam);
        m_shp.dy(iparam) = m_enr_derxy(1,iparam);
        m_shp.dz(iparam) = m_enr_derxy(2,iparam);
        m_shp.dxdx(iparam) = m_enr_derxy2(0,iparam);
        m_shp.dxdy(iparam) = m_enr_derxy2(3,iparam);
        m_shp.dxdz(iparam) = m_enr_derxy2(4,iparam);
        m_shp.dydx(iparam) = m_shp.dxdy(iparam);
        m_shp.dydy(iparam) = m_enr_derxy2(1,iparam);
        m_shp.dydz(iparam) = m_enr_derxy2(5,iparam);
        m_shp.dzdx(iparam) = m_shp.dxdz(iparam);
        m_shp.dzdy(iparam) = m_shp.dydz(iparam);
        m_shp.dzdz(iparam) = m_enr_derxy2(2,iparam);
      }
    }
    else // not xfem_assembly i.e. standard assembly
    {
      // fill approximation functions for standard FEM
      // remark: numparamvelx == numnode, for standard FEM
      for (size_t iparam = 0; iparam < ms_numnode; ++iparam)
      {
        m_shp.d0(iparam) = m_funct(iparam);
        m_shp.dx(iparam) = m_derxy(0,iparam);
        m_shp.dy(iparam) = m_derxy(1,iparam);
        m_shp.dz(iparam) = m_derxy(2,iparam);
        m_shp.dxdx(iparam) = m_derxy2(0,iparam);
        m_shp.dxdy(iparam) = m_derxy2(3,iparam);
        m_shp.dxdz(iparam) = m_derxy2(4,iparam);
        m_shp.dydx(iparam) = m_shp.dxdy(iparam);
        m_shp.dydy(iparam) = m_derxy2(1,iparam);
        m_shp.dydz(iparam) = m_derxy2(5,iparam);
        m_shp.dzdx(iparam) = m_shp.dxdz(iparam);
        m_shp.dzdy(iparam) = m_shp.dydz(iparam);
        m_shp.dzdz(iparam) = m_derxy2(2,iparam);
//        std::cout << "m_shp.d0(iparam)" << m_shp.d0(iparam) << std::endl;
      }
    }

    //-------------------------------------------------
    // evaluate enriched shape functions w.r.t. slave
    //-------------------------------------------------

    // (enriched) shape functions = approximation functions (P = N * \Psi)
    static XFEM::ApproxFunc<2,s_shpVecSize> s_shp;

    if (S_ASSTYPE == XFEM::xfem_assembly)
    {
      // temporary arrays holding enriched shape functions (N * \Psi)
      static LINALG::Matrix<s_shpVecSize,1> s_enr_funct;
      s_enr_funct.Clear();
      static LINALG::Matrix<3,s_shpVecSize> s_enr_derxy;
      s_enr_derxy.Clear();
      static LINALG::Matrix<6,s_shpVecSize> s_enr_derxy2;
      s_enr_derxy2.Clear();

      // shape functions and derivatives for nodal parameters (dofs)
      slave_enrvals.ComputeModifiedEnrichedNodalShapefunction(
          XFEM::PHYSICS::Velx,
          s_funct,
          s_derxy,
          s_derxy2,
          s_enr_funct,
          s_enr_derxy,
          s_enr_derxy2);

      // fill approximation functions for XFEM
      for (size_t iparam = 0; iparam != s_numparamvelx; ++iparam)
      {
        s_shp.d0(iparam) = s_enr_funct(iparam);
        s_shp.dx(iparam) = s_enr_derxy(0,iparam);
        s_shp.dy(iparam) = s_enr_derxy(1,iparam);
        s_shp.dz(iparam) = s_enr_derxy(2,iparam);
        s_shp.dxdx(iparam) = s_enr_derxy2(0,iparam);
        s_shp.dxdy(iparam) = s_enr_derxy2(3,iparam);
        s_shp.dxdz(iparam) = s_enr_derxy2(4,iparam);
        s_shp.dydx(iparam) = s_shp.dxdy(iparam);
        s_shp.dydy(iparam) = s_enr_derxy2(1,iparam);
        s_shp.dydz(iparam) = s_enr_derxy2(5,iparam);
        s_shp.dzdx(iparam) = s_shp.dxdz(iparam);
        s_shp.dzdy(iparam) = s_shp.dydz(iparam);
        s_shp.dzdz(iparam) = s_enr_derxy2(2,iparam);
      }
    }
    else // not xfem_assembly i.e. standard assembly
    {
      // fill approximation functions for standard FEM
      // remark: numparamvelx == numnode, for standard FEM
      for (size_t iparam = 0; iparam < ms_numnode; ++iparam)
      {
        s_shp.d0(iparam) = s_funct(iparam);
        s_shp.dx(iparam) = s_derxy(0,iparam);
        s_shp.dy(iparam) = s_derxy(1,iparam);
        s_shp.dz(iparam) = s_derxy(2,iparam);
        s_shp.dxdx(iparam) = s_derxy2(0,iparam);
        s_shp.dxdy(iparam) = s_derxy2(3,iparam);
        s_shp.dxdz(iparam) = s_derxy2(4,iparam);
        s_shp.dydx(iparam) = s_shp.dxdy(iparam);
        s_shp.dydy(iparam) = s_derxy2(1,iparam);
        s_shp.dydz(iparam) = s_derxy2(5,iparam);
        s_shp.dzdx(iparam) = s_shp.dxdz(iparam);
        s_shp.dzdy(iparam) = s_shp.dydz(iparam);
        s_shp.dzdz(iparam) = s_derxy2(2,iparam);
//        std::cout << "s_shp.d0(iparam)" << s_shp.d0(iparam) << std::endl;
      }
    }

    // TODO: brauche ich hier auch face?

    //--------------------------------------------------
    // get velocity and pressure at integration point
    //--------------------------------------------------
    // since functions are continuous at element boundaries, one can
    // evaluate the velocity using the master element, the slave element or the face
    // TODO: was passiert bei touched elements
    const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPoint(evelnp_m, m_shp.d0, m_numparamvelx);
    const double gppresnp = XFEM::interpolateScalarFieldToIntPoint(eprenp_m, m_shp.d0, m_numparampres);
//    double gppresnp = 0.0;
//    for (size_t iparam = 0; iparam != m_numparampres; ++iparam)
//        gppresnp += eprenp_m(iparam) * m_shp.d0(iparam);
//    std::cout << "gpvelnp  " << gpvelnp << std::endl;
//    std::cout << "gppresnp  " << gppresnp << std::endl;

    //----------------------------------------------
    // get velocity derivatives at integration point
    //----------------------------------------------
    // master
    static LINALG::Matrix<3,nsd> m_vderxy;
    m_vderxy.Clear();
    // vderxy = enr_derxy(j,k)*evelnp(i,k);
    for (size_t iparam = 0; iparam < m_numparamvelx; ++iparam)
    {
      for (size_t isd = 0; isd < nsd; ++isd)
      {
        m_vderxy(isd,0) += evelnp_m(isd,iparam) * m_shp.dx(iparam);
        m_vderxy(isd,1) += evelnp_m(isd,iparam) * m_shp.dy(iparam);
        m_vderxy(isd,2) += evelnp_m(isd,iparam) * m_shp.dz(iparam);
      }
    }
    // slave
    static LINALG::Matrix<3,nsd> s_vderxy;
    s_vderxy.Clear();
    // vderxy = enr_derxy(j,k)*evelnp(i,k);
    for (size_t iparam = 0; iparam < s_numparamvelx; ++iparam)
    {
      for (size_t isd = 0; isd < nsd; ++isd)
      {
        s_vderxy(isd,0) += evelnp_s(isd,iparam) * s_shp.dx(iparam);
        s_vderxy(isd,1) += evelnp_s(isd,iparam) * s_shp.dy(iparam);
        s_vderxy(isd,2) += evelnp_s(isd,iparam) * s_shp.dz(iparam);
      }
    }

    //----------------------------------------------
    // get divergence jump at interface
    //----------------------------------------------
    // master
    const double m_vdiv = m_vderxy(0,0) + m_vderxy(1,1) + m_vderxy(2,2);
    // salve
    const double s_vdiv = s_vderxy(0,0) + s_vderxy(1,1) + s_vderxy(2,2);
    // vdiv_slave - vdiv_master
    const double vdiv_jump = s_vdiv - m_vdiv;
//    std::cout << "vdiv_jump  " << vdiv_jump << std::endl;

    //----------------------------------------------
    // get gradient jump of velocity at interface
    //----------------------------------------------
    LINALG::Matrix<nsd,nsd> grad_vel_jump(true);
    grad_vel_jump.Update(1.0, s_vderxy, -1.0, m_vderxy);
//    std::cout << "grad_vel_jump  " << grad_vel_jump << std::endl;

    //----------------------------------------------
    // projection of velocity in normal direction of face
    //----------------------------------------------
    const double normal_vel = fabs(gpvelnp.Dot(normal));
//    std::cout << "normal_vel  " << normal_vel << std::endl;

    //-----------------------------------------------------
    // get second velocity derivatives at integration point
    //-----------------------------------------------------
    // master and slave
    static LINALG::Matrix<3,6> m_vderxy2;
    m_vderxy2.Clear();
    static LINALG::Matrix<3,6> s_vderxy2;
    s_vderxy2.Clear();
    //vderxy2 = evelnp(i,k)*enr_derxy2(j,k);
    if (higher_order_ele)
    {
      for (size_t iparam = 0; iparam < s_numparamvelx; ++iparam)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
        {
          s_vderxy2(isd,0) += evelnp_s(isd,iparam) * s_shp.dxdx(iparam);
          s_vderxy2(isd,1) += evelnp_s(isd,iparam) * s_shp.dydy(iparam);
          s_vderxy2(isd,2) += evelnp_s(isd,iparam) * s_shp.dzdz(iparam);
          s_vderxy2(isd,3) += evelnp_s(isd,iparam) * s_shp.dxdy(iparam);
          s_vderxy2(isd,4) += evelnp_s(isd,iparam) * s_shp.dxdz(iparam);
          s_vderxy2(isd,5) += evelnp_s(isd,iparam) * s_shp.dydz(iparam);
        }
      }
    }

    //--------------------------------------------
    // get pressure gradients at integration point
    //--------------------------------------------
    // master
    static LINALG::Matrix<nsd,1> m_gradp;
    m_gradp.Clear();
    // gradp = enr_derxy(i,j)*eprenp(j);
    for (size_t iparam = 0; iparam != m_numparampres; ++iparam)
    {
      m_gradp(0) += eprenp_m(iparam) * m_shp.dx(iparam);
      m_gradp(1) += eprenp_m(iparam) * m_shp.dy(iparam);
      m_gradp(2) += eprenp_m(iparam) * m_shp.dz(iparam);
    }
    // slave
    static LINALG::Matrix<nsd,1> s_gradp;
    s_gradp.Clear();
    // gradp = enr_derxy(i,j)*eprenp(j);
    for (size_t iparam = 0; iparam != s_numparampres; ++iparam)
    {
      s_gradp(0) += eprenp_s(iparam) * s_shp.dx(iparam);
      s_gradp(1) += eprenp_s(iparam) * s_shp.dy(iparam);
      s_gradp(2) += eprenp_s(iparam) * s_shp.dz(iparam);
    }

    //--------------------------------------------
    // get pressure jump at integration point
    //--------------------------------------------
    // grad(p_slave) - grad(p_master)
    LINALG::Matrix<nsd,1> grad_pres_jump(true);
    grad_pres_jump.Update(1.0, s_gradp, -1.0, m_gradp);
//    std::cout << "grad_pres_jump  " << grad_pres_jump << std::endl;

    //------------------------------------------
    // compute stabilization parameters (3 taus)
    //------------------------------------------

    double tau_conv = 0.0;   // convection-stabilization parameter
    double tau_div  = 0.0;   // divergence-stabilization parameter
    double tau_p    = 0.0;   // pressure-stabilization parameter

    // TODO: brauche ich den? wird jetzt nicht sinnvoll berechnet
    double tau_ghost = 0.0;  // ghost-penalty stabilization parameter due to Nitsche's method in the XFEM

    COMBUST::UTILS::computeStabilizationParamsEdgeBased(
            dynvisc, densaf, vel_norm, hk,
            tau_def, tau_conv, tau_div, tau_p, tau_ghost);
//    std::cout << "tau_conv  " << tau_conv << std::endl;
//    std::cout << "tau_div  " << tau_div << std::endl;
//    std::cout << "tau_p  " << tau_p << std::endl;

    //----------------------------------------------
    // build stiffness matrix stabilization entries
    //----------------------------------------------
    // that is master-master, master-slave, slave-master and slave-slave block
    COMBUST::BuildStandardFaceIntegrals<M_DISTYPE, S_DISTYPE, M_ASSTYPE, S_ASSTYPE, NUMDOF, m_shpVecSize, s_shpVecSize>(
            assembler, m_shp, s_shp, m_shp, s_shp,
            fac, timefac, rhstimefac, genalpha, instationary, densaf, dynvisc,
            gpvelnp, m_vderxy, s_vderxy,
            gppresnp, s_gradp, s_gradp,
            grad_pres_jump, vdiv_jump, normal_vel, grad_vel_jump,
            pres_stab, conv_stream_stab, conv_cross_stab, conti_stab,
            tau_conv, tau_div, tau_p);

   // dserror("ENDE!");
  } // end: loop all gaussian points

//  std::cout << "done" <<std::endl;

  //----------------------------------------
  // finish assembly
  //----------------------------------------
  // assemble master-master, master-slave, slave-master and slave-slave block
  // into patch matrix
  assembler.Finish();

//    assembler.PrintMat();
//    assembler.PrintRhs();
//
//  dserror("ENDE!");

  return;
}

} // end namespace COMBUST

#endif
