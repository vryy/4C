/*----------------------------------------------------------------------*/
/*!
\file combust3_facemat_twophaseflow.H

\brief system matrix formulation for two-phase flow problem / 3D combustion element

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_FACEMAT_TWOPHASEFLOW_H
#define COMBUST3_FACEMAT_TWOPHASEFLOW_H

#include "combust3.H"
#include "combust3_facemat_standard.H"
#include "combust3_stabilization.H"
#include "combust_time_integration_element.H"
#include "../drt_geometry/position_array.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_xfem/xfem_element_utils.H"
#include "../drt_xfem/xfem_utils_combust.H"
#include "../drt_inpar/inpar_xfem.H"

namespace COMBUST
{

template <DRT::Element::DiscretizationType DISTYPE,
          DRT::Element::DiscretizationType M_DISTYPE,
          DRT::Element::DiscretizationType S_DISTYPE,
          XFEM::AssemblyType M_ASSTYPE,
          XFEM::AssemblyType S_ASSTYPE,
          int NUMDOF,
          class M1, class M2, class V1, class V2, class V3, class V4>
void FaceMatTwoPhase(
     const DRT::ELEMENTS::Combust3IntFace*      ele,     ///< the element whose matrix is calculated
     const COMBUST::InterfaceHandleCombust*     ih,      ///< connection to the interface handler
     COMBUST::LocalPatchAssembler<M_DISTYPE, S_DISTYPE, M_ASSTYPE, S_ASSTYPE, NUMDOF>& assembler,
     const M1&                                  evelnp_m,
     const V1&                                  eprenp_m,
     const V2&                                  ephinp_m,
     const M2&                                  evelnp_s,
     const V3&                                  eprenp_s,
     const V4&                                  ephinp_s,
     Teuchos::RCP<const MAT::Material>          material,      ///< fluid material
     const INPAR::FLUID::TimeIntegrationScheme  timealgo,      ///< time discretization type
     const double                               time,          ///< current time
     const double                               dt,            ///< delta t (time step size)
     const double                               theta,         ///< factor for one step theta scheme
     const double                               ga_alphaF,
     const double                               ga_alphaM,
     const double                               ga_gamma,
     const bool                                 xfem_stab,       ///< flag to indicate whether ghost penalties should be used
     const INPAR::XFEM::FaceType                facetype,        ///< flag to define faces which have to be integrated twice
     const INPAR::FLUID::EOS_Pres               pres_stab,       ///< face term to evaluate: pressure stab
     const INPAR::FLUID::EOS_Conv_Stream        conv_stream_stab,///< face term to evaluate: conv stab stream
     const INPAR::FLUID::EOS_Conv_Cross         conv_cross_stab, ///< face term to evaluate: conv stab cross
     const INPAR::FLUID::EOS_Div                conti_stab,      ///< face term to evaluate: conti stab
     const INPAR::FLUID::EOS_ElementLength      hk_def,          ///< definition of characteristic element length
     const INPAR::FLUID::EOS_TauType            tau_def          ///< definition of stab parameter
)
{
  // currently merely working for pressure enrichments only

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

  // flag for higher order elements
  // currently not used
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  // get number of parameters (dofs) for each field on master and slave side
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t m_numparamvelx = XFEM::NumParam<ms_numnode,M_ASSTYPE>::get(*(master_ele->GetEleDofManager()), XFEM::PHYSICS::Velx);
  const size_t m_numparampres = XFEM::NumParam<ms_numnode,M_ASSTYPE>::get(*(master_ele->GetEleDofManager()), XFEM::PHYSICS::Pres);
  const size_t s_numparamvelx = XFEM::NumParam<ms_numnode,S_ASSTYPE>::get(*(slave_ele->GetEleDofManager()), XFEM::PHYSICS::Velx);
  const size_t s_numparampres = XFEM::NumParam<ms_numnode,S_ASSTYPE>::get(*(slave_ele->GetEleDofManager()), XFEM::PHYSICS::Pres);
  // this check is only valid if both fields are enriched
  // debug only
  dsassert(m_numparamvelx == m_numparampres, "assumption violation");
  dsassert(s_numparamvelx == s_numparampres, "assumption violation");

  //--------------------------------------------------
  // get material parameters
  //--------------------------------------------------
  // determine domain first
  bool indomainplus = true;
  // density \rho^{n+1}
  double densaf  = -1.0;
  // dynamic viscosity \mu^{n+1}
  double dynvisc = -1.0;
  // furthermore, we decide if ghost penalties have to be included for the current face
  bool add_ghost_penalties = false;
  // set density and dynamic viscosity
  {
    double densplus = -1.0;
    double densminus = -1.0;
    double dynviscplus = -1.0;
    double dynviscminus = -1.0;
    COMBUST::GetMaterialParams(material, true, densplus, dynviscplus);
    COMBUST::GetMaterialParams(material, false, densminus, dynviscminus);

    if (facetype == INPAR::XFEM::face_type_std)
    {
      //face is not intersected -> master and/or slave is also not intersected
      if ( (not master_ele->Bisected()) or (not slave_ele->Bisected()))
      {
        if (not master_ele->Bisected())
        {
          // information about domain integration cells
          const GEO::DomainIntCells&  domainIntCells(ih->ElementDomainIntCells(master_ele->Id()));
          // if element is uncut, we usually have only one cell, which is the first one
          // (exception if small cells have been removed in flamefront)
          GEO::DomainIntCells::const_iterator itcell = domainIntCells.begin();
          if (itcell == domainIntCells.end()) dserror("Domain integration cells expected for master!");
            indomainplus = itcell->getDomainPlus();

          // do not ask the level set field directly:
          // If the interface cuts off very, very small pieces, they may have been removed
          // in the cut algorithm or the flame front, meaning the element is uncut. However,
          // the level-set field is not apdated to this simplification. Therefore, we have to
          // ask the interface handler here not to run into this inconsistency.
          //          for (std::size_t inode = 0; inode < ms_numnode; inode++)
          //          {
          //              std::cout << ephinp_m(inode,0) << std::endl;
          //            if (XFEM::plusDomain(ephinp_m(inode,0)) == false)
          //            {
          //              indomainplus = false;
          //              break;
          //            }
          //          }
        }
        else
        {
          // information about domain from integration cells
          const GEO::DomainIntCells&  domainIntCells(ih->ElementDomainIntCells(slave_ele->Id()));
          // if element is uncut, we usually have only one cell, which is the first one
          // (exception if small cells have been removed in flamefront)
          GEO::DomainIntCells::const_iterator itcell = domainIntCells.begin();
          if (itcell == domainIntCells.end()) dserror("Domain integration cells expected for slave!");
              indomainplus = itcell->getDomainPlus();
        }

        // do we have an enriched element
        // meaning, that either master or slave have to be cut
        // in this case, we want to include ghost penalties
        if ( master_ele->Bisected() or slave_ele->Bisected())
          add_ghost_penalties = true;
      }
      else dserror("Neither master nor salve is uncut!");
    }
    else if (facetype == INPAR::XFEM::face_type_double_plus)
    {
      indomainplus = true;
      add_ghost_penalties = true;
    }
    else if (facetype == INPAR::XFEM::face_type_double_minus)
    {
      indomainplus = false;
      add_ghost_penalties = true;
    }
    else
      dserror("Unknown face type!");

    // undo setting for ghost penalties in case, we do not
    // want those terms but merely use a basic edge-based
    // stabilization also for the intersected terms
    if (not xfem_stab) add_ghost_penalties = false;

    // blend smoothly to minus domain if required
    COMBUST::BlendMaterial(ele, indomainplus, time, densaf, dynvisc, densplus, dynviscplus, densminus, dynviscminus);

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
  double timefac = 1.0; // stationary
  double rhstimefac = 1.0; // stationary
  if (instationary) // ost and bdf2
  {
    timefac = dt;
    rhstimefac = timefac;
  }

  if (genalpha) // overwritten for gen-alpha
  {
    timefac = COMBUST::TIMEINT::ComputeTimeFac(timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma);
    rhstimefac = timefac/ga_alphaF; // (gamma * dt ) / alpha_M
  }

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
      ele->FaceMasterNumber());

  // get the local gp coordinates w.r.t slave element
  DRT::UTILS::BoundaryGPToParentGP3(slave_xi_points,
      face_xi_points_slave,
      S_DISTYPE,
      DISTYPE,
      ele->FaceSlaveNumber());

  //--------------------------------------------
  // prepare evaluation of enrichment function
  //--------------------------------------------

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

  const XFEM::ElementEnrichmentValues master_enrvals(*master_ele,*(master_ele->GetEleDofManager()),master_cell,ephinp_m);
  const XFEM::ElementEnrichmentValues slave_enrvals(*slave_ele,*(slave_ele->GetEleDofManager()),slave_cell,ephinp_s);

  //----------------------------------------------
  // setup of stabilization parameter
  //----------------------------------------------
  // get characteristic element length
  const double hk = COMBUST::UTILS::computeHKfacestab<M_DISTYPE,S_DISTYPE,nsd>(master_ele, slave_ele, ele, m_xyze, s_xyze, hk_def);

  // get characteristic velocity
  double vel_norm = 0.0;
  // get the L_inf-norm of the parent's element velocity for stabilization
  // -> maximum value of all nodal values of the patch (including extensions for ghost penalties)

  if (M_ASSTYPE == XFEM::xfem_assembly)
  {
    const size_t m_shpVecSize_init = COMBUST::SizeFac<M_ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement;
    static LINALG::Matrix<m_shpVecSize_init,1> m_nodal_enr_funct(true);
    // shape function or enrichemnt function at node (either 1 (std dof), -2 or 2 (enr dof influencing current domain) or 0 (enr dof not influencing current domain)
    m_nodal_enr_funct.Clear();
    // number of field (i.e., velx) dofs per node (e.g., std + enr)
    std::vector<int> m_num_field_enr_per_node(ms_numnode);
    master_enrvals.GetNodalEnrichmentsFunctionValues(XFEM::PHYSICS::Velx,m_num_field_enr_per_node,m_nodal_enr_funct);

    // vector for nodal velocity
    // nodal velocity = std dof + enrcihment funct * enr dof
    LINALG::Matrix<nsd,ms_numnode> m_node_vel(true);

    std::size_t position = 0;
    for(std::size_t inode=0; inode<ms_numnode; inode++)
    {
      for (int rr=0; rr<m_num_field_enr_per_node[inode]; rr++)
      {
        for (std::size_t idim=0; idim<nsd; idim++)
        {
          m_node_vel(idim,inode) += m_nodal_enr_funct(position,0) * evelnp_m(idim,position);
          //std::cout << "inode  " << inode << "  pos  " << position << " m_nodal_enr_funct(position,0) " << m_nodal_enr_funct(position,0) << " evelnp_m(idim,position) " << evelnp_m(idim,position) << std::endl;
        }

        position++;
      }
    }
    if (position != m_numparamvelx) dserror("Unexpected number of nodal values");

    for(std::size_t idim=0; idim<nsd; idim++)
    {
      for(std::size_t inode=0; inode<ms_numnode; inode++)
      {
        if (fabs(m_node_vel(idim,inode)) > vel_norm)
          vel_norm = fabs(m_node_vel(idim,inode));
      }
    }
  }
  else
  {
    for(std::size_t idim=0; idim<nsd; idim++)
    {
      for(std::size_t inode=0; inode<ms_numnode; inode++)
      {
        if (fabs(evelnp_m(idim,inode)) > vel_norm)
          vel_norm = fabs(evelnp_m(idim,inode));
      }
    }
  }

  if (S_ASSTYPE == XFEM::xfem_assembly)
  {
    const size_t s_shpVecSize_init = COMBUST::SizeFac<S_ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement;

    static LINALG::Matrix<s_shpVecSize_init,1> s_nodal_enr_funct(true);
    s_nodal_enr_funct.Clear();
   // number of field (i.e., velx) dofs per node (e.g., std + enr)
    std::vector<int> s_num_field_enr_per_node(ms_numnode);
    slave_enrvals.GetNodalEnrichmentsFunctionValues(XFEM::PHYSICS::Velx,s_num_field_enr_per_node,s_nodal_enr_funct);

    // vector for nodal velocity
    LINALG::Matrix<nsd,ms_numnode> s_node_vel(true);

    std::size_t position = 0;
    for(std::size_t inode=0; inode<ms_numnode; inode++)
    {
      for (int rr=0; rr<s_num_field_enr_per_node[inode]; rr++)
      {
        for (std::size_t idim=0; idim<nsd; idim++)
          s_node_vel(idim,inode) += s_nodal_enr_funct(position,0) * evelnp_s(idim,position);

        position++;
      }
    }
    if (position != s_numparamvelx) dserror("Unexpected number of nodal values");

    for(std::size_t idim=0; idim<nsd; idim++)
    {
      for(std::size_t inode=0; inode<ms_numnode; inode++)
      {
        if (fabs(s_node_vel(idim,inode)) > vel_norm)
          vel_norm = fabs(s_node_vel(idim,inode));
      }
    }
  }
  else
  {
    for(std::size_t idim=0; idim<nsd; idim++)
    {
      for(std::size_t inode=0; inode<ms_numnode; inode++)
      {
        if (fabs(evelnp_s(idim,inode)) > vel_norm)
          vel_norm = fabs(evelnp_s(idim,inode));
      }
    }
  }

  //--------------------------------------
  // integration loop over Gaussian points
  //--------------------------------------
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    //--------------------------------------
    // evaluate shape functions w.r.t. master
    //--------------------------------------

    // shape functions and derivs of master element at gausspoint
    static LINALG::Matrix<ms_numnode,1> m_funct;
    static LINALG::Matrix<nsd,ms_numnode> m_deriv;
    DRT::UTILS::shape_function_3D(m_funct,master_xi_points(iquad,0),master_xi_points(iquad,1),master_xi_points(iquad,2),M_DISTYPE);
    DRT::UTILS::shape_function_3D_deriv1(m_deriv,master_xi_points(iquad,0),master_xi_points(iquad,1),master_xi_points(iquad,2),M_DISTYPE);

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
      DRT::UTILS::gder2<M_DISTYPE,ms_numnode>(m_xjm,m_derxy,deriv2,m_xyze,m_derxy2);
      deriv2.Clear();
      DRT::UTILS::shape_function_3D_deriv2(deriv2,slave_xi_points(iquad,0),slave_xi_points(iquad,1),slave_xi_points(iquad,2),S_DISTYPE);
      DRT::UTILS::gder2<S_DISTYPE,ms_numnode>(s_xjm,s_derxy,deriv2,s_xyze,s_derxy2);
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
    static XFEM::ApproxFunc<2,m_shpVecSize> m_shp_vel;
    static XFEM::ApproxFunc<2,m_shpVecSize> m_shp_pres;

    if (M_ASSTYPE == XFEM::xfem_assembly)
    {
      // temporary arrays holding enriched shape functions (N * \Psi)
      static LINALG::Matrix<m_shpVecSize,1> m_enr_funct_vel;
      m_enr_funct_vel.Clear();
      static LINALG::Matrix<3,m_shpVecSize> m_enr_derxy_vel;
      m_enr_derxy_vel.Clear();
      static LINALG::Matrix<6,m_shpVecSize> m_enr_derxy2_vel;
      m_enr_derxy2_vel.Clear();
      static LINALG::Matrix<m_shpVecSize,1> m_enr_funct_pres;
      m_enr_funct_pres.Clear();
      static LINALG::Matrix<3,m_shpVecSize> m_enr_derxy_pres;
      m_enr_derxy_pres.Clear();
      static LINALG::Matrix<6,m_shpVecSize> m_enr_derxy2_pres;
      m_enr_derxy2_pres.Clear();

      // shape functions and derivatives for nodal parameters (dofs)
      master_enrvals.ComputeModifiedEnrichedNodalShapefunction(
          XFEM::PHYSICS::Velx,
          m_funct,
          m_derxy,
          m_derxy2,
          m_enr_funct_vel,
          m_enr_derxy_vel,
          m_enr_derxy2_vel);

      // fill approximation functions for XFEM
      for (size_t iparam = 0; iparam != m_numparamvelx; ++iparam)
      {
        m_shp_vel.d0(iparam) = m_enr_funct_vel(iparam);
        m_shp_vel.dx(iparam) = m_enr_derxy_vel(0,iparam);
        m_shp_vel.dy(iparam) = m_enr_derxy_vel(1,iparam);
        m_shp_vel.dz(iparam) = m_enr_derxy_vel(2,iparam);
        m_shp_vel.dxdx(iparam) = m_enr_derxy2_vel(0,iparam);
        m_shp_vel.dxdy(iparam) = m_enr_derxy2_vel(3,iparam);
        m_shp_vel.dxdz(iparam) = m_enr_derxy2_vel(4,iparam);
        m_shp_vel.dydx(iparam) = m_shp_vel.dxdy(iparam);
        m_shp_vel.dydy(iparam) = m_enr_derxy2_vel(1,iparam);
        m_shp_vel.dydz(iparam) = m_enr_derxy2_vel(5,iparam);
        m_shp_vel.dzdx(iparam) = m_shp_vel.dxdz(iparam);
        m_shp_vel.dzdy(iparam) = m_shp_vel.dydz(iparam);
        m_shp_vel.dzdz(iparam) = m_enr_derxy2_vel(2,iparam);
      }

      // shape functions and derivatives for nodal parameters (dofs)
      master_enrvals.ComputeModifiedEnrichedNodalShapefunction(
          XFEM::PHYSICS::Pres,
          m_funct,
          m_derxy,
          m_derxy2,
          m_enr_funct_pres,
          m_enr_derxy_pres,
          m_enr_derxy2_pres);

      // fill approximation functions for XFEM
      for (size_t iparam = 0; iparam != m_numparampres; ++iparam)
      {
        m_shp_pres.d0(iparam) = m_enr_funct_pres(iparam);
        m_shp_pres.dx(iparam) = m_enr_derxy_pres(0,iparam);
        m_shp_pres.dy(iparam) = m_enr_derxy_pres(1,iparam);
        m_shp_pres.dz(iparam) = m_enr_derxy_pres(2,iparam);
        m_shp_pres.dxdx(iparam) = m_enr_derxy2_pres(0,iparam);
        m_shp_pres.dxdy(iparam) = m_enr_derxy2_pres(3,iparam);
        m_shp_pres.dxdz(iparam) = m_enr_derxy2_pres(4,iparam);
        m_shp_pres.dydx(iparam) = m_shp_pres.dxdy(iparam);
        m_shp_pres.dydy(iparam) = m_enr_derxy2_pres(1,iparam);
        m_shp_pres.dydz(iparam) = m_enr_derxy2_pres(5,iparam);
        m_shp_pres.dzdx(iparam) = m_shp_pres.dxdz(iparam);
        m_shp_pres.dzdy(iparam) = m_shp_pres.dydz(iparam);
        m_shp_pres.dzdz(iparam) = m_enr_derxy2_pres(2,iparam);
      }
    }
    else // not xfem_assembly i.e. standard assembly
    {
      // fill approximation functions for standard FEM
      // remark: numparamvelx == numnode, for standard FEM
      for (size_t iparam = 0; iparam < ms_numnode; ++iparam)
      {
        m_shp_vel.d0(iparam) = m_funct(iparam);
        m_shp_vel.dx(iparam) = m_derxy(0,iparam);
        m_shp_vel.dy(iparam) = m_derxy(1,iparam);
        m_shp_vel.dz(iparam) = m_derxy(2,iparam);
        m_shp_vel.dxdx(iparam) = m_derxy2(0,iparam);
        m_shp_vel.dxdy(iparam) = m_derxy2(3,iparam);
        m_shp_vel.dxdz(iparam) = m_derxy2(4,iparam);
        m_shp_vel.dydx(iparam) = m_shp_vel.dxdy(iparam);
        m_shp_vel.dydy(iparam) = m_derxy2(1,iparam);
        m_shp_vel.dydz(iparam) = m_derxy2(5,iparam);
        m_shp_vel.dzdx(iparam) = m_shp_vel.dxdz(iparam);
        m_shp_vel.dzdy(iparam) = m_shp_vel.dydz(iparam);
        m_shp_vel.dzdz(iparam) = m_derxy2(2,iparam);
      }

      for (size_t iparam = 0; iparam < ms_numnode; ++iparam)
      {
        m_shp_pres.d0(iparam) = m_funct(iparam);
        m_shp_pres.dx(iparam) = m_derxy(0,iparam);
        m_shp_pres.dy(iparam) = m_derxy(1,iparam);
        m_shp_pres.dz(iparam) = m_derxy(2,iparam);
        m_shp_pres.dxdx(iparam) = m_derxy2(0,iparam);
        m_shp_pres.dxdy(iparam) = m_derxy2(3,iparam);
        m_shp_pres.dxdz(iparam) = m_derxy2(4,iparam);
        m_shp_pres.dydx(iparam) = m_shp_pres.dxdy(iparam);
        m_shp_pres.dydy(iparam) = m_derxy2(1,iparam);
        m_shp_pres.dydz(iparam) = m_derxy2(5,iparam);
        m_shp_pres.dzdx(iparam) = m_shp_pres.dxdz(iparam);
        m_shp_pres.dzdy(iparam) = m_shp_pres.dydz(iparam);
        m_shp_pres.dzdz(iparam) = m_derxy2(2,iparam);
      }
    }

    //-------------------------------------------------
    // evaluate enriched shape functions w.r.t. slave
    //-------------------------------------------------

    // (enriched) shape functions = approximation functions (P = N * \Psi)
    static XFEM::ApproxFunc<2,s_shpVecSize> s_shp_vel;
    static XFEM::ApproxFunc<2,s_shpVecSize> s_shp_pres;

    if (S_ASSTYPE == XFEM::xfem_assembly)
    {
      // temporary arrays holding enriched shape functions (N * \Psi)
      static LINALG::Matrix<s_shpVecSize,1> s_enr_funct_vel;
      s_enr_funct_vel.Clear();
      static LINALG::Matrix<3,s_shpVecSize> s_enr_derxy_vel;
      s_enr_derxy_vel.Clear();
      static LINALG::Matrix<6,s_shpVecSize> s_enr_derxy2_vel;
      s_enr_derxy2_vel.Clear();
      static LINALG::Matrix<s_shpVecSize,1> s_enr_funct_pres;
      s_enr_funct_pres.Clear();
      static LINALG::Matrix<3,s_shpVecSize> s_enr_derxy_pres;
      s_enr_derxy_pres.Clear();
      static LINALG::Matrix<6,s_shpVecSize> s_enr_derxy2_pres;
      s_enr_derxy2_pres.Clear();

      // shape functions and derivatives for nodal parameters (dofs)
      slave_enrvals.ComputeModifiedEnrichedNodalShapefunction(
          XFEM::PHYSICS::Velx,
          s_funct,
          s_derxy,
          s_derxy2,
          s_enr_funct_vel,
          s_enr_derxy_vel,
          s_enr_derxy2_vel);

      // fill approximation functions for XFEM
      for (size_t iparam = 0; iparam != s_numparamvelx; ++iparam)
      {
        s_shp_vel.d0(iparam) = s_enr_funct_vel(iparam);
        s_shp_vel.dx(iparam) = s_enr_derxy_vel(0,iparam);
        s_shp_vel.dy(iparam) = s_enr_derxy_vel(1,iparam);
        s_shp_vel.dz(iparam) = s_enr_derxy_vel(2,iparam);
        s_shp_vel.dxdx(iparam) = s_enr_derxy2_vel(0,iparam);
        s_shp_vel.dxdy(iparam) = s_enr_derxy2_vel(3,iparam);
        s_shp_vel.dxdz(iparam) = s_enr_derxy2_vel(4,iparam);
        s_shp_vel.dydx(iparam) = s_shp_vel.dxdy(iparam);
        s_shp_vel.dydy(iparam) = s_enr_derxy2_vel(1,iparam);
        s_shp_vel.dydz(iparam) = s_enr_derxy2_vel(5,iparam);
        s_shp_vel.dzdx(iparam) = s_shp_vel.dxdz(iparam);
        s_shp_vel.dzdy(iparam) = s_shp_vel.dydz(iparam);
        s_shp_vel.dzdz(iparam) = s_enr_derxy2_vel(2,iparam);
      }

      // shape functions and derivatives for nodal parameters (dofs)
      slave_enrvals.ComputeModifiedEnrichedNodalShapefunction(
          XFEM::PHYSICS::Pres,
          s_funct,
          s_derxy,
          s_derxy2,
          s_enr_funct_pres,
          s_enr_derxy_pres,
          s_enr_derxy2_pres);

      // fill approximation functions for XFEM
      for (size_t iparam = 0; iparam != s_numparampres; ++iparam)
      {
        s_shp_pres.d0(iparam) = s_enr_funct_pres(iparam);
        s_shp_pres.dx(iparam) = s_enr_derxy_pres(0,iparam);
        s_shp_pres.dy(iparam) = s_enr_derxy_pres(1,iparam);
        s_shp_pres.dz(iparam) = s_enr_derxy_pres(2,iparam);
        s_shp_pres.dxdx(iparam) = s_enr_derxy2_pres(0,iparam);
        s_shp_pres.dxdy(iparam) = s_enr_derxy2_pres(3,iparam);
        s_shp_pres.dxdz(iparam) = s_enr_derxy2_pres(4,iparam);
        s_shp_pres.dydx(iparam) = s_shp_pres.dxdy(iparam);
        s_shp_pres.dydy(iparam) = s_enr_derxy2_pres(1,iparam);
        s_shp_pres.dydz(iparam) = s_enr_derxy2_pres(5,iparam);
        s_shp_pres.dzdx(iparam) = s_shp_pres.dxdz(iparam);
        s_shp_pres.dzdy(iparam) = s_shp_pres.dydz(iparam);
        s_shp_pres.dzdz(iparam) = s_enr_derxy2_pres(2,iparam);
      }
    }
    else // not xfem_assembly i.e. standard assembly
    {
      // fill approximation functions for standard FEM
      // remark: numparamvelx == numnode, for standard FEM
      for (size_t iparam = 0; iparam < ms_numnode; ++iparam)
      {
        s_shp_vel.d0(iparam) = s_funct(iparam);
        s_shp_vel.dx(iparam) = s_derxy(0,iparam);
        s_shp_vel.dy(iparam) = s_derxy(1,iparam);
        s_shp_vel.dz(iparam) = s_derxy(2,iparam);
        s_shp_vel.dxdx(iparam) = s_derxy2(0,iparam);
        s_shp_vel.dxdy(iparam) = s_derxy2(3,iparam);
        s_shp_vel.dxdz(iparam) = s_derxy2(4,iparam);
        s_shp_vel.dydx(iparam) = s_shp_vel.dxdy(iparam);
        s_shp_vel.dydy(iparam) = s_derxy2(1,iparam);
        s_shp_vel.dydz(iparam) = s_derxy2(5,iparam);
        s_shp_vel.dzdx(iparam) = s_shp_vel.dxdz(iparam);
        s_shp_vel.dzdy(iparam) = s_shp_vel.dydz(iparam);
        s_shp_vel.dzdz(iparam) = s_derxy2(2,iparam);
      }

      for (size_t iparam = 0; iparam < ms_numnode; ++iparam)
      {
        s_shp_pres.d0(iparam) = s_funct(iparam);
        s_shp_pres.dx(iparam) = s_derxy(0,iparam);
        s_shp_pres.dy(iparam) = s_derxy(1,iparam);
        s_shp_pres.dz(iparam) = s_derxy(2,iparam);
        s_shp_pres.dxdx(iparam) = s_derxy2(0,iparam);
        s_shp_pres.dxdy(iparam) = s_derxy2(3,iparam);
        s_shp_pres.dxdz(iparam) = s_derxy2(4,iparam);
        s_shp_pres.dydx(iparam) = s_shp_pres.dxdy(iparam);
        s_shp_pres.dydy(iparam) = s_derxy2(1,iparam);
        s_shp_pres.dydz(iparam) = s_derxy2(5,iparam);
        s_shp_pres.dzdx(iparam) = s_shp_pres.dxdz(iparam);
        s_shp_pres.dzdy(iparam) = s_shp_pres.dydz(iparam);
        s_shp_pres.dzdz(iparam) = s_derxy2(2,iparam);
      }
    }

    //--------------------------------------------------
    // get velocity and pressure at integration point
    //--------------------------------------------------
    // since functions are continuous at element boundaries, one can
    // evaluate the velocity using the master element, the slave element or the face
    // TODO: what happens for touched elements: Averaging of both sides at interface
    const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPoint(evelnp_m, m_shp_vel.d0, m_numparamvelx);
    const double gppresnp = XFEM::interpolateScalarFieldToIntPoint(eprenp_m, m_shp_pres.d0, m_numparampres);

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
        m_vderxy(isd,0) += evelnp_m(isd,iparam) * m_shp_vel.dx(iparam);
        m_vderxy(isd,1) += evelnp_m(isd,iparam) * m_shp_vel.dy(iparam);
        m_vderxy(isd,2) += evelnp_m(isd,iparam) * m_shp_vel.dz(iparam);
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
        s_vderxy(isd,0) += evelnp_s(isd,iparam) * s_shp_vel.dx(iparam);
        s_vderxy(isd,1) += evelnp_s(isd,iparam) * s_shp_vel.dy(iparam);
        s_vderxy(isd,2) += evelnp_s(isd,iparam) * s_shp_vel.dz(iparam);
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

    //----------------------------------------------
    // get gradient jump of velocity at interface
    //----------------------------------------------
    LINALG::Matrix<nsd,nsd> grad_vel_jump(true);
    grad_vel_jump.Update(1.0, s_vderxy, -1.0, m_vderxy);

    //----------------------------------------------
    // projection of velocity in normal direction of face
    //----------------------------------------------
    const double normal_vel = fabs(gpvelnp.Dot(normal));

    //-----------------------------------------------------
    // get second velocity derivatives at integration point
    //-----------------------------------------------------
    // master and slave
    static LINALG::Matrix<3,6> m_vderxy2;
    m_vderxy2.Clear();
    static LINALG::Matrix<3,6> s_vderxy2;
    s_vderxy2.Clear();

    if (higher_order_ele)
    {
      for (size_t iparam = 0; iparam < m_numparamvelx; ++iparam)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
        {
          m_vderxy2(isd,0) += evelnp_m(isd,iparam) * m_shp_vel.dxdx(iparam);
          m_vderxy2(isd,1) += evelnp_m(isd,iparam) * m_shp_vel.dydy(iparam);
          m_vderxy2(isd,2) += evelnp_m(isd,iparam) * m_shp_vel.dzdz(iparam);
          m_vderxy2(isd,3) += evelnp_m(isd,iparam) * m_shp_vel.dxdy(iparam);
          m_vderxy2(isd,4) += evelnp_m(isd,iparam) * m_shp_vel.dxdz(iparam);
          m_vderxy2(isd,5) += evelnp_m(isd,iparam) * m_shp_vel.dydz(iparam);
        }
      }

      for (size_t iparam = 0; iparam < s_numparamvelx; ++iparam)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
        {
          s_vderxy2(isd,0) += evelnp_s(isd,iparam) * s_shp_vel.dxdx(iparam);
          s_vderxy2(isd,1) += evelnp_s(isd,iparam) * s_shp_vel.dydy(iparam);
          s_vderxy2(isd,2) += evelnp_s(isd,iparam) * s_shp_vel.dzdz(iparam);
          s_vderxy2(isd,3) += evelnp_s(isd,iparam) * s_shp_vel.dxdy(iparam);
          s_vderxy2(isd,4) += evelnp_s(isd,iparam) * s_shp_vel.dxdz(iparam);
          s_vderxy2(isd,5) += evelnp_s(isd,iparam) * s_shp_vel.dydz(iparam);
        }
      }
    }

    //----------------------------------------------
    // get second derivative jump of velocity at interface
    //----------------------------------------------
    LINALG::Matrix<3,6> der2_vel_jump(true);
    if (higher_order_ele)
      der2_vel_jump.Update(1.0, s_vderxy2, -1.0, m_vderxy2);

    //--------------------------------------------
    // get pressure gradients at integration point
    //--------------------------------------------
    // master
    static LINALG::Matrix<nsd,1> m_gradp;
    m_gradp.Clear();
    // gradp = enr_derxy(i,j)*eprenp(j);
    for (size_t iparam = 0; iparam != m_numparampres; ++iparam)
    {
      m_gradp(0) += eprenp_m(iparam) * m_shp_pres.dx(iparam);
      m_gradp(1) += eprenp_m(iparam) * m_shp_pres.dy(iparam);
      m_gradp(2) += eprenp_m(iparam) * m_shp_pres.dz(iparam);
    }
    // slave
    static LINALG::Matrix<nsd,1> s_gradp;
    s_gradp.Clear();
    // gradp = enr_derxy(i,j)*eprenp(j);
    for (size_t iparam = 0; iparam != s_numparampres; ++iparam)
    {
      s_gradp(0) += eprenp_s(iparam) * s_shp_pres.dx(iparam);
      s_gradp(1) += eprenp_s(iparam) * s_shp_pres.dy(iparam);
      s_gradp(2) += eprenp_s(iparam) * s_shp_pres.dz(iparam);
    }

    //--------------------------------------------
    // get pressure jump at integration point
    //--------------------------------------------
    // grad(p_slave) - grad(p_master)
    LINALG::Matrix<nsd,1> grad_pres_jump(true);
    grad_pres_jump.Update(1.0, s_gradp, -1.0, m_gradp);

    //-----------------------------------------------------
    // get second pressure derivatives at integration point
    //-----------------------------------------------------
    // master
    static LINALG::Matrix<6,1> m_der2p;
    m_der2p.Clear();
    // slave
    static LINALG::Matrix<6,1> s_der2p;
    s_der2p.Clear();
    if (higher_order_ele)
    {
      for (size_t iparam = 0; iparam != m_numparampres; ++iparam)
      {
        m_der2p(0) += eprenp_m(iparam) * m_shp_pres.dxdx(iparam);
        m_der2p(1) += eprenp_m(iparam) * m_shp_pres.dydy(iparam);
        m_der2p(2) += eprenp_m(iparam) * m_shp_pres.dzdz(iparam);
        m_der2p(3) += eprenp_m(iparam) * m_shp_pres.dxdy(iparam);
        m_der2p(4) += eprenp_m(iparam) * m_shp_pres.dxdz(iparam);
        m_der2p(5) += eprenp_m(iparam) * m_shp_pres.dydz(iparam);
      }

      for (size_t iparam = 0; iparam != s_numparampres; ++iparam)
      {
        s_der2p(0) += eprenp_s(iparam) * s_shp_pres.dxdx(iparam);
        s_der2p(1) += eprenp_s(iparam) * s_shp_pres.dydy(iparam);
        s_der2p(2) += eprenp_s(iparam) * s_shp_pres.dzdz(iparam);
        s_der2p(3) += eprenp_s(iparam) * s_shp_pres.dxdy(iparam);
        s_der2p(4) += eprenp_s(iparam) * s_shp_pres.dxdz(iparam);
        s_der2p(5) += eprenp_s(iparam) * s_shp_pres.dydz(iparam);
      }
    }

    //-----------------------------------------------------
    // get jump of second pressure derivatives at integration point
    //-----------------------------------------------------
    // D^2(p_slave) - D^2(p_master)
    LINALG::Matrix<6,1> der2_pres_jump(true);
    if (higher_order_ele)
      der2_pres_jump.Update(1.0, s_der2p, -1.0, m_der2p);

    //------------------------------------------
    // compute stabilization parameters (3 taus)
    //------------------------------------------

    double tau_conv = 0.0;   // convection-stabilization parameter
    double tau_div  = 0.0;   // divergence-stabilization parameter
    double tau_p    = 0.0;   // pressure-stabilization parameter

    double tau_ghost_first_order = 0.0;   // ghost-penalty stabilization parameter due to Nitsche's method in the XFEM
    double tau_ghost_second_order = 0.0;  // ghost-penalty stabilization parameter due to Nitsche's method in the XFEM
    double tau_pres_second_order = 0.0;   // ghost-penalty stabilization parameter due to Nitsche's method in the XFEM

    COMBUST::UTILS::computeStabilizationParamsEdgeBased(
            dynvisc, densaf, vel_norm, normal_vel, hk,
            tau_def, tau_conv, tau_div, tau_p, instationary, timefac,
            add_ghost_penalties, tau_ghost_first_order, tau_ghost_second_order, tau_pres_second_order);

    //----------------------------------------------
    // build stiffness matrix stabilization entries
    //----------------------------------------------
    // that is master-master, master-slave, slave-master and slave-slave block
    COMBUST::BuildStandardFaceIntegrals<M_DISTYPE, S_DISTYPE, M_ASSTYPE, S_ASSTYPE, NUMDOF, m_shpVecSize, s_shpVecSize>(
            assembler, m_shp_vel, s_shp_vel, m_shp_pres, s_shp_pres,
            fac, timefac, rhstimefac, genalpha, instationary, densaf, dynvisc,
            gpvelnp, m_vderxy, s_vderxy,
            gppresnp, s_gradp, s_gradp,
            grad_pres_jump, vdiv_jump, grad_vel_jump,
            der2_vel_jump, der2_pres_jump,
            pres_stab, conv_stream_stab, conv_cross_stab, conti_stab,
            tau_conv, tau_div, tau_p, add_ghost_penalties, tau_ghost_first_order, tau_ghost_second_order, tau_pres_second_order,
            higher_order_ele);
  } // end: loop all gaussian points

  //----------------------------------------
  // finish assembly
  //----------------------------------------
  // assemble master-master, master-slave, slave-master and slave-slave block
  // into patch matrix
  assembler.Finish();
  // assembler.PrintMat();
  // assembler.PrintRhs();

  return;
}

} // end namespace COMBUST

#endif
