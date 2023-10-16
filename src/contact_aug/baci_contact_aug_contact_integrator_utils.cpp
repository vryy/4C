/*---------------------------------------------------------------------*/
/*! \file
\brief Utility methods for the contact integration.

\level 2

*/
/*---------------------------------------------------------------------*/

#include "baci_contact_aug_contact_integrator_utils.H"

#include "baci_contact_aug_element_utils.H"
#include "baci_contact_aug_integrator_policy.H"
#include "baci_contact_aug_projector.H"
#include "baci_contact_aug_utils.H"
#include "baci_contact_integrator.H"
#include "baci_contact_node.H"
#include "baci_io_pstream.H"
#include "baci_mortar_element.H"

#include <Teuchos_TimeMonitor.hpp>

// #define DEBUG_FIND_FEASIBLE_MASTER_ELEMENT

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::INTEGRATOR::FindFeasibleMasterElements(MORTAR::MortarElement& sele,
    const std::vector<MORTAR::MortarElement*>& meles, bool boundary_ele, CoIntegrator& wrapper,
    UniqueProjInfoPair& projInfo)
{
  TEUCHOS_FUNC_TIME_MONITOR(AUG::CONTACT_FUNC_NAME);

  const CORE::FE::CellType slavetype = sele.Shape();

  const unsigned msize = meles.size();
  const unsigned numGP = wrapper.nGP();
  const unsigned reserve_size = (msize ? numGP / msize : 0);

  const unsigned probdim = wrapper.Dim();

  CORE::GEN::reset(msize, projInfo);

  std::vector<double> found_alpha(msize, 0.0);
  std::vector<CORE::LINALG::Matrix<2, 1>> found_mxi(msize, CORE::LINALG::Matrix<2, 1>(true));
  std::vector<bool> is_on_meles(msize, false);

  std::vector<int> unique_ids(msize, -1);
  unsigned num_unique_ids = 0;
  double uniqueProjAlpha = std::numeric_limits<double>::max();

  double proj_tol = 0.0;

  for (unsigned gp = 0; gp < numGP; ++gp)
  {
    // get Gauss point in slave element coordinates
    const double eta[2] = {wrapper.Coordinate(gp, 0), wrapper.Coordinate(gp, 1)};
    const double sxi[2] = {eta[0], eta[1]};

    // reset is_on_meles flag
    std::fill(is_on_meles.begin(), is_on_meles.end(), false);

    // --> find master elements with an unique projection
    for (unsigned m = 0; m < msize; ++m)
    {
      double* mxi = found_mxi[m].A();
      double& projalpha = found_alpha[m];

      const CORE::FE::CellType mastertype = meles[m]->Shape();
      // project Gauss point onto master element
#ifdef DEBUG_FIND_FEASIBLE_MASTER_ELEMENT
      CONTACT::AUG::ProjectorBase* proj_ptr = nullptr;
      if (sele.Id() == -1 or sele.Id() == -1)
        proj_ptr = CONTACT::AUG::ProjectorBase::Get(probdim, slavetype, mastertype, true);
      else
        proj_ptr = CONTACT::AUG::ProjectorBase::Get(probdim, slavetype, mastertype);
#else
      CONTACT::AUG::ProjectorBase* proj_ptr =
          CONTACT::AUG::ProjectorBase::Get(probdim, slavetype, mastertype);
#endif

      CONTACT::AUG::ProjectorBase& proj = *proj_ptr;
      const bool conv = proj(sele, sxi, *meles[m], mxi, projalpha);

      // --- DEBUGGING --------------------------------------------------------
#ifdef DEBUG_FIND_FEASIBLE_MASTER_ELEMENT
      DRT::Node** snodes = sele.Nodes();
      if ((snodes[0]->Id() == 248 or snodes[1]->Id() == 251))
      //      if ( sele.Id() == 1261 )
      {
        std::cout << "Slave element #" << sele.Id() << std::endl;
        std::cout << "GP #" << gp << " | mele " << meles[m]->Id() << " [" << meles[m] << "]"
                  << std::setprecision(16) << ": sxi( " << sxi[0] << ") --> mxi( " << mxi[0]
                  << " ), alpha = " << projalpha << " projection "
                  << (conv ? "succeeded" : "failed") << std::endl;
      }
#endif

      // check GP projection
      proj_tol = 10.0 * proj.getRelativeSolutionTolerance();
      is_on_meles[m] = WithinBounds(mxi, mastertype, proj_tol);

      if (is_on_meles[m] and not conv)
        dserror(
            "The gp has a feasible projection, but the local Newton scheme "
            "did not converge!");
    }  // mele loop


    // reset
    num_unique_ids = 0;
    uniqueProjAlpha = std::numeric_limits<double>::max();

    for (unsigned m = 0; m < msize; ++m)
    {
      if (not is_on_meles[m]) continue;

      // check if the found projection point lies on a edge or a corner point
      if (std::abs(found_alpha[m] - uniqueProjAlpha) < proj_tol)
      {
        unique_ids[num_unique_ids++] = m;
      }
      // found a second master element with a feasible projection, but a smaller,
      // not necessarily shorter, distance value
      else if (found_alpha[m] < uniqueProjAlpha)
      {
        unique_ids.front() = m;
        num_unique_ids = 1;
        uniqueProjAlpha = found_alpha[m];
      }
    }

    if (num_unique_ids > 0)
    {
      // If the gp projects onto an edge or node the corresponding gp-weight
      // must be scaled accordingly.
      const double gp_scaling = 1.0 / num_unique_ids;
      for (unsigned i = 0; i < num_unique_ids; ++i)
      {
        const unsigned unique_id = unique_ids[i];

        MORTAR::MortarElement* mele = meles[unique_id];
        if (projInfo.find(mele) == projInfo.end())
        {
          projInfo[mele] = UniqueProjInfo(std::ceil(reserve_size));
        }

        projInfo[mele].Insert(gp, found_alpha[unique_id], found_mxi[unique_id].A(), gp_scaling);
      }
    }
    else if (not boundary_ele)
      std::cout << "*** warning *** Non-boundary element has non-projectable "
                   "GP \n";
  }

  // --- DEBUGGING --------------------------------------------------------
#ifdef DEBUG_FIND_FEASIBLE_MASTER_ELEMENT
  DRT::Node** snodes = sele.Nodes();
  const std::vector<int> desired_nids = {247, 248, 251};
  unsigned count = 0;

  for (unsigned i = 0; i < static_cast<unsigned>(sele.NumNode()); ++i)
  {
    const int sid = snodes[i]->Id();
    if (std::find(desired_nids.begin(), desired_nids.end(), sid) != desired_nids.end()) ++count;
  }

  if (count == 2)
  //  if ( ( snodes[0]->Id() == 7731 and snodes[1]->Id() == 7733 ) )
  //  if ( sele.Id() == 1260 or sele.Id() == 1261 )
  {
    for (auto& proj : projInfo)
    {
      std::cout << "\n---------------------------------------\n";
      std::cout << "sele #" << sele.Id() << ":\n";
      sele.Print(std::cout);
      std::cout << "\nSlave-nodes: ";
      for (unsigned i = 0; i < static_cast<unsigned>(sele.NumNode()); ++i)
      {
        std::cout << "#" << snodes[i]->Id() << " --> " << snodes[i]->X()[0] << ", "
                  << snodes[i]->X()[1] << ", " << snodes[i]->X()[2] << "\n";
      }
      std::cout << "\n";
      std::cout << "mele " << proj.first->Id() << "[" << proj.first << "]"
                << "\n";
      DRT::Node** mnodes = proj.first->Nodes();
      std::cout << "\nMaster-nodes: ";
      for (unsigned i = 0; i < static_cast<unsigned>(proj.first->NumNode()); ++i)
      {
        std::cout << "#" << mnodes[i]->Id() << " --> " << mnodes[i]->X()[0] << ", "
                  << mnodes[i]->X()[1] << ", " << mnodes[i]->X()[2] << "\n";
      }
      proj.second.Print(std::cout);
    }
  }
#endif


  return (projInfo.size() > 0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::INTEGRATOR::WithinBounds(
    const double* mxi, const CORE::FE::CellType type, const double tol)
{
  switch (type)
  {
    case CORE::FE::CellType::quad4:
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
    case CORE::FE::CellType::nurbs4:
    case CORE::FE::CellType::nurbs9:
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
        return false;

      break;
    }
    case CORE::FE::CellType::tri3:
    case CORE::FE::CellType::tri6:
    {
      if (mxi[0] < -tol or mxi[1] < -tol or mxi[0] > 1.0 + tol or mxi[1] > 1.0 + tol or
          mxi[0] + mxi[1] > 1.0 + 2 * tol)
        return false;

      break;
    }
    case CORE::FE::CellType::line2:
    case CORE::FE::CellType::line3:
    case CORE::FE::CellType::nurbs2:
    case CORE::FE::CellType::nurbs3:
    {
      if ((mxi[0] < -1.0 - tol) or (mxi[0] > 1.0 + tol)) return false;

      break;
    }
    default:
    {
      dserror("Unsupported element type %s!", DRT::DistypeToString(type).c_str());
      exit(EXIT_FAILURE);
    }
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::INTEGRATOR::BuildAveragedNormalAtSlaveNode(
    std::vector<ElementNormal>& adj_ele_normals, MORTAR::MortarNode& slavenode)
{
  // access the averaged nodal normal
  CORE::LINALG::Matrix<3, 1> avg_nodal_normal(slavenode.MoData().n(), true);

  // reset normal
  std::fill(avg_nodal_normal.A(), avg_nodal_normal.A() + 3, 0.0);

  const unsigned num_adj_eles = slavenode.NumElement();
  DRT::Element** adj_eles = slavenode.Elements();

  adj_ele_normals.clear();
  adj_ele_normals.resize(num_adj_eles, ElementNormal());

  for (unsigned e = 0; e < num_adj_eles; ++e)
  {
    MORTAR::MortarElement& adj_ele = dynamic_cast<MORTAR::MortarElement&>(*(adj_eles[e]));

    adj_ele_normals[e].ele_ = &adj_ele;

    // get parametric coordinates of the current slave node in the adjacent element
    const int ele_node_lid = adj_ele.GetLocalNodeId(slavenode.Id());
    double* xi = adj_ele_normals[e].xi_;
    adj_ele.LocalCoordinatesOfNode(ele_node_lid, xi);

    // evaluate the convective base vectors
    CORE::LINALG::Matrix<3, 2> tau;
    adj_ele.Metrics(xi, &tau(0, 0), &tau(0, 1));

    // evaluate the unit normal of the adjacent element at the current slave node
    CORE::LINALG::Matrix<3, 1>& unit_normal = adj_ele_normals[e].unit_n_;
    double& length_n_inv = adj_ele_normals[e].length_n_inv_;
    length_n_inv = UnitSlaveElementNormal(adj_ele, tau, unit_normal);

    // sum up all nodal unit element normals of the adjacent elements
    avg_nodal_normal.Update(1.0, unit_normal, 1.0);
  }

  // average the result
  const double smooth_nodal_normal_length = avg_nodal_normal.Norm2();

  if (smooth_nodal_normal_length == 0.0)
    dserror("Sum of slave unit normals at node %d has a length of zero!", slavenode.Id());

  avg_nodal_normal.Scale(1.0 / smooth_nodal_normal_length);

  return smooth_nodal_normal_length;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::INTEGRATOR::UnitSlaveElementNormal(const MORTAR::MortarElement& sele,
    const CORE::LINALG::Matrix<3, 2>& tau, CORE::LINALG::Matrix<3, 1>& unit_normal)
{
  const CORE::FE::CellType slavetype = sele.Shape();
  switch (slavetype)
  {
    case CORE::FE::CellType::line2:
    {
      const CONTACT::AUG::BaseSlaveIntPolicy<2, CORE::FE::CellType::line2> slave_policy;

      CORE::LINALG::Matrix<2, 1> unit_n(unit_normal.A(), true);
      return slave_policy.UnitSlaveElementNormal(sele, tau, unit_n);
    }
    case CORE::FE::CellType::nurbs2:
    {
      const CONTACT::AUG::BaseSlaveIntPolicy<2, CORE::FE::CellType::nurbs2> slave_policy;

      CORE::LINALG::Matrix<2, 1> unit_n(unit_normal.A(), true);
      return slave_policy.UnitSlaveElementNormal(sele, tau, unit_n);
    }
    case CORE::FE::CellType::nurbs3:
    {
      const CONTACT::AUG::BaseSlaveIntPolicy<2, CORE::FE::CellType::nurbs3> slave_policy;

      CORE::LINALG::Matrix<2, 1> unit_n(unit_normal.A(), true);
      return slave_policy.UnitSlaveElementNormal(sele, tau, unit_n);
    }
    case CORE::FE::CellType::quad4:
    {
      const CONTACT::AUG::BaseSlaveIntPolicy<3, CORE::FE::CellType::quad4> slave_policy;

      return slave_policy.UnitSlaveElementNormal(sele, tau, unit_normal);
    }
    case CORE::FE::CellType::tri3:
    {
      const CONTACT::AUG::BaseSlaveIntPolicy<3, CORE::FE::CellType::tri3> slave_policy;

      return slave_policy.UnitSlaveElementNormal(sele, tau, unit_normal);
    }
    case CORE::FE::CellType::nurbs4:
    {
      const CONTACT::AUG::BaseSlaveIntPolicy<3, CORE::FE::CellType::nurbs4> slave_policy;

      return slave_policy.UnitSlaveElementNormal(sele, tau, unit_normal);
    }
    case CORE::FE::CellType::nurbs9:
    {
      const CONTACT::AUG::BaseSlaveIntPolicy<3, CORE::FE::CellType::nurbs9> slave_policy;

      return slave_policy.UnitSlaveElementNormal(sele, tau, unit_normal);
    }
    default:
    {
      dserror("Unsupported slave element type! (enum = %d|\"%s\")", slavetype,
          DRT::DistypeToString(slavetype).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::INTEGRATOR::Deriv1st_AveragedSlaveNormal(CONTACT::CoNode& cnode,
    const std::vector<ElementNormal>& adj_ele_normals, const double avg_normal_length,
    Deriv1stVecMap& d_nodal_avg_normal)
{
  Deriv1stVecMap& d_avg_unit_normal = cnode.AugData().GetDeriv1st_N();
  CORE::GEN::reset(3, cnode.GetLinsize(), d_avg_unit_normal);
  CORE::GEN::reset(3, cnode.GetLinsize(), d_nodal_avg_normal);

  Deriv1stVecMap d_non_unit_normal;

  CORE::FE::CellType eletype = CORE::FE::CellType::dis_none;
  for (const ElementNormal& adj_ele_n : adj_ele_normals)
  {
    MORTAR::MortarElement& mo_ele = *adj_ele_n.ele_;

    eletype = mo_ele.Shape();

    const double* xi = adj_ele_n.xi_;

    Deriv1st_NonUnitSlaveNormal(xi, mo_ele, d_non_unit_normal);

    const CORE::LINALG::Matrix<3, 1>& adj_ele_unit_n = adj_ele_n.unit_n_;
    const double adj_ele_length_n_inv = adj_ele_n.length_n_inv_;

    Deriv1st_UnitSlaveNormal(eletype, adj_ele_unit_n, adj_ele_length_n_inv, d_non_unit_normal,
        d_nodal_avg_normal, false);
  }

  const CORE::LINALG::Matrix<3, 1> avg_unit_normal(cnode.MoData().n(), true);
  Deriv1st_UnitSlaveNormal(eletype, avg_unit_normal, 1.0 / avg_normal_length, d_nodal_avg_normal,
      d_avg_unit_normal, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::INTEGRATOR::Deriv1st_UnitSlaveNormal(const CORE::FE::CellType slavetype,
    const CORE::LINALG::Matrix<3, 1>& unit_normal, const double length_n_inv,
    const Deriv1stVecMap& d_non_unit_normal, Deriv1stVecMap& d_unit_normal, const bool reset)
{
  switch (slavetype)
  {
    case CORE::FE::CellType::line2:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::line2>::dim;
      const unsigned probdim = eledim + 1;

      const CORE::LINALG::Matrix<probdim, 1> unit_normal_red(unit_normal.A(), true);

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::line2> slave_policy;
      slave_policy.Deriv1st_UnitSlaveElementNormal(
          unit_normal_red, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    case CORE::FE::CellType::nurbs2:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::nurbs2>::dim;
      const unsigned probdim = eledim + 1;

      const CORE::LINALG::Matrix<probdim, 1> unit_normal_red(unit_normal.A(), true);

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::nurbs2> slave_policy;
      slave_policy.Deriv1st_UnitSlaveElementNormal(
          unit_normal_red, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    case CORE::FE::CellType::nurbs3:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::nurbs3>::dim;
      const unsigned probdim = eledim + 1;

      const CORE::LINALG::Matrix<probdim, 1> unit_normal_red(unit_normal.A(), true);

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::nurbs3> slave_policy;
      slave_policy.Deriv1st_UnitSlaveElementNormal(
          unit_normal_red, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    case CORE::FE::CellType::quad4:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::quad4>::dim;
      const unsigned probdim = eledim + 1;

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::quad4> slave_policy;
      slave_policy.Deriv1st_UnitSlaveElementNormal(
          unit_normal, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    case CORE::FE::CellType::tri3:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::tri3>::dim;
      const unsigned probdim = eledim + 1;

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::tri3> slave_policy;
      slave_policy.Deriv1st_UnitSlaveElementNormal(
          unit_normal, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    case CORE::FE::CellType::nurbs4:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::nurbs4>::dim;
      const unsigned probdim = eledim + 1;

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::nurbs4> slave_policy;
      slave_policy.Deriv1st_UnitSlaveElementNormal(
          unit_normal, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    case CORE::FE::CellType::nurbs9:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::nurbs9>::dim;
      const unsigned probdim = eledim + 1;

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::nurbs9> slave_policy;
      slave_policy.Deriv1st_UnitSlaveElementNormal(
          unit_normal, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    default:
    {
      dserror("Unsupported slave element type! (enum = %d|\"%s\")", slavetype,
          DRT::DistypeToString(slavetype).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::INTEGRATOR::Deriv1st_NonUnitSlaveNormal(
    const double* xi, MORTAR::MortarElement& sele, Deriv1stVecMap& d_non_unit_normal)
{
  const CORE::FE::CellType slavetype = sele.Shape();
  switch (slavetype)
  {
    case CORE::FE::CellType::line2:
    {
      Deriv1st_NonUnitSlaveNormal<CORE::FE::CellType::line2>(sele, xi, d_non_unit_normal);
      break;
    }
    case CORE::FE::CellType::nurbs2:
    {
      Deriv1st_NonUnitSlaveNormal<CORE::FE::CellType::nurbs2>(sele, xi, d_non_unit_normal);
      break;
    }
    case CORE::FE::CellType::nurbs3:
    {
      Deriv1st_NonUnitSlaveNormal<CORE::FE::CellType::nurbs3>(sele, xi, d_non_unit_normal);
      break;
    }
    case CORE::FE::CellType::quad4:
    {
      Deriv1st_NonUnitSlaveNormal<CORE::FE::CellType::quad4>(sele, xi, d_non_unit_normal);
      break;
    }
    case CORE::FE::CellType::tri3:
    {
      Deriv1st_NonUnitSlaveNormal<CORE::FE::CellType::tri3>(sele, xi, d_non_unit_normal);
      break;
    }
    case CORE::FE::CellType::nurbs4:
    {
      Deriv1st_NonUnitSlaveNormal<CORE::FE::CellType::nurbs4>(sele, xi, d_non_unit_normal);
      break;
    }
    case CORE::FE::CellType::nurbs9:
    {
      Deriv1st_NonUnitSlaveNormal<CORE::FE::CellType::nurbs9>(sele, xi, d_non_unit_normal);
      break;
    }
    default:
    {
      dserror("Unsupported slave element type! (enum = %d|\"%s\")", slavetype,
          DRT::DistypeToString(slavetype).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <CORE::FE::CellType slavetype>
void CONTACT::INTEGRATOR::Deriv1st_NonUnitSlaveNormal(
    MORTAR::MortarElement& sele, const double* xi, Deriv1stVecMap& d_non_unit_normal)
{
  const unsigned slavenumnode =
      CORE::DRT::UTILS::DisTypeToNumNodePerEle<slavetype>::numNodePerElement;
  const unsigned slavedim = CORE::DRT::UTILS::DisTypeToDim<slavetype>::dim;
  const unsigned probdim = slavedim + 1;

  CORE::LINALG::Matrix<probdim, slavenumnode, int> nodal_dofs;
  CONTACT::INTEGRATOR::GetElementNodalDofs(sele, nodal_dofs);

  // evaluate the convective base vectors
  CORE::LINALG::Matrix<3, 2> tau;
  sele.Metrics(xi, &tau(0, 0), &tau(0, 1));

  // derivatives with respect to xi^1 and xi^2
  CORE::LINALG::Matrix<slavedim, slavenumnode> deriv(true);

  // derivative of the convective base vectors with respect to the displacement
  const CORE::LINALG::Matrix<2, 1> xi_mat(xi, true);
  CONTACT::AUG::shape_function_deriv1<slavetype>(sele, xi_mat, deriv);

  const CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype> slave_policy;
  slave_policy.Deriv1st_NonUnitSlaveElementNormal(sele, nodal_dofs, deriv, tau, d_non_unit_normal);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::INTEGRATOR::Deriv2nd_AveragedSlaveNormal(CONTACT::CoNode& cnode,
    const std::vector<ElementNormal>& adj_ele_normals, const double avg_normal_length,
    const Deriv1stVecMap& d_nodal_avg_normal)
{
  Deriv2ndVecMap& dd_avg_unit_normal = cnode.AugData().GetDeriv2nd_N();
  CORE::GEN::reset(3, cnode.GetLinsize(), dd_avg_unit_normal);

  Deriv2ndVecMap dd_nodal_avg_normal;
  CORE::GEN::reset(3, cnode.GetLinsize(), dd_nodal_avg_normal);

  Deriv1stVecMap d_non_unit_normal;
  Deriv1stVecMap d_unit_normal;

  Deriv2ndVecMap dd_non_unit_normal(3, Deriv2ndMap());

  CORE::FE::CellType eletype = CORE::FE::CellType::dis_none;
  for (const ElementNormal& adj_ele_n : adj_ele_normals)
  {
    MORTAR::MortarElement& mo_ele = *adj_ele_n.ele_;

    eletype = mo_ele.Shape();

    const double* xi = adj_ele_n.xi_;

    Deriv1st_NonUnitSlaveNormal(xi, mo_ele, d_non_unit_normal);

    const CORE::LINALG::Matrix<3, 1>& adj_ele_unit_n = adj_ele_n.unit_n_;
    const double adj_ele_length_n_inv = adj_ele_n.length_n_inv_;

    Deriv1st_UnitSlaveNormal(
        eletype, adj_ele_unit_n, adj_ele_length_n_inv, d_non_unit_normal, d_unit_normal, true);

    Deriv2nd_NonUnitSlaveNormal(xi, mo_ele, dd_non_unit_normal);

    Deriv2nd_UnitSlaveNormal(eletype, adj_ele_unit_n, adj_ele_length_n_inv, d_non_unit_normal,
        d_unit_normal, dd_non_unit_normal, dd_nodal_avg_normal);
  }

  const CORE::LINALG::Matrix<3, 1> avg_unit_normal(cnode.MoData().n(), true);
  const Deriv1stVecMap& d_avg_unit_normal = cnode.AugData().GetDeriv1st_N();

  Deriv2nd_UnitSlaveNormal(eletype, avg_unit_normal, 1.0 / avg_normal_length, d_nodal_avg_normal,
      d_avg_unit_normal, dd_nodal_avg_normal, dd_avg_unit_normal);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::INTEGRATOR::Deriv2nd_NonUnitSlaveNormal(
    const double* xi, MORTAR::MortarElement& sele, Deriv2ndVecMap& dd_non_unit_normal)
{
  const CORE::FE::CellType slavetype = sele.Shape();
  switch (slavetype)
  {
    case CORE::FE::CellType::line2:
    {
      Deriv2nd_NonUnitSlaveNormal<CORE::FE::CellType::line2>(sele, xi, dd_non_unit_normal);
      break;
    }
    case CORE::FE::CellType::nurbs2:
    {
      Deriv2nd_NonUnitSlaveNormal<CORE::FE::CellType::nurbs2>(sele, xi, dd_non_unit_normal);
      break;
    }
    case CORE::FE::CellType::nurbs3:
    {
      Deriv2nd_NonUnitSlaveNormal<CORE::FE::CellType::nurbs3>(sele, xi, dd_non_unit_normal);
      break;
    }
    case CORE::FE::CellType::quad4:
    {
      Deriv2nd_NonUnitSlaveNormal<CORE::FE::CellType::quad4>(sele, xi, dd_non_unit_normal);
      break;
    }
    case CORE::FE::CellType::tri3:
    {
      Deriv2nd_NonUnitSlaveNormal<CORE::FE::CellType::tri3>(sele, xi, dd_non_unit_normal);
      break;
    }
    case CORE::FE::CellType::nurbs4:
    {
      Deriv2nd_NonUnitSlaveNormal<CORE::FE::CellType::nurbs4>(sele, xi, dd_non_unit_normal);
      break;
    }
    case CORE::FE::CellType::nurbs9:
    {
      Deriv2nd_NonUnitSlaveNormal<CORE::FE::CellType::nurbs9>(sele, xi, dd_non_unit_normal);
      break;
    }
    default:
    {
      dserror("Unsupported slave element type! (enum = %d|\"%s\")", slavetype,
          DRT::DistypeToString(slavetype).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <CORE::FE::CellType slavetype>
void CONTACT::INTEGRATOR::Deriv2nd_NonUnitSlaveNormal(
    MORTAR::MortarElement& sele, const double* xi, Deriv2ndVecMap& dd_non_unit_normal)
{
  const unsigned slavenumnode =
      CORE::DRT::UTILS::DisTypeToNumNodePerEle<slavetype>::numNodePerElement;
  const unsigned slavedim = CORE::DRT::UTILS::DisTypeToDim<slavetype>::dim;
  const unsigned probdim = slavedim + 1;

  CORE::LINALG::Matrix<probdim, slavenumnode, int> nodal_dofs;
  CONTACT::INTEGRATOR::GetElementNodalDofs(sele, nodal_dofs);

  // derivatives with respect to xi^1 and xi^2
  CORE::LINALG::Matrix<slavedim, slavenumnode> deriv(true);

  // derivative of the convective base vectors with respect to the displacement
  const CORE::LINALG::Matrix<2, 1> xi_mat(xi, true);
  CONTACT::AUG::shape_function_deriv1<slavetype>(sele, xi_mat, deriv);

  const CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype> slave_policy;
  slave_policy.Deriv2nd_NonUnitSlaveElementNormal(sele, nodal_dofs, deriv, dd_non_unit_normal);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::INTEGRATOR::Deriv2nd_UnitSlaveNormal(const CORE::FE::CellType slavetype,
    const CORE::LINALG::Matrix<3, 1>& unit_normal, const double length_n_inv,
    const Deriv1stVecMap& d_non_unit_normal, const Deriv1stVecMap& d_unit_normal,
    const Deriv2ndVecMap& dd_non_unit_normal, Deriv2ndVecMap& dd_unit_normal)
{
  switch (slavetype)
  {
    case CORE::FE::CellType::line2:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::line2>::dim;
      const unsigned probdim = eledim + 1;

      const CORE::LINALG::Matrix<probdim, 1> unit_normal_red(unit_normal.A(), true);

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::line2> slave_policy;
      slave_policy.Deriv2nd_UnitSlaveElementNormal(unit_normal_red, length_n_inv, d_non_unit_normal,
          d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    case CORE::FE::CellType::nurbs2:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::nurbs2>::dim;
      const unsigned probdim = eledim + 1;

      const CORE::LINALG::Matrix<probdim, 1> unit_normal_red(unit_normal.A(), true);

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::nurbs2> slave_policy;
      slave_policy.Deriv2nd_UnitSlaveElementNormal(unit_normal_red, length_n_inv, d_non_unit_normal,
          d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    case CORE::FE::CellType::nurbs3:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::nurbs3>::dim;
      const unsigned probdim = eledim + 1;

      const CORE::LINALG::Matrix<probdim, 1> unit_normal_red(unit_normal.A(), true);

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::nurbs3> slave_policy;
      slave_policy.Deriv2nd_UnitSlaveElementNormal(unit_normal_red, length_n_inv, d_non_unit_normal,
          d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    case CORE::FE::CellType::quad4:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::quad4>::dim;
      const unsigned probdim = eledim + 1;

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::quad4> slave_policy;
      slave_policy.Deriv2nd_UnitSlaveElementNormal(unit_normal, length_n_inv, d_non_unit_normal,
          d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    case CORE::FE::CellType::tri3:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::tri3>::dim;
      const unsigned probdim = eledim + 1;

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::tri3> slave_policy;
      slave_policy.Deriv2nd_UnitSlaveElementNormal(unit_normal, length_n_inv, d_non_unit_normal,
          d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    case CORE::FE::CellType::nurbs4:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::nurbs4>::dim;
      const unsigned probdim = eledim + 1;

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::nurbs4> slave_policy;
      slave_policy.Deriv2nd_UnitSlaveElementNormal(unit_normal, length_n_inv, d_non_unit_normal,
          d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    case CORE::FE::CellType::nurbs9:
    {
      const unsigned eledim = CORE::DRT::UTILS::DisTypeToDim<CORE::FE::CellType::nurbs9>::dim;
      const unsigned probdim = eledim + 1;

      const CONTACT::AUG::BaseSlaveIntPolicy<probdim, CORE::FE::CellType::nurbs9> slave_policy;
      slave_policy.Deriv2nd_UnitSlaveElementNormal(unit_normal, length_n_inv, d_non_unit_normal,
          d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    default:
    {
      dserror("Unsupported slave element type! (enum = %d|\"%s\")", slavetype,
          DRT::DistypeToString(slavetype).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, unsigned numnode>
void CONTACT::INTEGRATOR::GetElementNodalDofs(
    const MORTAR::MortarElement& ele, CORE::LINALG::Matrix<probdim, numnode, int>& nodal_dofs)
{
  const DRT::Node* const* mynodes = ele.Nodes();

  for (unsigned i = 0; i < numnode; ++i)
  {
    const auto& mynode = dynamic_cast<const MORTAR::MortarNode&>(*mynodes[i]);

    std::copy(mynode.Dofs().data(), mynode.Dofs().data() + probdim, &nodal_dofs(0, i));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::INTEGRATOR::LeviCivitaSymbol(const int i, const int j, const int k)
{
  const int unique_id = i + 4 * j + k * 16;
  switch (unique_id)
  {
    case 9:   // 1 2 0
    case 18:  // 2 0 1
    case 36:  // 0 1 2
    {
      return 1.0;
    }
    case 6:   // 2 1 0
    case 24:  // 0 2 1
    case 33:  // 1 0 2
    {
      return -1.0;
    }
    default:
      return 0.0;
  }
}

/*----------------------------------------------------------------------------*/
template void CONTACT::INTEGRATOR::GetElementNodalDofs(
    const MORTAR::MortarElement& ele, CORE::LINALG::Matrix<2, 2, int>& nodal_dofs);
template void CONTACT::INTEGRATOR::GetElementNodalDofs(
    const MORTAR::MortarElement& ele, CORE::LINALG::Matrix<2, 3, int>& nodal_dofs);
template void CONTACT::INTEGRATOR::GetElementNodalDofs(
    const MORTAR::MortarElement& ele, CORE::LINALG::Matrix<3, 3, int>& nodal_dofs);
template void CONTACT::INTEGRATOR::GetElementNodalDofs(
    const MORTAR::MortarElement& ele, CORE::LINALG::Matrix<3, 4, int>& nodal_dofs);
template void CONTACT::INTEGRATOR::GetElementNodalDofs(
    const MORTAR::MortarElement& ele, CORE::LINALG::Matrix<3, 9, int>& nodal_dofs);
