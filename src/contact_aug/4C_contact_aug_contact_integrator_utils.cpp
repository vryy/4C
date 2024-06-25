/*---------------------------------------------------------------------*/
/*! \file
\brief Utility methods for the contact integration.

\level 2

*/
/*---------------------------------------------------------------------*/

#include "4C_contact_aug_contact_integrator_utils.hpp"

#include "4C_contact_aug_element_utils.hpp"
#include "4C_contact_aug_integrator_policy.hpp"
#include "4C_contact_aug_projector.hpp"
#include "4C_contact_aug_utils.hpp"
#include "4C_contact_integrator.hpp"
#include "4C_contact_node.hpp"
#include "4C_io_pstream.hpp"
#include "4C_mortar_element.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

// #define DEBUG_FIND_FEASIBLE_MASTER_ELEMENT

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::INTEGRATOR::find_feasible_master_elements(Mortar::Element& sele,
    const std::vector<Mortar::Element*>& meles, bool boundary_ele, Integrator& wrapper,
    UniqueProjInfoPair& projInfo)
{
  TEUCHOS_FUNC_TIME_MONITOR(Aug::CONTACT_FUNC_NAME);

  const Core::FE::CellType slavetype = sele.Shape();

  const unsigned msize = meles.size();
  const unsigned numGP = wrapper.nGP();
  const unsigned reserve_size = (msize ? numGP / msize : 0);

  const unsigned probdim = wrapper.Dim();

  Core::Gen::reset(msize, projInfo);

  std::vector<double> found_alpha(msize, 0.0);
  std::vector<Core::LinAlg::Matrix<2, 1>> found_mxi(msize, Core::LinAlg::Matrix<2, 1>(true));
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
      double* mxi = found_mxi[m].data();
      double& projalpha = found_alpha[m];

      const Core::FE::CellType mastertype = meles[m]->Shape();
      // project Gauss point onto master element
#ifdef DEBUG_FIND_FEASIBLE_MASTER_ELEMENT
      CONTACT::Aug::ProjectorBase* proj_ptr = nullptr;
      if (sele.Id() == -1 or sele.Id() == -1)
        proj_ptr = CONTACT::Aug::ProjectorBase::get(probdim, slavetype, mastertype, true);
      else
        proj_ptr = CONTACT::Aug::ProjectorBase::get(probdim, slavetype, mastertype);
#else
      CONTACT::Aug::ProjectorBase* proj_ptr =
          CONTACT::Aug::ProjectorBase::get(probdim, slavetype, mastertype);
#endif

      CONTACT::Aug::ProjectorBase& proj = *proj_ptr;
      const bool conv = proj(sele, sxi, *meles[m], mxi, projalpha);

      // --- DEBUGGING --------------------------------------------------------
#ifdef DEBUG_FIND_FEASIBLE_MASTER_ELEMENT
      Core::Nodes::Node** snodes = sele.Nodes();
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
      proj_tol = 10.0 * proj.get_relative_solution_tolerance();
      is_on_meles[m] = WithinBounds(mxi, mastertype, proj_tol);

      if (is_on_meles[m] and not conv)
        FOUR_C_THROW(
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

        Mortar::Element* mele = meles[unique_id];
        if (projInfo.find(mele) == projInfo.end())
        {
          projInfo[mele] = UniqueProjInfo(std::ceil(reserve_size));
        }

        projInfo[mele].Insert(gp, found_alpha[unique_id], found_mxi[unique_id].data(), gp_scaling);
      }
    }
    else if (not boundary_ele)
      std::cout << "*** warning *** Non-boundary element has non-projectable "
                   "GP \n";
  }

  // --- DEBUGGING --------------------------------------------------------
#ifdef DEBUG_FIND_FEASIBLE_MASTER_ELEMENT
  Core::Nodes::Node** snodes = sele.Nodes();
  const std::vector<int> desired_nids = {247, 248, 251};
  unsigned count = 0;

  for (unsigned i = 0; i < static_cast<unsigned>(sele.num_node()); ++i)
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
      sele.print(std::cout);
      std::cout << "\nSlave-nodes: ";
      for (unsigned i = 0; i < static_cast<unsigned>(sele.num_node()); ++i)
      {
        std::cout << "#" << snodes[i]->Id() << " --> " << snodes[i]->X()[0] << ", "
                  << snodes[i]->X()[1] << ", " << snodes[i]->X()[2] << "\n";
      }
      std::cout << "\n";
      std::cout << "mele " << proj.first->Id() << "[" << proj.first << "]"
                << "\n";
      Core::Nodes::Node** mnodes = proj.first->Nodes();
      std::cout << "\nMaster-nodes: ";
      for (unsigned i = 0; i < static_cast<unsigned>(proj.first->num_node()); ++i)
      {
        std::cout << "#" << mnodes[i]->Id() << " --> " << mnodes[i]->X()[0] << ", "
                  << mnodes[i]->X()[1] << ", " << mnodes[i]->X()[2] << "\n";
      }
      proj.second.print(std::cout);
    }
  }
#endif


  return (projInfo.size() > 0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::INTEGRATOR::WithinBounds(
    const double* mxi, const Core::FE::CellType type, const double tol)
{
  switch (type)
  {
    case Core::FE::CellType::quad4:
    case Core::FE::CellType::quad8:
    case Core::FE::CellType::quad9:
    case Core::FE::CellType::nurbs4:
    case Core::FE::CellType::nurbs9:
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
        return false;

      break;
    }
    case Core::FE::CellType::tri3:
    case Core::FE::CellType::tri6:
    {
      if (mxi[0] < -tol or mxi[1] < -tol or mxi[0] > 1.0 + tol or mxi[1] > 1.0 + tol or
          mxi[0] + mxi[1] > 1.0 + 2 * tol)
        return false;

      break;
    }
    case Core::FE::CellType::line2:
    case Core::FE::CellType::line3:
    case Core::FE::CellType::nurbs2:
    case Core::FE::CellType::nurbs3:
    {
      if ((mxi[0] < -1.0 - tol) or (mxi[0] > 1.0 + tol)) return false;

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported element type %s!", Core::FE::CellTypeToString(type).c_str());
      exit(EXIT_FAILURE);
    }
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::INTEGRATOR::BuildAveragedNormalAtSlaveNode(
    std::vector<ElementNormal>& adj_ele_normals, Mortar::Node& slavenode)
{
  // access the averaged nodal normal
  Core::LinAlg::Matrix<3, 1> avg_nodal_normal(slavenode.MoData().n(), true);

  // reset normal
  std::fill(avg_nodal_normal.data(), avg_nodal_normal.data() + 3, 0.0);

  const unsigned num_adj_eles = slavenode.NumElement();
  Core::Elements::Element** adj_eles = slavenode.Elements();

  adj_ele_normals.clear();
  adj_ele_normals.resize(num_adj_eles, ElementNormal());

  for (unsigned e = 0; e < num_adj_eles; ++e)
  {
    Mortar::Element& adj_ele = dynamic_cast<Mortar::Element&>(*(adj_eles[e]));

    adj_ele_normals[e].ele_ = &adj_ele;

    // get parametric coordinates of the current slave node in the adjacent element
    const int ele_node_lid = adj_ele.GetLocalNodeId(slavenode.Id());
    double* xi = adj_ele_normals[e].xi_;
    adj_ele.local_coordinates_of_node(ele_node_lid, xi);

    // evaluate the convective base vectors
    Core::LinAlg::Matrix<3, 2> tau;
    adj_ele.Metrics(xi, &tau(0, 0), &tau(0, 1));

    // evaluate the unit normal of the adjacent element at the current slave node
    Core::LinAlg::Matrix<3, 1>& unit_normal = adj_ele_normals[e].unit_n_;
    double& length_n_inv = adj_ele_normals[e].length_n_inv_;
    length_n_inv = unit_slave_element_normal(adj_ele, tau, unit_normal);

    // sum up all nodal unit element normals of the adjacent elements
    avg_nodal_normal.update(1.0, unit_normal, 1.0);
  }

  // average the result
  const double smooth_nodal_normal_length = avg_nodal_normal.norm2();

  if (smooth_nodal_normal_length == 0.0)
    FOUR_C_THROW("Sum of slave unit normals at node %d has a length of zero!", slavenode.Id());

  avg_nodal_normal.scale(1.0 / smooth_nodal_normal_length);

  return smooth_nodal_normal_length;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::INTEGRATOR::unit_slave_element_normal(const Mortar::Element& sele,
    const Core::LinAlg::Matrix<3, 2>& tau, Core::LinAlg::Matrix<3, 1>& unit_normal)
{
  const Core::FE::CellType slavetype = sele.Shape();
  switch (slavetype)
  {
    case Core::FE::CellType::line2:
    {
      const CONTACT::Aug::BaseSlaveIntPolicy<2, Core::FE::CellType::line2> slave_policy;

      Core::LinAlg::Matrix<2, 1> unit_n(unit_normal.data(), true);
      return slave_policy.unit_slave_element_normal(sele, tau, unit_n);
    }
    case Core::FE::CellType::nurbs2:
    {
      const CONTACT::Aug::BaseSlaveIntPolicy<2, Core::FE::CellType::nurbs2> slave_policy;

      Core::LinAlg::Matrix<2, 1> unit_n(unit_normal.data(), true);
      return slave_policy.unit_slave_element_normal(sele, tau, unit_n);
    }
    case Core::FE::CellType::nurbs3:
    {
      const CONTACT::Aug::BaseSlaveIntPolicy<2, Core::FE::CellType::nurbs3> slave_policy;

      Core::LinAlg::Matrix<2, 1> unit_n(unit_normal.data(), true);
      return slave_policy.unit_slave_element_normal(sele, tau, unit_n);
    }
    case Core::FE::CellType::quad4:
    {
      const CONTACT::Aug::BaseSlaveIntPolicy<3, Core::FE::CellType::quad4> slave_policy;

      return slave_policy.unit_slave_element_normal(sele, tau, unit_normal);
    }
    case Core::FE::CellType::tri3:
    {
      const CONTACT::Aug::BaseSlaveIntPolicy<3, Core::FE::CellType::tri3> slave_policy;

      return slave_policy.unit_slave_element_normal(sele, tau, unit_normal);
    }
    case Core::FE::CellType::nurbs4:
    {
      const CONTACT::Aug::BaseSlaveIntPolicy<3, Core::FE::CellType::nurbs4> slave_policy;

      return slave_policy.unit_slave_element_normal(sele, tau, unit_normal);
    }
    case Core::FE::CellType::nurbs9:
    {
      const CONTACT::Aug::BaseSlaveIntPolicy<3, Core::FE::CellType::nurbs9> slave_policy;

      return slave_policy.unit_slave_element_normal(sele, tau, unit_normal);
    }
    default:
    {
      FOUR_C_THROW("Unsupported slave element type! (enum = %d|\"%s\")", slavetype,
          Core::FE::CellTypeToString(slavetype).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::INTEGRATOR::Deriv1st_AveragedSlaveNormal(CONTACT::Node& cnode,
    const std::vector<ElementNormal>& adj_ele_normals, const double avg_normal_length,
    Deriv1stVecMap& d_nodal_avg_normal)
{
  Deriv1stVecMap& d_avg_unit_normal = cnode.AugData().GetDeriv1st_N();
  Core::Gen::reset(3, cnode.GetLinsize(), d_avg_unit_normal);
  Core::Gen::reset(3, cnode.GetLinsize(), d_nodal_avg_normal);

  Deriv1stVecMap d_non_unit_normal;

  Core::FE::CellType eletype = Core::FE::CellType::dis_none;
  for (const ElementNormal& adj_ele_n : adj_ele_normals)
  {
    Mortar::Element& mo_ele = *adj_ele_n.ele_;

    eletype = mo_ele.Shape();

    const double* xi = adj_ele_n.xi_;

    Deriv1st_NonUnitSlaveNormal(xi, mo_ele, d_non_unit_normal);

    const Core::LinAlg::Matrix<3, 1>& adj_ele_unit_n = adj_ele_n.unit_n_;
    const double adj_ele_length_n_inv = adj_ele_n.length_n_inv_;

    Deriv1st_UnitSlaveNormal(eletype, adj_ele_unit_n, adj_ele_length_n_inv, d_non_unit_normal,
        d_nodal_avg_normal, false);
  }

  const Core::LinAlg::Matrix<3, 1> avg_unit_normal(cnode.MoData().n(), true);
  Deriv1st_UnitSlaveNormal(eletype, avg_unit_normal, 1.0 / avg_normal_length, d_nodal_avg_normal,
      d_avg_unit_normal, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::INTEGRATOR::Deriv1st_UnitSlaveNormal(const Core::FE::CellType slavetype,
    const Core::LinAlg::Matrix<3, 1>& unit_normal, const double length_n_inv,
    const Deriv1stVecMap& d_non_unit_normal, Deriv1stVecMap& d_unit_normal, const bool reset)
{
  switch (slavetype)
  {
    case Core::FE::CellType::line2:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::line2>;
      const unsigned probdim = eledim + 1;

      const Core::LinAlg::Matrix<probdim, 1> unit_normal_red(unit_normal.data(), true);

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::line2> slave_policy;
      slave_policy.deriv1st_unit_slave_element_normal(
          unit_normal_red, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    case Core::FE::CellType::nurbs2:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::nurbs2>;
      const unsigned probdim = eledim + 1;

      const Core::LinAlg::Matrix<probdim, 1> unit_normal_red(unit_normal.data(), true);

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::nurbs2> slave_policy;
      slave_policy.deriv1st_unit_slave_element_normal(
          unit_normal_red, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    case Core::FE::CellType::nurbs3:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::nurbs3>;
      const unsigned probdim = eledim + 1;

      const Core::LinAlg::Matrix<probdim, 1> unit_normal_red(unit_normal.data(), true);

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::nurbs3> slave_policy;
      slave_policy.deriv1st_unit_slave_element_normal(
          unit_normal_red, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    case Core::FE::CellType::quad4:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::quad4>;
      const unsigned probdim = eledim + 1;

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::quad4> slave_policy;
      slave_policy.deriv1st_unit_slave_element_normal(
          unit_normal, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    case Core::FE::CellType::tri3:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::tri3>;
      const unsigned probdim = eledim + 1;

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::tri3> slave_policy;
      slave_policy.deriv1st_unit_slave_element_normal(
          unit_normal, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    case Core::FE::CellType::nurbs4:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::nurbs4>;
      const unsigned probdim = eledim + 1;

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::nurbs4> slave_policy;
      slave_policy.deriv1st_unit_slave_element_normal(
          unit_normal, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    case Core::FE::CellType::nurbs9:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::nurbs9>;
      const unsigned probdim = eledim + 1;

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::nurbs9> slave_policy;
      slave_policy.deriv1st_unit_slave_element_normal(
          unit_normal, length_n_inv, d_non_unit_normal, d_unit_normal, reset);

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported slave element type! (enum = %d|\"%s\")", slavetype,
          Core::FE::CellTypeToString(slavetype).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::INTEGRATOR::Deriv1st_NonUnitSlaveNormal(
    const double* xi, Mortar::Element& sele, Deriv1stVecMap& d_non_unit_normal)
{
  const Core::FE::CellType slavetype = sele.Shape();
  switch (slavetype)
  {
    case Core::FE::CellType::line2:
    {
      Deriv1st_NonUnitSlaveNormal<Core::FE::CellType::line2>(sele, xi, d_non_unit_normal);
      break;
    }
    case Core::FE::CellType::nurbs2:
    {
      Deriv1st_NonUnitSlaveNormal<Core::FE::CellType::nurbs2>(sele, xi, d_non_unit_normal);
      break;
    }
    case Core::FE::CellType::nurbs3:
    {
      Deriv1st_NonUnitSlaveNormal<Core::FE::CellType::nurbs3>(sele, xi, d_non_unit_normal);
      break;
    }
    case Core::FE::CellType::quad4:
    {
      Deriv1st_NonUnitSlaveNormal<Core::FE::CellType::quad4>(sele, xi, d_non_unit_normal);
      break;
    }
    case Core::FE::CellType::tri3:
    {
      Deriv1st_NonUnitSlaveNormal<Core::FE::CellType::tri3>(sele, xi, d_non_unit_normal);
      break;
    }
    case Core::FE::CellType::nurbs4:
    {
      Deriv1st_NonUnitSlaveNormal<Core::FE::CellType::nurbs4>(sele, xi, d_non_unit_normal);
      break;
    }
    case Core::FE::CellType::nurbs9:
    {
      Deriv1st_NonUnitSlaveNormal<Core::FE::CellType::nurbs9>(sele, xi, d_non_unit_normal);
      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported slave element type! (enum = %d|\"%s\")", slavetype,
          Core::FE::CellTypeToString(slavetype).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <Core::FE::CellType slavetype>
void CONTACT::INTEGRATOR::Deriv1st_NonUnitSlaveNormal(
    Mortar::Element& sele, const double* xi, Deriv1stVecMap& d_non_unit_normal)
{
  const unsigned slavenumnode = Core::FE::num_nodes<slavetype>;
  const unsigned slavedim = Core::FE::dim<slavetype>;
  const unsigned probdim = slavedim + 1;

  Core::LinAlg::Matrix<probdim, slavenumnode, int> nodal_dofs;
  CONTACT::INTEGRATOR::GetElementNodalDofs(sele, nodal_dofs);

  // evaluate the convective base vectors
  Core::LinAlg::Matrix<3, 2> tau;
  sele.Metrics(xi, &tau(0, 0), &tau(0, 1));

  // derivatives with respect to xi^1 and xi^2
  Core::LinAlg::Matrix<slavedim, slavenumnode> deriv(true);

  // derivative of the convective base vectors with respect to the displacement
  const Core::LinAlg::Matrix<2, 1> xi_mat(xi, true);
  CONTACT::Aug::shape_function_deriv1<slavetype>(sele, xi_mat, deriv);

  const CONTACT::Aug::BaseSlaveIntPolicy<probdim, slavetype> slave_policy;
  slave_policy.deriv1st_non_unit_slave_element_normal(
      sele, nodal_dofs, deriv, tau, d_non_unit_normal);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::INTEGRATOR::Deriv2nd_AveragedSlaveNormal(CONTACT::Node& cnode,
    const std::vector<ElementNormal>& adj_ele_normals, const double avg_normal_length,
    const Deriv1stVecMap& d_nodal_avg_normal)
{
  Deriv2ndVecMap& dd_avg_unit_normal = cnode.AugData().GetDeriv2nd_N();
  Core::Gen::reset(3, cnode.GetLinsize(), dd_avg_unit_normal);

  Deriv2ndVecMap dd_nodal_avg_normal;
  Core::Gen::reset(3, cnode.GetLinsize(), dd_nodal_avg_normal);

  Deriv1stVecMap d_non_unit_normal;
  Deriv1stVecMap d_unit_normal;

  Deriv2ndVecMap dd_non_unit_normal(3, Deriv2ndMap());

  Core::FE::CellType eletype = Core::FE::CellType::dis_none;
  for (const ElementNormal& adj_ele_n : adj_ele_normals)
  {
    Mortar::Element& mo_ele = *adj_ele_n.ele_;

    eletype = mo_ele.Shape();

    const double* xi = adj_ele_n.xi_;

    Deriv1st_NonUnitSlaveNormal(xi, mo_ele, d_non_unit_normal);

    const Core::LinAlg::Matrix<3, 1>& adj_ele_unit_n = adj_ele_n.unit_n_;
    const double adj_ele_length_n_inv = adj_ele_n.length_n_inv_;

    Deriv1st_UnitSlaveNormal(
        eletype, adj_ele_unit_n, adj_ele_length_n_inv, d_non_unit_normal, d_unit_normal, true);

    Deriv2nd_NonUnitSlaveNormal(xi, mo_ele, dd_non_unit_normal);

    Deriv2nd_UnitSlaveNormal(eletype, adj_ele_unit_n, adj_ele_length_n_inv, d_non_unit_normal,
        d_unit_normal, dd_non_unit_normal, dd_nodal_avg_normal);
  }

  const Core::LinAlg::Matrix<3, 1> avg_unit_normal(cnode.MoData().n(), true);
  const Deriv1stVecMap& d_avg_unit_normal = cnode.AugData().GetDeriv1st_N();

  Deriv2nd_UnitSlaveNormal(eletype, avg_unit_normal, 1.0 / avg_normal_length, d_nodal_avg_normal,
      d_avg_unit_normal, dd_nodal_avg_normal, dd_avg_unit_normal);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::INTEGRATOR::Deriv2nd_NonUnitSlaveNormal(
    const double* xi, Mortar::Element& sele, Deriv2ndVecMap& dd_non_unit_normal)
{
  const Core::FE::CellType slavetype = sele.Shape();
  switch (slavetype)
  {
    case Core::FE::CellType::line2:
    {
      Deriv2nd_NonUnitSlaveNormal<Core::FE::CellType::line2>(sele, xi, dd_non_unit_normal);
      break;
    }
    case Core::FE::CellType::nurbs2:
    {
      Deriv2nd_NonUnitSlaveNormal<Core::FE::CellType::nurbs2>(sele, xi, dd_non_unit_normal);
      break;
    }
    case Core::FE::CellType::nurbs3:
    {
      Deriv2nd_NonUnitSlaveNormal<Core::FE::CellType::nurbs3>(sele, xi, dd_non_unit_normal);
      break;
    }
    case Core::FE::CellType::quad4:
    {
      Deriv2nd_NonUnitSlaveNormal<Core::FE::CellType::quad4>(sele, xi, dd_non_unit_normal);
      break;
    }
    case Core::FE::CellType::tri3:
    {
      Deriv2nd_NonUnitSlaveNormal<Core::FE::CellType::tri3>(sele, xi, dd_non_unit_normal);
      break;
    }
    case Core::FE::CellType::nurbs4:
    {
      Deriv2nd_NonUnitSlaveNormal<Core::FE::CellType::nurbs4>(sele, xi, dd_non_unit_normal);
      break;
    }
    case Core::FE::CellType::nurbs9:
    {
      Deriv2nd_NonUnitSlaveNormal<Core::FE::CellType::nurbs9>(sele, xi, dd_non_unit_normal);
      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported slave element type! (enum = %d|\"%s\")", slavetype,
          Core::FE::CellTypeToString(slavetype).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <Core::FE::CellType slavetype>
void CONTACT::INTEGRATOR::Deriv2nd_NonUnitSlaveNormal(
    Mortar::Element& sele, const double* xi, Deriv2ndVecMap& dd_non_unit_normal)
{
  const unsigned slavenumnode = Core::FE::num_nodes<slavetype>;
  const unsigned slavedim = Core::FE::dim<slavetype>;
  const unsigned probdim = slavedim + 1;

  Core::LinAlg::Matrix<probdim, slavenumnode, int> nodal_dofs;
  CONTACT::INTEGRATOR::GetElementNodalDofs(sele, nodal_dofs);

  // derivatives with respect to xi^1 and xi^2
  Core::LinAlg::Matrix<slavedim, slavenumnode> deriv(true);

  // derivative of the convective base vectors with respect to the displacement
  const Core::LinAlg::Matrix<2, 1> xi_mat(xi, true);
  CONTACT::Aug::shape_function_deriv1<slavetype>(sele, xi_mat, deriv);

  const CONTACT::Aug::BaseSlaveIntPolicy<probdim, slavetype> slave_policy;
  slave_policy.deriv2nd_non_unit_slave_element_normal(sele, nodal_dofs, deriv, dd_non_unit_normal);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::INTEGRATOR::Deriv2nd_UnitSlaveNormal(const Core::FE::CellType slavetype,
    const Core::LinAlg::Matrix<3, 1>& unit_normal, const double length_n_inv,
    const Deriv1stVecMap& d_non_unit_normal, const Deriv1stVecMap& d_unit_normal,
    const Deriv2ndVecMap& dd_non_unit_normal, Deriv2ndVecMap& dd_unit_normal)
{
  switch (slavetype)
  {
    case Core::FE::CellType::line2:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::line2>;
      const unsigned probdim = eledim + 1;

      const Core::LinAlg::Matrix<probdim, 1> unit_normal_red(unit_normal.data(), true);

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::line2> slave_policy;
      slave_policy.deriv2nd_unit_slave_element_normal(unit_normal_red, length_n_inv,
          d_non_unit_normal, d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    case Core::FE::CellType::nurbs2:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::nurbs2>;
      const unsigned probdim = eledim + 1;

      const Core::LinAlg::Matrix<probdim, 1> unit_normal_red(unit_normal.data(), true);

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::nurbs2> slave_policy;
      slave_policy.deriv2nd_unit_slave_element_normal(unit_normal_red, length_n_inv,
          d_non_unit_normal, d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    case Core::FE::CellType::nurbs3:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::nurbs3>;
      const unsigned probdim = eledim + 1;

      const Core::LinAlg::Matrix<probdim, 1> unit_normal_red(unit_normal.data(), true);

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::nurbs3> slave_policy;
      slave_policy.deriv2nd_unit_slave_element_normal(unit_normal_red, length_n_inv,
          d_non_unit_normal, d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    case Core::FE::CellType::quad4:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::quad4>;
      const unsigned probdim = eledim + 1;

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::quad4> slave_policy;
      slave_policy.deriv2nd_unit_slave_element_normal(unit_normal, length_n_inv, d_non_unit_normal,
          d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    case Core::FE::CellType::tri3:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::tri3>;
      const unsigned probdim = eledim + 1;

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::tri3> slave_policy;
      slave_policy.deriv2nd_unit_slave_element_normal(unit_normal, length_n_inv, d_non_unit_normal,
          d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    case Core::FE::CellType::nurbs4:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::nurbs4>;
      const unsigned probdim = eledim + 1;

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::nurbs4> slave_policy;
      slave_policy.deriv2nd_unit_slave_element_normal(unit_normal, length_n_inv, d_non_unit_normal,
          d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    case Core::FE::CellType::nurbs9:
    {
      const unsigned eledim = Core::FE::dim<Core::FE::CellType::nurbs9>;
      const unsigned probdim = eledim + 1;

      const CONTACT::Aug::BaseSlaveIntPolicy<probdim, Core::FE::CellType::nurbs9> slave_policy;
      slave_policy.deriv2nd_unit_slave_element_normal(unit_normal, length_n_inv, d_non_unit_normal,
          d_unit_normal, dd_non_unit_normal, dd_unit_normal);

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported slave element type! (enum = %d|\"%s\")", slavetype,
          Core::FE::CellTypeToString(slavetype).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, unsigned numnode>
void CONTACT::INTEGRATOR::GetElementNodalDofs(
    const Mortar::Element& ele, Core::LinAlg::Matrix<probdim, numnode, int>& nodal_dofs)
{
  const Core::Nodes::Node* const* mynodes = ele.Nodes();

  for (unsigned i = 0; i < numnode; ++i)
  {
    const auto& mynode = dynamic_cast<const Mortar::Node&>(*mynodes[i]);

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
    const Mortar::Element& ele, Core::LinAlg::Matrix<2, 2, int>& nodal_dofs);
template void CONTACT::INTEGRATOR::GetElementNodalDofs(
    const Mortar::Element& ele, Core::LinAlg::Matrix<2, 3, int>& nodal_dofs);
template void CONTACT::INTEGRATOR::GetElementNodalDofs(
    const Mortar::Element& ele, Core::LinAlg::Matrix<3, 3, int>& nodal_dofs);
template void CONTACT::INTEGRATOR::GetElementNodalDofs(
    const Mortar::Element& ele, Core::LinAlg::Matrix<3, 4, int>& nodal_dofs);
template void CONTACT::INTEGRATOR::GetElementNodalDofs(
    const Mortar::Element& ele, Core::LinAlg::Matrix<3, 9, int>& nodal_dofs);

FOUR_C_NAMESPACE_CLOSE
