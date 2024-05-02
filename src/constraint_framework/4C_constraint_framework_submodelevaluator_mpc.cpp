/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all multipoint constraint terms resulting
 from periodic displacement boundary conditions for representative volume elements (RVEs)
 and linear coupled degrees of freedom.

\level 3
*/
/*-----------------------------------------------------------*/

#include "4C_config.hpp"

#include "4C_constraint_framework_submodelevaluator_mpc.hpp"

#include "4C_beam3_base.hpp"
#include "4C_constraint_framework_equation_mpc.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_mpc_rve.hpp"
#include "4C_io.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_structure_new_timint_implicit.hpp"

#include <algorithm>
#include <iostream>
#include <vector>
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_OPEN
CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::RveMultiPointConstraintManager(
    Teuchos::RCP<const DRT::Discretization> disc_ptr, CORE::LINALG::SparseMatrix* st_ptr)
{
  discret_ptr_ = disc_ptr;
  stiff_ptr_ = st_ptr;

  CheckInput();

  // Maps of relevant boundary nodesets of the rve
  std::map<std::string, const std::vector<int>*> rveBoundaryNodeIdMap;

  // Map of the corner node ids
  std::map<std::string, int> rveCornerNodeIdMap;

  //  Map the Node IDs to the respective rve boundary --> rveBoundaryNodeIdMap
  BuildPeriodicRveBoundaryNodeMap(rveBoundaryNodeIdMap);

  // Map the Node ids to the respective corner of the rve --> rveCornerNodeIdMap
  switch (rve_ref_type_)
  {
    case INPAR::RVE_MPC::RveReferenceDeformationDefinition::automatic:
    {
      BuildPeriodicRveCornerNodeMap(rveBoundaryNodeIdMap, rveCornerNodeIdMap);
    }
    break;

    case INPAR::RVE_MPC::RveReferenceDeformationDefinition::manual:
    {
      if (rve_dim_ != INPAR::RVE_MPC::RveDimension::rve2d)
        FOUR_C_THROW("Manual Edge node definition is not implemented for 3D RVEs");

      // Read the reference points
      for (const auto& entry : point_periodic_rve_ref_conditions_)
      {
        const auto& str_id = entry->Get<std::string>("referenceNode");
        const auto* nodeInSet = entry->GetNodes();

        if (nodeInSet->size() > 1)
        {
          FOUR_C_THROW("There can only be a single node defined as a reference node");
        }
        rve_ref_node_map_[str_id.c_str()] = discret_ptr_->gNode(nodeInSet->data()[0]);

        IO::cout(IO::verbose) << "Map " << str_id.c_str() << " to node id " << (*nodeInSet)[0]
                              << IO::endl;
      }

      IO::cout(IO::verbose)
          << "Reference Points determined:"
          << "+--------------------------------------------------------------------+" << IO::endl;
      for (const auto& elem : rve_ref_node_map_)
      {
        IO::cout(IO::verbose) << elem.first << ": " << elem.second->Id() << ", ";
      }
      IO::cout(IO::verbose) << IO::endl;

      // calculate the Reference vectors between Ref. poitns
      r_xmxp_[0] = rve_ref_node_map_["N2"]->X()[0] - rve_ref_node_map_["N1L"]->X()[0];
      r_xmxp_[1] = rve_ref_node_map_["N2"]->X()[1] - rve_ref_node_map_["N1L"]->X()[1];
      IO::cout(IO::verbose) << "RVE reference vector (X- ---> X+ ) : [" << r_xmxp_[0] << ", "
                            << r_xmxp_[1] << "]" << IO::endl;

      r_ymyp_[0] = rve_ref_node_map_["N4"]->X()[0] - rve_ref_node_map_["N1B"]->X()[0];
      r_ymyp_[1] = rve_ref_node_map_["N4"]->X()[1] - rve_ref_node_map_["N1B"]->X()[1];
      IO::cout(IO::verbose) << "RVE reference vector (Y- ---> Y+ ) : [" << r_ymyp_[0] << ", "
                            << r_ymyp_[1] << "]" << IO::endl;
    }
    break;
  }

  // Create a vector with all MPCs describing the periodic BCs
  if (surface_periodic_rve_conditions_.size() != 0 || line_periodic_rve_conditions_.size() != 0)
    BuildPeriodicMPCs(rveBoundaryNodeIdMap, rveCornerNodeIdMap);


  // Add Linear Coupled Equation MPCs
  if (point_linear_coupled_equation_conditions_.size() != 0)
  {
    int nLinCe = BuildLinearMPCs();
    IO::cout(IO::verbose) << "Total number of linear coupled equations: " << nLinCe << IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::CheckInput()
{
  if (discret_ptr_->Comm().NumProc() > 1)
    FOUR_C_THROW("periodic boundary conditions for RVEs are not implemented in parallel.");


  mpc_parameter_list_ = GLOBAL::Problem::Instance()->RveMultiPointConstraintParams();

  rve_ref_type_ = CORE::UTILS::IntegralValue<INPAR::RVE_MPC::RveReferenceDeformationDefinition>(
      mpc_parameter_list_, "RVE_REFERENCE_POINTS");

  strategy_ = CORE::UTILS::IntegralValue<INPAR::RVE_MPC::EnforcementStrategy>(
      mpc_parameter_list_, "ENFORCEMENT");


  // Check the enforcement strategy
  switch (strategy_)
  {
    case INPAR::RVE_MPC::EnforcementStrategy::lagrangeMultiplier:
      FOUR_C_THROW("Constraint Enforcement via Lagrange Multiplier Methode is not impl.");

    case INPAR::RVE_MPC::EnforcementStrategy::penalty:
    {
      IO::cout(IO::minimal) << "Constraint enforcement strategy: Penalty method" << IO::endl;
      const Teuchos::ParameterList& MpcList =
          GLOBAL::Problem::Instance()->RveMultiPointConstraintParams();
      GetPenaltyParameterPtr() = MpcList.get<double>("PENALTY_PARAM");
      IO::cout(IO::verbose) << "Penalty weight used: " << GetPenaltyParameterPtr() << IO::endl;
    }
  }

  // Conditions definition
  discret_ptr_->GetCondition("LinePeriodicRve", line_periodic_rve_conditions_);
  discret_ptr_->GetCondition("SurfacePeriodicRve", surface_periodic_rve_conditions_);
  discret_ptr_->GetCondition("PointPeriodicRveReferenceNode", point_periodic_rve_ref_conditions_);
  discret_ptr_->GetCondition(
      "PointLinearCoupledEquation", point_linear_coupled_equation_conditions_);

  // Input Checks: Dimensions
  if (line_periodic_rve_conditions_.size() == 0 && surface_periodic_rve_conditions_.size() != 0)
  {
    IO::cout(IO::verbose) << "Rve dimension: 3d" << IO::endl;
    rve_dim_ = INPAR::RVE_MPC::RveDimension::rve3d;
  }
  else if (line_periodic_rve_conditions_.size() != 0 &&
           surface_periodic_rve_conditions_.size() == 0)
  {
    IO::cout(IO::verbose) << "Rve dimensions: 2d" << IO::endl;
    rve_dim_ = INPAR::RVE_MPC::RveDimension::rve2d;
  }
  else
  {
    FOUR_C_THROW("Periodic rve edge condition cannot be combinded with peridodic rve surf cond. ");
  }

  // Input Checks
  IO::cout(IO::verbose) << "There are " << line_periodic_rve_conditions_.size()
                        << " periodic rve edge conditions defined (2D)" << IO::endl;
  if (line_periodic_rve_conditions_.size() != 0)
  {
    if (line_periodic_rve_conditions_.size() != 4 && line_periodic_rve_conditions_.size() != 2)
    {
      IO::cout(IO::verbose) << "Number of Conditions: " << line_periodic_rve_conditions_.size()
                            << IO::endl;
      FOUR_C_THROW("For a 2D RVE either all or two opposing edges must be used for PBCs");
    }
  }

  IO::cout(IO::verbose) << "There are " << surface_periodic_rve_conditions_.size()
                        << " periodic rve surface conditions defined (3D)" << IO::endl;

  IO::cout(IO::verbose) << "There are " << point_linear_coupled_equation_conditions_.size()
                        << " linear coupled equations" << IO::endl;

  if (point_periodic_rve_ref_conditions_.size() == 0 &&
      rve_ref_type_ == INPAR::RVE_MPC::RveReferenceDeformationDefinition::manual)
  {
    FOUR_C_THROW(
        "A DESIGN POINT PERIODIC RVE 2D BOUNDARY REFERENCE CONDITIONS is req. for manual ref. "
        "point "
        "definition");
  }

  if (point_periodic_rve_ref_conditions_.size() != 0 &&
      rve_ref_type_ == INPAR::RVE_MPC::RveReferenceDeformationDefinition::automatic)
    FOUR_C_THROW("Set the RVE_REFERENCE_POINTS to manual");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::FindOppositeEdgeNode(
    const int nodeID, INPAR::RVE_MPC::RveEdgeIdentifiers edge,
    std::map<std::string, const std::vector<int>*>& rveBoundaryNodeIdMap)
{
  std::string newPos;
  std::array<double, 2> R_ipim;

  DRT::Node* nodeA = discret_ptr_->gNode(nodeID);

  switch (edge)
  {
    case INPAR::RVE_MPC::RveEdgeIdentifiers::Gamma_xm:
    {
      R_ipim = r_xmxp_;
      newPos = "x+";
      IO::cout(IO::debug) << " Find partner node of Node " << nodeID << " on Edge:  " << newPos
                          << IO::endl;
      break;
    }

    case INPAR::RVE_MPC::RveEdgeIdentifiers::Gamma_ym:
    {
      R_ipim = r_ymyp_;
      newPos = "y+";
      IO::cout(IO::debug) << " Find partner node of Node " << nodeID << " on Edge:  " << newPos
                          << IO::endl;
      break;
    }
    default:
    {
      FOUR_C_THROW("Specifiy the negative edge, 3D not implemented");
    }
  }

  // Calculate the Position of the opposing edge node
  std::vector<double> newPosition = {0.0, 0.0};
  for (int i = 0; i < 2; i++) newPosition[i] = nodeA->X()[i] + R_ipim[i];

  IO::cout(IO::debug) << "Calculated position of matching node: " << newPosition[0] << " / "
                      << newPosition[1] << IO::endl;

  // Loop all nodes of the relevant opposite edge line
  // ToDo: Switch to ArborX
  for (auto pairId : *rveBoundaryNodeIdMap[newPos])
  {
    DRT::Node* nodeB = discret_ptr_->gNode(pairId);

    if (std::abs(nodeB->X()[0] - newPosition[0]) < node_search_toler_)
    {
      if (std::abs(nodeB->X()[1] - newPosition[1]) < node_search_toler_)
      {
        IO::cout(IO::debug) << "Found Node Pair (IDs): A: " << nodeA->Id() << " B: " << pairId
                            << IO::endl;
        return pairId;
      }
    }
  }
  FOUR_C_THROW("No matching node found! - Is the mesh perodic?");
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::BuildPeriodicMPCs(
    std::map<std::string, const std::vector<int>*>& rveBoundaryNodeIdMap,
    std::map<std::string, int>& rveCornerNodeIdMap)
{
  std::vector<std::vector<DRT::Node*>> PBCs;
  std::vector<DRT::Node*> PBC;
  {
    IO::cout(IO::verbose)
        << "\nCreating Node Pairs for Periodic Boundary Conditions" << IO::endl
        << "+--------------------------------------------------------------------+" << IO::endl;
    IO::cout(IO::verbose) << "RVE Type: ";
  }
  switch (rve_ref_type_)
  {
    case INPAR::RVE_MPC::RveReferenceDeformationDefinition::automatic:
    {
      switch (rve_dim_)
      {
        case INPAR::RVE_MPC::RveDimension::rve3d:
        case INPAR::RVE_MPC::RveDimension::rve2d:
        {
          int numDim = 3;
          std::map<std::string, std::string> refEndNodeMap = {
              {"x", "N2"}, {"y", "N4"}, {"z", "N5"}};

          if (rve_dim_ == INPAR::RVE_MPC::rve3d)
          {
            IO::cout(IO::verbose) << "General 3D RVE" << IO::endl;
          }
          else if (rve_dim_ == INPAR::RVE_MPC::rve2d)
          {
            IO::cout(IO::verbose) << "General 2D RVE" << IO::endl;
            refEndNodeMap.erase("z");
            numDim = 2;
          }

          std::map<std::string, std::vector<double>> rveRefVecMap;
          std::vector<double> rveRefVector;

          for (const auto& surf : refEndNodeMap)
          {
            {
              rveRefVector.clear();
              for (int i = 0; i < numDim; ++i)
              {
                rveRefVector.push_back(
                    discret_ptr_->gNode(rveCornerNodeIdMap[surf.second])->X()[i] -
                    discret_ptr_->gNode(rveCornerNodeIdMap["N1"])->X()[i]);
              }
              rveRefVecMap[surf.first] = rveRefVector;

              IO::cout(IO::verbose) << "RVE reference vector " << surf.first << " dimension: "
                                    << "[ " << rveRefVector[0] << "; " << rveRefVector[1];
              if (rve_dim_ == INPAR::RVE_MPC::rve3d)
                IO::cout(IO::verbose) << "; " << rveRefVector[2];
              IO::cout(IO::verbose) << " ]" << IO::endl;
            }
          }
          // Create PBC Node Pairs:
          // #ToDo: Use ArborX for the search
          std::string surfId;
          std::vector<double> matchPosition = {0.0, 0.0, 0.0};

          for (const auto& surf : refEndNodeMap)
          {
            for (auto nodeP : *rveBoundaryNodeIdMap[surf.first + "+"])
            {
              IO::cout(IO::debug) << "Current SURF: " << surf.first << "Node+ ID: " << nodeP
                                  << IO::endl;

              for (int i = 0; i < numDim; ++i)
              {
                matchPosition[i] = discret_ptr_->gNode(nodeP)->X()[i] - rveRefVecMap[surf.first][i];
              }
              IO::cout(IO::debug) << "Node+ location: " << discret_ptr_->gNode(nodeP)->X()[0]
                                  << ", " << discret_ptr_->gNode(nodeP)->X()[1];
              if (rve_dim_ == INPAR::RVE_MPC::rve3d)
                IO::cout(IO::debug) << ", " << discret_ptr_->gNode(nodeP)->X()[2];
              IO::cout(IO::debug) << IO::endl;

              IO::cout(IO::debug) << "Position of matching Node: " << matchPosition[0] << ", "
                                  << matchPosition[1];
              if (rve_dim_ == INPAR::RVE_MPC::rve3d)
                IO::cout(IO::verbose) << ", " << matchPosition[2];
              IO::cout(IO::verbose) << IO::endl;
              for (auto nodeM : *rveBoundaryNodeIdMap[surf.first + "-"])
              {
                if (nodeP == rveCornerNodeIdMap[surf.second]) break;

                bool isMatch = true;
                for (int i = 0; i < numDim; ++i)
                {
                  if (std::abs(discret_ptr_->gNode(nodeM)->X()[i] - matchPosition[i]) >
                      node_search_toler_)
                  {
                    isMatch = false;
                    break;
                  }
                }
                if (isMatch)
                {
                  PBC.push_back(discret_ptr_->gNode(nodeP));
                  PBC.push_back(discret_ptr_->gNode(nodeM));
                  PBC.push_back(discret_ptr_->gNode(rveCornerNodeIdMap[surf.second]));
                  PBC.push_back(discret_ptr_->gNode(rveCornerNodeIdMap["N1"]));
                  PBCs.push_back(PBC);
                  PBC.clear();

                  {
                    IO::cout(IO::debug)
                        << "nodes coupled on " << surf.first << " +/- boundary: " << nodeP << " - "
                        << nodeM << " = " << rveCornerNodeIdMap[surf.second] << " - "
                        << rveCornerNodeIdMap["N1"] << IO::endl;
                  }
                }
              }
            }
          }
        }
        break;
      }
      break;
    }
    case INPAR::RVE_MPC::RveReferenceDeformationDefinition::manual:
    {
      IO::cout(IO::verbose) << "General 2D RVE" << IO::endl;
      for (const auto& elem : rve_ref_node_map_)
      {
        IO::cout(IO::verbose) << elem.first << " first " << elem.second->Id() << " sec "
                              << "\n";
      }

      /* Loop over X- Edge */
      for (auto nodeXm : *rveBoundaryNodeIdMap["x-"])
      {
        if (nodeXm != rve_ref_node_map_["N1L"]->Id())  // exlude N1 - N2 = N1 - N2
        {
          PBC.push_back(discret_ptr_->gNode(nodeXm));
          PBC.push_back(discret_ptr_->gNode(FindOppositeEdgeNode(
              nodeXm, INPAR::RVE_MPC::RveEdgeIdentifiers::Gamma_xm, rveBoundaryNodeIdMap)));
          PBC.push_back(rve_ref_node_map_["N1L"]);
          PBC.push_back(rve_ref_node_map_["N2"]);
          PBCs.push_back(PBC);

          IO::cout(IO::debug) << "PBC MPC Set created X- ---> X+ Edge:  " << IO::endl;
          for (auto* nnn : PBC)
          {
            IO::cout(IO::debug) << " ___ " << nnn->Id();
          }
          IO::cout(IO::debug) << IO::endl;
          PBC.clear();
        }
      }

      /* Loop over Y- Edge*/
      for (auto nodeYm : *rveBoundaryNodeIdMap["y-"])
      {
        if (nodeYm != rve_ref_node_map_["N1B"]->Id() && nodeYm != rve_ref_node_map_["N2"]->Id())
        {
          PBC.push_back(discret_ptr_->gNode(FindOppositeEdgeNode(
              nodeYm, INPAR::RVE_MPC::RveEdgeIdentifiers::Gamma_ym, rveBoundaryNodeIdMap)));
          PBC.push_back(discret_ptr_->gNode(nodeYm));
          PBC.push_back(rve_ref_node_map_["N4"]);
          PBC.push_back(rve_ref_node_map_["N1B"]);
          PBCs.push_back(PBC);
          PBC.clear();
        }
      }
      break;
    }
    default:
      FOUR_C_THROW("No ref def type defined");
  }

  // Ensure no constraint is enforced twice:
  int indx = 0;
  std::vector<int> ids;
  std::map<int, std::vector<int>> idListSet;

  for (const auto& pbc : PBCs)
  {
    for (auto* node : pbc)
    {
      ids.push_back(node->Id());
    }
    std::sort(ids.begin(), ids.end());
    idListSet[indx++] = ids;
    ids.clear();
  }


  // Print
  int nr = 0;
  IO::cout(IO::debug) << "Sorted Pair Ids:" << IO::endl;
  for (const auto& couple : idListSet)
  {
    IO::cout(IO::debug) << "Pair " << nr++ << ":";
    for (auto id : couple.second)
    {
      IO::cout(IO::debug) << id << ", ";
    }
    IO::cout(IO::debug) << IO::endl;
  }

  std::vector<int> idsToRemove;
  for (const auto& entryA : idListSet)
  {
    if (std::find(idsToRemove.begin(), idsToRemove.end(), entryA.first) == idsToRemove.end())
    {
      for (const auto& entryB : idListSet)
      {
        if (entryA.second == entryB.second && entryA.first != entryB.first)
        {
          idsToRemove.push_back(entryB.first);
          IO::cout(IO::debug) << "remove pair: " << entryB.first << IO::endl;
        }
      }
    }
  }

  // Remove duplicate constraints
  std::sort(idsToRemove.rbegin(), idsToRemove.rend());
  for (int id : idsToRemove)
  {
    PBCs.erase(PBCs.begin() + id);
  }
  IO::cout(IO::verbose) << "All Node Pairs found. Following Nodes are coupled:" << IO::endl;

  // Create the vector of MPC "Elements
  int mpcId = 0;
  std::vector<int> pbcDofs;
  std::vector<double> pbcCoefs = {1., -1., -1., 1.};

  int nDofCoupled = 3;
  if (rve_dim_ == INPAR::RVE_MPC::rve2d)
  {
    nDofCoupled = 2;
  }


  for (const auto& pbc : PBCs)
  {
    IO::cout(IO::debug) << "Node Pair: ";
    for (int dim = 0; dim < nDofCoupled; ++dim)
    {
      pbcDofs.clear();
      for (auto* node : pbc)
      {
        IO::cout(IO::debug) << node->Id() << ", ";

        // Create coupled equation dof list:
        pbcDofs.emplace_back(discret_ptr_->Dof(node)[dim]);
      }
      IO::cout(IO::debug) << "\ndofs used in linear coupled equation: \n ";
      for (auto dof : pbcDofs)
      {
        IO::cout(IO::debug) << dof << ", ";
      }
      IO::cout(IO::debug) << IO::endl;
      listMPCs_.emplace_back(Teuchos::rcp(new LinearCoupledEquation(mpcId++, pbcDofs, pbcCoefs)));
    }
    IO::cout(IO::debug) << "\n";
  }
  IO::cout(IO::verbose) << IO::endl
                        << "Total number of node pairs created for periodic boundary conditions: "
                        << listMPCs_.size() << IO::endl;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::BuildLinearMPCs()
{
  IO::cout(IO::verbose) << IO::endl << "-------------------------------------------" << IO::endl;
  IO::cout(IO::verbose) << "Reading linear coupled eq. from .dat file" << IO::endl;
  IO::cout(IO::verbose) << "linear MPC codition count: "
                        << point_linear_coupled_equation_conditions_.size() << IO::endl;


  int nEq = 0;
  for (const auto& ceTerm : point_linear_coupled_equation_conditions_)
  {
    nEq = std::max(nEq, (ceTerm->Get<int>("EQUATION_ID")) - 1);
  }
  IO::cout(IO::verbose) << "There are " << nEq + 1 << " linear MPC Equations defined" << IO::endl;

  int dofPos;
  int cond_num = 0;
  std::vector<std::vector<int>> constraintRowIds(nEq + 1);
  std::vector<std::vector<int>> constraintColIds(nEq + 1);
  std::vector<std::vector<double>> constraintCoeffs(nEq + 1);


  for (const auto& ceTerm : point_linear_coupled_equation_conditions_)
  {
    auto eq_id = (ceTerm->Get<int>("EQUATION_ID")) - 1;
    const auto* node_id = ceTerm->GetNodes();
    const auto& dofStr = ceTerm->IO::InputParameterContainer::Get<std::string>("DOF");
    auto coef = ceTerm->Get<double>("COEFFICIENT");
    auto* node = discret_ptr_->gNode(node_id->data()[0]);


    if (dofStr == "dispx")
    {
      dofPos = 0;
    }
    else if (dofStr == "dispy")
    {
      dofPos = 1;
    }
    else
    {
      FOUR_C_THROW(
          "Linear coupled equations (MPCs) are only implemented for 2D (dispx or dispy DOFs)");
    }
    auto dofID = discret_ptr_->Dof(node)[dofPos];


    IO::cout(IO::debug) << IO::endl;
    IO::cout(IO::debug) << "Condition Number " << cond_num++ << ": " << IO::endl;
    IO::cout(IO::debug) << "Eq.Id: " << eq_id << IO::endl;
    IO::cout(IO::debug) << "Node Id: " << node_id->data()[0] << IO::endl;
    IO::cout(IO::debug) << "Disp String: " << dofStr.c_str() << IO::endl;
    IO::cout(IO::debug) << "DOF ID: " << dofID << IO::endl;
    IO::cout(IO::debug) << "COEF: " << coef << IO::endl << IO::endl;

    // Save the linear MPCs
    constraintRowIds[eq_id].push_back(eq_id);
    constraintColIds[eq_id].push_back(dofID);
    constraintCoeffs[eq_id].push_back(coef);
    IO::cout(IO::debug) << "Added Term Equation with ID: " << eq_id << IO::endl;
    IO::cout(IO::debug) << "Current SIze constraintColIDs" << constraintColIds.size() << IO::endl;
  }
  // Get number of MPC already in the MPCs List
  int nMPC = 0;
  for (const auto& mpc : listMPCs_)
  {
    nMPC += mpc->GetNumberOfMPCs();
  }
  unsigned int i = 0;
  for (; i < constraintRowIds.size(); ++i)
  {
    listMPCs_.emplace_back(
        Teuchos::rcp(new LinearCoupledEquation(nMPC++, constraintColIds[i], constraintCoeffs[i])));

    IO::cout(IO::verbose) << "Linear MPC #" << i << "  Created: 0 = ";
    for (unsigned int o = 0; o < constraintColIds[i].size(); ++o)
    {
      IO::cout(IO::verbose) << " +" << constraintCoeffs[i][o] << "*d" << constraintColIds[i][o];
    }
    IO::cout(IO::verbose) << IO::endl;
  }
  IO::cout(IO::verbose) << "Number of Linear MPCs Created: " << i + 1 << IO::endl;
  IO::cout(IO::verbose) << "Number of Elements in the listMPC: " << nMPC << IO::endl;

  return i + 1;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::FindPeriodicRveCornerNodes(
    const std::vector<int>* edge1, const std::vector<int>* edge2)
{
  for (int nodeId : *edge2)
  {
    if (std::find(edge1->begin(), edge1->end(), nodeId) != edge1->end()) return nodeId;
  }
  return -1;
}

int CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::FindPeriodicRveCornerNodes(
    const std::vector<int>* surf1, const std::vector<int>* surf2, const std::vector<int>* surf3)
{
  std::vector<int> commonNodeIds12, commonNodeIds123;
  std::set_intersection(surf1->begin(), surf1->end(), surf2->begin(), surf2->end(),
      std::back_inserter(commonNodeIds12));

  std::set_intersection(surf3->begin(), surf3->end(), commonNodeIds12.begin(),
      commonNodeIds12.end(), std::back_inserter(commonNodeIds123));

  if (commonNodeIds123.size() < 1)
  {
    FOUR_C_THROW("No common node found");
  }

  else if (commonNodeIds123.size() > 1)
  {
    FOUR_C_THROW("More than one common node found");
  }
  return commonNodeIds123[0];
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::BuildPeriodicRveCornerNodeMap(
    std::map<std::string, const std::vector<int>*>& rveBoundaryNodeIdMap,
    std::map<std::string, int>& rveCornerNodeIdMap)
{
  switch (rve_dim_)
  {
    case INPAR::RVE_MPC::RveDimension::rve2d:
    {
      //* Get the Corner Node Ids */
      /*              N4 -- N
       *              |     |
       *      N   --- |     |---   N
       *      |                    |
       *    x-|                    |x+
       *      N1L ---|       | ----N2
       *             |       |
       *             N1B --  N
       *
       *      N1 ------ y- ------ N2
       */

      //* Get the Corner Node Ids */

      /*
       *      N4 ------y+------ N3
       *      |                    |
       *      x-                   |x+
       *      |                   |
       *      N1 ------ y- ------ N2
       */

      rveCornerNodeIdMap["N1"] =
          FindPeriodicRveCornerNodes(rveBoundaryNodeIdMap["x-"], rveBoundaryNodeIdMap["y-"]);
      rveCornerNodeIdMap["N2"] =
          FindPeriodicRveCornerNodes(rveBoundaryNodeIdMap["x+"], rveBoundaryNodeIdMap["y-"]);
      rveCornerNodeIdMap["N3"] =
          FindPeriodicRveCornerNodes(rveBoundaryNodeIdMap["x+"], rveBoundaryNodeIdMap["y+"]);
      rveCornerNodeIdMap["N4"] =
          FindPeriodicRveCornerNodes(rveBoundaryNodeIdMap["x-"], rveBoundaryNodeIdMap["y+"]);


      IO::cout(IO::verbose) << "\nAutomatically determined following reference Nodes: " << IO::endl;
      IO::cout(IO::verbose) << "---------------------------------------------------" << IO::endl;
      IO::cout(IO::verbose) << "N1: " << rveCornerNodeIdMap["N1"] << ", ";
      IO::cout(IO::verbose) << "N2: " << rveCornerNodeIdMap["N2"] << ", ";
      IO::cout(IO::verbose) << "N3: " << rveCornerNodeIdMap["N3"] << ", ";
      IO::cout(IO::verbose) << "N4: " << rveCornerNodeIdMap["N4"] << IO::endl;
    }
    break;
    case INPAR::RVE_MPC::RveDimension::rve3d:
    {
      //  z ^        N8 +  +   +   N7
      //    |      + .            + +
      //    |    +   .          +   +
      //    |  +     .        +     +
      //    N5 +  +  +  +  N3       +
      //    +       .      +        +
      //    +       N4 .  .+   . .. N3
      //    +      .       +     +
      //    +    .         +  +
      //    +  .           +
      //   N1 +  +  +  +  N2 ------> x
      //

      std::vector<std::string> boundaryNames = {"x+", "x-", "y+", "y-", "z+", "z-"};
      rveCornerNodeIdMap["N1"] = FindPeriodicRveCornerNodes(
          rveBoundaryNodeIdMap["x-"], rveBoundaryNodeIdMap["z-"], rveBoundaryNodeIdMap["y-"]);

      rveCornerNodeIdMap["N2"] = FindPeriodicRveCornerNodes(
          rveBoundaryNodeIdMap["x+"], rveBoundaryNodeIdMap["z-"], rveBoundaryNodeIdMap["y-"]);

      rveCornerNodeIdMap["N4"] = FindPeriodicRveCornerNodes(
          rveBoundaryNodeIdMap["x-"], rveBoundaryNodeIdMap["z-"], rveBoundaryNodeIdMap["y+"]);

      rveCornerNodeIdMap["N5"] = FindPeriodicRveCornerNodes(
          rveBoundaryNodeIdMap["x-"], rveBoundaryNodeIdMap["z+"], rveBoundaryNodeIdMap["y-"]);

      IO::cout(IO::verbose) << "\nAutomatically determined following reference Nodes: " << IO::endl;
      IO::cout(IO::verbose) << "---------------------------------------------------" << IO::endl;
      IO::cout(IO::verbose) << "N1: " << rveCornerNodeIdMap["N1"] << ", ";
      IO::cout(IO::verbose) << "N2: " << rveCornerNodeIdMap["N2"] << ", ";
      IO::cout(IO::verbose) << "N4: " << rveCornerNodeIdMap["N4"] << ", ";
      IO::cout(IO::verbose) << "N5: " << rveCornerNodeIdMap["N5"] << IO::endl;
      break;
    }
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::
    BuildPeriodicRveBoundaryNodeMap(
        std::map<std::string, const std::vector<int>*>& rveBoundaryNodeIdMap)
{
  switch (rve_dim_)
  {
    case INPAR::RVE_MPC::RveDimension::rve2d:
    {
      discret_ptr_->GetCondition("LinePeriodicRve", line_periodic_rve_conditions_);



      IO::cout(IO::verbose)
          << IO::endl
          << "Reading Line Conditions:  " << IO::endl
          << "+--------------------------------------------------------------------+" << IO::endl;
      for (const auto& conditionLine : line_periodic_rve_conditions_)
      {
        const auto& boundary =
            conditionLine->IO::InputParameterContainer::Get<std::string>("EdgeLineId");

        // Print the Edge Condition
        IO::cout(IO::verbose) << "EDGE: " << boundary.c_str() << " Node IDs: ";
        for (auto nodeId : *conditionLine->GetNodes())
        {
          IO::cout(IO::verbose) << nodeId << " ";
        }
        IO::cout(IO::verbose) << IO::endl;

        // Create EdgeNodeMap
        rveBoundaryNodeIdMap[boundary.c_str()] = conditionLine->GetNodes();
      }
    }
    break;

    case INPAR::RVE_MPC::RveDimension::rve3d:
    {
      discret_ptr_->GetCondition("SurfacePeriodicRve", surface_periodic_rve_conditions_);

      IO::cout(IO::verbose)
          << IO::endl
          << "Reading Surface Conditions:  " << IO::endl
          << "+--------------------------------------------------------------------+" << IO::endl;
      for (const auto& conditionLine : surface_periodic_rve_conditions_)
      {
        const auto& boundary =
            conditionLine->IO::InputParameterContainer::Get<std::string>("SurfId");

        IO::cout(IO::verbose) << "SURF: " << boundary.c_str() << " Node IDs: ";
        for (auto nodeId : *conditionLine->GetNodes())
        {
          IO::cout(IO::verbose) << nodeId << " ";
        }
        IO::cout(IO::verbose) << IO::endl;
        rveBoundaryNodeIdMap[boundary.c_str()] = conditionLine->GetNodes();
      }
      break;
    }
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::Reset() {}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE
