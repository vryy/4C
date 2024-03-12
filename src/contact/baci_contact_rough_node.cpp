/*-----------------------------------------------------------------------*/
/*! \file
\brief A class for a frictional contact node
\level 2
*/
/*-----------------------------------------------------------------------*/

#include "baci_contact_rough_node.hpp"

#include "baci_contact_defines.hpp"
#include "baci_contact_element.hpp"
#include "baci_global_data.hpp"
#include "baci_utils_function.hpp"

#include <mirco_topology.h>
#include <mirco_topologyutilities.h>

BACI_NAMESPACE_OPEN

CONTACT::RoughNodeType CONTACT::RoughNodeType::instance_;

CORE::COMM::ParObject* CONTACT::RoughNodeType::Create(const std::vector<char>& data)
{
  std::vector<double> x(3, 0.0);
  std::vector<int> dofs(0);

  CONTACT::RoughNode* node = new CONTACT::RoughNode(0, x, 0, dofs, false, false, 0, 0, 0, 0, 0, 0);
  node->Unpack(data);

  return node;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::RoughNode::RoughNode(int id, const std::vector<double>& coords, const int owner,
    const std::vector<int>& dofs, const bool isslave, const bool initactive,
    const int hurstexponentfunction, int initialtopologystddeviationfunction, int resolution,
    int randomtopologyflag, int randomseedflag, int randomgeneratorseed)
    : CONTACT::Node(id, coords, owner, dofs, isslave, initactive)
{
  if (isslave && hurstexponentfunction != 0)
  {
    hurstExponent_ = GLOBAL::Problem::Instance()
                         ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(hurstexponentfunction - 1)
                         .Evaluate(this->X().data(), 1, this->Dim());
    initialTopologyStdDeviation_ = GLOBAL::Problem::Instance()
                                       ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(
                                           initialtopologystddeviationfunction - 1)
                                       .Evaluate(this->X().data(), 1, this->Dim());

    const int N = pow(2, resolution);
    topology_ = Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix(N + 1, N + 1));

    std::string topologyFilePath = "";
    Teuchos::RCP<MIRCO::TopologyGeneration> surfacegenerator;
    // creating the correct surface object
    MIRCO::CreateSurfaceObject(resolution, initialTopologyStdDeviation_, hurstExponent_,
        randomseedflag, topologyFilePath, randomtopologyflag, randomgeneratorseed,
        surfacegenerator);
    surfacegenerator->GetSurface(*topology_);

    auto max_and_mean = MIRCO::ComputeMaxAndMean(*topology_);
    maxTopologyHeight_ = max_and_mean.max_;
  }
  return;
}

BACI_NAMESPACE_CLOSE
