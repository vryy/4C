/*-----------------------------------------------------------------------*/
/*! \file
\brief A class for a frictional contact node
\level 2
*/
/*-----------------------------------------------------------------------*/

#include "4C_contact_rough_node.hpp"

#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_global_data.hpp"
#include "4C_utils_function.hpp"

#ifdef FOUR_C_WITH_MIRCO

#include <mirco_topology.h>
#include <mirco_topologyutilities.h>

#endif

FOUR_C_NAMESPACE_OPEN

CONTACT::RoughNodeType CONTACT::RoughNodeType::instance_;

Core::Communication::ParObject* CONTACT::RoughNodeType::Create(const std::vector<char>& data)
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
    bool randomtopologyflag, bool randomseedflag, int randomgeneratorseed)
    : CONTACT::Node(id, coords, owner, dofs, isslave, initactive),
      hurstexponentfunction_(hurstexponentfunction),
      initialtopologystddeviationfunction_(initialtopologystddeviationfunction),
      resolution_(resolution),
      randomtopologyflag_(randomtopologyflag),
      randomseedflag_(randomseedflag),
      randomgeneratorseed_(randomgeneratorseed)
{
#ifdef FOUR_C_WITH_MIRCO
  if (isslave)
  {
    hurstExponent_ =
        Global::Problem::Instance()
            ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(hurstexponentfunction_ - 1)
            .Evaluate(this->X().data(), 1, this->Dim());
    initialTopologyStdDeviation_ = Global::Problem::Instance()
                                       ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(
                                           initialtopologystddeviationfunction_ - 1)
                                       .Evaluate(this->X().data(), 1, this->Dim());

    const int N = pow(2, resolution_);
    topology_.shape(N + 1, N + 1);

    std::string topologyFilePath = "";
    Teuchos::RCP<MIRCO::TopologyGeneration> surfacegenerator;
    // creating the correct surface object
    MIRCO::CreateSurfaceObject(resolution_, initialTopologyStdDeviation_, hurstExponent_,
        randomseedflag_, topologyFilePath, randomtopologyflag_, randomgeneratorseed_,
        surfacegenerator);
    surfacegenerator->GetSurface(topology_);

    auto max_and_mean = MIRCO::ComputeMaxAndMean(topology_);
    maxTopologyHeight_ = max_and_mean.max_;
  }
#else
  FOUR_C_THROW(
      "You are trying to create a RoughNode with FOUR_C_WITH_MIRCO flag turned off. Please enable "
      "this flag and build 4C again");
#endif
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 *----------------------------------------------------------------------*/
void CONTACT::RoughNode::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class Mortar::Node
  CONTACT::Node::Pack(data);

  AddtoPack(data, hurstexponentfunction_);
  AddtoPack(data, initialtopologystddeviationfunction_);
  AddtoPack(data, resolution_);
  AddtoPack(data, randomtopologyflag_);
  AddtoPack(data, randomseedflag_);
  AddtoPack(data, randomgeneratorseed_);

  AddtoPack(data, hurstExponent_);
  AddtoPack(data, initialTopologyStdDeviation_);
  AddtoPack(data, topology_);
  AddtoPack(data, maxTopologyHeight_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 *----------------------------------------------------------------------*/
void CONTACT::RoughNode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class CONTACT::Node
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  CONTACT::Node::Unpack(basedata);

  hurstexponentfunction_ = ExtractInt(position, data);
  initialtopologystddeviationfunction_ = ExtractInt(position, data);
  resolution_ = ExtractInt(position, data);
  randomtopologyflag_ = ExtractInt(position, data);
  randomseedflag_ = ExtractInt(position, data);
  randomgeneratorseed_ = ExtractInt(position, data);

  hurstExponent_ = ExtractDouble(position, data);
  initialTopologyStdDeviation_ = ExtractDouble(position, data);
  ExtractfromPack(position, data, topology_);
  maxTopologyHeight_ = ExtractDouble(position, data);

  // Check
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

FOUR_C_NAMESPACE_CLOSE
