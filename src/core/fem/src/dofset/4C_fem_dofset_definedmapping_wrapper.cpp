/*----------------------------------------------------------------------*/
/*! \file

 \brief  A dofset that does not rely on same GID/LID numbers but uses
         a defined node mapping instead (not implemented for element DOFs).

   \level 3


*/
/*----------------------------------------------------------------------*/


#include "4C_fem_dofset_definedmapping_wrapper.hpp"

#include "4C_coupling_matchingoctree.hpp"
#include "4C_fem_condition_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_base.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::DOFSets::DofSetDefinedMappingWrapper::DofSetDefinedMappingWrapper(
    Teuchos::RCP<DofSetInterface> sourcedofset,
    Teuchos::RCP<const Core::FE::Discretization> sourcedis, const std::string& couplingcond,
    const std::set<int> condids)
    : DofSetBase(),
      sourcedofset_(sourcedofset),
      targetlidtosourcegidmapping_(Teuchos::null),
      sourcedis_(sourcedis),
      couplingcond_(couplingcond),
      condids_(condids),
      filled_(false)
{
  if (sourcedofset_ == Teuchos::null) FOUR_C_THROW("Source dof set is null pointer.");
  if (sourcedis_ == Teuchos::null) FOUR_C_THROW("Source discretization is null pointer.");

  sourcedofset_->register_proxy(this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::DOFSets::DofSetDefinedMappingWrapper::~DofSetDefinedMappingWrapper()
{
  if (sourcedofset_ != Teuchos::null) sourcedofset_->unregister(this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::DOFSets::DofSetDefinedMappingWrapper::assign_degrees_of_freedom(
    const Core::FE::Discretization& dis, const unsigned dspos, const int start)
{
  if (sourcedofset_ == Teuchos::null)
    FOUR_C_THROW("No source dof set assigned to mapping dof set!");
  if (sourcedis_ == Teuchos::null)
    FOUR_C_THROW("No source discretization assigned to mapping dof set!");

  // get condition which defines the coupling on target discretization
  std::vector<Core::Conditions::Condition*> conds;
  dis.get_condition(couplingcond_, conds);

  // get condition which defines the coupling on source discretization
  std::vector<Core::Conditions::Condition*> conds_source;
  sourcedis_->get_condition(couplingcond_, conds_source);

  // get the respective nodes which are in the condition
  const bool use_coupling_id = condids_.size() != 1;
  std::map<int, Teuchos::RCP<std::vector<int>>> nodes;
  Core::Conditions::FindConditionedNodes(dis, conds, nodes, use_coupling_id);
  std::map<int, Teuchos::RCP<std::vector<int>>> nodes_source;
  Core::Conditions::FindConditionedNodes(*sourcedis_, conds_source, nodes_source, use_coupling_id);

  // map that will be filled with coupled nodes
  // mapping: target node gid to (source node gid, distance)
  std::map<int, std::pair<int, double>> coupling;

  // define iterators
  std::map<int, Teuchos::RCP<std::vector<int>>>::iterator iter_target;
  std::map<int, Teuchos::RCP<std::vector<int>>>::iterator iter_source;

  for (std::set<int>::iterator it = condids_.begin(); it != condids_.end(); ++it)
  {
    // find corresponding condition on source discretization
    iter_target = nodes.find(*it);
    // find corresponding condition on source discretization
    iter_source = nodes_source.find(*it);

    // get the nodes
    std::vector<int> sourcenodes;
    std::vector<int> targetnodes;
    if (iter_source != nodes_source.end()) sourcenodes = *iter_source->second;
    if (iter_target != nodes.end()) targetnodes = *iter_target->second;

    // initialize search tree for search
    Core::COUPLING::NodeMatchingOctree nodematchingtree;
    nodematchingtree.init(dis, targetnodes, 150, 1e-08);
    nodematchingtree.setup();

    // map that will be filled with coupled nodes for this condition
    // mapping: target node gid to (source node gid, distance)
    // note: FindMatch loops over all SOURCE (i.e. slave) nodes
    //       and finds corresponding target nodes.
    std::map<int, std::pair<int, double>> condcoupling;
    // match target and source nodes using octtree
    nodematchingtree.find_match(*sourcedis_, sourcenodes, condcoupling);

    // check if all nodes where matched for this condition ID
    if (targetnodes.size() != condcoupling.size())
      FOUR_C_THROW(
          "Did not get unique target to source spatial node coordinate mapping.\n"
          "targetnodes.size()=%d, coupling.size()=%d.\n"
          "The heterogeneous reaction strategy requires matching source and target meshes!",
          targetnodes.size(), condcoupling.size());

    // insert found coupling of this condition ID into map of all match nodes
    coupling.insert(condcoupling.begin(), condcoupling.end());

  }  // loop over all condition ids

  // clone communicator of target discretization
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(dis.get_comm().Clone());

  // extract permutation
  std::vector<int> targetnodes(dis.node_row_map()->MyGlobalElements(),
      dis.node_row_map()->MyGlobalElements() + dis.node_row_map()->NumMyElements());

  std::vector<int> patchedtargetnodes;
  patchedtargetnodes.reserve(coupling.size());
  std::vector<int> permsourcenodes;
  permsourcenodes.reserve(coupling.size());

  for (unsigned i = 0; i < targetnodes.size(); ++i)
  {
    const int gid = targetnodes[i];

    // We allow to hand in target nodes that do not take part in the
    // coupling. If this is undesired behaviour the user has to make
    // sure all nodes were used.
    if (coupling.find(gid) != coupling.end())
    {
      std::pair<int, double>& coupled = coupling[gid];
      patchedtargetnodes.push_back(gid);
      permsourcenodes.push_back(coupled.first);
    }
  }

  // Epetra maps
  Teuchos::RCP<Epetra_Map> targetnodemap = Teuchos::rcp(
      new Epetra_Map(-1, patchedtargetnodes.size(), patchedtargetnodes.data(), 0, *com));

  Teuchos::RCP<Epetra_Map> permsourcenodemap =
      Teuchos::rcp(new Epetra_Map(-1, permsourcenodes.size(), permsourcenodes.data(), 0, *com));

  // we expect to get maps of exactly the same shape
  if (not targetnodemap->PointSameAs(*permsourcenodemap))
    FOUR_C_THROW("target and permuted source node maps do not match");

  // export target nodes to source node distribution
  Teuchos::RCP<Epetra_IntVector> permsourcenodevec = Teuchos::rcp(
      new Epetra_IntVector(Copy, *targetnodemap, permsourcenodemap->MyGlobalElements()));

  // initialize the final mapping
  targetlidtosourcegidmapping_ = Teuchos::rcp(new Epetra_IntVector(*dis.node_col_map()));

  // default value -1
  targetlidtosourcegidmapping_->PutValue(-1);

  // export to column map
  Core::LinAlg::export_to(*permsourcenodevec, *targetlidtosourcegidmapping_);

  // filled.
  filled_ = true;

  // tell the proxies
  notify_assigned();

  return start;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetDefinedMappingWrapper::reset()
{
  targetlidtosourcegidmapping_ = Teuchos::null;
  filled_ = false;

  // tell the proxies
  notify_reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetDefinedMappingWrapper::disconnect(DofSetInterface* dofset)
{
  if (dofset == sourcedofset_.get())
  {
    sourcedofset_ = Teuchos::null;
    sourcedis_ = Teuchos::null;
  }
  else
    FOUR_C_THROW("cannot disconnect from non-connected DofSet");

  // clear my Teuchos::rcps.
  reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Core::Nodes::Node* Core::DOFSets::DofSetDefinedMappingWrapper::get_source_node(
    int targetLid) const
{
  // check
  FOUR_C_ASSERT(targetLid <= targetlidtosourcegidmapping_->MyLength(), "Target Lid out of range!");

  // get the gid of the source node
  int sourcegid = (*targetlidtosourcegidmapping_)[targetLid];

  // the target is not mapped -> return null pointer
  if (sourcegid == -1) return nullptr;
  // get the node from the source discretization
  return sourcedis_->g_node(sourcegid);
}

FOUR_C_NAMESPACE_CLOSE
