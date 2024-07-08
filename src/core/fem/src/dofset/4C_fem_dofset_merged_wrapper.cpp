/*---------------------------------------------------------------------*/
/*! \file

 \brief A proxy of a dofset that adds additional, existing degrees of freedom from the same
        discretization to nodes (not implemented for element DOFs).

 \level 2


*/
/*---------------------------------------------------------------------*/


#include "4C_fem_dofset_merged_wrapper.hpp"

#include "4C_coupling_matchingoctree.hpp"
#include "4C_fem_condition_utils.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <Epetra_Export.h>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::DOFSets::DofSetMergedWrapper::DofSetMergedWrapper(Teuchos::RCP<DofSetInterface> sourcedofset,
    Teuchos::RCP<const Core::FE::Discretization> sourcedis, const std::string& couplingcond_master,
    const std::string& couplingcond_slave)
    : DofSetBase(),
      master_nodegids_col_layout_(Teuchos::null),
      sourcedofset_(sourcedofset),
      sourcedis_(sourcedis),
      couplingcond_master_(couplingcond_master),
      couplingcond_slave_(couplingcond_slave),
      filled_(false)
{
  if (sourcedofset_ == Teuchos::null) FOUR_C_THROW("Source dof set is null pointer.");
  if (sourcedis_ == Teuchos::null) FOUR_C_THROW("Source discretization is null pointer.");

  sourcedofset_->register_proxy(this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::DOFSets::DofSetMergedWrapper::~DofSetMergedWrapper()
{
  if (sourcedofset_ != Teuchos::null) sourcedofset_->unregister(this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* Core::DOFSets::DofSetMergedWrapper::dof_row_map() const
{
  // the merged dofset does not add new dofs. So we can just return the
  // originial dof map here.
  return sourcedofset_->dof_row_map();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* Core::DOFSets::DofSetMergedWrapper::dof_col_map() const
{
  // the merged dofset does not add new dofs. So we can just return the
  // originial dof map here.
  return sourcedofset_->dof_col_map();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::DOFSets::DofSetMergedWrapper::assign_degrees_of_freedom(
    const Core::FE::Discretization& dis, const unsigned dspos, const int start)
{
  if (sourcedofset_ == Teuchos::null) FOUR_C_THROW("No source dof set assigned to merged dof set!");
  if (sourcedis_ == Teuchos::null)
    FOUR_C_THROW("No source discretization assigned to mapping dof set!");

  // get nodes to be coupled
  std::vector<int> masternodes;
  Core::Conditions::FindConditionedNodes(*sourcedis_, couplingcond_master_, masternodes);
  std::vector<int> slavenodes;
  Core::Conditions::FindConditionedNodes(dis, couplingcond_slave_, slavenodes);

  // initialize search tree
  auto tree = Core::COUPLING::NodeMatchingOctree();
  tree.init(*sourcedis_, masternodes, 150);
  tree.setup();

  // match master and slave nodes using octtree
  // master id -> slave id, distance
  std::map<int, std::pair<int, double>> coupling;
  tree.find_match(dis, slavenodes, coupling);

  // all nodes should be coupled
  if (masternodes.size() != coupling.size())
    FOUR_C_THROW(
        "Did not get 1:1 correspondence. \nmasternodes.size()=%d, coupling.size()=%d."
        "DofSetMergedWrapper requires matching slave and master meshes!",
        masternodes.size(), coupling.size());

  // initialize final mapping
  Teuchos::RCP<Epetra_IntVector> my_master_nodegids_row_layout =
      Teuchos::rcp(new Epetra_IntVector(*dis.node_row_map()));

  // loop over all coupled nodes
  for (unsigned i = 0; i < masternodes.size(); ++i)
  {
    // get master gid
    int gid = masternodes[i];

    // We allow to hand in master nodes that do not take part in the
    // coupling. If this is undesired behaviour the user has to make
    // sure all nodes were used.
    if (coupling.find(gid) != coupling.end())
    {
      std::pair<int, double>& coupled = coupling[gid];
      int slavegid = coupled.first;
      int slavelid = dis.node_row_map()->LID(slavegid);
      if (slavelid == -1) FOUR_C_THROW("slave gid %d was not found on this proc", slavegid);

      // save master gid at col lid of corresponding slave node
      (*my_master_nodegids_row_layout)[slavelid] = gid;
    }
  }

  // initialize final mapping
  master_nodegids_col_layout_ = Teuchos::rcp(new Epetra_IntVector(*dis.node_col_map()));

  // export to column map
  Core::LinAlg::Export(*my_master_nodegids_row_layout, *master_nodegids_col_layout_);


  ////////////////////////////////////////////////////
  // now we match source slave and aux dis
  ////////////////////////////////////////////////////

  // get nodes to be coupled
  masternodes.clear();
  Core::Conditions::FindConditionedNodes(*sourcedis_, couplingcond_slave_, masternodes);
  slavenodes.clear();
  Core::Conditions::FindConditionedNodes(dis, couplingcond_slave_, slavenodes);

  // initialize search tree
  tree.init(*sourcedis_, masternodes, 150);
  tree.setup();

  // match master and slave nodes using octtree
  // master id -> slave id, distance
  coupling.clear();
  tree.find_match(dis, slavenodes, coupling);

  // all nodes should be coupled
  if (masternodes.size() != coupling.size())
    FOUR_C_THROW(
        "Did not get 1:1 correspondence. \nmasternodes.size()=%d, coupling.size()=%d."
        "DofSetMergedWrapper requires matching slave and master meshes!",
        masternodes.size(), coupling.size());

  // initialize final mapping
  Teuchos::RCP<Epetra_IntVector> my_slave_nodegids_row_layout =
      Teuchos::rcp(new Epetra_IntVector(*dis.node_row_map()));

  // loop over all coupled nodes
  for (unsigned i = 0; i < masternodes.size(); ++i)
  {
    // get master gid
    int gid = masternodes[i];

    // We allow to hand in master nodes that do not take part in the
    // coupling. If this is undesired behaviour the user has to make
    // sure all nodes were used.
    if (coupling.find(gid) != coupling.end())
    {
      std::pair<int, double>& coupled = coupling[gid];
      int slavegid = coupled.first;
      int slavelid = dis.node_row_map()->LID(slavegid);
      if (slavelid == -1) FOUR_C_THROW("slave gid %d was not found on this proc", slavegid);

      // save master gid at col lid of corresponding slave node
      (*my_slave_nodegids_row_layout)[slavelid] = gid;
    }
  }

  // initialize final mapping
  slave_nodegids_col_layout_ = Teuchos::rcp(new Epetra_IntVector(*dis.node_col_map()));

  // export to column map
  Core::LinAlg::Export(*my_slave_nodegids_row_layout, *slave_nodegids_col_layout_);

  // set filled flag true
  filled_ = true;

  // tell the proxies
  notify_assigned();

  return start;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetMergedWrapper::reset()
{
  master_nodegids_col_layout_ = Teuchos::null;
  slave_nodegids_col_layout_ = Teuchos::null;

  // set filled flag
  filled_ = false;

  notify_reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetMergedWrapper::disconnect(DofSetInterface* dofset)
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

FOUR_C_NAMESPACE_CLOSE
