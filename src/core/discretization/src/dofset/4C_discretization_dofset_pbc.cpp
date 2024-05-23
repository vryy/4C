/*---------------------------------------------------------------------*/
/*! \file

\brief A modified set of degrees of freedom for periodic boundary
       conditions

\level 0


*/
/*---------------------------------------------------------------------*/

#include "4C_discretization_dofset_pbc.hpp"

#include "4C_lib_discret.hpp"
#include "4C_lib_element.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/07|
 *----------------------------------------------------------------------*/
CORE::Dofsets::PBCDofSet::PBCDofSet(Teuchos::RCP<std::map<int, std::vector<int>>> couplednodes)
    : DofSet(), perbndcouples_(Teuchos::null), myMaxGID_(-1)
{
  SetCoupledNodes(couplednodes);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::Dofsets::PBCDofSet::MaxAllGID() const { return myMaxGID_; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::Dofsets::PBCDofSet::MinAllGID() const { return myMinGID_; }


int CORE::Dofsets::PBCDofSet::assign_degrees_of_freedom(
    const DRT::Discretization& dis, const unsigned dspos, const int start)
{
  // temporarily store the slave node set
  Teuchos::RCP<std::set<int>> tempset = slavenodeids_;
  slavenodeids_ = Teuchos::rcp(new std::set<int>);

  // assign dofs using the empty slave node set. This way the dofrowmap_
  // contains exactly the entries as in a regular dofset
  DofSet::assign_degrees_of_freedom(dis, dspos, start);
  if (pccdofhandling_)
    FOUR_C_THROW("ERROR: Point coupling cinditions not yet implemented for PBCDofSet");

  myMaxGID_ = DofSet::MaxAllGID();
  myMinGID_ = DofSet::MinAllGID();

  // restore the slave node set
  slavenodeids_ = tempset;

  // assign dofs for the standard dofset, that is without periodic boundary
  // conditions and with the slave node set back in place
  int count = DofSet::assign_degrees_of_freedom(dis, dspos, start);


  // loop all master nodes and set the dofs of the slaves to the dofs of the master
  // remark: the previously assigned dofs of slave nodes are overwritten here
  for (std::map<int, std::vector<int>>::iterator master = perbndcouples_->begin();
       master != perbndcouples_->end(); ++master)
  {
    int master_lid = dis.NodeColMap()->LID(master->first);

    if (master_lid < 0)
    {
      FOUR_C_THROW("master gid %d not on proc %d, but required by slave %d", master->first,
          dis.Comm().MyPID(), master->second[0]);
    }

    for (std::vector<int>::iterator slave = master->second.begin(); slave != master->second.end();
         ++slave)
    {
      int slave_lid = dis.NodeColMap()->LID(*slave);

      if (slave_lid > -1)
      {
        (*numdfcolnodes_)[slave_lid] = (*numdfcolnodes_)[master_lid];
        (*idxcolnodes_)[slave_lid] = (*idxcolnodes_)[master_lid];
      }
      else
      {
#ifdef FOUR_C_ENABLE_ASSERTIONS
        if (dis.NodeRowMap()->MyGID(master->first))
        {
          FOUR_C_THROW("slave not on proc but master owned by proc\n");
        }
#endif
      }
    }
  }

  return count;
}


/*----------------------------------------------------------------------*
 |  update coupled nodes map                             rasthofer 07/11|
 |                                                       DA wichmann    |
 *----------------------------------------------------------------------*/
void CORE::Dofsets::PBCDofSet::SetCoupledNodes(
    Teuchos::RCP<std::map<int, std::vector<int>>> couplednodes)
{
  perbndcouples_ = couplednodes;
  slavenodeids_ = Teuchos::rcp(new std::set<int>);

  for (std::map<int, std::vector<int>>::iterator curr = perbndcouples_->begin();
       curr != perbndcouples_->end(); ++curr)
  {
    std::vector<int>& sids = curr->second;
    std::copy(sids.begin(), sids.end(), std::inserter(*slavenodeids_, slavenodeids_->begin()));
  }

  /// Build the connectivity between slave node and its master node
  build_slave_to_master_node_connectivity();

  return;
}

/*----------------------------------------------------------------------*
 |  Build the connectivity between slave node and its master node       |
 |                                                       schott 05/15   |
 *----------------------------------------------------------------------*/
void CORE::Dofsets::PBCDofSet::build_slave_to_master_node_connectivity()
{
  perbnd_slavetomaster_ = Teuchos::rcp(new std::map<int, int>);

  for (std::map<int, std::vector<int>>::const_iterator masterslavepair = perbndcouples_->begin();
       masterslavepair != perbndcouples_->end(); ++masterslavepair)
  {
    // loop slave nodes associated with master
    for (std::vector<int>::const_iterator iter = masterslavepair->second.begin();
         iter != masterslavepair->second.end(); ++iter)
    {
      const int slavegid = *iter;
      (*perbnd_slavetomaster_)[slavegid] = masterslavepair->first;
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
