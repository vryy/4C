/*----------------------------------------------------------------------*/
/*!
 \file drt_dofset_merged_proxy.cpp

 \brief A proxy of a dofset that adds additional, existing degrees of freedom from the same
        discretization to nodes (not implemented for element DOFs).

 \level 2

 \maintainer Anh-Tu Vuong
             vuong@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15251
 *----------------------------------------------------------------------*/


#include "drt_dofset_merged_proxy.H"

#include "../drt_lib/drt_matchingoctree.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../linalg/linalg_utils.H"

#include <Epetra_Export.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetMergedProxy::DofSetMergedProxy( Teuchos::RCP<DofSet>  dofset,
                                          const Teuchos::RCP<const DRT::Discretization>& sourcedis,
                                          const std::string& couplingcond_master,
                                          const std::string& couplingcond_slave)
  : DofSetProxy(&(*dofset)),
    master_nodegids_col_layout_(Teuchos::null),
    dofset_(dofset),
    sourcedis_(sourcedis),
    couplingcond_master_(couplingcond_master),
    couplingcond_slave_(couplingcond_slave)
{
  NotifyAssigned();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetMergedProxy::AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start)
{
  // make sure the real dofset gets the dofmaps
  DofSetProxy::AssignDegreesOfFreedom(dis,dspos,start);

  // copy communicator
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp( sourcedis_->Comm().Clone());

  // get nodes to be coupled
  std::vector<int> masternodes;
  DRT::UTILS::FindConditionedNodes(*sourcedis_,couplingcond_master_,masternodes);
  std::vector<int> slavenodes;
  DRT::UTILS::FindConditionedNodes(dis,couplingcond_slave_,slavenodes);

  // initialize search tree
  DRT::UTILS::NodeMatchingOctree tree = DRT::UTILS::NodeMatchingOctree();
  tree.Init(*sourcedis_,masternodes,150);
  tree.Setup();

  // match master and slave nodes using octtree
  // master id -> slave id, distance
  std::map<int,std::pair<int,double> > coupling;
  tree.FindMatch(dis, slavenodes, coupling);

  // all nodes should be coupled
  if (masternodes.size() != coupling.size())
    dserror("Did not get 1:1 correspondence. \nmasternodes.size()=%d, coupling.size()=%d."
        "DofSetMergedProxy requires matching slave and master meshes!",
            masternodes.size(), coupling.size());

  // initialize final mapping
  Teuchos::RCP<Epetra_IntVector>  my_master_nodegids_row_layout =
    Teuchos::rcp(new Epetra_IntVector(*dis.NodeRowMap()));

  // loop over all coupled nodes
  for (unsigned i=0; i<masternodes.size(); ++i)
  {
    // get master gid
    int gid = masternodes[i];

    // We allow to hand in master nodes that do not take part in the
    // coupling. If this is undesired behaviour the user has to make
    // sure all nodes were used.
    if (coupling.find(gid) != coupling.end())
    {
      std::pair<int,double>& coupled = coupling[gid];
      int slavegid = coupled.first;
      int slavelid=dis.NodeRowMap()->LID(slavegid);
      if(slavelid==-1)
        dserror("slave gid %d was not found on this proc",slavegid);

      // save master gid at col lid of corresponding slave node
      (*my_master_nodegids_row_layout)[slavelid] = gid;

    }
  }

  // initialize final mapping
  master_nodegids_col_layout_ =
    Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));

  // export to column map
  LINALG::Export(*my_master_nodegids_row_layout,*master_nodegids_col_layout_);


  ////////////////////////////////////////////////////
  // now we match source slave and aux dis
  ////////////////////////////////////////////////////

  // get nodes to be coupled
  masternodes.clear();
  DRT::UTILS::FindConditionedNodes(*sourcedis_,couplingcond_slave_,masternodes);
  slavenodes.clear();
  DRT::UTILS::FindConditionedNodes(dis,couplingcond_slave_,slavenodes);

  // initialize search tree
  tree.Init(*sourcedis_,masternodes,150);
  tree.Setup();

  // match master and slave nodes using octtree
  // master id -> slave id, distance
  coupling.clear();
  tree.FindMatch(dis, slavenodes, coupling);

  // all nodes should be coupled
  if (masternodes.size() != coupling.size())
    dserror("Did not get 1:1 correspondence. \nmasternodes.size()=%d, coupling.size()=%d."
        "DofSetMergedProxy requires matching slave and master meshes!",
            masternodes.size(), coupling.size());

  // initialize final mapping
  Teuchos::RCP<Epetra_IntVector>  my_slave_nodegids_row_layout =
    Teuchos::rcp(new Epetra_IntVector(*dis.NodeRowMap()));

  // loop over all coupled nodes
  for (unsigned i=0; i<masternodes.size(); ++i)
  {
    // get master gid
    int gid = masternodes[i];

    // We allow to hand in master nodes that do not take part in the
    // coupling. If this is undesired behaviour the user has to make
    // sure all nodes were used.
    if (coupling.find(gid) != coupling.end())
    {
      std::pair<int,double>& coupled = coupling[gid];
      int slavegid = coupled.first;
      int slavelid=dis.NodeRowMap()->LID(slavegid);
      if(slavelid==-1)
        dserror("slave gid %d was not found on this proc",slavegid);

      // save master gid at col lid of corresponding slave node
      (*my_slave_nodegids_row_layout)[slavelid] = gid;

    }
  }

  // initialize final mapping
  slave_nodegids_col_layout_ =
    Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));

  // export to column map
  LINALG::Export(*my_slave_nodegids_row_layout,*slave_nodegids_col_layout_);

  return start;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetMergedProxy::NotifyAssigned()
{
  // make sure the real dofset gets the dofmaps
  DofSetProxy::NotifyAssigned();
}
