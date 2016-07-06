/*----------------------------------------------------------------------*/
/*!
 \file drt_dofset_merged_proxy.cpp

 \brief A proxy of a dofset that adds additional, existing degrees of freedom from the same
        discretization to nodes (not implemented for element DOFs).

   \level 3

   \maintainer Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/


#include "drt_dofset_merged_proxy.H"

#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../linalg/linalg_utils.H"

#include <Epetra_Export.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetMergedProxy::DofSetMergedProxy( Teuchos::RCP<DofSet>  dofset,
                                          const Teuchos::RCP<const DRT::Discretization> sourcedis,
                                          const std::string& couplingcond_master,
                                          const std::string& couplingcond_slave)
  : DofSetProxy(&(*dofset)),
    slavetomasternodemapping_(Teuchos::null),
    dofset_(dofset),
    sourcedis_(sourcedis),
    couplingcond_master_(couplingcond_master),
    couplingcond_slave_(couplingcond_slave)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetMergedProxy::AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start)
{
  // make sure the real dofset gets the dofmaps
  DofSetProxy::AssignDegreesOfFreedom(dis,dspos,start);

  // copy communicator
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp( sourcedis_->Comm().Clone());

  // get nodes to be couplid
  std::vector<int> masternodes;
  DRT::UTILS::FindConditionedNodes(*sourcedis_,couplingcond_master_,masternodes);
  std::vector<int> slavenodes;
  DRT::UTILS::FindConditionedNodes(*sourcedis_,couplingcond_slave_,slavenodes);

  // initialize search tree
  DRT::UTILS::NodeMatchingOctree tree(*sourcedis_, masternodes);

  // match master and slave nodes using octtree
  std::map<int,std::pair<int,double> > coupling;
  tree.FindMatch(*sourcedis_, slavenodes, coupling);

  // all nodes should be coupled
  if (masternodes.size() != coupling.size())
    dserror("Did not get 1:1 correspondence. \nmasternodes.size()=%d, coupling.size()=%d."
        "DofSetMergedProxy requires matching slave and master meshes!",
            masternodes.size(), coupling.size());

    // extract permutation

  std::vector<int> patchedmasternodes;
  patchedmasternodes.reserve(coupling.size());
  std::vector<int> permslavenodes;
  permslavenodes.reserve(slavenodes.size());

  for (unsigned i=0; i<masternodes.size(); ++i)
  {
    const int gid = masternodes[i];

    // We allow to hand in master nodes that do not take part in the
    // coupling. If this is undesired behaviour the user has to make
    // sure all nodes were used.
    if (coupling.find(gid) != coupling.end())
    {
      std::pair<int,double>& coupled = coupling[gid];
      patchedmasternodes.push_back(gid);
      permslavenodes.push_back(coupled.first);
    }
  }

  // Epetra maps in original distribution

  Teuchos::RCP<Epetra_Map> masternodemap =
    Teuchos::rcp(new Epetra_Map(-1, patchedmasternodes.size(), &patchedmasternodes[0], 0, *com));

  Teuchos::RCP<Epetra_Map> slavenodemap =
    Teuchos::rcp(new Epetra_Map(-1, slavenodes.size(), &slavenodes[0], 0, *com));

  Teuchos::RCP<Epetra_Map> permslavenodemap =
    Teuchos::rcp(new Epetra_Map(-1, permslavenodes.size(), &permslavenodes[0], 0, *com));

  // we expect to get maps of exactly the same shape
  if (not masternodemap->PointSameAs(*permslavenodemap))
    dserror("master and permuted slave node maps do not match");

  // export master nodes to slave node distribution

  // To do so we create vectors that contain the values of the master
  // maps, assigned to the slave maps. On the master side we actually
  // create just a view on the map! This vector must not be changed!
  Teuchos::RCP<Epetra_IntVector> masternodevec =
    Teuchos::rcp(new Epetra_IntVector(View, *permslavenodemap, masternodemap->MyGlobalElements()));

  Teuchos::RCP<Epetra_IntVector> permmasternodevec =
    Teuchos::rcp(new Epetra_IntVector(*slavenodemap));

  // build exporter
  Epetra_Export masternodeexport(*permslavenodemap, *slavenodemap);
  const int err = permmasternodevec->Export(*masternodevec, masternodeexport, Insert);
  if (err)
    dserror("failed to export master nodes");

  // initialize final mapping
  slavetomasternodemapping_ =
    Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));

  // export to column map
  LINALG::Export(*permmasternodevec,*slavetomasternodemapping_);

  return start;
}
