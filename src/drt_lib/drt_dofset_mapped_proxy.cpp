/*----------------------------------------------------------------------*/
/*!
 \file drt_dofset_mapped_proxy.cpp

 \brief  A proxy of a dofset that does not rely on same GID/LID numbers but uses
         a defined node mapping instead (not implemented for element DOFs).

   \level 3

   \maintainer Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/


#include "drt_dofset_mapped_proxy.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_condition_utils.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetMappedProxy::DofSetMappedProxy( Teuchos::RCP<DofSet>  dofset,
                                          const Teuchos::RCP<const DRT::Discretization> sourcedis,
                                          const std::string& couplingcond)
  : DofSetProxy(&(*dofset)),
    sourcetotargetnodemapping_(Teuchos::null),
    dofset_(dofset),
    sourcedis_(sourcedis),
    couplingcond_(couplingcond)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetMappedProxy::AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start)
{
  // make sure the real dofset gets the dofmaps
  DofSetProxy::AssignDegreesOfFreedom(dis,dspos,start);

  //get condition which defines the coupling on target discretization
  std::vector<DRT::Condition*> conds;
  dis.GetCondition(couplingcond_, conds);

  //get condition which defines the coupling on source discretization
  std::vector<DRT::Condition*> conds_source;
  sourcedis_->GetCondition(couplingcond_, conds_source);

  // at least one condition needs to be defined on each discretization
  if(conds.size() == 0 or conds_source.size() == 0)
    dserror("No coupling condition defined on one or both discretizations!");
  if(conds.size() != conds_source.size())
    dserror("Unequal number of coupling conditions '%s' on discretizations!",couplingcond_.c_str());

  // get the respective nodes which are in the condition
  std::map<int, Teuchos::RCP<std::vector<int> > > nodes;
  DRT::UTILS::FindConditionedNodes(dis, conds, nodes);
  std::map<int, Teuchos::RCP<std::vector<int> > > nodes_source;
  DRT::UTILS::FindConditionedNodes(*sourcedis_, conds_source, nodes_source);

  // map that will be filled with coupled nodes
  // mapping: target node gid to (source node gid, distance)
  std::map<int,std::pair<int,double> > coupling;

  // define iterators
  std::map<int, Teuchos::RCP<std::vector<int> > >::iterator iter;
  std::map<int, Teuchos::RCP<std::vector<int> > >::iterator iter_source;

  // loop over nodal clouds of target condition
  for (iter = nodes.begin(); iter != nodes.end(); ++iter)
  {
    //get condition ID
    const int condID = iter->first;

    // find corresponding condition on source discretization
    iter_source = nodes_source.find(condID);
    if(iter_source == nodes_source.end())
      dserror("Condition ID %i not found in Coupling condition on source discretization!");

    // get the nodes
    Teuchos::RCP<std::vector<int> > sourcenodes = iter_source->second;

    // initialize search tree for search
    DRT::UTILS::NodeMatchingOctree tree(dis, *iter->second);

    // map that will be filled with coupled nodes for this condition
    // mapping: target node gid to (source node gid, distance)
    std::map<int,std::pair<int,double> > condcoupling;
    // match target and source nodes using octtree
    tree.FindMatch(*sourcedis_, *sourcenodes, condcoupling);

    // check if all nodes where matched for this condition ID
    if (sourcenodes->size() != condcoupling.size())
      dserror("Did not get 1:1 correspondence. \nsourcenodes.size()=%d, coupling.size()=%d."
          "The heterogeneous reaction strategy requires matching source and target meshes!",
          sourcenodes->size(), condcoupling.size());

    // insert found coupling of this condition ID into map of all match nodes
    coupling.insert(condcoupling.begin(), condcoupling.end());
  }

  // check if all nodes where matched for coupling
  if(dis.NodeRowMap()->NumMyElements() != (int) coupling.size())
    dserror("All nodes of the discretization need to be coupled. Number of nodes in discretization: %i, number of coupled nodes: %i",
        dis.NodeRowMap()->NumMyElements(), coupling.size());

  // clone communicator of target discretization
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp( dis.Comm().Clone());

  // extract permutation
  std::vector<int> targetnodes(dis.NodeRowMap()->MyGlobalElements(),
      dis.NodeRowMap()->MyGlobalElements() + dis.NodeRowMap()->NumMyElements());

  std::vector<int> patchedtargetnodes;
  patchedtargetnodes.reserve(coupling.size());
  std::vector<int> permsourcenodes;
  permsourcenodes.reserve(coupling.size());

  for (unsigned i=0; i<targetnodes.size(); ++i)
  {
    const int gid = targetnodes[i];

    // We allow to hand in target nodes that do not take part in the
    // coupling. If this is undesired behaviour the user has to make
    // sure all nodes were used.
    if (coupling.find(gid) != coupling.end())
    {
      std::pair<int,double>& coupled = coupling[gid];
      patchedtargetnodes.push_back(gid);
      permsourcenodes.push_back(coupled.first);
    }
  }

  // Epetra maps
  Teuchos::RCP<Epetra_Map> targetnodemap =
    Teuchos::rcp(new Epetra_Map(-1, patchedtargetnodes.size(), &patchedtargetnodes[0], 0, *com));

  Teuchos::RCP<Epetra_Map> permsourcenodemap =
    Teuchos::rcp(new Epetra_Map(-1, permsourcenodes.size(), &permsourcenodes[0], 0, *com));

  // we expect to get maps of exactly the same shape
  if (not targetnodemap->PointSameAs(*permsourcenodemap))
    dserror("target and permuted source node maps do not match");

  // export target nodes to source node distribution

  Teuchos::RCP<Epetra_IntVector> permsourcenodevec =
      Teuchos::rcp(new Epetra_IntVector(Copy, *dis.NodeRowMap(), permsourcenodemap->MyGlobalElements()));

  sourcetotargetnodemapping_ =
    Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));

  // export to column map
  LINALG::Export(*permsourcenodevec,*sourcetotargetnodemapping_);

  return start;
}
