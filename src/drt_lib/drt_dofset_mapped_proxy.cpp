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
#include "../drt_lib/drt_matchingoctree.H"
#include "../drt_lib/drt_condition_utils.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetMappedProxy::DofSetMappedProxy( Teuchos::RCP<DofSet>  dofset,
                                          const Teuchos::RCP<const DRT::Discretization> sourcedis,
                                          const std::string& couplingcond,
                                          const std::set<int> condids)
  : DofSetProxy(&(*dofset)),
    targetlidtosourcegidmapping_(Teuchos::null),
    dofset_(dofset),
    sourcedis_(sourcedis),
    couplingcond_(couplingcond),
    condids_(condids)
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

  // get the respective nodes which are in the condition
  std::map<int, Teuchos::RCP<std::vector<int> > > nodes;
  DRT::UTILS::FindConditionedNodes(dis, conds, nodes);
  std::map<int, Teuchos::RCP<std::vector<int> > > nodes_source;
  DRT::UTILS::FindConditionedNodes(*sourcedis_, conds_source, nodes_source);

  // map that will be filled with coupled nodes
  // mapping: target node gid to (source node gid, distance)
  std::map<int,std::pair<int,double> > coupling;

  // define iterators
  std::map<int, Teuchos::RCP<std::vector<int> > >::iterator iter_target;
  std::map<int, Teuchos::RCP<std::vector<int> > >::iterator iter_source;

  for (std::set<int>::iterator it=condids_.begin(); it!=condids_.end(); ++it)
  {
    // find corresponding condition on source discretization
    iter_target = nodes.find(*it);
    if(iter_target == nodes.end())
      dserror("Condition ID %i not found in Coupling condition on source discretization s!",*it,sourcedis_->Name().c_str());


    // find corresponding condition on source discretization
    iter_source = nodes_source.find(*it);
    if(iter_source == nodes_source.end())
      dserror("Condition ID %i not found in Coupling condition on source discretization %s!",*it,sourcedis_->Name().c_str());


    // get the nodes
    Teuchos::RCP<std::vector<int> > sourcenodes = iter_source->second;

    // initialize search tree for search
    DRT::UTILS::NodeMatchingOctree tree(dis, *iter_target->second);

    // map that will be filled with coupled nodes for this condition
    // mapping: target node gid to (source node gid, distance)
    // note: FindMatch loops over all SOURCE (i.e. slave) nodes
    //       and finds corresponding target nodes.
    std::map<int,std::pair<int,double> > condcoupling;
    // match target and source nodes using octtree
    tree.FindMatch(*sourcedis_, *sourcenodes, condcoupling);

    // check if all nodes where matched for this condition ID
    if (iter_target->second->size() != condcoupling.size())
      dserror("Did not get unique target to source spatial node coordinate mapping.\n"
          "targetnodes.size()=%d, coupling.size()=%d.\n"
          "The heterogeneous reaction strategy requires matching source and target meshes!",
          iter_target->second->size(), condcoupling.size());

    // insert found coupling of this condition ID into map of all match nodes
    coupling.insert(condcoupling.begin(), condcoupling.end());

  }// loop over all condition ids

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
      Teuchos::rcp(new Epetra_IntVector(Copy, *targetnodemap, permsourcenodemap->MyGlobalElements()));

  targetlidtosourcegidmapping_ =
    Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));

  targetlidtosourcegidmapping_->PutValue(-1);

  // export to column map
  LINALG::Export(*permsourcenodevec,*targetlidtosourcegidmapping_);


  return start;
}
