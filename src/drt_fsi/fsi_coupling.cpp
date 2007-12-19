
#ifdef CCADISCRET

#include "fsi_coupling.H"
#include "../drt_lib/drt_nodematchingoctree.H"

#include "fsi_conditiondofmap.H"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


using namespace std;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Coupling::Coupling()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Coupling::SetupConditionCoupling(const ConditionDofMap& master,
                                           const ConditionDofMap& slave)
{
  SetupCoupling(master.Discret(), slave.Discret(), master.Nodes(), slave.Nodes());

  // test for completeness
  if (static_cast<int>(master.Nodes().size())*genprob.ndim != masterdofmap_->NumMyElements())
    dserror("failed to setup master nodes properly");
  if (static_cast<int>(slave.Nodes().size())*genprob.ndim != slavedofmap_->NumMyElements())
    dserror("failed to setup slave nodes properly");

  // Now swap in the maps we already had.
  // So we did a little more work than required. But there are cases
  // where we have to do that work (fluid-ale coupling) and we want to
  // use just one setup implementation.
  //
  // The point is to make sure there is only one map for each
  // interface.

  if (not masterdofmap_->SameAs(*master.CondDofMap()))
    dserror("master dof map mismatch");

  if (not slavedofmap_->SameAs(*slave.CondDofMap()))
    dserror("master dof map mismatch");

  masterdofmap_ = master.CondDofMap();
  masterexport_ = rcp(new Epetra_Export(*permmasterdofmap_, *masterdofmap_));

  slavedofmap_ = slave.CondDofMap();
  slaveexport_ = rcp(new Epetra_Export(*permslavedofmap_, *slavedofmap_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Coupling::SetupCoupling(const DRT::Discretization& masterdis,
                                  const DRT::Discretization& slavedis,
                                  const std::vector<int>& masternodes,
                                  const std::vector<int>& slavenodes)
{
  std::vector<int> patchedmasternodes(masternodes);
  std::vector<int> permslavenodes;
  MatchNodes(masterdis, slavedis, patchedmasternodes, permslavenodes, slavenodes);

  // Epetra maps in original distribution

  Teuchos::RCP<Epetra_Map> masternodemap =
    rcp(new Epetra_Map(-1, patchedmasternodes.size(), &patchedmasternodes[0], 0, masterdis.Comm()));

  Teuchos::RCP<Epetra_Map> slavenodemap =
    rcp(new Epetra_Map(-1, slavenodes.size(), &slavenodes[0], 0, slavedis.Comm()));

  Teuchos::RCP<Epetra_Map> permslavenodemap =
    rcp(new Epetra_Map(-1, permslavenodes.size(), &permslavenodes[0], 0, slavedis.Comm()));

  FinishCoupling(masterdis, slavedis, masternodemap, slavenodemap, permslavenodemap);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Coupling::SetupCoupling(const DRT::Discretization& masterdis,
                                  const DRT::Discretization& slavedis,
                                  const Epetra_Map& masternodes,
                                  const Epetra_Map& slavenodes)
{
  vector<int> mastervect(masternodes.MyGlobalElements(),
                         masternodes.MyGlobalElements() + masternodes.NumMyElements());
  vector<int> slavevect(slavenodes.MyGlobalElements(),
                        slavenodes.MyGlobalElements() + slavenodes.NumMyElements());
  vector<int> permslavenodes;

  MatchNodes(masterdis, slavedis, mastervect, permslavenodes, slavevect);

  // Epetra maps in original distribution

  Teuchos::RCP<Epetra_Map> masternodemap =
    rcp(new Epetra_Map(-1, mastervect.size(), &mastervect[0], 0, masterdis.Comm()));

  Teuchos::RCP<Epetra_Map> slavenodemap =
    rcp(new Epetra_Map(slavenodes));

  Teuchos::RCP<Epetra_Map> permslavenodemap =
    rcp(new Epetra_Map(-1, permslavenodes.size(), &permslavenodes[0], 0, slavedis.Comm()));

  FinishCoupling(masterdis, slavedis, masternodemap, slavenodemap, permslavenodemap);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Coupling::MatchNodes(const DRT::Discretization& masterdis,
                               const DRT::Discretization& slavedis,
                               std::vector<int>& masternodes,
                               std::vector<int>& permslavenodes,
                               const std::vector<int>& slavenodes)
{
  // match master and slave nodes using Peter's octtree

  DRT::Utils::NodeMatchingOctree tree(masterdis, masternodes);

  map<int,pair<int,double> > coupling;
  tree.FindMatch(slavedis, slavenodes, coupling);

  if (masternodes.size() != coupling.size())
    dserror("Did not get 1:1 correspondence. masternodes.size()=%d, coupling.size()=%d",
            masternodes.size(), coupling.size());

  // extract permutation

  vector<int> patchedmasternodes;
  patchedmasternodes.reserve(coupling.size());
  permslavenodes.reserve(slavenodes.size());

  for (unsigned i=0; i<masternodes.size(); ++i)
  {
    int gid = masternodes[i];

    // We allow to hand in master nodes that do not take part in the
    // coupling. If this is undesired behaviour the user has to make
    // sure all nodes were used.
    if (coupling.find(gid) != coupling.end())
    {
      pair<int,double>& coupled = coupling[gid];
      if (coupled.second > 1e-7)
        dserror("Coupled nodes (%d,%d) do not match. difference=%e", gid, coupled.first, coupled.second);
      patchedmasternodes.push_back(gid);
      permslavenodes.push_back(coupled.first);
    }
  }

  // return new list of master nodes via reference
  swap(masternodes,patchedmasternodes);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Coupling::FinishCoupling(const DRT::Discretization& masterdis,
                                   const DRT::Discretization& slavedis,
                                   Teuchos::RCP<Epetra_Map> masternodemap,
                                   Teuchos::RCP<Epetra_Map> slavenodemap,
                                   Teuchos::RCP<Epetra_Map> permslavenodemap)
{
  // we expect to get maps of exactly the same shape
  if (not masternodemap->PointSameAs(*permslavenodemap))
    dserror("master and permutated slave node maps do not match");

  // export master nodes to slave node distribution

  // To do so we create vectors that contain the values of the master
  // maps, assigned to the slave maps. On the master side we actually
  // create just a view on the map! This vector must not be changed!
  Teuchos::RCP<Epetra_IntVector> masternodevec =
    rcp(new Epetra_IntVector(View, *permslavenodemap, masternodemap->MyGlobalElements()));

  Teuchos::RCP<Epetra_IntVector> permmasternodevec =
    rcp(new Epetra_IntVector(*slavenodemap));

  Epetra_Export masternodeexport(*permslavenodemap, *slavenodemap);
  int err = permmasternodevec->Export(*masternodevec, masternodeexport, Insert);
  if (err)
    dserror("failed to export master nodes");

  Teuchos::RCP<Epetra_Map> permmasternodemap =
    rcp(new Epetra_Map(-1, permmasternodevec->MyLength(), permmasternodevec->Values(), 0, masterdis.Comm()));

  if (not slavenodemap->PointSameAs(*permmasternodemap))
    dserror("slave and permutated master node maps do not match");

  masternodevec = Teuchos::null;
  permmasternodevec = Teuchos::null;

  BuildDofMaps(masterdis, masternodemap, permmasternodemap, masterdofmap_, permmasterdofmap_, masterexport_);
  BuildDofMaps(slavedis,  slavenodemap,  permslavenodemap,  slavedofmap_,  permslavedofmap_,  slaveexport_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Coupling::BuildDofMaps(const DRT::Discretization& dis,
                                 RCP<Epetra_Map> nodemap,
                                 RCP<Epetra_Map> permnodemap,
                                 RCP<Epetra_Map>& dofmap,
                                 RCP<Epetra_Map>& permdofmap,
                                 RCP<Epetra_Export>& exporter)
{
  // communicate dofs

  vector<int> dofmapvec;
  map<int, vector<int> > dofs;

  int* nodes = nodemap->MyGlobalElements();
  int numnode = nodemap->NumMyElements();

  for (int i=0; i<numnode; ++i)
  {
    DRT::Node* actnode = dis.gNode(nodes[i]);
    vector<int> dof = dis.Dof(actnode);
    copy(&dof[0], &dof[0]+genprob.ndim, back_inserter(dofs[nodes[i]]));
    copy(&dof[0], &dof[0]+genprob.ndim, back_inserter(dofmapvec));
  }

  // dof map is the original, unpermuted distribution of dofs
  dofmap = rcp(new Epetra_Map(-1, dofmapvec.size(), &dofmapvec[0], 0, dis.Comm()));

  dofmapvec.clear();

  DRT::Exporter exportdofs(*nodemap,*permnodemap,dis.Comm());
  exportdofs.Export(dofs);

  nodes = permnodemap->MyGlobalElements();
  numnode = permnodemap->NumMyElements();

  for (int i=0; i<numnode; ++i)
  {
    vector<int>& dof = dofs[nodes[i]];
    copy(dof.begin(), dof.end(), back_inserter(dofmapvec));
  }

  dofs.clear();

  // permutated dof map according to a given permutated node map
  permdofmap = rcp(new Epetra_Map(-1, dofmapvec.size(), &dofmapvec[0], 0, dis.Comm()));

  // prepare communication plan to create a dofmap out of a permutated
  // dof map
  exporter = rcp(new Epetra_Export(*permdofmap, *dofmap));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Coupling::MasterToSlave(Teuchos::RCP<const Epetra_Vector> mv) const
{
  Teuchos::RCP<Epetra_Vector> sv =
    Teuchos::rcp(new Epetra_Vector(*slavedofmap_));

  MasterToSlave(mv,sv);

  return sv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Coupling::SlaveToMaster(Teuchos::RCP<const Epetra_Vector> sv) const
{
  Teuchos::RCP<Epetra_Vector> mv =
    Teuchos::rcp(new Epetra_Vector(*masterdofmap_));

  SlaveToMaster(sv,mv);

  return mv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> FSI::Coupling::MasterToSlave(Teuchos::RCP<const Epetra_MultiVector> mv) const
{
  Teuchos::RCP<Epetra_MultiVector> sv =
    Teuchos::rcp(new Epetra_MultiVector(*slavedofmap_,mv->NumVectors()));

  MasterToSlave(mv,sv);

  return sv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> FSI::Coupling::SlaveToMaster(Teuchos::RCP<const Epetra_MultiVector> sv) const
{
  Teuchos::RCP<Epetra_MultiVector> mv =
    Teuchos::rcp(new Epetra_MultiVector(*masterdofmap_,sv->NumVectors()));

  SlaveToMaster(sv,mv);

  return mv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Coupling::MasterToSlave(Teuchos::RCP<const Epetra_MultiVector> mv, Teuchos::RCP<Epetra_MultiVector> sv) const
{
#ifdef DEBUG
  if (not mv->Map().SameAs(*masterdofmap_))
    dserror("master dof map vector expected");
  if (not sv->Map().SameAs(*slavedofmap_))
    dserror("slave dof map vector expected");
  if (sv->NumVectors()!=mv->NumVectors())
    dserror("column number mismatch %d!=%d",sv->NumVectors(),mv->NumVectors());
#endif

  Epetra_MultiVector perm(*permslavedofmap_,mv->NumVectors());
  copy(mv->Values(), mv->Values()+(mv->MyLength()*mv->NumVectors()), perm.Values());

  int err = sv->Export(perm,*slaveexport_,Insert);
  if (err)
    dserror("Export to slave distribution returned err=%d",err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Coupling::SlaveToMaster(Teuchos::RCP<const Epetra_MultiVector> sv, Teuchos::RCP<Epetra_MultiVector> mv) const
{
#ifdef DEBUG
  if (not mv->Map().SameAs(*masterdofmap_))
    dserror("master dof map vector expected");
  if (not sv->Map().SameAs(*slavedofmap_))
    dserror("slave dof map vector expected");
  if (sv->NumVectors()!=mv->NumVectors())
    dserror("column number mismatch %d!=%d",sv->NumVectors(),mv->NumVectors());
#endif

  Epetra_MultiVector perm(*permmasterdofmap_,sv->NumVectors());
  copy(sv->Values(), sv->Values()+(sv->MyLength()*sv->NumVectors()), perm.Values());

  int err = mv->Export(perm,*masterexport_,Insert);
  if (err)
    dserror("Export to master distribution returned err=%d",err);
}


#endif
