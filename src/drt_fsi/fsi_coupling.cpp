
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fsi_coupling.H"
#include "../drt_lib/drt_nodematchingoctree.H"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


using namespace std;


FSI::Coupling::Coupling()
{
}


void FSI::Coupling::SetupConditionCoupling(const DRT::Discretization& masterdis,
                                           const DRT::Discretization& slavedis,
                                           std::string coupname)
{
  set<int> masternodeset;
  set<int> slavenodeset;

  FindCondNodes(masterdis, coupname, masternodeset);
  FindCondNodes(slavedis,  coupname, slavenodeset);

  vector<int> masternodes(masternodeset.begin(), masternodeset.end());
  vector<int> slavenodes(slavenodeset.begin(), slavenodeset.end());

  SetupCoupling(masterdis, slavedis, masternodes, slavenodes);

  // test for completeness
  if (static_cast<int>(masternodes.size())*genprob.ndim != masterdofmap_->NumMyElements())
    dserror("failed to setup master nodes properly");
  if (static_cast<int>(slavenodes.size())*genprob.ndim != slavedofmap_->NumMyElements())
    dserror("failed to setup slave nodes properly");
}


void FSI::Coupling::SetupCoupling(const DRT::Discretization& masterdis,
                                  const DRT::Discretization& slavedis,
                                  const std::vector<int>& masternodes,
                                  const std::vector<int>& slavenodes)
{
  std::vector<int> patchedmasternodes(masternodes);
  std::vector<int> permslavenodes;
  MatchNodes(masterdis, slavedis, patchedmasternodes, permslavenodes, slavenodes);

  // Epetra maps in original distribution

  Teuchos::RefCountPtr<Epetra_Map> masternodemap =
    rcp(new Epetra_Map(-1, patchedmasternodes.size(), &patchedmasternodes[0], 0, masterdis.Comm()));

  Teuchos::RefCountPtr<Epetra_Map> slavenodemap =
    rcp(new Epetra_Map(-1, slavenodes.size(), &slavenodes[0], 0, slavedis.Comm()));

  Teuchos::RefCountPtr<Epetra_Map> permslavenodemap =
    rcp(new Epetra_Map(-1, permslavenodes.size(), &permslavenodes[0], 0, slavedis.Comm()));

  FinishCoupling(masterdis, slavedis, masternodemap, slavenodemap, permslavenodemap);
}


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

  Teuchos::RefCountPtr<Epetra_Map> masternodemap =
    rcp(new Epetra_Map(-1, mastervect.size(), &mastervect[0], 0, masterdis.Comm()));

  Teuchos::RefCountPtr<Epetra_Map> slavenodemap =
    rcp(new Epetra_Map(slavenodes));

  Teuchos::RefCountPtr<Epetra_Map> permslavenodemap =
    rcp(new Epetra_Map(-1, permslavenodes.size(), &permslavenodes[0], 0, slavedis.Comm()));

  FinishCoupling(masterdis, slavedis, masternodemap, slavenodemap, permslavenodemap);
}


void FSI::Coupling::MatchNodes(const DRT::Discretization& masterdis,
                               const DRT::Discretization& slavedis,
                               std::vector<int>& masternodes,
                               std::vector<int>& permslavenodes,
                               const std::vector<int>& slavenodes)
{
  // match master and slave nodes using Peter's octtree

  NodeMatchingOctree tree(masterdis, masternodes);

  map<int,pair<int,double> > coupling;
  tree.FindMatch(slavedis, slavenodes, coupling);

  if (masternodes.size() != coupling.size())
    dserror("did not get 1:1 correspondence");

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
      if (coupled.second > 1e-8)
        dserror("coupled nodes (%d,%d) do not match", gid, coupled.first);
      patchedmasternodes.push_back(gid);
      permslavenodes.push_back(coupled.first);
    }
  }

  // return new list of master nodes via reference
  swap(masternodes,patchedmasternodes);
}



void FSI::Coupling::FinishCoupling(const DRT::Discretization& masterdis,
                                   const DRT::Discretization& slavedis,
                                   Teuchos::RefCountPtr<Epetra_Map> masternodemap,
                                   Teuchos::RefCountPtr<Epetra_Map> slavenodemap,
                                   Teuchos::RefCountPtr<Epetra_Map> permslavenodemap)
{
  // we expect to get maps of exactly the same shape
  if (not masternodemap->PointSameAs(*permslavenodemap))
    dserror("master and permutated slave node maps do not match");

  // export master nodes to slave node distribution

  // To do so we create vectors that contain the values of the master
  // maps, assigned to the slave maps. On the master side we actually
  // create just a view on the map! This vector must not be changed!
  Teuchos::RefCountPtr<Epetra_IntVector> masternodevec =
    rcp(new Epetra_IntVector(View, *permslavenodemap, masternodemap->MyGlobalElements()));

  Teuchos::RefCountPtr<Epetra_IntVector> permmasternodevec =
    rcp(new Epetra_IntVector(*slavenodemap));

  Epetra_Export masternodeexport(*permslavenodemap, *slavenodemap);
  int err = permmasternodevec->Export(*masternodevec, masternodeexport, Insert);
  if (err)
    dserror("failed to export master nodes");

  Teuchos::RefCountPtr<Epetra_Map> permmasternodemap =
    rcp(new Epetra_Map(-1, permmasternodevec->MyLength(), permmasternodevec->Values(), 0, masterdis.Comm()));

  if (not slavenodemap->PointSameAs(*permmasternodemap))
    dserror("slave and permutated master node maps do not match");

  masternodevec = Teuchos::null;
  permmasternodevec = Teuchos::null;

  BuildDofMaps(masterdis, masternodemap, permmasternodemap, masterdofmap_, permmasterdofmap_, masterexport_);
  BuildDofMaps(slavedis,  slavenodemap,  permslavenodemap,  slavedofmap_,  permslavedofmap_,  slaveexport_);
}


void FSI::Coupling::BuildDofMaps(const DRT::Discretization& dis,
                                 RefCountPtr<Epetra_Map> nodemap,
                                 RefCountPtr<Epetra_Map> permnodemap,
                                 RefCountPtr<Epetra_Map>& dofmap,
                                 RefCountPtr<Epetra_Map>& permdofmap,
                                 RefCountPtr<Epetra_Export>& exporter)
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


Teuchos::RefCountPtr<Epetra_Vector> FSI::Coupling::MasterToSlave(Teuchos::RefCountPtr<Epetra_Vector> mv)
{
#ifdef DEBUG
  if (not mv->Map().SameAs(*masterdofmap_))
    dserror("master dof map vector expected");
#endif

  Epetra_Vector perm(*permslavedofmap_);
  copy(mv->Values(), mv->Values()+mv->MyLength(), perm.Values());

  Teuchos::RefCountPtr<Epetra_Vector> sv =
    rcp(new Epetra_Vector(*slavedofmap_));

  int err = sv->Export(perm,*slaveexport_,Insert);
  if (err)
    dserror("Export to slave distribution returned err=%d",err);

  return sv;
}


Teuchos::RefCountPtr<Epetra_Vector> FSI::Coupling::SlaveToMaster(Teuchos::RefCountPtr<Epetra_Vector> sv)
{
#ifdef DEBUG
  if (not sv->Map().SameAs(*slavedofmap_))
    dserror("slave dof map vector expected");
#endif

  Epetra_Vector perm(*permmasterdofmap_);
  copy(sv->Values(), sv->Values()+sv->MyLength(), perm.Values());

  Teuchos::RefCountPtr<Epetra_Vector> mv =
    rcp(new Epetra_Vector(*masterdofmap_));

  int err = mv->Export(perm,*masterexport_,Insert);
  if (err)
    dserror("Export to master distribution returned err=%d",err);

  return mv;
}


void FSI::Coupling::FindCondNodes(const DRT::Discretization& dis, std::string condname, std::set<int>& nodes)
{
  int myrank = dis.Comm().MyPID();
  vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);
  for (unsigned i=0; i<conds.size(); ++i)
  {
    const vector<int>* n = conds[i]->Get<vector<int> >("Node Ids");
    for (unsigned j=0; j<n->size(); ++j)
    {
      int gid = (*n)[j];
      if (dis.HaveGlobalNode(gid) and dis.gNode(gid)->Owner()==myrank)
      {
        nodes.insert(gid);
      }
    }
    //std::copy(n->begin(), n->end(), inserter(nodes, nodes.begin()));
  }
}


#endif
#endif
