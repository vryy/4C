
#ifdef CCADISCRET

#include "fsi_conditiondofmap.H"
#include "fsi_utils.H"

#include "../drt_lib/drt_discret.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C"
{
#include "../headers/standardtypes.h"
}


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::ConditionDofMap::ConditionDofMap(Teuchos::RCP<DRT::Discretization> dis)
  : dis_(dis)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConditionDofMap::SetupCondDofMap(std::string condname)
{
  std::set<int> nodeset;
  FindCondNodes(*dis_, condname, nodeset);

  nodes_.assign(nodeset.begin(), nodeset.end());

  std::vector<int> conddofmapvec;
  conddofmapvec.reserve(genprob.ndim*nodes_.size());

  std::vector<int> noneconddofmapvec;

  for (unsigned i=0; i<nodes_.size(); ++i)
  {
    DRT::Node* actnode = dis_->gNode(nodes_[i]);
    std::vector<int> dof = dis_->Dof(actnode);

    //copy(&dof[0], &dof[0]+genprob.ndim, back_inserter(dofs_[nodes_[i]]));
    copy(&dof[0], &dof[0]+genprob.ndim, back_inserter(conddofmapvec));
  }

  // dof map is the original, unpermuted distribution of dofs
  SetupCondDofMap(Teuchos::rcp(new Epetra_Map(-1, conddofmapvec.size(), &conddofmapvec[0], 0, dis_->Comm())));

  conddofmapvec.clear();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConditionDofMap::SetupCondDofMap(Teuchos::RCP<Epetra_Map> condmap)
{
  conddofmap_ = condmap;
  condimporter_ = Teuchos::rcp(new Epetra_Import(*conddofmap_, *dis_->DofRowMap()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConditionDofMap::SetupOtherDofMap()
{
  const Epetra_Map* dofmap = dis_->DofRowMap();

  int* gids = dofmap->MyGlobalElements();
  int elements = dofmap->NumMyElements();

  std::vector<int> otherdofmapvec;
  otherdofmapvec.reserve(elements);

  std::remove_copy_if(gids,
                      gids+elements,
                      back_inserter(otherdofmapvec),
                      FSI::Utils::MyGID(&*conddofmap_));

  otherdofmap_ = Teuchos::rcp(new Epetra_Map(-1,otherdofmapvec.size(),&otherdofmapvec[0],0,dofmap->Comm()));

  otherdofmapvec.clear();

  otherimporter_ = Teuchos::rcp(new Epetra_Import(*otherdofmap_, *dis_->DofRowMap()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::ConditionDofMap::ExtractCondVector(Teuchos::RCP<const Epetra_Vector> full) const
{
  Teuchos::RefCountPtr<Epetra_Vector> cond = Teuchos::rcp(new Epetra_Vector(*conddofmap_));
  ExtractCondVector(full,cond);
  return cond;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::ConditionDofMap::ExtractOtherVector(Teuchos::RCP<const Epetra_Vector> full) const
{
  Teuchos::RefCountPtr<Epetra_Vector> other = Teuchos::rcp(new Epetra_Vector(*otherdofmap_));
  ExtractOtherVector(full,other);
  return other;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConditionDofMap::ExtractCondVector(Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> cond) const
{
  int err = cond->Import(*full,*condimporter_,Insert);
  if (err)
    dserror("Import using importer returned err=%d",err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConditionDofMap::ExtractOtherVector(Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> other) const
{
  int err = other->Import(*full,*otherimporter_,Insert);
  if (err)
    dserror("Import using importer returned err=%d",err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::ConditionDofMap::InsertCondVector(Teuchos::RCP<const Epetra_Vector> cond) const
{
  Teuchos::RCP<Epetra_Vector> full = Teuchos::rcp(new Epetra_Vector(*dis_->DofRowMap()));
  InsertCondVector(cond,full);
  return full;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::ConditionDofMap::InsertOtherVector(Teuchos::RCP<const Epetra_Vector> other) const
{
  Teuchos::RCP<Epetra_Vector> full = Teuchos::rcp(new Epetra_Vector(*dis_->DofRowMap()));
  InsertOtherVector(other,full);
  return full;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConditionDofMap::InsertCondVector(Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const
{
  int err = full->Export(*cond,*condimporter_,Insert);
  if (err)
    dserror("Export using importer returned err=%d",err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConditionDofMap::InsertOtherVector(Teuchos::RCP<const Epetra_Vector> other, Teuchos::RCP<Epetra_Vector> full) const
{
  int err = full->Export(*other,*otherimporter_,Insert);
  if (err)
    dserror("Export using importer returned err=%d",err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConditionDofMap::FindCondNodes(const DRT::Discretization& dis,
                                         std::string condname,
                                         std::set<int>& nodes)
{
  int myrank = dis.Comm().MyPID();
  std::vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);
  for (unsigned i=0; i<conds.size(); ++i)
  {
    const vector<int>* n = conds[i]->Nodes();
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
