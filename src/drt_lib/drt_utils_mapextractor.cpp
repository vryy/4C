
#ifdef CCADISCRET

#include "drt_utils_mapextractor.H"
#include "linalg_utils.H"

#include <Teuchos_getConst.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C"
{
#include "../headers/standardtypes.h"
}


DRT::UTILS::MapExtractor::MapExtractor()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MapExtractor::SetupMaps(const Teuchos::RCP<DRT::Discretization> dis,
                                         const std::set<int>& conddofset, const std::set<int>& otherdofset)
{
  dofrowmap_ = Teuchos::rcp(dis->DofRowMap(),false);

  std::vector<int> conddofmapvec;
  conddofmapvec.reserve(conddofset.size());
  conddofmapvec.assign(conddofset.begin(), conddofset.end());
  conddofmap_ = Teuchos::rcp(new Epetra_Map(-1,
                                            conddofmapvec.size(),
                                            &conddofmapvec[0],
                                            0,
                                            dis->Comm()));
  condimporter_ = Teuchos::rcp(new Epetra_Import(*conddofmap_, *dofrowmap_));
  conddofmapvec.clear();

  std::vector<int> otherdofmapvec;
  otherdofmapvec.reserve(otherdofset.size());
  otherdofmapvec.assign(otherdofset.begin(), otherdofset.end());
  otherdofmap_ = Teuchos::rcp(new Epetra_Map(-1,
                                             otherdofmapvec.size(),
                                             &otherdofmapvec[0],
                                             0,
                                             dis->Comm()));
  otherimporter_ = Teuchos::rcp(new Epetra_Import(*otherdofmap_, *dofrowmap_));
  otherdofmapvec.clear();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MapExtractor::SetupMaps(Teuchos::RCP<const Epetra_Map> dofrowmap, Teuchos::RCP<Epetra_Map> condmap)
{
  dofrowmap_ = dofrowmap;
  conddofmap_ = condmap;
  otherdofmap_ = Teuchos::rcp(LINALG::SplitMap(*dofrowmap,*condmap));

  condimporter_ = Teuchos::rcp(new Epetra_Import(*conddofmap_, *dofrowmap_));
  otherimporter_ = Teuchos::rcp(new Epetra_Import(*otherdofmap_, *dofrowmap_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> DRT::UTILS::MapExtractor::ExtractCondVector(Teuchos::RCP<const Epetra_Vector> full) const
{
  Teuchos::RefCountPtr<Epetra_Vector> cond = Teuchos::rcp(new Epetra_Vector(*conddofmap_));
  ExtractCondVector(full,cond);
  return cond;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> DRT::UTILS::MapExtractor::ExtractOtherVector(Teuchos::RCP<const Epetra_Vector> full) const
{
  Teuchos::RefCountPtr<Epetra_Vector> other = Teuchos::rcp(new Epetra_Vector(*otherdofmap_));
  ExtractOtherVector(full,other);
  return other;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MapExtractor::ExtractCondVector(Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> cond) const
{
  int err = cond->Import(*full,*condimporter_,Insert);
  if (err)
    dserror("Import using importer returned err=%d",err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MapExtractor::ExtractOtherVector(Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> other) const
{
  int err = other->Import(*full,*otherimporter_,Insert);
  if (err)
    dserror("Import using importer returned err=%d",err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> DRT::UTILS::MapExtractor::InsertCondVector(Teuchos::RCP<const Epetra_Vector> cond) const
{
  Teuchos::RCP<Epetra_Vector> full = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  InsertCondVector(cond,full);
  return full;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> DRT::UTILS::MapExtractor::InsertOtherVector(Teuchos::RCP<const Epetra_Vector> other) const
{
  Teuchos::RCP<Epetra_Vector> full = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  InsertOtherVector(other,full);
  return full;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MapExtractor::InsertCondVector(Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const
{
  int err = full->Export(*cond,*condimporter_,Insert);
  if (err)
    dserror("Export using importer returned err=%d",err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MapExtractor::InsertOtherVector(Teuchos::RCP<const Epetra_Vector> other, Teuchos::RCP<Epetra_Vector> full) const
{
  int err = full->Export(*other,*otherimporter_,Insert);
  if (err)
    dserror("Export using importer returned err=%d",err);
}

#endif
