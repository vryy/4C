
#ifdef CCADISCRET

#include "linalg_mapextractor.H"
#include "linalg_utils.H"

#include <Teuchos_getConst.hpp>

#include <numeric>

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C"
{
#include "../headers/standardtypes.h"
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LINALG::MultiMapExtractor::MultiMapExtractor()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LINALG::MultiMapExtractor::MultiMapExtractor(const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<const Epetra_Map> >& maps)
{
  Setup(fullmap,maps);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::Setup(const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<const Epetra_Map> >& maps)
{
  fullmap_ = Teuchos::rcp(new Epetra_Map(fullmap));
  maps_ = maps;

  importer_.resize(maps_.size());
  for (unsigned i=0; i<importer_.size(); ++i)
  {
    importer_[i] = Teuchos::rcp(new Epetra_Import(*maps_[i], *fullmap_));
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::MultiMapExtractor::MergeMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps)
{
  if (maps.size()==0)
    dserror("no maps to merge");
  int maplength = 0;
  for (unsigned i=0; i<maps.size(); ++i)
  {
    if (not maps[i]->UniqueGIDs())
      dserror("map %d not unique", i);
    maplength += maps[i]->NumMyElements();
  }
  std::vector<int> mapentries;
  mapentries.reserve(maplength);
  for (unsigned i=0; i<maps.size(); ++i)
  {
    const Epetra_Map& map = *maps[i];
    std::copy(map.MyGlobalElements(),
              map.MyGlobalElements()+map.NumMyElements(),
              std::back_inserter(mapentries));
  }
  return Teuchos::rcp(new Epetra_Map(-1,mapentries.size(),&mapentries[0],0,maps[0]->Comm()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> LINALG::MultiMapExtractor::ExtractVector(const Epetra_Vector& full, int block) const
{
  Teuchos::RefCountPtr<Epetra_Vector> vec = Teuchos::rcp(new Epetra_Vector(*maps_[block]));
  ExtractVector(full,block,*vec);
  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> LINALG::MultiMapExtractor::ExtractVector(const Epetra_MultiVector& full, int block) const
{
  Teuchos::RefCountPtr<Epetra_MultiVector> vec = Teuchos::rcp(new Epetra_MultiVector(*maps_[block],full.NumVectors()));
  ExtractVector(full,block,*vec);
  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::ExtractVector(const Epetra_MultiVector& full, int block, Epetra_MultiVector& partial) const
{
  int err = partial.Import(full,*importer_[block],Insert);
  if (err)
    dserror("Import using importer returned err=%d",err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> LINALG::MultiMapExtractor::InsertVector(const Epetra_Vector& partial, int block) const
{
  Teuchos::RCP<Epetra_Vector> full = Teuchos::rcp(new Epetra_Vector(*fullmap_));
  InsertVector(partial,block,*full);
  return full;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> LINALG::MultiMapExtractor::InsertVector(const Epetra_MultiVector& partial, int block) const
{
  Teuchos::RCP<Epetra_MultiVector> full = Teuchos::rcp(new Epetra_MultiVector(*fullmap_,partial.NumVectors()));
  InsertVector(partial,block,*full);
  return full;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::InsertVector(const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full) const
{
  int err = full.Export(partial,*importer_[block],Insert);
  if (err)
    dserror("Export using importer returned err=%d",err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::AddVector(const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full, double scale) const
{
  Teuchos::RCP<Epetra_MultiVector> v = ExtractVector(full, block);
  v->Update(scale,partial,1.0);
  InsertVector(*v,block,full);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LINALG::MapExtractor::MapExtractor()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LINALG::MapExtractor::MapExtractor(const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> condmap, Teuchos::RCP<const Epetra_Map> othermap)
{
  Setup(fullmap, condmap, othermap);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MapExtractor::Setup(const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> condmap, Teuchos::RCP<const Epetra_Map> othermap)
{
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(othermap);
  maps.push_back(condmap);
  MultiMapExtractor::Setup(fullmap,maps);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> LINALG::MapExtractor::ExtractCondVector(Teuchos::RCP<const Epetra_Vector> full) const
{
  return ExtractVector(full,1);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> LINALG::MapExtractor::ExtractOtherVector(Teuchos::RCP<const Epetra_Vector> full) const
{
  return ExtractVector(full,0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MapExtractor::ExtractCondVector(Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> cond) const
{
  ExtractVector(full,1,cond);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MapExtractor::ExtractOtherVector(Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> other) const
{
  ExtractVector(full,0,other);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> LINALG::MapExtractor::InsertCondVector(Teuchos::RCP<const Epetra_Vector> cond) const
{
  return InsertVector(cond,1);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> LINALG::MapExtractor::InsertOtherVector(Teuchos::RCP<const Epetra_Vector> other) const
{
  return InsertVector(other,0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MapExtractor::InsertCondVector(Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const
{
  InsertVector(cond,1,full);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MapExtractor::InsertOtherVector(Teuchos::RCP<const Epetra_Vector> other, Teuchos::RCP<Epetra_Vector> full) const
{
  InsertVector(other,0,full);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MapExtractor::AddCondVector(Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const
{
  AddVector(cond,1,full);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MapExtractor::AddOtherVector(Teuchos::RCP<const Epetra_Vector> other, Teuchos::RCP<Epetra_Vector> full) const
{
  AddVector(other,0,full);
}


#endif
