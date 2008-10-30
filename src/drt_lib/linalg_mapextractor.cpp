/*----------------------------------------------------------------------*/
/*!
\file linalg_mapextractor.cpp

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de) 
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de                   

-------------------------------------------------------------------------
<\pre>
<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "linalg_mapextractor.H"
#include "linalg_utils.H"

#include <Teuchos_getConst.hpp>

#include <numeric>

#ifdef PARALLEL
#include <mpi.h>
#endif

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
    if (maps_[i]!=Teuchos::null)
    {
      importer_[i] = Teuchos::rcp(new Epetra_Import(*maps_[i], *fullmap_));
    }
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
    if (maps[i]==Teuchos::null)
      dserror("can not merge extractor with null maps");
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
  if (maps_[block]==Teuchos::null)
    dserror("null map at block %d",block);
  Teuchos::RefCountPtr<Epetra_Vector> vec = Teuchos::rcp(new Epetra_Vector(*maps_[block]));
  ExtractVector(full,block,*vec);
  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> LINALG::MultiMapExtractor::ExtractVector(const Epetra_MultiVector& full, int block) const
{
  if (maps_[block]==Teuchos::null)
    dserror("null map at block %d",block);
  Teuchos::RefCountPtr<Epetra_MultiVector> vec = Teuchos::rcp(new Epetra_MultiVector(*maps_[block],full.NumVectors()));
  ExtractVector(full,block,*vec);
  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::ExtractVector(const Epetra_MultiVector& full, int block, Epetra_MultiVector& partial) const
{
  if (maps_[block]==Teuchos::null)
    dserror("null map at block %d",block);
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
  if (maps_[block]==Teuchos::null)
    dserror("null map at block %d",block);
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
LINALG::MapExtractor::MapExtractor(const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> condmap)
{
  // initialise other DOFs by inserting all DOFs of full map
  std::set<int> othergids;
  const int* fullgids = fullmap.MyGlobalElements();
  for (int lid=0; lid<fullmap.NumMyElements(); ++lid)
    othergids.insert(fullgids[lid]);
  // throw away all DOFs which are in condmap
  const int* condgids = (*condmap).MyGlobalElements();
  for (int lid=0; lid<(*condmap).NumMyElements(); ++lid)
    othergids.erase(condgids[lid]);
  // create (non-overlapping) othermap for non-condmap DOFs
  Teuchos::RCP<Epetra_Map> othermap = Teuchos::null;
  if (othergids.size() > 0)
  {
    int* othergidsv = new int[othergids.size()];
    int lid = 0;
    for (std::set<int>::iterator gid=othergids.begin(); gid!=othergids.end(); ++gid)
    {
      othergidsv[lid] = *gid;
      lid += 1;
    }
    othermap = Teuchos::rcp(new Epetra_Map(-1, othergids.size(), othergidsv, fullmap.IndexBase(), fullmap.Comm()));
    delete [] othergidsv;
  }
  else
  {
    othermap = Teuchos::rcp(new Epetra_Map(0, fullmap.IndexBase(), fullmap.Comm()));
  }
  // create the extractor
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
void LINALG::MapExtractor::AddCondVector(double scale, Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const
{
  AddVector(cond,1,full,scale);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MapExtractor::AddOtherVector(Teuchos::RCP<const Epetra_Vector> other, Teuchos::RCP<Epetra_Vector> full) const
{
  AddVector(other,0,full);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MapExtractor::AddOtherVector(double scale, Teuchos::RCP<const Epetra_Vector> other, Teuchos::RCP<Epetra_Vector> full) const
{
  AddVector(other,0,full,scale);
}


#endif
