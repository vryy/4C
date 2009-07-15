/*!----------------------------------------------------------------------
\file drt_dofset.cpp
\brief A set of degrees of freedom

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
</pre>

<pre>
Maintainer: Ulrrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <iostream>
#include <algorithm>
#include <numeric>

#include "drt_dofset.H"
#include "drt_discret.H"
#include "drt_utils.H"

#include "linalg_utils.H"

#include "../headers/define_sizes.h"


/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSet::DofSet()
 : DRT::DofSetBase()
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSet::~DofSet()
{
  static_dofsets_.remove(this);
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                               ukue 04/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::DofSet& dofset)
{
  dofset.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this  (public)                                      ukue 04/07|
 *----------------------------------------------------------------------*/
void DRT::DofSet::Print(ostream& os) const
{
  for (int proc=0; proc < numdfcolelements_->Comm().NumProc(); ++proc)
  {
    if (proc == numdfcolelements_->Comm().MyPID())
    {
      if (numdfcolelements_->MyLength())
        os << "-------------------------- Proc " << proc << " :\n";
      for (int i=0; i<numdfcolelements_->MyLength(); ++i)
      {
        int numdf = (*numdfcolelements_)[i];
        int idx   = (*idxcolelements_)[i];
        os << i << ": ";
        for (int j=0; j<numdf; ++j)
          os << dofcolmap_->GID(idx+j) << " ";
        os << "\n";
      }
      os << endl;
    }
    numdfcolelements_->Comm().Barrier();
  }
  for (int proc=0; proc < numdfcolnodes_->Comm().NumProc(); ++proc)
  {
    if (proc == numdfcolnodes_->Comm().MyPID())
    {
      if (numdfcolnodes_->MyLength())
        os << "-------------------------- Proc " << proc << " :\n";
      for (int i=0; i<numdfcolnodes_->MyLength(); ++i)
      {
        int numdf = (*numdfcolnodes_)[i];
        int idx   = (*idxcolnodes_)[i];
        os << i << ": ";
        for (int j=0; j<numdf; ++j)
          os << dofcolmap_->GID(idx+j) << " ";
        os << "\n";
      }
      os << endl;
    }
    numdfcolnodes_->Comm().Barrier();
  }
}


/*----------------------------------------------------------------------*
 |  reset everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
void DRT::DofSet::Reset()
{
  dofrowmap_ = null;
  dofcolmap_ = null;
  numdfcolnodes_ = null;
  numdfcolelements_ = null;
  idxcolnodes_ = null;
  idxcolelements_ = null;
}


/*----------------------------------------------------------------------*
 |  setup everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
int DRT::DofSet::AssignDegreesOfFreedom(const Discretization& dis, const int start)
{
  if (!dis.Filled()) dserror("discretization Filled()==false");
  if (!dis.NodeRowMap()->UniqueGIDs()) dserror("Nodal row map is not unique");
  if (!dis.ElementRowMap()->UniqueGIDs()) dserror("Element row map is not unique");

  // A definite offset is currently not supported.
  if (start!=0)
    dserror("right now user specified dof offsets are not supported");

  // Add DofSets in order of assignment to list. Once it is there it has its
  // place and will get its starting id from the previous DofSet.
  AddDofSettoList();

  // We assume that all dof sets before this one have been set up. Otherwise
  // we'd have to reorder the list.
  //
  // There is no test anymore to make sure that all prior dof sets have been
  // assigned. It seems people like to manipulate dof sets. People do create
  // dof sets that do not contain any dofs (on its first assignment), people
  // even shift dof set numbers to create overlapping dof sets. This is
  // perfectly fine.
  //
  // However if you rely on non-overlapping dof sets, you have to
  // FillComplete() your discretizations in the order of their creation. This
  // is guaranteed for all discretizations read from the input file since the
  // input reader calls FillComplete(). If you create your own discretizations
  // try to understand what you do.

  // Get highest GID used so far
  int count = MaxGIDinList();

  // Now this is tricky. We have to care for nodes and elements, both
  // row and column maps. In general both nodes and elements can have
  // dofs. In both cases these dofs might be shared with other nodes
  // or elements. (The very general case. For elements we'd probably
  // don't need that.)
  //
  // The point is that we have to make sure the dof numbering of a
  // mesh is independent of its parallel distribution. Otherwise we
  // could not redistribute a mesh. We would not be able to use old
  // distributed vectors afterwards.
  //
  // Each object (node or element) could have a different number of
  // dofs. The parallel distribution is arbitrary. So we fall back to
  // two redundant vectors here to gather the number of dofs per node
  // or element.

  // numdf for all nodes and elements
  numdfcolnodes_ = rcp(new Epetra_IntVector(*dis.NodeColMap()));
  numdfcolelements_ = rcp(new Epetra_IntVector(*dis.ElementColMap()));

  // index of first dof for all nodes and elements
  idxcolnodes_ = rcp(new Epetra_IntVector(*dis.NodeColMap()));
  idxcolelements_ = rcp(new Epetra_IntVector(*dis.ElementColMap()));

  //////////////////////////////////////////////////////////////////

  // do the nodes first

  map<int,int> nidx;
  LINALG::AllreduceEMap(nidx, *dis.NodeRowMap());

  // build a redundant that holds all the node's numdof
  vector<int> sredundantnodes(dis.NumGlobalNodes());
  fill(sredundantnodes.begin(), sredundantnodes.end(), 0);

  int lnumdof = 0;

  vector<int> localnodenumdf(dis.NumMyRowNodes());

  // loop my row nodes and remember the number of degrees of freedom
  for (int i=0; i<dis.NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = dis.lRowNode(i);
    const int numele = actnode->NumElement();
    DRT::Element** myele = actnode->Elements();
    int numdf=0;
    for (int j=0; j<numele; ++j)
      numdf = max(numdf,myele[j]->NumDofPerNode(*actnode));
    const int gid = actnode->Id();
    sredundantnodes[nidx[gid]] = numdf;
    lnumdof += numdf;
  }

  // make all nodal numdfs globally known
  vector<int> rredundantnodes(dis.NumGlobalNodes());
  dis.Comm().SumAll(&sredundantnodes[0],&rredundantnodes[0],dis.NumGlobalNodes());

  int localcolpos=0;

  // Vectors to keep the dof gid of both nodes and elements. We need
  // both the row and column fashion to create the row and column map
  // later on.
  vector<int> localrowdofs;
  localrowdofs.reserve(lnumdof); // exact

  vector<int> localcoldofs;
  localcoldofs.reserve(lnumdof); // just a guess

  // We know all the numdfs of all the nodes and we even know the
  // nodes from our local row map. Use that information.
  // We have to loop the gids in order. This way we will get an
  // ordered set of dofs.
  for (map<int,int>::const_iterator i=nidx.begin(); i!=nidx.end(); ++i)
  {
    int numdf = sredundantnodes[i->second];

    // we know our row nodes from the redundant vector
    if (numdf!=0)
    {
      for (int j=0; j<numdf; ++j)
      {
        localrowdofs.push_back(count+j);
      }
    }

    // the col nodes we might have to look up
    if ((numdf!=0) || dis.NodeColMap()->MyGID(i->first))
    {
      numdf = rredundantnodes[i->second];

      // remember numdf and index for nodal based lookup
      int lid = dis.NodeColMap()->LID(i->first);
      (*numdfcolnodes_)[lid] = numdf;
      (*idxcolnodes_)[lid] = localcolpos;
      localcolpos += numdf;

      for (int j=0; j<numdf; ++j)
      {
        localcoldofs.push_back(count+j);
      }
    }
    else
    {
      // but in any case we need to increase our dof count
      numdf = rredundantnodes[i->second];
    }

    count += numdf;
  }

  sredundantnodes.clear();
  rredundantnodes.clear();
  nidx.clear();

  //////////////////////////////////////////////////////////////////

  // Now do all this fun again for the elements

  map<int,int> eidx;
  LINALG::AllreduceEMap(eidx, *dis.ElementRowMap());

  vector<int> sredundantelements(dis.NumGlobalElements());
  fill(sredundantelements.begin(), sredundantelements.end(), 0);

  vector<int> localelenumdf(dis.NumMyRowElements());

  // loop my row elements and set number of degrees of freedom
  for (int i=0; i<dis.NumMyRowElements(); ++i)
  {
    DRT::Element* actele = dis.lRowElement(i);
    const int gid = actele->Id();
    int numdf = actele->NumDofPerElement();
    sredundantelements[eidx[gid]] = numdf;
    localelenumdf[i] = numdf;
    lnumdof += numdf;
  }

  vector<int> rredundantelements(dis.NumGlobalElements());
  dis.Comm().SumAll(&sredundantelements[0],&rredundantelements[0],dis.NumGlobalElements());

  // enlarge the big dof vectors
  localrowdofs.reserve(lnumdof); // exact
  localcoldofs.reserve(lnumdof); // guess

  // We have to loop the gids in order. This way we will get an
  // ordered set of dofs.
  for (map<int,int>::const_iterator i=eidx.begin(); i!=eidx.end(); ++i)
  {
    int numdf = sredundantelements[i->second];

    // we know our row elements from the redundant vector
    if (numdf!=0)
    {
      for (int j=0; j<numdf; ++j)
      {
        localrowdofs.push_back(count+j);
      }
    }

    // the col elements we might have to look up
    if ((numdf!=0) || dis.ElementColMap()->MyGID(i->first))
    {
      numdf = rredundantelements[i->second];

      // remember numdf and index for elemental based lookup
      int lid = dis.ElementColMap()->LID(i->first);
      (*numdfcolelements_)[lid] = numdf;
      (*idxcolelements_)[lid] = localcolpos;
      localcolpos += numdf;

      for (int j=0; j<numdf; ++j)
      {
        localcoldofs.push_back(count+j);
      }
    }
    else
    {
      // but in any case we need to increase our dof count
      numdf = rredundantelements[i->second];
    }

    count += numdf;
  }

  sredundantelements.clear();
  rredundantelements.clear();
  eidx.clear();

  // Now finally we have everything in place to build the maps.

  dofrowmap_ = rcp(new Epetra_Map(-1,localrowdofs.size(),&localrowdofs[0],0,dis.Comm()));
  if (!dofrowmap_->UniqueGIDs()) dserror("Dof row map is not unique");

  dofcolmap_ = rcp(new Epetra_Map(-1,localcoldofs.size(),&localcoldofs[0],0,dis.Comm()));

  return count;
}


#endif  // #ifdef CCADISCRET
