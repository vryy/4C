/*!----------------------------------------------------------------------
\file drt_dofset.cpp
\brief A set of degrees of freedom

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <iostream>
#include <algorithm>
#include <numeric>

#include "drt_dofset.H"
#include "drt_discret.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::DofSet::DofSet()
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::DofSet::~DofSet()
{
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::DofSet& dofset)
{
  return os;
}


/*----------------------------------------------------------------------*
 |  print this  (public)                                     mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::DofSet::Print(ostream& os) const
{
  return;
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::DofSet::FindMyPos(int myrank, int numproc, const Epetra_Map& emap)
{
  vector<int> snum(numproc);
  vector<int> rnum(numproc);
  fill(snum.begin(), snum.end(), 0);
  snum[myrank] = emap.NumMyElements();

  emap.Comm().SumAll(&snum[0],&rnum[0],numproc);

  return std::accumulate(&rnum[0], &rnum[myrank], 0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::DofSet::AllreduceEMap(map<int,int>& idxmap, const Epetra_Map& emap)
{
  const int myrank  = emap.Comm().MyPID();
  const int numproc = emap.Comm().NumProc();

  idxmap.clear();

  int mynodepos = FindMyPos(myrank, numproc, emap);

  vector<int> sredundant(emap.NumGlobalElements());
  vector<int> rredundant(emap.NumGlobalElements());
  fill(sredundant.begin(), sredundant.end(), 0);

  int* gids = emap.MyGlobalElements();
  copy(gids, gids+emap.NumMyElements(), &sredundant[mynodepos]);

  emap.Comm().SumAll(&sredundant[0], &rredundant[0], emap.NumGlobalElements());
  sredundant.clear();

  for (unsigned i=0; i<rredundant.size(); ++i)
  {
    idxmap[rredundant[i]] = i;
  }
}


/*----------------------------------------------------------------------*
 |  setup everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
int DRT::DofSet::AssignDegreesOfFreedom(const Discretization& dis, const int start)
{
  if (!dis.Filled()) dserror("discretization Filled()==false");
  if (!dis.NodeRowMap()->UniqueGIDs()) dserror("Nodal row map is not unique");
  if (!dis.ElementRowMap()->UniqueGIDs()) dserror("Element row map is not unique");

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
  AllreduceEMap(nidx, *dis.NodeRowMap());

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

  int count=start;
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
  for (unsigned i=0; i<sredundantnodes.size(); ++i)
  {
    int numdf = sredundantnodes[nidx[i]];

    // we know our row nodes from the redundant vector
    if (numdf!=0)
    {
      for (int j=0; j<numdf; ++j)
      {
        localrowdofs.push_back(count+j);
      }
    }

    // the col nodes we might have to look up
    if ((numdf!=0) || dis.NodeColMap()->MyGID(i))
    {
      numdf = rredundantnodes[nidx[i]];

      // remember numdf and index for nodal based lookup
      int lid = dis.NodeColMap()->LID(i);
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
      numdf = rredundantnodes[nidx[i]];
    }

    count += numdf;
  }

  sredundantnodes.clear();
  rredundantnodes.clear();
  nidx.clear();

  //////////////////////////////////////////////////////////////////

  // Now do all this fun again for the elements

  map<int,int> eidx;
  AllreduceEMap(eidx, *dis.ElementRowMap());

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

  for (unsigned i=0; i<sredundantelements.size(); ++i)
  {
    int numdf = sredundantelements[eidx[i]];

    // we know our row elements from the redundant vector
    if (numdf!=0)
    {
      for (int j=0; j<numdf; ++j)
      {
        localrowdofs.push_back(count+j);
      }
    }

    // the col elements we might have to look up
    if ((numdf!=0) || dis.ElementColMap()->MyGID(i))
    {
      numdf = rredundantelements[eidx[i]];

      // remember numdf and index for elemental based lookup
      int lid = dis.ElementColMap()->LID(i);
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
      numdf = rredundantelements[eidx[i]];
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


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
