/*!----------------------------------------------------------------------
\file drt_xdofset.H

\brief dofset to manage the varying number and the varying meaning of 
       XFEM degrees of freedom

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifdef D_XFEM

#include "../drt_lib/drt_pbcdofset.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../headers/define_sizes.h"
#include <drt_xdofset.H>





/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/07|
 *----------------------------------------------------------------------*/
DRT::XDofSet::XDofSet(RefCountPtr<map<int,vector<int> > >  couplednodes)
{
  perbndcouples_=couplednodes;
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 05/07|
 *----------------------------------------------------------------------*/
DRT::XDofSet::~XDofSet()
{
  return;
}

/*----------------------------------------------------------------------*
 |  reset everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
void DRT::XDofSet::Reset()
{
  dofrowmap_ = null;
  dofcolmap_ = null;
  numdfcolnodes_ = null;
  numdfcolelements_ = null;
  idxcolnodes_ = null;
  idxcolelements_ = null;
}

/*----------------------------------------------------------------------*
 |  setup everything  (public)                               gammi 05/07|
 |  this function is a specialisation of the AssignDegreesOfFreedom of  |
 |  the base class DofSet                                               |
 *----------------------------------------------------------------------*/
int DRT::XDofSet::AssignDegreesOfFreedom(const Discretization& dis, const int start)
{
  dsassert(dis.Filled(), "discretization Filled()==false");
  dsassert(dis.NodeRowMap()->UniqueGIDs(), "Nodal row map is not unique");
  dsassert(dis.ElementRowMap()->UniqueGIDs(), "Element row map is not unique");

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
  Utils::AllreduceEMap(nidx, *dis.NodeRowMap());

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

  for (unsigned i=0; i<rredundantnodes.size(); ++i)
    if (rredundantnodes[i] > MAXDOFPERNODE)
      dserror("MAXDOFPERNODE=%d and numdf=%d found", MAXDOFPERNODE, rredundantnodes[i]);

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
  // We have to loop the gids in order. This way we will get an
  // ordered set of dofs.
  for (map<int,int>::iterator i=nidx.begin(); i!=nidx.end(); ++i)
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
  Utils::AllreduceEMap(eidx, *dis.ElementRowMap());

  vector<int> sredundantelements(dis.NumGlobalElements());
  fill(sredundantelements.begin(), sredundantelements.end(), 0);

  vector<int> localelenumdf(dis.NumMyRowElements());

  // loop my row elements and set number of degrees of freedom
  for (int i=0; i<dis.NumMyRowElements(); ++i)
  {
    const DRT::Element* actele = dis.lRowElement(i);
    const int gid = actele->Id();
    const int numdf = actele->NumDofPerElement();
    sredundantelements[eidx[gid]] = numdf;
    localelenumdf[i] = numdf;
    lnumdof += numdf;
  }

  vector<int> rredundantelements(dis.NumGlobalElements());
  dis.Comm().SumAll(&sredundantelements[0],&rredundantelements[0],dis.NumGlobalElements());

  for (unsigned i=0; i<rredundantelements.size(); ++i)
    if (rredundantelements[i] > MAXDOFPERNODE)
      dserror("MAXDOFPERNODE=%d and numdf=%d found", MAXDOFPERNODE, rredundantelements[i]);

  // enlarge the big dof vectors
  localrowdofs.reserve(lnumdof); // exact
  localcoldofs.reserve(lnumdof); // guess

  // We have to loop the gids in order. This way we will get an
  // ordered set of dofs.
  for (map<int,int>::iterator i=eidx.begin(); i!=eidx.end(); ++i)
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


#endif  // #ifdef D_XFEM
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
