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

#include <iostream>
#include <algorithm>
#include <numeric>

#include "drt_dofset.H"
#include "drt_dofset_proxy.H"
#include "drt_discret.H"
#include "drt_utils.H"

#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSet::DofSet()
  : DRT::DofSetBase(), filled_(false), dspos_( 0 )
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSet::~DofSet()
{
  for (std::list<DofSetProxy*>::iterator i=proxies_.begin(); i!=proxies_.end(); ++i)
  {
    (*i)->Disconnect(this);
  }
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
          os << (idx+j) << " ";
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
          os << (idx+j) << " ";
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

  filled_ = false;

  // tell all proxies
  NotifyReset();
}


/*----------------------------------------------------------------------*
 |  setup everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
int DRT::DofSet::AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start)
{
  if (!dis.Filled()) dserror("discretization Filled()==false");
  if (!dis.NodeRowMap()->UniqueGIDs()) dserror("Nodal row map is not unique");
  if (!dis.ElementRowMap()->UniqueGIDs()) dserror("Element row map is not unique");

  // A definite offset is currently not supported.
  if (start!=0)
    dserror("right now user specified dof offsets are not supported");

  dspos_ = dspos;

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

  // Get highest GID used so far and add one
  int count = MaxGIDinList(dis.Comm()) + 1;

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
  numdfcolnodes_ = Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));
  numdfcolelements_ = Teuchos::rcp(new Epetra_IntVector(*dis.ElementColMap()));

  // index of first dof for all nodes and elements
  idxcolnodes_ = Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));
  idxcolelements_ = Teuchos::rcp(new Epetra_IntVector(*dis.ElementColMap()));

  //////////////////////////////////////////////////////////////////

  // do the nodes first

  Epetra_IntVector numdfrownodes(*dis.NodeRowMap());
  Epetra_IntVector idxrownodes(*dis.NodeRowMap());

  int numrownodes = dis.NumMyRowNodes();
  for (int i=0; i<numrownodes; ++i)
  {
    DRT::Node* actnode = dis.lRowNode(i);
    //const int gid = actnode->Id();
    numdfrownodes[i] = NumDofPerNode(*actnode,dspos);
  }

  int minnodegid = dis.NodeRowMap()->MinAllGID();
  int maxnodenumdf = numdfrownodes.MaxValue();

  std::map<int,std::vector<int> > nodedofset;

  for (int i=0; i<numrownodes; ++i)
  {
    DRT::Node* actnode = dis.lRowNode(i);
    const int gid = actnode->Id();
    int numdf = numdfrownodes[i];
    int dof = count + ( gid-minnodegid )*maxnodenumdf;
    idxrownodes[i] = dof;
    std::vector<int> & dofs = nodedofset[gid];
    dofs.reserve( numdf );
    for ( int j=0; j<numdf; ++j )
    {
      dofs.push_back( dof+j );
    }
  }

  Epetra_Import nodeimporter( numdfcolnodes_->Map(), numdfrownodes.Map() );
  int err = numdfcolnodes_->Import( numdfrownodes, nodeimporter, Insert );
  if (err) dserror( "Import using importer returned err=%d", err );
  err = idxcolnodes_->Import( idxrownodes, nodeimporter, Insert );
  if (err) dserror( "Import using importer returned err=%d", err );

  count = idxrownodes.MaxValue() + maxnodenumdf + 1;

  //////////////////////////////////////////////////////////////////

  // Now do it again for the elements

  Epetra_IntVector numdfrowelements(*dis.ElementRowMap());
  Epetra_IntVector idxrowelements(*dis.ElementRowMap());

  int numrowelements = dis.NumMyRowElements();
  for (int i=0; i<numrowelements; ++i)
  {
    DRT::Element* actele = dis.lRowElement(i);
    //const int gid = actele->Id();
    int numdf = NumDofPerElement(*actele,dspos);
    numdfrowelements[i] = numdf;
  }

  int minelementgid = dis.ElementRowMap()->MinAllGID();
  int maxelementnumdf = numdfrowelements.MaxValue();

  std::map<int,std::vector<int> > elementdofset;

  for (int i=0; i<numrowelements; ++i)
  {
    DRT::Element* actelement = dis.lRowElement(i);
    const int gid = actelement->Id();
    int numdf = numdfrowelements[i];
    int dof = count + ( gid-minelementgid )*maxelementnumdf;
    idxrowelements[i] = dof;
    std::vector<int> & dofs = elementdofset[gid];
    dofs.reserve( numdf );
    for ( int j=0; j<numdf; ++j )
    {
      dofs.push_back( dof+j );
    }
  }

  Epetra_Import elementimporter( numdfcolelements_->Map(), numdfrowelements.Map() );
  err = numdfcolelements_->Import( numdfrowelements, elementimporter, Insert );
  if (err) dserror( "Import using importer returned err=%d", err );
  err = idxcolelements_->Import( idxrowelements, elementimporter, Insert );
  if (err) dserror( "Import using importer returned err=%d", err );

  // Now finally we have everything in place to build the maps.

  std::vector<int> localrowdofs;
  std::vector<int> localcoldofs;
  localrowdofs.reserve( numrownodes*maxnodenumdf + numrowelements*maxelementnumdf );
  localcoldofs.reserve( numrownodes*maxnodenumdf + numrowelements*maxelementnumdf );

  for ( std::map<int,std::vector<int> >::iterator i=nodedofset.begin();
        i!=nodedofset.end();
        ++i )
  {
    std::vector<int> & dofs = i->second;
    std::copy( dofs.begin(), dofs.end(), std::back_inserter( localrowdofs ) );
  }
  for ( std::map<int,std::vector<int> >::iterator i=elementdofset.begin();
        i!=elementdofset.end();
        ++i )
  {
    std::vector<int> & dofs = i->second;
    std::copy( dofs.begin(), dofs.end(), std::back_inserter( localrowdofs ) );
  }

  Exporter nodeexporter( *dis.NodeRowMap(), *dis.NodeColMap(), dis.Comm() );
  nodeexporter.Export( nodedofset );

  Exporter elementexporter( *dis.ElementRowMap(), *dis.ElementColMap(), dis.Comm() );
  elementexporter.Export( elementdofset );

  for ( std::map<int,std::vector<int> >::iterator i=nodedofset.begin();
        i!=nodedofset.end();
        ++i )
  {
    std::vector<int> & dofs = i->second;
    std::copy( dofs.begin(), dofs.end(), std::back_inserter( localcoldofs ) );
  }
  for ( std::map<int,std::vector<int> >::iterator i=elementdofset.begin();
        i!=elementdofset.end();
        ++i )
  {
    std::vector<int> & dofs = i->second;
    std::copy( dofs.begin(), dofs.end(), std::back_inserter( localcoldofs ) );
  }

  dofrowmap_ = Teuchos::rcp(new Epetra_Map(-1,localrowdofs.size(),&localrowdofs[0],0,dis.Comm()));
  if (!dofrowmap_->UniqueGIDs()) dserror("Dof row map is not unique");

  dofcolmap_ = Teuchos::rcp(new Epetra_Map(-1,localcoldofs.size(),&localcoldofs[0],0,dis.Comm()));

  filled_ = true;

  // tell all proxies
  NotifyAssigned();

  return count;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSet::RegisterProxy(DofSetProxy* proxy)
{
  proxies_.push_back(proxy);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSet::UnregisterProxy(DofSetProxy* proxy)
{
  proxies_.remove(proxy);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSet::NotifyAssigned()
{
  for (std::list<DofSetProxy*>::iterator i=proxies_.begin(); i!=proxies_.end(); ++i)
    (*i)->NotifyAssigned();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSet::NotifyReset()
{
  for (std::list<DofSetProxy*>::iterator i=proxies_.begin(); i!=proxies_.end(); ++i)
    (*i)->NotifyReset();
}


