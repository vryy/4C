/*!----------------------------------------------------------------------
\file drt_dofset_independent.cpp

\brief Implementation
\level 2
\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/

#include "drt_dofset_independent.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 *----------------------------------------------------------------------*/
DRT::IndependentDofSet::IndependentDofSet(bool ignoreminnodegid/*=false*/)
  : DRT::DofSet(),
    ignoreminnodegid_(ignoreminnodegid)
{
  return;
}


/*----------------------------------------------------------------------*
 |  cctor (public)                                                      |
 *----------------------------------------------------------------------*/
DRT::IndependentDofSet::IndependentDofSet(const IndependentDofSet& old)
  : DRT::DofSet(old),
    ignoreminnodegid_(old.ignoreminnodegid_)
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                                       |
 *----------------------------------------------------------------------*/
DRT::IndependentDofSet::~IndependentDofSet()
{
  return;
}


int DRT::IndependentDofSet::AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start)
{
  //========================================================================
  // This method has to be completely reimplemented in order to avoid
  // registering with the static_dofsets_ and starting the numbering at 0
  //========================================================================

  if (!dis.Filled()) dserror("discretization Filled()==false");
  if (!dis.NodeRowMap()->UniqueGIDs()) dserror("Nodal row map is not unique");
  if (!dis.ElementRowMap()->UniqueGIDs()) dserror("Element row map is not unique");

  // A definite offset is currently not supported.
  // TODO (kronbichler) find a better solution for this
  //if (start!=0)
  //  dserror("right now user specified dof offsets are not supported");

  dspos_ = dspos;

  // the starting GID is always 0
  int count = 0;

  // Now this is tricky. We have to care for nodes and elements, both in
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
    numdfrownodes[i] = NumDofPerNode(*actnode);
  }

  int minnodegid = dis.NodeRowMap()->MinAllGID();
  if(ignoreminnodegid_ == true)
    minnodegid = 0;
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
    int numdf = NumDofPerElement(*actele);
    //int numdf = dis.NumDof(dspos,actele);
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

