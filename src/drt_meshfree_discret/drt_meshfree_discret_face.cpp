/*!-------------------------------------------------------------------------*\
 * \file drt_meshfree_discret_face.cpp
 *
 * \brief definition of meshfree discretisation subclass Face
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*--------------------------------------------------------------------------*/

#include "drt_meshfree_discret.H"
#include "drt_meshfree_node.H"
#include "drt_meshfree_cell.H"
#include "../drt_fem_general/drt_utils_maxent_basisfunctions.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"

/*--------------------------------------------------------------------------*
 |  Face standardconstructor                                     (public) nis Jan14 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeDiscretization::Face::Face()
  : nodeids_(Teuchos::null),
    noderowptr_(0),
    nodecolmap_(),
    geometry_(Teuchos::null),
    discret_(NULL),
    dim_(-1),
    filled_(false),
    assigned_(false)
{
};

/*--------------------------------------------------------------------------*
 |  Face constructor                                     (public) nis Jan14 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeDiscretization::Face::Face(
  std::vector<int>& nodeids,
  int dim,
  DRT::MESHFREE::MeshfreeDiscretization* discret)
  : nodeids_(Teuchos::rcp(new std::vector<int>(nodeids))),
    noderowptr_(0),
    nodecolmap_(),
    geometry_(Teuchos::null),
    discret_(discret),
    dim_(dim),
    filled_(false),
    assigned_(false)
{
};

/*--------------------------------------------------------------------------*
 |  Building face node row and column pointers          (private) nis Jan14 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::Face::BuildFaceNodePointer()
{
  // get rank from discret
  const int myrank = discret_->Comm().MyPID();

  // auxiliary variables
  std::vector<int>::const_iterator it;

  // count row and column nodes to resize std::vectors
  int nummyrownodes     = 0;
  for (it=nodeids_->begin(); it!=nodeids_->end(); ++it)
    if (discret_->gNode(*it)->Owner()==myrank)
      ++nummyrownodes;

  // resize std::vectors
  noderowptr_.resize(nummyrownodes);

  // build row and column pointers
  int rowcount = 0;
  for (it=nodeids_->begin(); it!=nodeids_->end(); ++it)
  {
    DRT::Node* cnode = discret_->gNode(*it);
    // check whether node exists on proc...
    if (cnode==NULL)
      continue;
    // ... if it exists and ...
    else
    {
      // ... if node belongs to proc
      if (cnode->Owner()==myrank)
      {
        // ... add node as MeshfreeNode to row and column pointer
        MeshfreeNode* node = dynamic_cast<MeshfreeNode*>(cnode);
        if (node!=NULL)
        {
          noderowptr_[rowcount] = node;
          nodecolmap_[*it] = node;
          ++rowcount;
        }
        else
          dserror("Node %i could not be casted to MeshfreeNode or is not in discret of proc %i when building facenoderowptr!",cnode->Id(),myrank);
      }
      // ... else node is ghosted and will here be ghosted as well
      else
      {
        // ... add node as MeshfreeNode to column pointer only
        MeshfreeNode* node = dynamic_cast<MeshfreeNode*>(cnode);
        if (node!=NULL)
          nodecolmap_[*it] = node;
        else
          dserror("Node %i could not be casted to MeshfreeNode or is not in discret of proc %i when building facenodecolptr!",cnode->Id(),myrank);
      }
    }
  }

  // check if right number of nodes were added
  if (rowcount != nummyrownodes) dserror("%i row nodes added to pointer where %i were expected",rowcount,nummyrownodes);

  // we now have succesfully(?) filled to face node pointer/map
  filled_ = true;

  return;
}

/*--------------------------------------------------------------------------*
 |  Assigning column face nodes to row face nodes       (private) nis Jan14 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::Face::AssignNodesToNodes(
  bool setnodefacedim)
{
  // get range from discret
  const double range = discret_->GetSolutionApprox()->GetRange();

  // if vertex: check of ownership
  if (dim_==0)
  {
    // if own vertex
    if (noderowptr_.size())
    {
      // clear node and add node to itself
      noderowptr_[0]->ClearMyNodeTopology();
      noderowptr_[0]->AddNodePtr(noderowptr_[0]);
      if (setnodefacedim)
        noderowptr_[0]->SetFacedim(0);
      if (noderowptr_[0]->GetFaceDim()>dim_)
        dserror("Node face dim of node %i is %i on a %i-D face.",noderowptr_[0]->Id(),noderowptr_[0]->GetFaceDim(),0);
    }
  }
  // if line: brute force search on row nodes
  else if (dim_==1)
  {
    // get iterators of row and column nodes
    std::vector<MeshfreeNode*>::iterator rnode;
    std::map<int,MeshfreeNode*>::iterator cnode;
    // loop over row nodes
    for(rnode=noderowptr_.begin(); rnode!=noderowptr_.end(); ++rnode)
    {
      // check if node face dim is equal to face dim or, if larger, to be set
      if (((*rnode)->GetFaceDim()==dim_) or (((*rnode)->GetFaceDim()>dim_) and (setnodefacedim)))
      {
        // clear node topology
        (*rnode)->ClearMyNodeTopology();
        for(cnode=nodecolmap_.begin(); cnode!=nodecolmap_.end(); ++cnode)
        {
          // check distance to column nodes and add if in range
          if ((*rnode)->DistToPoint(cnode->second->X())<range)
            (*rnode)->AddNodePtr(cnode->first,cnode->second);
        }
        // set node face dim if required
        if (setnodefacedim)
          (*rnode)->SetFacedim(dim_);
      }
      // else, if node face dim is larger and set, something went wrong
      else if ((*rnode)->GetFaceDim()>dim_)
        dserror("Node face dim of node %i is %i on a %i-D face.", (*rnode)->Id(),(*rnode)->GetFaceDim(),dim_);
      // if node face dim is smaller, node is treated on lower dimensional face
    }
  }
  // if surface or volume: searchtree - prerequisite: correct ghosting
  else if (dim_==2 or dim_==3)
  {
    // get iterators of row and column nodes
    std::vector<MeshfreeNode*>::iterator rnode;
    std::map<int,MeshfreeNode*>::iterator cnode;

    // get node coordinates searchtree-style
    std::map<int,LINALG::Matrix<3,1> > nodepositions;
    for (cnode=nodecolmap_.begin(); cnode!=nodecolmap_.end(); ++cnode)
      nodepositions.insert( std::pair<int,LINALG::Matrix<3,1> >(cnode->first,LINALG::Matrix<3,1>(cnode->second->X())));

    // set up searchtree - always OCTREE although for surfaces a QUADTREE would be optimal
    GEO::SearchTree searchTree(8);
    searchTree.initializePointTree(GEO::getXAABBofPositions(nodepositions), nodepositions, GEO::TreeType(GEO::OCTTREE));

    // loop over row nodes
    for(rnode=noderowptr_.begin(); rnode!=noderowptr_.end(); ++rnode)
    {
      // check if node face dim is equal to face dim or, if larger, to be set
      if (((*rnode)->GetFaceDim()==dim_) or (((*rnode)->GetFaceDim()>dim_) and (setnodefacedim)))
      {
        // clear node topology
        (*rnode)->ClearMyNodeTopology();
        // find neighbours of current row node
        const LINALG::Matrix<3,1> rnodexyz((*rnode)->X());
        std::vector<int> neighbours = searchTree.searchPointsInRadius(nodepositions, rnodexyz, range);
        (*rnode)->SetMyNodeTopology(neighbours,nodecolmap_);
        // set node face dim if required
        if (setnodefacedim)
          (*rnode)->SetFacedim(dim_);
      }
      // else, if node face dim is larger and set, something went wrong
      else if ((*rnode)->GetFaceDim()>dim_)
        dserror("Node face dim of node %i is %i on a %i-D face.", (*rnode)->Id(), (*rnode)->GetFaceDim(), dim_);
      // if node face dim is smaller, node is treated on lower dimensional face
    }
  }
  else
    dserror("Face dimension is %i but has to be 0 (vertex), 1 (line), 2 (surface), or 3 (volume).",dim_);

  // we now have succesfully(?) assigned nodes to nodes
  assigned_ = true;

  return;
}

/*----------------------------------------------------------------------------*
 |  Assign nodes to cells on face                                   nis Jan14 |
 *----------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::Face::AssignNodesToCells()
{
  // get rank and range from discret
  const int myrank = discret_->Comm().MyPID();
  const double range = discret_->GetSolutionApprox()->GetRange();

  // initialize variables for neighbourhood-check
  LINALG::Matrix<3,1> center(true);
  double maxcellradius(0.0);

  if (geometry_->size()==0) dserror("No geometry exists on this face.");

  switch (geometry_->begin()->second->Shape())
  {
    case DRT::Element::line2:
    {
      // loop over all cells
      for (std::map<int,Teuchos::RCP<DRT::Element> >::iterator it = geometry_->begin(); it!=geometry_->end(); ++it)
      {
        // only deal with own cells
        if (it->second->Owner()!=myrank) continue;

        // create vectors to hold node ids and pointers
        std::vector<int> condnodeids(0);
        std::vector<DRT::Node*> condnodes(0);

        // get cell center and maximum radius
        DRT::MESHFREE::Cell* ccell = dynamic_cast<DRT::MESHFREE::Cell*>(it->second.get());
        if (ccell==NULL)
          dserror("Something went wrong when casting element to meshfree cell.");
        ccell->CenterAndMaxRadius(center, maxcellradius);
        const double truerange = range+maxcellradius;

        // brute force search for nodes in neighbourhood of lines
        std::map<int,MeshfreeNode*>::iterator cnode;
        for (cnode=nodecolmap_.begin(); cnode!=nodecolmap_.end(); ++cnode)
        {
          if (cnode->second->DistToPoint(center.A())<truerange)
          {
            condnodeids.push_back(cnode->first);
            condnodes.push_back(cnode->second);
          }
        }

        // set node ids and nodal pointers of cell
        ccell->DeleteNodes();
        ccell->SetNodeIds(condnodeids.size(), condnodeids.data());
        ccell->BuildNodalPointers(condnodes.data());
      }
      break;
    }
    case DRT::Element::tri3:
    case DRT::Element::quad4:
    case DRT::Element::tet4:
    case DRT::Element::hex8:
    {
      // get iterators of row and column nodes
      std::vector<MeshfreeNode*>::iterator rnode;
      std::map<int,MeshfreeNode*>::iterator cnode;

      // get node coordinates searchtree-style
      std::map<int,LINALG::Matrix<3,1> > nodepositions;
      for (cnode=nodecolmap_.begin(); cnode!=nodecolmap_.end(); ++cnode)
        nodepositions.insert( std::pair<int,LINALG::Matrix<3,1> >(cnode->first,LINALG::Matrix<3,1>(cnode->second->X())));

      // set up searchtree - always OCTREE although for surfaces a QUADTREE would be optimal
      GEO::SearchTree searchTree(8);
      searchTree.initializePointTree(GEO::getXAABBofPositions(nodepositions), nodepositions, GEO::TreeType(GEO::OCTTREE));

      // loop over all cells
      for (std::map<int,Teuchos::RCP<DRT::Element> >::iterator it = geometry_->begin(); it!=geometry_->end(); ++it)
      {
        // only deal with own cells
        if (it->second->Owner()!=myrank) continue;

        // get cell center and maximum radius
        DRT::MESHFREE::Cell* ccell = dynamic_cast<DRT::MESHFREE::Cell*>(it->second.get());
        if (ccell==NULL)
          dserror("Something went wrong when casting element to meshfree cell.");
        // find neighbours of current row node
        ccell->CenterAndMaxRadius(center, maxcellradius);
        std::vector<int> condnodeids = searchTree.searchPointsInRadius(nodepositions, center, range+maxcellradius);
        std::vector<DRT::Node*> condnodes(condnodeids.size());
        for (unsigned i=0; i<condnodeids.size(); ++i)
          condnodes[i] = nodecolmap_[condnodeids[i]];

        // set node ids and nodal pointers of cell
        ccell->DeleteNodes();
        ccell->SetNodeIds(condnodeids.size(), condnodeids.data());
        ccell->BuildNodalPointers(condnodes.data());
      }
      break;
    }
    default:
      dserror("Unknown discretization type for meshfree discretisation in neighbourhood search.");
  }

  return;
}
