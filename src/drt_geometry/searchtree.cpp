 /*!
\file searchtree.cpp

\brief provides a class with search tree

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 */

#ifdef CCADISCRET
#include "searchtree.H"
#include "intersection_service.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/standardtypes_cpp.H"
#include <Teuchos_TimeMonitor.hpp>



/*----------------------------------------------------------------------*
 | constructor SearchTree                                    u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::SearchTree::SearchTree(
    const int           max_depth): 
    max_depth_(max_depth),
    treeRoot_(NULL)      
{}
    


/*----------------------------------------------------------------------*
| destructor SearchTree                                     u.may 07/08|
*----------------------------------------------------------------------*/
GEO::SearchTree::~SearchTree(){}



/*----------------------------------------------------------------------*
 | initialize or rebuild tree with possibly new              u.may 07/08| 
 | discretization, elements are sorted according to the given map       |
 *----------------------------------------------------------------------*/
void GEO::SearchTree::initializeTree(
    const LINALG::Matrix<3,2>&            nodeBox,
    const std::map<int,std::set<int> >&   elementsByLabel,
    const TreeType                        treetype) 
{
  
  
  treeRoot_ = Teuchos::null;
  treeRoot_ = rcp(new TreeNode(NULL, max_depth_, nodeBox, treetype)); 

  // insert element map into tree root node
  if(elementsByLabel.size()>0)
    treeRoot_->setElementList(elementsByLabel);
  
}



/*----------------------------------------------------------------------*
 | initialize or rebuild tree with possibly new              u.may 07/08| 
 | discretization, elements are taken unsorted from discretization      |
 *----------------------------------------------------------------------*/
void GEO::SearchTree::initializeTree(
    const LINALG::Matrix<3,2>&  nodeBox,
    const DRT::Discretization&  dis,
    const TreeType              treetype) 
{

  treeRoot_ = Teuchos::null;
  treeRoot_ = rcp( new TreeNode(NULL, max_depth_, nodeBox, treetype)); 
 
  // inserts all elements in a map with key -1 and global id
  for (int i=0; i<dis.NumMyColElements(); i++)
    treeRoot_->insertElement(-1, dis.lColElement(i)->Id());
}



/*-----------------------------------------------------------------------*
 | initialize or rebuild tree with possibly new              cyron 04/09 | 
 | discretization, nodes are taken unsortedly from discretization; note: |
 | this initialization does not work with any specific data of the nodes;|
 | rather it just generally organizes a tree whose root node has as many |
 | tree elements  as the discretization nodes; any more specific action  |
 | is taken not before a specific query is applied                       |
 *----------------------------------------------------------------------*/
void GEO::SearchTree::initializePointTree(
    const LINALG::Matrix<3,2>&  nodeBox,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const TreeType              treetype) 
{

  treeRoot_ = Teuchos::null;
  treeRoot_ = rcp( new TreeNode(NULL, max_depth_, nodeBox, treetype)); 
    
  /*as each node of the tree the root node stores its elements in the class
   * variable elementList_, which is of type std::map<int, std::set<int> >;
   * thus this list allows not only for saving identification numbers for a list
   * of element assigned to the node; in addition to that it allows for saving 
   * several lists of such identification numbers separated by the label which is
   * stored in the int-component of the map; here we do not make use of this 
   * possibility and just apply the label -1 and one single related set which
   * contains all the identification numbers of all the elements assinged to
   * this tree node; the identification numbers of all points stores in
   * current positions are given by the first components of this map, respectively*/  
  for(std::map<int,LINALG::Matrix<3,1> >::const_iterator mapit = currentpositions.begin(); mapit != currentpositions.end(); mapit++)
    treeRoot_->insertElement(-1,mapit->first);

}



/*----------------------------------------------------------------------*
 | update tree                                               u.may 08/08|
 | only usefull for labeled structures,for networks use initializeTree  |
 *----------------------------------------------------------------------*/
void GEO::SearchTree::updateTree(
    const DRT::Discretization& 	dis,
    const std::map<int,LINALG::Matrix<3,1> >& 	currentpositions_old, 
    const std::map<int,LINALG::Matrix<3,1> >& 	currentpositions_new) 
{
  if(treeRoot_ == Teuchos::null)
    dserror("tree is not yet initialized !!!");
  
  std::vector< LINALG::Matrix<3,2> > AABBs_old = GEO::computeXAABBForLabeledStructures(dis, currentpositions_old, treeRoot_->getElementList());
  std::vector< LINALG::Matrix<3,2> > AABBs_new = GEO::computeXAABBForLabeledStructures(dis, currentpositions_new, treeRoot_->getElementList());

  for(unsigned int i = 0; i < AABBs_new.size(); i++)
    treeRoot_->updateTreeNode(AABBs_old[i], AABBs_new[i]);
}



/*----------------------------------------------------------------------*
 | returns xfem label and nearest object to point            u.may 07/08|
 *----------------------------------------------------------------------*/
int GEO::SearchTree::queryFSINearestObject(
    const DRT::Discretization&                  dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const std::map<int,LINALG::Matrix<3,2> >&   currentXAABBs,
    const LINALG::Matrix<3,1>&                  point,
    GEO::NearestObject&                         nearestobject)
{
  
  TEUCHOS_FUNC_TIME_MONITOR("GEO::SearchTree - queryTime");
  
  if(treeRoot_ == Teuchos::null)
    dserror("tree is not yet initialized !!!");

  if(!treeRoot_->getElementList().empty())
    return treeRoot_->queryFSINearestObject(dis, currentpositions, currentXAABBs, point, nearestobject);
  else 
    return 0;
}



/*----------------------------------------------------------------------*
 | returns xfem label of point                               u.may 07/08|
 *----------------------------------------------------------------------*/
int GEO::SearchTree::queryXFEMFSIPointType(
    const DRT::Discretization&                  dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const std::map<int,LINALG::Matrix<3,2> >&   currentXAABBs,
    const LINALG::Matrix<3,1>&                  point)
{
  
  TEUCHOS_FUNC_TIME_MONITOR("GEO::SearchTree - queryTime");
  
  if(treeRoot_ == Teuchos::null)
    dserror("tree is not yet initialized !!!");

  if(!treeRoot_->getElementList().empty())
    return treeRoot_->queryXFEMFSIPointType(dis, currentpositions, currentXAABBs, point);
  else 
    return 0;
}



/*----------------------------------------------------------------------*
 | fix intersection for contact with xfem FSI              u.may   02/09|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::moveContactNodes(
    const vector<vector<int> >&                   triangleList,
    vector< GEO::InterfacePoint >&                pointList, 
    const std::map<int,LINALG::Matrix<3,2> >&     triangleXAABBs,  
    const LINALG::Matrix<3,1>&                    querypoint,
    const int                                     querypointId)
{
  
  TEUCHOS_FUNC_TIME_MONITOR("GEO::SearchTree - queryTime");
  
  if(treeRoot_ == Teuchos::null)
    dserror("tree is not yet initialized !!!");

  if(!treeRoot_->getElementList().empty())
    treeRoot_->moveContactNodes( triangleList, pointList, triangleXAABBs,  querypoint, querypointId);

  return;
}



/*----------------------------------------------------------------------*
 | returns nodes in the radius of a given point              u.may 07/08|
 *----------------------------------------------------------------------*/
std::map<int,std::set<int> > GEO::SearchTree::searchElementsInRadius(
    const DRT::Discretization& 	                dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions, 
    const LINALG::Matrix<3,1>&                  point,
    const double                                radius, 
    const int                                   label) 
{
  TEUCHOS_FUNC_TIME_MONITOR("SearchTree - queryTime");

  std::map<int,std::set<int> > nodeset;

  if(treeRoot_ == Teuchos::null)
    dserror("tree is not yet initialized !!!");

  if(!(treeRoot_->getElementList().empty()))
    nodeset = treeRoot_->searchElementsInRadius(dis, currentpositions, point, radius, label);
  else
    dserror("element list is empty");

  return nodeset;
}


/*----------------------------------------------------------------------*
 | returns a vector of elements whose XAABB s are            u.may 09/08|
 | intersecting with the XAABB of a given volume element                |
 *----------------------------------------------------------------------*/
void GEO::SearchTree::queryIntersectionCandidates(
    const std::map<int,LINALG::Matrix<3,2> >&   currentXAABBs, 
    const std::vector< LINALG::Matrix<3,2> >&   structure_AABBs,
    const LINALG::Matrix<3,2>&                  eleXAABB,
    std::set<int>&                              elementset) 
{
  TEUCHOS_FUNC_TIME_MONITOR("SearchTree - queryTime");

  elementset.clear();
  
  if(treeRoot_ == Teuchos::null)
      dserror("tree is not yet initialized !!!");
  
  if(!intersectionOfXAABB<3>(treeRoot_->getNodeBox(), eleXAABB))
    return;

  bool intersection = false;
  for(unsigned int i = 0; i < structure_AABBs.size(); i++ )
    if(intersectionOfXAABB<3>(structure_AABBs[i], eleXAABB))
      intersection = true;
  
  if(!intersection)
    return;
  
  if(!(treeRoot_->getElementList().empty()))
    treeRoot_->queryIntersectionCandidates(currentXAABBs, eleXAABB, elementset);

  return;
}



/*----------------------------------------------------------------------*
 | returns contact elements for a given element              u.may 02/09|
 | for multibody contact                                                |                
 *----------------------------------------------------------------------*/
void GEO::SearchTree::searchMultibodyContactElements(
    const std::map<int,LINALG::Matrix<9,2> >&   currentKDOPs, 
    const LINALG::Matrix<9,2>&                  queryKDOP,
    const int                                   label,
    std::set<int>&                              contactEleIds)
{
  TEUCHOS_FUNC_TIME_MONITOR("SearchTree - queryTime");

  if(treeRoot_ == Teuchos::null)
    dserror("tree is not yet initialized !!!");

  if(!treeRoot_->getElementList().empty())
    treeRoot_->searchMultibodyContactElements(currentKDOPs, queryKDOP, label, contactEleIds);
  else
    dserror("element list is empty");

  return;
}



/*----------------------------------------------------------------------*
 | print tree to gmsh file                                 peder   07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::printTree(
    const string  prefix, 
    const int     step) const
{
  cout << endl << "writing... ";
  if (treeRoot_ == Teuchos::null) {
    cout << "nothing to write, tree not initialized yet -> done" << endl;
    return;
  }
  if (treeRoot_->getElementList().empty()){ 
    cout << "nothing to write, tree empty -> done" << endl;
    return;
  }
  
  std::stringstream filename;
  std::stringstream node_string;
  filename << prefix << "_octtree" << std::setw(5) << setfill('0') << step << ".pos";
  cout << " " << filename.str() << " ...";
  node_string << "View \" " << "XFEM_FSI_Octree \" {" << endl;  
  treeRoot_->printTreeNode(max_depth_,node_string);
  node_string << "};" << endl;
  std::ofstream f_system(filename.str().c_str());
  f_system << node_string.str();
  f_system.close();
  cout << " done" << endl;
}



/*----------------------------------------------------------------------*
 | evaluate tree metrics (METRICS)                         peder   07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::evaluateTreeMetrics(
    const int     step) const
{
  cout << "\t********************* TREE METRICS ******************" << endl;

  if (treeRoot_ == Teuchos::null)
  { 
    cout << "tree root was not initialized yet, nothing to print" << endl;
    return;
  }
  
  cout.precision(5);
  cout << "\tnumber tree nodes        : " << treeRoot_->getNumNodesInTree() << endl;
  cout << "\ttree depth               : " << max_depth_-treeRoot_->getDepth() << " (max: "<< max_depth_<<")"<< endl;
  cout << "\t***************************************************" << endl;
  return;
}



/*----------------------------------------------------------------------*
 | c-tor TreeNode                                            u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::SearchTree::TreeNode::TreeNode(
    const TreeNode* const                       parent,
    const int                                   depth, 
    const LINALG::Matrix<3,2> &                 nodeBox,
    const TreeType                              treeType):
    parent_(parent),
    treedepth_(depth),
    treeNodeType_(LEAF_NODE),
    treeType_(treeType),
    label_(-1),
    nodeBox_(nodeBox),
    xPlaneCoordinate_( (nodeBox_(0,0)+0.5*(nodeBox_(0,1)-nodeBox_(0,0))) ),
    yPlaneCoordinate_( (nodeBox_(1,0)+0.5*(nodeBox_(1,1)-nodeBox_(1,0))) ),
    zPlaneCoordinate_( (nodeBox_(2,0)+0.5*(nodeBox_(2,1)-nodeBox_(2,0))) )
{    
     
      children_.assign(getNumChildren(), Teuchos::null);

}
  
    

/*====================================================================*/
  /* octtree */
  /*--------------------------------------------------------------------*/
  /* numbering of children
   * child 0: 
   * child 1: 
   * child 2: 
   * child 3: 
   * child 4:
   * child 5: 
   * child 6: 
   * child 7: 
   * 
   *                      z
   *                      |           
   *             4========|================6¨
   *           //|        |               /||
   *          // |        |              //||
   *         //  |        |             // ||
   *        //   |        |            //  ||
   *       //    |        |           //   ||
   *      //     |        |          //    ||
   *     //      |        |         //     ||
   *     5=========================7       ||
   *    ||       |        |        ||      ||
   *    ||       |        o--------||---------y
   *    ||       |       /         ||      ||
   *    ||       0------/----------||------2
   *    ||      /      /           ||     //
   *    ||     /      /            ||    //
   *    ||    /      /             ||   //
   *    ||   /      /              ||  //
   *    ||  /      /               || //
   *    || /      x                ||//
   *    ||/                        ||/
   *     1=========================3
   *
   */
  /*====================================================================*/


/*  int index = 0;
    if (point(0) > xPlaneCoordinate_)
      index += 1;
    if (point(1) > yPlaneCoordinate_)
      octIdx += 2;
    if (point(2) > zPlaneCoordinate_)
      index += 4;
*/

    
    

/*----------------------------------------------------------------------*
 | d-tor TreeNode                                            u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::SearchTree::TreeNode::~TreeNode() {}



/*----------------------------------------------------------------------*
 | clears node and deletes chlidren                          u.may 08/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::clear()
{
  treeNodeType_ = LEAF_NODE;
  label_ = -1; 

  children_.clear();
  children_.assign(getNumChildren(),Teuchos::null);
}



/*----------------------------------------------------------------------*
 | sets element list                                         u.may 07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::setElementList(
   const std::map<int, std::set<int> >& elementsByLabel)
{
   elementList_ = elementsByLabel;
}



/*----------------------------------------------------------------------*
 | set label                                                 peder 07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::setLabel(
	const int label)
{
  label_ = label;
}



/*----------------------------------------------------------------------*
 | set nearest object                                        u.may 09/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::setNearestObject(
      const GEO::NearestObject&   nearestObject)
{
  nearestObject_ = nearestObject;
}



/*----------------------------------------------------------------------*
 | get center of treenode                                  peder   07/08|
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,1> GEO::SearchTree::TreeNode::getCenterCoord() const
{
  LINALG::Matrix<3,1> centerCoord;

  centerCoord(0) = this->xPlaneCoordinate_;
  centerCoord(1) = this->yPlaneCoordinate_;
  centerCoord(2) = this->zPlaneCoordinate_;
  return centerCoord;
}



/*----------------------------------------------------------------------*
 | get number of children                                  u.may   08/08|
 *----------------------------------------------------------------------*/
int GEO::SearchTree::TreeNode::getNumChildren(
    ) const
{
  if(treeType_ == OCTTREE)
    return 8;
  else if(treeType_ == QUADTREE)
    return 4;
  else
    dserror("treetype does not exist");
  
  return -1;
}



/*----------------------------------------------------------------------*
 | get child                                               peder   07/08|
 *----------------------------------------------------------------------*/
const Teuchos::RCP< GEO::SearchTree::TreeNode > GEO::SearchTree::TreeNode::getChild(
    const int index) const
{
  return children_[index];
}



/*----------------------------------------------------------------------*
 | get node box of a child specified by index              peder   07/08|
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,2> GEO::SearchTree::TreeNode::getChildNodeBox(
	const int index) const
{
  LINALG::Matrix<3,2> childNodeBox;
  
  // determine x-coordinates
  if ((index==1) || (index==3) || (index==5) || (index==7))
  {
    childNodeBox(0,0) = xPlaneCoordinate_;
    childNodeBox(0,1) = nodeBox_(0,1);
  }
  else 
  {
    childNodeBox(0,0) = nodeBox_(0,0);
    childNodeBox(0,1) = xPlaneCoordinate_;
  }  
  
  // determine y-coordinates
  if ((index==2) || (index==3) || (index==6) || (index==7))
  {
    childNodeBox(1,0) = yPlaneCoordinate_;
    childNodeBox(1,1) = nodeBox_(1,1);
  }
  else 
  {
    childNodeBox(1,0) = nodeBox_(1,0);
    childNodeBox(1,1) = yPlaneCoordinate_;
  }
  
  // determine z coordinates
  if (index>3)
  {
    childNodeBox(2,0) = zPlaneCoordinate_;
    childNodeBox(2,1) = nodeBox_(2,1);
  }
  else 
  {
    if(treeType_ == OCTTREE)
    {
      childNodeBox(2,0) = nodeBox_(2,0);
      childNodeBox(2,1) = zPlaneCoordinate_;
    }
    else
    {
      childNodeBox(2,0) = 0.0;
      childNodeBox(2,1) = 0.0;
    }
  }
  //  printf("created chldAABB(%f\t%f\t%f\t%f\t%f\t%f)\n", childNodeBox(0,0),childNodeBox(0,1),childNodeBox(1,0),childNodeBox(1,1),childNodeBox(2,0),childNodeBox(2,1));
  return childNodeBox;
}



/*----------------------------------------------------------------------*
 | get node box of a child specified by index              peder   07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::getChildNodeBox(
  const int             index,
  LINALG::Matrix<3,2>&  childNodeBox) const
{
  childNodeBox.Clear();
  
  // determine x-coordinates
  if ((index==1) || (index==3) || (index==5) || (index==7))
  {
    childNodeBox(0,0) = xPlaneCoordinate_;
    childNodeBox(0,1) = nodeBox_(0,1);
  }
  else 
  {
    childNodeBox(0,0) = nodeBox_(0,0);
    childNodeBox(0,1) = xPlaneCoordinate_;
  }  
  
  // determine y-coordinates
  if ((index==2) || (index==3) || (index==6) || (index==7))
  {
    childNodeBox(1,0) = yPlaneCoordinate_;
    childNodeBox(1,1) = nodeBox_(1,1);
  }
  else 
  {
    childNodeBox(1,0) = nodeBox_(1,0);
    childNodeBox(1,1) = yPlaneCoordinate_;
  }
  
  // determine z coordinates
  if (index>3)
  {
    childNodeBox(2,0) = zPlaneCoordinate_;
    childNodeBox(2,1) = nodeBox_(2,1);
  }
  else 
  {
    if(treeType_ == OCTTREE)
    {
      childNodeBox(2,0) = nodeBox_(2,0);
      childNodeBox(2,1) = zPlaneCoordinate_;
    }
    else
    {
      childNodeBox(2,0) = 0.0;
      childNodeBox(2,1) = 0.0;
    }
  }
  //  printf("created chldAABB(%f\t%f\t%f\t%f\t%f\t%f)\n", childNodeBox(0,0),childNodeBox(0,1),childNodeBox(1,0),childNodeBox(1,1),childNodeBox(2,0),childNodeBox(2,1));
  return;
}





/*----------------------------------------------------------------------*
 | insert element in tree node                              u.may  07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::insertElement(
    const int   label,
    const int   eleId) 
{
    elementList_[label].insert(eleId);
}



/*----------------------------------------------------------------------*
 | create children                                         u.may   08/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::createChildren(
    const DRT::Discretization&      	             dis,
    const std::map<int,LINALG::Matrix<3,1> >& 	   currentpositions)
{
  // create empty children
  for(int index = 0; index < getNumChildren(); index++)
    children_[index] = rcp(new TreeNode(this, (treedepth_-1), getChildNodeBox(index), treeType_));
   
  // insert elements into child node
  for (std::map<int, std::set<int> >::const_iterator labelIter = elementList_.begin(); labelIter != elementList_.end(); labelIter++)
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {
      std::vector<int> elementClassification = classifyElement(dis.gElement(*eleIter),currentpositions);
      for(unsigned int count = 0; count < elementClassification.size(); count++)
      {
        children_[elementClassification[count]]->insertElement(labelIter->first,*eleIter);
      }
    }

  // this node becomes an inner tree node
  treeNodeType_ = INNER_NODE;
}

/*----------------------------------------------------------------------*
 | create children for a search tree; this method was written especially|
 | for search trees which are concerned with searching points close to  |
 | other points                                            cyron   04/09|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::createChildren(
    const std::map<int,LINALG::Matrix<3,1> >&      currentpositions)
{
  /*create getNumChildren (for Octree e.g. 8) empty children; each child is a tree node itself; the child
   *node gets its volume (node box) at once; however, the elements situated in this node box are not yet 
   *assigned to the child; this happens only in the next step*/
  for(int index = 0; index < getNumChildren(); index++)
    children_[index] = rcp(new TreeNode(this, (treedepth_-1), getChildNodeBox(index), treeType_));
   
  /*insert tree elements into child node; we first loop throug all different lables in the variable elementList_;
   * in case that just points close to other points are looked for there is typically just one label so that the outer
   * loop is run through one time only; for each label the inner loop runs through all the elements in the set assigend
   * to this label*/
  
  
  for (std::map<int, std::set<int> >::const_iterator labelIter = elementList_.begin(); labelIter != elementList_.end(); labelIter++)
    for (std::set<int>::const_iterator pointIter = (labelIter->second).begin(); pointIter != (labelIter->second).end(); pointIter++)
    {
      int pointClassification = classifyPoint((currentpositions.find(*pointIter))->second);
      children_[pointClassification]->insertElement(labelIter->first,*pointIter);
    }

  // this node becomes an inner tree node, because now it has at least one child
  treeNodeType_ = INNER_NODE;
}



/*----------------------------------------------------------------------*
 | create children                                         u.may   08/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::createChildren(
    const std::map<int,LINALG::Matrix<3,2> >&      currentXAABBs)
{
  // create empty children
  static LINALG::Matrix<3,2> childNodeBox;
  for(int index = 0; index < getNumChildren(); index++)
  { 
    getChildNodeBox(index, childNodeBox );
    children_[index] = rcp(new TreeNode(this, (treedepth_-1), childNodeBox, treeType_));
  }
  
  static std::vector<int> elementClassification;
  // insert elements into child node
  for (std::map<int, std::set<int> >::const_iterator labelIter = elementList_.begin(); labelIter != elementList_.end(); labelIter++)
  {
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {    
      classifyXAABB(currentXAABBs.find(*eleIter)->second, elementClassification);
      for(unsigned int count = 0; count < elementClassification.size(); count++)
        children_[elementClassification[count]]->insertElement(labelIter->first,*eleIter);
    }
  }
  // this node becomes an inner tree node
  treeNodeType_ = INNER_NODE;
}



/*----------------------------------------------------------------------*
 | create children using KDOPS                             u.may   02/09|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::createChildren(
    const std::map<int,LINALG::Matrix<9,2> >&      currentKDOPs)
{
  // create empty children
  static LINALG::Matrix<3,2> childNodeBox;
  for(int index = 0; index < getNumChildren(); index++)
  { 
    getChildNodeBox(index, childNodeBox );
    children_[index] = rcp(new TreeNode(this, (treedepth_-1), childNodeBox, treeType_));
  }
  
  static std::vector<int> elementClassification;
  // insert elements into child node
  for (std::map<int, std::set<int> >::const_iterator labelIter = elementList_.begin(); labelIter != elementList_.end(); labelIter++)
  {
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {    
      classifyKDOP(currentKDOPs.find(*eleIter)->second, elementClassification);
      for(unsigned int count = 0; count < elementClassification.size(); count++)
        children_[elementClassification[count]]->insertElement(labelIter->first,*eleIter);
    }
  }
  // this node becomes an inner tree node
  treeNodeType_ = INNER_NODE;
}
  


/*----------------------------------------------------------------------*
 | set XFEM label of empty children                        u.may   08/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::setXFEMLabelOfEmptyChildren(
    const DRT::Discretization&                    dis,
    const std::map<int,LINALG::Matrix<3,1> >&     currentpositions)
{    
  for(int index = 0; index < getNumChildren(); index++)
  {
    // if one of the created children is empty, set label immediately
    if ((children_[index]->getElementList()).empty())
    {
      const LINALG::Matrix<3,1>  childNodeCenter(children_[index]->getCenterCoord());
      // xfem label has to be computed on this level because child is empty
      children_[index]->setLabel(getXFEMLabel(dis, currentpositions, childNodeCenter, elementList_));
    }
  }
}



/*----------------------------------------------------------------------*
 | set XFEM label of empty children                        u.may   02/09|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::setXFEMLabelOfEmptyChildren(
    const vector<vector<int> >&                   triangleList,
    const std::vector<GEO::InterfacePoint>&       pointList)
{    
  for(int index = 0; index < getNumChildren(); index++)
  {
    // if one of the created children is empty, set label immediately
    if ((children_[index]->getElementList()).empty())
    {
      const LINALG::Matrix<3,1>  childNodeCenter(children_[index]->getCenterCoord());
      // xfem label has to be computed on this level because child is empty
      children_[index]->setLabel(getXFEMLabel(triangleList, pointList, childNodeCenter, elementList_));
    }
  }
}



/*----------------------------------------------------------------------*
 | set XFEM label and nearest object of empty children     u.may   08/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::setXFEMLabelAndNearestObjectOfEmptyChildren(
    const DRT::Discretization&                    dis,
    const std::map<int,LINALG::Matrix<3,1> >&     currentpositions)
{    
  for(int index = 0; index < getNumChildren(); index++)
  {
    // if one of the created children is empty, set label immediately
    if ((children_[index]->getElementList()).empty())
    {
      const LINALG::Matrix<3,1> childNodeCenter(children_[index]->getCenterCoord());
      // xfem label has to be computed on this level because child is empty
      GEO::NearestObject nearestObject;
      int label = getXFEMLabelAndNearestObject(dis, currentpositions, childNodeCenter, elementList_, nearestObject);
      children_[index]->setLabel(label);
      children_[index]->setNearestObject(nearestObject);
    }
  }
}



/*----------------------------------------------------------------------*
 | classifies point, i.e. returns for each point the index of the child |
 | node in whose node box it is situated                   peder   07/08|
 *----------------------------------------------------------------------*/
inline int GEO::SearchTree::TreeNode::classifyPoint(
    const LINALG::Matrix<3,1>&   point) const
{
  
  int childIndex = 0;
  if (point(0) > xPlaneCoordinate_)
    childIndex += 1;
  if (point(1) > yPlaneCoordinate_)
    childIndex += 2;
  if(treeType_ == OCTTREE)
    if (point(2) > zPlaneCoordinate_)
      childIndex += 4;
  
  return childIndex;
}



/*----------------------------------------------------------------------*
 | classifiy AABB in node                                  u.may   07/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::TreeNode::classifyXAABB(
    const LINALG::Matrix<3,2>&   AABB) const 
{
  // collect all children the XAABB is lying in
  // use tolerances such that it is ensured that no child is left behind :-) 
  // XAABB s which are lying between plane-tol and plane + tol are collected in each of the two chlidren
  
  std::vector<int> octants;
  
  // check max_x greater than x-plane
  if (AABB(0, 1) > (xPlaneCoordinate_ - GEO::TOL7) ) 
  {
    // check max_y greater than y-plane
    if (AABB(1, 1) > (yPlaneCoordinate_ - GEO::TOL7) ) 
    {
      if(treeType_ == OCTTREE)
      {
        // check max_z greater than z-plane
        if (AABB(2, 1) > (zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(7);

        // check min_z less than z-plane
        if (AABB(2, 0) < (zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(3);
      }
      else if (treeType_ == QUADTREE)
        octants.push_back(3);
    }
    
    // check min_y less than y-plane
    if (AABB(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
      if(treeType_ == OCTTREE)
      {
        // check max_z greater than z-plane
        if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(5);
              
        // check min_z less than z-plane
        if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(1);
      }
      else if(treeType_ == QUADTREE)
        octants.push_back(1);
    }
  }

  // check min_x less than x-plane
  if (AABB(0, 0) < ( xPlaneCoordinate_ + GEO::TOL7) ) 
  {
    // check min_y less than y-plane
    if (AABB(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
      
      if(treeType_ == OCTTREE)
      {
        // check min_z less than z-plane
        if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(0);
      
        // check max_z greater than z-plane
        if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(4);
      }
      
      else if(treeType_ == QUADTREE)
        octants.push_back(0);
    }
    
    // check max_y greater than y-plane
    if (AABB(1, 1) > ( yPlaneCoordinate_ - GEO::TOL7) ) 
    {       
      if(treeType_ == OCTTREE)
      {
        // check max_z greater than z-plane
        if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(6);
      
        // check min_z less than z-plane
        if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(2);
      }
      
      else if(treeType_ == QUADTREE)
        octants.push_back(2);
    }
  }
  return octants;
}



/*----------------------------------------------------------------------*
 | classifiy AABB in node                                  u.may   07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::classifyXAABB(
    const LINALG::Matrix<3,2>&    AABB,
    std::vector<int>&             octants) const 
{
  
  // collect all children the XAABB is lying in
  // use tolerances such that it is ensured that no child is left behind :-) 
  // XAABB s which are lying between plane-tol and plane + tol are collected in each of the two chlidren
  
  octants.clear();
  octants.reserve(8);
  
  // check max_x greater than x-plane
  if (AABB(0, 1) > (xPlaneCoordinate_ - GEO::TOL7) ) 
  {
    // check max_y greater than y-plane
    if (AABB(1, 1) > (yPlaneCoordinate_ - GEO::TOL7) ) 
    {
        // check max_z greater than z-plane
        if (AABB(2, 1) > (zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(7);

        // check min_z less than z-plane
        if (AABB(2, 0) < (zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(3);

    }
    
    // check min_y less than y-plane
    if (AABB(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
        // check max_z greater than z-plane
        if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(5);
              
        // check min_z less than z-plane
        if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(1);
    }
  }

  // check min_x less than x-plane
  if (AABB(0, 0) < ( xPlaneCoordinate_ + GEO::TOL7) ) 
  {
    // check min_y less than y-plane
    if (AABB(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
        // check min_z less than z-plane
        if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(0);
      
        // check max_z greater than z-plane
        if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(4);

    }
    
    // check max_y greater than y-plane
    if (AABB(1, 1) > ( yPlaneCoordinate_ - GEO::TOL7) ) 
    {       
        // check max_z greater than z-plane
        if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(6);
      
        // check min_z less than z-plane
        if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(2);

    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | classifiy AABB in node                                  u.may   02/09|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::classifyKDOP(
    const LINALG::Matrix<9,2>&    KDOP,
    std::vector<int>&             octants) const 
{
  
  // collect all children the XAABB is lying in
  // use tolerances such that it is ensured that no child is left behind :-) 
  // XAABB s which are lying between plane-tol and plane + tol are collected in each of the two chlidren
  
  octants.clear();
  octants.reserve(8);
  
  // check max_x greater than x-plane
  if (KDOP(0, 1) > (xPlaneCoordinate_ - GEO::TOL7) ) 
  {
    // check max_y greater than y-plane
    if (KDOP(1, 1) > (yPlaneCoordinate_ - GEO::TOL7) ) 
    {
        // check max_z greater than z-plane
        if (KDOP(2, 1) > (zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(7);

        // check min_z less than z-plane
        if (KDOP(2, 0) < (zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(3);

    }
    
    // check min_y less than y-plane
    if (KDOP(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
        // check max_z greater than z-plane
        if (KDOP(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(5);
              
        // check min_z less than z-plane
        if (KDOP(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(1);
    }
  }

  // check min_x less than x-plane
  if (KDOP(0, 0) < ( xPlaneCoordinate_ + GEO::TOL7) ) 
  {
    // check min_y less than y-plane
    if (KDOP(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
        // check min_z less than z-plane
        if (KDOP(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(0);
      
        // check max_z greater than z-plane
        if (KDOP(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(4);

    }
    
    // check max_y greater than y-plane
    if (KDOP(1, 1) > ( yPlaneCoordinate_ - GEO::TOL7) ) 
    {       
        // check max_z greater than z-plane
        if (KDOP(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(6);
      
        // check min_z less than z-plane
        if (KDOP(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(2);

    }
  }
  return;
}




/*----------------------------------------------------------------------*
 | classifiy AABB in node                                  u.may   07/08|
 *----------------------------------------------------------------------*/
bool GEO::SearchTree::TreeNode::classifyXAABB(
    int&                         index,
    const LINALG::Matrix<3,2>&   AABB) const 
{
  
  // collect all children the XAABB is lying in
  // use tolerances such that it is ensured that no child is left behind :-) 
  // XAABB s which are lying between plane-tol and plane + tol are collected in each of the two chlidren
  
  bool oneIndex = true;
  index = -1;
  
  
  // check max_x greater than x-plane
  if (AABB(0, 1) > (xPlaneCoordinate_ - GEO::TOL7) ) 
  {
    // check max_y greater than y-plane
    if (AABB(1, 1) > (yPlaneCoordinate_ - GEO::TOL7) ) 
    {
      if(treeType_ == OCTTREE)
      {
        // check max_z greater than z-plane
        if (AABB(2, 1) > (zPlaneCoordinate_ - GEO::TOL7) )
        {
          if(index == -1)
            index = 7;
          else
          {
            index = -1;
            return false;
          }
        }
        // check min_z less than z-plane
        if (AABB(2, 0) < (zPlaneCoordinate_ + GEO::TOL7) )
        {
          if(index == -1)
            index = 3;
          else
          {
            index = -1;
            return false;
          }
        }
      }
      else if (treeType_ == QUADTREE)
        if(index == -1)
          index = 3;
        else
        {
          index = -1;
          return false;
        }
    }
    
    // check min_y less than y-plane
    if (AABB(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
      if(treeType_ == OCTTREE)
      {
        // check max_z greater than z-plane
        if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
        {
          if(index == -1)
            index = 5;
          else
          {
            index = -1;
            return false;
          }
        }
              
        // check min_z less than z-plane
        if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
        {
          if(index == -1)
            index = 1;
          else
          {
            index = -1;
            return false;
          }
        }
      }
      else if(treeType_ == QUADTREE)
        if(index == -1)
          index = 1;
        else
        {
          index = -1;
          return false;
        }
    }
  }

  // check min_x less than x-plane
  if (AABB(0, 0) < ( xPlaneCoordinate_ + GEO::TOL7) ) 
  {
    // check min_y less than y-plane
    if (AABB(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
      
      if(treeType_ == OCTTREE)
      {
        // check min_z less than z-plane
        if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
        {
          if(index == -1)
            index = 0;
          else
          {
            index = -1;
            return false;
          }
        }
        // check max_z greater than z-plane
        if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
        {
          if(index == -1)
            index = 4;
          else
          {
            index = -1;
            return false;
          }
        }
      }
      else if(treeType_ == QUADTREE)
        if(index == -1)
          index = 0;
        else
        {
          index = -1;
          return false;
        }
    }
    
    // check max_y greater than y-plane
    if (AABB(1, 1) > ( yPlaneCoordinate_ - GEO::TOL7) ) 
    {       
      if(treeType_ == OCTTREE)
      {
        // check max_z greater than z-plane
        if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
        {
          if(index == -1)
            index = 6;
          else
          {
            index = -1;
            return false;
          }
        }
      
        // check min_z less than z-plane
        if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
        {
          if(index == -1)
            index = 2;
          else
          {
            index = -1;
            return false;
          }
        }
      }
      else if(treeType_ == QUADTREE)
        if(index == -1)
          index = 2;
        else
        {
          index = -1;
          return false;
        }
    }
  }
  return oneIndex;
}



/*----------------------------------------------------------------------*
 | classifiy KDOP in node                                  u.may   02/09|
 *----------------------------------------------------------------------*/
bool GEO::SearchTree::TreeNode::classifyKDOP(
    int&                         index,
    const LINALG::Matrix<9,2>&   KDOP) const 
{
  // collect all children the XAABB is lying in
  // use tolerances such that it is ensured that no child is left behind :-) 
  // XAABB s which are lying between plane-tol and plane + tol are collected in each of the two chlidren
  bool oneIndex = true;
  index = -1;
  
  // check max_x greater than x-plane
  if (KDOP(0, 1) > (xPlaneCoordinate_ - GEO::TOL7) ) 
  {
    // check max_y greater than y-plane
    if (KDOP(1, 1) > (yPlaneCoordinate_ - GEO::TOL7) ) 
    {
      if(treeType_ == OCTTREE)
      {
        // check max_z greater than z-plane
        if (KDOP(2, 1) > (zPlaneCoordinate_ - GEO::TOL7) )
        {
          if(index == -1)
            index = 7;
          else
          {
            index = -1;
            return false;
          }
        }
        // check min_z less than z-plane
        if (KDOP(2, 0) < (zPlaneCoordinate_ + GEO::TOL7) )
        {
          if(index == -1)
            index = 3;
          else
          {
            index = -1;
            return false;
          }
        }
      }
      else if (treeType_ == QUADTREE)
        if(index == -1)
          index = 3;
        else
        {
          index = -1;
          return false;
        }
    }
    
    // check min_y less than y-plane
    if (KDOP(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
      if(treeType_ == OCTTREE)
      {
        // check max_z greater than z-plane
        if (KDOP(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
        {
          if(index == -1)
            index = 5;
          else
          {
            index = -1;
            return false;
          }
        }
              
        // check min_z less than z-plane
        if (KDOP(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
        {
          if(index == -1)
            index = 1;
          else
          {
            index = -1;
            return false;
          }
        }
      }
      else if(treeType_ == QUADTREE)
        if(index == -1)
          index = 1;
        else
        {
          index = -1;
          return false;
        }
    }
  }

  // check min_x less than x-plane
  if (KDOP(0, 0) < ( xPlaneCoordinate_ + GEO::TOL7) ) 
  {
    // check min_y less than y-plane
    if (KDOP(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
      
      if(treeType_ == OCTTREE)
      {
        // check min_z less than z-plane
        if (KDOP(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
        {
          if(index == -1)
            index = 0;
          else
          {
            index = -1;
            return false;
          }
        }
        // check max_z greater than z-plane
        if (KDOP(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
        {
          if(index == -1)
            index = 4;
          else
          {
            index = -1;
            return false;
          }
        }
      }
      else if(treeType_ == QUADTREE)
        if(index == -1)
          index = 0;
        else
        {
          index = -1;
          return false;
        }
    }
    
    // check max_y greater than y-plane
    if (KDOP(1, 1) > ( yPlaneCoordinate_ - GEO::TOL7) ) 
    {       
      if(treeType_ == OCTTREE)
      {
        // check max_z greater than z-plane
        if (KDOP(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
        {
          if(index == -1)
            index = 6;
          else
          {
            index = -1;
            return false;
          }
        }
      
        // check min_z less than z-plane
        if (KDOP(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
        {
          if(index == -1)
            index = 2;
          else
          {
            index = -1;
            return false;
          }
        }
      }
      else if(treeType_ == QUADTREE)
        if(index == -1)
          index = 2;
        else
        {
          index = -1;
          return false;
        }
    }
  }
  return oneIndex;
}



/*----------------------------------------------------------------------*
 | classifiy element in node                               peder   07/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::TreeNode::classifyElement(
    const DRT::Element*       		              element,
    const std::map<int,LINALG::Matrix<3,1> >& 	currentpositions
    ) const
{
  const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(element,currentpositions));
  GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
  GEO::checkRoughGeoType(element, xyze, eleGeoType);
  const LINALG::Matrix<3,2> elemXAABB(GEO::computeFastXAABB(element->Shape(), xyze, eleGeoType));  
  return classifyXAABB(elemXAABB);
}



/*----------------------------------------------------------------------*
 | classifiy element in tree node                          u.may   07/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::TreeNode::classifyElement(
        const DRT::Element*               element,
        const LINALG::SerialDenseMatrix&  xyze_element    
        ) const 
{
  GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
  GEO::checkRoughGeoType(element, xyze_element , eleGeoType);
  const LINALG::Matrix<3,2> elemXAABB(GEO::computeFastXAABB(element->Shape(), xyze_element , eleGeoType));  
  return classifyXAABB(elemXAABB);
}



/*----------------------------------------------------------------------*
 | classifiy element in node                               u.may   08/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::TreeNode::classifyRadius(
    const double                radius,
    const LINALG::Matrix<3,1>&	point
    ) const
{
  
  /*coordinates of axis-aligned boundary box around point, which stretechs over 
   * length radius in each coordinate direction*/
  LINALG::Matrix<3,2> radiusXAABB;  

  for(int dim = 0; dim < 3; dim++)
  {
    radiusXAABB(dim,0) = ( point(dim) - radius) - GEO::TOL7;
    radiusXAABB(dim,1) = ( point(dim) + radius) + GEO::TOL7;
  }
  
  if(treeType_ == QUADTREE)
  {
    radiusXAABB(2,0) = 0.0;
    radiusXAABB(2,1) = 0.0;
  }
   
  //return indices of childs which overlap with given axis-aligned boundary box of edge length 2*radius around point
  return classifyXAABB(radiusXAABB);
}



/*----------------------------------------------------------------------*
 | update tree node                                        u.may   08/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::updateTreeNode(
   const LINALG::Matrix<3,2>& AABB_old, 
   const LINALG::Matrix<3,2>& AABB_new)
{
  for(int i = 0; i < getNumChildren(); i++)
  {
    if(children_[i] != Teuchos::null)
    {
      if( GEO::inSameNodeBox(AABB_old, AABB_new, children_[i]->getNodeBox()) )
        children_[i]->updateTreeNode(AABB_old, AABB_new);
      else
        children_[i]->clear();
    }
  }
} 



/*----------------------------------------------------------------------*
 | return xfem label for point (interface method)          u.may   07/08|
 | and nearest object                                                   |
 *----------------------------------------------------------------------*/
int GEO::SearchTree::TreeNode::queryFSINearestObject(
    const DRT::Discretization&                      dis,
    const std::map<int,LINALG::Matrix<3,1> >&       currentpositions, 
    const std::map<int,LINALG::Matrix<3,2> >&       currentXAABBs,
    const LINALG::Matrix<3,1>&                      point,
    GEO::NearestObject&                             nearestObject
    ) 
{
  
  switch (treeNodeType_) 
  {
    case INNER_NODE:
    {       
      return children_[classifyPoint(point)]->queryFSINearestObject(dis, currentpositions, currentXAABBs, point, nearestObject);
      break;
    }
    case LEAF_NODE:   
    {
      if (elementList_.empty())
      {
        nearestObject = nearestObject_;
        return label_;
      }
      
      // max depth reached, counts reverse
      if (treedepth_ <= 0 || (elementList_.size()==1 && (elementList_.begin()->second).size() == 1) )
      {
        // nearest object refers only to the nearest object found in this particular tree node
        int xfemLabel = GEO::getXFEMLabelAndNearestObject(dis, currentpositions, point, elementList_, nearestObject);
        
        const TreeNode* workingNode = this;
        while(!GEO::pointInMinCircleInTreeNode(nearestObject.getPhysCoord(), point, workingNode->getNodeBox(), (!workingNode->hasParent()) ))
        {
          if(!workingNode->hasParent())
            dserror("this treenode has no parent");
          
          workingNode = workingNode->getParent();
          xfemLabel = GEO::getXFEMLabelAndNearestObject(dis, currentpositions, point, workingNode->getElementList(), nearestObject);
        }
        return xfemLabel;
      }

      // dynamically grow tree otherwise, create children and set label for empty children
      createChildren(currentXAABBs);
      setXFEMLabelOfEmptyChildren(dis, currentpositions);
      // search in apropriate child node
      return children_[classifyPoint(point)]->queryFSINearestObject(dis, currentpositions, currentXAABBs, point, nearestObject);
      break;
    }
    default:
      dserror("should not get here\n");
  }
  return -1;
}



/*----------------------------------------------------------------------*
 | return xfem label for point (interface method)          u.may   07/08|
 *----------------------------------------------------------------------*/
int GEO::SearchTree::TreeNode::queryXFEMFSIPointType(
    const DRT::Discretization&                    dis,
    const std::map<int,LINALG::Matrix<3,1> >&     currentpositions, 
    const std::map<int,LINALG::Matrix<3,2> >&     currentXAABBs,
    const LINALG::Matrix<3,1>&                    point
    ) 
{ 
  switch (treeNodeType_) 
  {
    case INNER_NODE:
    {       
      return children_[classifyPoint(point)]->queryXFEMFSIPointType(dis, currentpositions, currentXAABBs, point);
      break;
    }
    case LEAF_NODE:   
    {
      if (elementList_.empty())
        return label_;
      
      // max depth reached, counts reverse
      if (treedepth_ <= 0 || (elementList_.size()==1 && (elementList_.begin()->second).size() == 1) )
      {
        // nearest object refers only to the nearest object found in this particular tree node
        GEO::NearestObject nearestObjectInNode;
        int xfemLabel = GEO::getXFEMLabelAndNearestObject(dis, currentpositions, point, elementList_, nearestObjectInNode);
        
        const TreeNode* workingNode = this;
        //cout << "point = " << nearestObjectInNode.getPhysCoord() <<  endl;
        //cout << "nodebox = " << workingNode->nodeBox_ <<  endl;
        
        while(!GEO::pointInTreeNode(nearestObjectInNode.getPhysCoord(), workingNode->nodeBox_))
        {
          //cout << "point = " << nearestObjectInNode.getPhysCoord() <<  endl;
          //cout << "nodebox = " << workingNode->nodeBox_ <<  endl;
          
          
          if(!workingNode->hasParent())
            dserror("this treenode has no parent");
          
          workingNode = workingNode->getParent();
          xfemLabel = GEO::getXFEMLabelAndNearestObject(dis, currentpositions, point, workingNode->getElementList(), nearestObjectInNode);
        }
        return xfemLabel;
      }
      
      // dynamically grow tree otherwise, create children and set label for empty children
      createChildren(currentXAABBs);
      setXFEMLabelOfEmptyChildren(dis, currentpositions);
      // search in apropriate child node
      return children_[classifyPoint(point)]->queryXFEMFSIPointType(dis, currentpositions,currentXAABBs, point);

      break;
    }
    default:
      dserror("should not get here\n");
  }
  return -1;
}



/*----------------------------------------------------------------------*
 | fix intersection for contact with xfem FSI              u.may   02/09|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::moveContactNodes(
    const vector<vector<int> >&                   triangleList,
    vector< GEO::InterfacePoint >&                pointList,
    const std::map<int,LINALG::Matrix<3,2> >&     triangleXAABBs,  
    const LINALG::Matrix<3,1>&                    querypoint,
    const int                                     querypointId)
{ 
  switch (treeNodeType_) 
  {
    case INNER_NODE:
    {       
      return children_[classifyPoint(querypoint)]->moveContactNodes(triangleList, pointList, triangleXAABBs,  querypoint, querypointId);
      break; 
    }
    case LEAF_NODE:   
    {
      if(elementList_.empty())
      { 
        // point outside of a structure
        if(label_ == 0)
          return;
          
        // point inside a structure, step up and compute
        const TreeNode* workingNode = this;
        if(!workingNode->hasParent())
          dserror("this treenode has no parent");
          
        workingNode = workingNode->getParent();
        GEO::moveNodeOutOfStructure(triangleList, pointList, querypointId, elementList_);
        return;
      }
      // max depth reached, counts reverse
      if (treedepth_ <= 0 || (elementList_.size()==1 && (elementList_.begin()->second).size() == 1) )
      {
        // nearest object refers only to the nearest object found in this particular tree node
        GEO::NearestObject nearestObjectInNode;
        GEO::moveNodeOutOfStructure(triangleList, pointList, querypointId, elementList_);
    
        const TreeNode* workingNode = this;
        while(!GEO::pointInTreeNode(nearestObjectInNode.getPhysCoord(), workingNode->nodeBox_))
        {
          if(!workingNode->hasParent())
            dserror("this treenode has no parent");
          
          workingNode = workingNode->getParent();
          GEO::moveNodeOutOfStructure(triangleList, pointList, querypointId, elementList_);
        }
        return;
      }
      
      // dynamically grow tree otherwise, create children and set label for empty children
      createChildren(triangleXAABBs);
      setXFEMLabelOfEmptyChildren(triangleList, pointList);
      // search in apropriate child node
      children_[classifyPoint(querypoint)]->moveContactNodes(triangleList, pointList, triangleXAABBs,  querypoint, querypointId);
      return;
      break;
    }
    default:
      dserror("should not get here\n");
  }
  return;
}



/*----------------------------------------------------------------------*
 | returns nodes in the radius of a given point              u.may 08/08|
 *----------------------------------------------------------------------*/
std::map<int,std::set<int> > GEO::SearchTree::TreeNode::searchElementsInRadius(
    const DRT::Discretization& 	                  dis,
    const std::map<int,LINALG::Matrix<3,1> >&     currentpositions, 
    const LINALG::Matrix<3,1>&                    point,
    const double                                  radius, 
    const int                                     label) 
{  
  std::map<int,std::set<int> > eleMap;

  switch (treeNodeType_) 
  {
    case INNER_NODE:
    {       
      const vector<int> childindex = classifyRadius(radius, point);
      if(childindex.size() < 1)
        dserror("no child found\n");
      else if (childindex.size() ==1) // child node found which encloses AABB so step down
        return children_[childindex[0]]->searchElementsInRadius(dis, currentpositions, point, radius, label);
      else
        return GEO::getElementsInRadius(dis, currentpositions, point, radius, label, elementList_); 
      break;
    }
    case LEAF_NODE:   
    {
      if(elementList_.empty())
        return eleMap;

      // max depth reached, counts reverse
      if (treedepth_ <= 0 || (elementList_.size()==1 && (elementList_.begin()->second).size() == 1) )
        return GEO::getElementsInRadius(dis, currentpositions, point, radius, label, elementList_); 

      // dynamically grow tree otherwise, create children and set label for empty children
      // search in apropriate child node 
      const vector<int> childindex = classifyRadius(radius, point);
      if(childindex.size() < 1)
        dserror("no child found\n");
      else if (childindex.size() ==1) // child node found which encloses AABB so refine further
      { 
        createChildren(dis, currentpositions);
        return children_[childindex[0]]->searchElementsInRadius(dis, currentpositions, point, radius, label);
      }
      else // AABB does not fit into a single child node box anymore so don t refine any further
        return GEO::getElementsInRadius(dis, currentpositions, point, radius, label, elementList_); 
      break;
    }
    default:
      dserror("should not get here\n");
  }
  return eleMap;
}

/*-----------------------------------------------------------------------*
 |currentpositions assings a unique Id to a point, respectively. If the  |
 |points in currentpositions are nodes of a finite element discretization|
 |this unique Id is typically the GID of the finite elment nodes.        |
 |This method searches in the map currentpositions points at a distance  |
 |closer than radius to the querypoint and returns the Ids of these      |
 |points in a vector. Searching the tree is started by calling this      |
 |method for the tree itself; this method does not search itself. Rahter |
 |it starts a recursive search through all the nodes of the tree starting|
 |at the root node. The actual search is then carried out by the method  |
 |                                                            cyron 04/09|
 *-----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::searchPointsInRadius(
    const std::map<int,LINALG::Matrix<3,1> >&     currentpositions, 
    const LINALG::Matrix<3,1>&                    querypoint,
    const double                                  radius) 
{
  TEUCHOS_FUNC_TIME_MONITOR("SearchTree - queryTime");

  std::vector<int> nodes;

  /*treeRoot_ is the very node of the search tree which is no the highest level; 
   * chilren are nodes on a lower level than the root node; leaf nodes are nodes
   * on the lowest level; inner nodes are nodes neither on the highest nor lowest
   * level; stepping down the levels each node may split into several nodes on the
   * level below; these nodes below are called children of the node above;
   * 
   *                      0 (root node)
   *                    /   \
   *     (inner node)  0     0  (inner node & leaf node)
   *                 /  \   
   *    (leaf node) 0    0  (leaf node)
   *
   * each tree needs to have at least one node, which is the root node; each
   * node has certain elements; these elements are not necessarily finite
   * elements, but may be of any kind of entities, e.g. also nodes in a
   * finite element discretization; the elements of each node in the search
   * tree are administrated in the variable elementList_*/
  if(treeRoot_ == Teuchos::null)
    dserror("tree is not yet initialized !!!");

  if(!(treeRoot_->getElementList().empty()))
    nodes = treeRoot_->searchPointsInRadius(currentpositions,querypoint,radius);
  else
    dserror("tree element list is empty");

  return nodes;
}



/*-------------------------------------------------------------------------------*
 |returns vector of Ids of points closer radius to querypoint. This method       |
 |searches elements assigned to this tree node only; these elements represent a  |
 |subset of the points in currentpositions. In this subset the method searches   |
 |for points closer than radius to the querypoint.                    cyron 04/09|
 *------------------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::TreeNode::searchPointsInRadius(
    const std::map<int,LINALG::Matrix<3,1> >&     currentpositions, 
    const LINALG::Matrix<3,1>&                    querypoint,
    const double                                  radius) 
{  
  /*note: the tree is built dynamically upon query; first it needs to be initialized 
   * and has the root node only at this point; then once a query get in the root node
   * is identified as LEAF_NODE and in this method always the case LEAF_NODE applies;
   * in this case the tree is built up branch by branch as decribed in for the case
   * LEAF_NODE; at the same time for each node of the tree the variable elementList_
   * is defined; only as many new nodes are created as necessary to answer the very
   * query currently applied; as soon as new queries are applied for these queries the
   * already existing structure of children created for former queries may be applied;
   * then in this method the case INNER_NODE applies; however, it may also be that
   * the nodes of the tree created for former queries are not enough for the current
   * query; then again the case LEAF_NODE applies and the same procedure as throughout
   * the whole first query necessary is applied*/
  
  //vector for GIDs of nodes in the radius of a given querypoint
  std::vector<int> points;


  switch (treeNodeType_) 
  {
    /*if node is not at the lowest level of the tree; this case applies if current query has
     * not yet required any node structure of the tree which does not already exist from former
     * queries*/
    case INNER_NODE:
    {       
      /*get indices of all children nodes whose volume overlaps with an axis-aligned 
       * boundary box of edge length 2*radius around querypoint*/
      const vector<int> childindex = classifyRadius(radius,querypoint);
      
      //child node of tree found, which encloses AABB of querypoint completely; thus step down
      if (childindex.size() ==1) 
        return children_[childindex[0]]->searchPointsInRadius(currentpositions, querypoint, radius);
      
      /*no child overlaps AABB of querypoint; this has to be an error as stepping down 
       * the tree levels only happens if a child was found which comprises the AABB of
       * the querypoint completely; this child (i.e. the current tree node) has necessarily
       * one or more children itself which overlap with the querpoint's AABB*/     
      else if(childindex.size() < 1)
        dserror("no child found\n");
                
      /*several children found which overlap with AABB; so current tree node is the lowest one
       * which contains all the elements in the AABB of the querypoint; thus searching points
       * close to query point means searching among all the elements of the current tree node*/
      else
        return GEO::getPointsInRadius(currentpositions,querypoint,radius,elementList_); 
      break;
    }
    
    /*if node is at the lowest level of the tree; this means that a query requires a substructure of the
     * tree which does not yet exist; throughout the first query this case always applies*/   
    case LEAF_NODE:   
    {
      /*if volume of this tree node does not contain any tree elements then return points, which
       * is currently just an empty map*/ 
      if(elementList_.empty())
        return points;

      // max depth reached, counts reverse
      if (treedepth_ <= 0 || (elementList_.size()==1 && (elementList_.begin()->second).size() == 1) )
        return GEO::getPointsInRadius(currentpositions, querypoint, radius,elementList_); 

      /*in the following the tree grows dynamically first the vector childindex finds out with which
       * potential children of the current node the AABB of the querypoint would overlap; note that
       * these children do not exist yet. childindex just gives the indices of the childredn which
       * would matter if they already existed*/
      const vector<int> childindex = classifyRadius(radius, querypoint);
      
      /*if there is not even one potential child which would overlap with the AABB of the querypoint
       * there must be something wrong since stepping down to this node is only possible if at least
       * one of the elements of this node is relevant for the search; this element has to lie somewhere
       * in the volume of this node and this volume equals the union of all childrens' volumes*/
      if(childindex.size() < 1)
        dserror("no child found\n");
      
      /*there is only one potential child in whose volume there are elements relevant for this search;
       * create this child and step down to it to continue search within this child*/
      else if (childindex.size() ==1) // child node found which encloses AABB so refine further
      { 
        createChildren(currentpositions);
        return children_[childindex[0]]->searchPointsInRadius(currentpositions, querypoint, radius);
      }
      
      /*there are several potential children which are relevant for this search; thus there is no
       * need for further refinement of the tree; rather it's enough to search among the tree elements
       * of this very tree node*/ 
      else 
        return GEO::getPointsInRadius(currentpositions,querypoint,radius,elementList_); 
      break;
    }
    default:
      dserror("should not get here\n");
  }
  return points;
}



/*----------------------------------------------------------------------*
 | returns a set of elements whose XAABB s are               u.may 09/08|
 | intersecting with the XAABB of a given volume element                |
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::queryIntersectionCandidates(
    const std::map<int,LINALG::Matrix<3,2> >&   currentXAABBs, 
    const LINALG::Matrix<3,2>&                  eleXAABB,
    std::set<int>&                              elementset)
{
  switch(treeNodeType_) 
  {
    case INNER_NODE:
    {     
      int index = -1;
      if(classifyXAABB(index,eleXAABB))
      {
        children_[index]->queryIntersectionCandidates(currentXAABBs, eleXAABB, elementset);
        return;
      }
      else
      {
        GEO::getIntersectionCandidates(currentXAABBs, eleXAABB, elementList_, elementset); 
        return;
      }
      break;
    }
    case LEAF_NODE:   
    {
      if(elementList_.empty())
        return;

      // max depth reached, counts reverse
      if(treedepth_ <= 0 || (elementList_.begin()->second).size() == 1)
      {
        GEO::getIntersectionCandidates(currentXAABBs, eleXAABB, elementList_, elementset); 
        return;
      }
      // dynamically grow tree otherwise, create children and set label for empty children
      // search in apropriate child node  
      int index = -1;
      if(classifyXAABB(index,eleXAABB))
      {
        createChildren(currentXAABBs);
        children_[index]->queryIntersectionCandidates(currentXAABBs, eleXAABB, elementset);
        return;
      }
      else
      {
        GEO::getIntersectionCandidates(currentXAABBs, eleXAABB,  elementList_, elementset); 
        return;
      }
      break;
    }
    default:
      dserror("should not get here\n");
  }
  return;
}



/*----------------------------------------------------------------------*
 | returns contact elements based on kdops (CONTACT)         u.may 02/09|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::searchMultibodyContactElements(
    const std::map<int,LINALG::Matrix<9,2> >&   currentKDOPs, 
    const LINALG::Matrix<9,2>&                  queryKDOP,
    const int                                   label,
    std::set<int>&                              contactEleIds) 
{ 
  switch (treeNodeType_) 
  {
    case INNER_NODE:
    {     
      int index = -1;
      if(classifyKDOP(index, queryKDOP))
      {
        children_[index]->searchMultibodyContactElements(currentKDOPs, queryKDOP, label, contactEleIds);
        return;
      }
      else
      {
        GEO::searchMultibodyContactElements(currentKDOPs, queryKDOP, label, elementList_, contactEleIds);
        return;
      }
      break;
    }
    case LEAF_NODE:   
    {
      if(elementList_.empty())
        return;

      // max depth reached, counts reverse
      if (treedepth_ <= 0 || (elementList_.begin()->second).size() == 1)
      {
        GEO::searchMultibodyContactElements(currentKDOPs, queryKDOP, label, elementList_, contactEleIds);
        return;
      }
      // dynamically grow tree otherwise, create children and set label for empty children
      // search in apropriate child node  
      int index = -1;
      if(classifyKDOP(index, queryKDOP))
      {
        createChildren(currentKDOPs);
        children_[index]->searchMultibodyContactElements(currentKDOPs, queryKDOP, label, contactEleIds);
        return;
      }
      else
      {
        GEO::searchMultibodyContactElements(currentKDOPs, queryKDOP, label, elementList_, contactEleIds);
        return;
      }
      break;
    }
    default:
      dserror("should not get here\n");
  }
  return;
}



/*----------------------------------------------------------------------*
 | print tree node to gmsh file                            peder   07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::printTreeNode(
    const int       max_depth,
    stringstream&   fc) const
{
  if(treedepth_== max_depth)
  {
    LINALG::Matrix<3,8> printBox(true);
    printBox(0,0) = nodeBox_(0,0); printBox(1,0) = nodeBox_(1,0); printBox(2,0) = nodeBox_(2,0);
    printBox(0,1) = nodeBox_(0,0); printBox(1,1) = nodeBox_(1,1); printBox(2,1) = nodeBox_(2,0);
    printBox(0,2) = nodeBox_(0,0); printBox(1,2) = nodeBox_(1,1); printBox(2,2) = nodeBox_(2,1);
    printBox(0,3) = nodeBox_(0,0); printBox(1,3) = nodeBox_(1,0); printBox(2,3) = nodeBox_(2,1);
    printBox(0,4) = nodeBox_(0,1); printBox(1,4) = nodeBox_(1,0); printBox(2,4) = nodeBox_(2,0);
    printBox(0,5) = nodeBox_(0,1); printBox(1,5) = nodeBox_(1,1); printBox(2,5) = nodeBox_(2,0);
    printBox(0,6) = nodeBox_(0,1); printBox(1,6) = nodeBox_(1,1); printBox(2,6) = nodeBox_(2,1);
    printBox(0,7) = nodeBox_(0,1); printBox(1,7) = nodeBox_(1,0); printBox(2,7) = nodeBox_(2,1);
    fc << IO::GMSH::cellWithScalarToString(DRT::Element::hex8, 0, printBox)<< endl;
  }

  if(treeNodeType_== GEO::INNER_NODE)
  {
    for (int i=0; i < getNumChildren(); i++)
      if (children_[i] != Teuchos::null)
        children_[i]->printTreeNode(max_depth, fc);
   
  }
  else if (treeNodeType_== GEO::LEAF_NODE)
  {
    int factor = -1;
    
    if (label_ < 0)
      factor = 0; // more than one candidate in this octant
    else if (label_ == 0)
      factor = 1; // fluid octant
    else
      factor = 2; // solid octant
    LINALG::Matrix<3,8> printBox(true);
    printBox(0,0) = nodeBox_(0,0); printBox(1,0) = nodeBox_(1,0); printBox(2,0) = nodeBox_(2,0);
    printBox(0,1) = nodeBox_(0,0); printBox(1,1) = nodeBox_(1,1); printBox(2,1) = nodeBox_(2,0);
    printBox(0,2) = nodeBox_(0,0); printBox(1,2) = nodeBox_(1,1); printBox(2,2) = nodeBox_(2,1);
    printBox(0,3) = nodeBox_(0,0); printBox(1,3) = nodeBox_(1,0); printBox(2,3) = nodeBox_(2,1);
    printBox(0,4) = nodeBox_(0,1); printBox(1,4) = nodeBox_(1,0); printBox(2,4) = nodeBox_(2,0);
    printBox(0,5) = nodeBox_(0,1); printBox(1,5) = nodeBox_(1,1); printBox(2,5) = nodeBox_(2,0);
    printBox(0,6) = nodeBox_(0,1); printBox(1,6) = nodeBox_(1,1); printBox(2,6) = nodeBox_(2,1);
    printBox(0,7) = nodeBox_(0,1); printBox(1,7) = nodeBox_(1,0); printBox(2,7) = nodeBox_(2,1);
    fc << IO::GMSH::cellWithScalarToString(DRT::Element::hex8, factor + treedepth_ + max_depth, printBox)<< endl;
  }
}



/*----------------------------------------------------------------------*
 | get depth of tree node (METRICS)                        peder   07/08|
 *----------------------------------------------------------------------*/
int GEO::SearchTree::TreeNode::getDepth() const
{
  int depth = -1;
  
  if (treeNodeType_==LEAF_NODE)
    depth = treedepth_;
  else 
  {
    int tmp_depth = -1;
    depth = treedepth_;
    for (int i = 0; i < getNumChildren(); i++)
    {
      if (children_[i] != Teuchos::null)
        tmp_depth = children_[i]->getDepth();
      else 
        dserror("inner node has to have leave nodes");
      
      if(tmp_depth < depth)
        depth = tmp_depth;
    }
  }
  return depth;
}



/*----------------------------------------------------------------------*
 | get number of tree nodes in tree (METRICS)              u.may   08/08|
 *----------------------------------------------------------------------*/
int GEO::SearchTree::TreeNode::getNumNodesInTree() const
{
  int numTreeNodes = 0;
  
  switch(treeNodeType_)
  { 
    case LEAF_NODE:
    {
      numTreeNodes = 1;
      break;
    }
    case INNER_NODE:
    {
      for(int i = 0; i < getNumChildren(); i++)
        numTreeNodes += children_[i]->getNumNodesInTree();
  
      break;
    }
    default:
      dserror("wrong tree node type");
  }
  return numTreeNodes;
}


#endif  // #ifdef CCADISCRET

