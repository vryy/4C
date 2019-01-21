/*!----------------------------------------------------------------------
\file mortar_binarytree.cpp

\brief A class for performing mortar search in 2D/3D based on binarytrees

\level 1

\maintainer Alexander Popp
*----------------------------------------------------------------------*/

#include "mortar_binarytree.H"

#include "mortar_defines.H"
#include "mortar_element.H"
#include "mortar_node.H"

#include "../drt_contact/contact_defines.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_fixedsizematrix.H"


/*----------------------------------------------------------------------*
 |  ctor BinaryTreeNode (public)                              popp 10/08|
 *----------------------------------------------------------------------*/
MORTAR::BinaryTreeNode::BinaryTreeNode(MORTAR::BinaryTreeNodeType type,
    DRT::Discretization& discret, Teuchos::RCP<BinaryTreeNode> parent, std::vector<int> elelist,
    const Epetra_SerialDenseMatrix& dopnormals, const int& kdop, const int& dim,
    const bool& useauxpos, const int layer,
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& streenodesmap,
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& mtreenodesmap,
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& sleafsmap,
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& mleafsmap)
    : MORTAR::BaseBinaryTreeNode::BaseBinaryTreeNode(
          discret, elelist, dopnormals, kdop, dim, useauxpos, layer),
      type_(type),
      parent_(parent),
      streenodesmap_(streenodesmap),
      mtreenodesmap_(mtreenodesmap),
      sleafsmap_(sleafsmap),
      mleafsmap_(mleafsmap)
{
  // keep the constructor clean
  return;
}


/*----------------------------------------------------------------------*
 | get communicator (public)                                  popp 10/08|
 *----------------------------------------------------------------------*/
const Epetra_Comm& MORTAR::BinaryTreeNode::Comm() const { return Discret().Comm(); }

/*----------------------------------------------------------------------*
 | Initialize tree (public)                                   popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTreeNode::InitializeTree(double& enlarge)
{
  // return if proc. has no elements!
  if (Elelist().size() == 0) return;

  // calculate bounding volume
  CalculateSlabsDop();
  EnlargeGeometry(enlarge);

  // if current treenode is inner treenode
  if (type_ == SLAVE_INNER || type_ == MASTER_INNER)
  {
    // divide treenode
    DivideTreeNode();

    // check what to do with left child
    if (leftchild_->Elelist().size() == 0)
    {
      dserror("ERROR: InitializeTree:Processor has no leftchild elements-->return;");
      return;
    }
    else
    {
      // if leftchild is slave leaf
      if (leftchild_->Type() == SLAVE_LEAF) sleafsmap_[0].push_back(leftchild_);
      // if leaftchild is master leaf
      if (leftchild_->Type() == MASTER_LEAF) mleafsmap_[0].push_back(leftchild_);

      // recursively initialize the whole tree
      leftchild_->InitializeTree(enlarge);
    }

    // check what to do with right child
    if (rightchild_->Elelist().size() == 0)
    {
      dserror("ERROR: InitializeTree:Processor has no rightchild elements-->return;");
      return;
    }
    else
    {
      // if rightchild is slave leaf
      if (rightchild_->Type() == SLAVE_LEAF) sleafsmap_[1].push_back(rightchild_);
      // if rightchild is master leaf
      if (rightchild_->Type() == MASTER_LEAF) mleafsmap_[1].push_back(rightchild_);

      // recursively initialize the whole tree
      rightchild_->InitializeTree(enlarge);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Update slabs bottom up (public)                            popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTreeNode::UpdateSlabsBottomUp(double& enlarge)
{
  // if current treenode is inner node
  if (type_ == SLAVE_INNER || type_ == MASTER_INNER)
  {
    for (int k = 0; k < Kdop() / 2; k++)
    {
      // for minimum
      if (leftchild_->Slabs()(k, 0) <= rightchild_->Slabs()(k, 0))
        Slabs()(k, 0) = leftchild_->Slabs()(k, 0);
      else
        Slabs()(k, 0) = rightchild_->Slabs()(k, 0);

      // for maximum
      if (leftchild_->Slabs()(k, 1) >= rightchild_->Slabs()(k, 1))
        Slabs()(k, 1) = leftchild_->Slabs()(k, 1);
      else
        Slabs()(k, 1) = rightchild_->Slabs()(k, 1);
    }
  }

  // if current treenode is leafnode
  if (type_ == SLAVE_LEAF || type_ == MASTER_LEAF)
  {
    CalculateSlabsDop();

    EnlargeGeometry(enlarge);

    // Prints Slabs to std::cout
    // PrintSlabs();
  }  // current treenode is leaf

  return;
}
/*----------------------------------------------------------------------*
 | Divide treenode (public)                                   popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTreeNode::DivideTreeNode()
{
  // map of elements belonging to left / right child treenode
  std::vector<int> leftelements(0);
  std::vector<int> rightelements(0);

  // if only 2 elements in Treenode, create new treenodes out of them
  if (Elelist().size() == 2)
  {
    leftelements.push_back(Elelist()[0]);
    rightelements.push_back(Elelist()[1]);
  }
  // if more than 2 elements in Treenode
  else if (Elelist().size() > 2)
  {
    // calculate splitting area (split along longest side)
    double lmax = 0.0;                    // max. length of sides of DOP
    int splittingnormal = -1;             // defines side to split
    double xmedian[3] = {0.0, 0.0, 0.0};  // coordinates of centroid

    for (int i = 0; i < Kdop() / 2; ++i)
    {
      double lcurrent = abs(Slabs()(i, 1) - Slabs()(i, 0));
      if (lmax < lcurrent)
      {
        lmax = lcurrent;
        splittingnormal = i;
      }
    }

    // find median of centroid coordinates to divide area
    // coordinates of median
    for (int i = 0; i < Dim(); ++i)
      xmedian[i] = Slabs()(i, 1) - ((Slabs()(i, 1) - Slabs()(i, 0)) / 2);

    // compute d of ax+by+cz=d of splittingplane
    double d = xmedian[0] * Dopnormals()(splittingnormal, 0) +
               xmedian[1] * Dopnormals()(splittingnormal, 1) +
               xmedian[2] * Dopnormals()(splittingnormal, 2);

    // split treenode into two parts
    for (int i = 0; i < ((int)Elelist().size()); ++i)
    {
      bool isright = false;  // true, if element should be sorted into right treenode
      bool isleft = false;   // true, if element should be sorted into left treenode

      int gid = Elelist()[i];
      DRT::Element* element = Discret().gElement(gid);
      if (!element) dserror("ERROR: Cannot find element with gid %\n", gid);
      DRT::Node** nodes = element->Points();

      // vector of values of Hesse-Normalform of nodes of elements
      Epetra_SerialDenseVector axbycz;
      axbycz.Resize(element->NumPoint());

      for (int k = 0; k < element->NumPoint(); ++k)
      {
        MortarNode* mrtrnode = dynamic_cast<MortarNode*>(nodes[k]);
        if (!mrtrnode) dserror("ERROR: Null pointer!");
        const double* posnode = mrtrnode->X();

        // split along chosen area
        // ax+by+cz< or > d = criterion
        // compute ax+by+cz for chosen node
        if (Dim() == 2)
          axbycz[k] = posnode[0] * Dopnormals()(splittingnormal, 0) +
                      posnode[1] * Dopnormals()(splittingnormal, 1);
        else if (Dim() == 3)
          axbycz[k] = posnode[0] * Dopnormals()(splittingnormal, 0) +
                      posnode[1] * Dopnormals()(splittingnormal, 1) +
                      posnode[2] * Dopnormals()(splittingnormal, 2);
        else
          dserror("ERROR: Problem dimension must be 2D or 3D!");

        if (axbycz[k] >= d) isright = true;
        if (axbycz[k] < d) isleft = true;
      }

      if (isright == false && isleft == false)
        dserror("ERROR: Current element could neither be sorted into left- or right-child node!");

      // if element is split through, it is sorted into left treenode
      if (isright == true && isleft == true) isright = false;

      // sort elements into child treenodes
      if (isright) rightelements.push_back(gid);
      if (isleft) leftelements.push_back(gid);
    }
  }

  // if treenode splitting algorithm was not able to divide treenode
  // successfully (i.e all elements are in one child treenode),
  // then just put one element into the other treenode
  if (leftelements.size() == 0 && rightelements.size() > 1)
  {
    leftelements.push_back(rightelements[rightelements.size() - 1]);
    rightelements.pop_back();
  }
  if (rightelements.size() == 0 && leftelements.size() > 1)
  {
    rightelements.push_back(leftelements[leftelements.size() - 1]);
    leftelements.pop_back();
  }

  // define type of newly created children treenodes
  if (Elelist().size() >= 2)
  {
    // defines type of left and right TreeNode
    BinaryTreeNodeType lefttype = UNDEFINED;
    BinaryTreeNodeType righttype = UNDEFINED;

    // is the new left child treenode a leaf node?
    if (leftelements.size() == 1)
    {
      if (type_ == SLAVE_INNER)
        lefttype = SLAVE_LEAF;
      else if (type_ == MASTER_INNER)
        lefttype = MASTER_LEAF;
      else
        dserror("ERROR: Invalid TreeNodeType");
    }
    else
    {
      if (type_ == SLAVE_INNER)
        lefttype = SLAVE_INNER;
      else if (type_ == MASTER_INNER)
        lefttype = MASTER_INNER;
      else
        dserror("ERROR: Invalid TreeNodeType");
    }

    // is the new right child treenode a leaf node?
    if (rightelements.size() == 1)
    {
      if (type_ == SLAVE_INNER)
        righttype = SLAVE_LEAF;
      else if (type_ == MASTER_INNER)
        righttype = MASTER_LEAF;
      else
        dserror("ERROR: Invalid TreeNodeType");
    }
    else
    {
      if (type_ == SLAVE_INNER)
        righttype = SLAVE_INNER;
      else if (type_ == MASTER_INNER)
        righttype = MASTER_INNER;
      else
        dserror("ERROR: Invalid TreeNodeType");
    }

    // build left child treenode
    leftchild_ = Teuchos::rcp(new BinaryTreeNode(lefttype, Discret(), Teuchos::rcp(this, false),
        leftelements, Dopnormals(), Kdop(), Dim(), UseAuxPos(), (Layer() + 1), streenodesmap_,
        mtreenodesmap_, sleafsmap_, mleafsmap_));

    // build right child treenode
    rightchild_ = Teuchos::rcp(new BinaryTreeNode(righttype, Discret(), Teuchos::rcp(this, false),
        rightelements, Dopnormals(), Kdop(), Dim(), UseAuxPos(), (Layer() + 1), streenodesmap_,
        mtreenodesmap_, sleafsmap_, mleafsmap_));

    // update slave and mastertreenodes map
    // if parent treenode is slave
    if (type_ == SLAVE_INNER)
    {
      // if map of treenodes does not have enogh rows-->resize!
      if ((int)(streenodesmap_.size()) <= (Layer() + 1)) streenodesmap_.resize((Layer() + 2));

      // put new pointers to children into map
      streenodesmap_[(Layer() + 1)].push_back(leftchild_);
      streenodesmap_[(Layer() + 1)].push_back(rightchild_);
    }

    // if parent treenode is master
    if (type_ == MASTER_INNER)
    {
      // if map of treenodes does not have enogh rows-->resize!
      if ((int)(mtreenodesmap_.size()) <= (Layer() + 1)) mtreenodesmap_.resize((Layer() + 2));

      // put new pointers to children into map
      mtreenodesmap_[(Layer() + 1)].push_back(leftchild_);
      mtreenodesmap_[(Layer() + 1)].push_back(rightchild_);
    }
  }

  else
    dserror("ERROR: Only 1 or 0 elements in map-->TreeNode cannot be devided!!");

  return;
}

/*----------------------------------------------------------------------*
 | Print type of treenode to std::cout (public)               popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTreeNode::PrintType()
{
  switch (type_)
  {
    case SLAVE_INNER:
      std::cout << std::endl << "SLAVE_INNER ";
      break;
    case SLAVE_LEAF:
      std::cout << std::endl << "SLAVE_LEAF ";
      break;
    case MASTER_INNER:
      std::cout << std::endl << "MASTER_INNER ";
      break;
    case MASTER_LEAF:
      std::cout << std::endl << "MASTER_LEAF ";
      break;
    case NOSLAVE_ELEMENTS:
      std::cout << std::endl << "TreeNode contains no Slave-Elements=NO_SLAVEELEMENTS ";
      break;
    case NOMASTER_ELEMENTS:
      std::cout << std::endl << "TreeNode contains no Master-Elements=NO_MASTERELEMENTS ";
      break;
    case UNDEFINED:
      std::cout << std::endl << "UNDEFINED ";
      break;
    default:
      dserror("Unknown MORTAR::BinaryTreeNodeType detected.");
      break;
  }

  return;
}

/*----------------------------------------------------------------------*
 | Set slabs of current treenode with new slabs(public)       popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTreeNode::SetSlabs(Epetra_SerialDenseMatrix& newslabs)
{
  for (int i = 0; i < Kdop() / 2; ++i)
  {
    Slabs()(i, 0) = newslabs(i, 0);
    Slabs()(i, 1) = newslabs(i, 1);
  }
}

/*----------------------------------------------------------------------*
 |  ctor BinaryTree(public)                                   popp 10/08|
 *----------------------------------------------------------------------*/
MORTAR::BinaryTree::BinaryTree(DRT::Discretization& discret, Teuchos::RCP<Epetra_Map> selements,
    Teuchos::RCP<Epetra_Map> melements, int dim, double eps, bool useauxpos)
    : MORTAR::BaseBinaryTree(discret, dim, eps),
      selements_(selements),
      melements_(melements),
      useauxpos_(useauxpos)
{
  // keep the constructor clean
  return;
}

/*----------------------------------------------------------------------*
 | initialize the binary tree                              schmidt 12/18|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTree::Init()
{
  // call initialization method of base class
  MORTAR::BaseBinaryTree::Init();

  InitInternalVariables();

  // calculate minimal element length
  SetEnlarge();

  //**********************************************************************
  // initialize binary tree root nodes
  //**********************************************************************
  // create element lists
  std::vector<int> slist;
  std::vector<int> mlist;

  for (int i = 0; i < selements_->NumMyElements(); ++i)
  {
    int gid = selements_->GID(i);
    slist.push_back(gid);
  }

  for (int i = 0; i < melements_->NumMyElements(); ++i)
  {
    int gid = melements_->GID(i);
    mlist.push_back(gid);
  }

  // check slave root node case
  if (slist.size() >= 2)
  {
    sroot_ =
        Teuchos::rcp(new BinaryTreeNode(MORTAR::SLAVE_INNER, Discret(), sroot_, slist, DopNormals(),
            Kdop(), Dim(), useauxpos_, 0, streenodesmap_, mtreenodesmap_, sleafsmap_, mleafsmap_));

    // do initialization
    streenodesmap_[0].push_back(sroot_);
    sroot_->InitializeTree(Enlarge());
  }
  else if (slist.size() == 1)
  {
    sroot_ =
        Teuchos::rcp(new BinaryTreeNode(MORTAR::SLAVE_LEAF, Discret(), sroot_, slist, DopNormals(),
            Kdop(), Dim(), useauxpos_, 0, streenodesmap_, mtreenodesmap_, sleafsmap_, mleafsmap_));

    // trivial initialization
    streenodesmap_[0].push_back(sroot_);
    sleafsmap_[0].push_back(sroot_);
  }
  else
  {
    sroot_ = Teuchos::rcp(
        new BinaryTreeNode(MORTAR::NOSLAVE_ELEMENTS, Discret(), sroot_, slist, DopNormals(), Kdop(),
            Dim(), useauxpos_, 0, streenodesmap_, mtreenodesmap_, sleafsmap_, mleafsmap_));

    // trivial initialization
    streenodesmap_[0].push_back(sroot_);
  }

  // check master root node case
  if (mlist.size() >= 2)
  {
    mroot_ = Teuchos::rcp(
        new BinaryTreeNode(MORTAR::MASTER_INNER, Discret(), mroot_, mlist, DopNormals(), Kdop(),
            Dim(), useauxpos_, 0, streenodesmap_, mtreenodesmap_, sleafsmap_, mleafsmap_));

    // do initialization
    mtreenodesmap_[0].push_back(mroot_);
    mroot_->InitializeTree(Enlarge());
  }
  else if (mlist.size() == 1)
  {
    mroot_ =
        Teuchos::rcp(new BinaryTreeNode(MORTAR::MASTER_LEAF, Discret(), mroot_, mlist, DopNormals(),
            Kdop(), Dim(), useauxpos_, 0, streenodesmap_, mtreenodesmap_, sleafsmap_, mleafsmap_));

    // trivial initialization
    mtreenodesmap_[0].push_back(mroot_);
    mleafsmap_[0].push_back(mroot_);
  }
  else
  {
    mroot_ = Teuchos::rcp(
        new BinaryTreeNode(MORTAR::NOMASTER_ELEMENTS, Discret(), mroot_, mlist, DopNormals(),
            Kdop(), Dim(), useauxpos_, 0, streenodesmap_, mtreenodesmap_, sleafsmap_, mleafsmap_));

    // trivial initialization / error
    mtreenodesmap_[0].push_back(mroot_);

    // no error: round robin loop can handle this case!
    // dserror("No master element for Binarytree initialization on this processor");
  }

  /*
  // print binarytree to std::cout
  for (int k=0;k<Comm().NumProc();++k)
  {
    Comm().Barrier();
    if (Comm().MyPID()==k)
    {
      std::cout << "\n" << Comm().MyPID() << " Print tree with direct print function" << std::endl;
      std::cout <<"\n" <<Comm().MyPID()<< " Slave Tree:";
      PrintTree(sroot_);
      std::cout <<"\n" <<Comm().MyPID()<< " Master Tree:";
      PrintTree(mroot_);
    }
    Comm().Barrier();
  }

  for (int k=0;k<Comm().NumProc();++k)
  {
    Comm().Barrier();
    if (Comm().MyPID()==k)
    {
      std::cout << "\n" << Comm().MyPID() << " Print tree with print function of slave and master
  treemap" << std::endl; std::cout <<"\n" <<Comm().MyPID()<< " Slave Tree:";
      PrintTreeOfMap(streenodesmap_);
      std::cout <<"\n" <<Comm().MyPID()<< " Master Tree:";
      PrintTreeOfMap(mtreenodesmap_);
    }
    Comm().Barrier();
  }
  */

  return;
}

/*----------------------------------------------------------------------*
 |  Initialize internal variables (private)                schmidt 01/19|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTree::InitInternalVariables()
{
  // initialize sizes
  streenodesmap_.resize(1);
  mtreenodesmap_.resize(1);
  couplingmap_.resize(2);
  sleafsmap_.resize(2);
  mleafsmap_.resize(2);

  return;
}

/*----------------------------------------------------------------------*
 | clear found search elements from last eval.               farah 01/16|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTree::InitSearchElements()
{
  // loop over all elements to reset candidates / search lists
  // (use standard slave column map)
  for (int i = 0; i < selements_->NumMyElements(); ++i)
  {
    int gid = selements_->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i", gid);
    MortarElement* sele = dynamic_cast<MortarElement*>(ele);

    sele->MoData().SearchElements().resize(0);
  }
}


/*----------------------------------------------------------------------*
 | evaluate search to get sele to mele info (public)         farah 01/16|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTree::EvaluateSearch()
{
  // *********************************************************************
  // Possible versions for general mortar setting:
  // *********************************************************************
  //
  // 1) Combined Update and Search
  // -> In this case we only have to call SearchCombined(), which
  //    does buth top-down update (where necessary) and search.
  //
  // 2) Separate Update and Search
  // -> In this case we have to explicitly call and updating routine, i.e.
  //    UpdateTreeTopDown() or UpdateTreeBottomUp() before calling the
  //    search routine SearchSeparate(). Of course, the bottom-up
  //    update makes more sense here. For very large meshtying problems,
  //    this version is preferable and thus chosen as default.
  //
  // *********************************************************************

  // init search elements
  InitSearchElements();

  // calculate minimal element length
  SetEnlarge();

  // update tree in a top down way
  // binarytree_->UpdateTreeTopDown();

  // update tree in a bottom up way
  UpdateTreeBottomUp();

#ifdef MORTARGMSHCTN
  for (int i = 0; i < (int)(binarytree_->CouplingMap().size()); i++)
    binarytree_->CouplingMap()[i].clear();
  binarytree_->CouplingMap().clear();
  binarytree_->CouplingMap().resize(2);
#endif  // MORTARGMSHCTN
  // search with a separate algorithm
  SearchSeparate();

  // search with an combined algorithm
  // binarytree_->SearchCombined();
  return;
}


/*----------------------------------------------------------------------*
 | get communicator (public)                                  popp 10/08|
 *----------------------------------------------------------------------*/
const Epetra_Comm& MORTAR::BinaryTree::Comm() const { return Discret().Comm(); }

/*----------------------------------------------------------------------*
 | Find min. length of master and slave elements (public)     popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTree::SetEnlarge()
{
  double lmin = 1.0e12;

  // calculate mininmal length of slave elements
  for (int i = 0; i < selements_->NumMyElements(); ++i)
  {
    int gid = selements_->GID(i);
    DRT::Element* element = Discret().gElement(gid);
    if (!element) dserror("ERROR: Cannot find element with gid %\n", gid);
    MORTAR::MortarElement* mrtrelement = dynamic_cast<MORTAR::MortarElement*>(element);
    double mincurrent = mrtrelement->MinEdgeSize();
    if (mincurrent < lmin) lmin = mincurrent;
  }

  // calculate minimal length of master elements
  for (int i = 0; i < melements_->NumMyElements(); ++i)
  {
    int gid = melements_->GID(i);
    DRT::Element* element = Discret().gElement(gid);
    if (!element) dserror("ERROR: Cannot find element with gid %\n", gid);
    MORTAR::MortarElement* mrtrelement = dynamic_cast<MORTAR::MortarElement*>(element);
    double mincurrent = mrtrelement->MinEdgeSize();
    if (mincurrent < lmin) lmin = mincurrent;
  }

  if (lmin <= 0.0) dserror("ERROR: Minimal element length <= 0!");

  // set the class variables
  Enlarge() = Eps() * lmin;

  return;
}

/*----------------------------------------------------------------------*
 | Print tree (public)                                        popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTree::PrintTree(Teuchos::RCP<BinaryTreeNode> treenode)
{
  // if treenode has no elements (NOSLAVE_ELEMENTS,NOMASTER_ELEMENTS)
  if (treenode->Type() == NOSLAVE_ELEMENTS || treenode->Type() == NOMASTER_ELEMENTS)
  {
    std::cout << "\n" << Comm().MyPID() << " Tree has no element to print";
    return;
  }
  std::cout << "\n" << Comm().MyPID() << " Tree at layer: " << treenode->Layer() << " Elements: ";
  for (int i = 0; i < (int)(treenode->Elelist().size()); i++)
    std::cout << " " << treenode->Elelist()[i];

  // while treenode is inner node
  if (treenode->Type() == SLAVE_INNER || treenode->Type() == MASTER_INNER)
  {
    PrintTree(treenode->Leftchild());
    PrintTree(treenode->Rightchild());
  }

  return;
}

/*----------------------------------------------------------------------*
 | Print tree with treenodesmap_ (public)                     popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTree::PrintTreeOfMap(
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& treenodesmap)
{
  // print tree, elements listet in brackets (), belong to one treenode!
  for (int i = 0; i < (int)(treenodesmap.size()); i++)
  {
    std::cout << "\n" << Comm().MyPID() << " Tree at layer: " << i << " Elements: ";
    for (int k = 0; k < (int)(treenodesmap[i].size()); k++)
    {
      Teuchos::RCP<BinaryTreeNode> currentnode = treenodesmap[i][k];
      std::cout << " (";
      for (int l = 0; l < (int)(currentnode->Elelist().size()); l++)
      {
        std::cout << currentnode->Elelist()[l] << " ";
      }
      std::cout << ") ";
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Update tree topdown (public)                               popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTree::EvaluateUpdateTreeTopDown(Teuchos::RCP<BinaryTreeNode> treenode)
{
  // if no slave element on proc-->return
  if (treenode->Elelist().size() == 0) return;

  treenode->CalculateSlabsDop();
  treenode->EnlargeGeometry(Enlarge());

  if (treenode->Type() == SLAVE_INNER || treenode->Type() == MASTER_INNER)
  {
    EvaluateUpdateTreeTopDown(treenode->Leftchild());
    EvaluateUpdateTreeTopDown(treenode->Rightchild());
  }

  return;
}

/*----------------------------------------------------------------------*
 | Update tree bottom up based on list (public)               popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTree::EvaluateUpdateTreeBottomUp(
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& treenodesmap)
{
  // update tree bottom up (for every treelayer)
  for (int i = ((int)(treenodesmap.size() - 1)); i >= 0; i = i - 1)
  {
    for (int j = 0; j < (int)(treenodesmap[i].size()); j++)
      treenodesmap[i][j]->UpdateSlabsBottomUp(Enlarge());
  }

  return;
}

/*----------------------------------------------------------------------*
 | Search for contact (public)                                popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTree::EvaluateSearchSeparate(
    Teuchos::RCP<BinaryTreeNode> streenode, Teuchos::RCP<BinaryTreeNode> mtreenode)
{
  // tree needs to be updated before running contact search!

  // if there are no elements
  if (streenode->Type() == NOSLAVE_ELEMENTS || mtreenode->Type() == NOMASTER_ELEMENTS) return;

  // check if treenodes intercept
  // (they only intercept if ALL slabs intercept!)
  int nintercepts = 0;

  for (int i = 0; i < Kdop() / 2; ++i)
  {
    if (streenode->Slabs()(i, 0) <= mtreenode->Slabs()(i, 0))
    {
      if (streenode->Slabs()(i, 1) >= mtreenode->Slabs()(i, 0))
        nintercepts++;
      else if (streenode->Slabs()(i, 1) >= mtreenode->Slabs()(i, 1))
        nintercepts++;
    }
    else if (streenode->Slabs()(i, 0) >= mtreenode->Slabs()(i, 0))
    {
      if (mtreenode->Slabs()(i, 1) >= streenode->Slabs()(i, 1))
        nintercepts++;
      else if (mtreenode->Slabs()(i, 1) >= streenode->Slabs()(i, 0))
        nintercepts++;
    }
  }

  // treenodes intercept
  if (nintercepts == Kdop() / 2)
  {
    // slave and master treenodes are inner nodes
    if (streenode->Type() == SLAVE_INNER && mtreenode->Type() == MASTER_INNER)
    {
      // std::cout <<"\n"<< Comm().MyPID() << " 2 inner nodes!";
      EvaluateSearchSeparate(streenode->Leftchild(), mtreenode->Leftchild());
      EvaluateSearchSeparate(streenode->Leftchild(), mtreenode->Rightchild());
      EvaluateSearchSeparate(streenode->Rightchild(), mtreenode->Leftchild());
      EvaluateSearchSeparate(streenode->Rightchild(), mtreenode->Rightchild());
    }

    // slave treenode is inner, master treenode is leaf
    if (streenode->Type() == SLAVE_INNER && mtreenode->Type() == MASTER_LEAF)
    {
      // std::cout <<"\n"<< Comm().MyPID() << " slafe inner, master leaf!";
      EvaluateSearchSeparate(streenode->Leftchild(), mtreenode);
      EvaluateSearchSeparate(streenode->Rightchild(), mtreenode);
    }

    // slave treenode is leaf,  master treenode is inner
    if (streenode->Type() == SLAVE_LEAF && mtreenode->Type() == MASTER_INNER)
    {
      // std::cout <<"\n"<< Comm().MyPID() << " slave leaf, master inner!";
      EvaluateSearchSeparate(streenode, mtreenode->Leftchild());
      EvaluateSearchSeparate(streenode, mtreenode->Rightchild());
    }

    // both treenodes are leaf --> feasible pair
    if (streenode->Type() == SLAVE_LEAF && mtreenode->Type() == MASTER_LEAF)
    {
      int sgid = (int)streenode->Elelist()[0];  // global id of slave element
      int mgid = (int)mtreenode->Elelist()[0];  // global id of masterelement
      // std::cout <<"\n"<< Comm().MyPID() << "TreeDividedContact found between slave-Element: "
      //     << sgid <<"and master-Element: "<< mgid;
      DRT::Element* element = Discret().gElement(sgid);
      MORTAR::MortarElement* selement = dynamic_cast<MORTAR::MortarElement*>(element);
      selement->AddSearchElements(mgid);
    }
  }

#ifdef MORTARGMSHCTN  // for plotting contacting treenodes
  if (streenode->Type() == SLAVE_LEAF && mtreenode->Type() == MASTER_LEAF &&
      nintercepts == kdop_ / 2)
  {
    couplingmap_[0].push_back(streenode);
    couplingmap_[1].push_back(mtreenode);
  }
#endif  // #ifdef MORTARGMSHCTN

  return;
}

/*----------------------------------------------------------------------*
 | Search for contact (public)                                popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BinaryTree::EvaluateSearchCombined(
    Teuchos::RCP<BinaryTreeNode> streenode, Teuchos::RCP<BinaryTreeNode> mtreenode)
{
  // root nodes need to be updated before running combined contact search!

  // if there are no elements
  if (streenode->Type() == NOSLAVE_ELEMENTS || mtreenode->Type() == NOMASTER_ELEMENTS) return;

  // check if treenodes intercept
  // (they only intercept if ALL slabs intercept!)
  int nintercepts = 0;

  for (int i = 0; i < Kdop() / 2; ++i)
  {
    if (streenode->Slabs()(i, 0) <= mtreenode->Slabs()(i, 0))
    {
      if (streenode->Slabs()(i, 1) >= mtreenode->Slabs()(i, 0))
        nintercepts++;
      else if (streenode->Slabs()(i, 1) >= mtreenode->Slabs()(i, 1))
        nintercepts++;
    }
    else if (streenode->Slabs()(i, 0) >= mtreenode->Slabs()(i, 0))
    {
      if (mtreenode->Slabs()(i, 1) >= streenode->Slabs()(i, 1))
        nintercepts++;
      else if (mtreenode->Slabs()(i, 1) >= streenode->Slabs()(i, 0))
        nintercepts++;
    }
  }

  // treenodes intercept
  if (nintercepts == Kdop() / 2)
  {
    // slave and master treenodes are inner nodes
    if (streenode->Type() == SLAVE_INNER && mtreenode->Type() == MASTER_INNER)
    {
      // std::cout <<"\n"<< Comm().MyPID() << " 2 inner nodes!";
      streenode->Leftchild()->CalculateSlabsDop();
      streenode->Leftchild()->EnlargeGeometry(Enlarge());
      streenode->Rightchild()->CalculateSlabsDop();
      streenode->Rightchild()->EnlargeGeometry(Enlarge());
      mtreenode->Leftchild()->CalculateSlabsDop();
      mtreenode->Leftchild()->EnlargeGeometry(Enlarge());
      mtreenode->Rightchild()->CalculateSlabsDop();
      mtreenode->Rightchild()->EnlargeGeometry(Enlarge());

      EvaluateSearchCombined(streenode->Leftchild(), mtreenode->Leftchild());
      EvaluateSearchCombined(streenode->Leftchild(), mtreenode->Rightchild());
      EvaluateSearchCombined(streenode->Rightchild(), mtreenode->Leftchild());
      EvaluateSearchCombined(streenode->Rightchild(), mtreenode->Rightchild());
    }

    // slave treenode is inner,  master treenode is leaf
    if (streenode->Type() == SLAVE_INNER && mtreenode->Type() == MASTER_LEAF)
    {
      // std::cout <<"\n"<< Comm().MyPID() << " slafe inner, master leaf!";
      streenode->Leftchild()->CalculateSlabsDop();
      streenode->Leftchild()->EnlargeGeometry(Enlarge());
      streenode->Rightchild()->CalculateSlabsDop();
      streenode->Rightchild()->EnlargeGeometry(Enlarge());

      EvaluateSearchCombined(streenode->Leftchild(), mtreenode);
      EvaluateSearchCombined(streenode->Rightchild(), mtreenode);
    }

    // slave treenode is leaf,  master treenode is inner
    if (streenode->Type() == SLAVE_LEAF && mtreenode->Type() == MASTER_INNER)
    {
      // std::cout <<"\n"<< Comm().MyPID() << " slave leaf, master inner!";
      mtreenode->Leftchild()->CalculateSlabsDop();
      mtreenode->Leftchild()->EnlargeGeometry(Enlarge());
      mtreenode->Rightchild()->CalculateSlabsDop();
      mtreenode->Rightchild()->EnlargeGeometry(Enlarge());

      EvaluateSearchCombined(streenode, mtreenode->Leftchild());
      EvaluateSearchCombined(streenode, mtreenode->Rightchild());
    }

    // both treenodes are leaf --> feasible pair
    if (streenode->Type() == SLAVE_LEAF && mtreenode->Type() == MASTER_LEAF)
    {
      // std::cout <<"\n"<< Comm().MyPID() << " 2 leaf nodes!";
      int sgid = (int)streenode->Elelist()[0];  // global id of slave element
      int mgid = (int)mtreenode->Elelist()[0];  // global id of master element
      // std::cout <<"\n"<< Comm().MyPID() << "TreeCombinedContact found between slave-Element: "
      //     << sgid <<"and master-Element: "<< mgid;
      DRT::Element* element = Discret().gElement(sgid);
      MORTAR::MortarElement* selement = dynamic_cast<MORTAR::MortarElement*>(element);
      selement->AddSearchElements(mgid);
    }
  }

#ifdef MORTARGMSHCTN  // for plotting contacting treenodes
  if (streenode->Type() == SLAVE_LEAF && mtreenode->Type() == MASTER_LEAF &&
      nintercepts == kdop_ / 2)
  {
    couplingmap_[0].push_back(streenode);
    couplingmap_[1].push_back(mtreenode);
  }
#endif  // #ifdef MORTARGMSHCTN

  return;
}
