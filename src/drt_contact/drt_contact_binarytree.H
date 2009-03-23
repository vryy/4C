/*!----------------------------------------------------------------------
\file drt_contact_binarytree.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

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
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CONTACT_BINARYTREE_H
#define DRT_CONTACT_BINARYTREE_H

#include "../drt_lib/drt_discret.H"
#include "drt_celement.H"

using namespace std;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

//! identifies node types in a tree
enum BinaryTreeNodeType {
    SLAVE_INNER,   			///< indicates a slave inner node (has children)
    SLAVE_LEAF,   	 		///< indicates a slave leaf node (no further children)
    MASTER_INNER,  			///< indicates a master inner node (has children)
    MASTER_LEAF,	    	///< indicates a master leaf node (no further children)
    NOSLAVE_ELEMENTS,		///< inidicates that there are no slave elements on this (root) treenode
    NOMASTER_ELEMENTS,	///< inidicates that there are no master elements on this (root) treenode
    UNDEFINED           ///< indiactes an undefined tree node
};
  
// forward declarations

/*!
\brief A class for performing contact search in 2D / 3D based
       on a binary search tree

\author popp (popp@lnm.mw.tum.de)
*/

class BinaryTree
{
  private: class BinaryTreeNode{
     
  private:
    
    //! type of BinaryTreeNode
    BinaryTreeNodeType type_;
    
    //! contact interface discretization
    DRT::Discretization& idiscret_;
    
    //! pointer to the parent BinaryTreeNode
    RCP<BinaryTreeNode> parent_; 
  
    //! pointer to the left child TreeNode
    RCP<BinaryTreeNode> leftchild_;
    
    //! pointer to the right child TreeNode
    RCP<BinaryTreeNode> rightchild_;
    
    //! geometry slabs of TreeNode, saved as Min|Max
    Epetra_SerialDenseMatrix slabs_;
   
    //! list with all elements of BinaryTreeNode
    vector<int> elelist_;
    
    //! reference to DOP normals
    const Epetra_SerialDenseMatrix& dopnormals_;
    
    //! reference to no. of vertices
    const int& kdop_;
    
    //! reference to dim. of problem
    const int& dim_;
    
    //! layer of TreeNode in Tree (0=rootnode!)
    const int layer_;
    
    //! reference to map of all slave treenodes, sorted by layer
    vector<vector<RCP<BinaryTreeNode> > > & streenodesmap_;	
 
    //! reference to map of all master treenodes, sorted by layer
    vector<vector<RCP<BinaryTreeNode> > > &  mtreenodesmap_;
    
    //! reference to map of all slave leaf treenodes
    vector<vector<RCP<BinaryTreeNode> > > & sleafsmap_;	
 
    //! reference to map of all master leaf treenodes
    vector<vector<RCP<BinaryTreeNode> > > &  mleafsmap_;


  public:
    
    /*!
    \brief constructor of a tree node
    
    \param type        type of BinaryTreeNode
    \param discret     contact interface discretization 
    \param parent      points to parent tree node
    \param elelist     list of all elements in BinaryTreeNode
    \param dopnormals  reference to DOP normals
    \param kdop        refernece to no. of vertices
    \param dim				 dimension of problem    
    \param layer			 current layer of treenode
    \param ...map			 references to maps
    
    */
    BinaryTreeNode(BinaryTreeNodeType type, DRT::Discretization& discret, 
                   RCP<BinaryTreeNode> parent, vector<int> elelist,
                   const Epetra_SerialDenseMatrix& dopnormals, const int& kdop,
                   const int& dim, const int layer,
                   vector<vector<RCP<BinaryTreeNode> > > & streenodesmap,
                   vector<vector<RCP<BinaryTreeNode> > > & mtreenodesmap,
                   vector<vector<RCP<BinaryTreeNode> > > & sleafsmap,
                   vector<vector<RCP<BinaryTreeNode> > > & mleafsmap);
        
    /*!
    \brief Destructor

    */
    virtual ~BinaryTreeNode() {};
    
    //! @name Evaluation methods
    
    /*!
    \brief Calculate slabs of dop in initialization

    */
    void CalculateSlabsDop(bool isinit);
    
    /*!
    \brief Update slabs of current treenode in bottom up way

    */
    void UpdateSlabsBottomUp(double & eps);

    /*!
    \brief Initialize Tree

    */    
    void InitializeTree(double& eps);
    
    /*!
    \brief Enlarge geometry of a Treenode by an offset, dependent on size

    */    
    void EnlargeGeometry(double& eps);
    
    /*!
    \brief Divide a TreeNode into two child nodes

    */ 
    void DivideTreeNode();
    
    /*!
    \brief Print type of TreeNode to std::cout

    */
    void PrintType();
    
    /*!
    \brief Print slabs to std::cout

    */
    void PrintSlabs();

    /*!
    \brief Set slabs of current treenode

    */
    void SetSlabs(Epetra_SerialDenseMatrix& newslabs);
    
    //@}  
    
    //! @name Access methods
    
    /*!
    \brief Get communicator
    
    */
    virtual const Epetra_Comm& Comm() const { return idiscret_.Comm();}
    
    /*!
    \brief Get discretization of the interface

    */
    virtual inline DRT::Discretization& Discret() const { return idiscret_; }
    
    /*!
    \brief Return pointer to type of treenode

    */
    virtual inline BinaryTreeNodeType Type() const { return type_; }
    
    /*!
    \brief Return pointer to right child

    */
    virtual inline RCP<BinaryTreeNode> Rightchild() const { return rightchild_; }
    
    /*!
    \brief Return pointer to left child

    */
    virtual inline RCP<BinaryTreeNode> Leftchild() const { return leftchild_; }
  
    /*!
    \brief Return pointer to parent

    */
    virtual inline RCP<BinaryTreeNode> Parent() const { return parent_; }
    
    /*!
    \brief Return pointer to slabs of DOP

    */
    virtual inline Epetra_SerialDenseMatrix Slabs() const { return slabs_; }
    
    /*!
    \brief Return layer of current TreeNode

    */
    virtual inline const int Layer() const { return layer_; }
    
    /*!
    \brief Return pointer to element list of TreeNode

    */
    virtual inline vector<int> Elelist() const { return elelist_; }
    
    //@}  
    
    //! @name Visualization methods (GMSH)
    
    /*!
    \brief Print slabs of DOP to file for GMSH output
    \param filename     filename to which treenodes are plotted

    */
    void PrintDopsForGmsh(std::string filename);

    /*!
    \brief Plot a point in GMSH to given file

    */
    void PlotGmshPoint(std::string filename, double* position0, int nr);
   
    /*!
    \brief Plot n quadrangle in GMSH to given file

    */
    void PlotGmshQuadrangle(std::string filename, double* position0, double* position1,
                                                  double* position2, double* position3);

    
    /*!
    \brief Plot a triangle in GMSH to given file

    */
    void PlotGmshTriangle(std::string filename, double* position0, double* position1,
                                                                   double* position2);
 
    //@}  
    
  }; // class BinaryTreeNode
  
public:
  
  /*!
  \brief Standard constructor
  
  Constructs an instance of this class.<br>
  
  \param discret (in): 		The contact interface discretization
  \param selements (in): 	All slave elements (column map)
  \param melements (in): 	All master elements (fully overlapping map)
  \param dim (in): 				The problem dimension 
  
  */
  BinaryTree(DRT::Discretization& discret,
             RCP<Epetra_Map> selements,
             RCP<Epetra_Map> melements,
             int dim,
             double eps);

  /*!
  \brief Destructor

  */
  virtual ~BinaryTree() {};
  
  //! @name Query methods
  
  /*!
  \brief Print full tree

	*/
	void PrintTree(RCP<BinaryTreeNode> treenode);
	
	/*!
	\brief Print full tree out of map of treenodes
	
	*/
	void PrintTreeOfMap(vector<vector<RCP<BinaryTreeNode> > > & treenodesmap);
  
  //@}
	
  //! @name Access methods
  
  /*!
  \brief Get communicator

  */
  virtual const Epetra_Comm& Comm() const { return idiscret_.Comm();}
  
  /*!
  \brief Get discretization of the interface

  */
  virtual inline DRT::Discretization& Discret() const { return idiscret_; }
  
  /*!
  \brief Get matrix of DOP normals

  */
  virtual inline const Epetra_SerialDenseMatrix& DopNormals() { return dopnormals_; }
  
  /*!
  \brief Get number of vertices of DOP

  */
  virtual inline const int& Kdop() { return kdop_; }
  
  /*!
  \brief Return dim of Problem

  */
  virtual inline const int& Dim() { return dim_; }
  
  /*!
  \brief Return reference to slave treenodesmap

  */
  virtual inline vector<vector<RCP<BinaryTreeNode> > >& Streenodesmap() { return streenodesmap_; }
  
  /*!
  \brief Return reference to master treenodesmap

  */
  virtual inline vector<vector<RCP<BinaryTreeNode> > >& Mtreenodesmap() { return mtreenodesmap_; }
 
  /*!
  \brief Return reference to contact-treenodesmap

  */
  virtual inline vector<vector<RCP<BinaryTreeNode> > >& ContactMap() { return contactmap_; }
  
  /*!
  \brief Return pointer to sroot-treenode

  */
  virtual inline RCP<BinaryTreeNode>& Sroot() { return sroot_; }
  
  /*!
  \brief Return pointer to mroot-treenode

  */
  virtual inline RCP<BinaryTreeNode>& Mroot() { return mroot_; }
  
  /*!
  \brief Return minimal element length

  */
  virtual inline double & MinLengthEle() { return minlengthele_; }
 
  /*!
  \brief Return factor "eps" to enlarge dops

  */
  virtual inline double & Eps() { return eps_; }
  
  /*!
  \brief Return factor "enlarge" to enlarge dops

  */
  virtual inline double & Enlarge() { return enlarge_; }
    
  //@}
  
  //! @name Evaluation methods
 
  /*!
  \brief Calculate minimal element length / inflation factor "enlarge"

  */
  void SetEnlarge(bool isinit);
  
  /*!
  \brief Update master and slave tree in a topdown way

  */     
  void UpdateTreeTopDown(){EvaluateUpdateTreeTopDown(sroot_);
                           EvaluateUpdateTreeTopDown(mroot_);
                           return;}
   
  /*!
  \brief Evaluate update of master and slave tree in a topdown way

  */     
  void EvaluateUpdateTreeTopDown(RCP<BinaryTreeNode> treenode);
  
  /*!
  \brief Updates master and slave tree in a bottom up way

  */     
  void UpdateTreeBottomUp(){EvaluateUpdateTreeBottomUp(streenodesmap_);
                            EvaluateUpdateTreeBottomUp(mtreenodesmap_);
                            return;}
 
  /*!
  \brief Evaluate update of master and slave tree in a bottom up way

  */
  void EvaluateUpdateTreeBottomUp(vector<vector<RCP<BinaryTreeNode> > >& treenodesmap);

  /*!
  \brief Search for contact between treenodes for initialization of active nodes

  */   
  void SearchContactInit(RCP<BinaryTreeNode> streenode,RCP<BinaryTreeNode> mtreenode);

  /*!
  \brief Search for contact between TreeNodes

  */
  void SearchContactSeparate(){EvaluateSearchContactSeparate(sroot_, mroot_);
                               return;}
    
  /*!
  \brief Evaluate Binary search tree for separate update and search

  */     
  void EvaluateSearchContactSeparate(RCP<BinaryTreeNode> streenode, RCP<BinaryTreeNode> mtreenode);

  /*!
  \brief Search for contact between TreeNodes and automatically update treeEvaluate Binary search tree for combined search and update

  */
  void SearchContactCombined(){// if proc has no slave elements --> return;
                               if (sroot_->Type()==4) return;
                               sroot_->CalculateSlabsDop(false);
                               sroot_->EnlargeGeometry(enlarge_);
                               mroot_->CalculateSlabsDop(false);
                               mroot_->EnlargeGeometry(enlarge_);
                               EvaluateSearchContactCombined(sroot_, mroot_);
                               return;}
  
  /*!
  \brief Evaluate Binary search tree for combined search and update

  */     
  void EvaluateSearchContactCombined(RCP<BinaryTreeNode> streenode, RCP<BinaryTreeNode> mtreenode);
    
  //@}
  
      
protected:
    
  // don't want = operator and cctor
  BinaryTree operator = (const BinaryTree& old);
  BinaryTree(const BinaryTree& old);

  DRT::Discretization&     						   idiscret_;       // contact interface discretization
  RCP<Epetra_Map>           					   selements_;      // all slave elements on surface (column map)
  RCP<Epetra_Map>           					   melements_;      // all master elements on surface (full map)
  vector<vector<RCP<BinaryTreeNode> > >  streenodesmap_;	// map of all slave treenodes, sorted by layers
  vector<vector<RCP<BinaryTreeNode> > >  mtreenodesmap_;	// map of all slave treenodes, sorted by layers
  vector<vector<RCP<BinaryTreeNode> > >  contactmap_;			// map of all treenodes, that are in contact, st/mt
  vector<vector<RCP<BinaryTreeNode> > >  sleafsmap_;			// map of all slave leaf treenodes, [0]=leftchild,[1]=rightchild 
  vector<vector<RCP<BinaryTreeNode> > >  mleafsmap_;			// map of all slave leaf treenodes, [0]=leftchild,[1]=rightlchild
  RCP<BinaryTreeNode>       					   sroot_;          // slave root treenode
  RCP<BinaryTreeNode>       					   mroot_;          // master root treenode
  double															   minlengthele_;		// minimum of length of fe-elements
  const int                 					   dim_;            // problem dimension (2D or 3D)
  double 															   eps_;						// epsilon for enlarging dops (of user)
  double															   enlarge_;				// to enlarge dops
  Epetra_SerialDenseMatrix  					   dopnormals_;     // defining vertices of DOP
  int                       					   kdop_;           // sets k for DOP (8 for 2D, 18 for 3D)   
     
}; // class BinaryTree
}  // namespace CONTACT


#endif  // #ifndef DRT_CONTACT_BINARYTREE_H
#endif  // #ifdef CCADISCRET
