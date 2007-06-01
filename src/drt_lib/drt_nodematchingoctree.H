/*!----------------------------------------------------------------------
\file drt_nodematchingoctree.H
\brief connect nodes from two nodesets by their distance


<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef NODEMATCHINGOCTREE_H
#define NODEMATCHINGOCTREE_H

#include "drt_utils.H"
#include "drt_discret.H"
#include "drt_parobject.H"
#include "drt_exporter.H"
#include "Teuchos_RefCountPtr.hpp"

extern "C"
{
#include "../headers/standardtypes.h"
}
class OctreeElement;


class NodeMatchingOctree
{
public:
  

  /*!
   
  \brief Standard Constructor (public)
  
  <pre> 
  
  Set up processor local octree                               gammi 04/07


  On input a list of so called "master"-nodes and a discretisation which
  contains these "master"-nodes. All nodeids from masternodes owned by
  this proc are put into a octree structure, i.e. are attached to leafs
  of an octree.

  The number maxnodeperleaf specifies the upper bound for the number of
  nodeids attached to a octree leaf. If the number would be higher, the
  octree element would be split into two children (which could be leafs
  on their own)
  
  If a nodeid is associated with a child of an octree element, it is
  checked whether its coordinate is inside the bounding box of the
  octree element. This is done only with a tolerance tol to ensure that
  each node is guaranteed to be owned by one of the children an could
  be found later on.
  
  </pre>
  
  \param   actdis          (i) discretisation
  \param   masternodeids   (i) list of masternodeids
  \param   maxnodeperleaf  (i) parameter for octree
  \param   tol             (i) tolerance for octree
  
  \return void
  
  
  */
  NodeMatchingOctree(
    RefCountPtr<DRT::Discretization> actdis,
    const vector <int> &             masternodeids,
    int                              maxnodeperleaf,
    double                           tol
    );

  /*!
  \brief Destructor

  */
  virtual ~NodeMatchingOctree();



  //! @name Methods to create global node matching

  /*!
   
  \brief Search for closest (slave) nodes on all procs
  
  <pre> 

  Search for closest (slave) nodes on all processors to given (master)
  nodeset (only in bounding box of masternodes)               gammi 05/07

   1) each proc generates a list of his slavenodes (out of slavenodeids)

   2) the list is communicated in a round robin pattern to all the
      other procs.
  
   3) the proc checks the package from each proc and calcs the min
      distance on each --- the result is kept if distance is smaller
      than on the preceding processors
  
  Again, dofsforpbcplane is used to overwrite the coordinate of the
  slave perpendicular to the plane of the periodic boundary condition
  by the coordinate of the masterplane (so to speak, the closest node
  matching is "exact")

  On output, we get a map from masterids on this proc to slavenodeids
  (on arbitrary procs).
  
  </pre>
  
  \param  slavenodeids     (i) list of slavenodeids
  \param  dofsforpbcplane  (i) specification of the plane containing
                               parallel to the boundary condition
                               (used to "shift" nodes in the normal
                                direction) 
  \param  midtosid         (o) map from master to slavenodes
  
  \return void
  
  */  
  void CreateGlobalNodeMatching(
    const vector <int> slavenodeids,
    vector <int>       dofsforpbcplane,
    map<int,int> &     midtosid
    );

  //@}

  
protected:

   //! @name Methods to create processor local node matching
 
  /*!
   
  \brief 
  
  <pre> 

  Search closest node to given coordinate x in local octree.
  
  returns false if node is not in bounding box of local octree

  gammi 05/07
  
  </pre>
  
  \param  x                  (i) coordinate of slavenode
  \param  idofclosestpoint   (o) the node id of the closest point
  \param  distofclosestpoint (o) the distance to the closest point
  
  \return  bool        false if node is not in bounding box of
                       local octree
  
  */
  bool SearchClosestNodeOnThisProc(
    vector <double> x,
    int           & idofclosestpoint,
    double        & distofclosestpoint
  );
  
  //@}
  
private:

  //! \brief all nodes in leafs are nodes of this discretistion
  RefCountPtr<DRT::Discretization> discret_;
  //! \brief root of the local search tree
  RefCountPtr<OctreeElement>       octreeroot_;
  //! \brief coordinate of one point in the master plane
  vector <double>                  masterplanecoords_;
}; // class NodeMatchingOctree





class OctreeElement
{

public:
  

  /*!
  \brief Standard Constructor

  <pre>
  
  Create one element in octree                              gammi 05/07

  nodeids is a list of nodes out of discretisation actdis belonging to
  the element or its children. bounding box specifies the geometry of
  the element. layer is the number of subdivisions after the octree
  root. maxnodeperleaf specifies the max number of nodes contained by
  one leaf element. tol is a tolerance used to assign the nodes to the
  elements by checking if they are in the bounding box or not.
  
  </pre> 

  \param  actdis          (i) the discretisation
  \param  nodeids         (i) nodes attached to element
  \param  boundingbox     (i) bounding box of element
  \param  layer           (i) depth in tree
  \param  maxnodeperleaf  (i) how many nodes does
                              one leaf contain?
  \param  tol             (i) Tolerance for octree

  \return  void
  
  */

  OctreeElement(
    RefCountPtr<DRT::Discretization> actdis,
    vector <int> &                   nodeids,
    Epetra_SerialDenseMatrix&        boundingbox,
    int                              layer,
    int                              maxnodeperleaf,
    double                           tol
    );

  /*!
  \brief Destructor

  */
  virtual ~OctreeElement();


  //! @name Octree element functions

  /*!
  \brief Question if coordinate is in bounding box

  <pre>
  
  Is x in the bounding box of the element?              gammi 05/07
  
  </pre> 

  \param x         (i) coordinate of point

  \return  bool true if point in bounding box
  
  */
  bool IsPointInBoundingBox(
    vector <double> &x
    );

  /*!
  \brief Question if octree element is leaf

  <pre>
  
  is the octree element a leaf, i.e. does it posess a vector
  of nodes attached to it?                         gammi 05/07
  
  </pre> 

  \return  bool true if element is leaf
  
  */
  bool IsLeaf();

  /*!
  \brief Return child containing point

  <pre>
  
  Determine which child contains the coordinate x. In some cases
  the return value may not be the unique answer (overlap!)

                                                   gammi 05/07
  
  </pre>

  \param   x   (i) coordinate

  \return  RefCountPtr<OctreeElement> child
  
  */
  RefCountPtr<OctreeElement> ReturnChildContainingPoint(
    vector <double> &x
    );

  /*!
  \brief Find closest point in leaf

  <pre>

  leaf contains still maxnodeperleaf nodes. This function
  determines the closest point of all nodes attached to this leaf.

                                                   gammi 05/07
  
  </pre>

  \param   x                  (i) coordinate of point
  \param   idofclosestpoint   (o) global id of closest point
  \param   distofclosestpoint (o) distance of closest point

  \return void
  
  */
  
  void SearchClosestNodeInLeaf(
    vector <double> & x,
    int             & idofclosestpoint,
    double          & distofclosestpoint
  );


  /*!
    Print some information on the octree leaf
  
  <pre>

  cout the node coordinates associated to a leaf
                                                   gammi 05/07
  
  </pre>

  \param  os (i) 

  \return void
  
  */
  void Print(ostream& os) const;
  
  //@}  
protected:

private:
  
  //! \brief all nodes belong to this discretisation
  RefCountPtr<DRT::Discretization> discret_;
  //! \brief the bounding box of the element
  Epetra_SerialDenseMatrix         boundingbox_;
  //! \brief the nodeids
  vector <int>                     nodeids_;
  //! \brief depth of the octree element in the tree
  int                              layer_;
  //! \brief  pointer to first child
  RefCountPtr<OctreeElement>       octreechild1_;
  //! \brief  pointer to second child
  RefCountPtr<OctreeElement>       octreechild2_;
  
}; // class OctreeElement


#endif  // not NODEMATCHINGOCTREE_H

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
