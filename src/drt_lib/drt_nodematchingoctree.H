/*!----------------------------------------------------------------------
\file drt_nodematchingoctree.H

\brief connect nodes from two nodesets by their distance

<pre>
\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef NODEMATCHINGOCTREE_H
#define NODEMATCHINGOCTREE_H

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include <vector>

namespace DRT {

class Discretization;

namespace UTILS {

class OctreeElement;

/// Parallel octree to establish neighborhood relations between two sets of nodes
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
  be found later on. The same tolerance is used for finding the smallest
  distance between two nodes. Thus, it should be of the order or magnitude
  of the smallest element size.

  </pre>

  \param   actdis          (i) discretisation
  \param   masternodeids   (i) list of masternodeids
  \param   maxnodeperleaf  (i) parameter for octree
  \param   tol             (i) tolerance for octree

  \return void


  */
  NodeMatchingOctree(
    const DRT::Discretization&       actdis,
    const std::vector <int> &        masternodeids,
    int                              maxnodeperleaf=150,
    double                           tol=1e-8
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
  \param  rotangle         (i) angle (RAD) for rotation of slave plane
  \param  midtosid         (o) map from master to slavenodes

  \return void

  */
  void CreateGlobalNodeMatching(
    const std::vector<int>  & slavenodeids,
    const std::vector<int>  & dofsforpbcplane,
    const double            rotangle,
    std::map<int,std::vector<int> > & midtosid
    );

  /// find pairs of nearest nodes
  /*!
    Find the pairs of master node and slave node with the smallest
    distance. The sets of nodes might come from different
    discretisations, however the communicators of both discretisations
    need to span the same processors. For communication the master
    discretisation communicator is used.

    Internally the search is based on the octtree spanned by the
    master nodes. The slave nodes are gathered on their owning
    processors and send in a round robin fashion to each processor
    once.

    \param slavedis     (i) discretization the slave nodes belong to
    \param slavenodeids (i) gids of nodes to match
    \param coupling     (o) master node gid to (slave node gid, distance)

   */
  void FindMatch(const DRT::Discretization& slavedis,
                 const std::vector<int>& slavenodeids,
                 std::map<int,std::pair<int,double> >& coupling);

  //@}


protected:

   //! @name Methods to create processor local node matching

  /*!

  \brief Declaration

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
    const std::vector<double>& x,
    int                      & idofclosestpoint,
    double                   & distofclosestpoint,
    bool                       searchsecond = false
  );

  //@}

private:

  //! \brief all nodes in leafs are nodes of this discretization
  const DRT::Discretization&       discret_;
  //! \brief order of magnitude of smallest element size (used for tolerances)
  const double                     tol_;
  //! \brief root of the local search tree
  Teuchos::RCP<OctreeElement>      octreeroot_;
  //! \brief coordinate of one point in the master plane
  std::vector <double>             masterplanecoords_;
}; // class NodeMatchingOctree



/// Leaf in parallel octree
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
    const DRT::Discretization&       actdis,
    std::vector<int> &               nodeids,
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
    const std::vector <double> &x
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

  \return  Teuchos::RCP<OctreeElement> child

  */
  Teuchos::RCP<OctreeElement> ReturnChildContainingPoint(
    const std::vector<double> &x
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
  \param   elesize            (i) tolerance for node matching
  \param   x                  (i) flag for search of second match

  \return void

  */

  void SearchClosestNodeInLeaf(
    const std::vector<double> & x,
    int                       & idofclosestpoint,
    double                    & distofclosestpoint,
    const double              & elesize,
    bool                        searchsecond
  );


  /*!
    Print some information on the octree leaf

  <pre>

  std::cout the node coordinates associated to a leaf
                                                   gammi 05/07

  </pre>

  \param  os (i)

  \return void

  */
  void Print(std::ostream& os) const;

  //@}
protected:

private:

  //! \brief all nodes belong to this discretisation
  const DRT::Discretization&       discret_;
  //! \brief the bounding box of the element
  Epetra_SerialDenseMatrix         boundingbox_;
  //! \brief the nodeids
  std::vector <int>                nodeids_;
  //! \brief depth of the octree element in the tree
  int                              layer_;
  //! \brief  pointer to first child
  Teuchos::RCP<OctreeElement>       octreechild1_;
  //! \brief  pointer to second child
  Teuchos::RCP<OctreeElement>       octreechild2_;

}; // class OctreeElement

}
}

#endif  // not NODEMATCHINGOCTREE_H
