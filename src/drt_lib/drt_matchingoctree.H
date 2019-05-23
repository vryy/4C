/*---------------------------------------------------------------------*/
/*!

\brief connect nodes from two nodesets by their distance

\level 1

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef MATCHINGOCTREE_H
#define MATCHINGOCTREE_H

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include <vector>
#include "../drt_lib/drt_dserror.H"

namespace DRT
{
  class Discretization;
  class Element;
  class Node;
  class PackBuffer;
  class ParObject;

  namespace UTILS
  {
    class OctreeElement;

    //! Base class for parallel octree to establish neighborhood relations between entities
    class MatchingOctree
    {
     public:
      //! Constructor
      MatchingOctree();

      //! Destructor
      virtual ~MatchingOctree();

      /*! \brief initialize this class (public)


      Initialize processor local octree                               rauch 09/16

      On input a list of so called "master"-entities and a discretisation which
      contains these "master"-entities. Make sure, that you provide the octree
      only with those entity gids, you want in the processor local octree.
      You may want to match only entities, which are owned, ghosted, or both.

      The number maxeleperleaf specifies the upper bound for the number of
      eleids attached to a octree leaf. If the number would be higher, the
      octree element would be split into two children (which could be leafs
      on their own)

      If an entity id is associated with a child of an octree element, it is
      checked whether its coordinate is inside the bounding box of the
      octree element. This is done only with a tolerance tol to ensure that
      each node is guaranteed to be owned by one of the children and could
      be found later on. The same tolerance is used for finding the smallest
      distance between two nodes. Thus, it should be of the order or magnitude
      of the smallest element size.


      \param   actdis         (i) discretisation
      \param   mastereleids   (i) list of master entity ids
      \param   maxeleperleaf  (i) parameter for octree
      \param   tol            (i) tolerance for octree

      \return void  */
      virtual int Init(const DRT::Discretization& actdis, const std::vector<int>& masternodeids,
          const int maxnodeperleaf = 150, const double tol = 1e-08);

      //! setup this class
      virtual int Setup();

     protected:
      //! @name Methods to create processor local node matching

      /*! \brief Declaration


      Search closest entity to given coordinate x in local octree.

      returns false if entity is not in bounding box of local octree

      gammi 05/07


      \param  x                  (i) coordinate of slavenode
      \param  idofclosestpoint   (o) the node id of the closest point
      \param  distofclosestpoint (o) the distance to the closest point

      \return  bool        false if node is not in bounding box of
                           local octree  */
      virtual bool SearchClosestEntityOnThisProc(const std::vector<double>& x,
          int& idofclosestpoint, double& distofclosestpoint, bool searchsecond = false);

      //@}


      //! @name Pure virtual methods which have to be overridden by specialized octrees

      //! calc unique coordinate of entity
      virtual void CalcPointCoordinate(
          const DRT::Discretization* dis, const int id, double* coord) = 0;

      //! calc unique coordinate of entity
      virtual void CalcPointCoordinate(DRT::ParObject* entity, double* coord) = 0;

      //! check if entity is on calling proc
      virtual bool CheckHaveEntity(const DRT::Discretization* dis, const int id) = 0;

      //! returns true if entity is owned by calling proc
      virtual bool CheckEntityOwner(const DRT::Discretization* dis, const int id) = 0;

      //! pack entity to PackPuffer
      virtual void PackEntity(PackBuffer& data, const DRT::Discretization* dis, const int id) = 0;

      //! unpack entity to PackPuffer
      virtual void UnPackEntity(std::vector<char>::size_type& index,
          std::vector<char>& rblockofnodes, std::vector<char>& data) = 0;

      //! check if unpacked type is correct
      virtual int CheckValidEntityType(Teuchos::RCP<DRT::ParObject> o) = 0;

      //! create an octree element
      virtual Teuchos::RCP<OctreeElement> CreateOctreeElement(std::vector<int>& nodeidstoadd,
          Epetra_SerialDenseMatrix& boundingboxtoadd, int layer) = 0;

      //@}

     public:
      //! @name Methods to create global node matching

      /*! \brief Search for closest (slave) nodes on all procs


      Search for closest (slave) entities on all processors to given (master)
      entity set (only in bounding box of master entities)        gammi 05/07

       1) each proc gets a list of his slave entities (slavenodeids)

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


      \param  slavenodeids     (i) list of slavenodeids
      \param  dofsforpbcplane  (i) specification of the plane containing
                                   parallel to the boundary condition
                                   (used to "shift" nodes in the normal
                                    direction)
      \param  rotangle         (i) angle (RAD) for rotation of slave plane
      \param  midtosid         (o) map from master to slavenodes

      \return void  */
      virtual void CreateGlobalEntityMatching(const std::vector<int>& slavenodeids,
          const std::vector<int>& dofsforpbcplane, const double rotangle,
          std::map<int, std::vector<int>>& midtosid);

      /*! \brief find pairs of nearest nodes

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
        \param coupling     (o) master node gid to (slave node gid, distance) */
      virtual void FindMatch(const DRT::Discretization& slavedis,
          const std::vector<int>& slavenodeids, std::map<int, std::pair<int, double>>& coupling);

      /*! \brief find pairs of nearest nodes

        Basically, the algorithm performed in this method is equal to the one
        described in the documentation to \ref FindMatch .

        \note On exit, the map 'coupling' maps the slave entity gids (key) to
              the matching master entity gids, distance, and the information,
              whether the current proc owns the matching master entity or not.
              This is different from the method \ref FindMatch.

        \param slavedis     (i) discretization the slave entity belongs to
        \param slavenodeids (i) gids of entity to match
        \param coupling     (o) slave entity gid to (slave entity gid, distance, owned) */
      virtual void FillSlaveToMasterGIDMapping(const DRT::Discretization& slavedis,
          const std::vector<int>& slavenodeids, std::map<int, std::vector<double>>& coupling);

      //@}

     protected:
      //! \brief all nodes in leafs are nodes of this discretization
      const DRT::Discretization* discret_;
      //! \brief order of magnitude of smallest element size (used for tolerances)
      double tol_;
      //! \brief root of the local search tree
      Teuchos::RCP<OctreeElement> octreeroot_;
      //! \brief coordinate of one point in the master plane
      std::vector<double> masterplanecoords_;
      //! \brief ids of entities to be coupled (e.g. nodes in \ref NodeMatchingOctree )
      const std::vector<int>* masterentityids_;
      //! maximum number of tree nodes per leaf
      int maxtreenodesperleaf_;

     private:
      //! flag indicating if class is setup
      bool issetup_;

      //! flag indicating if class is initialized
      bool isinit_;

     protected:
      //! returns true if Setup() was called and is still valid
      bool IsSetup() { return issetup_; };

      //! returns true if Init(..) was called and is still valid
      bool IsInit() { return isinit_; };

      //! check if \ref Setup() was called
      void CheckIsSetup()
      {
        if (not IsSetup()) dserror("Setup() was not called.");
      };

      //! check if \ref Init() was called
      void CheckIsInit()
      {
        if (not IsInit()) dserror("Init(...) was not called.");
      };

     private:
      //! set flag true after setup or false if setup became invalid
      void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

      //! set flag true after init or false if init became invalid
      void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };

    };  // MatchingOctree


    //! Parallel octree to establish neighborhood relations between two sets of nodes
    class NodeMatchingOctree : public MatchingOctree
    {
     public:
      //! Constructor
      NodeMatchingOctree();

      //! Destructor
      virtual ~NodeMatchingOctree();

     protected:
      //! calc unique coordinate, associated to node
      virtual void CalcPointCoordinate(const DRT::Discretization* dis, const int id, double* coord);

      //! calc unique coordinate of Node
      virtual void CalcPointCoordinate(DRT::ParObject* entity, double* coord);

      //! check if node with gid = id is on calling proc
      virtual bool CheckHaveEntity(const DRT::Discretization* dis, const int id);

      //! returns true if node with gid = id is owned by calling proc
      virtual bool CheckEntityOwner(const DRT::Discretization* dis, const int id);

      //! pack node to PackPuffer
      virtual void PackEntity(PackBuffer& data, const DRT::Discretization* dis, const int id);

      //! unpack node from PackPuffer
      virtual void UnPackEntity(std::vector<char>::size_type& index,
          std::vector<char>& rblockofnodes, std::vector<char>& data);

      //! check if unpacked type is correct
      virtual int CheckValidEntityType(Teuchos::RCP<DRT::ParObject> o);

      //! create an octree element
      virtual Teuchos::RCP<OctreeElement> CreateOctreeElement(
          std::vector<int>& nodeidstoadd, Epetra_SerialDenseMatrix& boundingboxtoadd, int layer);

    };  // class NodeMatchingOctree


    //! Parallel octree to establish neighborhood relations between two sets of elements
    class ElementMatchingOctree : public MatchingOctree
    {
     public:
      //! Constructor
      ElementMatchingOctree();

      //! Destructor
      virtual ~ElementMatchingOctree();

     protected:
      //! calc unique coordinate, associated to element
      virtual void CalcPointCoordinate(const DRT::Discretization* dis, const int id, double* coord);

      //! calc unique coordinate of entity
      virtual void CalcPointCoordinate(DRT::ParObject* entity, double* coord);

      //! check if element with gid = id is on calling proc
      virtual bool CheckHaveEntity(const DRT::Discretization* dis, const int id);

      //! returns true if element is owned by calling proc
      virtual bool CheckEntityOwner(const DRT::Discretization* dis, const int id);

      //! pack element to PackPuffer
      virtual void PackEntity(PackBuffer& data, const DRT::Discretization* dis, const int id);

      //! unpack element from PackPuffer
      virtual void UnPackEntity(std::vector<char>::size_type& index,
          std::vector<char>& rblockofnodes, std::vector<char>& data);

      //! check if unpacked type is correct
      virtual int CheckValidEntityType(Teuchos::RCP<DRT::ParObject> o);

      //! create an octree element
      virtual Teuchos::RCP<OctreeElement> CreateOctreeElement(
          std::vector<int>& nodeidstoadd, Epetra_SerialDenseMatrix& boundingboxtoadd, int layer);

     private:
      //! \brief Map that stores id and RCP to flying nodes
      //!
      //! The flying nodes are sent around in the round robin loop
      //! together with the corresponding element.
      //!
      //! \note We need this map as member, because in \ref UnPackEntity
      //!       we Extract the nodes and fill this map. Then, we set the
      //!       nodal pointers in the communicated element by providing
      //!       the element with this map. The RCPs stored in this map have
      //!       to survive until they re used in \ref CalcPointCoordinate .
      //!       Therefore, we need to have a global temporary storage for
      //!       the nodal pointers, which is provided by this member.
      //!
      //! \author Andreas Rauch
      //! \date   10/16
      std::map<int, Teuchos::RCP<DRT::Node>> nodes_;

    };  // class ElementMatchingOctree


    //! Leaf in parallel octree
    class OctreeElement
    {
     public:
      //! Standard Constructor
      OctreeElement();

      //! Destructor
      virtual ~OctreeElement();

      /*! \brief initialize this class


      Initialize one element in octree                              rauch 09/16

      nodeids is a list of tree-nodes out of discretisation actdis belonging to
      the tree-element or its children. bounding box specifies the geometry of
      the tree-element. layer is the number of subdivisions after the octree
      root. maxnodeperleaf specifies the max number of nodes contained by
      one leaf element. tol is a tolerance used to assign the tree-nodes to the
      tree-elements by checking if they are in the bounding box or not.


      \param  actdis          (i) the discretisation
      \param  nodeids         (i) tree-nodes attached to tree-element
      \param  boundingbox     (i) bounding box of element
      \param  layer           (i) depth in tree
      \param  maxnodeperleaf  (i) how many nodes does
                                  one leaf contain?
      \param  tol             (i) Tolerance for octree

      \return  int  */
      virtual int Init(const DRT::Discretization& actdis, std::vector<int>& nodeids,
          const Epetra_SerialDenseMatrix& boundingbox, const int layer, const int maxnodeperleaf,
          const double tol);

      //! setup this class
      virtual int Setup();

      /*! \brief Find closest point in leaf


      leaf contains still maxnodeperleaf nodes. This function
      determines the closest point of all nodes attached to this leaf.

                                                       gammi 05/07


      \param   x                  (i) coordinate of point
      \param   idofclosestpoint   (o) global id of closest point
      \param   distofclosestpoint (o) distance of closest point
      \param   elesize            (i) tolerance for node matching
      \param   searchsecond       (i) flag for search of second match

      \return void  */
      virtual void SearchClosestNodeInLeaf(const std::vector<double>& x, int& idofclosestpoint,
          double& distofclosestpoint, const double& elesize, bool searchsecond);

      //! @name Octree element functions

      /*! \brief Question if coordinate is in bounding box


      Is x in the bounding box of the element?              gammi 05/07


      \param x         (i) coordinate of point

      \return  bool true if point in bounding box  */
      bool IsPointInBoundingBox(const std::vector<double>& x);

      /*! \brief Question if octree element is leaf


      is the octree element a leaf, i.e. does it posess a vector
      of nodes attached to it?                         gammi 05/07


      \return  bool true if element is leaf  */
      bool IsLeaf();

      /*! \brief Return child containing point


      Determine which child contains the coordinate x. In some cases
      the return value may not be the unique answer (overlap!)

                                                       gammi 05/07


      \param   x   (i) coordinate

      \return  Teuchos::RCP<OctreeElement> child  */
      Teuchos::RCP<OctreeElement> ReturnChildContainingPoint(const std::vector<double>& x);

      /*! \brief Print some information on the octree leaf


      std::cout the node coordinates associated to a leaf
                                                       gammi 05/07


      \param  os (i)

      \return void  */
      void Print(std::ostream& os) const;

      //@}

     protected:
      //! calc unique coordinate of entity
      virtual void CalcPointCoordinate(
          const DRT::Discretization* dis, const int id, double* coord) = 0;

      //! create an octree element
      virtual Teuchos::RCP<OctreeElement> CreateOctreeElement(std::vector<int>& nodeidstoadd,
          Epetra_SerialDenseMatrix& boundingboxtoadd, int layer) = 0;

     protected:
      //! all nodes belong to this discretisation
      const DRT::Discretization* discret_;
      //! the bounding box of the element
      Epetra_SerialDenseMatrix boundingbox_;
      //! the nodeids (needs to be a copy !)
      std::vector<int> nodeids_;
      //! \brief depth of the octree element in the tree
      int layer_;
      //! maximum number of tree nodes per leaf
      int maxtreenodesperleaf_;
      //! tolerance for octree
      double tol_;
      //! pointer to first child
      Teuchos::RCP<OctreeElement> octreechild1_;
      //! pointer to second child
      Teuchos::RCP<OctreeElement> octreechild2_;

     private:
      //! flag indicating if class is setup
      bool issetup_;

      //! flag indicating if class is initialized
      bool isinit_;

     protected:
      //! returns true if Setup() was called and is still valid
      bool IsSetup() { return issetup_; };

      //! returns true if Init(..) was called and is still valid
      bool IsInit() { return isinit_; };

      //! check if \ref Setup() was called
      void CheckIsSetup()
      {
        if (not IsSetup()) dserror("Setup() was not called.");
      };

      //! check if \ref Init() was called
      void CheckIsInit()
      {
        if (not IsInit()) dserror("Init(...) was not called.");
      };

     private:
      //! set flag true after setup or false if setup became invalid
      void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

      //! set flag true after init or false if init became invalid
      void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };

    };  // class OctreeElement


    //! Leaf in parallel octree for FEM nodes
    class OctreeNodalElement : public OctreeElement
    {
     public:
      //! Constructor
      OctreeNodalElement();

      //! Destructor
      virtual ~OctreeNodalElement() { return; };

     protected:
      //! calc unique coordinate of node
      virtual void CalcPointCoordinate(const DRT::Discretization* dis, const int id, double* coord);

      //! create an octree element
      virtual Teuchos::RCP<OctreeElement> CreateOctreeElement(
          std::vector<int>& nodeidstoadd, Epetra_SerialDenseMatrix& boundingboxtoadd, int layer);

    };  // OctreeNodalElement


    //! Leaf in parallel octree for FEM elements
    class OctreeElementElement : public OctreeElement
    {
     public:
      //! Constructor
      OctreeElementElement();

      //! Destructor
      virtual ~OctreeElementElement() { return; };

     protected:
      //! calc unique coordinate of element
      virtual void CalcPointCoordinate(const DRT::Discretization* dis, const int id, double* coord);

      //! create an octree element
      virtual Teuchos::RCP<OctreeElement> CreateOctreeElement(
          std::vector<int>& nodeidstoadd, Epetra_SerialDenseMatrix& boundingboxtoadd, int layer);

    };  // OctreeNodalElement

  }  // namespace UTILS
}  // namespace DRT

#endif  // not MATCHINGOCTREE_H
