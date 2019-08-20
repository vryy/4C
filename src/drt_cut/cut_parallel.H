/*----------------------------------------------------------------------*/
/*! \file

\brief provides the basic parallel cut classes "Parallel"


\level 1

\maintainer Christoph Ager
 *------------------------------------------------------------------------------------------------*/


#ifndef CUT_PARALLEL_H_
#define CUT_PARALLEL_H_

#include "cut_meshintersection.H"


namespace DRT
{
  class Discretization;
  class PackBuffer;
}  // namespace DRT


namespace GEO
{
  namespace CUT
  {
    class Node;
    class Edge;
    class Side;
    class Element;
    class ElementHandle;
    class Mesh;
    class ParentIntersection;
    //  class MeshIntersection;


    /*!
    \brief this class is the basic TIMEINT class for the projection, adaption or
           something else in XFEM-problems between consecutive time steps
     */
    class Parallel
    {
     public:
      /*!
      \brief Basic CUT parallel constructor
       */
      explicit Parallel(const Teuchos::RCP<DRT::Discretization>& discret, GEO::CUT::Mesh& mesh,
          GEO::CUT::ParentIntersection& parentintersection);

      /*!
      \brief Destructor
      */
      virtual ~Parallel(){};

      /*!
      \brief Communicates the node positions
       */
      void CommunicateNodePositions();

      /*!
      \brief Communicate the node dofset number for single volumecells
       */
      void CommunicateNodeDofSetNumbers(bool include_inner);


     protected:
     private:
      /*!
      \brief Export data whether proc has finished to set node positions
       */
      void exportCommunicationFinished(bool& procDone);

      /*!
      \brief Export position data to neighbor proc and receive data from previous proc
       */
      void exportNodePositionData();

      /*!
      \brief Distribute received node positions on my processor
       */
      void distributeMyReceivedNodePositionData();

      /*!
      \brief set received node positions for node and distribute it to facets, vcs ...
       */
      void setPositionForNode(const Node* n, const Point::PointPosition& pos);

      /*!
      \brief Export dofset data for volumecells to neighbor proc and receive data from previous proc
       */
      void exportDofSetData(bool include_inner);

      /*!
      \brief find the volumecell on myrank for which we received data stored in vc_data
       */
      VolumeCell* findVolumeCell(MeshIntersection::DofSetData&
              vc_data  ///< volumecell data which have to be identified on myrank
      );

      /*!
      \brief distribute received dofset data for volumecells my processor
       */
      void distributeDofSetData();

      /*!
      \brief Pack the point coordinates
       */
      void packPoints(
          DRT::PackBuffer& dataSend, std::vector<LINALG::Matrix<3, 1>>& points_coords) const;

      /*!
      \brief Unpack the point coordinates
       */
      void unpackPoints(std::vector<char>::size_type& posinData,  //!< position in data
          std::vector<char>& dataRecv,                            //!< received data
          std::vector<LINALG::Matrix<3, 1>>& points_coords        //!< point coordinates
          ) const;

      /*!
      \brief Basic function sending data to destination and receiving data from source
       */
      void sendData(DRT::PackBuffer& dataSend,  //!< pack buffer
          int& dest,                            //!< destination proc
          int& source,                          //!< source proc
          std::vector<char>& dataRecv           //!< received data
          ) const;

      /*!
      \brief Print current dofset data
       */
      void printDofSetData();

      /*!
      \brief Get the index of nid in the vector of element's (eid) node Ids
       */
      int getDofSetVecIndex(int nid,  //!< nid for which index has to be found
          int eid                     //!< element id
      );

      //! replace nds (nodal dofset number) vectors for a volumecell set with set_index
      void ReplaceNdsVectors(ElementHandle* e,                   //!< element handle
          const std::vector<plain_volumecell_set>& ele_vc_sets,  //!< sets of volumecells
          std::vector<std::vector<int>>&
              nodaldofset_vc_sets,  //!< corresponding nodal dofsets for sets of volumecells
          int set_index,            //!< number of set
          std::map<int, int>& node_dofsetnumber_map  //!< map of nid and corresponding dofset number
                                                     //!< for all nodes in element
      );

      // data accessing

      //! discretization
      Teuchos::RCP<DRT::Discretization> discret_;

      //! current processor id
      const int myrank_;

      //! number of processors
      const int numproc_;

      GEO::CUT::Mesh& mesh_;  // mesh that carries the cut information of the mesh

      GEO::CUT::ParentIntersection& parentintersection_;

      //! map of node-Id and current node-Pos index
      std::map<int, int> curr_undecidedNodePos_;

      /*!
      \brief map of shadow-nodes (identified via boundary nodes of quad8 sides and the 20 nodes of
      the hex20 element for inner center node) and current node-Pos index
       */
      std::map<plain_int_set, int> curr_undecidedNodePos_shadow_;

      std::vector<Teuchos::RCP<GEO::CUT::MeshIntersection::DofSetData>> dofSetData_;


    };  // class PARALLEL
  }     // namespace CUT
}  // namespace GEO



#endif /* CUT_PARALLEL_H_ */
