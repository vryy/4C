/*----------------------------------------------------------------------*/
/*! \file

\brief provides the basic parallel cut classes "Parallel"


\level 1

 *------------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_CUT_PARALLEL_HPP
#define FOUR_C_CUT_PARALLEL_HPP

#include "4C_config.hpp"

#include "4C_cut_meshintersection.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{
  class PackBuffer;
}

namespace Discret
{
  class Discretization;
}  // namespace Discret


namespace Core::Geo
{
  namespace Cut
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
      explicit Parallel(const Teuchos::RCP<Discret::Discretization>& discret,
          Core::Geo::Cut::Mesh& mesh, Core::Geo::Cut::ParentIntersection& parentintersection);

      /*!
      \brief Destructor
      */
      virtual ~Parallel() = default;

      /*!
      \brief Communicates the node positions
       */
      void communicate_node_positions();

      /*!
      \brief Communicate the node dofset number for single volumecells
       */
      void communicate_node_dof_set_numbers(bool include_inner);


     protected:
     private:
      /*!
      \brief Export data whether proc has finished to set node positions
       */
      void export_communication_finished(bool& procDone);

      /*!
      \brief Export position data to neighbor proc and receive data from previous proc
       */
      void export_node_position_data();

      /*!
      \brief Distribute received node positions on my processor
       */
      void distribute_my_received_node_position_data();

      /*!
      \brief set received node positions for node and distribute it to facets, vcs ...
       */
      void set_position_for_node(const Node* n, const Point::PointPosition& pos);

      /*!
      \brief Export dofset data for volumecells to neighbor proc and receive data from previous proc
       */
      void export_dof_set_data(bool include_inner);

      /*!
      \brief find the volumecell on myrank for which we received data stored in vc_data
       */
      VolumeCell* find_volume_cell(
          MeshIntersection::DofSetData&
              vc_data,  ///< volumecell data which have to be identified on myrank
          double tol    ///< geometric tolerance to fit the points of volumecells on two different
                        ///< processors
      );

      /*!
      \brief distribute received dofset data for volumecells my processor
       */
      void distribute_dof_set_data();

      /*!
      \brief Pack the point coordinates
       */
      void pack_points(Core::Communication::PackBuffer& dataSend,
          std::vector<Core::LinAlg::Matrix<3, 1>>& points_coords) const;

      /*!
      \brief Unpack the point coordinates
       */
      void unpack_points(std::vector<char>::size_type& posinData,  //!< position in data
          std::vector<char>& dataRecv,                             //!< received data
          std::vector<Core::LinAlg::Matrix<3, 1>>& points_coords   //!< point coordinates
      ) const;

      /*!
      \brief Basic function sending data to destination and receiving data from source
       */
      void send_data(Core::Communication::PackBuffer& dataSend,  //!< pack buffer
          int& dest,                                             //!< destination proc
          int& source,                                           //!< source proc
          std::vector<char>& dataRecv                            //!< received data
      ) const;

      /*!
      \brief Print current dofset data
       */
      void print_dof_set_data();

      /*!
      \brief Get the index of nid in the vector of element's (eid) node Ids
       */
      int get_dof_set_vec_index(int nid,  //!< nid for which index has to be found
          int eid                         //!< element id
      );

      //! replace nds (nodal dofset number) vectors for a volumecell set with set_index
      void replace_nds_vectors(ElementHandle* e,                 //!< element handle
          const std::vector<plain_volumecell_set>& ele_vc_sets,  //!< sets of volumecells
          std::vector<std::vector<int>>&
              nodaldofset_vc_sets,  //!< corresponding nodal dofsets for sets of volumecells
          int set_index,            //!< number of set
          std::map<int, int>& node_dofsetnumber_map  //!< map of nid and corresponding dofset number
                                                     //!< for all nodes in element
      );

      // data accessing

      //! discretization
      Teuchos::RCP<Discret::Discretization> discret_;

      //! current processor id
      const int myrank_;

      //! number of processors
      const int numproc_;

      Core::Geo::Cut::Mesh& mesh_;  // mesh that carries the cut information of the mesh

      Core::Geo::Cut::ParentIntersection& parentintersection_;

      //! map of node-Id and current node-Pos index
      std::map<int, int> curr_undecided_node_pos_;

      /*!
      \brief map of shadow-nodes (identified via boundary nodes of quad8 sides and the 20 nodes of
      the hex20 element for inner center node) and current node-Pos index
       */
      std::map<plain_int_set, int> curr_undecided_node_pos_shadow_;

      std::vector<Teuchos::RCP<Core::Geo::Cut::MeshIntersection::DofSetData>> dof_set_data_;


    };  // class PARALLEL
  }     // namespace Cut
}  // namespace Core::Geo



FOUR_C_NAMESPACE_CLOSE

#endif
