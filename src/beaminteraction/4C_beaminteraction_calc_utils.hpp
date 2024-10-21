// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_CALC_UTILS_HPP
#define FOUR_C_BEAMINTERACTION_CALC_UTILS_HPP

#include "4C_config.hpp"

#include "4C_binstrategy_utils.hpp"
#include "4C_inpar_beaminteraction.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_mapextractor.hpp"

#include <Epetra_FEVector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class SerialDenseVector;
  class SerialDenseMatrix;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace Core::Nodes
{
  class Node;
}


namespace Core::Geo::MeshFree
{
  class BoundingBox;
}


namespace BEAMINTERACTION
{
  class CrosslinkingParams;
  class BeamLink;
  namespace Utils
  {
    /// specific MultiMapExtractor to handle different types of element during beam interaction
    class MapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        beam = 0,
        sphere = 1,
        solid = 2
      };

      MAP_EXTRACTOR_VECTOR_METHODS(beam, beam)
      MAP_EXTRACTOR_VECTOR_METHODS(sphere, sphere)
      MAP_EXTRACTOR_VECTOR_METHODS(solid, solid)
    };

    /// class for comparing Core::Elements::Element* (and Core::Nodes::Node*) in a std::set
    /*! -------------------------------------------------------------------------
     * overwrites standard < for pointers, this is necessary to ensure same order
     * of neighboring elements for crosslink check and therefore for random numbers
     * independent of pointer addresses. Without this,
     * simulation with crosslinker is not wrong, but depends randomly on memory
     * allocation, i.e. pointer addresses. Without random numbers, everything is fine
     * with default compare operator
    *  \author J. Eichinger March 2017
     -------------------------------------------------------------------------*/
    class Less
    {
     public:
      template <typename ELEMENT>
      bool operator()(ELEMENT const* first, ELEMENT const* second) const
      {
        return first->id() < second->id();
      }
    };

    /*! -------------------------------------------------------------------------
     * class for comparing std::set< std::pair < int, int > >
     *  \author J. Eichinger March 2017
     -------------------------------------------------------------------------*/
    class StdPairComparatorOrderCounts
    {
     public:
      bool operator()(std::pair<int, int> const& lhs, std::pair<int, int> const& rhs) const
      {
        return (lhs.first == rhs.first) ? (lhs.second < rhs.second) : (lhs.first < rhs.first);
      }
    };


    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool is_beam_element(Core::Elements::Element const& element);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool is_rigid_sphere_element(Core::Elements::Element const& element);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool is_beam_node(Core::Nodes::Node const& node);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool is_rigid_sphere_node(Core::Nodes::Node const& node);

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    bool is_beam_centerline_node(Core::Nodes::Node const& node);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void periodic_boundary_consistent_dis_vector(Core::LinAlg::Vector<double>& dis,
        const Core::Geo::MeshFree::BoundingBox& pbb, const Core::FE::Discretization& discret);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    inline int calculate_number_of_beam_elements_from_number_of_nodes_on_filament(
        int const numnodes, int const numnodesperele)
    {
      // from: nodesperfil = nodesperele + ( numele - 1 ) * ( nodesperele - 1 )
      return ((numnodes - numnodesperele) / (numnodesperele - 1.0)) + 1.0;
    }

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    std::vector<int> permutation(int number);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void get_current_element_dis(Core::FE::Discretization const& discret,
        Core::Elements::Element const* ele, const Core::LinAlg::Vector<double>& ia_discolnp,
        std::vector<double>& eledisp);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void get_current_unshifted_element_dis(Core::FE::Discretization const& discret,
        Core::Elements::Element const* ele, const Core::LinAlg::Vector<double>& ia_discolnp,
        Core::Geo::MeshFree::BoundingBox const& pbb, std::vector<double>& eledisp);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    template <typename T>
    void set_filament_binding_spot_positions(
        Teuchos::RCP<Core::FE::Discretization> discret, T& params);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void extend_ghosting_for_filament_bspot_setup(
        std::set<int>& relevantfilaments, Core::FE::Discretization& discret);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void determine_off_my_rank_nodes_with_relevant_ele_cloud_for_filament_bspot_setup(
        std::set<int>& relevantfilaments, std::set<int>& setofrequirednodes,
        Core::FE::Discretization& discret);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void compute_filament_length_and_sort_its_elements(
        std::vector<Core::Elements::Element*>& sortedfilamenteles, std::vector<int> const* nodeids,
        double& filreflength, Core::FE::Discretization& discret);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void set_binding_spots_positions_on_filament(
        std::vector<Core::Elements::Element*>& sortedfilamenteles, double start,
        Inpar::BEAMINTERACTION::CrosslinkerType linkertype, int numbspot,
        double filamentbspotinterval, double tol);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void get_pos_and_triad_of_binding_spot(Core::Elements::Element* ele,
        Core::Geo::MeshFree::BoundingBox& pbb, Inpar::BEAMINTERACTION::CrosslinkerType linkertype,
        int locbspotnum, Core::LinAlg::Matrix<3, 1>& bspotpos,
        Core::LinAlg::Matrix<3, 3>& bspottriad, std::vector<double>& eledisp);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void get_pos_and_triad_of_binding_spot(Core::FE::Discretization const& discret,
        Core::Elements::Element* ele, Core::LinAlg::Vector<double>& ia_discolnp,
        Core::Geo::MeshFree::BoundingBox& pbb, Inpar::BEAMINTERACTION::CrosslinkerType linkertype,
        int locbspotnum, Core::LinAlg::Matrix<3, 1>& bspotpos,
        Core::LinAlg::Matrix<3, 3>& bspottriad);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    bool is_distance_out_of_range(Core::LinAlg::Matrix<3, 1> const& pos1,
        Core::LinAlg::Matrix<3, 1> const& pos2, double const lowerbound, double const upperbound);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    bool is_enclosed_angle_out_of_range(Core::LinAlg::Matrix<3, 1> const& direction1,
        Core::LinAlg::Matrix<3, 1> const& direction2, double const lowerbound,
        double const upperbound);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    bool do_beam_elements_share_nodes(
        Core::Elements::Element const* const beam, Core::Elements::Element const* const nbbeam);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void fe_assemble_ele_force_stiff_into_system_vector_matrix(
        const Core::FE::Discretization& discret, std::vector<int> const& elegid,
        std::vector<Core::LinAlg::SerialDenseVector> const& elevec,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<Epetra_FEVector> fe_sysvec,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> fe_sysmat);


    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     * \brief Get the number of centerline DOFs for a given beam element.
     * @param ele (in) Pointer to the element.
     */
    unsigned int get_number_of_element_centerline_dof(const Core::Elements::Element* elerline_gid);

    /**
     * \brief Get the global indices of the centerline DOFs of a beam element.
     * @param discret (in) Pointer to the discretization.
     * @param ele (in) Pointer to the element.
     * @param ele_centerline_dof_indices (out) Vector with global indices of centerline DOFs in
     * the element.
     */
    template <unsigned int n_centerline_dof>
    void get_element_centerline_gid_indices(Core::FE::Discretization const& discret,
        const Core::Elements::Element* ele,
        Core::LinAlg::Matrix<n_centerline_dof, 1, int>& centerline_gid);

    /**
     * \brief Get the local indices of the centerline DOFs of an element.
     * @param discret (in) Reference to the discretization.
     * @param ele (in) Pointer to the element.
     * @param ele_centerline_dof_indices (out) Vector with local indices of centerline DOFs in the
     * element.
     * @param num_dof (out) Number total DOFs on the element.
     */
    void get_element_centerline_dof_indices(Core::FE::Discretization const& discret,
        const Core::Elements::Element* ele, std::vector<unsigned int>& ele_centerline_dof_indices,
        unsigned int& num_dof);

    /**
     * \brief Get the GID of the rotational DOFs of an element
     * @param discret (in) Reference to the discretization.
     * @param ele (in) Pointer to the element.
     * @return element_rot_gid_indices (out) Vector with GID of rotational DOFs in the
     * element.
     */
    std::vector<int> get_element_rot_gid_indices(
        const Core::FE::Discretization& discret, const Core::Elements::Element* element);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void assemble_centerline_dof_force_stiff_into_element_force_stiff(
        Core::FE::Discretization const& discret, std::vector<int> const& elegid,
        std::vector<Core::LinAlg::SerialDenseVector> const& eleforce_centerlineDOFs,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> const& elestiff_centerlineDOFs,
        std::vector<Core::LinAlg::SerialDenseVector>* eleforce,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>>* elestiff);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     * \brief Assemble a matrix with columns based on centerline DOFs of an element into a matrix
     * with columns based on all DOFs of the element. Example: Mortar coupling matrices as the
     * rows correspond the Lagrange multipliers and the columns correspond to the centerline DOFs.
     * @param discret (in) Pointer to the discretization.
     * @param element (in) Pointer to the element.
     * @param row_matrix_centerlineDOFs (in) Matrix where the columns correspond to the centerline
     * DOFs.
     * @param row_matrix_elementDOFs (out) Matrix where the columns correspond to all Element DOFs
     * (the rest will be 0).
     */
    void assemble_centerline_dof_col_matrix_into_element_col_matrix(
        Core::FE::Discretization const& discret, const Core::Elements::Element* element,
        Core::LinAlg::SerialDenseMatrix const& row_matrix_centerlineDOFs,
        Core::LinAlg::SerialDenseMatrix& row_matrix_elementDOFs);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void extract_pos_dof_vec_absolute_values(Core::FE::Discretization const& discret,
        Core::Elements::Element const* ele, const Core::LinAlg::Vector<double>& ia_discolnp,
        std::vector<double>& element_posdofvec_absolutevalues);
    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void extract_pos_dof_vec_values(Core::FE::Discretization const& discret,
        Core::Elements::Element const* ele, const Core::LinAlg::Vector<double>& ia_discolnp,
        std::vector<double>& element_posdofvec_values);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    template <class T1, class T2>
    void apply_binding_spot_force_to_parent_elements(Core::FE::Discretization const& discret,
        Core::Geo::MeshFree::BoundingBox& pbb, Core::LinAlg::Vector<double>& disp_np_col,
        BEAMINTERACTION::BeamLink& elepairptr,
        std::vector<Core::LinAlg::SerialDenseVector> const& bspotforce,
        std::vector<Core::LinAlg::SerialDenseVector>& eleforce);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    template <class T1, class T2>
    void apply_binding_spot_stiff_to_parent_elements(Core::FE::Discretization const& discret,
        Core::Geo::MeshFree::BoundingBox& pbb, Core::LinAlg::Vector<double>& disp_np_col,
        BEAMINTERACTION::BeamLink& elepairptr,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> const& bspotstiff,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>>& elestiff);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    template <class T1, class T2>
    void apply_binding_spot_force_stiff_to_parent_elements(Core::FE::Discretization const& discret,
        Core::Geo::MeshFree::BoundingBox& pbb, Core::LinAlg::Vector<double>& disp_np_col,
        BEAMINTERACTION::BeamLink& elepairptr,
        std::vector<Core::LinAlg::SerialDenseVector> const& bspotforce,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> const& bspotstiff,
        std::vector<Core::LinAlg::SerialDenseVector>& eleforce,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>>& elestiff);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void setup_ele_type_map_extractor(Teuchos::RCP<const Core::FE::Discretization> const& discret,
        Teuchos::RCP<Core::LinAlg::MultiMapExtractor>& eletypeextractor);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void update_dof_map_of_vector(Core::FE::Discretization& discret,
        Teuchos::RCP<Core::LinAlg::Vector<double>>& dofmapvec,
        Teuchos::RCP<Core::LinAlg::Vector<double>> old = Teuchos::null);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    long long cantor_pairing(std::pair<int, int> const& pair);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    std::pair<int, int> cantor_de_pairing(long long z);

    //! convert element @p ele to bin content type
    Core::Binstrategy::Utils::BinContentType convert_element_to_bin_content_type(
        const Core::Elements::Element* ele);

  }  // namespace Utils
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
