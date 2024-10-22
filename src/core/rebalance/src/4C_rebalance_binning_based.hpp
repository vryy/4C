// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_REBALANCE_BINNING_BASED_HPP
#define FOUR_C_REBALANCE_BINNING_BASED_HPP

#include "4C_config.hpp"

#include "4C_binstrategy_utils.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

#include <functional>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::IO
{
  class OutputControl;
}

namespace Core::Rebalance
{
  /*! \brief Rebalance discretizations in input vector using BinningStrategy
     *
    \param vector_of_discretizations (in) : vector containing RCPs to discretizations */
  void rebalance_discretizations_by_binning(const Teuchos::ParameterList& binning_params,
      Teuchos::RCP<Core::IO::OutputControl> output_control,
      const std::vector<Teuchos::RCP<Core::FE::Discretization>>& vector_of_discretizations,
      std::function<const Core::Nodes::Node&(const Core::Nodes::Node& node)> correct_node = nullptr,
      std::function<std::vector<std::array<double, 3>>(const Core::FE::Discretization&,
          const Core::Elements::Element&, Teuchos::RCP<const Core::LinAlg::Vector<double>> disnp)>
          determine_relevant_points = nullptr,
      bool revertextendedghosting = false);

  /*!
    \brief Ghost the discretization handed in to this method on all procs

    \param distobeghosted (in) : RCP to the discretisation which is ghosted after execution of
    this method

    \return void */
  void ghost_discretization_on_all_procs(Core::FE::Discretization& distobeghosted);

  /*! \brief Rebalance elements matching to another discretization.
   *
   *    A equal distribution of nodes may lead to unequal distribution
   *    of elements, depending on the algorithm, determining which elements are
   *    owned and which elements are ghosted. This method makes sure, that element
   *    ownerships match and nodes are distributed accordingly.
   *
   * \note The discretization serving as template and the discretization
   *       supposed to be rebalanced have to be matching! That means,
   *       nodal positions have to coincide.
   *
  \param dis_template        (in) : discretization with parallel distr., serving as template
  \param dis_to_rebalance    (in) : discretization which is rebalanced matching to dis_template
*/
  void match_element_distribution_of_matching_discretizations(
      Core::FE::Discretization& dis_template, Core::FE::Discretization& dis_to_rebalance);

  /*! \brief Rebalance conditioned elements matching other conditioned elements.
   *
   * Rebalance the elements of discretization dis_to_rebalance, wearing the
   * condition condname_rebalance to match the parallel distribution of the
   * elements of discretization dis_template, wearing the condition condname_rebalance.
   *
   * \note The elements serving as template and the elements supposed to be rebalanced
   * have to be matching! That means, nodal positions have to coincide.
   *
  \param dis_template        (i) : discretization with parallel distr., serving as template
  \param dis_to_rebalance    (i) : discretization which is rebalanced matching to dis_template
  \param condname_template   (i) : condition on template dis with template distribution
  \param condname_rebalance  (i) : condition on elements to be rebalanced
   \param print  (i) : print elemental and nodal redistribution */
  void match_element_distribution_of_matching_conditioned_elements(
      Core::FE::Discretization& dis_template, Core::FE::Discretization& dis_to_rebalance,
      const std::string& condname_template, const std::string& condname_rebalance,
      bool print = false);

  /*! \brief Get a column vector made of a row vector.

  Using this method, a reference to a vector is returned.
  If the vector is supplied in dis.DofColMap() the vector itself will be returned.
  If the vector is NOT supplied in dis.DofColMap(), but in dis.dof_row_map(),
  a new vector with column map is allocated and the supplied vector is exported to it and
  returned.

  \note The very same functionality is used in discretization::set_state()!

  \param name (in): discretization
  \param state (in): vector of some data  */
  Teuchos::RCP<const Core::LinAlg::Vector<double>> get_col_version_of_row_vector(
      const Core::FE::Discretization& dis,
      const Teuchos::RCP<const Core::LinAlg::Vector<double>> state, const int nds = 0);


  /// recompute nodecolmap of standard discretization to include all nodes as of subdicretization
  Teuchos::RCP<Epetra_Map> compute_node_col_map(
      const Core::FE::Discretization& sourcedis,  ///< standard discretization we want to rebalance
      const Core::FE::Discretization& subdis      ///< subdiscretization prescribing ghosting
  );

  /*! \brief Fill processor local row and col vectors with element ids fitting the desired
  parallel distribution.

  \param dis_template        (i): template for parallel distribution
  \param dis_to_rebalance    (i): discretization supposed to be distributed matching template
  \param row_id_vec_to_fill  (o): on exit this vector contains the matched element row gids
  \param col_id_vec_to_fill  (o): on exit this vector contains the matched element col gids  */
  void match_element_row_col_distribution(const Core::FE::Discretization& dis_template,
      const Core::FE::Discretization& dis_to_rebalance, std::vector<int>& row_id_vec_to_fill,
      std::vector<int>& col_id_vec_to_fill);

  /*! \brief Fill processor local row and col vectors with node ids fitting the desired parallel
  distribution.

  \param dis_template        (i): template for parallel distribution
  \param dis_to_rebalance    (i): discretization supposed to be distributed matching template
  \param row_id_vec_to_fill  (o): on exit this vector contains the matched node row gids
  \param col_id_vec_to_fill  (o): on exit this vector contains the matched node col gids  */
  void match_nodal_row_col_distribution(const Core::FE::Discretization& dis_template,
      const Core::FE::Discretization& dis_to_rebalance, std::vector<int>& row_id_vec_to_fill,
      std::vector<int>& col_id_vec_to_fill);
}  // namespace Core::Rebalance


FOUR_C_NAMESPACE_CLOSE

#endif
