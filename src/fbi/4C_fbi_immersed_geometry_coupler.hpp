// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FBI_IMMERSED_GEOMETRY_COUPLER_HPP
#define FOUR_C_FBI_IMMERSED_GEOMETRY_COUPLER_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_vector.hpp"

#include <map>
#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class ConstraintEnforcerFactory;
}

namespace Core::Binstrategy
{
  class BinningStrategy;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Geo
{
  class SearchTree;
}

namespace FBI
{
  class GeometryCouplerFactory;

  /**
   * \brief Class to wrap all generic geometric and parallelism related functionality in order
   * to simulate a beam immersed in fluid
   */
  class FBIGeometryCoupler
  {
    friend GeometryCouplerFactory;

   public:
    /// empty destructor
    virtual ~FBIGeometryCoupler() = default;

    /**
     * \brief Sets the binning strategy in the binning coupler
     *
     * \param[in] binning binning strategy
     */
    virtual void set_binning(std::shared_ptr<Core::Binstrategy::BinningStrategy> binning);

    /** \brief Setup the Geometry object
     *
     * \param[in] discretizations vector containing the structure and fluid discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    virtual void setup(std::vector<std::shared_ptr<Core::FE::Discretization>>&,
        std::shared_ptr<const Core::LinAlg::Vector<double>> structure_displacement);

    /**
     * \brief Performs the search to find possible beam-fluid element pairs
     *
     * Each fluid processor checks for every beam element if it is embedded in one of its fluid
     * elements. The search is implemented via a searchtree.
     *
     * \param[in] discretizations vector containing the structure and fluid discretization
     * respectively
     * \param[in] column_structure_displacement vector containing the structure displacements for
     * all column nodes
     *
     * \returns map relating the beam element IDs to a vector of nearby fluid element IDs
     */
    virtual std::shared_ptr<std::map<int, std::vector<int>>> search(
        std::vector<std::shared_ptr<Core::FE::Discretization>>& discretizations,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& column_structure_displacement);

    /**
     * \brief Ghosts ALL beams to all processors in order to do the search
     *
     *\param[in, out] discretization structure discretization
     *
     * \note For now we assume that we only have beam elements in the structure discretization
     */

    virtual void extend_beam_ghosting(Core::FE::Discretization& discretization);

    /**
     * \brief Handles the parallel communication necessary to create the beam-fluid pairs
     *
     * For now the pairs have to be evaluated on the proc owning the beam element to ensure that
     * every pair is only evaluated once. This function handles all necessary communication to send
     * the data from the owner of the fluid element ( which did the search) to the owner of the beam
     * element ( who ensures uniqueness of the beam segments). This becomes obsolete as soon as the
     * geometry pair is adapted to handle the communication of the tracking data instead.
     *
     * \param[in] discretizations vector containing the structure and fluid discretizations
     * respectively
     * \param[in] pairids  a map containing a map relating all beam element ids to a set
     * of fluid elements ids which they potentially cut
     */
    virtual void prepare_pair_creation(
        std::vector<std::shared_ptr<Core::FE::Discretization>>& discretizations,
        std::shared_ptr<std::map<int, std::vector<int>>> pairids);

    /** \brief Update distribution of elements to bins
     *
     * \param[in] structure_discretization structure discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    virtual void update_binning(std::shared_ptr<Core::FE::Discretization>& structure_discretization,
        std::shared_ptr<const Core::LinAlg::Vector<double>> structure_column_displacement) {};

   protected:
    /**
     * \brief Please use FBI::GeometryCouplerFactory::create_geometry_coupler to create an instance
     * of this class
     */
    FBIGeometryCoupler();

    /**
     * \brief Computes the reference nodal positions needed for the search
     *
     * \param[in] dis field dicretization from which to extract the positions
     * \param[in, out] positions map relating the node IDs to reference positions
     *
     */
    virtual void compute_fixed_positions(Core::FE::Discretization& dis,
        std::shared_ptr<std::map<int, Core::LinAlg::Matrix<3, 1>>> positions) const;

    /**
     * \brief Computes the reference current positions needed for the search
     *
     * \param[in] dis field dicretization from which to extract the positions
     * \param[in, out] positions map relating the node IDs to reference positions
     * \param[in] disp current displacements
     */
    virtual void compute_current_positions(Core::FE::Discretization& dis,
        std::shared_ptr<std::map<int, Core::LinAlg::Matrix<3, 1>>> positions,
        std::shared_ptr<const Core::LinAlg::Vector<double>> disp) const;

    /// Get function for the fludi positions
    virtual std::shared_ptr<const std::map<int, Core::LinAlg::Matrix<3, 1>>> get_fluid_positions()
        const final
    {
      return fluidpositions_;
    };
    /// Get function for the beam positions
    virtual std::shared_ptr<const std::map<int, Core::LinAlg::Matrix<3, 1>>> get_beam_positions()
        const final
    {
      return beampositions_;
    };


   private:
    /// Map storing the nodal positions of the fluid for the search
    std::shared_ptr<std::map<int, Core::LinAlg::Matrix<3, 1>>> fluidpositions_;

    /// Map storing the centerline positions of the beam for the search
    std::shared_ptr<std::map<int, Core::LinAlg::Matrix<3, 1>>> beampositions_;

    /// 3D search tree for embedded discretization
    std::shared_ptr<Core::Geo::SearchTree> searchtree_;

    /** \brief The search radius is used for the pair search and describes the maximum distance of a
     * beam node and a respective fluid node in order for the respective element pair to be found
     * and used for the integration
     *
     */
    const double searchradius_;

    /// Flag determining the fluid stabilization type
    bool edgebased_fluidstabilization_;
  };
}  // namespace FBI

FOUR_C_NAMESPACE_CLOSE

#endif
