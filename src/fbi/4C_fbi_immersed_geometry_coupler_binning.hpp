// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FBI_IMMERSED_GEOMETRY_COUPLER_BINNING_HPP
#define FOUR_C_FBI_IMMERSED_GEOMETRY_COUPLER_BINNING_HPP

#include "4C_config.hpp"

#include "4C_binstrategy.hpp"
#include "4C_fbi_immersed_geometry_coupler.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_vector.hpp"

#include <Epetra_Map.h>

#include <map>
#include <memory>
#include <set>
#include <vector>

FOUR_C_NAMESPACE_OPEN

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
  class GeometrycouplerFactory;
  /**
   * \brief Class to wrap all generic geometric and parallelism related functionality in order
   * to simulate a beam immersed in fluid using binning as a presorting strategy
   */
  class FBIBinningGeometryCoupler : public FBIGeometryCoupler
  {
    friend FBI::GeometryCouplerFactory;

   public:
    /**
     * \brief Sets the binning strategy in the binning coupler
     *
     * \param[in] binning binning strategy
     */
    void set_binning(std::shared_ptr<Core::Binstrategy::BinningStrategy> binning) override;

    /** \brief Setup the Geometry object
     *
     * \param[in] discretizations vector containing the structure and fluid discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    void setup(std::vector<std::shared_ptr<Core::FE::Discretization>>&,
        std::shared_ptr<const Core::LinAlg::Vector<double>> structure_displacement) override;

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
    std::shared_ptr<std::map<int, std::vector<int>>> search(
        std::vector<std::shared_ptr<Core::FE::Discretization>>& discretizations,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& column_structure_displacement)
        override;

    /** \brief Update distribution of elements to bins
     *
     * \param[in] structure_discretization structure discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    void update_binning(std::shared_ptr<Core::FE::Discretization>& structure_discretization,
        std::shared_ptr<const Core::LinAlg::Vector<double>> structure_column_displacement) override;

   protected:
    /**
     * \brief Please use FBI::GeometryCouplerFactory::create_geometry_coupler to create an instance
     * of this class
     */
    FBIBinningGeometryCoupler();

    /**
     * \brief Computes the reference current positions needed for the search
     *
     * \param[in] dis field dicretization from which to extract the positions
     * \param[in, out] positions map relating the node IDs to reference positions
     * \param[in] disp current displacements
     */
    void compute_current_positions(Core::FE::Discretization& dis,
        std::shared_ptr<std::map<int, Core::LinAlg::Matrix<3, 1>>> positions,
        std::shared_ptr<const Core::LinAlg::Vector<double>> disp) const override;

    /** \brief Setup the Binning object
     *
     * \param[in] discretizations vector containing the structure and fluid discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    virtual void setup_binning(std::vector<std::shared_ptr<Core::FE::Discretization>>&,
        std::shared_ptr<const Core::LinAlg::Vector<double>> structure_displacement);

    /** \brief Partition the Problem into bins
     *
     * \param[in] discretizations vector containing the structure and fluid discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    virtual void partition_geometry(std::vector<std::shared_ptr<Core::FE::Discretization>>&,
        std::shared_ptr<const Core::LinAlg::Vector<double>> structure_displacement);

   private:
    /// binning strategy
    std::shared_ptr<Core::Binstrategy::BinningStrategy> binstrategy_;

    /// Map relating bins to elements they contain
    std::map<int, std::set<int>> bintoelemap_;

    /// Row map of the bin discretization
    std::shared_ptr<Epetra_Map> binrowmap_;
  };
}  // namespace FBI

FOUR_C_NAMESPACE_CLOSE

#endif
