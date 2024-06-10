/*----------------------------------------------------------------------*/
/*! \file

\brief Class containing geometric operations usually needed for the coupling of an embedded
body using a binning strategy to pre-sort the beam elements, for which an octree search needs to be
performed afterwards

\level 3

*----------------------------------------------------------------------*/
#ifndef FOUR_C_FBI_IMMERSED_GEOMETRY_COUPLER_BINNING_HPP
#define FOUR_C_FBI_IMMERSED_GEOMETRY_COUPLER_BINNING_HPP

#include "4C_config.hpp"

#include "4C_binstrategy.hpp"
#include "4C_fbi_immersed_geometry_coupler.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <map>
#include <set>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace BINSTRATEGY
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
    void SetBinning(Teuchos::RCP<BINSTRATEGY::BinningStrategy> binning) override;

    /** \brief Setup the Geoemtry object
     *
     * \param[in] discretizations vector containing the structure and fluid discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    void Setup(std::vector<Teuchos::RCP<Core::FE::Discretization>>&,
        Teuchos::RCP<const Epetra_Vector> structure_displacement) override;

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
    Teuchos::RCP<std::map<int, std::vector<int>>> Search(
        std::vector<Teuchos::RCP<Core::FE::Discretization>>& discretizations,
        Teuchos::RCP<const Epetra_Vector>& column_structure_displacement) override;

    /** \brief Update distribution of elements to bins
     *
     * \param[in] structure_discretization structure discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    void UpdateBinning(Teuchos::RCP<Core::FE::Discretization>& structure_discretization,
        Teuchos::RCP<const Epetra_Vector> structure_column_displacement) override;

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
        Teuchos::RCP<std::map<int, Core::LinAlg::Matrix<3, 1>>> positions,
        Teuchos::RCP<const Epetra_Vector> disp) const override;

    /** \brief Setup the Binning object
     *
     * \param[in] discretizations vector containing the structure and fluid discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    virtual void setup_binning(std::vector<Teuchos::RCP<Core::FE::Discretization>>&,
        Teuchos::RCP<const Epetra_Vector> structure_displacement);

    /** \brief Partition the Problem into bins
     *
     * \param[in] discretizations vector containing the structure and fluid discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    virtual void partition_geometry(std::vector<Teuchos::RCP<Core::FE::Discretization>>&,
        Teuchos::RCP<const Epetra_Vector> structure_displacement);

   private:
    /// binning strategy
    Teuchos::RCP<BINSTRATEGY::BinningStrategy> binstrategy_;

    /// Map relating bins to elements they contain
    std::map<int, std::set<int>> bintoelemap_;

    /// Row map of the bin discretization
    Teuchos::RCP<Epetra_Map> binrowmap_;
  };
}  // namespace FBI

FOUR_C_NAMESPACE_CLOSE

#endif
