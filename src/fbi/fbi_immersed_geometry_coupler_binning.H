/*----------------------------------------------------------------------*/
/*! \file

\brief Class containing geometric operations usually needed for the coupling of an embedded
body using a binning strategy to pre-sort the beam elements, for which an octree search needs to be
performed afterwards

\level 3

*----------------------------------------------------------------------*/
#ifndef FBI_IMMERSED_GEOMETRY_COUPLER_BINNING_H
#define FBI_IMMERSED_GEOMETRY_COUPLER_BINNING_H

#include "linalg_fixedsizematrix.H"
#include <Teuchos_RCP.hpp>
#include <map>
#include <vector>
#include <set>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include "fbi_immersed_geometry_coupler.H"
#include "binstrategy.H"

namespace BINSTRATEGY
{
  class BinningStrategy;
}

namespace DRT
{
  class Discretization;
}

namespace CORE::GEO
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
    /// empty destructor
    virtual ~FBIBinningGeometryCoupler() = default;

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
    virtual void Setup(std::vector<Teuchos::RCP<DRT::Discretization>>&,
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
    virtual Teuchos::RCP<std::map<int, std::vector<int>>> Search(
        std::vector<Teuchos::RCP<DRT::Discretization>>& discretizations,
        Teuchos::RCP<const Epetra_Vector>& column_structure_displacement);

    /** \brief Update distribution of elements to bins
     *
     * \param[in] structure_discretization structure discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    void UpdateBinning(Teuchos::RCP<DRT::Discretization>& structure_discretization,
        Teuchos::RCP<const Epetra_Vector> structure_column_displacement) override;

   protected:
    /**
     * \brief Please use FBI::GeometryCouplerFactory::CreateGeometryCoupler to create an instance of
     * this class
     */
    FBIBinningGeometryCoupler();

    /**
     * \brief Computes the reference current positions needed for the search
     *
     * \param[in] dis field dicretization from which to extract the positions
     * \param[in, out] positions map relating the node IDs to reference positions
     * \param[in] disp current displacements
     */
    void ComputeCurrentPositions(DRT::Discretization& dis,
        Teuchos::RCP<std::map<int, LINALG::Matrix<3, 1>>> positions,
        Teuchos::RCP<const Epetra_Vector> disp) const override;

    /** \brief Setup the Binning object
     *
     * \param[in] discretizations vector containing the structure and fluid discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    virtual void SetupBinning(std::vector<Teuchos::RCP<DRT::Discretization>>&,
        Teuchos::RCP<const Epetra_Vector> structure_displacement);

    /** \brief Partition the Problem into bins
     *
     * \param[in] discretizations vector containing the structure and fluid discretization
     * \param[in] structure_displacement vector containing the column structure displacement
     */
    virtual void PartitionGeometry(std::vector<Teuchos::RCP<DRT::Discretization>>&,
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

#endif
