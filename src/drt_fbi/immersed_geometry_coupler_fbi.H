/*----------------------------------------------------------------------*/
/*! \file

\brief Class containing geometric operations usually needed for the coupling of an embedded body

\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/
#ifndef GEOMETRYCOUPLER_FBI_H
#define GEOMETRYCOUPLER_FBI_H

#include "../linalg/linalg_fixedsizematrix.H"
#include <Teuchos_RCP.hpp>
#include <map>
#include <vector>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

namespace ADAPTER
{
  class ConstraintEnforcerFactory;
}

namespace DRT
{
  class Discretization;
}

namespace GEO
{
  class SearchTree;
}

namespace FBI
{
  /**
   * \brief Class to wrap all generic geometric and parallelism related functionality in order
   * to simulate a beam immersed in fluid
   */
  class FBIGeometryCoupler
  {
    friend ADAPTER::ConstraintEnforcerFactory;

   public:
    /// empty destructor
    virtual ~FBIGeometryCoupler(){};

    /** \brief Setup the Geoemtry object
     *
     * \param[in] discretizations vector containing the structure and fluid discretization
     */
    virtual void Setup(std::vector<Teuchos::RCP<DRT::Discretization>>&);

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

    /**
     * \brief Ghosts ALL beams to all processors in order to do the search
     *
     *\param[in, out] discretization structure discretization
     *
     * \note For now we assume that we only have beam elements in the structure discretization
     */

    virtual void ExtendBeamGhosting(DRT::Discretization& discretization);

    /**
     * \brief Handles the parallel communication necessary to create the beam-fluid pairs
     *
     * For now the pairs have to be evaluated on the proc owning the beam element to ensure that
     * every pair is only evaluated once. This function handles all necessary communication to send
     * the data from the owner of the fluid element ( which did the search) to the owner of the beam
     * element ( who ensures uniqueness of the beam segments). This becomes obsolete as soon as the
     * geoemtry pair is adapted to handle the communication of the tracking data instead.
     *
     * \param[in] discretizations vector containing the structure and fluid discretizations
     * respectively
     * \param[in] pairids  a map containing a map relating all beam element ids to a set
     * of fluid elements ids which they potentially cut
     */
    virtual void PreparePairCreation(
        std::vector<Teuchos::RCP<DRT::Discretization>>& discretizations,
        Teuchos::RCP<std::map<int, std::vector<int>>> pairids);

   protected:
    FBIGeometryCoupler();
    /**
     * \brief Computes the reference nodal positions needed for the search
     *
     * \param[in] dis field dicretization from which to extract the positions
     * \param[in, out] positions map relating the node IDs to reference positions
     *
     */
    virtual void ComputeFixedPositions(DRT::Discretization& dis,
        Teuchos::RCP<std::map<int, LINALG::Matrix<3, 1>>> positions) const;

    /**
     * \brief Computes the reference current positions needed for the search
     *
     * \param[in] dis field dicretization from which to extract the positions
     * \param[in, out] positions map relating the node IDs to reference positions
     * \param[in] disp current displacements
     */
    virtual void ComputeCurrentPositions(DRT::Discretization& dis,
        Teuchos::RCP<std::map<int, LINALG::Matrix<3, 1>>> positions,
        Teuchos::RCP<const Epetra_Vector> disp) const;

    /// Get function for the fludi positions
    virtual Teuchos::RCP<const std::map<int, LINALG::Matrix<3, 1>>> GetFluidPositions() const final
    {
      return fluidpositions_;
    };
    /// Get function for the beam positions
    virtual Teuchos::RCP<const std::map<int, LINALG::Matrix<3, 1>>> GetBeamPositions() const final
    {
      return beampositions_;
    };

    /// Get function for the search tree
    virtual Teuchos::RCP<const GEO::SearchTree> GetSearchtree() const final { return searchtree_; };

   private:
    /// Map storing the nodal positions of the fluid for the search
    Teuchos::RCP<std::map<int, LINALG::Matrix<3, 1>>> fluidpositions_;

    /// Map storing the centerline positions of the beam for the search
    Teuchos::RCP<std::map<int, LINALG::Matrix<3, 1>>> beampositions_;

    /// 3D search tree for embedded discretization
    Teuchos::RCP<GEO::SearchTree> searchtree_;

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

#endif
