/*----------------------------------------------------------------------*/
/*! \file

\brief Class containing geometric operations usually needed for the coupling of an embedded body

\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/
#ifndef GEOMETRYCOUPLER_FBI_H
#define GEOMETRYCOUPLER_FBI_H

#include <Teuchos_RCP.hpp>
#include <map>
#include <vector>
#include <Epetra_Vector.h>

namespace DRT
{
  class Discretization;
}

namespace GEO
{
  class SearchTree;
}

namespace ADAPTER
{
  class ConstraintEnforcerFactory;
}
namespace FBI
{
  /**
   * \brief Class to encapsule all generic geometric and parallelism related functionality in order
   * to simulate a beam immersed in fluid
   */
  class FBIGeometryCoupler
  {
    friend ADAPTER::ConstraintEnforcerFactory;

   public:
    /// empty destructor
    virtual ~FBIGeometryCoupler(){};

    /// Setup the Geoemtry object
    virtual void Setup(Teuchos::RCP<std::vector<Teuchos::RCP<DRT::Discretization>>>);

    virtual Teuchos::RCP<std::map<int, std::vector<int>>> Search(
        Teuchos::RCP<std::vector<Teuchos::RCP<DRT::Discretization>>> discretizations,
        Teuchos::RCP<const Epetra_Vector>& column_structure_displacement,
        Teuchos::RCP<const Epetra_Vector>& column_structure_velocity,
        Teuchos::RCP<const Epetra_Vector>& column_fluid_velocity);

    /**
     * \brief Ghosts ALL beams to all processors in order to do the search
     *
     * \note For now we assume that we only have beam elements in the structure discretization
     */

    virtual void ExtendBeamGhosting(Teuchos::RCP<DRT::Discretization> discretizations);

    /**
     * \brief Handles the parallel communication necessary to create the beam-fluid pairs
     *
     * For now the pairs have to be evaluated on the proc owning the beam element to ensure that
     * every pair is only evaluated once. This function handles all necessary communication to send
     * the data from the owner of the fluid element ( which did the search) to the owner of the beam
     * element ( who ensures uniqueness of the beam segments). This becomes obsolete as soon as the
     * geoemtry pair is adapted to handle the communication of the tracking data instead.
     *
     * \params[in] pairids  a map containing a map relating all beam element ids to a set of fluid
     * elements ids which they potentially cut
     */
    virtual void PreparePairCreation(
        Teuchos::RCP<std::vector<Teuchos::RCP<DRT::Discretization>>> discretizations,
        Teuchos::RCP<std::map<int, std::vector<int>>> pairids);

   protected:
    FBIGeometryCoupler();
    /**
     * \brief Computes the reference nodal positions needed for the search
     *
     */
    virtual void ComputeFixedPositions(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<std::map<int, LINALG::Matrix<3, 1>>> positions) const;

    /**
     * \brief Computes the reference current positions needed for the search
     */
    virtual void ComputeCurrentPositions(Teuchos::RCP<DRT::Discretization> dis,
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
