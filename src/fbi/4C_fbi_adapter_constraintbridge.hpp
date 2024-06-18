/*----------------------------------------------------------------------*/
/*! \file

\brief Abstract class to be overloaded by different adapter implementations connecting a constraint
enforcement technique with a discretization approach for Fluid-beam interaction.

\level 2

*----------------------------------------------------------------------*/

#ifndef FOUR_C_FBI_ADAPTER_CONSTRAINTBRIDGE_HPP
#define FOUR_C_FBI_ADAPTER_CONSTRAINTBRIDGE_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

#include <Epetra_FEVector.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <vector>


FOUR_C_NAMESPACE_OPEN

// Forward declaration

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace Core::LinAlg
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace Core::LinAlg
namespace BEAMINTERACTION
{
  class BeamContactPair;
}

namespace GEOMETRYPAIR
{
  class LineTo3DEvaluationData;
}

namespace FBI
{
  class BeamToFluidMeshtyingParams;
  namespace UTILS
  {
    class FBIAssemblyStrategy;
  }
}  // namespace FBI

namespace Adapter
{
  class ConstraintEnforcerFactory;

  /**
   *   \brief Abstract class to be overloaded by different adapter implementations connecting a
   * constraint enforcement technique with a discretization approach for Fluid-beam interaction.
   *
   * The idea is, that this method should act a bit like a mini
   * BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact and manage the fluid-beam interaction pairs and
   * their assembly into global contribution matrices
   */
  class FBIConstraintBridge
  {
    friend ConstraintEnforcerFactory;

   public:
    /// empty destructor
    virtual ~FBIConstraintBridge() = default;

    /**
     * \brief Initializes all members of the class
     *
     * \params[in] beam_map Row Map of the structure discretization
     * \params[in} fluid_map Row Map of the fluid discretization
     * \params[in] fluidmatrix system matrix of the fluid matrix in correct from
     * \params[in] fluidmeshtying bool indicating if fluid meshtying is included
     */
    virtual void setup(const Epetra_Map* beam_map, const Epetra_Map* fluid_map,
        Teuchos::RCP<Core::LinAlg::SparseOperator> fluidmatrix, bool fluidmeshtying);

    /**
     * \brief Computes the coupling matrices
     *
     * This is where the magic happens. The global meshtying contributions are integrated using
     * information of the beam elements, the fluid elements and their position relative to each
     * other.
     *
     */
    virtual void evaluate(Teuchos::RCP<const Core::FE::Discretization> discretization1,
        Teuchos::RCP<const Core::FE::Discretization> discretization2,
        Teuchos::RCP<const Epetra_Vector> fluid_vel,
        Teuchos::RCP<const Epetra_Vector> beam_vel) = 0;

    /**
     * \brief Wraps the ResetState function of the pair
     *
     * Here, the current values of the single fields lie positions and velocities are handed
     * to/updated in the pair
     */
    virtual void ResetPair(const std::vector<double> beam_centerline_dofvec,
        const std::vector<double> fluid_nodal_dofvec,
        Teuchos::RCP<BEAMINTERACTION::BeamContactPair> interactionpair);

    /// Creates a fluid_beam_meshtying pair
    virtual void CreatePair(std::vector<Core::Elements::Element const*> elements,
        std::vector<double> beam_centerline_dofvec, std::vector<double> fluid_nodal_dofvec);

    // Get function for the meshtying pairs meshtying_pairs_
    virtual Teuchos::RCP<std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>> GetPairs()
        const final
    {
      return meshtying_pairs_;
    };

    /// returns data container holding all beam interaction related parameters
    virtual Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> GetParams() const final
    {
      return beam_interaction_params_;
    };

    /// returns data container geometry_evaluation-data_ holding all geometry related evaluation
    /// data
    virtual Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData> GetGeometryData() const final
    {
      return geometry_evaluation_data_;
    };

    /// Clears the pair vector and segmentation information
    virtual void clear();

    /// Resets class members of the bridge
    virtual void ResetBridge() = 0;

    /// Sets the fluid solve flag
    virtual void PrepareFluidSolve() = 0;

    /// Matrix containing only structure side contributions \f$C_{ss}\f$
    virtual Teuchos::RCP<const Core::LinAlg::SparseMatrix> GetCss() const = 0;

    /// Matrix containing only fluid side contributions \f$C_{ff}\f$
    virtual Teuchos::RCP<const Core::LinAlg::SparseOperator> GetCff() const = 0;

    /// Matrix containing mixed fluid side contributions \f$C_{fs}\f$
    virtual Teuchos::RCP<const Core::LinAlg::SparseMatrix> GetCfs() const = 0;

    /// Matrix containing mixed structure side contributions \f$C_{sf}\f$
    virtual Teuchos::RCP<const Core::LinAlg::SparseMatrix> GetCsf() const = 0;

    /// Force vector acting on the fluid side \f$f_f\f$
    virtual Teuchos::RCP<const Epetra_FEVector> get_fluid_coupling_residual() const = 0;

    /// Force vector acting on the structure side \f$f_s\f$
    virtual Teuchos::RCP<const Epetra_FEVector> get_structure_coupling_residual() const = 0;

   protected:
    /** \brief You will have to use the Adapter::ConstraintEnforcerFactory
     *
     */
    FBIConstraintBridge();

    /// data container holding all beam interaction related parameters
    Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> beam_interaction_params_;

    /// Store the assembly strategy here to hand into the assembly manager
    Teuchos::RCP<FBI::UTILS::FBIAssemblyStrategy> assemblystrategy_;

   private:
    /// meshtying pairs
    Teuchos::RCP<std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>> meshtying_pairs_;

    /// data container holding all geometry related evaluation data
    Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData> geometry_evaluation_data_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
