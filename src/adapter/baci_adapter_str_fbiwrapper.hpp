/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for Fluid-beam interaction problems


\level 2
*/

#ifndef FOUR_C_ADAPTER_STR_FBIWRAPPER_HPP
#define FOUR_C_ADAPTER_STR_FBIWRAPPER_HPP

#include "baci_config.hpp"

#include "baci_adapter_str_fsiwrapper.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace STR
{
  namespace TIMINT
  {
    class ParamsRuntimeOutput;
  }
}  // namespace STR
namespace CORE::LINALG
{
  class MultiMapExtractor;
}  // namespace CORE::LINALG

namespace BEAMINTERACTION
{
  class BeamToSolidVolumeMeshtyingVisualizationOutputParams;
}
namespace ADAPTER
{
  class FBIStructureWrapper : public FSIStructureWrapper
  {
   public:
    /// constructor
    explicit FBIStructureWrapper(Teuchos::RCP<Structure> structure);

    /// extracts interface velocities at \f$t_{n}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();

    /// extracts interface velocities at \f$t_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp();

    /// Predictor for interface velocities
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceVelnp();

    /** \brief linear structure solve with just a interface load
     *
     * Overloads RelaxationSolve of base class with an error message, because it is not implemented
     * for fluid-beam interaction yet
     */
    Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce) override;

    /// switch structure field to block matrix in fsi simulations
    void UseBlockMatrix() override { FOUR_C_THROW("Not yet implemented\n"); };

    /// extract interface displacements at \f$t_{n}\f$
    Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn() override;

    /// extract interface displacements at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp() override;

    /// Predictor for interface displacements
    Teuchos::RCP<Epetra_Vector> PredictInterfaceDispnp() override;

    /** \brief Apply interface forces to structural solver
     *
     * This prepares a new solve of the structural field within one time
     * step. This implementation overloads the base class function because we are not using a
     * condition for the interface, since we expect all beam elements to be immersed.
     */
    void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce) override;

    /// rebuild FSI interface from structure side
    void RebuildInterface() override;

    /// Setup the multi map extractor after ghosting of the structure discretization
    virtual void SetupMultiMapExtractor();

    /// Get Runtime Output data
    virtual Teuchos::RCP<const STR::TIMINT::ParamsRuntimeOutput> GetIOData();

   private:
    /// Map extractor seperating the beam elements from the structure elements
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> eletypeextractor_;

  };  // class FSIStructureWrapper
}  // namespace ADAPTER

FOUR_C_NAMESPACE_CLOSE

#endif
