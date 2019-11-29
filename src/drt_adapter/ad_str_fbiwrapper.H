/*----------------------------------------------------------------------*/
/*! \file
\file ad_str_fbiwrapper.H

\brief Structural adapter for Fluid-beam interaction problems containing the interface to the
structure and methods dependent on the interface

\maintainer Nora Hagmeyer

\level 3
*/

#ifndef FBI_STRUCTURE_WRAPPER_H
#define FBI_STRUCTURE_WRAPPER_H

#include "ad_str_fsiwrapper.H"
#include "../drt_fsi/fsi_str_model_evaluator_partitioned.H"

// forward declarations

namespace LINALG
{
  class MultiMapExtractor;
}  // namespace LINALG

namespace ADAPTER
{
  class FBIStructureWrapper : public FSIStructureWrapper
  {
   public:
    /// constructor
    explicit FBIStructureWrapper(Teuchos::RCP<Structure> structure);

    /// Destructor
    virtual ~FBIStructureWrapper() { ; }

    /// extracts interface velocities at \f$t_{n}\f$
    Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();

    /// extracts interface velocities at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp();

    /// Predictor for interface velocities
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceVelnp();

    /** \brief linear structure solve with just a interface load
     *
     * Overloads RelaxationSolve of base class with an error message, because it is not implemented
     * for fluid-beam interaction yet
     */
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(
        Teuchos::RCP<Epetra_Vector> iforce) override;

    /// switch structure field to block matrix in fsi simulations
    virtual void UseBlockMatrix() override { dserror("Not yet implemented\n"); };

    /// extract interface displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn() override;

    /// extract interface displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp() override;

    /// Predictor for interface displacements
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceDispnp() override;

    /** \brief Apply interface forces to structural solver
     *
     * This prepares a new solve of the structural field within one time
     * step. This implementation overloads the base class function because we are not using a
     * condition for the interface, since we expect all beam elements to be immersed.
     */
    virtual void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce) override;

    /// rebuild FSI interface from structure side
    virtual void RebuildInterface() override;

   private:
    /// Map extractor seperating the beam elements from the structure elements
    Teuchos::RCP<LINALG::MultiMapExtractor> eletypeextractor_;

  };  // class FSIStructureWrapper
}  // namespace ADAPTER

#endif
