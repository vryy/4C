/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
be directly assembled into the global matrices.


\level 1

*/
/*-----------------------------------------------------------*/


#ifndef PARTITIONED_PENALTYCOUPLING_ASSEMBLY_DIRECT_H
#define PARTITIONED_PENALTYCOUPLING_ASSEMBLY_DIRECT_H

#include <Epetra_FEVector.h>
#include "partitioned_penaltycoupling_assembly_manager.H"

namespace FBI
{
  namespace UTILS
  {
    class FBIAssemblyStrategy;
  }
}  // namespace FBI

namespace BEAMINTERACTION
{
  namespace SUBMODELEVALUATOR
  {
    /**
     * \brief This class collects local force and stiffness terms of the pairs and adds them
     * directly into the global force vector and stiffness matrix.
     */
    class PartitionedBeamInteractionAssemblyManagerDirect
        : public PartitionedBeamInteractionAssemblyManager
    {
     public:
      /**
       * \brief Constructor.
       * \param[in] assembly_contact_elepairs Vector with element pairs to be evaluated by this
       * class.
       * \param[in] assemblystrategy Object determining how the local matrices are assembled into
       * the global one
       */
      PartitionedBeamInteractionAssemblyManagerDirect(
          const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>
              assembly_contact_elepairs,
          Teuchos::RCP<FBI::UTILS::FBIAssemblyStrategy> assemblystrategy);

      /**
       * \brief empty Destructor
       */
      ~PartitionedBeamInteractionAssemblyManagerDirect() override = default;

      /**
       * \brief Evaluate all force and stiffness terms and add them to the global matrices.
       * \param[in] fluid_dis (in) Pointer to the fluid disretization
       * \param[in] solid_dis (in) Pointer to the solid disretization
       * \params[inout] ff Global force vector acting on the fluid
       * \params[inout] fb Global force vector acting on the beam
       * \params[inout] cff  Global stiffness matrix coupling fluid to fluid DOFs
       * \params[inout] cbb  Global stiffness matrix coupling beam to beam DOFs
       * \params[inout] cfb  Global stiffness matrix coupling beam to fluid DOFs
       * \params[inout] cbf  Global stiffness matrix coupling fluid to beam DOFs
       */
      void EvaluateForceStiff(const DRT::Discretization& discretization1,
          const DRT::Discretization& discretization2, Teuchos::RCP<Epetra_FEVector>& ff,
          Teuchos::RCP<Epetra_FEVector>& fb, Teuchos::RCP<LINALG::SparseOperator> cff,
          Teuchos::RCP<LINALG::SparseMatrix>& cbb, Teuchos::RCP<LINALG::SparseMatrix>& cfb,
          Teuchos::RCP<LINALG::SparseMatrix>& cbf, Teuchos::RCP<const Epetra_Vector> fluid_vel,
          Teuchos::RCP<const Epetra_Vector> beam_vel) override;

     protected:
      /// Object determining how the local matrices are assembled into the global one
      Teuchos::RCP<FBI::UTILS::FBIAssemblyStrategy> assemblystrategy_;
    };
  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

#endif
