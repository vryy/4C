/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
be directly assembled into the global matrices.

\maintainer Nora Hagmeyer

\level 1

*/
/*-----------------------------------------------------------*/


#ifndef PARTITIONED_PENALTYCOUPLING_ASSEMBLY_DIRECT_H
#define PARTITIONED_PENALTYCOUPLING_ASSEMBLY_DIRECT_H

#include <Epetra_FEVector.h>
#include "partitioned_penaltycoupling_assembly_manager.H"


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
       * \params[in] assembly_contact_elepairs Vector with element pairs to be evaluated by this
       * class.
       */
      PartitionedBeamInteractionAssemblyManagerDirect(
          const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>
              assembly_contact_elepairs);

      /**
       * \brief empty Destructor
       */
      virtual ~PartitionedBeamInteractionAssemblyManagerDirect(){};

      /**
       * \brief Evaluate all force and stiffness terms and add them to the global matrices.
       * \params[in] fluid_dis (in) Pointer to the fluid disretization
       * \params[in] solid_dis (in) Pointer to the solid disretization
       * \params[inout] ff Global force vector acting on the fluid
       * \params[inout] fb Global force vector acting on the beam
       * \params[inout] cff  Global stiffness matrix coupling fluid to fluid DOFs
       * \params[inout] cbb  Global stiffness matrix coupling beam to beam DOFs
       * \params[inout] cfb  Global stiffness matrix coupling beam to fluid DOFs
       * \params[inout] cbf  Global stiffness matrix coupling fluid to beam DOFs
       */
      void EvaluateForceStiff(
          const Teuchos::RCP<const std::vector<Teuchos::RCP<DRT::Discretization>>>& discretizations,
          Teuchos::RCP<Epetra_FEVector>& ff, Teuchos::RCP<Epetra_FEVector>& fb,
          Teuchos::RCP<LINALG::SparseMatrix>& cff, Teuchos::RCP<LINALG::SparseMatrix>& cbb,
          Teuchos::RCP<LINALG::SparseMatrix>& cfb,
          Teuchos::RCP<LINALG::SparseMatrix>& cbf) override;
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

#endif
