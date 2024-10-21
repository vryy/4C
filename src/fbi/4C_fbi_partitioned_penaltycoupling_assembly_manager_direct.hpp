// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FBI_PARTITIONED_PENALTYCOUPLING_ASSEMBLY_MANAGER_DIRECT_HPP
#define FOUR_C_FBI_PARTITIONED_PENALTYCOUPLING_ASSEMBLY_MANAGER_DIRECT_HPP

#include "4C_config.hpp"

#include "4C_fbi_partitioned_penaltycoupling_assembly_manager.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN

namespace FBI
{
  namespace Utils
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
          Teuchos::RCP<FBI::Utils::FBIAssemblyStrategy> assemblystrategy);

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
      void evaluate_force_stiff(const Core::FE::Discretization& discretization1,
          const Core::FE::Discretization& discretization2, Teuchos::RCP<Epetra_FEVector>& ff,
          Teuchos::RCP<Epetra_FEVector>& fb, Teuchos::RCP<Core::LinAlg::SparseOperator> cff,
          Teuchos::RCP<Core::LinAlg::SparseMatrix>& cbb,
          Teuchos::RCP<Core::LinAlg::SparseMatrix>& cfb,
          Teuchos::RCP<Core::LinAlg::SparseMatrix>& cbf,
          Teuchos::RCP<const Core::LinAlg::Vector<double>> fluid_vel,
          Teuchos::RCP<const Core::LinAlg::Vector<double>> beam_vel) override;

     protected:
      /// Object determining how the local matrices are assembled into the global one
      Teuchos::RCP<FBI::Utils::FBIAssemblyStrategy> assemblystrategy_;
    };
  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
