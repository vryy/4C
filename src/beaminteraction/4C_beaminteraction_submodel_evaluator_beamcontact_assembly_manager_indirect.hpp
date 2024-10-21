// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_ASSEMBLY_MANAGER_INDIRECT_HPP
#define FOUR_C_BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_ASSEMBLY_MANAGER_INDIRECT_HPP


#include "4C_config.hpp"

#include "4C_beaminteraction_submodel_evaluator_beamcontact_assembly_manager.hpp"

FOUR_C_NAMESPACE_OPEN


// Forward declaration.
namespace BEAMINTERACTION
{
  class BeamToSolidMortarManager;
  class BeamToSolidParamsBase;
}  // namespace BEAMINTERACTION


namespace BEAMINTERACTION
{
  namespace SUBMODELEVALUATOR
  {
    /**
     * \brief This class collects local coupling terms of the pairs (D and M) and assembles them
     * into the global coupling matrices. Those global coupling matrices are then multiplied with
     * each other and added to the global force vector and stiffness matrix.
     */
    class BeamContactAssemblyManagerInDirect : public BeamContactAssemblyManager
    {
     public:
      /**
       * \brief Constructor.
       */
      BeamContactAssemblyManagerInDirect(
          const Teuchos::RCP<BEAMINTERACTION::BeamToSolidMortarManager>& mortar_manager)
          : BeamContactAssemblyManager(), mortar_manager_(mortar_manager){};

      /**
       * \brief Evaluate all force and stiffness terms and add them to the global matrices.
       * @param discret (in) Pointer to the disretization.
       * @param data_state (in) Beam interaction data state.
       * @param fe_sysvec (out) Global force vector.
       * @param fe_sysmat (out) Global stiffness matrix.
       */
      void evaluate_force_stiff(Teuchos::RCP<Core::FE::Discretization> discret,
          const Teuchos::RCP<const Solid::ModelEvaluator::BeamInteractionDataState>& data_state,
          Teuchos::RCP<Epetra_FEVector> fe_sysvec,
          Teuchos::RCP<Core::LinAlg::SparseMatrix> fe_sysmat) override;

      /**
       * \brief Return a const pointer to the mortar manager.
       */
      inline Teuchos::RCP<const BEAMINTERACTION::BeamToSolidMortarManager> get_mortar_manager()
          const
      {
        return mortar_manager_;
      }

      double get_energy(
          const Teuchos::RCP<const Core::LinAlg::Vector<double>>& disp) const override;

     private:
      //! Pointer to the mortar manager. This object stores the relevant mortar matrices.
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidMortarManager> mortar_manager_;
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
