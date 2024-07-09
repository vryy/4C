/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global mortar matrices. The pairs in this
class can not be directly assembled into the global matrices. They have to be assembled into the
global coupling matrices M and D first.


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_FBI_PARTITIONED_PENALTYCOUPLING_ASSEMBLY_MANAGER_INDIRECT_HPP
#define FOUR_C_FBI_PARTITIONED_PENALTYCOUPLING_ASSEMBLY_MANAGER_INDIRECT_HPP


#include "4C_config.hpp"

#include "4C_fbi_partitioned_penaltycoupling_assembly_manager.hpp"

FOUR_C_NAMESPACE_OPEN


// Forward declaration.
namespace BEAMINTERACTION
{
  class BeamToFluidMortarManager;
}  // namespace BEAMINTERACTION

namespace FBI
{
  class BeamToFluidMeshtyingParams;
}
namespace BEAMINTERACTION
{
  namespace SUBMODELEVALUATOR
  {
    /**
     * \brief This class collects local mortar coupling terms of the pairs (D and M) and assembles
     * them into the global coupling matrices. Those global coupling matrices are then multiplied
     * with each other and added to the global force vector and stiffness matrix.
     */
    class PartitionedBeamInteractionAssemblyManagerIndirect
        : public PartitionedBeamInteractionAssemblyManager
    {
     public:
      /**
       * \brief Constructor.
       * @param assembly_contact_elepairs (in) Vector with element pairs to be evaluated by this
       * class.
       */
      PartitionedBeamInteractionAssemblyManagerIndirect(
          std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>& assembly_contact_elepairs,
          Teuchos::RCP<const Core::FE::Discretization>& discretization1,
          Teuchos::RCP<const Core::FE::Discretization>& discretization2,
          Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> beam_contact_params_ptr);

      /**
       * \brief Evaluate all force and stiffness terms and add them to the global matrices.
       * @param discret (in) Pointer to the disretization.
       * @param fe_sysvec (out) Global force vector.
       * @param fe_sysmat (out) Global stiffness matrix.
       * @param disp (in) Current displacement vector.
       */
      void evaluate_force_stiff(const Core::FE::Discretization& discretization1,
          const Core::FE::Discretization& discretization2, Teuchos::RCP<Epetra_FEVector>& ff,
          Teuchos::RCP<Epetra_FEVector>& fb, Teuchos::RCP<Core::LinAlg::SparseOperator> cff,
          Teuchos::RCP<Core::LinAlg::SparseMatrix>& cbb,
          Teuchos::RCP<Core::LinAlg::SparseMatrix>& cfb,
          Teuchos::RCP<Core::LinAlg::SparseMatrix>& cbf,
          Teuchos::RCP<const Epetra_Vector> fluid_vel,
          Teuchos::RCP<const Epetra_Vector> beam_vel) override;

      /**
       * \brief Return a const pointer to the mortar manager.
       */
      inline Teuchos::RCP<const BEAMINTERACTION::BeamToFluidMortarManager> get_mortar_manager()
          const
      {
        return mortar_manager_;
      }

     private:
      //! Pointer to the mortar manager. This object stores the relevant mortar matrices.
      Teuchos::RCP<BEAMINTERACTION::BeamToFluidMortarManager> mortar_manager_;
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
