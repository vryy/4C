/*-----------------------------------------------------------*/
/*! \file
\file partitioned_penaltycoupling_assembly_manager_direct.cpp

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
be directly assembled into the global matrices.


\level 1

*/
/*-----------------------------------------------------------*/

#include "4C_fbi_partitioned_penaltycoupling_assembly_manager_direct.hpp"

#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_fbi_calc_utils.hpp"
#include "4C_fbi_fluid_assembly_strategy.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_sparseoperator.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerDirect::
    PartitionedBeamInteractionAssemblyManagerDirect(
        std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> assembly_contact_elepairs,
        Teuchos::RCP<FBI::UTILS::FBIAssemblyStrategy> assemblystrategy)
    : PartitionedBeamInteractionAssemblyManager(assembly_contact_elepairs),
      assemblystrategy_(assemblystrategy)
{
}


/**
 *
 */
void BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerDirect::
    evaluate_force_stiff(const Core::FE::Discretization& discretization1,
        const Core::FE::Discretization& discretization2, Teuchos::RCP<Epetra_FEVector>& ff,
        Teuchos::RCP<Epetra_FEVector>& fb, Teuchos::RCP<Core::LinAlg::SparseOperator> cff,
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& cbb,
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& cfb,
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& cbf, Teuchos::RCP<const Epetra_Vector> fluid_vel,
        Teuchos::RCP<const Epetra_Vector> beam_vel)
{
  // resulting discrete element force vectors of the two interacting elements
  std::vector<Core::LinAlg::SerialDenseVector> eleforce(2);

  // resulting discrete force vectors (centerline DOFs only!) of the two
  // interacting elements
  std::vector<Core::LinAlg::SerialDenseVector> eleforce_centerlineDOFs(2);

  // linearizations
  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> elestiff(
      2, std::vector<Core::LinAlg::SerialDenseMatrix>(2));

  // linearizations (centerline DOFs only!)
  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> elestiff_centerlineDOFs(
      2, std::vector<Core::LinAlg::SerialDenseMatrix>(2));

  // element gids of interacting elements
  std::vector<int> elegids(2);

  // are non-zero stiffness values returned which need assembly?
  bool pair_is_active = false;


  for (auto& elepairptr : assembly_contact_elepairs_)
  {
    // pre_evaluate the pair
    elepairptr->pre_evaluate();
  }

  for (auto& elepairptr : assembly_contact_elepairs_)
  {
    // Evaluate the pair and check if there is active contact
    pair_is_active =
        elepairptr->Evaluate(&(eleforce_centerlineDOFs[0]), &(eleforce_centerlineDOFs[1]),
            &(elestiff_centerlineDOFs[0][0]), &(elestiff_centerlineDOFs[0][1]),
            &(elestiff_centerlineDOFs[1][0]), &(elestiff_centerlineDOFs[1][1]));

    if (pair_is_active)
    {
      elegids[0] = elepairptr->Element1()->Id();
      elegids[1] = elepairptr->Element2()->Id();

      // assemble force vector and stiffness matrix affecting the centerline DoFs only
      // into element force vector and stiffness matrix ('all DoFs' format, as usual)
      FBI::UTILS::AssembleCenterlineDofForceStiffIntoFBIElementForceStiff(discretization1,
          discretization2, elegids, eleforce_centerlineDOFs, elestiff_centerlineDOFs, &eleforce,
          &elestiff);

      // assemble the contributions into force and stiffness matrices
      assemblystrategy_->Assemble(discretization1, discretization2, elegids, eleforce, elestiff, fb,
          ff, cbb, cff, cbf, cfb);
    }
  }
  int err = fb->GlobalAssemble();
  if (err) printf("Global assembly failed with error %i", err);
  err = ff->GlobalAssemble();
  if (err) printf("Global assembly failed with error %i", err);
}

FOUR_C_NAMESPACE_CLOSE
