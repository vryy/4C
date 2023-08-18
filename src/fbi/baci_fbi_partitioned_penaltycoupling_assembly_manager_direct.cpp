/*-----------------------------------------------------------*/
/*! \file
\file partitioned_penaltycoupling_assembly_manager_direct.cpp

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
be directly assembled into the global matrices.


\level 1

*/
/*-----------------------------------------------------------*/

#include "baci_fbi_partitioned_penaltycoupling_assembly_manager_direct.H"

#include "baci_beaminteraction_calc_utils.H"
#include "baci_beaminteraction_contact_pair.H"
#include "baci_fbi_calc_utils.H"
#include "baci_fbi_fluid_assembly_strategy.H"
#include "baci_lib_element.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_linalg_sparsematrix.H"
#include "baci_linalg_sparseoperator.H"


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
    EvaluateForceStiff(const DRT::Discretization& discretization1,
        const DRT::Discretization& discretization2, Teuchos::RCP<Epetra_FEVector>& ff,
        Teuchos::RCP<Epetra_FEVector>& fb, Teuchos::RCP<CORE::LINALG::SparseOperator> cff,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& cbb,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& cfb,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& cbf, Teuchos::RCP<const Epetra_Vector> fluid_vel,
        Teuchos::RCP<const Epetra_Vector> beam_vel)
{
  // resulting discrete element force vectors of the two interacting elements
  std::vector<CORE::LINALG::SerialDenseVector> eleforce(2);

  // resulting discrete force vectors (centerline DOFs only!) of the two
  // interacting elements
  std::vector<CORE::LINALG::SerialDenseVector> eleforce_centerlineDOFs(2);

  // linearizations
  std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> elestiff(
      2, std::vector<CORE::LINALG::SerialDenseMatrix>(2));

  // linearizations (centerline DOFs only!)
  std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> elestiff_centerlineDOFs(
      2, std::vector<CORE::LINALG::SerialDenseMatrix>(2));

  // element gids of interacting elements
  std::vector<int> elegids(2);

  // are non-zero stiffness values returned which need assembly?
  bool pair_is_active = false;


  for (auto& elepairptr : assembly_contact_elepairs_)
  {
    // PreEvaluate the pair
    elepairptr->PreEvaluate();
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
