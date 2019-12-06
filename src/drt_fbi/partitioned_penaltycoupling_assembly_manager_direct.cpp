/*-----------------------------------------------------------*/
/*! \file
\file partitioned_penaltycoupling_assembly_manager_direct.cpp

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
be directly assembled into the global matrices.

\maintainer Nora Hagmeyer

\level 1

*/
/*-----------------------------------------------------------*/

#include "partitioned_penaltycoupling_assembly_manager_direct.H"

#include "../drt_beaminteraction/beam_contact_pair.H"
#include "../drt_beaminteraction/beaminteraction_calc_utils.H"
#include "../drt_fbi/fbi_calc_utils.H"

#include "../drt_lib/drt_element.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

/**
 *
 */
BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerDirect::
    PartitionedBeamInteractionAssemblyManagerDirect(
        std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> assembly_contact_elepairs)
    : PartitionedBeamInteractionAssemblyManager(assembly_contact_elepairs)
{
}


/**
 *
 */
void BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerDirect::
    EvaluateForceStiff(const DRT::Discretization& discretization1,
        const DRT::Discretization& discretization2, Teuchos::RCP<Epetra_FEVector>& ff,
        Teuchos::RCP<Epetra_FEVector>& fb, Teuchos::RCP<LINALG::SparseMatrix>& cff,
        Teuchos::RCP<LINALG::SparseMatrix>& cbb, Teuchos::RCP<LINALG::SparseMatrix>& cfb,
        Teuchos::RCP<LINALG::SparseMatrix>& cbf)
{
  // resulting discrete element force vectors of the two interacting elements
  std::vector<LINALG::SerialDenseVector> eleforce(2);

  // resulting discrete force vectors (centerline DOFs only!) of the two
  // interacting elements
  std::vector<LINALG::SerialDenseVector> eleforce_centerlineDOFs(2);

  // linearizations
  std::vector<std::vector<LINALG::SerialDenseMatrix>> elestiff(
      2, std::vector<LINALG::SerialDenseMatrix>(2));

  // linearizations (centerline DOFs only!)
  std::vector<std::vector<LINALG::SerialDenseMatrix>> elestiff_centerlineDOFs(
      2, std::vector<LINALG::SerialDenseMatrix>(2));

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
    pair_is_active = elepairptr->Evaluate(&eleforce_centerlineDOFs[0], &eleforce_centerlineDOFs[1],
        &elestiff_centerlineDOFs[0][0], &elestiff_centerlineDOFs[0][1],
        &elestiff_centerlineDOFs[1][0], &elestiff_centerlineDOFs[1][1]);

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
      BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrices(discretization1,
          discretization2, elegids, eleforce, elestiff, ff, fb, cff, cbb, cfb, cbf);
    }
  }
  int err = fb->GlobalAssemble();
  if (err) printf("Global assembly failed with error %i", err);
  err = ff->GlobalAssemble();
  if (err) printf("Global assembly failed with error %i", err);
}
