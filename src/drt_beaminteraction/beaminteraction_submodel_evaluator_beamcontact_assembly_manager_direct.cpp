/*-----------------------------------------------------------*/
/*!

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
be directly assembled into the global matrices.

\maintainer Ivo Steinbrecher

\level 3

*/


#include "beaminteraction_submodel_evaluator_beamcontact_assembly_manager_direct.H"

#include "beam_contact_pair.H"
#include "beaminteraction_calc_utils.H"

#include "../drt_lib/drt_element.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

/**
 *
 */
BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerDirect::
    BeamContactAssemblyManagerDirect(
        std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> assembly_contact_elepairs)
    : BeamContactAssemblyManager(assembly_contact_elepairs)
{
}


/**
 *
 */
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerDirect::EvaluateForceStiff(
    Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<Epetra_FEVector> fe_sysvec,
    Teuchos::RCP<LINALG::SparseMatrix> fe_sysmat, Teuchos::RCP<const Epetra_Vector> disp)
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
      BEAMINTERACTION::UTILS::AssembleCenterlineDofForceStiffIntoElementForceStiff(*discret,
          elegids, eleforce_centerlineDOFs, elestiff_centerlineDOFs, &eleforce, &elestiff);


      // Fixme
      eleforce[0].Scale(-1.0);
      eleforce[1].Scale(-1.0);

      // assemble the contributions into force vector class variable
      // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
      BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix(
          *discret, elegids, eleforce, elestiff, fe_sysvec, fe_sysmat);
    }
  }
}
