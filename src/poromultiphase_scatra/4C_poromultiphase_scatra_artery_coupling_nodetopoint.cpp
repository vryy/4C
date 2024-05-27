/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for coupling 1D node to point in 3D (non-conforming) between
        poromultiphase_scatra-framework and flow in artery networks
        including scalar transport

\level 3

    *----------------------------------------------------------------------*/

#include "4C_poromultiphase_scatra_artery_coupling_nodetopoint.hpp"

#include "4C_discretization_condition_selector.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_pair.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeToPoint::
    PoroMultiPhaseScaTraArtCouplNodeToPoint(Teuchos::RCP<DRT::Discretization> arterydis,
        Teuchos::RCP<DRT::Discretization> contdis, const Teuchos::ParameterList& couplingparams,
        const std::string& condname, const std::string& artcoupleddofname,
        const std::string& contcoupleddofname)
    : PoroMultiPhaseScaTraArtCouplNonConforming(
          arterydis, contdis, couplingparams, condname, artcoupleddofname, contcoupleddofname)
{
  // user info
  if (myrank_ == 0)
  {
    std::cout << "<                                                  >" << std::endl;
    print_out_coupling_method();
    std::cout << "<                                                  >" << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "\n";
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeToPoint::Setup()
{
  // call base class
  POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::Setup();


  // preevaluate coupling pairs
  pre_evaluate_coupling_pairs();

  // print out summary of pairs
  if (contdis_->Name() == "porofluid" &&
      (CORE::UTILS::IntegralValue<int>(couplingparams_, "PRINT_OUT_SUMMARY_PAIRS")))
    output_coupling_pairs();

  // error-checks
  if (has_varying_diam_)
    FOUR_C_THROW("Varying diameter not yet possible for node-to-point coupling");
  if (!evaluate_in_ref_config_)
    FOUR_C_THROW("Evaluation in current configuration not yet possible for node-to-point coupling");

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeToPoint::pre_evaluate_coupling_pairs()
{
  // pre-evaluate
  for (auto& coupl_elepair : coupl_elepairs_) coupl_elepair->pre_evaluate(Teuchos::null);

  // delete the inactive pairs
  coupl_elepairs_.erase(
      std::remove_if(coupl_elepairs_.begin(), coupl_elepairs_.end(),
          [](const Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>
                  coupling_pair) { return not coupling_pair->is_active(); }),
      coupl_elepairs_.end());

  // output
  int total_numactive_pairs = 0;
  int numactive_pairs = static_cast<int>(coupl_elepairs_.size());
  comm().SumAll(&numactive_pairs, &total_numactive_pairs, 1);
  if (myrank_ == 0)
  {
    std::cout << total_numactive_pairs
              << " Artery-to-PoroMultiphaseScatra coupling pairs are active" << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeToPoint::Evaluate(
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
{
  if (!issetup_) FOUR_C_THROW("Setup() has not been called");


  // call base class
  POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::Evaluate(sysmat, rhs);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeToPoint::SetupSystem(
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_cont,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_art, Teuchos::RCP<const Epetra_Vector> rhs_cont,
    Teuchos::RCP<const Epetra_Vector> rhs_art,
    Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_cont,
    Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_art)
{
  // call base class
  POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetupSystem(sysmat, rhs,
      sysmat_cont, sysmat_art, rhs_cont, rhs_art, dbcmap_cont, dbcmap_art->CondMap(),
      dbcmap_art->CondMap());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeToPoint::ApplyMeshMovement()
{
  if (!evaluate_in_ref_config_)
    FOUR_C_THROW("Evaluation in current configuration not possible for node-to-point coupling");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeToPoint::blood_vessel_volume_fraction()
{
  FOUR_C_THROW("Output of vessel volume fraction not possible for node-to-point coupling");

  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeToPoint::print_out_coupling_method()
    const
{
  std::cout << "<Coupling-Method: 1D node to coincident point in 3D>" << std::endl;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeToPoint::output_coupling_pairs() const
{
  if (myrank_ == 0)
  {
    std::cout << "\nSummary of coupling pairs (segments):" << std::endl;
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
  }
  comm().Barrier();
  for (const auto& coupl_elepair : coupl_elepairs_)
  {
    std::cout << "Proc " << std::right << std::setw(2) << myrank_ << ": Artery-ele " << std::right
              << std::setw(5) << coupl_elepair->Ele1GID() << ": <---> continuous-ele " << std::right
              << std::setw(7) << coupl_elepair->Ele2GID() << std::endl;
  }
  comm().Barrier();
  if (myrank_ == 0) std::cout << "\n";
}

FOUR_C_NAMESPACE_CLOSE
