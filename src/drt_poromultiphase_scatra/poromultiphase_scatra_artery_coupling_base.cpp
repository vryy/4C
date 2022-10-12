/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for coupling between poromultiphase_scatra-
        framework and flow in artery networks including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_artery_coupling_base.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_ParameterListExceptions.hpp>

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase::PoroMultiPhaseScaTraArtCouplBase(
    Teuchos::RCP<DRT::Discretization> arterydis, Teuchos::RCP<DRT::Discretization> contdis,
    const Teuchos::ParameterList& couplingparams, const std::string& condname,
    const std::string& artcoupleddofname, const std::string& contcoupleddofname)
    : arterydis_(arterydis),
      contdis_(contdis),
      myrank_(arterydis->Comm().MyPID()),
      evaluate_in_ref_config_(DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->PoroFluidMultiPhaseDynamicParams().sublist("ARTERY COUPLING"),
          "EVALUATE_IN_REF_CONFIG")),
      comm_(arterydis->Comm())
{
  // safety check
  if (arterydis_->NumGlobalNodes() == 0)
    dserror("artery discretization does not seem to have any nodes");

  // get the actual coupled DOFs
  // 1) 1D artery discretization
  int word1;
  int dummy = 0;
  std::istringstream coupled_art_dof_stream(
      Teuchos::getNumericStringParameter(couplingparams, artcoupleddofname));
  while (coupled_art_dof_stream >> word1)
  {
    // check ascending order
    if (dummy > 0)
      if ((int)(word1 - 1) <= coupleddofs_art_[dummy - 1])
        dserror("DOFs have to be ordered in ascending order");
    coupleddofs_art_.push_back((int)(word1 - 1));
    dummy++;
  }

  // 2) 2D, 3D continuous field discretization
  dummy = 0;
  std::istringstream coupled_poro_dof_stream(
      Teuchos::getNumericStringParameter(couplingparams, contcoupleddofname));
  while (coupled_poro_dof_stream >> word1)
  {
    // check ascending order
    if (dummy > 0)
      if ((int)(word1 - 1) <= coupleddofs_cont_[dummy - 1])
        dserror("DOFs have to be ordered in ascending order");
    coupleddofs_cont_.push_back((int)(word1 - 1));
    dummy++;
  }

  // no coupling selected by user
  if (coupleddofs_cont_.size() == 1 and coupleddofs_art_.size() == 1 and
      coupleddofs_cont_[0] < 0 and coupleddofs_art_[0] < 0)
  {
    coupleddofs_cont_.resize(0);
    coupleddofs_art_.resize(0);
  }

  if (coupleddofs_cont_.size() != coupleddofs_art_.size())
    dserror("size mismatch between COUPLEDDOFS_ART and COUPLEDDOFS_PORO");

  num_coupled_dofs_ = coupleddofs_cont_.size();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase::RecomputeCoupledDOFsForNTP(
    std::vector<DRT::Condition*> coupcond, unsigned int couplingnode)
{
  coupleddofs_art_ = {coupcond[couplingnode]->GetInt("COUPLEDDOF_REDUCED") - 1};
  coupleddofs_cont_ = {coupcond[couplingnode]->GetInt("COUPLEDDOF_PORO") - 1};
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>&
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase::FullMap() const
{
  return globalex_->FullMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<LINALG::MultiMapExtractor>&
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase::GlobalExtractor() const
{
  return globalex_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase::SetSolutionVectors(
    Teuchos::RCP<const Epetra_Vector> phinp_cont, Teuchos::RCP<const Epetra_Vector> phin_cont,
    Teuchos::RCP<const Epetra_Vector> phinp_art)
{
  // do nothing
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase::SetNearbyElePairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  // do nothing
}
