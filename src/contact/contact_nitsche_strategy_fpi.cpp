/*---------------------------------------------------------------------*/
/*! \file

\brief Nitsche contact solving strategy for problems with FPI

\level 3


*/
/*---------------------------------------------------------------------*/

#include "contact_nitsche_strategy_fpi.H"

#include "contact_interface.H"
#include "contact_nitsche_strategy_fsi.H"


void CONTACT::CoNitscheStrategyFpi::SetState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  CONTACT::CoNitscheStrategyPoro::SetState(statename, vec);
  if (statename == MORTAR::state_new_displacement)
  {
    DoContactSearch();
  }
}

void CONTACT::CoNitscheStrategyFpi::DoContactSearch()
{
  for (auto& interface : interface_)
  {
    interface->Initialize();
    interface->EvaluateSearchBinarytree();
    interface->EvaluateNodalNormals();
    interface->ExportNodalNormals();
  }
}

bool CONTACT::CoNitscheStrategyFpi::CheckNitscheContactState(
    CONTACT::CoElement* cele,         // the contact element
    const LINALG::Matrix<2, 1>& xsi,  // local coord on the ele element
    const double& full_fsi_traction,  // stressfluid + penalty
    double& gap                       // gap
)
{
  return CONTACT::UTILS::CheckNitscheContactState(
      *ContactInterfaces()[0], pen_n_, weighting_, cele, xsi, full_fsi_traction, gap);
}
