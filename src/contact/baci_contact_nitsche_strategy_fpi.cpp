/*---------------------------------------------------------------------*/
/*! \file

\brief Nitsche contact solving strategy for problems with FPI

\level 3


*/
/*---------------------------------------------------------------------*/

#include "baci_contact_nitsche_strategy_fpi.H"

#include "baci_contact_interface.H"
#include "baci_contact_nitsche_strategy_fsi.H"

BACI_NAMESPACE_OPEN


void CONTACT::NitscheStrategyFpi::SetState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  CONTACT::NitscheStrategyPoro::SetState(statename, vec);
  if (statename == MORTAR::state_new_displacement)
  {
    DoContactSearch();
  }
}

void CONTACT::NitscheStrategyFpi::DoContactSearch()
{
  for (auto& interface : interface_)
  {
    interface->Initialize();
    interface->EvaluateSearchBinarytree();
    interface->EvaluateNodalNormals();
    interface->ExportNodalNormals();
  }
}

bool CONTACT::NitscheStrategyFpi::CheckNitscheContactState(
    CONTACT::Element* cele,                 // the contact element
    const CORE::LINALG::Matrix<2, 1>& xsi,  // local coord on the ele element
    const double& full_fsi_traction,        // stressfluid + penalty
    double& gap                             // gap
)
{
  return CONTACT::UTILS::CheckNitscheContactState(
      *ContactInterfaces()[0], pen_n_, weighting_, cele, xsi, full_fsi_traction, gap);
}

BACI_NAMESPACE_CLOSE
