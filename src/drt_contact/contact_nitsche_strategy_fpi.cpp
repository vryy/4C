/*---------------------------------------------------------------------*/
/*! \file

\brief Nitsche contact solving strategy for problems with FPI

\level 3

\maintainer Christoph Ager

*/
/*---------------------------------------------------------------------*/

#include "contact_nitsche_strategy_fpi.H"
#include "contact_nitsche_strategy_fsi.H"

//#include "contact_nitsche_integrator_fpi.H"
//#include "contact_element.H"
#include "contact_interface.H"

//#include "../drt_mortar/mortar_projector.H"

//#include "../drt_lib/drt_discret.H"

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
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->Initialize();
    interface_[i]->EvaluateSearchBinarytree();
    interface_[i]->EvaluateNodalNormals();
    interface_[i]->ExportNodalNormals();
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
