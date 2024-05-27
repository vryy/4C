/*---------------------------------------------------------------------*/
/*! \file

\brief Nitsche contact solving strategy for problems with FPI

\level 3


*/
/*---------------------------------------------------------------------*/

#include "4C_contact_nitsche_strategy_fpi.hpp"

#include "4C_contact_interface.hpp"
#include "4C_contact_nitsche_strategy_fsi.hpp"

FOUR_C_NAMESPACE_OPEN


void CONTACT::NitscheStrategyFpi::set_state(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  CONTACT::NitscheStrategyPoro::set_state(statename, vec);
  if (statename == MORTAR::state_new_displacement)
  {
    do_contact_search();
  }
}

void CONTACT::NitscheStrategyFpi::do_contact_search()
{
  for (auto& interface : interface_)
  {
    interface->Initialize();
    interface->evaluate_search_binarytree();
    interface->evaluate_nodal_normals();
    interface->export_nodal_normals();
  }
}

bool CONTACT::NitscheStrategyFpi::check_nitsche_contact_state(
    CONTACT::Element* cele,                 // the contact element
    const CORE::LINALG::Matrix<2, 1>& xsi,  // local coord on the ele element
    const double& full_fsi_traction,        // stressfluid + penalty
    double& gap                             // gap
)
{
  return CONTACT::UTILS::check_nitsche_contact_state(
      *ContactInterfaces()[0], pen_n_, weighting_, cele, xsi, full_fsi_traction, gap);
}

FOUR_C_NAMESPACE_CLOSE
