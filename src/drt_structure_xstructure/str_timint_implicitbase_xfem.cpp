/*----------------------------------------------------------------------------*/
/*!

\brief This file contains the adaptions of the implicit structural time
       integration due to XFEM

\level 2

\maintainer Anh-Tu Vuong

*/
/*----------------------------------------------------------------------------*/

#include "../drt_structure_new/str_timint_implicitbase.H"
#include "../solver_nonlin_nox/nox_nln_group.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::ImplicitBase::DestroyNoxState()
{
  Teuchos::RCP<NOX::Abstract::Group> solgrp = SolutionGroupPtr();
  Teuchos::RCP<NOX::NLN::Group> nlngrp = Teuchos::rcp_dynamic_cast<NOX::NLN::Group>(solgrp, true);

  nlngrp->DestroyJacobian();

  return true;
}
