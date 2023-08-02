/*-----------------------------------------------------------*/
/*! \file

\brief Extension of Trilinos's NOX::Solver::SingleStep to include Baci's specific inner/outer
tests

\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_solver_nonlin_nox_solver_singlestep.H"  // class definition

#include "baci_solver_nonlin_nox_aux.H"
#include "baci_solver_nonlin_nox_group.H"
#include "baci_solver_nonlin_nox_inner_statustest_generic.H"

// templated status tests
#include "baci_solver_nonlin_nox_statustest_activeset.H"
#include "baci_solver_nonlin_nox_statustest_normf.H"
#include "baci_solver_nonlin_nox_statustest_normupdate.H"
#include "baci_solver_nonlin_nox_statustest_normwrms.H"

#include <NOX_Abstract_Group.H>
#include <NOX_Solver_SolverUtils.H>
#include <NOX_Utils.H>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Solver::SingleStep::SingleStep(const Teuchos::RCP<NOX::Abstract::Group>& grp,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
    const Teuchos::RCP<Teuchos::ParameterList>& params)
    : NOX::Solver::SingleStep(grp, params)
{
  // call own init() after base init() was called.
  init(innerTests);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::SingleStep::init(
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests)
{
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::StatusTest::StatusType NOX::NLN::Solver::SingleStep::getStatus() const { return status; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::Utils& NOX::NLN::Solver::SingleStep::GetUtils() const { return *utilsPtr; }
