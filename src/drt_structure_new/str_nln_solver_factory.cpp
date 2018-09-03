/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_factory.cpp

\maintainer Michael Hiermeier

\date Sep 2, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_nln_solver_factory.H"

// supported nonlinear solvers
#include "str_nln_solver_fullnewton.H"
#include "str_nln_solver_nonlin.H"
#include "str_nln_solver_nox.H"
#include "str_nln_solver_ptc.H"
#include "str_nln_solver_uzawa.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::NLN::SOLVER::Factory::Factory()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::NLN::SOLVER::Generic> STR::NLN::SOLVER::Factory::BuildNlnSolver(
    const enum INPAR::STR::NonlinSolTech& nlnSolType) const
{
  Teuchos::RCP<STR::NLN::SOLVER::Generic> nlnSolver = Teuchos::null;

  switch (nlnSolType)
  {
    case INPAR::STR::soltech_newtonfull:
      nlnSolver = Teuchos::rcp(new STR::NLN::SOLVER::FullNewton());
      break;
    case INPAR::STR::soltech_nox_nln:
      nlnSolver = Teuchos::rcp(new STR::NLN::SOLVER::Nox());
      break;
    case INPAR::STR::soltech_ptc:
      nlnSolver = Teuchos::rcp(new STR::NLN::SOLVER::PseudoTransient());
      break;
    case INPAR::STR::soltech_nlnsol:
      //      nlnSolver = Teuchos::rcp(new STR::NLN::SOLVER::Nonlin());
      //      break;
    case INPAR::STR::soltech_newtonuzawanonlin:
    case INPAR::STR::soltech_newtonuzawalin:
      //      nlnSolver = Teuchos::rcp(new STR::NLN::SOLVER::Uzawa());
      //      break;
    default:
      dserror("Solution technique \"%s\" is not implemented.",
          INPAR::STR::NonlinSolTechString(nlnSolType).c_str());
      break;
  }

  return nlnSolver;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::NLN::SOLVER::Generic> STR::NLN::SOLVER::BuildNlnSolver(
    const enum INPAR::STR::NonlinSolTech& nlnSolType)
{
  Factory factory;
  return factory.BuildNlnSolver(nlnSolType);
}
