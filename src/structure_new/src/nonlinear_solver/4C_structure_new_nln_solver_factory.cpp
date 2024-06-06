/*-----------------------------------------------------------*/
/*! \file

\brief Factory for nonlinear solvers in structural dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_nln_solver_factory.hpp"

#include "4C_structure_new_nln_solver_fullnewton.hpp"
#include "4C_structure_new_nln_solver_nox.hpp"
#include "4C_structure_new_nln_solver_ptc.hpp"
#include "4C_structure_new_nln_solver_singlestep.hpp"
#include "4C_structure_new_nln_solver_uzawa.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Nln::SOLVER::Factory::Factory()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::Nln::SOLVER::Generic> STR::Nln::SOLVER::Factory::build_nln_solver(
    const enum Inpar::STR::NonlinSolTech& nlnSolType) const
{
  Teuchos::RCP<STR::Nln::SOLVER::Generic> nlnSolver = Teuchos::null;

  switch (nlnSolType)
  {
    case Inpar::STR::soltech_newtonfull:
      nlnSolver = Teuchos::rcp(new STR::Nln::SOLVER::FullNewton());
      break;
    case Inpar::STR::soltech_nox_nln:
      nlnSolver = Teuchos::rcp(new STR::Nln::SOLVER::Nox());
      break;
    case Inpar::STR::soltech_ptc:
      nlnSolver = Teuchos::rcp(new STR::Nln::SOLVER::PseudoTransient());
      break;
    case Inpar::STR::soltech_singlestep:
      nlnSolver = Teuchos::rcp(new STR::Nln::SOLVER::SingleStep());
      break;
    case Inpar::STR::soltech_newtonuzawanonlin:
    case Inpar::STR::soltech_newtonuzawalin:
      //      nlnSolver = Teuchos::rcp(new STR::Nln::SOLVER::Uzawa());
      //      break;
    default:
      FOUR_C_THROW("Solution technique \"%s\" is not implemented.",
          Inpar::STR::NonlinSolTechString(nlnSolType).c_str());
      break;
  }

  return nlnSolver;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::Nln::SOLVER::Generic> STR::Nln::SOLVER::build_nln_solver(
    const enum Inpar::STR::NonlinSolTech& nlnSolType)
{
  Factory factory;
  return factory.build_nln_solver(nlnSolType);
}

FOUR_C_NAMESPACE_CLOSE
