/*-----------------------------------------------------------*/
/*!
\file nox_nln_linearsystem_factory.cpp

\maintainer Michael Hiermeier

\date Jul 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_linearsystem_factory.H"

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparseoperator.H"

#include <Teuchos_ParameterList.hpp>

#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Epetra_Scaling.H>

//// header files for different linearSystems
#include "../drt_contact_aug/nox_nln_contact_linearsystem.H"
#include "../drt_structure_new/nox_nln_str_linearsystem.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LinSystem::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> NOX::NLN::LinSystem::Factory::BuildLinearSystem(
    Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& linSolvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& J,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<LINALG::SparseOperator>& M,
    const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject)
{
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys = Teuchos::null;

  // pure structural case
  if (linSolvers.size()==1 and
      linSolvers.find(NOX::NLN::sol_structure)!=linSolvers.end())
  {
    linSys = Teuchos::rcp(new NOX::NLN::STR::LinearSystem(printParams,linearSolverParams,linSolvers,iReq,iJac,J,iPrec,M,cloneVector,scalingObject));
  }
  // structural/contact case
  else if (linSolvers.size()==2 and
      linSolvers.find(NOX::NLN::sol_structure)!=linSolvers.end() and
      linSolvers.find(NOX::NLN::sol_contact)!=linSolvers.end())
  {
    linSys = Teuchos::rcp(new NOX::NLN::CONTACT::LinearSystem(printParams,linearSolverParams,linSolvers,iReq,iJac,J,iPrec,M,cloneVector,scalingObject));
  }
  else
  {
    dserror("ERROR - NOX::NLN::LinSystem::Factory::BuildLinearSystem - No capable LinearSystem constructor was found!");
  }

  // return the linear system
  return linSys;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> NOX::NLN::LinSystem::BuildLinearSystem(
    Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& linSolvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& J,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<LINALG::SparseOperator>& M,
    const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject)
{
  Factory factory;

  return factory.BuildLinearSystem(printParams,linearSolverParams,
      linSolvers,iReq,iJac,J,iPrec,M,cloneVector,scalingObject);
}
