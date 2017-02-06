/*-----------------------------------------------------------*/
/*!
\file nox_nln_linearsystem_factory.cpp

\brief %NOX::NLN factory to create a %NOX::Epetra::LinearSystem.

\maintainer Michael Hiermeier

\date Jul 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_linearsystem_factory.H"
#include "nox_nln_globaldata.H"

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
#include "../drt_cardiovascular0d/nox_nln_cardiovascular0d_linearsystem.H"
#include "../drt_constraint/nox_nln_lagpenconstraint_linearsystem.H"
#include "../drt_structure_new/nox_nln_str_linearsystem.H"
#include "../drt_scatra/nox_nln_scatra_linearsystem.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LinSystem::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> NOX::NLN::LinSystem::Factory::
    BuildLinearSystem(
    const NOX::NLN::LinSystem::LinearSystemType& linsystype,
    NOX::NLN::GlobalData& noxNlnGlobalData,
    const Teuchos::RCP<LINALG::SparseOperator>& jac,
    const Teuchos::RCP<NOX::Epetra::Vector>& cloneVector,
    const Teuchos::RCP<LINALG::SparseOperator>& precMat,
    const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject) const
{
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys = Teuchos::null;

  // extract some stuff from the NOX::NLN::GlobalData object
  const NOX::NLN::LinearSystem::SolverMap& linSolvers =
      noxNlnGlobalData.GetLinSolvers();
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq =
      noxNlnGlobalData.GetRequiredInterface();
  const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac =
      noxNlnGlobalData.GetJacobianInterface();
  const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec =
      noxNlnGlobalData.GetPreconditionerInterface();

  Teuchos::ParameterList& params = noxNlnGlobalData.GetNlnParameterList();
  // printing parameters
  Teuchos::ParameterList& printParams = params.sublist("Printing",true);
  // linear solver parameters
  Teuchos::ParameterList& lsParams    = params.sublist("Direction",true).
      sublist("Newton",true).sublist("Linear Solver",true);

  switch (linsystype)
  {
    // pure structural case
    case NOX::NLN::LinSystem::linear_system_structure:
    {
      linSys = Teuchos::rcp(new NOX::NLN::STR::LinearSystem(
          printParams,lsParams,linSolvers,iReq,iJac,jac,iPrec,precMat,
          *cloneVector,scalingObject));
      break;
    }
    // structural/contact case
    case NOX::NLN::LinSystem::linear_system_structure_contact:
    {
      const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr =
          noxNlnGlobalData.GetConstraintInterfaces();
      const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec =
          noxNlnGlobalData.GetConstraintPrecInterfaces();

      linSys = Teuchos::rcp(new NOX::NLN::CONTACT::LinearSystem(
          printParams,lsParams,linSolvers,iReq,iJac,iConstr,jac,iPrec,
          iConstrPrec,precMat,*cloneVector,scalingObject));
      break;
    }
    // structural/cardiovascular0d case
    case NOX::NLN::LinSystem::linear_system_structure_cardiovascular0d:
    {
      linSys = Teuchos::rcp(new NOX::NLN::CARDIOVASCULAR0D::LinearSystem(
          printParams,lsParams,linSolvers,iReq,iJac,jac,iPrec,precMat,
          *cloneVector,scalingObject));
      break;
    }
    // structural/constraint case
    case NOX::NLN::LinSystem::linear_system_structure_lag_pen_constraint:
    {
      const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr =
          noxNlnGlobalData.GetConstraintInterfaces();
      const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec =
          noxNlnGlobalData.GetConstraintPrecInterfaces();

      linSys = Teuchos::rcp(new NOX::NLN::LAGPENCONSTRAINT::LinearSystem(
          printParams,lsParams,linSolvers,iReq,iJac,iConstr,jac,iPrec,
          iConstrPrec,precMat,*cloneVector,scalingObject));

      break;
    }
    // scalar transport case
    case NOX::NLN::LinSystem::linear_system_scatra:
    {
      linSys = Teuchos::rcp(new NOX::NLN::SCATRA::LinearSystem(
          printParams,lsParams,linSolvers,iReq,iJac,jac,iPrec,precMat,
          *cloneVector,scalingObject));
      break;
    }

    // default case
    default:
    {
      dserror("ERROR - NOX::NLN::LinSystem::Factory::BuildLinearSystem - "
          "No capable LinearSystem constructor was found (enum = %s|%i)!",
          NOX::NLN::LinSystem::LinearSystemType2String(linsystype).c_str(),
          linsystype);
      break;
    }
  } // end switch

  // return the linear system
  return linSys;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> NOX::NLN::LinSystem::BuildLinearSystem(
    const NOX::NLN::LinSystem::LinearSystemType& linsystype,
    NOX::NLN::GlobalData& noxNlnGlobalData,
    const Teuchos::RCP<LINALG::SparseOperator>& jac,
    const Teuchos::RCP<NOX::Epetra::Vector>& cloneVector,
    const Teuchos::RCP<LINALG::SparseOperator>& precMat,
    const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject
    )
{
  Factory factory;
  return factory.BuildLinearSystem(linsystype,noxNlnGlobalData,jac,cloneVector,
      precMat,scalingObject);
}
