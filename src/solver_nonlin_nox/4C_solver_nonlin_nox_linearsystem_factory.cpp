/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN factory to create a %::NOX::Epetra::LinearSystem.



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_linearsystem_factory.hpp"

#include "4C_cardiovascular0d_nox_nln_linearsystem.hpp"
#include "4C_constraint_nox_nln_lagpenconstraint_linearsystem.hpp"
#include "4C_contact_aug_nox_nln_contact_linearsystem.hpp"
#include "4C_contact_nox_nln_meshtying_linearsystem.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_solver_nonlin_nox_globaldata.hpp"
#include "4C_structure_new_nox_nln_str_linearsystem.hpp"
#include "4C_utils_exceptions.hpp"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Scaling.H>
#include <NOX_Epetra_Vector.H>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LinSystem::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Epetra::LinearSystem> NOX::NLN::LinSystem::Factory::BuildLinearSystem(
    const NOX::NLN::LinSystem::LinearSystemType& linsystype, NOX::NLN::GlobalData& noxNlnGlobalData,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& jac,
    const Teuchos::RCP<::NOX::Epetra::Vector>& cloneVector,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& precMat,
    const Teuchos::RCP<::NOX::Epetra::Scaling>& scalingObject) const
{
  Teuchos::RCP<::NOX::Epetra::LinearSystem> linSys = Teuchos::null;

  // extract some stuff from the NOX::NLN::GlobalData object
  const NOX::NLN::LinearSystem::SolverMap& linSolvers = noxNlnGlobalData.GetLinSolvers();
  const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq =
      noxNlnGlobalData.get_required_interface();
  const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac =
      noxNlnGlobalData.get_jacobian_interface();
  const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec =
      noxNlnGlobalData.get_preconditioner_interface();

  Teuchos::ParameterList& params = noxNlnGlobalData.GetNlnParameterList();
  // printing parameters
  Teuchos::ParameterList& printParams = params.sublist("Printing", true);
  // linear solver parameters
  Teuchos::ParameterList& lsParams =
      params.sublist("Direction", true).sublist("Newton", true).sublist("Linear Solver", true);

  switch (linsystype)
  {
    // pure structural case
    case NOX::NLN::LinSystem::linear_system_structure:
    {
      linSys = Teuchos::rcp(new NOX::NLN::STR::LinearSystem(printParams, lsParams, linSolvers, iReq,
          iJac, jac, iPrec, precMat, *cloneVector, scalingObject));
      break;
    }
    // structural/contact case
    case NOX::NLN::LinSystem::linear_system_structure_contact:
    {
      const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr =
          noxNlnGlobalData.get_constraint_interfaces();
      const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec =
          noxNlnGlobalData.get_constraint_prec_interfaces();

      linSys = Teuchos::rcp(new NOX::NLN::CONTACT::LinearSystem(printParams, lsParams, linSolvers,
          iReq, iJac, iConstr, jac, iPrec, iConstrPrec, precMat, *cloneVector, scalingObject));
      break;
    }
    case NOX::NLN::LinSystem::linear_system_structure_meshtying:
    {
      const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr =
          noxNlnGlobalData.get_constraint_interfaces();
      const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec =
          noxNlnGlobalData.get_constraint_prec_interfaces();

      linSys = Teuchos::rcp(new NOX::NLN::MESHTYING::LinearSystem(printParams, lsParams, linSolvers,
          iReq, iJac, iConstr, jac, iPrec, iConstrPrec, precMat, *cloneVector, scalingObject));
      break;
    }
    // structural/cardiovascular0d case
    case NOX::NLN::LinSystem::linear_system_structure_cardiovascular0d:
    {
      linSys = Teuchos::rcp(new NOX::NLN::CARDIOVASCULAR0D::LinearSystem(printParams, lsParams,
          linSolvers, iReq, iJac, jac, iPrec, precMat, *cloneVector, scalingObject));
      break;
    }
    // structural/constraint case
    case NOX::NLN::LinSystem::linear_system_structure_lag_pen_constraint:
    {
      const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr =
          noxNlnGlobalData.get_constraint_interfaces();
      const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec =
          noxNlnGlobalData.get_constraint_prec_interfaces();

      linSys = Teuchos::rcp(
          new NOX::NLN::LAGPENCONSTRAINT::LinearSystem(printParams, lsParams, linSolvers, iReq,
              iJac, iConstr, jac, iPrec, iConstrPrec, precMat, *cloneVector, scalingObject));

      break;
    }

    // default case
    default:
    {
      FOUR_C_THROW(
          "ERROR - NOX::NLN::LinSystem::Factory::BuildLinearSystem - "
          "No capable LinearSystem constructor was found (enum = %s|%i)!",
          NOX::NLN::LinSystem::LinearSystemType2String(linsystype).c_str(), linsystype);
      break;
    }
  }  // end switch

  // return the linear system
  return linSys;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Epetra::LinearSystem> NOX::NLN::LinSystem::BuildLinearSystem(
    const NOX::NLN::LinSystem::LinearSystemType& linsystype, NOX::NLN::GlobalData& noxNlnGlobalData,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& jac,
    const Teuchos::RCP<::NOX::Epetra::Vector>& cloneVector,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& precMat,
    const Teuchos::RCP<::NOX::Epetra::Scaling>& scalingObject)
{
  Factory factory;
  return factory.BuildLinearSystem(
      linsystype, noxNlnGlobalData, jac, cloneVector, precMat, scalingObject);
}

FOUR_C_NAMESPACE_CLOSE
