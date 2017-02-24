/*-----------------------------------------------------------*/
/*!
\file loca_nln_problem.cpp

\maintainer Michael Hiermeier

\date Nov 20, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "loca_nln_problem.H"
#include "loca_nln_statustest_factory.H"
#include "loca_nln_interface_required.H"
#include "loca_nln_group.H"

#include "../solver_nonlin_nox/nox_nln_globaldata.H"

#include "../drt_lib/drt_dserror.H"

#include <Epetra_Vector.h>

#include <LOCA_GlobalData.H>
#include <LOCA_Parameter_Vector.H>
#include <LOCA_Parameter_SublistParser.H>

#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LOCA::NLN::Problem::Problem(
    const Teuchos::RCP<NOX::NLN::GlobalData>& nox_nln_global_data_ptr,
    const Teuchos::RCP<LOCA::GlobalData>& loca_global_data_ptr,
    const Teuchos::RCP<NOX::Epetra::Vector>& x_ptr,
    const Teuchos::RCP<LINALG::SparseOperator>& jac_ptr,
    const Teuchos::RCP<LOCA::ParameterVector>& loca_param_vec_ptr )
    : NOX::NLN::Problem(nox_nln_global_data_ptr,x_ptr,jac_ptr),
      loca_global_data_ptr_(loca_global_data_ptr),
      loca_param_vec_ptr_(loca_param_vec_ptr),
      isLocaStatusTest_(false)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Abstract::Group> LOCA::NLN::Problem::CreateGroup(
    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys) const
{
  CheckInit();

  Teuchos::RCP<NOX::Abstract::Group> locagrp = Teuchos::null;
  Teuchos::ParameterList& pnox = noxNlnGlobalData_->GetNlnParameterList();
  Teuchos::RCP<LOCA::NLN::Interface::Required> ireq =
      Teuchos::rcp_dynamic_cast<LOCA::NLN::Interface::Required>(
          noxNlnGlobalData_->GetRequiredInterface());
  if (ireq.is_null())
    dserror("Dynamic cast to \"LOCA::NLN::Interface::Required\" failed!");

  if (noxNlnGlobalData_->GetIsConstrained())
  {
    dserror("LOCA can not handle constrained problems at the moment!");
  }
  else
  {
    locagrp = Teuchos::rcp(new LOCA::NLN::Group(loca_global_data_ptr_,
        pnox.sublist("Group Options"),pnox.sublist("Printing"),ireq,**xVector_,linsys,
        *loca_param_vec_ptr_));
  }

  return locagrp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::NLN::Problem::CreateStatusTests(
    Teuchos::RCP<NOX::StatusTest::Generic>& outerTest,
    Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTest) const
{
  if (!isLocaStatusTest_)
    dserror("Create the LOCA::StatusTest first!");

  NOX::NLN::Problem::CreateStatusTests(outerTest,innerTest);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::NLN::Problem::CreateStatusTests(
    Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
    Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
    Teuchos::RCP<LOCA::StatusTest::Abstract>& locaTests) const
{
  Teuchos::RCP<Teuchos::ParameterList> p = loca_global_data_ptr_->parsedParams->
      getSublist("Status Test");

  Teuchos::ParameterList& p_loca_st = p->sublist("LOCA Status Test",true);
  locaTests = LOCA::NLN::StatusTest::buildStatusTests(p_loca_st,loca_global_data_ptr_);
  if (!locaTests.is_null())


  CreateStatusTests(outerTests,innerTests);

  return;
}

