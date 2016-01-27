/*-----------------------------------------------------------*/
/*!
\file str_timint_implicitbase.cpp

\maintainer Michael Hiermeier

\date Dec 16, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_timint_implicitbase.H"
#include "str_integrator.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_blocksparsematrix.H"

#include <NOX_Epetra_Vector.H>
#include <NOX_Abstract_Group.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::ImplicitBase::ImplicitBase()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::TIMINT::ImplicitBase::RHS()
{
  const NOX::Abstract::Group& solgrp = GetSolutionGroup();
  const NOX::Epetra::Vector& F =
      dynamic_cast<const NOX::Epetra::Vector&>(solgrp.getF());
  return DataGlobalState().ExportDisplEntries(F.getEpetraVector());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::TIMINT::ImplicitBase::Freact()
{
  CheckInitSetup();
  return DataGlobalState().GetMutableFreactNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> STR::TIMINT::ImplicitBase::SystemMatrix()
{
  CheckInitSetup();
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>
      (DataGlobalState().GetMutableJacobian());
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> STR::TIMINT::ImplicitBase::BlockSystemMatrix()
{
  CheckInitSetup();
  return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>
      (DataGlobalState().GetMutableJacobian());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::ImplicitBase::UseBlockMatrix(
        Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
        Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps)
{
  dserror("Currently disabled!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::STC_Scale STR::TIMINT::ImplicitBase::GetSTCAlgo()
{
  return DataSDyn().GetSTCAlgoType();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> STR::TIMINT::ImplicitBase::GetSTCMat()
{
  dserror("Not yet implemented!");
  /* See the scaling object in the NOX::NLN::Epetra::LinearSystem class.
   * The STC matrix has to be implemented as a scaling object or as a
   * preconditioner. Both are part of the linear system. */
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::ImplicitBase::Evaluate(
    Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  CheckInitSetup();
  dserror("Not yet implemented!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::TIMINT::ImplicitBase::InitialGuess()
{
  CheckInitSetup();
  dserror("Not yet imlemented!");
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::ImplicitBase::Update(double endtime)
{
  CheckInitSetup();
  PreUpdate();
  Integrator().UpdateStepState();
  SetTimeNp(endtime);
  UpdateStepTime();
  Integrator().UpdateStepElement();
  PostUpdate();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::ImplicitBase::Output(bool forced_writerestart)
{
  CheckInitSetup();
  Integrator().OutputStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::ImplicitBase::PrintStep()
{
  CheckInitSetup();
  // FixMe
  if (DataGlobalState().GetMyRank() == 0)
    std::cout << "FixMe: The PrintStep() routine is not yet implemented!" << std::endl;
}
