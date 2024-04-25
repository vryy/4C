/*-----------------------------------------------------------*/
/*! \file

\brief This class summarizes the functionality which all
       implicit time integration strategies share and have in
       common.


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_timint_implicitbase.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_structure_new_integrator.hpp"

#include <NOX_Abstract_Group.H>
#include <NOX_Epetra_Vector.H>
#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::ImplicitBase::ImplicitBase()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::TIMINT::ImplicitBase::GetF() const
{
  const ::NOX::Abstract::Group& solgrp = GetSolutionGroup();
  const ::NOX::Epetra::Vector& F = dynamic_cast<const ::NOX::Epetra::Vector&>(solgrp.getF());
  return GetDataGlobalState().ExtractDisplEntries(F.getEpetraVector());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::TIMINT::ImplicitBase::Freact()
{
  CheckInitSetup();
  return DataGlobalState().GetFreactNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> STR::TIMINT::ImplicitBase::SystemMatrix()
{
  CheckInitSetup();
  return Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(DataGlobalState().GetJacobian());
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> STR::TIMINT::ImplicitBase::BlockSystemMatrix()
{
  CheckInitSetup();
  return Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(
      DataGlobalState().GetJacobian());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::ImplicitBase::UseBlockMatrix(
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> domainmaps,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> rangemaps)
{
  FOUR_C_THROW("Currently disabled!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::StcScale STR::TIMINT::ImplicitBase::GetSTCAlgo() { return DataSDyn().GetSTCAlgoType(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> STR::TIMINT::ImplicitBase::GetSTCMat()
{
  FOUR_C_THROW("Not yet implemented!");
  /* See the scaling object in the NOX::NLN::Epetra::LinearSystem class.
   * The STC matrix has to be implemented as a scaling object or as a
   * preconditioner. Both are part of the linear system. */
  // group->linearsystem->scalingobject
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::TIMINT::ImplicitBase::InitialGuess()
{
  CheckInitSetup();
  FOUR_C_THROW("Not yet implemented!");
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
void STR::TIMINT::ImplicitBase::PrintStep()
{
  CheckInitSetup();

  if (DataGlobalState().GetMyRank() != 0 or GroupId() != 0) return;

  const int stepmax = DataSDyn().GetStepMax();
  const int stepn = DataGlobalState().GetStepN();
  const double& timen = DataGlobalState().GetTimeN();
  const double& dt = (*DataGlobalState().GetDeltaTime())[0];
  const int newtoniter = DataGlobalState().GetNlnIterationNumber(stepn);
  double wct = DataGlobalState().GetTimer()->totalElapsedTime(true);

  // open outstd::stringstream
  std::ostringstream oss;

  /* Output of the following quantities
   * time   : total simulated time
   * dt     : used time step
   * nlniter: number of nonlinear solver steps
   * wct    : wall clock time */
  oss << "Finalised step " << std::setw(1) << stepn;
  oss << " / " << std::setw(1) << stepmax;
  oss << " | time " << std::setw(9) << std::setprecision(3) << std::scientific << timen;
  oss << " | dt " << std::setw(9) << std::setprecision(3) << std::scientific << dt;
  oss << " | nlniter " << std::setw(1) << newtoniter;
  oss << " | wct " << std::setw(8) << std::setprecision(2) << std::scientific << wct;
  oss << "\n--------------------------------------------------------------------------------\n";

  // print to ofile (could be done differently...)
  fprintf(stdout, "%s\n", oss.str().c_str());

  // print it, now
  fflush(stdout);
}

FOUR_C_NAMESPACE_CLOSE
