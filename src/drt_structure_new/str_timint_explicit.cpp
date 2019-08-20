/*-----------------------------------------------------------*/
/*! \file

\brief explicit structural time integration

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_timint_explicit.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Explicit::Explicit() : STR::TIMINT::Base()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::Setup()
{
  // safety check
  CheckInit();

  Base::Setup();

  // set isSetup flag
  issetup_ = true;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::PrepareTimeStep()
{
  CheckInitSetup();
  // things that need to be done before Predict
  PrePredict();

  // ToDo prepare contact for new time step
  // PrepareStepContact();

  // things that need to be done after Predict
  PostPredict();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  CheckInitSetup();
  dserror(
      "All monolithically coupled problems work with implicit time "
      "integration schemes. Thus, calling Evaluate() in an explicit scheme "
      "is not possible.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::SetState(const Teuchos::RCP<Epetra_Vector>& x)
{
  dserror(
      "All coupled problems work with implicit time "
      "integration schemes. Thus, calling SetState() in an explicit scheme "
      "is not considered, yet.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::ConvergenceStatus STR::TIMINT::Explicit::Solve()
{
  CheckInitSetup();
  IntegrateStep();
  return INPAR::STR::conv_success;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::PreparePartitionStep()
{
  // do nothing for explicit time integrators
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::Update(double endtime)
{
  CheckInitSetup();
  dserror("Not implemented. No time adaptivity available for explicit time integration.");
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::PrintStep()
{
  CheckInitSetup();
  // FixMe
  if (DataGlobalState().GetMyRank() == 0)
    std::cout << "FixMe: The PrintStep() routine is not yet implemented!" << std::endl;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::STC_Scale STR::TIMINT::Explicit::GetSTCAlgo()
{
  CheckInitSetup();
  dserror("GetSTCAlgo() has not been tested for explicit time integration.");
  return INPAR::STR::stc_none;
};


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> STR::TIMINT::Explicit::GetSTCMat()
{
  CheckInitSetup();
  dserror("GetSTCMat() has not been tested for explicit time integration.");
  return Teuchos::null;
};


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Explicit::Integrate()
{
  dserror(
      "The function is unused since the ADAPTER::StructureTimeLoop "
      "wrapper gives you all the flexibility you need.");
  return 0;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Explicit::IntegrateStep()
{
  CheckInitSetup();
  dserror("Not yet implemented!");
  return 0;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::TIMINT::Explicit::InitialGuess()
{
  dserror("InitialGuess() is not available for explicit time integration");
  return Teuchos::null;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::TIMINT::Explicit::GetF() const
{
  dserror("RHS() is not available for explicit time integration");
  return Teuchos::null;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::TIMINT::Explicit::Freact()
{
  CheckInitSetup();
  dserror("Not implemented!");
  return Teuchos::null;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> STR::TIMINT::Explicit::SystemMatrix()
{
  dserror("SystemMatrix() is not available for explicit time integration");
  return Teuchos::null;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> STR::TIMINT::Explicit::BlockSystemMatrix()
{
  dserror("BlockSystemMatrix() is not available for explicit time integration");
  return Teuchos::null;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
    Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps)
{
  dserror("UseBlockMatrix() is not available for explicit time integration");
}
///@}
